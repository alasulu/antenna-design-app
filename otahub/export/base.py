"""A small geometry representation shared by every simulator backend.

Archetype specs carry electrical design formulas, not solid models, so an
exporter has to know how to turn primary dimensions into shapes. That
knowledge lives here once, in a neutral form, and each backend renders it.
Writing it twice - once in VBA and once in IronPython - would guarantee the
two drift apart.

Honesty rule: an archetype with no builder exports its PARAMETERS and says
plainly that no solid geometry was generated. A half-built model that looks
complete is worse than none.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable

from ..core.archetype import DesignResult


@dataclass(frozen=True, slots=True)
class Brick:
    name: str
    material: str
    x: tuple[float, float]
    y: tuple[float, float]
    z: tuple[float, float]


@dataclass(frozen=True, slots=True)
class Cylinder:
    name: str
    material: str
    axis: str                     # "x" | "y" | "z"
    radius: float
    span: tuple[float, float]     # start, end along `axis`
    centre: tuple[float, float] = (0.0, 0.0)   # the other two coordinates


@dataclass(frozen=True, slots=True)
class Cone:
    """Truncated cone. Set one radius to zero for a true cone."""
    name: str
    material: str
    axis: str                     # "x" | "y" | "z"
    radius_start: float           # at span[0]
    radius_end: float             # at span[1]
    span: tuple[float, float]
    centre: tuple[float, float] = (0.0, 0.0)


@dataclass(frozen=True, slots=True)
class Torus:
    """A ring. `major_radius` is axis-to-tube-centre, `minor_radius` the tube.

    Both backends have this as a primitive, which is why the loop family can be
    built at all - unlike the horns, which would need a loft.
    """
    name: str
    material: str
    axis: str                     # normal to the plane of the ring
    major_radius: float
    minor_radius: float
    centre: tuple[float, float, float] = (0.0, 0.0, 0.0)


@dataclass(frozen=True, slots=True)
class Sphere:
    name: str
    material: str
    radius: float
    centre: tuple[float, float, float] = (0.0, 0.0, 0.0)


@dataclass(frozen=True, slots=True)
class Subtract:
    """Boolean difference, applied after every solid exists.

    Kept out of the solid list deliberately: it is an operation on named
    solids, not a shape, and both backends want it emitted last.
    """
    target: str
    tools: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class DiscretePort:
    """A lumped/discrete port between two points, for wire and probe feeds."""
    name: str
    start: tuple[float, float, float]
    end: tuple[float, float, float]
    impedance: float = 50.0


@dataclass(slots=True)
class Model:
    """A complete exportable model."""
    archetype: str
    title: str
    parameters: dict[str, float]              # name -> value, SI
    #: name -> the unit the spec declares for it. Exporters need this to know
    #: what is a length: guessing from magnitude wrote a 319 ohm resistance
    #: out as "319105 mm".
    units: dict[str, str] = field(default_factory=dict)
    solids: list = field(default_factory=list)
    #: Boolean operations, applied in order after the solids are created.
    operations: list = field(default_factory=list)
    ports: list = field(default_factory=list)
    frequency_hz: float = 1e9
    notes: list[str] = field(default_factory=list)
    built_geometry: bool = True

    @property
    def has_geometry(self) -> bool:
        return bool(self.solids)


#: archetype key -> builder(design) -> Model
BUILDERS: dict[str, Callable[[DesignResult], Model]] = {}


def builder(key: str):
    """Register a geometry builder for an archetype."""
    def wrap(fn):
        BUILDERS[key] = fn
        return fn
    return wrap


def _param(design: DesignResult, *names: str, default: float | None = None) -> float:
    for name in names:
        value = design.get(name)
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            return float(value)
    if default is not None:
        return default
    raise KeyError(f"design {design.archetype!r} has none of {names}")


def _base_model(design: DesignResult, title: str) -> Model:
    """Collect every numeric value the design involved.

    Requirements are included alongside derived parameters: a substrate
    thickness the user supplied is just as much part of the exported model as
    a patch length the engine computed, and leaving it out of the simulator's
    parameter list would make the geometry undrivable.
    """
    merged = {**design.requirements, **design.parameters}
    numeric = {k: float(v) for k, v in merged.items()
               if isinstance(v, (int, float)) and not isinstance(v, bool)
               and k not in ("k0",)}
    units = dict(design.units)
    units.setdefault("f0", "Hz")
    return Model(archetype=design.archetype, title=title, parameters=numeric,
                 units=units,
                 frequency_hz=float(design.requirements.get("f0", 1e9)))


# --------------------------------------------------------------- builders

@builder("half_wave_dipole")
@builder("resonant_dipole")
@builder("short_dipole")
@builder("dipole_arbitrary_length")
def _dipole(design: DesignResult) -> Model:
    """Two collinear wire arms along z with a discrete port in the gap."""
    model = _base_model(design, "Centre-fed dipole")
    total = _param(design, "L")
    radius = _param(design, "aw", default=total / 2000.0)
    gap = max(total / 200.0, radius * 2.0)
    half = total / 2.0
    model.solids += [
        Cylinder("arm_upper", "PEC", "z", radius, (gap / 2.0, half)),
        Cylinder("arm_lower", "PEC", "z", radius, (-half, -gap / 2.0)),
    ]
    model.ports.append(DiscretePort("port1", (0.0, 0.0, -gap / 2.0), (0.0, 0.0, gap / 2.0)))
    model.notes += [
        f"Feed gap set to {gap * 1e3:.4g} mm; the spec does not define one, and the "
        "gap materially affects the computed input reactance.",
        "Wire modelled as a solid PEC cylinder. Thin-wire theory assumes a "
        "perfect conductor with no end caps.",
    ]
    return model


@builder("quarter_wave_monopole")
def _monopole(design: DesignResult) -> Model:
    model = _base_model(design, "Quarter-wave monopole over a ground plane")
    height = _param(design, "h")
    radius = _param(design, "aw", default=height / 1000.0)
    gap = max(height / 100.0, radius * 2.0)
    ground = 4.0 * height
    model.solids += [
        Cylinder("monopole", "PEC", "z", radius, (gap, height + gap)),
        Cylinder("ground_plane", "PEC", "z", ground / 2.0, (-ground / 200.0, 0.0)),
    ]
    model.ports.append(DiscretePort("port1", (0.0, 0.0, 0.0), (0.0, 0.0, gap)))
    model.notes += [
        f"Ground plane rendered as a finite disc {ground * 1e3:.4g} mm across "
        "(2 wavelengths). The spec's impedance assumes an INFINITE plane; a "
        "finite one raises the input resistance and tilts the pattern upward.",
        f"Feed gap {gap * 1e3:.4g} mm between the plane and the monopole base.",
    ]
    return model


@builder("rectangular_patch")
@builder("rectangular_patch_inset")
def _patch(design: DesignResult) -> Model:
    """Substrate, ground, patch, and an inset notch where one is designed."""
    model = _base_model(design, "Rectangular microstrip patch")
    w = _param(design, "W")
    length = _param(design, "L")
    h = _param(design, "h")
    eps_r = _param(design, "eps_r", default=4.4)
    margin = max(w, length) * 0.6
    sub_w, sub_l = w + 2 * margin, length + 2 * margin
    model.solids += [
        Brick("substrate", f"eps_r={eps_r:g}", (-sub_w / 2, sub_w / 2),
              (-sub_l / 2, sub_l / 2), (0.0, h)),
        Brick("ground", "PEC", (-sub_w / 2, sub_w / 2), (-sub_l / 2, sub_l / 2),
              (0.0, 0.0)),
        Brick("patch", "PEC", (-w / 2, w / 2), (-length / 2, length / 2), (h, h)),
    ]
    y0 = design.get("y0")
    if isinstance(y0, (int, float)) and y0 > 0:
        feed_w = h * 2.0
        model.solids.append(
            Brick("inset_notch", "VOID", (-feed_w * 1.5, feed_w * 1.5),
                  (-length / 2, -length / 2 + float(y0)), (h, h)))
        model.notes.append(
            f"Inset notch cut {float(y0) * 1e3:.4g} mm deep. Notch WIDTH is not "
            "given by the transmission-line model and is set here to twice the "
            "substrate thickness; it affects the achieved impedance and should "
            "be tuned in the solver.")
    model.ports.append(DiscretePort(
        "port1", (0.0, -length / 2 + float(y0 or 0.0), 0.0),
        (0.0, -length / 2 + float(y0 or 0.0), h)))
    model.notes += [
        f"Substrate extended {margin * 1e3:.4g} mm beyond the patch on each side; "
        "the model assumes an infinite substrate and ground, so a finite board "
        "shifts resonance slightly and raises back radiation.",
        "Feed rendered as a vertical probe. An edge-fed or aperture-coupled "
        "design needs a different feed structure than this.",
    ]
    return model


@builder("circular_patch")
def _circular_patch(design: DesignResult) -> Model:
    model = _base_model(design, "Circular microstrip patch")
    a = _param(design, "a")
    h = _param(design, "h")
    eps_r = _param(design, "eps_r", default=2.2)
    sub = a * 3.0
    model.solids += [
        Brick("substrate", f"eps_r={eps_r:g}", (-sub, sub), (-sub, sub), (0.0, h)),
        Brick("ground", "PEC", (-sub, sub), (-sub, sub), (0.0, 0.0)),
        Cylinder("patch", "PEC", "z", a, (h, h)),
    ]
    feed_r = a * 0.3
    model.ports.append(DiscretePort("port1", (feed_r, 0.0, 0.0), (feed_r, 0.0, h)))
    model.notes += [
        f"Probe placed at r = {feed_r * 1e3:.4g} mm (0.3a) as a starting point. "
        "The spec gives no feed position for the circular patch, so this must "
        "be tuned for the target input impedance.",
    ]
    return model


@builder("open_ended_waveguide")
def _oewg(design: DesignResult) -> Model:
    model = _base_model(design, "Open-ended rectangular waveguide")
    a = _param(design, "a_wg")
    b = _param(design, "b_wg")
    length = 2.0 * a
    model.solids.append(
        Brick("guide_interior", "VACUUM", (-a / 2, a / 2), (-b / 2, b / 2), (0.0, length)))
    model.notes += [
        f"Interior volume only, {length * 1e3:.4g} mm long. Assign a waveguide "
        "port to the z = 0 face and radiation boundaries around the aperture.",
        "Walls are implied by the surrounding PEC boundary, not modelled as solids.",
    ]
    return model


@builder("folded_dipole")
def _folded_dipole(design: DesignResult) -> Model:
    """Two parallel arms joined at the ends, fed at the centre of one."""
    model = _base_model(design, "Folded dipole")
    total = _param(design, "L")
    sep = _param(design, "d_sep", default=total / 50.0)
    radius = _param(design, "aw", default=total / 2000.0)
    gap = max(total / 200.0, radius * 2.0)
    half = total / 2.0
    y_fed, y_par = -sep / 2.0, sep / 2.0
    model.solids += [
        Cylinder("fed_upper", "PEC", "z", radius, (gap / 2.0, half), (0.0, y_fed)),
        Cylinder("fed_lower", "PEC", "z", radius, (-half, -gap / 2.0), (0.0, y_fed)),
        Cylinder("parasitic", "PEC", "z", radius, (-half, half), (0.0, y_par)),
        Cylinder("end_top", "PEC", "y", radius, (y_fed, y_par), (0.0, half)),
        Cylinder("end_bottom", "PEC", "y", radius, (y_fed, y_par), (0.0, -half)),
    ]
    model.ports.append(DiscretePort(
        "port1", (0.0, y_fed, -gap / 2.0), (0.0, y_fed, gap / 2.0)))
    model.notes += [
        f"Conductor spacing {sep * 1e3:.4g} mm, from the spec. The N^2 impedance "
        "step assumes EQUAL conductor diameters; unequal ones give a step "
        "between 1 and N^2, which this model does not offer.",
        f"Feed gap {gap * 1e3:.4g} mm, chosen here - the spec does not define one.",
        "End bars modelled as straight cylinders. A real folded dipole has "
        "rounded ends, worth a fraction of a percent in resonant length.",
    ]
    return model


@builder("dipole_over_ground")
def _dipole_over_ground(design: DesignResult) -> Model:
    """Horizontal dipole along x, at height h above a finite ground plane."""
    model = _base_model(design, "Horizontal dipole over a ground plane")
    total = _param(design, "L")
    height = _param(design, "h")
    radius = _param(design, "aw", default=total / 2000.0)
    gap = max(total / 200.0, radius * 2.0)
    half = total / 2.0
    ground = max(4.0 * total, 4.0 * height)
    model.solids += [
        Cylinder("arm_plus", "PEC", "x", radius, (gap / 2.0, half), (0.0, height)),
        Cylinder("arm_minus", "PEC", "x", radius, (-half, -gap / 2.0), (0.0, height)),
        Brick("ground_plane", "PEC", (-ground / 2, ground / 2),
              (-ground / 2, ground / 2), (0.0, 0.0)),
    ]
    model.ports.append(DiscretePort(
        "port1", (-gap / 2.0, 0.0, height), (gap / 2.0, 0.0, height)))
    model.notes += [
        f"Ground plane rendered {ground * 1e3:.4g} mm square. The spec's image "
        "theory assumes an INFINITE perfect conductor; a finite plane raises "
        "low-angle radiation and shifts the feed resistance.",
        f"Feed gap {gap * 1e3:.4g} mm at the dipole centre.",
    ]
    return model


@builder("turnstile_dipole")
def _turnstile(design: DesignResult) -> Model:
    """Crossed dipoles along x and y, each with its own port."""
    model = _base_model(design, "Turnstile: crossed dipoles in quadrature")
    total = _param(design, "L")
    radius = _param(design, "aw", default=total / 2000.0)
    gap = max(total / 200.0, radius * 2.0)
    half = total / 2.0
    offset = radius * 4.0
    model.solids += [
        Cylinder("x_arm_plus", "PEC", "x", radius, (gap / 2.0, half), (0.0, -offset)),
        Cylinder("x_arm_minus", "PEC", "x", radius, (-half, -gap / 2.0), (0.0, -offset)),
        Cylinder("y_arm_plus", "PEC", "y", radius, (gap / 2.0, half), (0.0, offset)),
        Cylinder("y_arm_minus", "PEC", "y", radius, (-half, -gap / 2.0), (0.0, offset)),
    ]
    model.ports += [
        DiscretePort("port1", (-gap / 2.0, 0.0, -offset), (gap / 2.0, 0.0, -offset)),
        DiscretePort("port2", (0.0, -gap / 2.0, offset), (0.0, gap / 2.0, offset)),
    ]
    model.notes += [
        f"The two dipoles are offset {2 * offset * 1e3:.4g} mm in z so they do "
        "not intersect. A real turnstile crosses them in one plane with a "
        "bonded or insulated junction; the offset perturbs the axial ratio "
        "slightly and is a modelling compromise, not part of the design.",
        "TWO ports. Excite them with equal amplitude and 90 degrees of phase "
        "difference; the spec's quarter-wave line is not modelled, because "
        "driving the ports directly is both easier and more accurate.",
        f"Feed gap {gap * 1e3:.4g} mm on each dipole.",
    ]
    return model


@builder("quarter_wave_shorted_patch")
def _shorted_patch(design: DesignResult) -> Model:
    model = _base_model(design, "Quarter-wave shorted patch")
    w = _param(design, "W")
    length = _param(design, "L")
    h = _param(design, "h")
    eps_r = _param(design, "eps_r", default=4.4)
    margin = max(w, length) * 0.6
    sub_w, sub_l = w + 2 * margin, length + 2 * margin
    model.solids += [
        Brick("substrate", f"eps_r={eps_r:g}", (-sub_w / 2, sub_w / 2),
              (-sub_l / 2, sub_l / 2), (0.0, h)),
        Brick("ground", "PEC", (-sub_w / 2, sub_w / 2), (-sub_l / 2, sub_l / 2),
              (0.0, 0.0)),
        Brick("patch", "PEC", (-w / 2, w / 2), (-length / 2, length / 2), (h, h)),
        Brick("shorting_wall", "PEC", (-w / 2, w / 2),
              (-length / 2, -length / 2), (0.0, h)),
    ]
    feed_y = -length / 2 + length * 0.3
    model.ports.append(DiscretePort("port1", (0.0, feed_y, 0.0), (0.0, feed_y, h)))
    model.notes += [
        "Shorting wall spans the FULL patch width, which is what the spec's "
        "resonant length assumes. Narrowing it toward a single pin lowers the "
        "resonance and turns this into a PIFA.",
        f"Probe placed {length * 0.3 * 1e3:.4g} mm from the shorted edge (0.3L) "
        "as a starting point; the spec gives no feed position, so tune it.",
    ]
    return model


@builder("rectangular_dra")
def _rect_dra(design: DesignResult) -> Model:
    model = _base_model(design, "Rectangular dielectric resonator antenna")
    w = _param(design, "w")
    length = _param(design, "Lr")
    d = _param(design, "d")
    eps_r = _param(design, "eps_r", default=10.0)
    ground = max(w, length) * 4.0
    probe_r = min(w, length) / 40.0
    model.solids += [
        Brick("ground", "PEC", (-ground / 2, ground / 2), (-ground / 2, ground / 2),
              (0.0, 0.0)),
        Brick("resonator", f"eps_r={eps_r:g}", (-w / 2, w / 2),
              (-length / 2, length / 2), (0.0, d)),
        Cylinder("probe", "PEC", "z", probe_r, (0.0, d * 0.6),
                 (w / 2 + probe_r * 2, 0.0)),
    ]
    model.ports.append(DiscretePort(
        "port1", (w / 2 + probe_r * 2, 0.0, -d * 0.05),
        (w / 2 + probe_r * 2, 0.0, 0.0)))
    model.notes += [
        "The probe stands at the centre of the face normal to x, where the mode "
        "whose magnetic dipole runs along the length Lr (y) has its vertical E "
        "field - the mode the spec designs. The mode along w is not excited "
        "from there. Slot and microstrip feeds behave differently and are not "
        "modelled; probe height and offset both need tuning in the solver.",
        "The spec's resonance is fitted to FDTD ringdowns of the isolated brick "
        "on an infinite ground (within about 0.2%), so a shift in the solver "
        "comes from the probe and the finite ground plane, not the sizing.",
        f"Ground plane {ground * 1e3:.4g} mm square, standing in for an infinite one.",
    ]
    return model


@builder("cylindrical_dra")
def _cyl_dra(design: DesignResult) -> Model:
    model = _base_model(design, "Cylindrical dielectric resonator antenna")
    a = _param(design, "a")
    height = _param(design, "h")
    eps_r = _param(design, "eps_r", default=10.0)
    ground = a * 8.0
    probe_r = a / 20.0
    model.solids += [
        Cylinder("ground", "PEC", "z", ground / 2.0, (0.0, 0.0)),
        Cylinder("resonator", f"eps_r={eps_r:g}", "z", a, (0.0, height)),
        Cylinder("probe", "PEC", "z", probe_r, (0.0, height * 0.6),
                 (a + probe_r * 2, 0.0)),
    ]
    model.ports.append(DiscretePort(
        "port1", (a + probe_r * 2, 0.0, -height * 0.05),
        (a + probe_r * 2, 0.0, 0.0)))
    model.notes += [
        "Probe placed just outside the puck to excite HE11. Its height and "
        "radial offset set the coupling and must be tuned; the spec models neither.",
        f"Ground plane rendered as a disc {ground * 1e3:.4g} mm across.",
    ]
    return model


@builder("hemispherical_dra")
def _hemi_dra(design: DesignResult) -> Model:
    """A sphere with everything below the ground plane cut away."""
    model = _base_model(design, "Hemispherical dielectric resonator antenna")
    a = _param(design, "a")
    eps_r = _param(design, "eps_r", default=10.0)
    ground = a * 8.0
    probe_r = a / 20.0
    model.solids += [
        Cylinder("ground", "PEC", "z", ground / 2.0, (0.0, 0.0)),
        Sphere("resonator_sphere", f"eps_r={eps_r:g}", a),
        Brick("lower_half", "VOID", (-a * 1.1, a * 1.1), (-a * 1.1, a * 1.1),
              (-a * 1.1, 0.0)),
        Cylinder("probe", "PEC", "z", probe_r, (0.0, a * 0.5), (a * 0.65, 0.0)),
    ]
    model.operations.append(Subtract("resonator_sphere", ("lower_half",)))
    model.ports.append(DiscretePort(
        "port1", (a * 0.65, 0.0, -a * 0.05), (a * 0.65, 0.0, 0.0)))
    model.notes += [
        "Built as a full sphere with the lower half subtracted, because neither "
        "backend has a hemisphere primitive.",
        "Probe inside the dielectric at 0.65a, a common starting point. Its "
        "depth and offset set the coupling and are not part of the spec, which "
        "solves the ISOLATED resonator.",
        f"Ground plane rendered as a disc {ground * 1e3:.4g} mm across. The "
        "exact solution assumes it is infinite - that is what makes the "
        "hemisphere equivalent to a whole sphere.",
    ]
    return model


@builder("conical_monopole")
def _conical_monopole(design: DesignResult) -> Model:
    model = _base_model(design, "Conical monopole over a ground plane")
    height = _param(design, "height")
    base_d = _param(design, "base_diameter")
    ground = max(base_d * 3.0, height * 3.0)
    gap = height / 50.0
    model.solids += [
        Cylinder("ground", "PEC", "z", ground / 2.0, (0.0, 0.0)),
        Cone("cone", "PEC", "z", 0.0, base_d / 2.0, (gap, gap + height)),
    ]
    model.ports.append(DiscretePort("port1", (0.0, 0.0, 0.0), (0.0, 0.0, gap)))
    model.notes += [
        f"Feed gap {gap * 1e3:.4g} mm between the apex and the plane. The gap "
        "sets the high-frequency limit and is not given by the spec.",
        f"Ground plane rendered as a disc {ground * 1e3:.4g} mm across. The "
        "spec's impedance assumes an infinite plane AND an infinite cone; a "
        "truncated cone's impedance oscillates about that value with frequency.",
        "Solid cone. A skeletal wire cone behaves similarly at the low end and "
        "worse at the high end.",
    ]
    return model


@builder("biconical")
def _biconical(design: DesignResult) -> Model:
    model = _base_model(design, "Biconical antenna")
    slant = _param(design, "Lc")
    theta_full = _param(design, "theta_h", default=0.5236)
    half_angle = theta_full / 2.0
    base_r = slant * math.sin(half_angle)
    height = slant * math.cos(half_angle)
    gap = height / 50.0
    model.solids += [
        Cone("cone_upper", "PEC", "z", 0.0, base_r, (gap / 2.0, gap / 2.0 + height)),
        Cone("cone_lower", "PEC", "z", base_r, 0.0, (-gap / 2.0 - height, -gap / 2.0)),
    ]
    model.ports.append(DiscretePort(
        "port1", (0.0, 0.0, -gap / 2.0), (0.0, 0.0, gap / 2.0)))
    model.notes += [
        f"Half angle {math.degrees(half_angle):.4g} deg, from the spec's full "
        f"cone angle of {math.degrees(theta_full):.4g} deg.",
        f"Feed gap {gap * 1e3:.4g} mm between the apexes - the spec assumes "
        "apexes meeting at a point, which cannot be meshed.",
        "The spec's characteristic impedance is the INFINITE-cone value. A "
        "truncated bicone's input impedance oscillates about it, converging as "
        "the cones grow electrically long.",
    ]
    return model


@builder("discone")
def _discone(design: DesignResult) -> Model:
    """Disc on top, cone below with its apex at the feed."""
    model = _base_model(design, "Discone")
    cone_h = _param(design, "cone_height")
    base_d = _param(design, "cone_base_diameter")
    disc_d = _param(design, "disc_diameter")
    gap = _param(design, "feed_gap")
    disc_t = disc_d / 200.0
    model.solids += [
        Cone("cone", "PEC", "z", base_d / 2.0, 0.0, (-cone_h, 0.0)),
        Cylinder("disc", "PEC", "z", disc_d / 2.0, (gap, gap + disc_t)),
    ]
    model.ports.append(DiscretePort("port1", (0.0, 0.0, 0.0), (0.0, 0.0, gap)))
    model.notes += [
        f"Cone apex at the origin opening downward, disc {gap * 1e3:.4g} mm "
        "above it. The coax inner connects to the disc and the shield to the "
        "cone, which this discrete port stands in for.",
        f"Disc rendered {disc_t * 1e3:.4g} mm thick; a real one is sheet metal "
        "and its thickness is not a design parameter.",
        "The feed gap is the main control on the high-frequency limit and is "
        "the least certain dimension in the spec - expect to tune it.",
    ]
    return model



@builder("small_circular_loop")
@builder("one_wavelength_circular_loop")
def _circular_loop(design: DesignResult) -> Model:
    """A wire ring in the x-y plane, fed across a gap cut out of it."""
    model = _base_model(design, "Circular loop")
    radius = _param(design, "a")
    wire = _param(design, "b", "aw", default=radius / 100.0)
    gap = max(radius / 50.0, wire * 3.0)
    model.solids += [
        Torus("loop", "PEC", "z", radius, wire),
        Brick("feed_gap_cut", "VOID", (radius - 2 * wire, radius + 2 * wire),
              (-gap / 2, gap / 2), (-2 * wire, 2 * wire)),
    ]
    model.operations.append(Subtract("loop", ("feed_gap_cut",)))
    model.ports.append(DiscretePort(
        "port1", (radius, -gap / 2, 0.0), (radius, gap / 2, 0.0)))
    model.notes += [
        f"Feed gap {gap * 1e3:.4g} mm, cut from the ring with a boolean. The "
        "spec assumes an unbroken loop, so the gap is a modelling necessity "
        "rather than part of the design - and for an electrically small loop it "
        "adds capacitance that shifts the tuning.",
        f"Wire radius {wire * 1e3:.4g} mm as a solid PEC tube.",
        "A small loop's radiation resistance is milliohms, so the simulated "
        "input impedance will be dominated by whatever loss the solver models. "
        "Set the conductor to a real metal rather than PEC if efficiency is the "
        "question.",
    ]
    return model


@builder("small_square_loop")
@builder("quad_loop_square")
@builder("alford_loop")
def _square_loop(design: DesignResult) -> Model:
    """Four wires round a square, fed across a gap in the bottom side."""
    model = _base_model(design, "Square loop")
    side = _param(design, "s")
    wire = _param(design, "b", "aw", default=side / 200.0)
    gap = max(side / 100.0, wire * 3.0)
    h = side / 2.0
    model.solids += [
        Cylinder("side_top", "PEC", "x", wire, (-h, h), (h, 0.0)),
        Cylinder("side_left", "PEC", "y", wire, (-h, h), (-h, 0.0)),
        Cylinder("side_right", "PEC", "y", wire, (-h, h), (h, 0.0)),
        Cylinder("side_bottom_a", "PEC", "x", wire, (-h, -gap / 2), (-h, 0.0)),
        Cylinder("side_bottom_b", "PEC", "x", wire, (gap / 2, h), (-h, 0.0)),
    ]
    model.ports.append(DiscretePort(
        "port1", (-gap / 2, -h, 0.0), (gap / 2, -h, 0.0)))
    model.notes += [
        f"Feed gap {gap * 1e3:.4g} mm in the centre of the bottom side; the spec "
        "does not define one.",
        "Corners are butt joints between straight tubes rather than mitred or "
        "radiused, which is a small perturbation at the current maximum.",
        "The Alford loop's corner capacitors, if this is one, are NOT modelled - "
        "without them the current is not uniform and the azimuthal pattern will "
        "not be the near-circle the spec predicts.",
    ]
    return model


@builder("halo_loop")
def _halo_loop(design: DesignResult) -> Model:
    """A half-wave dipole bent into a ring, fed across the gap between its tips."""
    model = _base_model(design, "Halo loop")
    diameter = _param(design, "Dm")
    gap = _param(design, "g")
    radius = diameter / 2.0
    wire = _param(design, "b", "aw", default=radius / 60.0)
    model.solids += [
        Torus("halo", "PEC", "z", radius, wire),
        Brick("tip_gap_cut", "VOID", (radius - 2 * wire, radius + 2 * wire),
              (-gap / 2, gap / 2), (-2 * wire, 2 * wire)),
    ]
    model.operations.append(Subtract("halo", ("tip_gap_cut",)))
    model.ports.append(DiscretePort(
        "port1", (radius, -gap / 2, 0.0), (radius, gap / 2, 0.0)))
    model.notes += [
        f"Tip gap {gap * 1e3:.4g} mm, from the spec - unlike the other loops "
        "here, the halo's gap IS a design parameter, because the capacitance "
        "across it sets the resonance.",
        "Modelled as a full ring with the gap cut out. A real halo is usually a "
        "ring with flattened or capacitor-hatted tips, which increases that "
        "capacitance and lowers the resonant frequency.",
        "This archetype is marked low confidence in the catalogue; treat the "
        "simulated resonance as the answer and the spec's as a starting point.",
    ]
    return model


@builder("half_wave_slot")
def _ground_plane_slot(design: DesignResult) -> Model:
    """A slot cut through a finite ground plane, fed across its middle."""
    model = _base_model(design, "Slot in a ground plane")
    length = _param(design, "L")
    width = _param(design, "w", default=length / 20.0)
    plate = 3.0 * length
    thick = plate / 2000.0
    model.solids += [
        Brick("ground_plane", "PEC", (-plate / 2, plate / 2),
              (-plate / 2, plate / 2), (-thick, 0.0)),
        Brick("slot_cut", "VOID", (-length / 2, length / 2),
              (-width / 2, width / 2), (-thick * 1.5, thick * 0.5)),
    ]
    model.operations.append(Subtract("ground_plane", ("slot_cut",)))
    model.ports.append(DiscretePort(
        "port1", (0.0, -width / 2, -thick / 2), (0.0, width / 2, -thick / 2)))
    model.notes += [
        f"Ground plane rendered {plate * 1e3:.4g} mm square and "
        f"{thick * 1e3:.4g} mm thick. Babinet's principle assumes an INFINITE, "
        "zero-thickness, perfectly conducting screen; a finite plate adds edge "
        "diffraction and a real thickness lowers the resonance slightly.",
        "The slot radiates from BOTH faces, so compare against the spec's "
        "two-sided impedance unless you back it with a cavity - see "
        "cavity_backed_slot, which halves the radiation and doubles the "
        "resistance.",
        "Port across the slot at its centre, where the slot voltage is highest.",
    ]
    return model


@builder("folded_slot")
def _folded_slot(design: DesignResult) -> Model:
    """N parallel slots sharing a ground plane, fed across one of them."""
    model = _base_model(design, "Folded slot")
    length = _param(design, "L")
    width = _param(design, "w", default=length / 20.0)
    n = int(round(_param(design, "N", default=2.0)))
    pitch = 3.0 * width
    plate = 3.0 * length
    thick = plate / 2000.0
    offsets = [(i - (n - 1) / 2.0) * pitch for i in range(n)]
    model.solids.append(
        Brick("ground_plane", "PEC", (-plate / 2, plate / 2),
              (-plate / 2, plate / 2), (-thick, 0.0)))
    for i, y in enumerate(offsets):
        model.solids.append(
            Brick(f"slot_cut_{i + 1}", "VOID", (-length / 2, length / 2),
                  (y - width / 2, y + width / 2), (-thick * 1.5, thick * 0.5)))
    model.operations.append(
        Subtract("ground_plane", tuple(f"slot_cut_{i + 1}" for i in range(n))))
    y0 = offsets[0]
    model.ports.append(DiscretePort(
        "port1", (0.0, y0 - width / 2, -thick / 2),
        (0.0, y0 + width / 2, -thick / 2)))
    model.notes += [
        f"{n} slots on a {pitch * 1e3:.4g} mm pitch, three slot widths apart. "
        "The N^2 impedance division assumes the slots are CLOSE compared with a "
        "wavelength; widening this pitch moves the result away from the spec.",
        "Only the first slot is driven. The others are shorted at their ends by "
        "the surrounding conductor, which is what makes the folding work.",
        "Infinite zero-thickness screen assumed, as for the plain slot.",
    ]
    return model


@builder("cavity_backed_slot")
def _cavity_backed_slot(design: DesignResult) -> Model:
    """A slot in a ground plane with a box behind it, so it radiates one way."""
    model = _base_model(design, "Cavity-backed slot")
    length = _param(design, "slot_length")
    width = _param(design, "slot_width", default=length / 20.0)
    depth = _param(design, "cavity_depth")
    plate = 3.0 * length
    thick = plate / 2000.0
    cav_l, cav_w = length * 1.2, max(width * 6.0, length * 0.3)
    model.solids += [
        Brick("ground_plane", "PEC", (-plate / 2, plate / 2),
              (-plate / 2, plate / 2), (-thick, 0.0)),
        Brick("slot_cut", "VOID", (-length / 2, length / 2),
              (-width / 2, width / 2), (-thick * 1.5, thick * 0.5)),
        Brick("cavity", "VACUUM", (-cav_l / 2, cav_l / 2),
              (-cav_w / 2, cav_w / 2), (-depth - thick, -thick)),
    ]
    model.operations.append(Subtract("ground_plane", ("slot_cut",)))
    model.ports.append(DiscretePort(
        "port1", (0.0, -width / 2, -thick / 2), (0.0, width / 2, -thick / 2)))
    model.notes += [
        f"Cavity {cav_l * 1e3:.4g} x {cav_w * 1e3:.4g} x {depth * 1e3:.4g} mm. "
        "Only the DEPTH comes from the spec - the lateral size is chosen here, "
        "and it matters: a cavity close to the slot in width loads it and "
        "shifts the resonance.",
        "The cavity is drawn as a vacuum volume; enclose it in PEC walls in the "
        "solver, or the slot radiates backwards and the point is lost.",
        "Backing the slot halves the radiated power and doubles the input "
        "resistance against the two-sided case.",
    ]
    return model


def _slotted_guide(design: DesignResult, title: str, count: int,
                   spacing: float, offset: float, alternate: bool) -> Model:
    """Shared body for the three waveguide slot archetypes."""
    model = _base_model(design, title)
    a = _param(design, "a_wg")
    b = _param(design, "b_wg")
    slot_l = _param(design, "slot_length")
    slot_w = slot_l / 16.0
    wall = min(a, b) / 20.0
    run = max((count - 1) * spacing + 4 * slot_l, 4 * slot_l)
    model.solids += [
        Brick("guide_interior", "VACUUM", (-a / 2, a / 2), (-b / 2, b / 2),
              (0.0, run)),
        Brick("broad_wall", "PEC", (-a / 2, a / 2), (b / 2, b / 2 + wall),
              (0.0, run)),
    ]
    tools = []
    for i in range(count):
        z = (i - (count - 1) / 2.0) * spacing + run / 2.0
        x = offset * (-1 if (alternate and i % 2) else 1)
        name = f"slot_cut_{i + 1}"
        tools.append(name)
        model.solids.append(
            Brick(name, "VOID", (x - slot_w / 2, x + slot_w / 2),
                  (b / 2 - wall, b / 2 + 2 * wall),
                  (z - slot_l / 2, z + slot_l / 2)))
    model.operations.append(Subtract("broad_wall", tuple(tools)))
    model.notes += [
        f"Broad wall drawn {wall * 1e3:.4g} mm thick with the slots cut through "
        "it. The other three walls are left to the surrounding PEC boundary, as "
        "for open_ended_waveguide.",
        f"Slot width set to a sixteenth of its length ({slot_w * 1e3:.4g} mm); "
        "the spec gives only the length. Width affects bandwidth more than "
        "resonance, but it is a choice made here and not by the design.",
        "Assign a waveguide port to the z = 0 face. A resonant array needs a "
        "short a quarter guide wavelength beyond the last slot; a "
        "travelling-wave array needs a matched load there instead.",
    ]
    return model


@builder("waveguide_longitudinal_slot")
def _single_guide_slot(design: DesignResult) -> Model:
    model = _slotted_guide(design, "Longitudinal slot in a waveguide broad wall",
                           1, 0.0, _param(design, "x1"), False)
    model.notes.append(
        "One slot at the offset the spec gives. Its conductance goes as "
        "sin^2(pi*x/a), so the offset is the whole design variable here.")
    return model


@builder("waveguide_slot_array_resonant")
def _resonant_slot_array(design: DesignResult) -> Model:
    n = int(round(_param(design, "N", default=12.0)))
    model = _slotted_guide(design, "Resonant longitudinal slot array", n,
                           _param(design, "spacing"), _param(design, "offset"),
                           True)
    model.notes.append(
        f"{n} slots at half a GUIDE wavelength, offsets ALTERNATING about the "
        "centreline. That alternation undoes the 180 degrees of propagation "
        "phase between neighbours; put them all on one side and the array "
        "splits into two beams off broadside instead of one on it.")
    return model


@builder("waveguide_slot_array_travelling_wave")
def _travelling_slot_array(design: DesignResult) -> Model:
    n = int(round(_param(design, "N", default=20.0)))
    a = _param(design, "a_wg")
    model = _slotted_guide(design, "Travelling-wave longitudinal slot array", n,
                           _param(design, "spacing"), a * 0.13, True)
    model.notes.append(
        "Slot offsets are UNIFORM here at 13% of the broad wall. A real "
        "travelling-wave array tapers them along the guide to hold the aperture "
        "distribution as power drains away; the spec gives no taper, so none is "
        "applied.")
    return model


@builder("long_wire_travelling")
def _long_wire(design: DesignResult) -> Model:
    """One straight wire, fed at one end and terminated at the other."""
    model = _base_model(design, "Terminated travelling-wave long wire")
    length = _param(design, "L")
    wire = _param(design, "aw", default=length / 4000.0)
    gap = max(length / 2000.0, wire * 3.0)
    model.solids.append(Cylinder("wire", "PEC", "x", wire, (gap, length + gap)))
    model.ports += [
        DiscretePort("port1", (0.0, 0.0, 0.0), (gap, 0.0, 0.0)),
        DiscretePort("port2", (length + gap, 0.0, 0.0),
                     (length + 2 * gap, 0.0, 0.0), impedance=600.0),
    ]
    model.notes += [
        "TWO ports. Port 1 is the feed; port 2 stands in for the TERMINATION and "
        "must be a real absorbing load - the travelling wave is the whole "
        "premise. Leave it open and the reflected wave restores a standing-wave "
        "pattern and the beam splits.",
        "600 ohm is a placeholder. The right value is the wire's characteristic "
        "impedance against ground, which depends on a height this spec does not "
        "model.",
        "Free space, no ground plane. A real long wire runs above earth, which "
        "adds an image and tilts the beam upward.",
    ]
    return model


@builder("leaky_wave_line_source")
def _leaky_wave(design: DesignResult) -> Model:
    """A waveguide with a continuous longitudinal slit in its broad wall."""
    model = _base_model(design, "Uniform leaky-wave line source")
    a = _param(design, "a")
    length = _param(design, "L")
    b = a * 0.45
    wall = a / 20.0
    slit = a / 12.0
    model.solids += [
        Brick("guide_interior", "VACUUM", (-a / 2, a / 2), (-b / 2, b / 2),
              (0.0, length)),
        Brick("broad_wall", "PEC", (-a / 2, a / 2), (b / 2, b / 2 + wall),
              (0.0, length)),
        Brick("slit_cut", "VOID", (-slit / 2, slit / 2),
              (b / 2 - wall, b / 2 + 2 * wall), (length * 0.05, length * 0.95)),
    ]
    model.operations.append(Subtract("broad_wall", ("slit_cut",)))
    model.notes += [
        f"Narrow wall set to {b * 1e3:.4g} mm, 0.45 of the broad wall, because "
        "the spec gives only a. Feed a waveguide port at z = 0 and a matched "
        "load at the far end.",
        f"Slit width {slit * 1e3:.4g} mm. This is the leakage control and the "
        "spec parameterises it only as alpha_norm, never as a dimension, so "
        "expect to tune it - the slit sets how much power reaches the load "
        "rather than radiating.",
        "The slit stops short of both ends so the port faces stay solid.",
    ]
    return model


@builder("planar_monopole_rectangular")
def _planar_monopole_rect(design: DesignResult) -> Model:
    """A flat plate standing on edge above a ground plane."""
    model = _base_model(design, "Rectangular planar monopole")
    height = _param(design, "Lp")
    width = _param(design, "Wp")
    gap = _param(design, "p_gap")
    thick = max(width, height) / 500.0
    ground = 4.0 * max(width, height)
    model.solids += [
        Brick("ground_plane", "PEC", (-ground / 2, ground / 2),
              (-ground / 2, ground / 2), (-thick, 0.0)),
        Brick("plate", "PEC", (-width / 2, width / 2), (-thick / 2, thick / 2),
              (gap, gap + height)),
    ]
    model.ports.append(DiscretePort("port1", (0.0, 0.0, 0.0), (0.0, 0.0, gap)))
    model.notes += [
        f"Feed gap {gap * 1e3:.4g} mm, from the spec - a real design parameter "
        "here, since the lower cut-off goes as 1/(L + W/2 + p).",
        f"Plate rendered {thick * 1e3:.4g} mm thick; a real one is sheet metal or "
        "copper on a substrate, and its thickness is not a design variable.",
        f"Ground plane {ground * 1e3:.4g} mm square. A small ground plane raises "
        "the cut-off and spoils the pattern; the spec assumes a large one.",
    ]
    return model


@builder("planar_monopole_circular")
def _planar_monopole_disc(design: DesignResult) -> Model:
    """A disc standing on edge above a ground plane."""
    model = _base_model(design, "Circular disc planar monopole")
    radius = _param(design, "r_disc")
    gap = _param(design, "p_gap")
    thick = radius / 250.0
    ground = 8.0 * radius
    model.solids += [
        Brick("ground_plane", "PEC", (-ground / 2, ground / 2),
              (-ground / 2, ground / 2), (-thick, 0.0)),
        Cylinder("disc", "PEC", "y", radius, (-thick / 2, thick / 2),
                 (0.0, gap + radius)),
    ]
    model.ports.append(DiscretePort("port1", (0.0, 0.0, 0.0), (0.0, 0.0, gap)))
    model.notes += [
        f"Disc of radius {radius * 1e3:.4g} mm standing VERTICALLY, its lowest "
        f"point {gap * 1e3:.4g} mm above the plane. The gap is a design "
        "parameter: the cut-off goes as 1/(2r + r + p).",
        f"Rendered {thick * 1e3:.4g} mm thick.",
        "Marked low confidence in the catalogue - the equivalent-cylinder "
        "cut-off rule is empirical and good to perhaps 10%.",
    ]
    return model


@builder("annular_ring_patch")
def _annular_ring(design: DesignResult) -> Model:
    """A ring of copper on a substrate: outer disc with the inner one removed."""
    model = _base_model(design, "Annular ring patch")
    a_in = _param(design, "a_in")
    b_out = _param(design, "b_out")
    h = _param(design, "h")
    eps_r = _param(design, "eps_r", default=2.2)
    sub = b_out * 3.0
    model.solids += [
        Brick("substrate", f"eps_r={eps_r:g}", (-sub, sub), (-sub, sub), (0.0, h)),
        Brick("ground", "PEC", (-sub, sub), (-sub, sub), (0.0, 0.0)),
        Cylinder("ring", "PEC", "z", b_out, (h, h)),
        Cylinder("ring_hole", "VOID", "z", a_in, (h, h)),
    ]
    model.operations.append(Subtract("ring", ("ring_hole",)))
    feed_r = (a_in + b_out) / 2.0
    model.ports.append(DiscretePort("port1", (feed_r, 0.0, 0.0), (feed_r, 0.0, h)))
    model.notes += [
        f"Probe at the MEAN radius {feed_r * 1e3:.4g} mm as a starting point. The "
        "spec gives no feed position, and a ring is awkward to feed - probe, "
        "coupled microstrip and aperture each perturb the TM11 mode differently.",
        "Drawn as a zero-thickness annulus. Real etched copper has thickness, "
        "which lowers the resonance by a fraction of a percent.",
        "The spec's resonance is fitted to the exact Bessel cross-product and "
        "carries no fringing correction on the radii, so expect the simulated "
        "resonance a few percent low.",
    ]
    return model


@builder("pifa")
def _pifa(design: DesignResult) -> Model:
    """Top plate over a ground plane, shorted along one edge, fed beside it."""
    model = _base_model(design, "Planar inverted-F antenna")
    length = _param(design, "L")
    width = _param(design, "W")
    height = _param(design, "h")
    short_w = _param(design, "Ws", default=width)
    ground = max(4.0 * length, 4.0 * width)
    thick = height / 20.0
    model.solids += [
        Brick("ground_plane", "PEC", (-ground / 2, ground / 2),
              (-ground / 2, ground / 2), (-thick, 0.0)),
        Brick("top_plate", "PEC", (-width / 2, width / 2),
              (-length / 2, length / 2), (height, height)),
        Brick("shorting_wall", "PEC", (-short_w / 2, short_w / 2),
              (-length / 2, -length / 2), (0.0, height)),
    ]
    feed_y = -length / 2 + length * 0.15
    model.ports.append(DiscretePort("port1", (0.0, feed_y, 0.0), (0.0, feed_y, height)))
    model.notes += [
        f"Shorting wall {short_w * 1e3:.4g} mm wide against a {width * 1e3:.4g} mm "
        "plate. That ratio is the PIFA's main tuning control: a full-width short "
        "resonates near L + h = lambda/4, a narrow one nearer L + W + h = lambda/4.",
        f"Probe {length * 0.15 * 1e3:.4g} mm from the shorted edge. The spec gives "
        "no feed position and the impedance is very sensitive to it.",
        "Air between plate and ground; a dielectric there lowers the resonance.",
        "Marked low confidence. A PIFA is dominated by the ground plane it sits "
        "on, rendered here as a plain rectangle rather than the handset it would "
        "really be.",
    ]
    return model


@builder("stacked_patch")
def _stacked_patch(design: DesignResult) -> Model:
    """Driven patch on its substrate, parasitic patch suspended above it."""
    model = _base_model(design, "Stacked patch")
    w = _param(design, "W")
    length = _param(design, "L")
    w2 = _param(design, "W2")
    l2 = _param(design, "L2")
    h = _param(design, "h")
    h2 = _param(design, "h2")
    eps_r = _param(design, "eps_r", default=2.2)
    margin = max(w, length) * 0.6
    sub_w, sub_l = max(w, w2) + 2 * margin, max(length, l2) + 2 * margin
    model.solids += [
        Brick("substrate", f"eps_r={eps_r:g}", (-sub_w / 2, sub_w / 2),
              (-sub_l / 2, sub_l / 2), (0.0, h)),
        Brick("ground", "PEC", (-sub_w / 2, sub_w / 2), (-sub_l / 2, sub_l / 2),
              (0.0, 0.0)),
        Brick("driven_patch", "PEC", (-w / 2, w / 2), (-length / 2, length / 2),
              (h, h)),
        Brick("parasitic_patch", "PEC", (-w2 / 2, w2 / 2), (-l2 / 2, l2 / 2),
              (h + h2, h + h2)),
    ]
    model.ports.append(DiscretePort(
        "port1", (0.0, -length / 4, 0.0), (0.0, -length / 4, h)))
    model.notes += [
        f"Parasitic patch suspended {h2 * 1e3:.4g} mm above the driven one, in "
        "air. The spacer is not modelled: foam or honeycomb is close to air, a "
        "dielectric one is not.",
        "Probe at a quarter of the driven patch's length, which is a guess. The "
        "spec gives no feed position, and on a stack the probe inductance is "
        "often large enough to need a series capacitor.",
        "Marked LOW CONFIDENCE: the bandwidth multiplier is an expectation drawn "
        "from published designs, not a computed result. The geometry here is "
        "sound; the predicted bandwidth is the part to distrust.",
    ]
    return model


@builder("fresnel_zone_plate")
def _zone_plate(design: DesignResult) -> Model:
    """Concentric metal rings on a flat plane: a Soret amplitude zone plate.

    Zone m runs from r(m-1) to r(m) with r(m) = sqrt(m*lambda*F + (m*lambda/2)^2).
    Odd zones are left open and even zones filled with metal, so the blocked
    contributions are the ones that would have arrived out of phase.
    """
    model = _base_model(design, "Fresnel zone plate")
    f0 = float(design.requirements.get("f0", 1e9))
    focal = _param(design, "F")
    zones = int(round(_param(design, "M", default=4.0)))
    lam = 2.99792458e8 / f0
    radii = [math.sqrt(m * lam * focal + (m * lam / 2.0) ** 2)
             for m in range(0, zones + 1)]
    thick = radii[-1] / 200.0
    rings = 0
    for m in range(2, zones + 1, 2):          # even zones carry the metal
        rings += 1
        outer, inner = radii[m], radii[m - 1]
        model.solids += [
            Cylinder(f"ring_{rings}", "PEC", "z", outer, (0.0, thick)),
            Cylinder(f"ring_{rings}_hole", "VOID", "z", inner,
                     (-thick, 2 * thick)),
        ]
        model.operations.append(
            Subtract(f"ring_{rings}", (f"ring_{rings}_hole",)))
    model.notes += [
        f"{rings} metal rings from {zones} zones, outermost radius "
        f"{radii[-1] * 1e3:.4g} mm. Zone edges from "
        "r(m) = sqrt(m*lambda*F + (m*lambda/2)^2) - the exact form, not the "
        "sqrt(m*lambda*F) approximation, which matters when F is only a few "
        "wavelengths.",
        f"Rendered {thick * 1e3:.4g} mm thick. A real plate is etched foil on a "
        "thin dielectric, which is not modelled - the dielectric shifts the "
        "focus slightly.",
        "This is the OPAQUE (Soret) plate, central zone open. Its efficiency "
        "under a real feed is the spec's aperture_efficiency (0.107 for the "
        "default design), not the uniform-illumination 1/pi^2. The phase-reversing "
        "and stepped versions replace these rings with dielectric steps and are "
        "NOT built here.",
        "No feed is included. Illuminate it from the focal point with a separate "
        "source; there is no port to drive.",
    ]
    model.built_geometry = rings > 0
    return model


@builder("luneburg_lens")
def _luneburg(design: DesignResult) -> Model:
    """A Luneburg lens as the concentric shells one is actually built from.

    The ideal gradient n(r) = sqrt(2 - (r/R)^2) is continuous and unmakeable, so
    real lenses are stepped. Each shell here gets the index at its own mid-radius.
    """
    model = _base_model(design, "Luneburg lens")
    radius = _param(design, "R")
    shells = 8
    edges = [radius * (i / shells) ** 0.5 for i in range(shells + 1)]
    for i in range(shells, 0, -1):
        outer, inner = edges[i], edges[i - 1]
        mid = 0.5 * (outer + inner) / radius
        eps = 2.0 - mid ** 2
        model.solids.append(
            Sphere(f"shell_{i}", f"eps_r={eps:.4g}", outer))
        if inner > 0:
            model.solids.append(Sphere(f"shell_{i}_core", "VOID", inner))
            model.operations.append(Subtract(f"shell_{i}", (f"shell_{i}_core",)))
    model.notes += [
        f"{shells} shells, boundaries placed at equal steps in r^2 rather than "
        "in r, because the permittivity is exactly linear in r^2 - that spacing "
        "makes each shell span the same permittivity range.",
        f"Permittivity runs from {2.0 - (0.5 * edges[1] / radius) ** 2:.4g} at "
        f"the core to {2.0 - (0.5 * (edges[-1] + edges[-2]) / radius) ** 2:.4g} "
        "at the rim, against the ideal 2.0 and 1.0. The stepping is the whole "
        "difference between this and the ideal lens: the steps scatter, costing "
        "a few tenths of a dB and raising the sidelobes.",
        "No feed is included. A Luneburg lens is fed from a point ON its surface, "
        "and any number of feeds can share it - which is the reason to build one.",
        "Materials are named by permittivity; define them in the solver before "
        "running the macro.",
    ]
    return model


@builder("top_loaded_monopole")
def _top_loaded(design: DesignResult) -> Model:
    """Short vertical with a disc hat, over a ground plane."""
    model = _base_model(design, "Capacitively top-loaded monopole")
    height = _param(design, "h")
    hat_r = _param(design, "a_hat")
    wire = _param(design, "aw", default=height / 500.0)
    gap = max(height / 200.0, wire * 3.0)
    ground = max(6.0 * hat_r, 3.0 * height)
    hat_t = max(hat_r / 200.0, wire)
    model.solids += [
        Cylinder("ground_plane", "PEC", "z", ground / 2.0, (-ground / 400.0, 0.0)),
        Cylinder("rod", "PEC", "z", wire, (gap, gap + height)),
        Cylinder("top_hat", "PEC", "z", hat_r, (gap + height, gap + height + hat_t)),
    ]
    model.ports.append(DiscretePort("port1", (0.0, 0.0, 0.0), (0.0, 0.0, gap)))
    model.notes += [
        f"Hat radius {hat_r * 1e3:.4g} mm comes from the spec's "
        "hat_radius_over_lambda, which is an INPUT there and not derived from "
        "the loading it achieves. The spec's beta_top - how nearly uniform the "
        "hat makes the rod current - is likewise an input, so the two are not "
        "tied together here. Simulate to find the beta this hat really gives.",
        f"Feed gap {gap * 1e3:.4g} mm between the plane and the rod base.",
        f"Ground plane rendered as a disc {ground * 1e3:.4g} mm across. The "
        "spec's directivity of exactly 3.0 assumes it is infinite; a real "
        "radial ground system also adds loss that often exceeds the radiation "
        "resistance entirely.",
    ]
    return model


@builder("inductively_loaded_monopole")
def _loaded_monopole(design: DesignResult) -> Model:
    """Short whip with a gap where the loading coil goes."""
    model = _base_model(design, "Inductively loaded monopole")
    height = _param(design, "h")
    wire = _param(design, "aw", default=height / 500.0)
    gap = max(height / 200.0, wire * 3.0)
    coil_gap = height / 20.0
    ground = 3.0 * height
    model.solids += [
        Cylinder("ground_plane", "PEC", "z", ground / 2.0, (-ground / 400.0, 0.0)),
        Cylinder("rod_lower", "PEC", "z", wire, (gap, gap + coil_gap)),
        Cylinder("rod_upper", "PEC", "z", wire,
                 (gap + 2 * coil_gap, gap + height + coil_gap)),
    ]
    model.ports += [
        DiscretePort("port1", (0.0, 0.0, 0.0), (0.0, 0.0, gap)),
        DiscretePort("port2", (0.0, 0.0, gap + coil_gap),
                     (0.0, 0.0, gap + 2 * coil_gap)),
    ]
    model.notes += [
        "TWO ports. Port 1 is the feed; port 2 is the gap where the LOADING "
        "COIL goes, and it must be given the spec's loading_inductance_H with "
        "the spec's coil_q as a series resistance. Leave it open and the whip "
        "is just a short whip.",
        f"Coil gap {coil_gap * 1e3:.4g} mm at a twentieth of the height, which "
        "puts it near the base. The spec captures coil POSITION only through "
        "beta_top, so moving this gap up the rod will not match the spec's "
        "numbers unless beta_top is changed to suit.",
        "A real coil is a wound solenoid with its own self-capacitance and "
        "radiation; a lumped port is the idealisation the spec assumes.",
        f"Ground plane rendered as a disc {ground * 1e3:.4g} mm across; the "
        "spec's ground loss is an input, not modelled geometrically.",
    ]
    return model


@builder("multiturn_small_loop")
def _multiturn_loop(design: DesignResult) -> Model:
    """A stack of coaxial rings standing in for a closely wound coil."""
    model = _base_model(design, "Multi-turn small loop")
    radius = _param(design, "a")
    turns = int(round(_param(design, "N", default=10.0)))
    wire = _param(design, "b", "aw", default=radius / 200.0)
    pitch = 3.0 * wire
    gap = max(radius / 50.0, wire * 3.0)
    for i in range(turns):
        z = (i - (turns - 1) / 2.0) * pitch
        model.solids.append(
            Torus(f"turn_{i + 1}", "PEC", "z", radius, wire, (0.0, 0.0, z)))
    z0 = -(turns - 1) / 2.0 * pitch
    model.solids.append(
        Brick("feed_gap_cut", "VOID", (radius - 2 * wire, radius + 2 * wire),
              (-gap / 2, gap / 2), (z0 - 2 * wire, z0 + 2 * wire)))
    model.operations.append(Subtract("turn_1", ("feed_gap_cut",)))
    model.ports.append(DiscretePort(
        "port1", (radius, -gap / 2, z0), (radius, gap / 2, z0)))
    model.notes += [
        f"{turns} turns as SEPARATE coaxial rings on a {pitch * 1e3:.4g} mm "
        "pitch, not as a helix. For a closely wound coil the difference is "
        "small, and the rings are far easier to mesh - but they are also not "
        "electrically connected to each other here.",
        "CONNECT THE TURNS before solving, or this models one driven ring and "
        f"{turns - 1} parasitic ones rather than an N-turn loop. The spec's "
        "N-squared radiation resistance depends entirely on the turns carrying "
        "the same current in series.",
        f"Feed gap {gap * 1e3:.4g} mm cut from the bottom turn.",
        "The spec's coil length l_coil is not used here; the pitch is set from "
        "the wire radius instead, so a loosely wound coil will not match.",
    ]
    return model


def build(design: DesignResult) -> Model:
    """Build a model for a design, or a parameters-only model if we cannot.

    A missing builder is reported, never faked.
    """
    fn = BUILDERS.get(design.archetype)
    if fn is None:
        model = _base_model(design, design.archetype.replace("_", " ").title())
        model.built_geometry = False
        model.notes.append(
            f"NO SOLID GEOMETRY GENERATED for '{design.archetype}'. This toolkit "
            "only builds models for archetypes whose construction is unambiguous "
            "from their primary dimensions. The synthesised parameters are "
            "exported below as named variables so you can drive your own "
            "geometry from them.")
        return model
    try:
        return fn(design)
    except KeyError as exc:
        model = _base_model(design, design.archetype)
        model.built_geometry = False
        model.notes.append(
            f"Geometry could not be built: {exc}. Supply the missing requirement "
            "and re-synthesise. Parameters exported below.")
        return model
