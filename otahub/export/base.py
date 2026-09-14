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
        "The probe beside the resonator is ONE of several ways to excite TE111 "
        "- slot and microstrip feeds behave differently and are not modelled. "
        "Probe height and offset both need tuning in the solver.",
        "The spec's magnetic-wall resonance runs 10-20% HIGH, so expect the "
        "simulated resonance BELOW the design frequency and scale the "
        "resonator up accordingly.",
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
