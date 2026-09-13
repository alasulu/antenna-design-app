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
