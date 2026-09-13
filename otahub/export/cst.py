"""CST Studio Suite VBA macro generation.

Emits a .bas macro that builds the model when run from CST's macro editor.
Dimensions become named CST parameters so the geometry stays drivable from
the parameter list rather than being frozen numbers.
"""
from __future__ import annotations

import math

from .base import Brick, Cylinder, DiscretePort, Model

_MM = 1e3          # the macro works in millimetres


def _p(name: str) -> str:
    return f"{name}_mm"


def _expr(value: float) -> str:
    return f"{value * _MM:.6f}"


def render(model: Model) -> str:
    """Render a complete VBA macro for `model`."""
    out: list[str] = []
    add = out.append

    add("' " + "=" * 72)
    add(f"' OTA Hub Antenna Toolkit - CST Studio Suite macro")
    add(f"' archetype : {model.archetype}")
    add(f"' title     : {model.title}")
    add(f"' frequency : {model.frequency_hz / 1e9:.6g} GHz")
    add("' " + "=" * 72)
    for note in model.notes:
        for line in _wrap(note):
            add(f"' NOTE: {line}")
    if not model.built_geometry:
        add("'")
        add("' This macro defines PARAMETERS ONLY - it builds no solids.")
    add("' " + "=" * 72)
    add("")
    add("Option Explicit")
    add("")
    add("Sub Main ()")
    add("")
    add("    ' ---- units -------------------------------------------------------")
    add("    With Units")
    add('        .Geometry "mm"')
    add('        .Frequency "GHz"')
    add('        .Time "ns"')
    add("    End With")
    add("")
    add("    ' ---- parameters --------------------------------------------------")
    add("    ' Lengths are in millimetres; everything else is dimensionless.")
    for name, value in sorted(model.parameters.items()):
        if name in ("lambda0",):
            continue
        kind, rendered = classify(name, value, model.units)
        label = _p(name) if kind == "length" else (
            f"{name}_GHz" if kind == "frequency" else name)
        unit_note = {"length": "mm", "frequency": "GHz", "angle": "deg"}.get(kind, "")
        suffix = f"    ' [{unit_note}]" if unit_note else ""
        add(f'    StoreParameter "{label}", {rendered}{suffix}')
    add("")
    freq = model.frequency_hz / 1e9
    add("    ' ---- frequency range ---------------------------------------------")
    add(f"    Solver.FrequencyRange \"{freq * 0.7:.6g}\", \"{freq * 1.3:.6g}\"")
    add("")

    if not model.solids:
        add("    ' No geometry is generated for this archetype - see the notes above.")
        add("")
        add("End Sub")
        return "\n".join(out) + "\n"

    add("    ' ---- materials ---------------------------------------------------")
    for material in sorted({s.material for s in model.solids}):
        if material.startswith("eps_r="):
            eps = material.split("=", 1)[1]
            add("    With Material")
            add("        .Reset")
            add(f'        .Name "substrate_eps{eps.replace(".", "p")}"')
            add('        .Type "Normal"')
            add(f'        .Epsilon "{eps}"')
            add('        .Mu "1.0"')
            add('        .TanD "0.02"')
            add('        .TanDModel "ConstTanD"')
            add('        .Colour "0.8", "0.8", "0.4"')
            add('        .Transparency "50"')
            add("        .Create")
            add("    End With")
    add("")

    add("    ' ---- geometry ----------------------------------------------------")
    for solid in model.solids:
        add(_render_solid(solid))
    add("")

    if model.ports:
        add("    ' ---- ports -------------------------------------------------------")
        for port in model.ports:
            add(_render_port(port))
        add("")

    add("    ' ---- boundaries --------------------------------------------------")
    add("    With Boundary")
    add('        .Xmin "expanded open" : .Xmax "expanded open"')
    add('        .Ymin "expanded open" : .Ymax "expanded open"')
    add('        .Zmin "expanded open" : .Zmax "expanded open"')
    add("    End With")
    add("")
    add("    ' ---- far field monitor -------------------------------------------")
    add("    With Monitor")
    add("        .Reset")
    add(f'        .Name "farfield_{freq:.6g}GHz"')
    add('        .Domain "Frequency"')
    add('        .FieldType "Farfield"')
    add(f'        .Frequency "{freq:.6g}"')
    add("        .Create")
    add("    End With")
    add("")
    add("End Sub")
    return "\n".join(out) + "\n"


def _render_solid(solid) -> str:
    material = _material_name(solid.material)
    if isinstance(solid, Brick):
        return "\n".join([
            "    With Brick",
            "        .Reset",
            f'        .Name "{solid.name}"',
            '        .Component "component1"',
            f'        .Material "{material}"',
            f'        .Xrange "{_expr(solid.x[0])}", "{_expr(solid.x[1])}"',
            f'        .Yrange "{_expr(solid.y[0])}", "{_expr(solid.y[1])}"',
            f'        .Zrange "{_expr(solid.z[0])}", "{_expr(solid.z[1])}"',
            "        .Create",
            "    End With",
        ])
    if isinstance(solid, Cylinder):
        axis = solid.axis.lower()
        other = {"x": ("y", "z"), "y": ("x", "z"), "z": ("x", "y")}[axis]
        return "\n".join([
            "    With Cylinder",
            "        .Reset",
            f'        .Name "{solid.name}"',
            '        .Component "component1"',
            f'        .Material "{material}"',
            f'        .Axis "{axis}"',
            f'        .Outerradius "{_expr(solid.radius)}"',
            '        .Innerradius "0"',
            f'        .{axis.upper()}range "{_expr(solid.span[0])}", "{_expr(solid.span[1])}"',
            f'        .{other[0].upper()}center "{_expr(solid.centre[0])}"',
            f'        .{other[1].upper()}center "{_expr(solid.centre[1])}"',
            "        .Segments \"0\"",
            "        .Create",
            "    End With",
        ])
    return f"    ' unsupported solid type {type(solid).__name__}"


def _render_port(port: DiscretePort) -> str:
    return "\n".join([
        "    With DiscretePort",
        "        .Reset",
        f'        .PortNumber "1"',
        '        .Type "SParameter"',
        f'        .Impedance "{port.impedance:.6g}"',
        f'        .SetP1 "False", "{_expr(port.start[0])}", "{_expr(port.start[1])}", '
        f'"{_expr(port.start[2])}"',
        f'        .SetP2 "False", "{_expr(port.end[0])}", "{_expr(port.end[1])}", '
        f'"{_expr(port.end[2])}"',
        '        .LocalCoordinates "False"',
        "        .Create",
        "    End With",
    ])


def _material_name(material: str) -> str:
    if material == "PEC":
        return "PEC"
    if material == "VACUUM":
        return "Vacuum"
    if material == "VOID":
        return "Vacuum"
    if material.startswith("eps_r="):
        return "substrate_eps" + material.split("=", 1)[1].replace(".", "p")
    return material


def classify(name: str, value: float, units: dict[str, str]) -> tuple[str, str]:
    """Return (kind, rendered_value) for a parameter, driven by its spec unit.

    Using the declared unit rather than guessing from magnitude: an earlier
    heuristic treated anything under 1000 as a length and wrote a 319 ohm
    input resistance into the macro as "319105 mm".
    """
    unit = (units.get(name) or "").strip()
    if unit == "m":
        return "length", f"{value * _MM:.6f}"
    if unit == "Hz":
        return "frequency", f"{value / 1e9:.9g}"
    if unit in ("deg",):
        return "angle", f"{value:.9g}"
    if unit in ("rad",):
        return "angle", f"{math.degrees(value):.9g}"
    return "scalar", f"{value:.10g}"


def _wrap(text: str, width: int = 74) -> list[str]:
    import textwrap
    return textwrap.wrap(text, width) or [""]
