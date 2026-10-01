"""Ansys HFSS script generation.

Emits an IronPython script for the HFSS scripting interface. Run it from
Tools > Run Script inside Electronics Desktop, or with ansysedt -RunScript.
Dimensions become design variables so the model stays parametric.
"""
from __future__ import annotations

from .base import (Brick, Cone, Cylinder, DiscretePort, Model, Sphere,
                   Subtract, Torus, Unite, dielectric_name, parse_dielectric)
from .cst import classify

_MM = 1e3


def _var(name: str) -> str:
    return name if name.isidentifier() else name.replace("-", "_")


def render(model: Model) -> str:
    """Render a complete HFSS IronPython script for `model`."""
    out: list[str] = []
    add = out.append

    add("# " + "=" * 72)
    add("# OTA Hub Antenna Toolkit - Ansys HFSS script")
    add(f"# archetype : {model.archetype}")
    add(f"# title     : {model.title}")
    add(f"# frequency : {model.frequency_hz / 1e9:.6g} GHz")
    add("#")
    add("# Run from Electronics Desktop: Tools > Run Script")
    add("# " + "=" * 72)
    for note in model.notes:
        for line in _wrap(note):
            add(f"# NOTE: {line}")
    if not model.built_geometry:
        add("#")
        add("# This script defines VARIABLES ONLY - it builds no geometry.")
    add("# " + "=" * 72)
    add("")
    add("import ScriptEnv")
    add('ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")')
    add("oDesktop.RestoreWindow()")
    add("oProject = oDesktop.NewProject()")
    add('oProject.InsertDesign("HFSS", "OTAHub_%s", "DrivenModal", "")'
        % model.archetype[:24])
    add('oDesign = oProject.SetActiveDesign("OTAHub_%s")' % model.archetype[:24])
    add('oEditor = oDesign.SetActiveEditor("3D Modeler")')
    add("")
    add("# ---- design variables ---------------------------------------------")
    add("# Lengths in millimetres; everything else dimensionless.")
    add("variables = [")
    for name, value in sorted(model.parameters.items()):
        if name in ("lambda0",):
            continue
        kind, rendered = classify(name, value, model.units)
        suffix = {"length": "mm", "frequency": "GHz", "angle": "deg"}.get(kind, "")
        label = _var(name) + ("_GHz" if kind == "frequency" else "")
        add(f'    ("{label}", "{rendered}{suffix}"),')
    add("]")
    add("for vname, vvalue in variables:")
    add("    oDesign.ChangeProperty([")
    add('        "NAME:AllTabs",')
    add('        ["NAME:LocalVariableTab", ["NAME:PropServers", "LocalVariables"],')
    add('         ["NAME:NewProps", ["NAME:" + vname, "PropType:=", "VariableProp",')
    add('          "UserDef:=", True, "Value:=", vvalue]]]])')
    add("")

    if not model.solids:
        add("# No geometry is generated for this archetype - see the notes above.")
        add('oDesktop.AddMessage("", "", 0, "OTA Hub: variables only, no geometry.")')
        return "\n".join(out) + "\n"

    add("# ---- materials ----------------------------------------------------")
    materials = sorted({s.material for s in model.solids
                        if s.material.startswith("eps_r=")})
    if materials:
        add("oDefinitionManager = oProject.GetDefinitionManager()")
        for material in materials:
            eps, tand = parse_dielectric(material)      # the design's loss, not an imposed 0.02
            add("oDefinitionManager.AddMaterial([")
            add(f'    "NAME:{dielectric_name(material)}", "CoordinateSystemType:=", "Cartesian",')
            add(f'    "permittivity:=", "{eps}",')
            add(f'    "dielectric_loss_tangent:=", "{tand}"])')
        add("")

    add("# ---- geometry ------------------------------------------------------")
    for solid in model.solids:
        add(_render_solid(solid))
        add("")

    if model.operations:
        add("# ---- boolean operations ---------------------------------------------")
        for op in model.operations:
            add(_render_operation(op))
        add("")

    sheets = _pec_sheets(model)
    if sheets:
        add("# ---- PEC sheets ----------------------------------------------------")
        add("# A sheet has no volume, so a material on it means nothing: a zero-thickness")
        add("# conductor (patch, ground, shorting wall) needs a Perfect E boundary.")
        add('oModule = oDesign.GetModule("BoundarySetup")')
        add(f'oModule.AssignPerfectE(["NAME:PerfE_sheets", "Objects:=", {sheets!r},')
        add('    "InfGroundPlane:=", False])')
        add("")

    if model.ports:
        add("# ---- excitation ----------------------------------------------------")
        for index, port in enumerate(model.ports, start=1):
            add(_render_port(port, index))
        add("")

    freq = model.frequency_hz / 1e9
    lo, hi = (f / 1e9 for f in model.band_hz)
    # a wideband design meshes at the top of its band, the usual HFSS practice
    mesh = hi if hi > 1.35 * freq else freq
    add("# ---- radiation boundary --------------------------------------------")
    add("# An airbox a quarter wavelength clear of the structure on every side.")
    add('oEditor.CreateRegion([')
    add('    "NAME:RegionParameters",')
    add('    "+XPaddingType:=", "Percentage Offset", "+XPadding:=", "25",')
    add('    "-XPaddingType:=", "Percentage Offset", "-XPadding:=", "25",')
    add('    "+YPaddingType:=", "Percentage Offset", "+YPadding:=", "25",')
    add('    "-YPaddingType:=", "Percentage Offset", "-YPadding:=", "25",')
    add('    "+ZPaddingType:=", "Percentage Offset", "+ZPadding:=", "25",')
    add('    "-ZPaddingType:=", "Percentage Offset", "-ZPadding:=", "25"],')
    add('    ["NAME:Attributes", "Name:=", "Region", "Flags:=", "Wireframe#",')
    add('     "MaterialValue:=", "\\"vacuum\\"", "SolveInside:=", True])')
    add('oModule = oDesign.GetModule("BoundarySetup")')
    add('oModule.AssignRadiation(["NAME:Rad1", "Objects:=", ["Region"],')
    add('    "IsFssReference:=", False, "IsForPML:=", False])')
    add("")
    add("# ---- solution setup -------------------------------------------------")
    add('oModule = oDesign.GetModule("AnalysisSetup")')
    add('oModule.InsertSetup("HfssDriven", ["NAME:Setup1",')
    add(f'    "Frequency:=", "{mesh:.6g}GHz", "MaxDeltaS:=", 0.02,')
    add('    "MaximumPasses:=", 12, "MinimumPasses:=", 2])')
    add('oModule.InsertFrequencySweep("Setup1", ["NAME:Sweep",')
    add('    "IsEnabled:=", True, "RangeType:=", "LinearCount",')
    add(f'    "RangeStart:=", "{lo:.6g}GHz", "RangeEnd:=", "{hi:.6g}GHz",')
    add('    "RangeCount:=", 401, "Type:=", "Interpolating", "SaveFields:=", False])')
    add("")
    add("# ---- far field infinite sphere --------------------------------------")
    add('oModule = oDesign.GetModule("RadField")')
    add('oModule.InsertInfiniteSphereSetup(["NAME:Infinite Sphere1",')
    add('    "UseCustomRadiationSurface:=", False,')
    add('    "ThetaStart:=", "0deg", "ThetaStop:=", "180deg", "ThetaStep:=", "1deg",')
    add('    "PhiStart:=", "0deg", "PhiStop:=", "360deg", "PhiStep:=", "2deg"])')
    add("")
    add('oDesktop.AddMessage("", "", 0, "OTA Hub: model built. Review notes in the header.")')
    return "\n".join(out) + "\n"


def _render_solid(solid) -> str:
    material = _material_name(solid.material)
    solve_inside = "False" if material == '"pec"' else "True"
    if isinstance(solid, Brick):
        dx = solid.x[1] - solid.x[0]
        dy = solid.y[1] - solid.y[0]
        dz = solid.z[1] - solid.z[0]
        flat = [ax for ax, d in (("X", dx), ("Y", dy), ("Z", dz)) if abs(d) < 1e-15]
        if flat:                # a zero-thickness sheet, normal to that axis
            # AEDT's rectangle spans the next two axes cyclically: Z -> (X, Y),
            # X -> (Y, Z), Y -> (Z, X). Every sheet used to be drawn normal to Z,
            # which gave a shorting wall (normal to Y) a height of zero.
            normal = flat[0]
            width, height = {"Z": (dx, dy), "X": (dy, dz), "Y": (dz, dx)}[normal]
            return "\n".join([
                "oEditor.CreateRectangle([",
                '    "NAME:RectangleParameters", "IsCovered:=", True,',
                f'    "XStart:=", "{_mm(solid.x[0])}", "YStart:=", "{_mm(solid.y[0])}",',
                f'    "ZStart:=", "{_mm(solid.z[0])}",',
                f'    "Width:=", "{_mm(width)}", "Height:=", "{_mm(height)}",',
                f'    "WhichAxis:=", "{normal}"],',
                f'    ["NAME:Attributes", "Name:=", "{solid.name}",',
                f'     "MaterialValue:=", {material}, "SolveInside:=", False])',
            ])
        return "\n".join([
            "oEditor.CreateBox([",
            '    "NAME:BoxParameters",',
            f'    "XPosition:=", "{_mm(solid.x[0])}", "YPosition:=", "{_mm(solid.y[0])}",',
            f'    "ZPosition:=", "{_mm(solid.z[0])}",',
            f'    "XSize:=", "{_mm(dx)}", "YSize:=", "{_mm(dy)}", "ZSize:=", "{_mm(dz)}"],',
            f'    ["NAME:Attributes", "Name:=", "{solid.name}",',
            f'     "MaterialValue:=", {material}, "SolveInside:=", {solve_inside}])',
        ])
    if isinstance(solid, Cylinder):
        axis = solid.axis.upper()
        pos = {"Z": (solid.centre[0], solid.centre[1], solid.span[0]),
               "X": (solid.span[0], solid.centre[0], solid.centre[1]),
               "Y": (solid.centre[0], solid.span[0], solid.centre[1])}[axis]
        height = solid.span[1] - solid.span[0]
        lines = [
            "oEditor.CreateCylinder([",
            '    "NAME:CylinderParameters",',
            f'    "XCenter:=", "{_mm(pos[0])}", "YCenter:=", "{_mm(pos[1])}",',
            f'    "ZCenter:=", "{_mm(pos[2])}",',
            f'    "Radius:=", "{_mm(solid.radius)}", "Height:=", "{_mm(height)}",',
            f'    "WhichAxis:=", "{axis}", "NumSides:=", "0"],',
            f'    ["NAME:Attributes", "Name:=", "{solid.name}",',
            f'     "MaterialValue:=", {material}, "SolveInside:=", {solve_inside}])',
        ]
        if solid.rotate_z:
            lines += [
                "oEditor.Rotate([",
                f'    "NAME:Selections", "Selections:=", "{solid.name}", "NewPartsModelFlag:=", "Model"],',
                f'    ["NAME:RotateParameters", "RotateAxis:=", "Z", "RotateAngle:=", "{solid.rotate_z:.9g}deg"])',
            ]
        return "\n".join(lines)
    if isinstance(solid, Cone):
        axis = solid.axis.upper()
        pos = {"Z": (solid.centre[0], solid.centre[1], solid.span[0]),
               "X": (solid.span[0], solid.centre[0], solid.centre[1]),
               "Y": (solid.centre[0], solid.span[0], solid.centre[1])}[axis]
        height = solid.span[1] - solid.span[0]
        return "\n".join([
            "oEditor.CreateCone([",
            '    "NAME:ConeParameters",',
            f'    "XCenter:=", "{_mm(pos[0])}", "YCenter:=", "{_mm(pos[1])}",',
            f'    "ZCenter:=", "{_mm(pos[2])}",',
            f'    "BottomRadius:=", "{_mm(solid.radius_start)}",',
            f'    "TopRadius:=", "{_mm(solid.radius_end)}",',
            f'    "Height:=", "{_mm(height)}", "WhichAxis:=", "{axis}"],',
            f'    ["NAME:Attributes", "Name:=", "{solid.name}",',
            f'     "MaterialValue:=", {material}, "SolveInside:=", {solve_inside}])',
        ])
    if isinstance(solid, Torus):
        return "\n".join([
            "oEditor.CreateTorus([",
            '    "NAME:TorusParameters",',
            f'    "XCenter:=", "{_mm(solid.centre[0])}", '
            f'"YCenter:=", "{_mm(solid.centre[1])}",',
            f'    "ZCenter:=", "{_mm(solid.centre[2])}",',
            f'    "MajorRadius:=", "{_mm(solid.major_radius)}",',
            f'    "MinorRadius:=", "{_mm(solid.minor_radius)}",',
            f'    "WhichAxis:=", "{solid.axis.upper()}"],',
            f'    ["NAME:Attributes", "Name:=", "{solid.name}",',
            f'     "MaterialValue:=", {material}, "SolveInside:=", {solve_inside}])',
        ])
    if isinstance(solid, Sphere):
        return "\n".join([
            "oEditor.CreateSphere([",
            '    "NAME:SphereParameters",',
            f'    "XCenter:=", "{_mm(solid.centre[0])}", '
            f'"YCenter:=", "{_mm(solid.centre[1])}",',
            f'    "ZCenter:=", "{_mm(solid.centre[2])}",',
            f'    "Radius:=", "{_mm(solid.radius)}"],',
            f'    ["NAME:Attributes", "Name:=", "{solid.name}",',
            f'     "MaterialValue:=", {material}, "SolveInside:=", {solve_inside}])',
        ])
    return f"# unsupported solid type {type(solid).__name__}"


def _pec_sheets(model: Model) -> list[str]:
    """The PEC solids of zero thickness that survive the booleans."""
    gone = set()
    for op in model.operations:
        if isinstance(op, Unite) or (isinstance(op, Subtract) and not op.keep_tools):
            gone.update(op.tools)
    out = []
    for s in model.solids:
        if s.material != "PEC" or s.name in gone:
            continue
        if isinstance(s, Brick) and min(abs(s.x[1] - s.x[0]), abs(s.y[1] - s.y[0]), abs(s.z[1] - s.z[0])) < 1e-15:
            out.append(s.name)
        elif isinstance(s, Cylinder) and abs(s.span[1] - s.span[0]) < 1e-15:
            out.append(s.name)
    return out


def _render_operation(op) -> str:
    if isinstance(op, Subtract):
        tools = ",".join(op.tools)
        return "\n".join([
            "# boolean: remove the tool solids from the target"
            + (" (the tools stay)" if op.keep_tools else ""),
            "oEditor.Subtract([",
            '    "NAME:Selections",',
            f'    "Blank Parts:=", "{op.target}", "Tool Parts:=", "{tools}"],',
            f'    ["NAME:SubtractParameters", "KeepOriginals:=", {op.keep_tools}])',
        ])
    if isinstance(op, Unite):
        return "\n".join([
            "# boolean: join the touching conductors into one",
            "oEditor.Unite([",
            f'    "NAME:Selections", "Selections:=", "{",".join((op.target,) + op.tools)}"],',
            '    ["NAME:UniteParameters", "KeepOriginals:=", False])',
        ])
    return f"# unsupported operation {type(op).__name__}"


def _port_sheet(port: DiscretePort) -> tuple[list[float], float, str]:
    """(minimum corner, side, normal axis) of a square sheet that holds the port's
    integration line: the line's own axis and one beside it, the side as long as
    the gap. Square, so it does not matter which in-plane axis AEDT runs Width
    along. The sheet used to be normal to Y whatever the feed: a horizontal feed
    got zero height, a y-directed one a sheet the line did not lie in."""
    d = [e - s for s, e in zip(port.start, port.end)]
    along = max(range(3), key=lambda i: abs(d[i]))
    g = abs(d[along])
    beside = 0 if along != 0 else 1                 # x beside a z or y line, y beside an x line
    normal = ({0, 1, 2} - {along, beside}).pop()
    corner = list(port.start)
    corner[along] = min(port.start[along], port.end[along])
    corner[beside] = port.start[beside] - g / 2
    return corner, g, "XYZ"[normal]


def _render_port(port: DiscretePort, index: int) -> str:
    corner, side, normal = _port_sheet(port)
    return "\n".join([
        f'# lumped port {index}: {port.impedance:.6g} ohm',
        "oEditor.CreateRectangle([",
        '    "NAME:RectangleParameters", "IsCovered:=", True,',
        f'    "XStart:=", "{_mm(corner[0])}", "YStart:=", "{_mm(corner[1])}",',
        f'    "ZStart:=", "{_mm(corner[2])}",',
        f'    "Width:=", "{_mm(side)}", "Height:=", "{_mm(side)}",',
        f'    "WhichAxis:=", "{normal}"],',
        f'    ["NAME:Attributes", "Name:=", "{port.name}_sheet",',
        '     "MaterialValue:=", "\\"vacuum\\"", "SolveInside:=", True])',
        'oModule = oDesign.GetModule("BoundarySetup")',
        f'oModule.AssignLumpedPort(["NAME:{port.name}",',
        f'    "Objects:=", ["{port.name}_sheet"], "RenormalizeAllTerminals:=", True,',
        '    "DoDeembed:=", False,',
        f'    ["NAME:Modes", ["NAME:Mode1", "ModeNum:=", 1, "UseIntLine:=", True,',
        f'     ["NAME:IntLine", "Start:=", ["{_mm(port.start[0])}", '
        f'"{_mm(port.start[1])}", "{_mm(port.start[2])}"],',
        f'      "End:=", ["{_mm(port.end[0])}", "{_mm(port.end[1])}", '
        f'"{_mm(port.end[2])}"]],',
        f'     "CharImp:=", "Zpi", "RenormImp:=", "{port.impedance:.6g}ohm"]],',
        # the port's own impedance; RenormImp only renormalises the reported S
        f'    "Impedance:=", "{port.impedance:.6g}ohm"])',
    ])


def _material_name(material: str) -> str:
    if material == "PEC":
        return '"pec"'
    if material in ("VACUUM", "VOID"):
        return '"vacuum"'
    if material.startswith("eps_r="):
        return '"%s"' % dielectric_name(material)
    return f'"{material}"'


def _mm(value: float) -> str:
    return f"{value * _MM:.6f}mm"


def _wrap(text: str, width: int = 74) -> list[str]:
    import textwrap
    return textwrap.wrap(text, width) or [""]
