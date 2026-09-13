"""Ansys HFSS script generation.

Emits an IronPython script for the HFSS scripting interface. Run it from
Tools > Run Script inside Electronics Desktop, or with ansysedt -RunScript.
Dimensions become design variables so the model stays parametric.
"""
from __future__ import annotations

from .base import Brick, Cylinder, DiscretePort, Model
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
            eps = material.split("=", 1)[1]
            name = "substrate_eps" + eps.replace(".", "p")
            add("oDefinitionManager.AddMaterial([")
            add(f'    "NAME:{name}", "CoordinateSystemType:=", "Cartesian",')
            add(f'    "permittivity:=", "{eps}",')
            add('    "dielectric_loss_tangent:=", "0.02"])')
        add("")

    add("# ---- geometry ------------------------------------------------------")
    for solid in model.solids:
        add(_render_solid(solid))
        add("")

    if model.ports:
        add("# ---- excitation ----------------------------------------------------")
        for index, port in enumerate(model.ports, start=1):
            add(_render_port(port, index))
        add("")

    freq = model.frequency_hz / 1e9
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
    add(f'    "Frequency:=", "{freq:.6g}GHz", "MaxDeltaS:=", 0.02,')
    add('    "MaximumPasses:=", 12, "MinimumPasses:=", 2])')
    add('oModule.InsertFrequencySweep("Setup1", ["NAME:Sweep",')
    add('    "IsEnabled:=", True, "Type:=", "Interpolating",')
    add(f'    "StartValue:=", "{freq * 0.7:.6g}GHz", "StopValue:=", "{freq * 1.3:.6g}GHz",')
    add('    "Count:=", 401, "SaveFields:=", False])')
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
        if abs(dz) < 1e-15:     # a zero-thickness sheet
            return "\n".join([
                "oEditor.CreateRectangle([",
                '    "NAME:RectangleParameters", "IsCovered:=", True,',
                f'    "XStart:=", "{_mm(solid.x[0])}", "YStart:=", "{_mm(solid.y[0])}",',
                f'    "ZStart:=", "{_mm(solid.z[0])}",',
                f'    "Width:=", "{_mm(dx)}", "Height:=", "{_mm(dy)}",',
                '    "WhichAxis:=", "Z"],',
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
        return "\n".join([
            "oEditor.CreateCylinder([",
            '    "NAME:CylinderParameters",',
            f'    "XCenter:=", "{_mm(pos[0])}", "YCenter:=", "{_mm(pos[1])}",',
            f'    "ZCenter:=", "{_mm(pos[2])}",',
            f'    "Radius:=", "{_mm(solid.radius)}", "Height:=", "{_mm(height)}",',
            f'    "WhichAxis:=", "{axis}", "NumSides:=", "0"],',
            f'    ["NAME:Attributes", "Name:=", "{solid.name}",',
            f'     "MaterialValue:=", {material}, "SolveInside:=", {solve_inside}])',
        ])
    return f"# unsupported solid type {type(solid).__name__}"


def _render_port(port: DiscretePort, index: int) -> str:
    return "\n".join([
        f'# lumped port {index}: {port.impedance:.6g} ohm',
        "oEditor.CreateRectangle([",
        '    "NAME:RectangleParameters", "IsCovered:=", True,',
        f'    "XStart:=", "{_mm(port.start[0])}", "YStart:=", "{_mm(port.start[1])}",',
        f'    "ZStart:=", "{_mm(port.start[2])}",',
        f'    "Width:=", "{_mm(max(abs(port.end[0] - port.start[0]), 1e-4))}",',
        f'    "Height:=", "{_mm(port.end[2] - port.start[2] or port.end[1] - port.start[1])}",',
        '    "WhichAxis:=", "Y"],',
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
        f'     "CharImp:=", "Zpi", "RenormImp:=", "{port.impedance:.6g}ohm"]]])',
    ])


def _material_name(material: str) -> str:
    if material == "PEC":
        return '"pec"'
    if material in ("VACUUM", "VOID"):
        return '"vacuum"'
    if material.startswith("eps_r="):
        return '"substrate_eps%s"' % material.split("=", 1)[1].replace(".", "p")
    return f'"{material}"'


def _mm(value: float) -> str:
    return f"{value * _MM:.6f}mm"


def _wrap(text: str, width: int = 74) -> list[str]:
    import textwrap
    return textwrap.wrap(text, width) or [""]
