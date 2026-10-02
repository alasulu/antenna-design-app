"""Parametric export: the CST and HFSS scripts write the geometry as expressions in
the variables they declare, so that editing a variable in the simulator - by hand,
or by its optimiser - moves the model the way re-running the builder would.

Checked on the scripts themselves, read back by a small evaluator of this file's
own (Python's ast under a whitelist: numbers, names, + - * /, unary minus and the
four functions the scripts use, with HFSS's unit suffixes turned into SI factors).
Nothing of the exporter's is used to read its output.
"""
from __future__ import annotations

import ast
import dataclasses
import math
import operator
import re

import pytest

from otahub.export import build, cst, hfss
from otahub.export.base import BUILDERS
from tests.test_export import CASES

_OPS = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul, ast.Div: operator.truediv}
_CALLS = {"Sqr": math.sqrt, "sqrt": math.sqrt, "Sin": math.sin, "sin": math.sin,
          "Cos": math.cos, "cos": math.cos, "Abs": abs, "abs": abs}
#: HFSS's unit suffixes as SI factors, under names no variable has
_UNITS = {"__mm": 1e-3, "__GHz": 1e9, "__deg": math.pi / 180.0}


def evaluate(expr: str, variables: dict[str, float]) -> float:
    """One field of a script. Anything outside the scripts' arithmetic fails, and a
    name the script never declared is a KeyError."""
    def ev(node):
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return float(node.value)
        if isinstance(node, ast.Name):
            return variables[node.id]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return -ev(node.operand)
        if isinstance(node, ast.BinOp) and type(node.op) in _OPS:
            return _OPS[type(node.op)](ev(node.left), ev(node.right))
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id in _CALLS and len(node.args) == 1 and not node.keywords):
            return _CALLS[node.func.id](ev(node.args[0]))
        raise ValueError(f"{expr!r}: {ast.dump(node)} is not script arithmetic")
    return ev(ast.parse(expr, mode="eval").body)


def names(expr: str) -> set[str]:
    return {n.id for n in ast.walk(ast.parse(expr, mode="eval")) if isinstance(n, ast.Name)} - set(_UNITS)


_CST_FIELD = re.compile(r"^\s+\.(Xrange|Yrange|Zrange|Outerradius|Innerradius|Xcenter|Ycenter|Zcenter"
                        r"|Bottomradius|Topradius|CenterRadius|Center|SetP1|SetP2|Angle) (.*)$")


def cst_script(text: str):
    """(declared parameters, geometry fields) of a CST macro: lengths in mm, the
    fields as (keyword, expression) in the order written."""
    declared = {m[1]: float(m[2]) for m in re.finditer(r'StoreParameter "([^"]+)", (\S+)', text)}
    fields = []
    for line in text.splitlines():
        match = _CST_FIELD.match(line)
        if match:
            fields += [(match[1], f) for f in re.findall(r'"([^"]*)"', match[2]) if f not in ("False", "True")]
    return declared, fields


_HFSS_FIELD = re.compile(
    r'"(XPosition|YPosition|ZPosition|XSize|YSize|ZSize|XStart|YStart|ZStart|Width|Height|XCenter'
    r'|YCenter|ZCenter|Radius|BottomRadius|TopRadius|MajorRadius|MinorRadius|RotateAngle):=", "([^"]*)"'
    r'|"(Start|End):=", \[([^\]]*)\]')


def _si(text: str) -> str:
    """HFSS's 1.5mm, 2.4GHz, 30deg as Python: (1.5*__mm)."""
    return re.sub(r"(?<![\w.])(\d+(?:\.\d+)?)(mm|GHz|deg)\b", r"(\1*__\2)", text)


def hfss_script(text: str):
    """(declared variables, geometry fields) of an HFSS script, in SI."""
    block = text[text.index("variables = ["):text.index("for vname, vvalue in variables")]
    declared = {name: evaluate(_si(value), _UNITS)
                for name, value in re.findall(r'\("([^"]+)", "([^"]+)"\)', block)}
    fields = []
    for m in _HFSS_FIELD.finditer(text[text.index("# ---- design variables"):]):
        if m[1]:
            fields.append((m[1], _si(m[2])))
        else:                                   # an integration line's end point
            fields += [(m[3], _si(p)) for p in re.findall(r'"([^"]*)"', m[4])]
    return declared, fields


def cst_value(text: str, expr: str) -> float:
    """An expression from a CST macro, with the parameters the macro declares (mm)."""
    return evaluate(expr, cst_script(text)[0])


def hfss_value(text: str, expr: str) -> float:
    """An expression from an HFSS script, with the variables it declares (SI)."""
    return evaluate(_si(expr), {**_UNITS, **hfss_script(text)[0]})


#: backend -> (render, read, absolute tolerance in its length unit: 1e-12 m)
BACKENDS = {"cst": (cst.render, cst_script, 1e-9), "hfss": (hfss.render, hfss_script, 1e-12)}


def numeric(model):
    """The same model with every expression dropped, as the scripts were written
    before: the builder's numbers."""
    def plain(v):
        return tuple(plain(x) for x in v) if isinstance(v, tuple) else float(v) if isinstance(v, float) else v

    def strip(item):
        return dataclasses.replace(item, **{f.name: plain(getattr(item, f.name)) for f in dataclasses.fields(item)})
    return dataclasses.replace(model, solids=[strip(s) for s in model.solids],
                               ports=[strip(p) for p in model.ports])


def values(fields, declared):
    return [(key, evaluate(expr, {**_UNITS, **declared})) for key, expr in fields]


def mismatches(got, want, absolute):
    assert [k for k, _ in got] == [k for k, _ in want], "the two scripts differ in structure"
    return [(i, k, g, w) for i, ((k, g), (_, w)) in enumerate(zip(got, want))
            if abs(g - w) > 1e-9 * max(abs(g), abs(w)) + absolute]


def _design(registry, key):
    """The first of the spec's known cases that builds geometry, else test_export's case."""
    archetype = registry[key]
    for case in archetype.spec.known_cases:
        design = archetype.synthesize(**case.given)
        if build(design).built_geometry:
            return design
    return archetype.synthesize(**CASES[key])


@pytest.fixture(scope="module")
def designs(registry):
    return {key: _design(registry, key) for key in sorted(BUILDERS)}


@pytest.mark.parametrize("backend", sorted(BACKENDS))
def test_every_expression_reproduces_the_builders_numbers(designs, backend):
    """Every coordinate, radius, span, port end and rotation the script writes,
    evaluated with the values it declares, is the builder's number - to 1e-9, or
    1e-12 m for a zero."""
    render, read, absolute = BACKENDS[backend]
    bad, parametric = {}, {}
    for key, design in designs.items():
        model = build(design)
        assert model.built_geometry, key
        declared, fields = read(render(model))
        _, plain = read(render(numeric(model)))
        assert not any(names(e) for _, e in plain), f"{key}: the numeric script names a variable"
        wrong = mismatches(values(fields, declared), values(plain, {}), absolute)
        if wrong:
            bad[key] = wrong[:3]
        parametric[key] = sum(1 for _, e in fields if names(e))
    assert not bad
    # and the expressions are really there: every builder writes some of its geometry in variables
    assert all(parametric.values()), [k for k, n in parametric.items() if not n]


def _with(design, name, value):
    """The design with one exported value changed and nothing re-derived from it."""
    def edit(values):
        return {**values, name: value} if name in values else dict(values)
    return dataclasses.replace(design, requirements=edit(design.requirements),
                               parameters=edit(design.parameters))


def _structure(model):
    return ([(type(s).__name__, s.name) for s in model.solids], [p.name for p in model.ports],
            model.operations)


@pytest.mark.parametrize("backend", sorted(BACKENDS))
def test_editing_a_variable_moves_the_geometry_as_the_builder_would(designs, backend):
    """Change one declared variable in the script and every expression re-evaluated
    with it must give the builder's own geometry for that change.

    Re-synthesising from a changed requirement would re-derive other parameters
    too (a new W moves L), which the script cannot know about - in the simulator
    the other variables stay put. So the reference is the BUILDER run on the
    design with that one exported value changed by 2% and every other held at its
    exported value: the editing is the same on both sides, and only the geometry
    is compared. Every exported value is tried - a value the builder reads but the
    expressions do not follow shows up as a mismatch. A change that alters the
    structure (a count rounding to another integer, a branch taken differently
    such as a dielectric spacer appearing) cannot be followed by editing a
    number, and is skipped."""
    render, read, absolute = BACKENDS[backend]
    bad, moved = {}, {}
    for key, design in designs.items():
        model = build(design)
        declared, fields = read(render(model))
        moved[key] = []
        for name, value in sorted(model.parameters.items()):
            if name == "lambda0" or value == 0:           # not declared; scaling 0 changes nothing
                continue
            try:
                edited = build(_with(design, name, value * 1.02))
            except ValueError:                             # the builder refuses it (a probe off its patch)
                continue
            if not edited.built_geometry or _structure(edited) != _structure(model):
                continue
            changed, _ = read(render(edited))
            differ = [k for k in declared if declared[k] != changed[k]]   # the one edited variable
            assert len(differ) <= 1 and all(k.startswith(name) for k in differ), (key, name, differ)
            _, plain = read(render(numeric(edited)))
            got, want = values(fields, changed), values(plain, {})
            wrong = mismatches(got, want, absolute)
            if wrong:
                bad[f"{key}:{name}"] = wrong[:3]
            if mismatches(got, values(fields, declared), absolute):
                moved[key].append(name)
    assert not bad
    # the check is not vacuous: in every model some variable moves something
    assert all(moved.values()), [k for k, v in moved.items() if not v]


def _edit_cst(text: str, label: str, factor: float) -> str:
    """The macro with one StoreParameter's value scaled, as a user would edit it."""
    def scale(m):
        return f'StoreParameter "{label}", {float(m[1]) * factor!r}'
    return re.sub(rf'StoreParameter "{re.escape(label)}", (\S+)', scale, text, count=1)


def _cst_field(text: str, solid: str, keyword: str, declared) -> list[float]:
    block = text[text.index(f'.Name "{solid}"'):]
    block = block[:block.index("End With")]
    line = next(l for l in block.splitlines() if l.strip().startswith(f".{keyword} "))
    return [evaluate(f, declared) for f in re.findall(r'"([^"]*)"', line)]


def test_a_patch_follows_its_width_in_cst(registry):
    """Edit W_mm in the macro: the patch's x range is +-W/2 of the new width, and the
    board, which the builder sizes from the larger of W and L, grows with it."""
    design = registry["rectangular_patch"].synthesize(f0=2.4e9, eps_r=4.4, h=1.6e-3)
    text = _edit_cst(cst.render(build(design)), "W_mm", 1.5)
    declared, _ = cst_script(text)
    w, length = design.get("W") * 1.5 * 1e3, design.get("L") * 1e3
    assert _cst_field(text, "patch", "Xrange", declared) == pytest.approx([-w / 2, w / 2], rel=1e-12)
    margin = max(w, length) * 0.6
    assert _cst_field(text, "substrate", "Xrange", declared) == pytest.approx(
        [-(w / 2 + margin), w / 2 + margin], rel=1e-12)


def test_a_dipole_follows_its_length_in_hfss(registry):
    """Edit L in the HFSS variables: the upper arm's tip is at the new L/2."""
    design = registry["half_wave_dipole"].synthesize(f0=300e6, aw=1e-3)
    text = hfss.render(build(design))
    declared, _ = hfss_script(text)
    declared["L"] *= 0.9
    block = text[:text.index('"Name:=", "arm_upper"')]
    block = block[block.rindex("CreateCylinder"):]
    field = dict(re.findall(r'"(ZCenter|Height):=", "([^"]*)"', block))
    tip = sum(evaluate(_si(field[k]), {**_UNITS, **declared}) for k in ("ZCenter", "Height"))
    assert tip == pytest.approx(design.get("L") * 0.9 / 2, rel=1e-12)


def test_a_loop_and_a_guide_follow_their_sizes(registry):
    """A circular loop's ring follows its radius a; an open-ended guide's aperture
    and length follow a_wg."""
    loop = cst.render(build(registry["one_wavelength_circular_loop"].synthesize(f0=300e6, b=0.001)))
    loop = _edit_cst(loop, "a_mm", 1.25)
    declared, _ = cst_script(loop)
    outer, inner = _cst_field(loop, "loop", "Outerradius", declared)[0], _cst_field(loop, "loop", "Innerradius", declared)[0]
    assert (outer + inner) / 2 == pytest.approx(declared["a_mm"], rel=1e-12)
    guide = cst.render(build(registry["open_ended_waveguide"].synthesize(f0=10e9, a_wg=0.02286, b_wg=0.01016)))
    guide = _edit_cst(guide, "a_wg_mm", 1.1)
    declared, _ = cst_script(guide)
    a = declared["a_wg_mm"]
    assert _cst_field(guide, "guide_interior", "Xrange", declared) == pytest.approx([-a / 2, a / 2], rel=1e-12)
    assert _cst_field(guide, "guide_interior", "Zrange", declared) == pytest.approx([0.0, 2 * a], rel=1e-12)
