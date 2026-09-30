"""Synthesis-solver behaviour, independent of any particular spec file."""
from __future__ import annotations

import pytest

from otahub.core.archetype import Archetype, SynthesisError
from otahub.core.spec import ArchetypeSpec


def _spec(**over) -> ArchetypeSpec:
    base = dict(
        key="demo", name="Demo", family="test",
        freq_range_hz=(1e6, 1e12),
        parameters=(),
        synthesis=(), analysis=(), known_cases=(), references=("test",),
    )
    base.update(over)
    return ArchetypeSpec(**base)


def _from_json(d) -> Archetype:
    return Archetype(ArchetypeSpec.from_json(d, "test"))


DEMO = {
    "key": "demo", "name": "Demo", "freq_range_hz": [1e6, 1e12],
    "parameters": [
        {"symbol": "f0", "role": "requirement", "unit": "Hz"},
        {"symbol": "sigma", "role": "material", "unit": "S/m", "typical": "5.8e7"},
        {"symbol": "eta_target", "role": "requirement", "typical": "0.5"},
    ],
    "synthesis": [
        {"output": "L", "expr": "0.5 * c / f0", "units_out": "m"},
        {"output": "b", "expr": "L * eta_target / 100", "units_out": "m"},
    ],
    "analysis": [{"metric": "rs", "expr": "np.sqrt(pi * f0 * mu0 / sigma)", "units_out": "ohm"}],
    "references": ["unit test"],
}


def test_partial_synthesis_keeps_what_it_can_compute():
    """An unresolvable optional rule must not discard the whole design.

    This is the failure that made 79 of 100 spec known-cases fail: a case
    asking only for length was thrown away because an unrelated wire-radius
    rule lacked an efficiency target.
    """
    d = _from_json(DEMO).synthesize(f0=300e6)
    assert d.parameters["L"] == pytest.approx(0.49965, rel=1e-3)
    assert "b" in d.unresolved and "eta_target" in d.unresolved["b"]
    assert any("eta_target" in w for w in d.warnings)


def test_strict_mode_refuses_an_incomplete_design():
    with pytest.raises(SynthesisError, match="incomplete design"):
        _from_json(DEMO).synthesize(strict=True, f0=300e6)


def test_supplying_the_missing_requirement_completes_the_design():
    d = _from_json(DEMO).synthesize(f0=300e6, eta_target=0.5)
    assert d.unresolved == {}
    assert d.parameters["b"] == pytest.approx(0.49965 * 0.5 / 100, rel=1e-3)


def test_material_parameters_default_but_requirements_never_do():
    """Copper is a defensible default; a design requirement is not."""
    d = _from_json(DEMO).synthesize(f0=300e6)
    assert d.parameters["sigma"] == pytest.approx(5.8e7)   # material defaulted
    assert "b" in d.unresolved                              # requirement did not


def test_explicit_material_value_overrides_the_default():
    d = _from_json(DEMO).synthesize(f0=300e6, sigma=3.54e7)
    assert d.metrics["rs"] == pytest.approx(
        (3.14159265358979 * 300e6 * 1.25663706212e-6 / 3.54e7) ** 0.5, rel=1e-6)


def test_frequency_outside_the_stated_band_warns_rather_than_silently_computing():
    a = _from_json({**DEMO, "freq_range_hz": [1e9, 2e9]})
    d = a.synthesize(f0=300e6, eta_target=0.5)
    assert any("outside this archetype's stated validity band" in w for w in d.warnings)


def test_lambda_and_k_are_derived_automatically():
    d = _from_json(DEMO).synthesize(f0=300e6, eta_target=0.5)
    assert d.parameters["lambda0"] == pytest.approx(0.99931, rel=1e-4)
    assert d.parameters["k0"] == pytest.approx(6.2872, rel=1e-4)


def test_expectation_lookup_tolerates_unit_tagged_keys():
    from otahub.core.archetype import _lookup
    d = _from_json(DEMO).synthesize(f0=300e6, eta_target=0.5)
    assert _lookup(d, "L_m") == pytest.approx(d.parameters["L"])
    assert _lookup(d, "L") == pytest.approx(d.parameters["L"])
    assert _lookup(d, "nonexistent_quantity_xyz") is None


def test_absolute_tolerance_is_used_for_expectations_of_zero():
    """Relative error is undefined against zero; tol_abs is the only sane test.

    A resonant dipole left with 3 ohm of residual reactance against an
    expected 0 is an excellent result, but relative comparison scores it as a
    300% failure.
    """
    import math
    from otahub.core.archetype import _compare
    assert _compare(3.0, 0.0, tol_pct=1.0) == (False, math.inf)
    passed, _ = _compare(3.0, 0.0, tol_pct=1.0, tol_abs=6.0)
    assert passed
    passed, _ = _compare(9.0, 0.0, tol_pct=1.0, tol_abs=6.0)
    assert not passed


def test_absolute_tolerance_overrides_a_tight_relative_one():
    from otahub.core.archetype import _compare
    assert _compare(101.0, 100.0, tol_pct=0.1)[0] is False
    assert _compare(101.0, 100.0, tol_pct=0.1, tol_abs=2.0)[0] is True


def test_complex_impedance_uses_magnitude_of_the_difference():
    from otahub.core.archetype import _compare
    passed, err = _compare(complex(73.0, 42.0), complex(73.1, 42.5), tol_pct=2.0)
    assert passed and err < 2.0


def test_a_supplied_value_outranks_the_specs_own_default(registry):
    """The quietest kind of wrong: asking for one antenna and being answered
    about another.

    Many specs synthesise a nominal geometry - a halo's gap, a loop's
    circumference - so the archetype is usable from a frequency alone. The
    resolver used to run those rules unconditionally, which meant a caller who
    supplied the value watched it be overwritten by the default and got a
    design sheet for a different antenna, with no warning anywhere. Every
    number on the sheet was self-consistent; it just answered a question
    nobody asked.
    """
    lam = 2.99792458e8 / 1.46e8
    b = 0.002 * lam          # a wire radius, so the design resolves fully
    default = registry["halo_loop"].synthesize(f0=1.46e8, b=b)
    asked = registry["halo_loop"].synthesize(f0=1.46e8, b=b, g=0.03 * lam)
    assert default.get("g") == pytest.approx(0.015 * lam, rel=1e-6)
    assert asked.get("g") == pytest.approx(0.03 * lam, rel=1e-12)
    # and it must actually propagate, not merely be recorded
    assert asked.get("resonant_circumference_m") != pytest.approx(
        default.get("resonant_circumference_m"), rel=1e-4)


def test_overriding_still_leaves_everything_else_derived(registry):
    """Overriding one output must not strand the rules that depend on it."""
    lam = 2.99792458e8 / 1.46e8
    d = registry["halo_loop"].synthesize(f0=1.46e8, g=0.02 * lam, b=0.002 * lam)
    assert not d.unresolved, d.unresolved
    assert d.get("Lc") + d.get("g") == pytest.approx(
        d.get("resonant_circumference_m"), rel=1e-12)


def test_a_tiny_expectation_is_still_compared_relatively():
    """Below 1e-15 the comparison used to divide by 1, so 1e-12 passed a 10%
    check against an expected 1e-16 - ten thousand times too large."""
    from otahub.core.archetype import _compare
    passed, err = _compare(1e-12, 1e-16, tol_pct=10.0)
    assert not passed and err == pytest.approx(999900.0)
    assert _compare(1.05e-16, 1e-16, tol_pct=10.0)[0]


def test_strict_synthesis_refuses_a_rule_that_failed():
    """A rule that raises (1/0 here) used to drop out of the design, so strict
    mode returned L = None with nothing unresolved."""
    a = _from_json({**DEMO, "synthesis": [{"output": "L", "expr": "1/0", "units_out": "m"}]})
    loose = a.synthesize(f0=1e9)
    assert loose.get("L") is None and any("'L' failed" in w for w in loose.warnings)
    with pytest.raises(SynthesisError, match="L failed"):
        a.synthesize(strict=True, f0=1e9)


def test_a_band_design_outside_the_stated_band_warns_too():
    """Only f0 was checked; a design asked for with f_low alone never warned."""
    a = _from_json({**DEMO, "parameters": DEMO["parameters"] + [
        {"symbol": "f_low", "role": "requirement", "unit": "Hz"}],
        "synthesis": [{"output": "L", "expr": "0.5 * c / f_low", "units_out": "m"}],
        "analysis": []})
    assert any("f_low=1000 Hz is outside" in w for w in a.synthesize(f_low=1e3).warnings)
    assert not any("outside" in w for w in a.synthesize(f_low=1e9).warnings)


def test_malformed_formula_is_reported_not_raised():
    """`doctor` crashed on a syntax error: problems() caught ExprError only."""
    spec = ArchetypeSpec.from_json({**DEMO, "synthesis": [{"output": "L", "expr": "1 +"}],
                                    "known_cases": [{"given": {"f0": 1e9}, "expect": {"L": 1.0}}]},
                                   "test")
    assert any("cannot parse" in p for p in spec.problems())


@pytest.mark.parametrize("expect,why", [
    ({}, "expects nothing"),
    ({"L": "long"}, "not a number"),
    ({"L": 0.0}, "without tol_abs"),
])
def test_a_known_case_that_cannot_fail_is_a_problem(expect, why):
    spec = ArchetypeSpec.from_json({**DEMO, "known_cases": [{"given": {"f0": 1e9}, "expect": expect}]},
                                   "test")
    assert any(why in p for p in spec.problems())
