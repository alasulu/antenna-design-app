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
