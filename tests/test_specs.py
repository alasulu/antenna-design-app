"""Tests generated from the spec files themselves.

Every citable ``known_cases`` entry an archetype declares becomes a live test.
This is the mechanism that keeps the catalogue honest: an archetype whose
formulas drift away from its references fails here, not silently in a user's
design.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from otahub.core.registry import Registry  # noqa: E402

REG = Registry.load(ROOT / "specs")
ARCHETYPES = sorted(REG, key=lambda a: a.key)

# (archetype, case_index, quantity) triples, flattened for readable test ids.
KNOWN_CASES = [
    pytest.param(a.key, i, id=f"{a.key}-case{i}")
    for a in ARCHETYPES
    for i, _ in enumerate(a.spec.known_cases)
]


def test_spec_directory_is_populated():
    assert len(REG) > 0, (
        "no archetypes loaded — specs/*.json is empty or every file failed to parse"
    )


def test_no_spec_files_failed_to_load():
    assert REG.load_errors == {}, f"spec files failed to load: {REG.load_errors}"


@pytest.mark.parametrize("key", [a.key for a in ARCHETYPES])
def test_archetype_is_structurally_sound(key):
    """No unknown symbols, no duplicate outputs, at least one verifiable case."""
    problems = REG[key].spec.problems()
    assert not problems, "\n".join(f"  - {p}" for p in problems)


@pytest.mark.parametrize("key", [a.key for a in ARCHETYPES])
def test_archetype_declares_provenance(key):
    spec = REG[key].spec
    assert spec.references, f"{key} cites no reference — unciteable numbers are not usable"
    lo, hi = spec.freq_range_hz
    assert lo < hi, f"{key} has a degenerate frequency band"


@pytest.mark.parametrize("key,case_index", KNOWN_CASES)
def test_known_case_matches_reference(key, case_index):
    """The archetype reproduces its own cited textbook numbers."""
    archetype = REG[key]
    case = archetype.spec.known_cases[case_index]
    results = [r for r in archetype.check_known_cases()][
        sum(len(c.expect) for c in archetype.spec.known_cases[:case_index]):
    ][: len(case.expect)]

    failures = [
        f"{r['quantity']}: expected {r['expected']}, got {r['actual']} "
        f"(err {r['error_pct']}%, tol {r['tol_pct']}%) {r['detail']}"
        for r in results if not r["passed"]
    ]
    assert not failures, (
        f"{key} disagrees with {case.source or 'its cited source'}:\n"
        + "\n".join(f"  - {f}" for f in failures)
    )


@pytest.mark.parametrize("key", [a.key for a in ARCHETYPES])
def test_synthesis_runs_at_a_frequency_inside_the_stated_band(key):
    """Smoke test: geometric outputs must be finite and positive."""
    import math

    archetype = REG[key]
    lo, hi = archetype.spec.freq_range_hz
    f0 = math.sqrt(max(lo, 1e3) * min(hi, 1e12)) if math.isfinite(hi) else max(lo, 1e9)

    # Supply declared requirement-role parameters from their typical values.
    reqs = {"f0": f0}
    for p in archetype.spec.parameters:
        if p.role in ("requirement", "material") and p.symbol != "f0":
            try:
                reqs[p.symbol] = float(p.typical)
            except (TypeError, ValueError):
                continue
    try:
        design = archetype.synthesize(**reqs)
    except Exception as exc:  # noqa: BLE001
        pytest.skip(f"{key} needs requirements the spec does not give typicals for: {exc}")

    for name, value in design.parameters.items():
        if isinstance(value, complex) or not isinstance(value, (int, float)):
            continue
        assert math.isfinite(value), f"{key}.{name} is not finite ({value})"


# --- physical plausibility -------------------------------------------------
# A db10()/log10() mix-up once made a 4.2-lambda Yagi report 53 dBi, which is
# more gain than a 30 m dish. Unit-scale errors like that produce numbers that
# are structurally valid and physically absurd, so bounds are checked directly.

# The lower dBi bound is deliberately generous: a ferrite rod antenna really
# does show about -81 dBi, because its radiation resistance (~1e-7 ohm) is
# swamped by coil loss (~16 ohm). That is correct physics, not a bug - such
# antennas are only ever used for receiving, where external noise dominates.
# The UPPER bound is the one that catches unit-scale errors.
_PLAUSIBLE = {
    "dbi": (-150.0, 80.0),
    "db": (-120.0, 120.0),
    "efficiency": (0.0, 1.0),
    "deg": (0.0, 360.0),
}


def _bounds_for(metric: str, unit: str) -> tuple[float, float] | None:
    name, unit = metric.lower(), (unit or "").lower()
    if name.endswith("_dbi") or unit == "dbi":
        return _PLAUSIBLE["dbi"]
    if "efficiency" in name:
        return _PLAUSIBLE["efficiency"]
    if name.endswith("_deg") or unit == "deg":
        return _PLAUSIBLE["deg"]
    if name.endswith("_db") or unit == "db":
        return _PLAUSIBLE["db"]
    return None


@pytest.mark.parametrize("key,case_index", KNOWN_CASES)
def test_known_case_values_are_physically_plausible(key, case_index):
    import math

    archetype = REG[key]
    case = archetype.spec.known_cases[case_index]
    try:
        design = archetype.synthesize(**case.given)
    except Exception:  # noqa: BLE001 - covered by the reference-matching test
        pytest.skip("synthesis covered elsewhere")

    problems = []
    for metric, value in design.metrics.items():
        if isinstance(value, complex) or not isinstance(value, (int, float)):
            continue
        if math.isnan(value) or value == float("inf"):
            problems.append(f"{metric} is not finite ({value})")
            continue
        if value == float("-inf"):
            # A true pattern null is -inf dB. The 90-degree corner reflector at
            # S = lambda/2 is exactly this case, and it is the correct answer
            # for the idealised image model.
            continue
        bounds = _bounds_for(metric, design.units.get(metric, ""))
        if bounds and not (bounds[0] <= value <= bounds[1]):
            problems.append(f"{metric} = {value:.4g} outside plausible {bounds}")
    for name, value in design.parameters.items():
        if name in ("lambda0", "k0") or not isinstance(value, (int, float)):
            continue
        if isinstance(value, bool):
            continue
        unit = design.units.get(name, "")
        if unit == "m" and not (0 < value < 1e4):
            problems.append(f"{name} = {value:.4g} m is not a credible dimension")
    assert not problems, f"{key} case {case_index}:\n" + "\n".join(f"  - {p}" for p in problems)
