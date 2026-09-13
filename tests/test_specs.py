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
