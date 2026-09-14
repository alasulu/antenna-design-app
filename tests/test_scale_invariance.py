"""Dimensional audit of every archetype, using no reference data whatsoever.

Maxwell's equations have no preferred length. Take any antenna, multiply every
frequency by S, divide every length by S, and multiply conductivity by S, and
you have a different antenna with *identical* electrical behaviour: the same
directivity, the same beamwidths, the same impedances, the same efficiency.
Lengths shrink by S, frequencies grow by S, areas by S squared, and so on.

That makes a powerful check. Every quantity a spec produces must follow the
power of S its declared unit implies, and nothing else. A formula with a hidden
dimensional constant, a length that is secretly in centimetres, or a unit label
that does not match what the expression computes will all break it.

The value of this test is that it needs no textbook. The `known_cases` harness
checks each archetype against numbers a human chose; this one checks the
formulas against physics, and it caught a radius carried in centimetres while
labelled dimensionless, plus a unit typo, on specs that had been passing their
own cited cases since the first session.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from otahub.core.registry import Registry  # noqa: E402

REG = Registry.load(ROOT / "specs")
ARCHETYPES = sorted(REG, key=lambda a: a.key)

#: Scale factor. Large enough that a dimensional error cannot hide in rounding.
S = 10.0

# Case matters: "S" is siemens, invariant like its reciprocal the ohm, while
# "s" would be seconds and scales as 1/S. Folding their case together silently
# gives conductance the time rule.
_POWER_EXACT = {"S": 0, "s": -1, "S/m": 1, "H": -1, "F": -1}
_POWER = {
    "hz": 1, "m": -1, "m^2": -2, "m^3": -3, "mm^3": -3,
    "rad/m": 1, "1/m": 1,
    # Angles carry no length, so a per-radian quantity scales like its numerator:
    # the Archimedean spiral's growth is metres per radian, the equiangular
    # spiral's is dimensionless per radian.
    "m/rad": -1, "1/rad": 0,
    "-": 0, "dbi": 0, "db": 0, "deg": 0, "ohm": 0, "ohm^2": 0, "rad": 0, "": 0,
}


def power_of(unit: str) -> int | None:
    """Power of S a quantity in `unit` follows. None means 'no rule declared'."""
    u = (unit or "-").strip()
    if u in _POWER_EXACT:
        return _POWER_EXACT[u]
    return _POWER.get(u.lower())


def _inputs(archetype, scale: float) -> dict[str, float]:
    """Supplied values at the given scale, built from each parameter's typical."""
    out: dict[str, float] = {}
    for p in archetype.spec.parameters:
        if p.role not in ("requirement", "material", "assumption"):
            continue
        if p.unit == "Hz":
            lo, hi = archetype.spec.freq_range_hz
            out[p.symbol] = math.sqrt(max(lo, 1e6) * min(hi, 1e12)) * scale
            continue
        try:
            value = float(p.typical)
        except (TypeError, ValueError):
            continue
        k = power_of(p.unit)
        out[p.symbol] = value * (scale ** k) if k else value
    return out


def _pairs(archetype):
    try:
        base = archetype.synthesize(**_inputs(archetype, 1.0))
        scaled = archetype.synthesize(**_inputs(archetype, S))
    except Exception as exc:  # noqa: BLE001
        pytest.skip(f"cannot synthesise from typicals: {type(exc).__name__}: {exc}")
    return base, scaled


@pytest.mark.parametrize("key", [a.key for a in ARCHETYPES])
def test_every_declared_unit_has_a_scaling_rule(key):
    """A unit this test does not recognise is usually a typo ('m2' for 'm^2')."""
    archetype = REG[key]
    declared = {r.units_out for r in archetype.spec.synthesis}
    declared |= {r.units_out for r in archetype.spec.analysis}
    declared |= {p.unit for p in archetype.spec.parameters}
    unknown = sorted(u for u in declared if u and power_of(u) is None)
    assert not unknown, (
        f"{key} declares unit(s) {unknown} with no scaling rule. Either it is a "
        f"typo, or the rule belongs in this test's table."
    )


@pytest.mark.parametrize("key", [a.key for a in ARCHETYPES])
def test_archetype_is_dimensionally_consistent(key):
    """Every quantity must scale as the power of S its unit implies."""
    archetype = REG[key]
    base, scaled = _pairs(archetype)

    problems = []
    for name, v1 in list(base.parameters.items()) + list(base.metrics.items()):
        v2 = scaled.parameters.get(name, scaled.metrics.get(name))
        if v2 is None:
            continue
        if isinstance(v1, complex) or isinstance(v2, complex):
            v1, v2 = abs(v1), abs(v2)
        if isinstance(v1, bool) or isinstance(v2, bool):
            continue
        if not isinstance(v1, (int, float)) or not isinstance(v2, (int, float)):
            continue
        if not (math.isfinite(v1) and math.isfinite(v2)):
            continue
        if abs(v1) < 1e-30:
            continue
        k = power_of(base.units.get(name, ""))
        if k is None:
            continue                      # reported by the unit test above
        want = v1 * S ** k
        err = abs(v2 - want) / max(abs(want), 1e-30)
        if err > 1e-6:
            problems.append(
                f"{name} [{base.units.get(name, '-') or '-'}]: {v1:.6g} at 1x, "
                f"{v2:.6g} at {S:.0f}x, expected {want:.6g} ({err*100:.3g}% off)"
            )
    assert not problems, (
        f"{key} is not dimensionally consistent under a {S:.0f}x scale change:\n"
        + "\n".join(f"  - {p}" for p in problems)
    )


def test_the_audit_actually_covers_the_catalogue():
    """Guard against the check silently degrading into nothing."""
    covered = 0
    for archetype in ARCHETYPES:
        try:
            base = archetype.synthesize(**_inputs(archetype, 1.0))
        except Exception:  # noqa: BLE001
            continue
        covered += sum(
            1
            for name in list(base.parameters) + list(base.metrics)
            if power_of(base.units.get(name, "")) is not None
        )
    assert covered > 700, f"only {covered} quantities audited; the harness has lost coverage"
