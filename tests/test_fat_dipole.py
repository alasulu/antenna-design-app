"""Fat dipoles, and therefore wide slots: where the feed starts to decide.

A delta gap on a wire fatter than about 0.01 wavelengths has no converged
answer - its capacitance diverges as the mesh refines - which capped the slot
family's checks at w/L = 0.05. `otahub.num.bor` solves the conductor as an open
tube with a gap of physical length, its mesh graded toward both rims (the rim
current is edge-singular; an ungraded mesh left every length 0.15-0.48% long).
The answer it gives is the finding: past a = 0.005 wavelengths the resonance
depends on the gap. resonant_dipole's length law sits inside that spread from
0.0075 wavelengths on, and 0.2-0.35% above it at 0.005 and below - thin-wire
theory has no rim; its resistance law, fitted to delta-gap solutions, sits at
the small-gap end; and a fat enough tube with a small enough gap has no
resonance at all.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from otahub.num import bor
from otahub.num.strip import tube_profile

DATA = json.loads((Path(__file__).parent / "data" / "fat_dipole_bor.json").read_text())


def _rows(a):
    return [r for r in DATA["rows"] if r["a"] == a and r["L"] is not None]


def _law(registry, a):
    return registry["resonant_dipole"].synthesize(f0=3e8, aw=a * 2.99792458e8 / 3e8)


@pytest.mark.parametrize("a", [0.0075, 0.01, 0.015, 0.02, 0.025])
def test_the_length_law_is_inside_the_feeds_own_spread(registry, a):
    Ls = [r["L"] for r in _rows(a)]
    law = _law(registry, a).metrics["length_over_lambda"]
    assert min(Ls) <= law <= max(Ls), (law, min(Ls), max(Ls))


def test_the_spread_grows_with_the_radius():
    spread = {a: (max(r["L"] for r in _rows(a)) - min(r["L"] for r in _rows(a))) / min(r["L"] for r in _rows(a))
              for a in (0.002, 0.005, 0.01, 0.02)}
    assert spread[0.002] < 0.002 and spread[0.005] < 0.007 and spread[0.01] < 0.02
    assert spread[0.02] > 0.05


def test_the_resistance_law_sits_at_the_small_gap_end(registry):
    for a in (0.005, 0.0075, 0.01):
        Rs = [r["R"] for r in _rows(a)]
        law = _law(registry, a).metrics["input_resistance_ohm"]
        assert law >= max(Rs) - 0.5


@pytest.mark.parametrize("a", [0.002, 0.005])
def test_on_thin_wires_the_length_law_sits_just_above_the_tube(registry, a):
    """The law has no rim; the rim-resolved tube resonates 0.2-0.35% shorter even
    at its smallest gap - inside the law's own declared 0.42%."""
    law = _law(registry, a).metrics["length_over_lambda"]
    top = max(r["L"] for r in _rows(a))
    assert 0.001 < law / top - 1 < 0.0042


def _tube(L, a, gap, seg):
    return tube_profile(L, a, seg, gap)


@pytest.mark.slow
def test_a_fat_tube_with_a_small_gap_never_resonates():
    """a = 0.025 wavelengths, gap a/2: the reactance tops out below zero."""
    a = 0.025
    xs = [0.40 + 0.02 * i for i in range(10)]
    peak = max(bor.solve(_tube(L, a, a / 2, 0.005)).input_impedance.imag for L in xs)
    assert peak < 0


@pytest.mark.slow
def test_recorded_resonance_reproduces_live():
    r = [r for r in DATA["rows"] if r["a"] == 0.01 and r["gap_over_a"] == 2.0][0]
    z = bor.solve(_tube(r["L"], 0.01, 0.02, 0.005)).input_impedance
    assert abs(z.imag) < 0.1
    assert z.real == pytest.approx(r["R"], rel=1e-4)
