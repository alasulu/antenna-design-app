"""Biconical and conical monopole: an angle-convention bug, and solved directivity.

The biconical spec documented theta_h as the HALF angle and then used
cot(theta_h/4) - the formula written for the FULL cone angle - so it reported
the impedance of a cone half as wide: 188 ohm for the classic 47-degree cone
that presents 100. Its sibling specs used the half angle correctly, and two
cross-consistency tests had been written with doubled angles to make them
agree. Those tests now compare like with like (see test_cross_consistency.py).

The directivities were flat constants - the half-wave dipole's 1.6409 for the
biconical, the short monopole's 3.0 for the conical monopole. Both are now
solved at the lowest design frequency, from a wire cage whose directivity
converges to 0.2% in the wire count. Its impedance does not converge, so it is
used below only for an ordering that no plausible convergence could reverse.
"""
from __future__ import annotations

import functools
import math

import numpy as np
import pytest

from otahub.num import mom

pytestmark = pytest.mark.slow   # full-wave solves; skip with -m 'not slow'

ETA = 376.730313412


@functools.lru_cache(maxsize=32)
def _solve(slant, half, wires=16):
    m, f = mom.bicone_cage(slant, half, wires)
    return mom.solve(m, f)


def test_the_solver_sides_with_the_half_angle_impedance():
    """Even unconverged, a 47-degree cage oscillates about 100 ohm, nowhere
    near the 188 the old full-angle formula gave for the same cone."""
    rs = [_solve(s, 47.0).input_impedance.real for s in (0.4, 0.5, 0.6, 0.75)]
    mean = float(np.mean(rs))
    half_form = ETA / math.pi * math.log(1 / math.tan(math.radians(47.0) / 2))
    full_form = ETA / math.pi * math.log(1 / math.tan(math.radians(47.0) / 4))
    assert abs(mean - half_form) < abs(mean - full_form), (mean, half_form, full_form)


def test_biconical_impedance_is_the_half_angle_form(registry):
    d = registry["biconical"].synthesize(f0=1e9, theta_h=math.radians(47.0))
    assert d.metrics["characteristic_impedance_ohm"] == pytest.approx(99.9, abs=0.2)


def test_directivity_converges_in_the_wire_count():
    a = mom.directivity(_solve(0.25, 47.0, 16), 50, 60)
    b = mom.directivity(_solve(0.25, 47.0, 24), 50, 60)
    assert a == pytest.approx(b, rel=5e-3)


@pytest.mark.parametrize("half", [10.0, 30.0, 47.0, 60.0])
def test_spec_directivity_matches_the_cage(half, registry):
    want = mom.directivity(_solve(0.25, half), 50, 60)
    got = registry["biconical"].synthesize(
        f0=1e9, theta_h=math.radians(half)).metrics["directivity_linear"]
    assert got == pytest.approx(want, rel=5e-3)


def test_narrow_cones_are_dipoles_and_wide_ones_are_not(registry):
    """The old note's qualitative claim, now with numbers: a narrow cone
    slightly exceeds the half-wave dipole, a wide one falls well below."""
    def d(half):
        return registry["biconical"].synthesize(
            f0=1e9, theta_h=math.radians(half)).metrics["directivity_linear"]
    assert d(5.0) > 1.6409 > d(30.0) > d(47.0) > d(65.0)
    assert 1.6409 / d(65.0) > 1.2


def test_directivity_moves_across_the_band_in_opposite_directions():
    """Not a constant, and not even monotone in the same sense for every flare:
    a 30-degree cone gains directivity up the band, a 47-degree one loses it
    as its beam lifts off the horizon."""
    narrow = [mom.directivity(_solve(s, 30.0), 50, 60) for s in (0.25, 0.5)]
    wide = [mom.directivity(_solve(s, 47.0), 50, 60) for s in (0.25, 0.5)]
    assert narrow[1] > narrow[0]
    assert wide[1] < wide[0]


def test_the_conical_monopoles_three_was_the_short_monopole_figure(registry):
    """3.0 is right for a SHORT monopole. A quarter-wave cone of 47 degrees
    lands near it by coincidence; narrow and wide cones do not."""
    def d(half):
        return registry["conical_monopole"].synthesize(
            f_low=1e9, cone_half_angle_deg=half).metrics["directivity_linear"]
    assert d(47.0) == pytest.approx(2.9, abs=0.05)
    assert d(5.0) > 3.3 and d(65.0) < 2.75
