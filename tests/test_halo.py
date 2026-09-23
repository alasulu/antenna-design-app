"""The halo: a half-wave dipole bent round until its ends nearly meet.

Verified by continuity rather than by a second formula, the same way the quad
loop was. At zero bend `mom.arc` IS `mom.dipole`, and the dipole is already
validated two ways; following the bend from straight to almost-closed therefore
never leaves validated ground. The gap is modelled as the real air gap it is,
and a lumped load at the same place - a different code path entirely - agrees
on the resistance to four figures.

The spec's own reasoning was sound and two of its numbers were close. What it
got wrong was the azimuth ripple, which is about double what it claimed, and
its geometry, which subtracted the gap twice.
"""
from __future__ import annotations

import functools
import math

import numpy as np
import pytest

from otahub.num import mom

NSEG = 72
LAMBDA = 2.99792458e8 / 1.46e8
CASES = [(0.010, 0.002), (0.015, 0.002), (0.020, 0.003), (0.030, 0.001)]


@functools.lru_cache(maxsize=64)
def _resonance(gap: float, a: float, nseg: int = NSEG):
    lo, hi = 0.38, 0.60
    for _ in range(28):
        mid = 0.5 * (lo + hi)
        m = mom.halo(mid, gap, a, nseg)
        if mom.solve(m, m.n_basis // 2).input_impedance.imag < 0:
            lo = mid
        else:
            hi = mid
    m = mom.halo(lo, gap, a, nseg)
    return lo, mom.solve(m, m.n_basis // 2)


def _horizon(sol):
    _, total, _ = mom.pattern_power(sol, 100, 100)
    phi = np.linspace(0, 2 * math.pi, 721)
    e_th, e_ph = mom.far_field(sol, np.full(721, math.pi / 2), phi)
    u = 4 * math.pi * (np.abs(e_th) ** 2 + np.abs(e_ph) ** 2) / total
    pol = float((np.abs(e_ph) ** 2
                 / (np.abs(e_th) ** 2 + np.abs(e_ph) ** 2)).min())
    return float(u.max()), float(u.min()), pol


def _design(registry, gap, a):
    return registry["halo_loop"].synthesize(
        f0=1.46e8, g=gap * LAMBDA, b=a * LAMBDA)


# -------------------------------------------------- the geometry that was wrong

@pytest.mark.parametrize("gap,a", CASES)
def test_conductor_plus_gap_equals_the_ring(gap, a, registry):
    """The identity the spec broke. It synthesised a 0.47 lambda conductor
    beside a 0.5 lambda ring and a 0.015 lambda gap: 0.47 + 0.015 = 0.485, not
    0.5, because the gap had been subtracted twice."""
    d = _design(registry, gap, a)
    assert d.get("Lc") + d.get("g") == pytest.approx(
        d.get("resonant_circumference_m"), rel=1e-12)


def test_the_nominal_ring_is_only_a_starting_point(registry):
    """0.5 lambda is where you start, not where you end. It sits near the top
    of the resonant range but not outside it: a narrow gap resonates smaller,
    a wide gap on thin wire very slightly larger. Asserting it is always
    smaller would be an overreach, and was - this test caught that."""
    ratios = []
    for gap, a in CASES:
        d = _design(registry, gap, a)
        ratios.append(d.get("resonant_circumference_m") / d.get("circumference_m"))
    assert all(0.90 < r < 1.03 for r in ratios), ratios
    narrow = _design(registry, 0.010, 0.002)
    assert narrow.get("resonant_circumference_m") < narrow.get("circumference_m")


# ------------------------------------------------------ the spec against the solver

@pytest.mark.parametrize("gap,a", CASES)
def test_spec_matches_the_solver(gap, a, registry):
    cres, sol = _resonance(gap, a)
    d = _design(registry, gap, a)
    assert d.get("resonant_circumference_m") / LAMBDA == pytest.approx(cres, rel=5e-3)
    assert d.metrics["input_resistance_ohm"] == pytest.approx(
        sol.input_impedance.real, rel=1e-2)
    hmax, hmin, _ = _horizon(sol)
    assert d.metrics["directivity_linear"] == pytest.approx(hmax, rel=5e-3)
    assert d.metrics["azimuth_ripple_db"] == pytest.approx(
        10 * math.log10(hmax / hmin), rel=2e-2)


# ------------------------------------------------ what a halo is sold as

def test_polarisation_round_the_horizon_is_exactly_horizontal():
    """Not approximately: the horizon field is 100% E_phi at every azimuth.
    This is the whole reason to build a halo, and it is the one claim the spec
    made that is exactly right."""
    _, sol = _resonance(0.015, 0.002)
    _, _, pol = _horizon(sol)
    assert pol == pytest.approx(1.0, abs=1e-6)


def test_the_peak_really_is_in_the_plane_of_the_ring():
    """A halo is a horizon antenna, and it earns that: the maximum round the
    horizon IS the global maximum, not merely a local one."""
    for gap, a in CASES:
        _, sol = _resonance(gap, a)
        hmax, _, _ = _horizon(sol)
        assert mom.directivity(sol, 100, 100) == pytest.approx(hmax, rel=1e-3)


def test_the_azimuth_ripple_is_about_double_what_the_spec_claimed(registry):
    """The substantive error. 1.5 dB was asserted; it is 2.7 to 3.3 dB, so the
    worst-case azimuth gain is 1 to 1.6 dB worse than the spec promised - in
    exactly the number an omnidirectional link budget runs on."""
    for gap, a in CASES:
        d = _design(registry, gap, a)
        assert d.metrics["azimuth_ripple_db"] > 2.5
        assert d.metrics["gain_min_azimuth_dbi"] < -1.0


def test_resistance_is_about_a_fifth_of_a_dipoles(registry):
    """The spec's reasoning, checked: bending the dipole round makes the two
    halves' currents partly cancel. Following the bend continuously shows it."""
    straight = mom.input_impedance(mom.dipole(0.47, 1e-3, 60)).real
    bent = []
    for sub in (0.0, 2.0, 4.0, 6.0):
        m = mom.arc(0.47, sub, 1e-3, 60)
        bent.append(mom.solve(m, m.n_basis // 2).input_impedance.real)
    assert bent[0] == pytest.approx(straight, rel=1e-9)
    assert bent[0] > bent[1] > bent[2] > bent[3]
    assert 0.1 < bent[-1] / straight < 0.3


# --------------------------------------------- the two new solver capabilities

def test_an_arc_at_zero_bend_is_a_dipole():
    """The anchor for every bent-wire result here."""
    for length in (0.2, 0.47, 0.75):
        m = mom.arc(length, 0.0, 1e-3, 60)
        got = mom.solve(m, m.n_basis // 2).input_impedance
        assert got == pytest.approx(
            mom.input_impedance(mom.dipole(length, 1e-3, 60)), rel=1e-12)


def test_a_huge_lumped_load_is_the_same_as_cutting_the_wire():
    """Two unrelated ways of breaking a ring at one point. The resistances must
    agree; the reactances need not, because a cut has no gap capacitance and a
    real gap does."""
    a = 1e-3
    cut = mom.solve(mom.loop(0.5, a, 60), 0, loads={30: 1e9}).input_impedance
    m = mom.halo(0.5, 0.004, a, 60)
    opened = mom.solve(m, m.n_basis // 2).input_impedance
    assert cut.real == pytest.approx(opened.real, rel=2e-3)
    assert cut.imag > opened.imag > 0


def test_a_zero_load_changes_nothing():
    plain = mom.solve(mom.loop(0.5, 1e-3, 40), 0).input_impedance
    loaded = mom.solve(mom.loop(0.5, 1e-3, 40), 0, loads={20: 0.0}).input_impedance
    assert loaded == pytest.approx(plain, rel=1e-12)


def test_a_gap_capacitor_makes_a_halo_smaller_not_larger():
    """Worth pinning because it is the opposite of the intuition that a
    capacitor 'tunes it down'. End-loading a dipole shortens it, so a
    capacitively tuned halo sits below the self-resonant ring size - and a
    0.5 lambda ring cannot be series-resonated by a gap capacitor at all: the
    only zero crossing there is a parallel antiresonance in the kilohms."""
    a = 1e-3
    below = [mom.solve(mom.loop(0.5, a, 60), 0,
                       loads={30: -1j * xc}).input_impedance for xc in (200, 800, 3200)]
    assert all(z.imag > 100.0 for z in below), below
    wide, _ = _resonance(0.030, a)
    narrow, _ = _resonance(0.010, a)
    assert narrow < wide, "a smaller gap means more end capacitance, so a smaller ring"
