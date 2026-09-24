"""The folded dipole, and bandwidth that is computed rather than assumed.

Bandwidth is the metric this catalogue has most often carried as "indicative",
because the closed forms give Q only through an assumed value. The solver can
do better, and does it twice by routes sharing no algebra: directly, by finding
where VSWR crosses 2 either side of resonance, and through the antenna Q of
Yaghjian and Best, Q = (w0/2R0)|dZ/dw|. The two agreeing to 2% is what makes
either worth putting in a spec, and the way they part company - widening as the
bandwidth grows - is exactly how a narrowband approximation should fail.

The headline result is a correction to a correction. The classic 4:1 step-up
holds, and holds tightly (3.949 to 4.021), but only between structures AT THEIR
OWN RESONANCE. At a fixed frequency it is not 4 and not even constant: it runs
2.2 to 4.9 across a 20% span of length, because the transmission-line mode does
not scale with the radiating one.
"""
from __future__ import annotations

import functools
import math

import numpy as np
import pytest

from otahub.num import mom

pytestmark = pytest.mark.slow   # full-wave solves; skip with -m 'not slow'

NOMINAL_L = 0.48
NSEG = 40          # R moves 0.2% and bandwidth 0.03% from here to 96
LAMBDA = 1.0            # geometry is already in wavelengths
CASES = [(0.01, 1e-4), (0.02, 1e-4), (0.02, 1e-3), (0.04, 3e-4)]


def _folded_z(sep, aw, nseg=NSEG):
    def f(scale):
        m = mom.folded_dipole_wire(NOMINAL_L * scale, sep * scale,
                                   aw * scale, nseg)
        nodes = np.array([m.node_of(n) for n in range(m.n_basis)])
        feed = int(np.argmin(np.abs(nodes[:, 0])
                             + np.abs(nodes[:, 1] - 0.5 * sep * scale)))
        return mom.solve(m, feed).input_impedance
    return f


def _plain_z(aw, nseg=NSEG):
    return lambda s: mom.input_impedance(
        mom.dipole(NOMINAL_L * s, aw * s, nseg))


@functools.lru_cache(maxsize=32)
def _solved(sep, aw):
    f = _folded_z(sep, aw)
    s0 = mom.resonant_scale(f)
    z = f(s0)
    fbw, _ = mom.vswr_bandwidth(f, s0, z.real, step=1e-2)
    return s0, z, fbw, mom.antenna_q(f, s0)


@functools.lru_cache(maxsize=32)
def _solved_plain(aw):
    f = _plain_z(aw)
    s0 = mom.resonant_scale(f)
    z = f(s0)
    fbw, _ = mom.vswr_bandwidth(f, s0, z.real, step=1e-2)
    return s0, z, fbw


def _design(registry, sep, aw):
    return registry["folded_dipole"].synthesize(
        f0=3e8, d_sep=sep * 2.99792458e8 / 3e8, aw=aw * 2.99792458e8 / 3e8, N=2.0)


# ------------------------------------------------ the two bandwidth routes

@pytest.mark.parametrize("sep,aw", CASES)
def test_the_two_bandwidth_routes_agree(sep, aw):
    """The check that licenses a computed bandwidth at all."""
    _, _, fbw, q = _solved(sep, aw)
    assert 0.7071 / q == pytest.approx(fbw, rel=0.03)


def test_a_plain_dipole_gets_wider_as_the_wire_thickens():
    """The oldest rule in wire antennas, and a sanity check on the sweep."""
    widths = [_solved_plain(a)[2] for a in (1e-4, 3e-4, 1e-3)]
    assert widths[0] < widths[1] < widths[2], widths
    assert 0.06 < widths[0] < 0.09, widths[0]


# ------------------------------------------------------- the 4:1, qualified

@pytest.mark.parametrize("sep,aw", CASES)
def test_four_to_one_holds_between_resonances(sep, aw):
    _, z, _, _ = _solved(sep, aw)
    _, pz, _ = _solved_plain(aw)
    assert z.real / pz.real == pytest.approx(4.0, rel=0.02)


def test_four_to_one_does_not_hold_at_a_fixed_frequency():
    """The distinction the spec now draws. If someone later 'simplifies' the
    step-up to a frequency-independent 4, this fails."""
    f, p = _folded_z(0.02, 1e-4), _plain_z(1e-4)
    ratios = [abs(f(s) / p(s)) for s in (0.90, 1.00, 1.10)]
    assert max(ratios) / min(ratios) > 1.5, ratios


# ------------------------------------------------------ the spec's numbers

@pytest.mark.parametrize("sep,aw", CASES)
def test_spec_matches_the_solver(sep, aw, registry):
    s0, z, fbw, _ = _solved(sep, aw)
    d = _design(registry, sep, aw)
    lam = 2.99792458e8 / 3e8
    assert d.get("resonant_length_m") / lam == pytest.approx(
        NOMINAL_L * s0, rel=5e-3)
    assert d.metrics["input_resistance_ohm"] == pytest.approx(z.real, rel=5e-3)
    assert d.metrics["fractional_bandwidth_vswr2"] == pytest.approx(fbw, rel=2e-2)


def test_the_nominal_length_is_longer_than_the_resonant_one(registry):
    """0.48 lambda is a starting point; wider spacing shortens it markedly."""
    lam = 2.99792458e8 / 3e8
    got = []
    for sep, aw in CASES:
        d = _design(registry, sep, aw)
        assert d.get("resonant_length_m") < d.get("L")
        got.append(d.get("resonant_length_m") / lam)
    assert min(got) > 0.44 and max(got) < 0.48, got


# --------------------------------------------- claims that were overstated

def test_it_is_not_twice_a_plain_dipoles_bandwidth(registry):
    """The spec said 'roughly twice'. It is 1.36 to 1.72 times, never 2."""
    for sep, aw in CASES:
        _, _, fbw, _ = _solved(sep, aw)
        _, _, pfbw = _solved_plain(aw)
        assert 1.3 < fbw / pfbw < 1.8, (sep, aw, fbw / pfbw)


def test_the_pattern_is_not_quite_a_plain_dipoles(registry):
    """The spec said 'indistinguishable'. The two conductors are a short
    end-fire pair, so the directivity is systematically up, and more so the
    wider they are spaced - small, but not noise."""
    lam = 2.99792458e8 / 3e8
    close = _design(registry, 0.01, 1e-4).metrics["directivity_linear"]
    wide = _design(registry, 0.04, 3e-4).metrics["directivity_linear"]
    assert close > 1.64093 * 0.995
    assert wide > close, (close, wide)
    assert wide / close > 1.01


def test_resonant_resistance_is_below_the_specs_old_flat_292(registry):
    """4 x 73.079 used the induced-EMF value where a real dipole's resonant
    resistance is about 72, so the old figure ran 1 to 3% high."""
    for sep, aw in CASES:
        r = _design(registry, sep, aw).metrics["input_resistance_ohm"]
        assert 282.0 < r < 292.0, (sep, aw, r)
