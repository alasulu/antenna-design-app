"""The single-feed CP patch, checked against a cavity model of its real outline.

`truncated_corner_cp_patch` carried three errors that would each have spoiled a
built antenna, and nothing had checked them because it asserted 3 of the 16
quantities it produced:

- its mode split was f0/(2Q), which is each mode's offset from the CP centre,
  not the separation the two modes need (f0/Q);
- it said to feed on a diagonal - the NEARLY-SQUARE patch's rule. On a diagonal
  of a corner-cut square one of the two modes vanishes, so that feed radiates
  linear polarisation;
- it sized the square to resonate at f0, but cutting the corners lifts one mode
  and leaves the other, so the CP centre landed 0.79% high - outside its own
  axial-ratio band, with 7.9 dB of axial ratio left at f0.

The arbiter is `otahub.num.patch_cavity`: Neumann modes of the truncated square
by finite elements, a probe feed, the broadside polarisation summed over 30
modes. It assumes nothing the spec's rules assume.
"""
from __future__ import annotations

import functools
import math

import numpy as np
import pytest
from scipy.optimize import brentq

from otahub.num import patch_cavity as pc

GIVEN = dict(f0=2.4e9, eps_r=2.2, h=0.0016)
CENTRELINE = (0.2, 0.5)


@pytest.fixture(scope="module")
def cp(registry):
    return registry["truncated_corner_cp_patch"]


def _two_mode(x, d=1.0):
    """Mode amplitudes, modes at x = -d and +d in x = 2Q(f - fc)/fc."""
    return 1 / (1 + 1j * (x + d)), 1 / (1 + 1j * (x - d))


def _ar_db(ea, eb):
    total = abs(ea) ** 2 + abs(eb) ** 2
    circ = 2 * np.imag(ea * np.conj(eb))
    lin = math.sqrt(max(total ** 2 - circ ** 2, 0.0))
    return 20 * math.log10(math.sqrt((total + lin) / max(total - lin, 1e-300)))


# ------------------------------------------------- the two-mode model (quick)

def test_perfect_cp_needs_the_modes_one_q_apart_not_half(cp):
    """Modes at fc(1 -+ 1/(2Q)) give 0 dB; at half that split, not even close."""
    assert _ar_db(*_two_mode(0.0, 1.0)) < 1e-6
    assert _ar_db(*_two_mode(0.0, 0.5)) > 3.0
    d = cp.synthesize(**GIVEN)
    assert d.metrics["mode_split_hz"] == pytest.approx(
        GIVEN["f0"] / d.metrics["quality_factor"], rel=1e-12)


def test_the_axial_ratio_band_is_the_two_mode_models(cp):
    x3 = brentq(lambda x: _ar_db(*_two_mode(x)) - 3.0, 0.0, 1.0)
    d = cp.synthesize(**GIVEN)
    assert d.metrics["axial_ratio_bandwidth_3db"] * d.metrics["quality_factor"] == pytest.approx(
        x3, rel=1e-5)


def test_the_impedance_band_is_twice_one_modes(cp):
    """Two staggered resonances flatten the impedance: |Gamma| = x^2/(4 + x^2)."""
    def swr(x):
        ea, eb = _two_mode(x)
        z = ea + eb                                   # series modes; z(0) = 1 exactly
        g = abs((z - 1) / (z + 1))
        return (1 + g) / (1 - g)
    assert swr(0.0) == pytest.approx(1.0, abs=1e-12)
    xv = brentq(lambda x: swr(x) - 2.0, 0.5, 5.0)
    d = cp.synthesize(**GIVEN)
    q = d.metrics["quality_factor"]
    assert d.metrics["fractional_bandwidth_vswr2"] == pytest.approx(xv / q, rel=1e-6)
    assert d.metrics["fractional_bandwidth_vswr2"] == pytest.approx(
        2 * d.metrics["fractional_bandwidth_vswr2_one_mode"], rel=1e-9)


def test_the_square_is_sized_below_f0(cp):
    d = cp.synthesize(**GIVEN)
    assert d.get("f_sq") < GIVEN["f0"]
    # 1/(2 Q0) to first order; the cut's second-order term, about -6 u^2 of it,
    # is 5% at the full-wave Q0 of 54 (u = 0.096)
    assert GIVEN["f0"] / d.get("f_sq") - 1 == pytest.approx(0.5 / d.get("Q0"), rel=0.07)
    assert GIVEN["f0"] / d.get("f_sq") - 1 < 0.5 / d.get("Q0")


# ------------------------------------------------ the cavity model (slow)

@functools.lru_cache(maxsize=8)
def _cavity(n, m):
    return pc.truncated_square(n, m)


@pytest.mark.slow
def test_the_uncut_square_stays_degenerate_on_the_mesh():
    s = _cavity(80, 0)
    assert s.lam[2] == pytest.approx(s.lam[1], rel=1e-9)
    assert s.lam[1] == pytest.approx(math.pi ** 2, rel=2e-4)


@pytest.mark.slow
def test_the_cut_splits_by_twice_its_area_to_first_order_and_less_beyond():
    small, large = _cavity(240, 14), _cavity(240, 34)
    r_small = small.split / small.cut_over_side ** 2
    r_large = large.split / large.cut_over_side ** 2
    assert 1.95 < r_small < 2.0
    assert r_large < r_small - 0.02


@pytest.mark.slow
def test_a_diagonal_feed_radiates_linear_polarisation():
    s = _cavity(240, 22)
    q, _ = s.circular_q(CENTRELINE)
    assert s.best_axial_ratio(q, CENTRELINE)[0] < 0.01
    assert s.best_axial_ratio(q, (0.5, 0.2))[0] < 0.01          # the other centreline
    assert s.best_axial_ratio(q, (0.2, 0.2))[0] > 100.0         # the cut diagonal
    assert s.best_axial_ratio(q, (0.2, 0.8))[0] > 100.0         # the other diagonal


@pytest.mark.slow
def test_the_cavity_model_reproduces_both_bandwidths():
    s = _cavity(240, 22)
    q, kc = s.circular_q(CENTRELINE)
    f = lambda k: s.axial_ratio_db(k, q, CENTRELINE) - 3.0
    ar_band = (brentq(f, kc, kc + 0.05) - brentq(f, kc - 0.05, kc)) / kc
    assert ar_band * q == pytest.approx(0.347107, rel=2e-3)
    r0 = s.pair_impedance(kc, q, CENTRELINE).real
    def swr(k):
        z = s.pair_impedance(k, q, CENTRELINE)
        g = abs((z - r0) / (z + r0))
        return (1 + g) / (1 - g) - 2.0
    vswr_band = (brentq(swr, kc, kc + 0.1) - brentq(swr, kc - 0.1, kc)) / kc
    assert vswr_band * q == pytest.approx(math.sqrt(2), rel=3e-3)


@pytest.mark.slow
@pytest.mark.parametrize("n,m", [(270, 11), (300, 28), (330, 60)])
def test_the_spec_cuts_what_the_cavity_model_makes_circular(n, m, cp):
    """On meshes the fits never saw: hand the spec the Q a cut is perfectly
    circular at, and it must ask for that cut and put the CP centre there."""
    s = _cavity(n, m)
    q, kc = s.circular_q(CENTRELINE)
    d = cp.synthesize(**GIVEN, Q0=q)
    assert d.metrics["truncation_over_side"] == pytest.approx(s.cut_over_side, rel=1e-3)
    assert GIVEN["f0"] / d.get("f_sq") - 1 == pytest.approx(kc / s.k_square - 1, rel=5e-3)


@pytest.mark.slow
def test_the_design_is_circular_at_f0_and_the_old_one_was_not(cp):
    d = cp.synthesize(**GIVEN)
    q, f_sq = d.get("Q0"), d.get("f_sq")
    # meshes re-matched when Q0 became the full-wave radiation-and-surface-wave Q
    # (54.06; it was the cavity-current radiation Q, 62.3, and before that 59.0):
    new = _cavity(246, 24)                 # c/L 0.097561 against the spec's 0.097555
    assert new.axial_ratio_db(new.k_square * GIVEN["f0"] / f_sq, q, CENTRELINE) < 0.1
    old = _cavity(260, 25)                 # the classic cut, c/L = sqrt(1/(2Q)) = 0.0962
    assert d.get("u_cp") == pytest.approx(25 / 260, rel=1e-3)
    assert old.axial_ratio_db(old.k_square, q, CENTRELINE) > 7.0   # square tuned to f0



def test_the_uncut_squares_mesh_split_is_as_documented():
    """The anti-diagonal mesh keeps the reflection, not the rotation: the uncut
    pair splits by 4.0e-5 at 4 cells, 3.0e-7 at 12 (frequency, as `split`)."""
    assert pc.truncated_square(4, 0, count=5).split == pytest.approx(4.03e-5, rel=0.01)
    assert pc.truncated_square(12, 0, count=5).split == pytest.approx(3.04e-7, rel=0.02)
