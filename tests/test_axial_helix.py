"""The axial-mode helix: Kraus measured, and a correction that went the wrong way.

The spec knew Kraus's gain formula overestimates, and told designers to use a
"corrected" formula instead: 8.3 + 10 log10(C^2 N S) + 20 log10(1 + N/10). Its
note called that sub-linear in N. It is super-linear, and it exceeds Kraus for
every N >= 5 - by 2.5 dB at ten turns and 8.5 at thirty - in exactly the regime
where Kraus is already too high. That needed no solver to see, only arithmetic.

The solver says what the right answer is. It models the helix over an infinite
ground plane by image theory, and is cross-checked here two ways that do not
share the pattern integral: power balance, and directivity by reciprocity from
the open-circuit voltage under plane-wave illumination.
"""
from __future__ import annotations

import functools
import math

import numpy as np
import pytest

from otahub.core.constants import ETA0
from otahub.num import mom
from otahub.num.mom import _half_weights

pytestmark = pytest.mark.slow   # full-wave solves; skip with -m 'not slow'

K = 2.0 * math.pi


@functools.lru_cache(maxsize=16)
def _solve(turns, C=1.0, pitch=13.0, spt=12):
    m, feed = mom.helix_over_ground(C, pitch, turns, seg_per_turn=spt)
    return m, feed, mom.solve(m, feed)


def _axial_directivity(m, sol):
    _, total, _ = mom.pattern_power(sol, 70, 70)
    eth, eph = mom.far_field(sol, np.array([1e-6]), np.array([0.0]))
    return 2.0 * 4.0 * math.pi * (abs(eth[0]) ** 2 + abs(eph[0]) ** 2) / total


def _design(registry, N, C=1.0, pitch=13.0):
    return registry["axial_mode_helix"].synthesize(
        f0=1e9, N=N, C_over_lambda=C, pitch_deg=pitch)


# --------------------------------------------------- the correction's direction

def test_the_old_correction_would_have_exceeded_kraus_where_kraus_is_high():
    """Pure arithmetic on the removed formula, kept so the reasoning survives."""
    S = math.tan(math.radians(13))
    for N in (5, 10, 20, 30):
        base = 10 * math.log10(N * S)
        kraus = 11.8 + base
        corrected = 8.3 + base + 20 * math.log10(1 + N / 10)
        assert corrected > kraus - 0.05, N
    assert 8.3 + 20 * math.log10(4) - 11.8 > 8.0     # N = 30: +8.5 dB


def test_the_wrong_direction_formula_is_gone(registry):
    metrics = {r.metric for r in registry["axial_mode_helix"].spec.analysis}
    assert "gain_corrected_dbi" not in metrics
    assert {"gain_dbi", "gain_kraus_dbi", "hpbw_deg", "hpbw_kraus_deg"} <= metrics


# ---------------------------------------------------------- the image model

def test_the_image_structure_is_mirror_symmetric():
    """What makes 'twice the full structure's directivity' legitimate."""
    m, _, sol = _solve(6)
    th = np.linspace(0.1, 1.4, 9)
    ph = np.full_like(th, 0.7)
    a1, b1 = mom.far_field(sol, th, ph)
    a2, b2 = mom.far_field(sol, math.pi - th, ph)
    u1 = np.abs(a1) ** 2 + np.abs(b1) ** 2
    u2 = np.abs(a2) ** 2 + np.abs(b2) ** 2
    assert np.max(np.abs(u1 - u2)) / u1.max() < 1e-5


def test_power_balance_holds_on_the_helix():
    _, _, sol = _solve(8)
    assert mom.radiated_power(sol, 70, 70) == pytest.approx(sol.circuit_power, rel=1e-3)


def test_receive_mode_reciprocity_gives_the_same_directivity():
    """A second route to D that uses the terminal voltage and input resistance
    rather than the radiation integral."""
    m, feed, sol = _solve(8)
    tx = _axial_directivity(m, sol)
    Zi = np.linalg.inv(mom.impedance_matrix(m))
    z_in = 1.0 / Zi[feed, feed]
    x, w = np.polynomial.legendre.leggauss(10)
    best = 0.0
    for pol in (np.array([1, 1j, 0]), np.array([1, -1j, 0])):
        pol = pol / math.sqrt(2)
        V = np.zeros(m.n_basis, dtype=complex)
        for n in range(m.n_basis):
            for p, (al, be), _ in _half_weights(m, n):
                A, u, L = m.seg_a[p], m.seg_t[p], m.seg_len[p]
                s = 0.5 * L * (x + 1)
                pts = A[None, :] + s[:, None] * u[None, :]
                V[n] += 0.5 * L * ((al + be * s) * (pol @ u)
                                   * np.exp(1j * K * pts[:, 2]) * w).sum()
        v_oc = (Zi @ V)[feed] * z_in
        ae = abs(v_oc) ** 2 / (8 * z_in.real / (2 * ETA0))
        best = max(best, 2 * 4 * math.pi * ae)
    assert 10 * math.log10(best / tx) == pytest.approx(0.0, abs=0.1)


# ------------------------------------------------------ the spec's numbers

@pytest.mark.parametrize("N,C,pitch", [(5, 1.0, 13.0), (8, 0.9, 12.0), (10, 1.0, 13.0)])
def test_spec_gain_matches_the_solver(N, C, pitch, registry):
    """At a coarser mesh than the fit used, so allow its ~0.15 dB offset."""
    m, feed = mom.helix_over_ground(C, pitch, N, seg_per_turn=16)
    want = 10 * math.log10(_axial_directivity(m, mom.solve(m, feed)))
    assert _design(registry, N, C, pitch).metrics["gain_dbi"] == pytest.approx(want, abs=0.4)


def test_kraus_grows_more_optimistic_with_length(registry):
    """Roughly right only for a short, small helix; several dB high for a long
    one - the progressive divergence the literature reports, with a size on it.
    The error grows with circumference as well as with turns."""
    corner = _design(registry, 3, C=0.9).metrics          # the one place it holds
    assert abs(corner["gain_kraus_dbi"] - corner["gain_dbi"]) < 1.0
    for N in (8, 10, 15, 20):
        for C in (0.9, 1.0, 1.1):
            m = _design(registry, N, C).metrics
            assert 1.6 < m["gain_kraus_dbi"] - m["gain_dbi"] < 5.6, (N, C)
    gaps = [(_design(registry, N).metrics["gain_kraus_dbi"]
             - _design(registry, N).metrics["gain_dbi"]) for N in (3, 6, 10, 20)]
    assert gaps == sorted(gaps), gaps


def test_steeper_pitch_lowers_the_real_gain_and_raises_kraus(registry):
    """The two disagree on the SIGN of the pitch dependence, not just its size."""
    lo, hi = _design(registry, 10, pitch=12.0).metrics, _design(registry, 10, pitch=14.0).metrics
    assert hi["gain_dbi"] < lo["gain_dbi"] - 0.5
    assert hi["gain_kraus_dbi"] > lo["gain_kraus_dbi"]


def test_kraus_beamwidth_is_too_narrow_by_the_same_token(registry):
    m = _design(registry, 10).metrics
    assert m["hpbw_deg"] > m["hpbw_kraus_deg"] * 1.2


def test_the_band_edge_is_in_the_table(registry):
    """The cliff at C/lambda = 1.2 on long helices, which the fit leaves out."""
    tab = registry["axial_mode_helix"].spec.tables["mom_gain"]
    cols = tab["columns"]
    rows = {(r[0], r[1], r[2]): r[cols.index("gain_dbi")] for r in tab["rows"]}
    assert rows[(20, 1.2, 12.0)] < rows[(20, 1.1, 12.0)] - 2.0
