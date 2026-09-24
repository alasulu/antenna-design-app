"""The Taylor taper for discrete arrays, checked by the pattern it produces.

The taper used to be Taylor's continuous line-source distribution sampled at
the element positions. Over 357 designs that overshot the design sidelobe level
by more than 0.5 dB in 156, by up to 2.3 dB. It is now Villeneuve's discrete
distribution, built by placing the array polynomial's zeros. Every check here
judges the result by something the construction does not use: the array factor
evaluated densely, the separate Dolph-Chebyshev implementation, and Taylor's
continuous distribution re-derived from its textbook formula.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from otahub.arrays.tapers import dolph_chebyshev, taylor_nbar


def _af(w, psi):
    k = np.arange(len(w)) - (len(w) - 1) / 2
    return np.cos(np.outer(psi, k)) @ w                  # symmetric real weights


def _peak_sidelobe_db(w, n_pts=40001):
    psi = np.linspace(0.0, math.pi, n_pts)
    db = 20 * np.log10(np.abs(_af(w, psi)) / w.sum() + 1e-300)
    i = 1
    while i < len(db) - 1 and not (db[i] <= db[i - 1] and db[i] <= db[i + 1]):
        i += 1
    return float(db[i:].max())


def _designs():
    for sll in (-25.0, -30.0, -35.0, -40.0):
        a = math.acosh(10 ** (abs(sll) / 20)) / math.pi
        for n in (10, 16, 20, 33, 64):
            for nbar in range(3, 9):
                if nbar >= 2 * a * a + 0.5:
                    yield sll, n, nbar


@pytest.mark.parametrize("sll,n,nbar", list(_designs()))
def test_the_realised_sidelobes_sit_at_the_design_level(sll, n, nbar):
    got = _peak_sidelobe_db(taylor_nbar(n, sll, nbar))
    assert sll - 0.7 < got < sll + 0.06


def test_the_small_array_the_line_source_failed_on():
    """N = 6 at -25 dB: the sampled line source realised -22.7 dB."""
    assert _peak_sidelobe_db(taylor_nbar(6, -25.0, 5)) < -25.0 + 0.06


@pytest.mark.parametrize("n,nbar", [(20, 5), (33, 4), (64, 6)])
def test_every_zero_past_nbar_is_the_uniform_arrays(n, nbar):
    """The defining property: from the nbar-th on, the pattern nulls exactly
    where a uniform array's do, psi = 2 pi m / N."""
    w = taylor_nbar(n, -30.0, nbar)
    m = np.arange(nbar, (n - 1) // 2 + 1)
    assert np.abs(_af(w, 2 * math.pi * m / n)).max() < 1e-12 * w.sum()


@pytest.mark.parametrize("n", [7, 10, 15, 16])
def test_past_the_last_zero_pair_it_is_dolph_chebyshev(n):
    assert taylor_nbar(n, -30.0, 9) == pytest.approx(dolph_chebyshev(n, -30.0), abs=1e-12)


def test_large_arrays_converge_to_taylors_line_source():
    """Taylor's continuous distribution from the textbook, sampled at the
    elements: g(p) = 1 + 2 sum_m F_m cos(2 pi m p), with F_m the classic
    product over the dilated zeros."""
    n, sll, nbar = 400, -30.0, 5
    a = math.acosh(10 ** (abs(sll) / 20)) / math.pi
    sigma = nbar / math.sqrt(a * a + (nbar - 0.5) ** 2)
    zn = [sigma * math.sqrt(a * a + (i - 0.5) ** 2) for i in range(1, nbar)]
    p = (np.arange(n) - (n - 1) / 2) / n
    g = np.ones(n)
    for m in range(1, nbar):
        num = np.prod([1 - (m / z) ** 2 for z in zn])
        den = np.prod([1 - (m / i) ** 2 for i in range(1, nbar) if i != m])
        g += 2 * ((-1) ** (m + 1)) * num / (2 * den) * np.cos(2 * math.pi * m * p)
    g /= g.max()
    assert np.abs(taylor_nbar(n, sll, nbar) - g).max() < 2e-5


def test_a_thousand_elements_do_not_overflow():
    w = taylor_nbar(1000, -30.0, 5)
    assert np.all(np.isfinite(w)) and np.all(w > 0)
    assert _peak_sidelobe_db(w, 200001) < -30.0 + 0.06
