"""The annular ring's TM11 mode in `patch_q`, and the spec's fits to it.

On a thin board the ring radiates as the two edge walls' magnetic ring
currents in free space over the ground; that model is written out here from
the Bessel-function far field, independently of `patch_q`'s current-sheet
integral through the slab, and the two must meet as the board vanishes. The
spec's radiation Q and substrate factor are fits to `patch_q`; they are checked
against fresh solves at designs that were not in the survey.
"""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import jv, jvp, yv, yvp

from otahub.num import patch_q as pq

C = 2.99792458e8


def _ring_k(a, b):
    """Lowest root of J1'(ka) Y1'(kb) = J1'(kb) Y1'(ka), found by its own scan."""
    g = lambda k: jvp(1, k * a) * yvp(1, k * b) - jvp(1, k * b) * yvp(1, k * a)
    ks = np.linspace(0.3, 1.7, 700) * 2 / (a + b)
    v = np.sign([g(k) for k in ks])
    i = int(np.flatnonzero(v[1:] != v[:-1])[0])
    return brentq(g, ks[i], ks[i + 1], xtol=1e-13 / b)


def _ring_currents(k0, k, a, b, eps_r, h):
    """(directivity, radiation Q) of the two edge walls' magnetic ring currents."""
    def R(r):
        return jv(1, k * r) * yvp(1, k * a) - yv(1, k * r) * jvp(1, k * a)

    wa, wb = -a * R(a), b * R(b)

    def pattern(t):
        sa, sb = k0 * a * math.sin(t), k0 * b * math.sin(t)
        e_th = wb * (jv(0, sb) - jv(2, sb)) + wa * (jv(0, sa) - jv(2, sa))
        e_ph = math.cos(t) * (wb * (jv(0, sb) + jv(2, sb)) + wa * (jv(0, sa) + jv(2, sa)))
        return (e_th ** 2 + e_ph ** 2) * math.sin(t)

    power = quad(pattern, 0, math.pi / 2, limit=400)[0]
    energy = quad(lambda r: R(r) ** 2 * r, a, b, limit=400)[0]
    return 4 * (wb + wa) ** 2 / power, 4 * eps_r * energy / (k0 * h * power)


@pytest.mark.parametrize("eps_r,ratio", [(2.2, 1.5), (2.2, 3.0), (10.2, 2.0), (4.4, 1.2)])
def test_patch_q_ring_meets_the_edge_currents_on_thin_board(eps_r, ratio, registry):
    f = 1e10
    k0 = 2 * math.pi * f / C
    h = 1e-4 * C / f
    d = registry["annular_ring_patch"].synthesize(f0=f, eps_r=eps_r, h=h, ratio=ratio)
    a, b = d.get("a_in"), d.get("b_out")
    mode = pq.annulus(a, b)
    d_ring, q_ring = _ring_currents(k0, _ring_k(a, b), a, b, eps_r, h)
    assert pq.radiation_q(*mode, eps_r, h, f) == pytest.approx(q_ring, rel=3e-3)
    assert pq.directivity(*mode, eps_r, h, f) == pytest.approx(d_ring, rel=3e-3)


@pytest.mark.parametrize("g", [
    dict(f0=3.5e9, eps_r=3.0, h=0.0017, ratio=1.8),      # h/lambda0 0.020
    dict(f0=2.0e9, eps_r=8.0, h=0.0015, ratio=2.7),      # h/lambda0 0.010
])
def test_spec_ring_q_and_directivity_match_a_fresh_solve(g, registry):
    d = registry["annular_ring_patch"].synthesize(**g)
    mode = pq.annulus(d.get("a_in"), d.get("b_out"))
    q = pq.radiation_q(*mode, g["eps_r"], g["h"], g["f0"])
    assert d.metrics["radiation_q"] == pytest.approx(q, rel=0.01)
    assert d.metrics["fractional_bandwidth_vswr2"] == pytest.approx(1 / (math.sqrt(2) * q), rel=0.01)
    assert d.metrics["directivity_linear"] == pytest.approx(pq.directivity(*mode, g["eps_r"], g["h"], g["f0"]), rel=0.01)


def test_substrate_factor_is_one_on_thin_board_and_directivity_stops_at_its_domain(registry):
    thin = registry["annular_ring_patch"].synthesize(f0=1e10, eps_r=4.4, h=1e-6, ratio=2.0)
    assert thin.metrics["substrate_directivity_factor"] == pytest.approx(1.0, rel=1e-3)
    assert thin.metrics["directivity_linear"] == pytest.approx(
        thin.metrics["directivity_thin_substrate_linear"], rel=1e-3)
    thick = registry["annular_ring_patch"].synthesize(f0=1e10, eps_r=10.2, h=0.001, ratio=2.0)
    assert math.isnan(thick.metrics["directivity_linear"])         # h sqrt(eps_r)/lambda0 = 0.106
