"""Patch bandwidths from a solved radiation Q, not a borrowed formula.

The rectangular patch carried the simplified Jackson-Alexopoulos bandwidth,
3.771 (eps_r-1)/eps_r^2 (h/lambda0)(W/L); the circular and triangular patches
borrowed it, and the corner-truncated CP patch took its Q0 - and so its cut -
from it. `otahub.num.patch_q` computes Q = omega W / P from the cavity mode's
stored energy and the space wave of the patch's own current on its grounded
substrate. It is checked before it is used: a vanishing patch reproduces
Jackson's substrate factor c1, and on air its directivity reproduces the
cavity model's edge-current integrals, which share nothing with it.
"""
from __future__ import annotations

import math

import pytest

from otahub.num import patch_q as pq

C0 = 2.99792458e8


@pytest.mark.parametrize("eps_r", [1.0, 2.2, 4.4, 10.2])
def test_a_vanishing_patch_reproduces_jacksons_substrate_factor(eps_r):
    f = 1e10
    lam = C0 / f
    h, L, W = 0.0005 * lam, 0.01 * lam, 0.013 * lam
    q = pq.radiation_q(*pq.rectangle(L, W, 20), eps_r, h, f, n_theta=64, n_phi=128)
    c1 = 1 - 1 / eps_r + 2 / (5 * eps_r ** 2)
    assert q == pytest.approx(3 / 16 * eps_r / c1 * (lam / h) * (L / W), rel=4e-3)


@pytest.mark.parametrize("key,geom,tol", [
    ("circular_patch", lambda d: pq.disc(d.get("a_eff")), 1e-3),
    ("triangular_patch", lambda d: pq.triangle(d.get("a_eff")), 1e-3)])
def test_on_air_the_space_wave_is_the_edge_current_pattern(registry, key, geom, tol):
    """Electric current on the patch versus magnetic current on its edge: the
    same field by equivalence when there is no dielectric."""
    f = 1e10
    h = 0.001 * C0 / f
    d = registry[key].synthesize(f0=f, eps_r=1.0, h=h)
    assert pq.directivity(*geom(d), 1.0, h, f) == pytest.approx(d.metrics["directivity_linear"], rel=tol)


@pytest.mark.slow
@pytest.mark.parametrize("key,geom", [
    # the rectangular patch's Q is full-wave now: see test_patch_sdm
    ("circular_patch", lambda d: pq.disc(d.get("a_eff"))),
    ("triangular_patch", lambda d: pq.triangle(d.get("a_eff")))])
@pytest.mark.parametrize("f0,eps_r,h", [(3.1e9, 1.7, 0.004), (6.2e9, 3.66, 0.000508), (1.2e9, 9.8, 0.00254)])
def test_spec_q_is_the_solved_q_off_the_fitting_grid(registry, key, geom, f0, eps_r, h):
    d = registry[key].synthesize(f0=f0, eps_r=eps_r, h=h)
    q = pq.radiation_q(*geom(d), eps_r, h, f0)
    assert d.metrics["radiation_q"] == pytest.approx(q, rel=0.006)
    assert d.metrics["fractional_bandwidth_vswr2"] == pytest.approx(1 / (math.sqrt(2) * q), rel=0.006)


def test_an_air_patch_has_bandwidth(registry):
    """The simplified formula's (eps_r - 1) factor gave an air patch zero."""
    for key in ("rectangular_patch", "circular_patch", "triangular_patch"):
        d = registry[key].synthesize(f0=2.4e9, eps_r=1.0, h=0.005)
        assert d.metrics["fractional_bandwidth_vswr2"] > 0.02


def test_the_triangle_is_narrower_than_the_rectangle_it_borrowed_from(registry):
    for eps_r in (2.2, 4.4, 10.2):
        tri = registry["triangular_patch"].synthesize(f0=2e9, eps_r=eps_r, h=0.0016)
        old = 3.771 * (eps_r - 1) / eps_r ** 2 * 0.0016 * 2e9 / C0
        assert tri.metrics["fractional_bandwidth_vswr2"] < old / 1.35
