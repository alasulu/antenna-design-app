"""The quarter-wave shorted patch's bandwidth, solved rather than halved.

The spec carried half the full patch's bandwidth, "INDICATIVE - this halving has
not been computed". `otahub.num.patch_q` now carries the current a shorting wall
takes down to the ground, radiating through the grounded slab as a vertical
current. The wall's slab factor is checked against a plane-wave boundary-value
solve it does not share code with, and on air the whole shorted patch against
the cavity model's magnetic currents on its three open walls, written out here.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from otahub.num import patch_q as pq

C0 = 2.99792458e8
EPS0 = 8.8541878128e-12
ETA0 = 376.730313668


# ------------------------------------------------------------------ independent references

def _slab_ez_over_ex(theta, eps_r, k0, h, z):
    """E_z(z) / E_x(h) inside a grounded slab under an incident TM plane wave,
    by matching the boundary conditions numerically."""
    w = k0 * C0
    kx = -k0 * math.sin(theta)
    kz0 = k0 * math.cos(theta)
    kz1 = k0 * np.sqrt(eps_r - math.sin(theta) ** 2 + 0j)
    e1, e0 = EPS0 * eps_r, EPS0
    M = np.array([
        [1, -1, 0],
        [np.exp(-1j * kz1 * h), np.exp(1j * kz1 * h), -np.exp(-1j * kz0 * h)],
        [kz1 / (w * e1) * np.exp(-1j * kz1 * h), -kz1 / (w * e1) * np.exp(1j * kz1 * h), -kz0 / (w * e0) * np.exp(-1j * kz0 * h)],
    ], dtype=complex)
    rhs = np.array([0, np.exp(1j * kz0 * h), -kz0 / (w * e0) * np.exp(1j * kz0 * h)], dtype=complex)
    A, B, _ = np.linalg.solve(M, rhs)
    ex_h = kz1 / (w * e1) * (A * np.exp(-1j * kz1 * h) - B * np.exp(1j * kz1 * h))
    hy = A * np.exp(-1j * kz1 * z) + B * np.exp(1j * kz1 * z)
    return (-kx * hy / (w * e1)) / ex_h


def _magnetic_wall_q(walls, energy, k0, omega, nt=96, nph=192):
    """Q from magnetic line currents M = 2 psi (z x n) on the open walls, in free space."""
    pts, wts, mx, my = [], [], [], []
    x, w = np.polynomial.legendre.leggauss(64)
    t = 0.5 * (x + 1)
    for p0, p1, psi, n_out in walls:
        p = np.outer(1 - t, p0) + np.outer(t, p1)
        ps = psi(p[:, 0])
        pts.append(p)
        wts.append(0.5 * math.dist(p0, p1) * w)
        mx.append(-2 * ps * n_out[1])
        my.append(2 * ps * n_out[0])
    P, Wt, Mx, My = np.vstack(pts), np.concatenate(wts), np.concatenate(mx), np.concatenate(my)
    xt, wt = np.polynomial.legendre.leggauss(nt)
    th, wth = 0.25 * math.pi * (xt + 1), 0.25 * math.pi * wt
    ph = 2 * math.pi * np.arange(nph) / nph
    tot = 0.0
    for i, t_ in enumerate(th):
        st, ct = math.sin(t_), math.cos(t_)
        E = np.exp(1j * k0 * st * (np.outer(np.cos(ph), P[:, 0]) + np.outer(np.sin(ph), P[:, 1])))
        fx, fy = E @ (Mx * Wt), E @ (My * Wt)
        mth = ct * (fx * np.cos(ph) + fy * np.sin(ph))
        mph = -fx * np.sin(ph) + fy * np.cos(ph)
        tot += wth[i] * st * float(np.sum(np.abs(mth) ** 2 + np.abs(mph) ** 2)) * 2 * math.pi / nph
    return omega * energy / (k0 ** 2 / (32 * math.pi ** 2 * ETA0) * tot)


# ------------------------------------------------------------------ the arbiter

@pytest.mark.parametrize("eps_r", [1.0, 2.2, 10.2])
@pytest.mark.parametrize("theta", [0.3, 1.0, 1.45])
def test_the_wall_slab_factor_is_the_plane_wave_solution(eps_r, theta):
    k0, h = 2 * math.pi / 0.03, 0.002
    z = np.linspace(0, h, 4001)
    num = np.trapezoid(_slab_ez_over_ex(theta, eps_r, k0, h, z), z)
    assert pq.slab_vertical_factor(theta, eps_r, k0) == pytest.approx(num, rel=1e-6)


def test_on_air_the_wall_factor_is_a_vertical_current_and_its_image():
    k0, h = 2 * math.pi / 0.03, 0.004
    for th in (0.2, 0.8, 1.3):
        c, s = math.cos(th), math.sin(th)
        t = math.tan(k0 * h * c)
        G = 2j * c * t * c / (1j * c * t + c)
        image = -2 * math.sin(k0 * h * c) * np.exp(-1j * k0 * h * c) * s / (k0 * c)
        assert G * pq.slab_vertical_factor(th, 1.0, k0) == pytest.approx(image, rel=1e-12)


@pytest.mark.parametrize("h_lam", [0.002, 0.01])
def test_on_air_the_shorted_patch_is_its_three_open_walls(h_lam):
    f = 1e10
    lam = C0 / f
    k0, om, h, W, L = 2 * math.pi / lam, 2 * math.pi * f, h_lam * lam, 0.6 * lam, 0.25 * lam
    psi = lambda x: np.sin(math.pi * x / (2 * L))
    walls = [((L, 0), (L, W), psi, (1, 0)), ((0, 0), (L, 0), psi, (0, -1)), ((0, W), (L, W), psi, (0, 1))]
    ref = _magnetic_wall_q(walls, EPS0 / (2 * h) * (L * W / 2), k0, om)
    g = pq.shorted_rectangle(L, W)
    assert pq.radiation_q(*g[:4], 1.0, h, f, wall=g[4]) == pytest.approx(ref, rel=1.5e-3)
    # the open edge alone is not the antenna: the side walls radiate a third of the power
    alone = _magnetic_wall_q(walls[:1], EPS0 / (2 * h) * (L * W / 2), k0, om)
    assert alone > 1.25 * ref
    # and without the wall's current the space wave is simply wrong
    assert pq.radiation_q(*g[:4], 1.0, h, f) > 3 * ref


# ------------------------------------------------------------------ the spec

def test_the_full_patch_comparison_uses_the_same_cavity_model(registry):
    """Like with like: the full patch's cavity-current Q on the same board, not
    the rectangular spec's full-wave Q, which runs 1.6-28% lower."""
    for g in (dict(f0=2.4e9, eps_r=4.4, h=0.0016), dict(f0=1e10, eps_r=2.2, h=0.000787)):
        s = registry["quarter_wave_shorted_patch"].synthesize(**g).metrics
        r = registry["rectangular_patch"].synthesize(**g)
        q = pq.radiation_q(*pq.rectangle(r.get("L_textbook") + 2 * r.get("dL"), r.get("W")), g["eps_r"], g["h"], g["f0"])
        assert s["full_patch_radiation_q"] == pytest.approx(q, rel=6e-3)
        assert r.metrics["radiation_q"] < s["full_patch_radiation_q"]


def test_shorting_is_not_halving(registry):
    rel = {}
    for er in (1.0, 2.0, 4.4, 12.0):
        s = registry["quarter_wave_shorted_patch"].synthesize(f0=1e10, eps_r=er, h=0.01 * C0 / 1e10 / math.sqrt(er))
        rel[er] = s.metrics["bandwidth_relative_to_full_patch"]
    assert rel[1.0] > 1.8                       # broader than the full patch on air
    assert rel[2.0] == pytest.approx(1.0, abs=0.1)
    assert 0.6 < rel[4.4] < 0.75
    assert 0.55 < rel[12.0] < 0.6               # approaching, never reaching, a half
    assert rel[1.0] > rel[2.0] > rel[4.4] > rel[12.0]


def test_an_air_shorted_patch_has_a_bandwidth(registry):
    s = registry["quarter_wave_shorted_patch"].synthesize(f0=2.4e9, eps_r=1.0, h=0.005)
    assert s.metrics["fractional_bandwidth_vswr2"] > 0.08


@pytest.mark.slow
@pytest.mark.parametrize("eps_r,h_lam", [(1.4, 0.02), (2.9, 0.008), (6.6, 0.03), (11.5, 0.012)])
def test_the_fit_against_a_fresh_solve(registry, eps_r, h_lam):
    f = 1e10
    h = h_lam * C0 / f
    d = registry["quarter_wave_shorted_patch"].synthesize(f0=f, eps_r=eps_r, h=h)
    g = pq.shorted_rectangle(d.get("L") + d.get("dL"), d.get("W"))
    q = pq.radiation_q(*g[:4], eps_r, h, f, wall=g[4])
    assert d.metrics["radiation_q"] == pytest.approx(q, rel=5e-3)
