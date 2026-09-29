"""Sectoral horns across their unflared, waveguide-sized dimension.

The aperture-power form (Balanis), 4 pi |int E|^2 / (lambda^2 int |E|^2), is
exact for large apertures; the open-ended waveguide showed it reads a small one
2 dB low. A sectoral horn keeps one dimension at waveguide size, so its far field
is integrated here under two aperture models - a Huygens aperture radiating into
all space, and the same aperture in a ground plane (checked against the open-ended
guide's verified value). The two disagree by about a decibel; the FDTD of the whole
horn (tests/test_fdtd_farfield.py) decided between them: the ground plane for the
E-plane horn, whose gain_dbi now carries its difference from the aperture-power
form, and the aperture-power form itself for the H-plane horn.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

K = 2 * math.pi
LAM = 0.02998                    # 10 GHz, metres


def _directivities(ax, by, tx, py_edge, px_edge, nt=200, nph=400, n=56):
    """(aperture-power, Huygens, ground-plane) directivity of an E_y aperture ax x by
    (wavelengths): cosine across x if tx, quadratic edge phases px_edge (x) and py_edge (y)."""
    gx, wx = np.polynomial.legendre.leggauss(n)
    gy, wy = np.polynomial.legendre.leggauss(n)
    x, y = 0.5 * ax * gx, 0.5 * by * gy
    fx = (np.cos(math.pi * x / ax) if tx else np.ones_like(x)) * np.exp(-1j * K * px_edge * (2 * x / ax) ** 2)
    fy = np.exp(-1j * K * py_edge * (2 * y / by) ** 2)
    ax_w, by_w = 0.5 * ax * wx, 0.5 * by * wy
    U0 = abs(np.sum(fx * ax_w) * np.sum(fy * by_w)) ** 2
    d_ap = 4 * math.pi * U0 / (np.sum(np.abs(fx) ** 2 * ax_w) * np.sum(np.abs(fy) ** 2 * by_w))
    out = [d_ap]
    for full in (True, False):
        tmax = math.pi if full else math.pi / 2
        g, w = np.polynomial.legendre.leggauss(nt)
        th, wt = 0.5 * tmax * (g + 1), 0.5 * tmax * w
        ph = 2 * math.pi * (np.arange(nph) + 0.5) / nph
        T, P = np.meshgrid(th, ph, indexing="ij")
        u, v = np.sin(T) * np.cos(P), np.sin(T) * np.sin(P)
        Fx = np.tensordot(np.exp(1j * K * np.multiply.outer(u, x)), fx * ax_w, axes=([-1], [0]))
        Fy = np.tensordot(np.exp(1j * K * np.multiply.outer(v, y)), fy * by_w, axes=([-1], [0]))
        el = ((1 + np.cos(T)) / 2) ** 2 if full else np.sin(P) ** 2 + np.cos(T) ** 2 * np.cos(P) ** 2
        Prad = np.sum(el * np.abs(Fx * Fy) ** 2 * np.sin(T) * wt[:, None]) * 2 * math.pi / nph
        out.append(4 * math.pi * U0 / Prad)
    return [10 * math.log10(v) for v in out]


def test_the_ground_plane_branch_is_the_open_ended_guides():
    _, _, ground = _directivities(0.7625275, 0.3389011, True, 0.0, 0.0)
    assert ground == pytest.approx(6.309, abs=0.01)


def test_the_h_plane_horn_keeps_the_aperture_power_form(registry):
    d = registry["h_plane_sectoral_horn"].synthesize(f0=10e9, b_wg=0.01016, rho=0.3)
    a1, b = d.get("a1") / LAM, 0.01016 / LAM
    ap, huy, gnd = _directivities(a1, b, True, 0.0, (a1 ** 2) / (8 * 0.3 / LAM))
    assert d.metrics["gain_dbi"] == pytest.approx(ap, abs=0.01)
    assert 0.5 < gnd - ap < 0.7 and 1.55 < huy - ap < 1.75
    assert huy - gnd > 0.9                       # the aperture models themselves disagree by a decibel


def test_the_e_plane_horn_gain_is_the_ground_plane_aperture(registry):
    d = registry["e_plane_sectoral_horn"].synthesize(f0=10e9, a_wg=0.02286, rho=0.3)
    aw, b1 = 0.02286 / LAM, d.get("b1") / LAM
    ap, huy, gnd = _directivities(aw, b1, True, (b1 ** 2) / (8 * 0.3 / LAM), 0.0)
    assert d.metrics["gain_aperture_power_dbi"] == pytest.approx(ap, abs=0.01)
    assert 0.1 < huy - ap < 0.35 and 1.0 < gnd - ap < 1.25
    assert d.metrics["gain_dbi"] == pytest.approx(gnd, abs=0.05)   # the fitted correction


@pytest.mark.parametrize("a, rho, flare", [(0.62, 7.0, 0.8), (0.95, 15.0, 1.15), (1.4, 5.0, 0.7), (1.9, 20.0, 1.25)])
def test_the_e_plane_correction_holds_between_its_fitting_points(registry, a, rho, flare):
    d = registry["e_plane_sectoral_horn"].synthesize(f0=10e9, a_wg=a * LAM, rho=rho * LAM, flare=flare)
    b1 = d.get("b1") / LAM
    ap, _, gnd = _directivities(a, b1, True, b1 ** 2 / (8 * rho), 0.0)
    assert d.metrics["gain_aperture_power_dbi"] == pytest.approx(ap, abs=0.01)
    assert d.metrics["gain_dbi"] == pytest.approx(gnd, abs=0.05)


def test_the_e_plane_correction_is_nan_outside_its_range(registry):
    for a, rho in ((0.5, 10.0), (2.5, 10.0), (0.8, 2.0)):        # the last has b1 = 2 wavelengths
        d = registry["e_plane_sectoral_horn"].synthesize(f0=10e9, a_wg=a * LAM, rho=rho * LAM)
        assert math.isnan(d.metrics["gain_dbi"]) and not math.isnan(d.metrics["gain_aperture_power_dbi"])
