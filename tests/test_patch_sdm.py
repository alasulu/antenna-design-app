"""A full-wave arbiter for patches: the spectral-domain method of moments.

`otahub.num.patch_sdm` solves for a rectangular patch's current on its grounded
slab instead of assuming the cavity mode. It is checked here on problems whose
answers it does not share code with, and then used to settle what a patch's
directivity is: the cavity current radiated through the slab (`patch_q`) holds
to about 1%, and the two-slot formula the rectangular specs used does not.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from otahub.num import patch_q as pq
from otahub.num import patch_sdm as sdm

C0 = 2.99792458e8
DATA = json.loads((Path(__file__).parent / "data" / "patch_sdm.json").read_text())


# ------------------------------------------------------------------ against closed forms

@pytest.mark.parametrize("eps_r", [2.2, 4.4, 10.2])
def test_a_vanishing_dipole_radiates_and_launches_surface_waves_as_jackson_says(eps_r):
    """Jackson & Alexopoulos, thin slab: P_space = eta k0^4 h^2 |Il|^2 c1 / (6 pi),
    P_surface = eta k0^5 h^3 |Il|^2 (1 - 1/eps_r)^3 / 8."""
    f = 1e10
    lam = C0 / f
    k0, h, L = 2 * math.pi / lam, 0.001 * lam, 0.002 * lam
    p = sdm.RectPatch(eps_r, h, L, L)
    il = L * L                                      # moment of the one basis function
    c1 = 1 - 1 / eps_r + 0.4 / eps_r ** 2
    eta = 376.730313668
    space = 0.5 * p.R_space(f)[0, 0]
    total = 0.5 * p.Z_trunc(f, 60.0)[0, 0].real
    assert space == pytest.approx(eta * k0 ** 4 * h ** 2 * il ** 2 * c1 / (6 * math.pi), rel=6e-3)
    assert total - space == pytest.approx(eta * k0 ** 5 * h ** 3 * il ** 2 * (1 - 1 / eps_r) ** 3 / 8, rel=6e-3)


def _eeff_kirschning_jansen(f, w, h, er):
    u = w / h
    a = 1 + math.log((u ** 4 + (u / 52) ** 2) / (u ** 4 + 0.432)) / 49 + math.log(1 + (u / 18.1) ** 3) / 18.7
    b = 0.564 * ((er - 0.9) / (er + 3)) ** 0.053
    e0 = (er + 1) / 2 + (er - 1) / 2 * (1 + 10 / u) ** (-a * b)
    fn = f * 1e-9 * h * 1e3
    P1 = 0.27488 + (0.6315 + 0.525 / (1 + 0.0157 * fn) ** 20) * u - 0.065683 * math.exp(-8.7513 * u)
    P2 = 0.33622 * (1 - math.exp(-0.03442 * er))
    P3 = 0.0363 * math.exp(-4.6 * u) * (1 - math.exp(-(fn / 38.7) ** 4.97))
    P4 = 1 + 2.751 * (1 - math.exp(-(er / 15.916) ** 8))
    P = P1 * P2 * ((0.1844 + P3 * P4) * fn) ** 1.5763
    return er - (er - e0) / (1 + P)


@pytest.mark.parametrize("er,u,h", [(2.2, 1.0, 0.5e-3), (4.4, 3.0, 0.8e-3), (10.2, 1.0, 0.635e-3)])
def test_the_same_greens_function_gives_microstrip_its_effective_permittivity(er, u, h):
    """An infinite strip: J_x Maxwell across it, J_y odd; det Z(beta) = 0."""
    from scipy.optimize import brentq
    f = 1e10
    k0 = 2 * math.pi * f / C0
    w = u * h
    p = sdm.RectPatch(er, h, 1.0, w)
    x, wq = np.polynomial.legendre.leggauss(64)
    edges = np.concatenate([[0], k0 * np.geomspace(0.05, 800.0, 120)])
    ky = np.concatenate([0.5 * (b - a) * (x + 1) + a for a, b in zip(edges[:-1], edges[1:])])
    wk = np.concatenate([0.5 * (b - a) * wq for a, b in zip(edges[:-1], edges[1:])])
    Yx = sdm._cheb_ft(ky, w, 0).real
    Yy = sdm._sine_ft(ky, w, 2)

    def det(s):
        beta = s * k0
        kr = np.sqrt(beta ** 2 + ky ** 2)
        ztm, zte = p.znode(kr + 0j, k0)
        ca, sa = beta / kr, ky / kr
        z11 = 2 * np.sum(wk * Yx * Yx * (ztm * ca ** 2 + zte * sa ** 2))
        z22 = 2 * np.sum(wk * np.conj(Yy) * Yy * (ztm * sa ** 2 + zte * ca ** 2))
        z12 = 2 * np.sum(wk * Yx * Yy * (ztm - zte) * ca * sa)
        return (z11 * z22 + z12 * z12).real

    ss = np.linspace(1.0, math.sqrt(er) * 0.999, 60)
    v = [det(s) for s in ss]
    i = next(k for k in range(1, len(ss)) if np.sign(v[k]) != np.sign(v[k - 1]))
    eeff = brentq(det, ss[i - 1], ss[i], xtol=1e-12) ** 2
    assert eeff == pytest.approx(_eeff_kirschning_jansen(f, w, h, er), rel=8e-3)


def test_space_wave_power_from_the_spectrum_equals_the_far_field():
    f = 1e10
    lam = C0 / f
    p = sdm.RectPatch(4.4, 0.0128 * lam, 0.22 * lam, 0.29 * lam, **sdm.BASIS_SMALL)
    v = np.array([1.0, -0.8, 0.05])
    spectral = 0.5 * v @ p.R_space(f) @ v
    far, _ = p.far_field(v, f)
    assert far == pytest.approx(spectral, rel=1e-4)


# ------------------------------------------------------------------ what it settles

@pytest.mark.slow
def test_full_wave_directivity_of_a_thick_high_permittivity_patch(registry):
    f0 = 1e10
    lam = C0 / f0
    er, h = 10.2, 0.0245 * lam
    d = registry["rectangular_patch"].synthesize(f0=f0, eps_r=er, h=h)
    p = sdm.RectPatch(er, h, d.get("L"), d.get("W"), **sdm.BASIS_SMALL)
    m = sdm.mode_metrics(p, sdm.resonance(p, f0))
    cavity = pq.directivity(*pq.rectangle(d.get("L") + 2 * d.get("dL"), d.get("W")), er, h, f0)
    assert m["far_over_spectral"] == pytest.approx(1.0, rel=1e-4)
    assert m["directivity"] == pytest.approx(cavity, rel=0.02)


@pytest.mark.parametrize("case", DATA["cases"], ids=lambda c: f"eps{c['eps_r']}-h{c['h_lam']}")
def test_the_recorded_full_wave_solutions_back_the_cavity_current(case):
    """Largest basis: patch_q's directivity within 1% everywhere, where the
    two-slot formula was 2-13% low."""
    full = case["sets"][-1]["D"]
    assert case["patch_q_D"] == pytest.approx(full, rel=0.01)
    assert case["two_slot_D"] < 0.95 * full


@pytest.mark.parametrize("case", DATA["cases"], ids=lambda c: f"eps{c['eps_r']}-h{c['h_lam']}")
def test_the_rectangular_spec_directivity_is_the_full_wave_one(registry, case):
    lam = C0 / 1e10
    d = registry["rectangular_patch"].synthesize(f0=1e10, eps_r=case["eps_r"], h=case["h_lam"] * lam)
    assert d.metrics["directivity_linear"] == pytest.approx(case["sets"][-1]["D"], rel=0.015)
