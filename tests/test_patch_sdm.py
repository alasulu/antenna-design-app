"""A full-wave arbiter for patches: the spectral-domain method of moments.

`otahub.num.patch_sdm` solves for a rectangular patch's current on its grounded
slab instead of assuming the cavity mode. It is checked here on problems whose
answers it does not share code with - a vanishing dipole, an infinite
microstrip line, and an FDTD ringdown of a whole patch (`patch_fdtd`) - and then
used to settle what a patch's directivity, resonance and Q are.
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
    il = (math.pi * L / 4) * (math.pi * L / 2)      # moment of U0 sqrt(1-t^2) x T0/sqrt(1-t^2)
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
    p = sdm.RectPatch(4.4, 0.0128 * lam, 0.22 * lam, 0.29 * lam, **sdm.BASIS_FULL)
    v = np.array([1.0, -0.4, 0.05, 0.02, -0.01, 0.005, 0.03, -0.01, 0.01])
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
    p = sdm.RectPatch(er, h, d.get("L"), d.get("W"), **sdm.BASIS_FULL)
    m = sdm.mode_metrics(p, sdm.resonance(p, 0.95 * f0))
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


# ------------------------------------------------------------------ against the FDTD

FULL = json.loads((Path(__file__).parent / "data" / "patch_fullwave.json").read_text())


@pytest.mark.parametrize("which", ["fdtd_check", "fdtd_check_fr4"])
def test_fdtd_and_sdm_agree_on_whole_patches(which):
    """The same patch two unrelated ways: FDTD at three cell sizes taken to
    zero, and the spectral MoM with its default basis - thick eps_r 10.2 and FR-4."""
    from otahub.num import patch_fdtd
    ref = FULL[which]
    runs = ref["fdtd_runs"]
    f_inf = patch_fdtd.extrapolate([r["nh"] for r in runs], [r["f_over_f0"] for r in runs])
    q_inf = patch_fdtd.extrapolate([r["nh"] for r in runs], [r["Q"] for r in runs])
    assert ref["sdm"]["f_over_f0"] == pytest.approx(f_inf, rel=1e-3)
    assert ref["sdm"]["Q_total"] == pytest.approx(q_inf, rel=0.01)


@pytest.mark.slow
def test_the_coarsest_fdtd_run_reproduces():
    from otahub.num import patch_fdtd
    ref = FULL["fdtd_check"]
    r = ref["fdtd_runs"][0]
    cell = ref["h_lam"] / r["nh"]
    f, q = patch_fdtd.ringdown(ref["eps_r"], r["nh"], r["nl"], r["nw"], cell)
    assert f / cell == pytest.approx(r["f_over_f0"], rel=2e-3)
    assert q == pytest.approx(r["Q"], rel=0.02)


@pytest.mark.slow
def test_the_default_basis_reproduces_the_recorded_check():
    ref = FULL["fdtd_check"]
    lam = C0 / 1e10
    p = sdm.RectPatch(ref["eps_r"], ref["h_lam"] * lam, ref["L_lam"] * lam, ref["W_lam"] * lam, **sdm.BASIS_FULL)
    f = sdm.resonance(p, 0.96e10)
    assert f / 1e10 == pytest.approx(ref["sdm"]["f_over_f0"], rel=2e-4)


# ------------------------------------------------------------------ what the rectangular spec now carries

def test_the_textbook_patch_resonates_low_and_more_so_on_thick_board(registry):
    lam = C0 / 1e10
    ratios = {}
    for er, hl in ((2.2, 0.003), (2.2, 0.01), (2.2, 0.04), (10.2, 0.02)):
        d = registry["rectangular_patch"].synthesize(f0=1e10, eps_r=er, h=hl * lam)
        ratios[(er, hl)] = d.metrics["full_wave_resonance_hz"] / 1e10
        assert d.metrics["L_full_wave_m"] < d.get("L")
    assert 0.99 < ratios[(2.2, 0.003)] < 0.996
    assert ratios[(2.2, 0.003)] > ratios[(2.2, 0.01)] > ratios[(2.2, 0.04)] > 0.93
    assert ratios[(10.2, 0.02)] < 0.96


def test_surface_wave_efficiency_is_one_on_air_and_jacksons_on_thin_board(registry):
    lam = C0 / 1e10
    air = registry["rectangular_patch"].synthesize(f0=1e10, eps_r=1.0, h=0.02 * lam)
    assert air.metrics["surface_wave_efficiency"] == pytest.approx(1.0, abs=2e-3)
    for er in (2.2, 4.4, 10.2):
        hl = 0.004
        d = registry["rectangular_patch"].synthesize(f0=1e10, eps_r=er, h=hl * lam)
        c1 = 1 - 1 / er + 0.4 / er ** 2
        jackson = 1 / (1 + 0.75 * math.pi * 2 * math.pi * hl / c1 * (1 - 1 / er) ** 3)
        assert d.metrics["surface_wave_efficiency"] == pytest.approx(jackson, rel=5e-3)
        assert d.metrics["fractional_bandwidth_vswr2_with_surface_waves"] > d.metrics["fractional_bandwidth_vswr2"]


@pytest.mark.slow
def test_the_rectangular_spec_against_a_fresh_full_wave_solve(registry):
    """Off the fitting grid: resonance of the textbook length, the length that
    resonates, and the radiation Q and surface-wave efficiency there."""
    g = dict(f0=3.1e9, eps_r=3.66, h=0.0022)
    d = registry["rectangular_patch"].synthesize(**g)
    L, W = d.get("L"), d.get("W")
    fr = sdm.resonance(sdm.RectPatch(g["eps_r"], g["h"], L, W, **sdm.BASIS_FULL), 0.97 * g["f0"])
    Lfw = sdm.resonant_length(g["eps_r"], g["h"], W, g["f0"], L * fr / g["f0"])
    m = sdm.mode_metrics(sdm.RectPatch(g["eps_r"], g["h"], Lfw, W, **sdm.BASIS_FULL), g["f0"])
    assert d.metrics["full_wave_resonance_hz"] == pytest.approx(fr, rel=3e-3)
    assert d.metrics["L_full_wave_m"] == pytest.approx(Lfw, rel=3e-3)
    assert d.metrics["radiation_q"] == pytest.approx(m["q_radiation"], rel=0.015)
    assert d.metrics["surface_wave_efficiency"] == pytest.approx(m["efficiency"], rel=3e-3)
