"""The FDTD's far field: a closed near-to-far-field box, and what it settled.

`horn_fdtd.BoxTransform` DFTs the tangential fields on a box round the radiator
and radiates the equivalent currents with their images through the grid's
symmetry walls. It is checked live against sources whose directivity is known
exactly, and a patch in air against the spectral MoM. With it, the FDTD of the
whole horn decided the sectoral horns' narrow dimension, which the aperture
models had left a decibel apart (section 4 of the handover), and gave the PIFA
a directivity where the spec had a placeholder. Those runs take minutes to
hours, so they are kept in tests/data/fdtd_farfield.json and the models they are
compared with are recomputed here.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from otahub.num import horn_fdtd as hf
from tests.test_sectoral_narrow_dimension import _directivities

DATA = json.loads((Path(__file__).parent / "data" / "fdtd_farfield.json").read_text())
LAM = 2.99792458e8 / 10e9


@pytest.mark.parametrize("height", [0.25, 0.6])
def test_a_dipole_over_the_plane_returns_the_array_factor(height):
    d40 = hf.dipole_over_plane(height, cells_per_lambda=40)
    exact = hf.dipole_over_plane_exact(height)
    assert d40 == pytest.approx(exact, rel=0.008)
    # second order: the error falls about fourfold as the cells halve
    err20 = hf.dipole_over_plane(height, cells_per_lambda=20) / exact - 1
    assert 3.0 < err20 / (d40 / exact - 1) < 5.0


def test_a_lone_dipole_on_the_free_grid_is_one_and_a_half():
    assert hf.dipole_over_plane(0.0, cells_per_lambda=40, free=True) == pytest.approx(1.5, rel=0.003)


@pytest.mark.parametrize("run", DATA["air_patch"], ids=lambda r: f"nh{r['nh']}")
def test_an_air_patch_agrees_with_the_spectral_mom(run):
    from otahub.num import patch_sdm as sdm
    cell = run["h"] / run["nh"]
    nl, nw = int(round(run["L"] / 2 / cell)), int(round(run["W"] / 2 / cell))
    p = sdm.RectPatch(1.0, run["h"], 2 * nl * cell, 2 * nw * cell, **sdm.BASIS_FULL)
    mom = sdm.mode_metrics(p, run["f_dft"] * 2.99792458e8)["directivity"]
    assert run["D_fdtd"] == pytest.approx(mom, rel=0.007)


def test_sheet_directivity_live_on_a_coarse_air_patch():
    from otahub.num import patch_fdtd as pf
    from otahub.num import patch_sdm as sdm
    nh, cell = 2, 0.02
    nl, nw = 11, 12
    out = pf.sheet_directivity(1.0, nh, lambda x, y: (x <= nw) & (y <= nl) & (x >= 0) & (y >= 0), nw, 0.0, nl,
                               0.97 * cell / 0.44 * 0.5, (0, nl - 2), box_gap=4, air=14)
    assert out["f_dft"] == pytest.approx(out["f"], rel=0.005)
    assert out["peak"] == pytest.approx(out["directivity"], rel=1e-4)        # a patch peaks at broadside
    p = sdm.RectPatch(1.0, nh * cell, 2 * nl * cell, 2 * nw * cell, **sdm.BASIS_FULL)
    mom = sdm.mode_metrics(p, out["f_dft"] / cell * 2.99792458e8)["directivity"]
    assert out["directivity"] == pytest.approx(mom, rel=0.01)


def _models(h):
    a0, b0, A, B = h["dims"]
    if h["kind"] == "H":
        return _directivities(A, B, True, 0.0, A ** 2 / (8 * h["rho"]))
    return _directivities(A, B, True, B ** 2 / (8 * h["rho"]), 0.0)


@pytest.mark.parametrize("h", [h for h in DATA["horns"] if h["kind"] == "E"], ids=lambda h: f"a{h['dims'][0]:.3f}-N{h['N']}")
def test_the_e_plane_horn_is_the_aperture_in_a_ground_plane(h, registry):
    ap, huy, gnd = _models(h)
    assert h["fdtd_dbi"] == pytest.approx(gnd, abs=0.15)
    if h["dims"][0] < 1.0:                                  # a waveguide-sized a_wg: the aperture-power form reads low
        assert h["fdtd_dbi"] - ap > 0.6 and h["fdtd_dbi"] - huy > 0.5
    a0, _, _, B = h["dims"]
    d = registry["e_plane_sectoral_horn"].synthesize(f0=10e9, a_wg=a0 * LAM, rho=h["rho"] * LAM,
                                                     flare=B / math.sqrt(2 * h["rho"]))
    assert d.metrics["gain_dbi"] == pytest.approx(h["fdtd_dbi"], abs=0.2)


def test_the_h_plane_horn_keeps_the_aperture_power_form():
    default = [h for h in DATA["horns"] if h["kind"] == "H" and h["dims"][1] < 0.4]
    assert sorted(h["N"] for h in default) == [32, 40, 48]    # three grids, each its own staircase of b_wg
    for h in default:
        ap, huy, gnd = _models(h)
        assert -0.35 < h["fdtd_dbi"] - ap < -0.15
        assert gnd - h["fdtd_dbi"] > 0.7 and huy - h["fdtd_dbi"] > 1.7
    for h in (h for h in DATA["horns"] if h["kind"] == "H" and h["dims"][1] >= 0.5):
        ap, _, _ = _models(h)
        assert abs(h["fdtd_dbi"] - ap) < 0.3


def test_the_flange_transform_reproduces_the_te10_aperture_integral(registry):
    """Fed an ideal TE10 field, the aperture-plane transform (M = -z x E imaged through
    the flange) must give what the spec integrates independently, whatever the box."""
    a, b, N = 0.75, 0.35, 40
    xa, yb = int(round(a / 2 * N)), int(round(b / 2 * N))
    spec = registry["open_ended_waveguide"].synthesize(f0=10e9, a_wg=a * LAM, b_wg=b * LAM)
    for m in (4, 10):
        box = hf.BoxTransform(xa + m, yb + m, 5, 1.0 / N, 0.5, k_plane=0)
        X, Y = np.meshgrid(np.arange(box.i0) + 0.5, np.arange(box.j0) + 0.5, indexing="ij")
        ey = np.where((X < xa) & (Y < yb), np.cos(math.pi * X / (2 * xa)), 0.0).astype(complex)
        zero = np.zeros_like(ey)
        box.acc = {"zlo_Ex": zero, "zlo_Ey": ey, "zlo_Hx": zero, "zlo_Hy": zero}
        got = 10 * math.log10(box.directivity(n_theta=80, n_phi=96))
        assert got == pytest.approx(spec.metrics["directivity_dbi"], abs=0.005), m


def test_a_flanged_guide_reads_a_sixth_of_a_decibel_above_the_te10_aperture(registry):
    """Live at 40 cells a wavelength, and the finer grids from the data: the guide's own
    aperture field (edge fields, the evanescent modes the flange excites) adds 0.17 dB
    to the pure-TE10 aperture on three grids."""
    runs = list(DATA["oewg"])
    a, b = 2 * round(0.02286 / LAM / 2 * 40) / 40, 2 * round(0.01016 / LAM / 2 * 40) / 40
    live = hf.sectoral_horn(a, b, a, b, 0.0, cells_per_lambda=40, guide_len=1.5, flange=True, periods=25)
    assert 10 * math.log10(live["directivity"]) == pytest.approx(next(r for r in runs if r["N"] == 40)["fdtd_dbi"], abs=1e-6)
    for r in runs:
        model = registry["open_ended_waveguide"].synthesize(f0=10e9, a_wg=r["a"] * LAM, b_wg=r["b"] * LAM)
        assert 0.15 < r["fdtd_dbi"] - model.metrics["directivity_dbi"] < 0.19, r["N"]


def test_the_pifa_directivity_is_the_fdtds(registry):
    runs = DATA["pifa"]
    for r in runs:
        assert r["peak_theta_deg"] == pytest.approx(90.0, abs=2.0)        # at the horizon
        if r["nh"] != 4:
            continue
        lam = 2.99792458e8 / 1e9
        d = registry["pifa"].synthesize(f0=1e9, h=r["h"] * lam, W=r["W"] * lam, Ws=r["W"] * lam)
        assert d.metrics["directivity_dbi"] == pytest.approx(10 * math.log10(r["D_peak"]), abs=0.05)
        assert d.metrics["directivity_broadside_dbi"] == pytest.approx(10 * math.log10(r["D_broadside"]), abs=0.2)
    for fine in (r for r in runs if r["nh"] == 6):
        coarse = next(r for r in runs if r["nh"] == 4 and r["h"] == fine["h"] and r["W"] == fine["W"])
        assert 10 * math.log10(fine["D_peak"] / coarse["D_peak"]) == pytest.approx(0.0, abs=0.04)


def test_a_partial_short_has_no_pifa_directivity(registry):
    d = registry["pifa"].synthesize(f0=1e9, h=0.0105, W=0.036, Ws=0.018)
    assert math.isnan(d.metrics["directivity_dbi"]) and math.isnan(d.metrics["directivity_broadside_dbi"])


def test_a_far_field_off_the_resonance_is_refused():
    """The ring-down accepts 0.5-1.6 f_guess but the DFT covers 0.8-1.2: a coarse air
    patch resonating 42% above its guess got its far field 15.5% below the resonance,
    under a docstring promising 0.5%."""
    from otahub.num import patch_fdtd as pf
    with pytest.raises(ValueError, match="nearest DFT frequency"):
        pf.sheet_directivity(1.0, 2, lambda x, y: (x <= 12) & (y <= 11) & (x >= 0) & (y >= 0), 12, 0.0, 11,
                             0.013227272727272726, (0, 9), box_gap=4, air=14, npml=8, periods=15)
