"""The dual-mode (Potter) horn's step, solved - and what its flare does.

`conical_horn_dual_mode` took its TM11 share as an assumption because the step
that launches it was not solved. `otahub.num.waveguide_step` solves it by mode
matching, and `otahub.num.bor_fdtd` - an FDTD at azimuthal order 1 that shares
nothing with it - checks it. The same mode matching cascaded over a stepped cone
models a whole horn, and the FDTD checks that too; it shows the first-order
phasing (local propagation constants, no mode conversion) is not good enough
to set the phasing length.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from otahub.num import waveguide_step as ws

DATA = json.loads((Path(__file__).parent / "data" / "potter_fullwave.json").read_text())


def _extrapolate(N, v):
    """f = f_inf + C N^-p through three points, the order fitted."""
    from scipy.optimize import brentq
    N, v = np.asarray(N, float), np.asarray(v, float)
    r = (v[0] - v[1]) / (v[1] - v[2])
    p = brentq(lambda q: (N[0] ** -q - N[1] ** -q) / (N[1] ** -q - N[2] ** -q) - r, 0.3, 4.0)
    C = (v[1] - v[2]) / (N[1] ** -p - N[2] ** -p)
    return v[2] - C * N[2] ** -p


# ------------------------------------------------------------------ the step

def test_mode_matching_conserves_power_and_converges():
    lam = 2.99792458e8 / 1e10
    R, T, M1, M2, P = ws.step(0.5 * 1.1 * lam, 0.5 * 1.4 * lam, 1e10, 12)
    assert (P["reflected"].sum() + P["transmitted"].sum()) / P["incident"] == pytest.approx(1.0, abs=1e-8)
    a = ws.tm11_launch(1.1, 1.4, 12)
    b = ws.tm11_launch(1.1, 1.4, 20)
    assert a[0] == pytest.approx(b[0], rel=2e-3) and a[1] == pytest.approx(b[1], abs=0.1)


def test_the_fdtd_of_the_step_converges_on_mode_matching():
    runs = DATA["step"]["fdtd"]
    a, b = DATA["step"]["radii_cells_at_60"]
    share = _extrapolate([r["N"] for r in runs], [r["share"] for r in runs])
    phase = _extrapolate([r["N"] for r in runs], [r["phase"] for r in runs])
    mm = ws.tm11_launch(2 * a / 60, 2 * b / 60, 20)
    assert share == pytest.approx(mm[0], rel=3e-3)
    assert phase == pytest.approx(mm[1], abs=0.2)


@pytest.mark.slow
def test_the_coarsest_fdtd_step_reproduces():
    from otahub.num import bor_fdtd
    a, b = DATA["step"]["radii_cells_at_60"]
    share, phase = bor_fdtd.step_launch(a, b, 60)
    assert share == pytest.approx(DATA["step"]["fdtd"][0]["share"], rel=1e-5)
    assert phase == pytest.approx(DATA["step"]["fdtd"][0]["phase"], abs=1e-3)


def test_a_narrow_input_guide_puts_the_step_near_tm11_cutoff():
    """A 1.0-wavelength input guide reaches a 0.15 share only 3% above TM11
    cutoff; 1.1 leaves room below f0."""
    from scipy.optimize import brentq
    near = brentq(lambda d: ws.tm11_launch(1.0, d)[0] - 0.15, 1.2199, 1.69)
    clear = brentq(lambda d: ws.tm11_launch(1.1, d)[0] - 0.15, 1.2199, 1.69)
    assert near / 1.2197 - 1 < 0.03 < 0.1 < clear / 1.2197 - 1


# ------------------------------------------------------------------ a whole horn

def _fdtd_staircase(a_in, a_st, ell, a_m, Lz, N):
    """The FDTD's own wall profile, as (radius, length) sections in wavelengths."""
    prof = [(a_in, 10), (a_st, ell)]
    radii = np.round(a_st + (a_m - a_st) * np.arange(Lz + 1) / Lz).astype(int)
    runs, r0, n = [], radii[0], 0
    for r in radii:
        if r == r0:
            n += 1
        else:
            runs.append((r0, n)); r0, n = r, 1
    runs.append((r0, n))
    prof[-1] = (prof[-1][0], prof[-1][1] + runs[0][1])
    prof += runs[1:] + [(a_m, 2)]
    return [(r / N, l / N) for r, l in prof]


@pytest.mark.slow
def test_the_stepped_cone_cascade_is_what_the_fdtd_sees():
    g = DATA["horn"]["geometry_cells_at_60"]
    lam = 1.0
    prof = _fdtd_staircase(*g, 60)
    S, first, last = ws.cascade(prof, 2.99792458e8 / lam, per=18)
    a_m = g[3] / 60
    L = a_m / math.sin(math.atan((g[3] - g[1]) / g[4]))
    share, phase = ws.aperture_share(S[2][:, 0], last, a_m, math.sqrt(L * L - a_m * a_m), 2.99792458e8)
    fine = DATA["horn"]["fdtd_corrected"][-1]
    assert share == pytest.approx(fine["share"], rel=0.05)          # the FDTD still converging upward
    assert phase == pytest.approx(fine["phase"], abs=1.0)


def test_the_first_order_phasing_misses_by_twenty_degrees():
    """Both full-wave routes put the aperture phase near 15 degrees where the
    first-order chain says -7: the flare converts between the modes."""
    h = DATA["horn"]
    fdtd = h["fdtd_corrected"][-1]["phase"]
    assert h["cascade"]["corrected_phase"] == pytest.approx(fdtd, abs=1.0)
    assert abs(h["first_order"]["phase"] - fdtd) > 15
    assert h["plain_cone_corrected"][0]["share"] < 0.01               # a plain cone barely converts
    assert abs(h["cascade"]["slope_deg_per_percent"]) > 1.3 * abs(h["first_order"]["slope_deg_per_percent"])


# ------------------------------------------------------------------ the spec

@pytest.mark.parametrize("d_in,p", [(1.1, 0.15), (1.05, 0.13), (0.95, 0.3), (1.18, 0.08)])
def test_the_spec_step_launches_the_share_it_is_asked_for(registry, d_in, p):
    d = registry["conical_horn_dual_mode"].synthesize(f0=1e10, L=0.3, tm11_fraction=p, d_in_over_lambda=d_in)
    ds = d.metrics["step_diameter_over_lambda"]
    share, phase, _ = ws.tm11_launch(d_in, ds)
    assert share == pytest.approx(p, abs=3e-3)
    assert d.metrics["step_phase_deg"] == pytest.approx(phase, abs=0.4)


def test_an_unreachable_share_has_no_step(registry):
    d = registry["conical_horn_dual_mode"].synthesize(f0=1e10, L=0.3, tm11_fraction=0.35, d_in_over_lambda=1.2)
    assert math.isnan(d.get("d_step"))


def test_the_band_is_narrow_and_shrinks_with_length(registry):
    short = registry["conical_horn_dual_mode"].synthesize(f0=1e10, L=0.3, tm11_fraction=0.13, d_in_over_lambda=1.1).metrics
    long_ = registry["conical_horn_dual_mode"].synthesize(f0=1e10, L=1.0, tm11_fraction=0.13, d_in_over_lambda=1.1).metrics
    assert 0.0 < long_["fractional_bandwidth_xpol_30db"] < short["fractional_bandwidth_xpol_30db"] < 0.06
    above = registry["conical_horn_dual_mode"].synthesize(f0=1e10, L=0.3, tm11_fraction=0.2, d_in_over_lambda=1.1).metrics
    assert above["fractional_bandwidth_xpol_30db"] == pytest.approx(0.0, abs=1e-12)
