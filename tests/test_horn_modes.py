"""Two-mode horns by aperture integration: the diagonal horn and Potter's
dual-mode conical horn.

Both carried their selling points as assertions - the dual-mode horn 0.62
efficiency, +0.85 dB over a smooth horn, 68 lambda/D "equal in E and H", -30 dB
cross-polar; the diagonal horn cross-polar lobes "around -19 dB, not modelled".
`otahub.num.horn_pattern` now transforms both field components of each
aperture, so co- and cross-polar patterns (Ludwig-3) come out of the same
integral. It is checked against the scalar TE11 routine it does not share code
with, and against the diagonal horn's efficiency in Fresnel integrals.
"""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.special import fresnel

from otahub.num import horn_pattern as hp

C0 = 2.99792458e8
LAM = C0 / 1e10


def _fresnel_diag_efficiency(s):
    """E-plane sectoral phase factor times H-plane sectoral efficiency."""
    q = 2 * math.sqrt(s)
    sq, cq = fresnel(q)
    e = (cq ** 2 + sq ** 2) / q ** 2
    r = math.sqrt(8 * s)
    u, v = (1 / r + r) / math.sqrt(2), (1 / r - r) / math.sqrt(2)
    su, cu = fresnel(u)
    sv, cv = fresnel(v)
    return e * ((cu - cv) ** 2 + (su - sv) ** 2) / (8 * s)


# ------------------------------------------------------------------ the arbiter

def test_the_vector_transform_reproduces_the_scalar_te11_routine():
    """dual_mode(0) is TE11 alone; its co-polar E- and H-plane fields must equal
    te11(), which integrates a separately written field."""
    f = hp.dual_mode(0.0)
    for th in (0.0, 0.1, 0.3):
        assert hp.disc_pattern(f, 6.0, 0.3, th, 0.0)[0] == pytest.approx(hp.te11(6.0, 0.3, th, "E"), rel=1e-9)
        assert hp.disc_pattern(f, 6.0, 0.3, th, math.pi / 2)[0] == pytest.approx(hp.te11(6.0, 0.3, th, "H"), rel=1e-9)


def test_te11_alone_has_the_textbook_efficiency_and_cross_pol():
    f = hp.dual_mode(0.0)
    assert hp.disc_efficiency(f, 0.0) == pytest.approx(0.8368, abs=2e-4)
    x = hp.peak_cross(lambda t, phi: hp.disc_pattern(f, 12.0, 0.0, t, phi), 12.0, math.pi / 4)
    assert x == pytest.approx(-18.3, abs=0.2)          # -17.6 at the 3/8 phase error of the spec
    # no cross-polar field at all in the principal planes
    assert abs(hp.disc_pattern(f, 12.0, 0.0, 0.1, 0.0)[1]) < 1e-9 * abs(hp.disc_pattern(f, 12.0, 0.0, 0.0, 0.0)[0])


@pytest.mark.parametrize("s", [0.02, 0.075, 0.25, 0.5])
def test_diagonal_efficiency_matches_the_fresnel_product(s):
    assert hp.diagonal_efficiency(s) == pytest.approx(_fresnel_diag_efficiency(s), rel=1e-6)


def test_diagonal_efficiency_with_no_phase_error_is_8_over_pi_squared():
    assert hp.diagonal_efficiency(0.0) == pytest.approx(8 / math.pi ** 2, rel=1e-6)


def test_diagonal_beams_are_equal_by_symmetry_and_cross_pol_vanishes_between_them():
    for th in (0.05, 0.2, 0.4):
        e = hp.diagonal(5.0, 0.2, th, math.pi / 4)
        h = hp.diagonal(5.0, 0.2, th, 3 * math.pi / 4)
        assert abs(e[0]) == pytest.approx(abs(h[0]), rel=1e-9)
        assert abs(e[1]) < 1e-9 * abs(e[0]) and abs(h[1]) < 1e-9 * abs(h[0])


# ------------------------------------------------------------------ the specs

def test_dual_mode_efficiency_formula_is_the_integral(registry):
    for p in (0.0, 0.05, 0.15, 0.3):
        d = registry["conical_horn_dual_mode"].synthesize(f0=1e10, L=0.3, tm11_fraction=p)
        assert d.metrics["aperture_efficiency"] == pytest.approx(hp.disc_efficiency(hp.dual_mode(p), 0.375), rel=1e-6)


def test_dual_mode_costs_gain_rather_than_adding_it(registry):
    for p in (0.05, 0.1, 0.15, 0.25):
        d = registry["conical_horn_dual_mode"].synthesize(f0=1e10, L=0.3, tm11_fraction=p)
        assert d.metrics["gain_advantage_over_smooth_db"] < 0
    # and with no TM11 it IS the smooth horn
    d0 = registry["conical_horn_dual_mode"].synthesize(f0=1e10, L=0.3, tm11_fraction=0.0)
    assert d0.metrics["gain_advantage_over_smooth_db"] == pytest.approx(0.0, abs=1e-4)


def test_dual_mode_beams_meet_at_the_quoted_fraction(registry):
    d = registry["conical_horn_dual_mode"].synthesize(f0=1e10, L=1.0, tm11_fraction=0.097)
    assert d.metrics["tm11_fraction_for_equal_beams"] == pytest.approx(0.097)
    assert d.metrics["hpbw_e_deg"] == pytest.approx(d.metrics["hpbw_h_deg"], rel=1e-3)
    dm = d.get("dm") / LAM
    assert d.metrics["hpbw_e_deg"] * dm == pytest.approx(78.1, rel=0.01)       # not 68


def test_dual_mode_cross_pol_null_sits_at_the_quoted_fraction(registry):
    reg = registry["conical_horn_dual_mode"]
    xs = {p: reg.synthesize(f0=1e10, L=0.3, tm11_fraction=p).metrics["cross_pol_db"] for p in np.arange(0.1, 0.16, 0.0025)}
    best = min(xs, key=xs.get)
    assert best == pytest.approx(reg.synthesize(f0=1e10, L=0.3, tm11_fraction=0.1).metrics["tm11_fraction_for_min_cross_pol"], abs=0.003)
    assert xs[best] < -34.0


def test_dual_mode_outside_the_surveyed_domain_is_nan(registry):
    reg = registry["conical_horn_dual_mode"]
    assert math.isnan(reg.synthesize(f0=1e10, L=0.3, tm11_fraction=0.45).metrics["cross_pol_db"])
    assert math.isnan(reg.synthesize(f0=1e10, L=0.3, tm11_fraction=0.45).metrics["hpbw_e_deg"])
    assert math.isnan(reg.synthesize(f0=1e10, L=0.05, tm11_fraction=0.1).metrics["cross_pol_db"])     # dm = 2.2 lambda


def test_diagonal_spec_efficiency_falls_with_phase_error(registry):
    reg = registry["diagonal_horn"]
    long_ = reg.synthesize(f0=1e10, a_ap=0.06, R_axial=5.0).metrics
    short = reg.synthesize(f0=1e10, a_ap=0.06, R_axial=0.05).metrics
    assert long_["aperture_efficiency_with_phase_error"] == pytest.approx(8 / math.pi ** 2, rel=2e-3)
    assert short["gain_dbi"] < long_["gain_dbi"] - 2.0
    assert long_["gain_dbi"] == pytest.approx(long_["gain_ideal_dbi"], abs=0.01)


def test_diagonal_cross_pol_is_worse_than_the_old_19_db(registry):
    for R in (0.2, 0.5, 5.0):
        d = registry["diagonal_horn"].synthesize(f0=1e10, a_ap=0.12, R_axial=R)
        assert -16.0 < d.metrics["cross_pol_peak_db"] < -12.0


# ------------------------------------------------------------------ fits vs arbiter, held out

@pytest.mark.slow
@pytest.mark.parametrize("L,p", [(0.2, 0.07), (0.6, 0.2), (2.5, 0.12), (0.15, 0.33)])
def test_dual_mode_spec_against_a_fresh_integration(registry, L, p):
    d = registry["conical_horn_dual_mode"].synthesize(f0=1e10, L=L, tm11_fraction=p)
    D = d.get("dm") / LAM
    f = hp.dual_mode(p)
    he = hp.hpbw(lambda t: hp.disc_pattern(f, D, 0.375, t, 0.0)[0], D)
    hh = hp.hpbw(lambda t: hp.disc_pattern(f, D, 0.375, t, math.pi / 2)[0], D)
    xp = hp.peak_cross(lambda t, phi: hp.disc_pattern(f, D, 0.375, t, phi), D, math.pi / 4, n=400)
    assert d.metrics["hpbw_e_deg"] == pytest.approx(he, rel=2e-3)
    assert d.metrics["hpbw_h_deg"] == pytest.approx(hh, rel=2e-3)
    assert d.metrics["cross_pol_db"] == pytest.approx(xp, abs=0.1)


@pytest.mark.slow
@pytest.mark.parametrize("Dap,s", [(2.5, 0.15), (4.5, 0.33), (9.0, 0.05), (20.0, 0.45)])
def test_diagonal_spec_against_a_fresh_integration(registry, Dap, s):
    a = Dap * LAM
    R = a * a / (8 * s * LAM)
    d = registry["diagonal_horn"].synthesize(f0=1e10, a_ap=a, R_axial=R)
    h = hp.hpbw(lambda t: hp.diagonal(Dap, s, t, math.pi / 4)[0], Dap)
    assert d.metrics["hpbw_deg"] == pytest.approx(h, rel=2e-3)
    if Dap <= 12:
        x = hp.peak_cross(lambda t, phi: hp.diagonal(Dap, s, t, phi), Dap, 0.0, n=400)
        assert d.metrics["cross_pol_peak_db"] == pytest.approx(x, abs=0.02)


def test_all_the_power_in_tm11_is_tm11():
    """dual_mode(1) returned pure TE11 (its TM11 coefficient fell to zero). TM11
    alone has no boresight co-polar field: efficiency zero, as a direct
    integration of its field gives."""
    assert hp.disc_efficiency(hp.dual_mode(1.0), 0.0) < 1e-20
    assert hp.disc_efficiency(hp.dual_mode(0.999999), 0.0) == pytest.approx(8.368e-7, rel=1e-3)
    with pytest.raises(ValueError):
        hp.dual_mode(1.2)


def test_cross_polar_level_is_against_the_co_polar_peak():
    """A long diagonal horn's co-polar maximum leaves boresight (8.6 deg at an
    edge error of 0.8); the level must be taken against that maximum."""
    th = np.linspace(0.0, math.pi / 2 * 0.999, 200)
    co = max(abs(hp.diagonal(6.0, 0.8, t, 0.0)[0]) for t in th)
    x = max(abs(hp.diagonal(6.0, 0.8, t, 0.0)[1]) for t in th)
    got = hp.peak_cross(lambda t, phi: hp.diagonal(6.0, 0.8, t, phi), 6.0, 0.0)
    assert got == pytest.approx(20 * math.log10(x / co), abs=1e-12)
    assert co > 1.2 * abs(hp.diagonal(6.0, 0.8, 0.0, 0.0)[0])
