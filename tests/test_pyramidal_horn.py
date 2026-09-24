"""The pyramidal horn, made buildable - checked by routes the spec does not use.

The spec used to give the horn a single apex (rho1 = rho2), which fits only a
feed guide with the aperture's own sqrt(3/2) aspect ratio. Real guides are
about 2:1, so its flares met a real guide at two different lengths - 15 mm
apart on WR-90 at 20 dBi. It now solves for the two apex distances so the
flares meet the guide together, with the optimum proportions written in the
AXIAL distances the aperture phase actually depends on.

Every check here is independent of the spec's algebra: the flare geometry is
rebuilt from the output dimensions, the gain comes from direct 2-D aperture
integration, the quartic's root is compared with bisection on the unsquared
condition, and the beamwidths are integrated too.
"""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.optimize import brentq

C0 = 2.99792458e8
WR90 = dict(a_wg=0.02286, b_wg=0.01016)
WR28 = dict(a_wg=0.007112, b_wg=0.003556)
WR187 = dict(a_wg=0.04755, b_wg=0.02215)
DESIGNS = [(10e9, 20.0, WR90), (11e9, 22.6, WR90), (10e9, 15.0, WR90), (10e9, 12.0, WR90),
           (30e9, 25.0, WR28), (30e9, 30.0, WR28), (5e9, 18.0, WR187)]


@pytest.fixture(scope="module")
def horn(registry):
    return registry["pyramidal_horn"]


def _aperture(a1, b1, r1, r2, lam, n=1201):
    k = 2 * math.pi / lam
    x, w = np.polynomial.legendre.leggauss(n)
    xs, ys = 0.5 * a1 * x, 0.5 * b1 * x
    fx = np.cos(math.pi * xs / a1) * np.exp(-1j * k * xs ** 2 / (2 * r2))
    fy = np.exp(-1j * k * ys ** 2 / (2 * r1))
    return (xs, w, fx), (ys, w, fy)


def _gain_db(a1, b1, r1, r2, lam):
    (xs, w, fx), (ys, _, fy) = _aperture(a1, b1, r1, r2, lam)
    num = abs((fx * w).sum() * 0.5 * a1) ** 2 * abs((fy * w).sum() * 0.5 * b1) ** 2
    den = ((abs(fx) ** 2 * w).sum() * 0.5 * a1) * ((abs(fy) ** 2 * w).sum() * 0.5 * b1)
    return 10 * math.log10(4 * math.pi * num / (lam ** 2 * den))


def _hpbw_times_size(coord, w, f, size, lam):
    k = 2 * math.pi / lam
    u = lambda th: abs((f * np.exp(1j * k * coord * math.sin(th)) * w).sum()) ** 2
    th = brentq(lambda t: u(t) - 0.5 * u(0.0), 1e-6, 1.5)
    return 2 * math.degrees(th) * size / lam


def _design(horn, f0, g, guide):
    d = horn.synthesize(f0=f0, G_target=g, **guide)
    return d, {k: float(d.get(k)) for k in ("a1", "b1", "rho_1", "rho_2", "p_len")}


@pytest.mark.parametrize("f0,g,guide", DESIGNS)
def test_both_flares_meet_the_feed_guide_together(horn, f0, g, guide):
    """Rebuild each flare as a straight line from its apex to the aperture edge
    and find where it is as wide as the guide."""
    _, h = _design(horn, f0, g, guide)
    z_e = -h["rho_1"] + h["rho_1"] * guide["b_wg"] / h["b1"]   # aperture at z = 0
    z_h = -h["rho_2"] + h["rho_2"] * guide["a_wg"] / h["a1"]
    assert z_e == pytest.approx(z_h, abs=1e-9)
    assert -z_e == pytest.approx(h["p_len"], rel=1e-9)
    assert h["a1"] > guide["a_wg"] and h["b1"] > guide["b_wg"]


@pytest.mark.parametrize("f0,g,guide", DESIGNS)
def test_the_gain_is_on_target_by_aperture_integration(horn, f0, g, guide):
    d, h = _design(horn, f0, g, guide)
    got = _gain_db(h["a1"], h["b1"], h["rho_1"], h["rho_2"], C0 / f0)
    assert got == pytest.approx(g, abs=1e-3)
    assert d.metrics["gain_dbi"] == pytest.approx(g, abs=1e-6)


@pytest.mark.parametrize("f0,g,guide", DESIGNS)
def test_the_quartic_root_is_the_bisected_one(horn, f0, g, guide):
    d, h = _design(horn, f0, g, guide)
    lam = C0 / f0
    P = d.get("P_opt")

    def mismatch(chi):
        r1, r2 = chi * lam, P * P / chi * lam
        a1, b1 = math.sqrt(3 * lam * r2), math.sqrt(2 * lam * r1)
        return r1 * (1 - guide["b_wg"] / b1) - r2 * (1 - guide["a_wg"] / a1)
    chi = brentq(mismatch, 0.3 * P, 3 * P, xtol=1e-15)
    assert h["rho_1"] / lam == pytest.approx(chi, rel=1e-9)


def test_the_old_single_apex_horn_could_not_be_built(horn):
    """rho1 = rho2 on WR-90 at 20 dBi: the E-plane flare reaches the guide
    15 mm further back than the H-plane one."""
    lam = C0 / 10e9
    rho = 100.0 * lam / (0.51 * 4 * math.pi * math.sqrt(6))
    a1, b1 = math.sqrt(3 * lam * rho), math.sqrt(2 * lam * rho)
    gap = rho * (1 - WR90["b_wg"] / b1) - rho * (1 - WR90["a_wg"] / a1)
    assert gap > 0.014
    assert horn.synthesize(f0=10e9, G_target=20.0, **WR90).metrics[
        "realizability_residual_m"] == pytest.approx(0.0, abs=1e-12)


def test_a_single_apex_is_what_a_guide_of_the_apertures_shape_gets(horn):
    """Continuity with the old model: feed it a guide with a/b = sqrt(3/2) and
    the two apex distances coincide."""
    d = horn.synthesize(f0=10e9, G_target=20.0, a_wg=0.01 * math.sqrt(1.5), b_wg=0.01)
    assert d.get("rho_1") == pytest.approx(d.get("rho_2"), rel=1e-9)


def test_the_textbook_slant_form_undershoots(horn):
    """Balanis's procedure writes the optimum proportions in SLANT lengths and
    sizes with an implied efficiency of 0.5105. Solved exactly and integrated,
    that horn falls short of its target - which is why this spec does not use it."""
    shortfalls = []
    for f0, g in ((10e9, 20.0), (11e9, 22.6), (10e9, 15.0)):
        lam = C0 / f0
        G0 = 10 ** (g / 10)
        K = G0 / (2 * math.pi) * math.sqrt(3 / (2 * math.pi))

        def geom(chi):
            return (K / math.sqrt(chi) * lam, math.sqrt(2 * chi) * lam, chi * lam,
                    G0 ** 2 / (8 * math.pi ** 3 * chi) * lam)

        def mismatch(chi):
            a1, b1, re, rh = geom(chi)
            return ((b1 - WR90["b_wg"]) * math.sqrt((re / b1) ** 2 - 0.25)
                    - (a1 - WR90["a_wg"]) * math.sqrt((rh / a1) ** 2 - 0.25))
        chi1 = G0 / (2 * math.pi * math.sqrt(2 * math.pi))
        a1, b1, re, rh = geom(brentq(mismatch, 0.5 * chi1, 1.5 * chi1))
        r1, r2 = math.sqrt(re ** 2 - b1 ** 2 / 4), math.sqrt(rh ** 2 - a1 ** 2 / 4)
        shortfalls.append(g - _gain_db(a1, b1, r1, r2, lam))
    assert 0.05 < shortfalls[0] < 0.25 and 0.05 < shortfalls[1] < 0.25
    assert shortfalls[2] > 0.6


@pytest.mark.parametrize("f0,g,guide", [d for d in DESIGNS if 15.0 <= d[1] <= 25.0])
def test_the_beamwidth_constants_hold_near_the_optimum(horn, f0, g, guide):
    d, h = _design(horn, f0, g, guide)
    lam = C0 / f0
    (xs, wx, fx), (ys, wy, fy) = _aperture(h["a1"], h["b1"], h["rho_1"], h["rho_2"], lam)
    assert _hpbw_times_size(ys, wy, fy, h["b1"], lam) == pytest.approx(54.1, rel=0.007)
    assert _hpbw_times_size(xs, wx, fx, h["a1"], lam) == pytest.approx(78.1, rel=0.01)
