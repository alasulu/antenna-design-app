"""Horn beamwidths by aperture integration, not flat coefficients.

Every horn here carried a beamwidth of the form K lambda / D with K fixed -
54 and 78 for the sectoral horns, 60/70 for the smooth conical, 66 for the
corrugated - "indicative, near the optimum flare only", while each spec lets
the flare vary. `otahub.num.horn_pattern` integrates the aperture field with its
quadratic phase error. It is checked first against the closed form the E-plane
pattern has in Fresnel integrals, which it does not use.
"""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.special import fresnel

from otahub.num import horn_pattern as hp

C0 = 2.99792458e8


def _fresnel_uniform(D, s, theta):
    a = 8 * math.pi * s / D ** 2
    b = 2 * math.pi * math.sin(theta)
    x0 = b / (2 * a)
    k = math.sqrt(2 * a / math.pi)
    s2, c2 = fresnel(k * (D / 2 - x0))
    s1, c1 = fresnel(k * (-D / 2 - x0))
    return (1 + math.cos(theta)) / 2 * ((c2 - c1) - 1j * (s2 - s1)) / k * np.exp(1j * b * b / (4 * a))


@pytest.mark.parametrize("D,s", [(5.0, 0.25), (12.0, 0.1), (3.0, 0.45)])
def test_quadrature_matches_the_fresnel_closed_form(D, s):
    for th in (0.0, 0.05, 0.2, 0.5):
        assert hp.slit_uniform(D, s, th) == pytest.approx(_fresnel_uniform(D, s, th), rel=1e-10)


def test_the_old_coefficients_are_the_optimum_flare_values():
    """Right where they were derived, and only there."""
    e = hp.hpbw(lambda t: hp.slit_uniform(20.0, 0.25, t), 20.0) * 20.0
    h = hp.hpbw(lambda t: hp.slit_cosine(20.0, 0.375, t), 20.0) * 20.0
    assert e == pytest.approx(54.0, rel=0.005) and h == pytest.approx(78.0, rel=0.005)
    assert hp.hpbw(lambda t: hp.slit_uniform(20.0, 0.4, t), 20.0) * 20.0 > 62.0


def test_a_zero_rim_j0_aperture_is_wider_than_66():
    """The corrugated horn's HE11 field, with no phase error: 75.7 lambda/D."""
    assert hp.hpbw(lambda t: hp.he11(20.0, 0.0, t), 20.0) * 20.0 == pytest.approx(75.8, abs=0.3)


def test_an_over_flared_e_plane_horn_has_no_beamwidth_to_quote(registry):
    d = registry["e_plane_sectoral_horn"].synthesize(f0=1e10, a_wg=0.02286, rho=0.3, flare=1.4)
    assert math.isnan(d.metrics["hpbw_e_deg"])
    ok = registry["e_plane_sectoral_horn"].synthesize(f0=1e10, a_wg=0.02286, rho=0.3, flare=1.3)
    assert math.isfinite(ok.metrics["hpbw_e_deg"])


@pytest.mark.slow
@pytest.mark.parametrize("key,metric,given,D_of,s_of,pattern", [
    ("e_plane_sectoral_horn", "hpbw_e_deg", dict(f0=9e9, a_wg=0.02286, rho=0.45, flare=0.8),
     lambda d, lam: d.get("b1") / lam, lambda d: d.metrics["max_phase_error_wavelengths"], hp.slit_uniform),
    ("h_plane_sectoral_horn", "hpbw_h_deg", dict(f0=11e9, b_wg=0.01016, rho=0.2, flare=1.15),
     lambda d, lam: d.get("a1") / lam, lambda d: d.metrics["max_phase_error_wavelengths"], hp.slit_cosine),
    ("conical_horn", "hpbw_e_deg", dict(f0=12e9, L=0.25, flare=1.1),
     lambda d, lam: d.get("dm") / lam, lambda d: d.get("s_phase"), lambda D, s, t: hp.te11(D, s, t, "E")),
    ("conical_horn", "hpbw_h_deg", dict(f0=12e9, L=0.25, flare=1.1),
     lambda d, lam: d.get("dm") / lam, lambda d: d.get("s_phase"), lambda D, s, t: hp.te11(D, s, t, "H")),
    ("corrugated_conical_horn", "hpbw_deg", dict(f0=20e9, L=0.4, flare=1.2),
     lambda d, lam: d.get("dm") / lam, lambda d: d.get("s_phase"), hp.he11)])
def test_spec_beamwidth_is_the_aperture_integral(registry, key, metric, given, D_of, s_of, pattern):
    d = registry[key].synthesize(**given)
    lam = C0 / given["f0"]
    D, s = D_of(d, lam), s_of(d)
    assert d.metrics[metric] == pytest.approx(hp.hpbw(lambda t: pattern(D, s, t), D), rel=0.006)
