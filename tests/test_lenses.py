"""Collimating lenses, derived from the feed through a ray-traced lens.

The hyperbolic dielectric lens asserted eta_ill = 0.80 and the metal-plate lens
0.75, each silently including spillover, and both quoted 70 lambda/D. A lens
maps feed angle to aperture radius through its own face, and that mapping
decides the illumination: a hyperbola TAPERS the aperture on top of the feed,
an ellipse (n < 1) BRIGHTENS the rim. Every check here comes from
`otahub.num.lens`, which traces rays with Snell's law at the actual surface.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from otahub.num import lens

REXOLITE, SILICON, PTFE = math.sqrt(2.53), math.sqrt(11.7), math.sqrt(2.1)
PLATES_06 = math.sqrt(1 - (1 / 1.2) ** 2)
PLATES_08 = math.sqrt(1 - (1 / 1.6) ** 2)


def _nf(ft, trim_deg):
    return 2 * ft / (20 * math.log10(math.cos(math.radians(trim_deg))))


def _hyper(registry, eps, trim, ft):
    return registry["hyperbolic_dielectric_lens"].synthesize(f0=1e11, eps_r=eps, F=0.05,
                                                              theta_rim_deg=trim, feed_taper_db=ft)


def _metal(registry, s_over_l, trim, ft):
    return registry["metal_plate_lens"].synthesize(f0=1e10, F=0.3, plate_spacing_over_lambda=s_over_l,
                                                    theta_rim_deg=trim, feed_taper_db=ft)


@pytest.mark.parametrize("n", [PTFE, REXOLITE, SILICON, PLATES_06, PLATES_08])
def test_traced_rays_leave_parallel_where_the_formula_puts_them(n):
    for t in np.linspace(0.02, 0.85 * lens.rim_limit(n), 12):
        direction, rho = lens.trace(n, t)
        assert abs(direction[0]) < 1e-8
        assert rho == pytest.approx(abs(n - 1) * math.sin(t) / abs(n * math.cos(t) - 1), rel=1e-12)


@pytest.mark.parametrize("eps,trim,ft", [(2.53, 32.0, -7.0), (11.7, 55.0, -9.5), (2.1, 20.0, -12.0)])
def test_dielectric_lens_efficiencies_are_the_traced_ones(registry, eps, trim, ft):
    d = _hyper(registry, eps, trim, ft)
    taper, spill = lens.efficiencies(math.sqrt(eps), math.radians(trim), _nf(ft, trim))
    assert d.get("eta_ill") == pytest.approx(taper, rel=2e-4)
    assert d.get("eta_spill") == pytest.approx(spill, rel=2e-4)


@pytest.mark.parametrize("s_over_l,trim,ft", [(0.6, 38.0, -11.0), (0.8, 25.0, -8.5), (0.55, 20.0, -15.0)])
def test_metal_plate_lens_efficiencies_are_the_traced_ones(registry, s_over_l, trim, ft):
    d = _metal(registry, s_over_l, trim, ft)
    n = math.sqrt(1 - (1 / (2 * s_over_l)) ** 2)
    taper, spill = lens.efficiencies(n, math.radians(trim), _nf(ft, trim))
    assert d.get("eta_ill") == pytest.approx(taper, rel=2e-4)
    assert d.get("eta_spill") == pytest.approx(spill, rel=2e-4)


@pytest.mark.slow
@pytest.mark.parametrize("kind,a,trim,ft", [("hyper", 2.53, 32.0, -7.0), ("hyper", 11.7, 55.0, -9.5),
                                            ("metal", 0.6, 38.0, -11.0), ("metal", 0.8, 25.0, -8.5)])
def test_beamwidth_and_edge_illumination_are_the_traced_ones(registry, kind, a, trim, ft):
    d = _hyper(registry, a, trim, ft) if kind == "hyper" else _metal(registry, a, trim, ft)
    n = math.sqrt(a) if kind == "hyper" else math.sqrt(1 - (1 / (2 * a)) ** 2)
    nf = _nf(ft, trim)
    lam = 2.99792458e8 / (1e11 if kind == "hyper" else 1e10)
    hp = lens.hpbw(n, math.radians(trim), nf)
    assert d.metrics["hpbw_deg"] == pytest.approx(hp * lam / d.get("D_ap"), abs=0.3 * lam / d.get("D_ap"))
    _, _, _, A = lens.aperture(n, math.radians(trim), nf)
    assert d.metrics["aperture_edge_illumination_db"] == pytest.approx(20 * math.log10(A[-1] / A[0]), abs=0.02)


def test_a_hyperbola_tapers_and_an_ellipse_brightens(registry):
    h = _hyper(registry, 2.53, 40.0, -8.0)
    m = _metal(registry, 0.6, 30.0, -12.0)
    assert h.metrics["aperture_edge_illumination_db"] < -8.0 - 10.0
    assert m.metrics["aperture_edge_illumination_db"] > -12.0 + 3.0


def test_the_flat_efficiencies_were_wrong_both_ways(registry):
    """0.80 flattered the default Rexolite lens; 0.75 undersold the metal plates."""
    h = _hyper(registry, 2.53, 40.0, -8.0)
    m = _metal(registry, 0.6, 30.0, -12.0)
    assert h.get("eta_ill") * h.get("eta_spill") < 0.55
    assert m.get("eta_ill") * m.get("eta_spill") > 0.85


def test_an_elliptical_face_folds_back_past_acos_n():
    """The limit the metal-plate spec never stated: beyond acos(n) the aperture
    radius shrinks again, so rays would cross."""
    lim = lens.rim_limit(PLATES_06)
    assert lens.trace(PLATES_06, lim + 0.05)[1] < lens.trace(PLATES_06, lim)[1]
    assert lens.trace(PLATES_06, lim - 0.05)[1] < lens.trace(PLATES_06, lim)[1]
