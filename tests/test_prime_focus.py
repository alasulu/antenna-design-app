"""The prime-focus dish, derived from its feed and checked by aperture integration.

`prime_focus_parabolic` asserted four numbers that all follow from one feed
pattern: 0.82 illumination efficiency, 0.85 spillover, 70 lambda/D and a -24 dB
sidelobe. The first two double-counted spillover - 0.82 is the classic optimum
of the taper-spillover PRODUCT - and cost the budget 0.74 dB. They are now
derived from the edge taper and f/D through the Silver cos^n model, and every
check here comes from `otahub.num.paraboloid`, which integrates the aperture
field directly rather than the feed-angle integral the spec evaluates.
"""
from __future__ import annotations

import math

import pytest
from scipy.integrate import quad

from otahub.num import paraboloid as pb

F0 = 1e10
LAM = 2.99792458e8 / F0


@pytest.fixture(scope="module")
def dish(registry):
    return registry["prime_focus_parabolic"]


def _design(dish, fd, et):
    return dish.synthesize(f0=F0, D=1.0, f_over_D=fd, edge_taper_db=et, eps_rms=0.0, d_blockage=0.0)


@pytest.mark.parametrize("fd,et", [(0.3, -8.0), (0.33, -9.5), (0.4, -11.0), (0.42, -12.3),
                                   (0.55, -6.2), (0.6, -14.0), (0.75, -17.7), (0.95, -4.4)])
def test_taper_and_spillover_are_the_aperture_integrals(dish, fd, et):
    d = _design(dish, fd, et)
    taper, spill = pb.efficiencies(pb.feed_exponent(et, fd), fd)
    assert d.get("eta_ill") == pytest.approx(taper, rel=1e-5)
    assert d.get("eta_spill") == pytest.approx(spill, rel=1e-6)


@pytest.mark.slow
@pytest.mark.parametrize("fd,et,db", [(0.33, -9.5, 0.0), (0.42, -12.3, 0.063), (0.55, -6.2, 0.18),
                                      (0.75, -17.7, 0.0), (0.95, -4.4, 0.11), (0.38, -16.3, 0.137)])
def test_beamwidth_and_peak_sidelobe_are_the_blocked_hankel_transforms(dish, fd, et, db):
    """Off the fit's grid in edge taper, f/D and blockage."""
    d = dish.synthesize(f0=F0, D=1.0, f_over_D=fd, edge_taper_db=et, eps_rms=0.0, d_blockage=db)
    hp, peak = pb.blocked_beam(et, fd, db)
    assert d.metrics["hpbw_deg"] == pytest.approx(hp * LAM, abs=0.012 * LAM)
    assert d.metrics["peak_sidelobe_db"] == pytest.approx(peak, abs=0.2)


def test_with_no_blockage_the_peak_is_the_first_sidelobe_at_the_default_taper(dish):
    """Two transforms: the unblocked `beam` and the blocked one at zero blockage."""
    hp, sll = pb.beam(pb.feed_exponent(-11.0, 0.4), 0.4)
    d = _design(dish, 0.4, -11.0)
    assert d.metrics["peak_sidelobe_db"] == pytest.approx(sll, abs=0.05)
    assert d.metrics["hpbw_deg"] == pytest.approx(hp * LAM, abs=0.005 * LAM)


def _blockage(et, fd, db):
    """By quadrature of the aperture field over rho - not the feed-angle integral the spec uses."""
    n = pb.feed_exponent(et, fd)
    def a(r):
        t = 2 * math.atan(r / (2 * fd))
        return math.cos(t) ** (n / 2) * (1 + math.cos(t)) / 2 * r
    return (quad(a, db / 2, 0.5, epsrel=1e-12)[0] / quad(a, 0, 0.5, epsrel=1e-12)[0]) ** 2


@pytest.mark.parametrize("fd,et,db", [(0.4, -11.0, 0.1), (0.3, -4.0, 0.05), (0.8, -18.0, 0.2), (0.55, -9.0, 0.13)])
def test_blockage_is_weighted_by_the_aperture_field(dish, fd, et, db):
    """The central disc shadows the brightest part of a tapered aperture, so it costs
    more than its area: (1 - (d/D)^2)^2 understated the loss."""
    d = dish.synthesize(f0=F0, D=1.0, f_over_D=fd, edge_taper_db=et, eps_rms=0.0, d_blockage=db)
    assert d.metrics["eta_blockage"] == pytest.approx(_blockage(et, fd, db), abs=2e-5)
    assert d.metrics["eta_blockage"] < (1 - db ** 2) ** 2
    unblocked = _design(dish, fd, et).metrics["aperture_efficiency"]
    assert d.metrics["aperture_efficiency"] == pytest.approx(unblocked * d.metrics["eta_blockage"], rel=1e-12)


def test_a_cos_squared_feed_reaches_the_classic_optimum(dish):
    """Silver's result: 0.829 at a 66-degree rim angle and -10.9 dB edge taper."""
    fd = 1 / (4 * math.tan(math.radians(33.0)))
    et = 20 * math.log10(math.cos(math.radians(66.0)) * (1 + math.cos(math.radians(66.0))) / 2)
    d = _design(dish, fd, et)
    assert d.get("feed_n") == pytest.approx(2.0, rel=1e-9)
    assert d.get("eta_ill") * d.get("eta_spill") == pytest.approx(0.829, abs=1e-3)


def test_the_old_defaults_counted_spillover_twice(dish):
    d = _design(dish, 0.4, -11.0)
    product = d.get("eta_ill") * d.get("eta_spill")
    assert 10 * math.log10(product / (0.82 * 0.85)) > 0.7
    assert product == pytest.approx(0.82, abs=0.015)          # 0.82 was the PRODUCT


@pytest.mark.parametrize("fd", [0.4, 0.5, 0.6])
def test_the_best_edge_taper_is_near_eleven_db(dish, fd):
    tapers = [-8.0 - 0.25 * i for i in range(29)]              # -8 to -15 dB
    best = max(tapers, key=lambda et: _design(dish, fd, et).metrics["aperture_efficiency"])
    assert -11.5 <= best <= -10.0
