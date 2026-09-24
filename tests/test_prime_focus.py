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
@pytest.mark.parametrize("fd,et", [(0.33, -9.5), (0.42, -12.3), (0.55, -6.2), (0.75, -17.7),
                                   (0.95, -4.4), (0.38, -19.3)])
def test_beamwidth_and_first_sidelobe_are_the_hankel_transforms(dish, fd, et):
    """Off the fit's grid in both edge taper and f/D."""
    d = _design(dish, fd, et)
    hp, sll = pb.beam(pb.feed_exponent(et, fd), fd)
    assert d.metrics["hpbw_deg"] == pytest.approx(hp * LAM, abs=0.03 * LAM)
    assert d.metrics["first_sidelobe_db"] == pytest.approx(sll, abs=0.1)


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
