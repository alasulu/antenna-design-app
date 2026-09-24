"""The offset dish, derived from its feed and checked by 2-D aperture integration.

`offset_parabolic` carried the prime-focus defaults - 0.82 x 0.85, with
spillover counted twice - and a +0.09 dB "gain advantage" that counted only
the blockage it removes. Its efficiencies and beamwidth are now derived from
the feed taper and geometry, and every check here comes from
`otahub.num.paraboloid`, which integrates the offset aperture directly.
"""
from __future__ import annotations

import math

import pytest
from scipy.integrate import quad

from otahub.num import paraboloid as pb

F0 = 12e9
LAM = 2.99792458e8 / F0


@pytest.fixture(scope="module")
def dish(registry):
    return registry["offset_parabolic"]


def _design(dish, fd, h, ft):
    return dish.synthesize(f0=F0, D=1.0, f_over_D_parent=fd, h0_over_D=h, feed_taper_db=ft, eps_rms=0.0)


def _n(ft, ts):
    return 2 * ft / (20 * math.log10(math.cos(ts)))


@pytest.mark.parametrize("fd,h", [(0.5, 0.6), (0.3, 0.55), (0.8, 1.0), (0.4, 0.75), (1.0, 1.2)])
def test_the_rim_is_a_circular_cone_about_the_bisector(dish, fd, h):
    """What makes the spillover exact - checked ray by ray, not assumed."""
    angles = pb.offset_rim_angles(fd, h)
    assert angles.max() - angles.min() < 1e-12
    assert _design(dish, fd, h, -10.0).get("theta_star") == pytest.approx(angles.mean(), rel=1e-12)


@pytest.mark.parametrize("fd", [0.3, 0.45, 0.7])
def test_with_no_offset_it_is_the_prime_focus_dish(fd):
    """Two integrators written differently must agree where the geometries coincide."""
    n = pb.feed_exponent(-11.0, fd)
    assert pb.offset_taper(fd, 0.0, n) == pytest.approx(pb.efficiencies(n, fd)[0], rel=1e-9)


@pytest.mark.parametrize("fd,h,ft", [(0.45, 0.65, -9.3), (0.33, 0.9, -13.7), (0.75, 1.1, -5.5),
                                     (0.9, 0.58, -17.2), (0.55, 0.75, -11.0)])
def test_taper_spillover_and_beamwidth_are_the_aperture_integrals(dish, fd, h, ft):
    """Off the fit's grid in all three variables."""
    d = _design(dish, fd, h, ft)
    _, ts = pb.offset_geometry(fd, h)
    n = _n(ft, ts)
    assert d.get("eta_ill") == pytest.approx(pb.offset_taper(fd, h, n), abs=0.003)
    g = lambda t: 2 * (n + 1) * math.cos(t) ** n * math.sin(t)
    assert d.get("eta_spill") == pytest.approx(quad(g, 0, ts)[0] / quad(g, 0, math.pi / 2)[0], rel=1e-6)
    w1, w2 = pb.offset_beamwidths(fd, h, n)
    assert abs(w1 / w2 - 1) < 0.02
    assert d.metrics["hpbw_deg"] == pytest.approx(0.5 * (w1 + w2) * LAM, abs=0.1 * LAM)


def test_the_symmetric_formula_would_have_overstated_the_taper(dish):
    """The tilted spreading loss is real: at a short parent focal length the
    symmetric formula with theta* is 15-20% optimistic."""
    fd, h, ft = 0.3, 0.6, -10.0
    _, ts = pb.offset_geometry(fd, h)
    n = _n(ft, ts)
    i = quad(lambda t: math.sqrt(2 * (n + 1) * math.cos(t) ** n) * math.tan(t / 2), 0, ts)[0]
    symmetric = (1 / math.tan(ts / 2)) ** 2 * i * i / (1 - math.cos(ts) ** (n + 1))
    assert symmetric / pb.offset_taper(fd, h, n) > 1.15


def test_the_offset_wins_on_gain_only_at_long_focal_lengths(dish):
    """The blockage it removes is worth 0.09 dB. At short and moderate focal
    lengths the tilted illumination costs more; only long, gently offset dishes
    come out ahead, and then by hundredths of a dB."""
    for fd, h, ft in ((0.5, 0.6, -10.0), (0.4, 0.8, -12.0), (0.7, 1.0, -8.0), (0.3, 0.55, -10.0)):
        assert _design(dish, fd, h, ft).metrics["gain_advantage_over_prime_focus_db"] < 0.0
    ahead = _design(dish, 1.0, 0.6, -10.0).metrics["gain_advantage_over_prime_focus_db"]
    assert 0.0 < ahead < 0.08


def test_the_old_default_budget_counted_spillover_twice(dish):
    d = _design(dish, 0.5, 0.6, -10.0)
    assert 10 * math.log10(d.get("eta_ill") * d.get("eta_spill") / (0.82 * 0.85)) > 0.4
