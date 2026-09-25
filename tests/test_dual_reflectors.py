"""Cassegrain and Gregorian, derived through an equivalent paraboloid that was traced.

Both dual-reflector specs asserted taper 0.85 and spillover 0.90, took blockage
from the uniform-illumination (1 - (Ds/D)^2)^2, and quoted 70 lambda/D. The
equivalent-paraboloid theorem says a dual reflector illuminates its aperture as
a single paraboloid of focal length M*F would; here it is checked by tracing
rays through the real hyperboloid and ellipsoid, and every efficiency the specs
now derive is compared with one computed from those traced rays.
"""
from __future__ import annotations

import math

import pytest

from otahub.num import paraboloid as pb

F0 = 1e10
LAM = 2.99792458e8 / F0
KINDS = [("cassegrain", "cassegrain"), ("gregorian_dual_reflector", "gregorian")]


def _n(et, M, fd):
    tf = 2 * math.atan(1 / (4 * M * fd))
    return 2 * (et / 20 - math.log10((1 + math.cos(tf)) / 2)) / math.log10(math.cos(tf))


def _design(registry, key, M, fd, dsd, et):
    return registry[key].synthesize(f0=F0, D=3.0, Ds=3.0 * dsd, f_over_D=fd, magnification=M,
                                    edge_taper_db=et, eps_rms=0.0)


@pytest.mark.parametrize("kind", ["cassegrain", "gregorian"])
@pytest.mark.parametrize("M", [1.5, 2.5, 4.0, 9.0])
def test_the_traced_rays_are_the_equivalent_paraboloids(kind, M):
    """tan(psi/2) = M tan(theta/2) ray by ray - which is what licenses the rest."""
    for th in (0.02, 0.1, 0.25, 0.4):
        psi = pb.dual_trace(kind, M, th)
        assert math.tan(psi / 2) / math.tan(th / 2) == pytest.approx(M, rel=1e-9)


@pytest.mark.parametrize("key,kind", KINDS)
@pytest.mark.parametrize("M,fd,dsd,et", [(3.0, 0.38, 0.08, -10.2), (5.0, 0.32, 0.12, -12.6),
                                         (2.2, 0.45, 0.17, -7.5)])
def test_the_efficiencies_are_the_traced_rays(registry, key, kind, M, fd, dsd, et):
    d = _design(registry, key, M, fd, dsd, et)
    taper, spill, block = pb.dual_efficiencies(kind, M, fd, dsd, _n(et, M, fd))
    assert d.get("eta_ill") == pytest.approx(taper, rel=1e-5)
    assert d.get("eta_spill") == pytest.approx(spill, rel=1e-5)
    assert d.metrics["eta_blockage"] == pytest.approx(block, abs=3e-4)


@pytest.mark.slow
@pytest.mark.parametrize("fdeq,dsd,et", [(1.1, 0.07, -9.5), (1.6, 0.13, -12.3), (2.2, 0.03, -6.2),
                                         (2.8, 0.18, -17.7), (0.9, 0.11, -4.4)])
def test_beamwidth_and_peak_sidelobe_are_the_blocked_transforms(registry, fdeq, dsd, et):
    """Off the fit's grid in all three variables."""
    d = _design(registry, "cassegrain", 4.0, fdeq / 4.0, dsd, et)
    hp, peak = pb.blocked_beam(et, fdeq, dsd)
    assert d.metrics["hpbw_deg"] == pytest.approx(hp * LAM / 3.0, abs=0.02 * LAM / 3.0)
    assert d.metrics["peak_sidelobe_db"] == pytest.approx(peak, abs=0.3)


@pytest.mark.parametrize("dsd", [0.05, 0.1, 0.15, 0.2])
def test_the_uniform_blockage_rule_understated_the_loss(registry, dsd):
    """The subreflector shadows the brightest part of a tapered aperture."""
    d = _design(registry, "cassegrain", 4.0, 0.35, dsd, -11.0)
    assert d.metrics["eta_blockage"] < (1 - dsd ** 2) ** 2


def test_with_no_blockage_the_beam_is_the_prime_focus_dishs():
    """Two transforms written separately, meeting where the geometries coincide."""
    fd = 0.8
    assert pb.blocked_beam(-11.0, fd, 0.0) == pytest.approx(pb.beam(pb.feed_exponent(-11.0, fd), fd), abs=1e-3)


def test_cassegrain_and_gregorian_are_the_same_antenna_to_geometric_optics(registry):
    c = _design(registry, "cassegrain", 4.0, 0.35, 0.1, -11.0)
    g = _design(registry, "gregorian_dual_reflector", 4.0, 0.35, 0.1, -11.0)
    for m in ("aperture_efficiency", "hpbw_deg", "peak_sidelobe_db"):
        assert c.metrics[m] == pytest.approx(g.metrics[m], rel=1e-12)
