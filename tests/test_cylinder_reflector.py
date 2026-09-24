"""The parabolic cylinder, derived from its line feed and checked on its aperture.

`cylindrical_parabolic` asserted 0.81 x 0.88 for taper and spillover and borrowed
the dish's 70 lambda/W beamwidth. Under a line feed the wave is cylindrical and
the aperture one-dimensional: the product reaches 0.90 at f/W = 0.4, the best
edge taper is a lighter -8 dB rather than the dish's -11, and the beam is
58 lambda/W. Every check comes from `otahub.num.paraboloid`, which integrates
over the aperture coordinate rather than the feed angle the spec uses.
"""
from __future__ import annotations

import math

import pytest

from otahub.num import paraboloid as pb

F0 = 3e9
LAM = 2.99792458e8 / F0


@pytest.fixture(scope="module")
def cyl(registry):
    return registry["cylindrical_parabolic"]


def _design(cyl, fw, et):
    return cyl.synthesize(f0=F0, W=1.0, Lc=3.0, f_over_W=fw, edge_taper_db=et, eps_rms=0.0)


@pytest.mark.parametrize("fw,et", [(0.3, -6.0), (0.33, -9.5), (0.4, -8.0), (0.42, -12.3),
                                   (0.55, -6.2), (0.75, -17.7), (0.95, -4.4)])
def test_taper_and_spillover_are_the_aperture_integrals(cyl, fw, et):
    d = _design(cyl, fw, et)
    taper, spill = pb.cylinder_efficiencies(pb.cylinder_feed_exponent(et, fw), fw)
    assert d.get("eta_ill") == pytest.approx(taper, rel=1e-5)
    assert d.get("eta_spill") == pytest.approx(spill, rel=2e-4)


@pytest.mark.slow
@pytest.mark.parametrize("fw,et", [(0.33, -9.5), (0.42, -12.3), (0.55, -6.2), (0.75, -17.7),
                                   (0.95, -4.4), (0.38, -19.3)])
def test_beamwidth_and_peak_sidelobe_are_the_fourier_transforms(cyl, fw, et):
    """Off the fit's grid in both variables."""
    d = _design(cyl, fw, et)
    hp, peak = pb.cylinder_beam(pb.cylinder_feed_exponent(et, fw), fw)
    assert d.metrics["hpbw_focusing_plane_deg"] == pytest.approx(hp * LAM, abs=0.05 * LAM)
    assert d.metrics["peak_sidelobe_focusing_plane_db"] == pytest.approx(peak, abs=0.15)


def test_a_nearly_uniform_aperture_is_the_uniform_line_source():
    """n = 0 on a long focal length: 50.8 lambda/W and -13.3 dB, as it must be."""
    hp, peak = pb.cylinder_beam(0.0, 2.0)
    assert hp == pytest.approx(50.8, abs=0.1)
    assert peak == pytest.approx(-13.26, abs=0.1)


@pytest.mark.parametrize("fw", [0.3, 0.4, 0.5, 0.7])
def test_a_line_aperture_wants_a_lighter_taper_than_a_dish(cyl, fw):
    tapers = [-4.0 - 0.25 * i for i in range(45)]              # -4 to -15 dB
    best = max(tapers, key=lambda et: _design(cyl, fw, et).metrics["aperture_efficiency"])
    assert -9.0 <= best <= -6.0


def test_the_old_budget_understated_the_gain_by_a_decibel(cyl):
    d = _design(cyl, 0.4, -8.0)
    assert 10 * math.log10(d.get("eta_ill") * d.get("eta_spill") / (0.81 * 0.88)) > 0.95
