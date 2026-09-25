"""The Fresnel zone plate, from the feed through the plate by Kirchhoff diffraction.

Its gain was the first-order grating efficiency (1/pi^2, 4/pi^2, 8/pi^2) times a
perfectly illuminated aperture: no feed taper, no spillover. For a phase plate
that left the budget 0.9-1.4 dB optimistic; for an opaque plate the error could
go either way, because its unfocused zero order passes the open zones too.
Every reference here is `otahub.num.zone_plate`, which integrates over a
Cartesian grid on the plate with each point's own path length - not the
feed-angle integral the spec evaluates.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from otahub.num import zone_plate as zp

C0 = 2.99792458e8


def _design(registry, f0, F, M, levels, taper):
    return registry["fresnel_zone_plate"].synthesize(f0=f0, F=F, M=M, phase_levels=levels,
                                                      feed_taper_db=taper)


def _arbiter_args(f0, F, M, taper):
    Fl = F * f0 / C0
    rim = math.atan(zp.zone_radius(M, Fl) / Fl)
    return Fl, taper / (10 * math.log10(math.cos(rim)))


@pytest.mark.parametrize("f0,F,M,levels,taper", [
    (30e9, 0.15, 4, 1, -10.0), (30e9, 0.15, 4, 2, -13.0), (30e9, 0.15, 5, 4, -7.0),
    (24e9, 0.2, 11, 1, -15.0), (77e9, 0.05, 6, 2, -9.0), (10e9, 0.6, 7, 4, -12.0)])
def test_aperture_efficiency_is_the_plate_integral(registry, f0, F, M, levels, taper):
    d = _design(registry, f0, F, M, levels, taper)
    Fl, n = _arbiter_args(f0, F, M, taper)
    D = 2 * zp.zone_radius(M, Fl)
    eta = zp.gain(Fl, M, levels, n, cells=1200) / (math.pi * D) ** 2
    assert d.metrics["aperture_efficiency"] == pytest.approx(eta, rel=1e-3)


def test_uniform_illumination_recovers_the_grating_efficiencies():
    """A distant isotropic-ish feed on many zones: the plate integral must fall
    back to 1/pi^2, 4/pi^2 and 8/pi^2, which is what the spec used to quote."""
    for levels, grating in ((1, 1 / math.pi ** 2), (2, 4 / math.pi ** 2), (4, 8 / math.pi ** 2)):
        F, M = 2000.0, 40
        D = 2 * zp.zone_radius(M, F)
        eta = zp.gain(F, M, levels, 0.0, cells=1500) / (math.pi * D) ** 2
        # an isotropic feed also spills almost everything past so small a rim
        spill = 1 - math.cos(math.atan(D / 2 / F))
        assert eta / spill == pytest.approx(grating, rel=0.01)


def test_phase_plates_factorise_and_opaque_plates_do_not(registry):
    """For a phase plate efficiency = grating x (the feed's taper and
    spillover); for an opaque one the zero order adds through the bright
    central zone, so it beats 1/pi^2 under a heavy taper with few zones."""
    base = _design(registry, 30e9, 0.15, 4, 4, -10.0)
    reversal = _design(registry, 30e9, 0.15, 4, 2, -10.0)
    ideal = base.metrics["aperture_efficiency"] / (8 / math.pi ** 2)
    assert reversal.metrics["aperture_efficiency"] == pytest.approx(4 / math.pi ** 2 * ideal, rel=0.01)
    opaque = _design(registry, 30e9, 0.15, 4, 1, -17.0)
    assert opaque.metrics["aperture_efficiency"] > 1.1 / math.pi ** 2


def test_the_old_gain_left_the_feed_out(registry):
    for levels in (2, 4):
        d = _design(registry, 30e9, 0.15, 4, levels, -10.0)
        old = d.metrics["focusing_efficiency"] * (math.pi * d.metrics["aperture_diameter_m"] * 30e9 / C0) ** 2
        assert 0.9 < 10 * math.log10(old) - d.metrics["gain_dbi"] < 1.1


@pytest.mark.slow
@pytest.mark.parametrize("f0,F,M,levels,taper", [
    (30e9, 0.15, 4, 1, -10.0), (30e9, 0.15, 5, 1, -10.0), (40e9, 0.12, 8, 2, -14.0),
    (20e9, 0.4, 13, 4, -6.0)])
def test_beamwidth_and_bandwidth_are_the_plate_integral(registry, f0, F, M, levels, taper):
    """Off the fitting grid in taper and geometry, both parities of the opaque plate."""
    d = _design(registry, f0, F, M, levels, taper)
    Fl, n = _arbiter_args(f0, F, M, taper)
    assert d.metrics["hpbw_deg"] == pytest.approx(zp.hpbw(Fl, M, levels, n, cells=800), rel=0.01)
    assert d.metrics["gain_bandwidth_1db"] == pytest.approx(
        zp.gain_bandwidth(Fl, M, levels, n, 1.0, cells=500), rel=0.01)


@pytest.mark.slow
def test_an_opaque_plate_beam_depends_on_the_parity_of_its_zone_count(registry):
    """An even M ends on a blocked zone, an odd one on an open zone; the
    zero order makes the beams differ. Both come from the 2-D integral."""
    Fl, n4 = _arbiter_args(30e9, 0.15, 4, -10.0)
    _, n5 = _arbiter_args(30e9, 0.15, 5, -10.0)
    hp4 = zp.hpbw(Fl, 4, 1, n4, cells=800) * 2 * zp.zone_radius(4, Fl)
    hp5 = zp.hpbw(Fl, 5, 1, n5, cells=800) * 2 * zp.zone_radius(5, Fl)
    assert hp4 / hp5 > 1.08
    rev4 = zp.hpbw(Fl, 4, 2, n4, cells=800) * 2 * zp.zone_radius(4, Fl)
    rev5 = zp.hpbw(Fl, 5, 2, n5, cells=800) * 2 * zp.zone_radius(5, Fl)
    assert abs(rev4 / rev5 - 1) < 0.03
