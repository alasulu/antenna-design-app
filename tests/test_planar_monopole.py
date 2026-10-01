"""The planar monopoles' sizing rule against the planar MoM.

f_low = 7.2e7/(L + r_eq + p) with the equivalent cylinder's radius from equal
AREA (Agrawall, Kumar and Ray 1998; Ray 2008): W/(2 pi) for a plate, r/4 for a
disc. The specs used W/2 and r, and drew every plate about a quarter too small.
"""
from __future__ import annotations

import math

import pytest
from scipy.optimize import brentq

from otahub.num import rwg

C = 2.99792458e8


def _vswr(z, z0=50.0):
    g = abs((z - z0) / (z + z0))
    return (1 + g) / (1 - g)


@pytest.mark.slow
def test_a_square_plate_designed_for_f_low_starts_its_band_there(registry):
    d = registry["planar_monopole_rectangular"].synthesize(f_low=1.5e9, aspect=1.0, p_gap=0.001)
    L, W, p = d.get("Lp"), d.get("Wp"), d.get("p_gap")
    edge = brentq(lambda f: _vswr(rwg.plate_monopole(W * f / C, L * f / C, p * f / C, 0.002 * f / C)) - 2.0,
                  1.2e9, 2.0e9, xtol=2e6)
    assert edge == pytest.approx(1.5e9, rel=0.05)
    # the old half-width reading would have put this plate's edge at
    old = 7.2e7 / (L + W / 2 + p)
    assert abs(old / edge - 1) > 0.15


def test_the_disc_and_the_plate_use_the_same_equal_area_rule(registry):
    disc = registry["planar_monopole_circular"].synthesize(f_low=1e9, p_gap=0.001)
    r = disc.get("r_disc")
    assert disc.get("equivalent_cylinder_radius_m") == pytest.approx(r / 4)
    assert 2 * math.pi * (r / 4) * (2 * r) == pytest.approx(math.pi * r * r)
    plate = registry["planar_monopole_rectangular"].synthesize(f_low=1e9, aspect=1.5, p_gap=0.001)
    W, L = plate.get("Wp"), plate.get("Lp")
    assert 2 * math.pi * plate.get("equivalent_cylinder_radius_m") * L == pytest.approx(W * L)
