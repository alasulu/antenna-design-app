"""Corner reflectors against an exact image-theory method of moments.

A dipole in an infinite PEC corner is exactly its image set, and solving that
set as one MoM problem gives the real current rather than an assumed sinusoid
(`otahub.num.corner`). The specs' directivity - array factor times dipole
directivity times a resistance ratio from induced-EMF mutual terms - is
checked against it here and holds to a few hundredths of a dB. The input
IMPEDANCE is another matter: the induced-EMF figure is the resistance of an
assumed sinusoid on an exactly half-wave element, and the real driving point
of that element differs by up to 40% in resistance and carries a reactance the
old spec did not report.
"""
from __future__ import annotations

import functools
import math

import pytest

from otahub.num import corner

pytestmark = pytest.mark.slow   # full-wave solves; skip with -m 'not slow'

F0 = 2.99792458e8               # lambda = 1 m, so metres are wavelengths


@functools.lru_cache(maxsize=64)
def _solve(spacing, deg, radius=1e-3):
    return corner.solve_corner(spacing, deg, radius)


@pytest.mark.parametrize("key,deg,spacing", [
    ("corner_reflector_90", 90, s) for s in (0.1, 0.25, 0.35, 0.5, 0.6, 0.75)] + [
    ("corner_reflector_60", 60, s) for s in (0.35, 0.5, 0.7, 0.8)])
def test_the_directivity_is_the_exact_image_solutions(registry, key, deg, spacing):
    d = registry[key].synthesize(f0=F0, S_over_lambda=spacing)
    exact = 10 * math.log10(_solve(spacing, deg).directivity_boresight())
    assert d.metrics["directivity_dbi"] == pytest.approx(exact, abs=0.06)


@pytest.mark.parametrize("deg,spacing", [(90, 0.25), (90, 0.5), (60, 0.7)])
def test_the_images_conserve_power(deg, spacing):
    """Every wedge of the image system radiates what the real feed delivers."""
    assert _solve(spacing, deg).power_balance() == pytest.approx(1.0, abs=1e-4)


def test_the_directivity_is_nearly_flat_below_half_a_wavelength(registry):
    """The validity note's claim, 12.6 dBi near the vertex falling to 11.8 at
    lambda/2 - with the real current as well as the idealised one."""
    d = [10 * math.log10(_solve(s, 90).directivity_boresight()) for s in (0.1, 0.25, 0.5)]
    assert d[0] == pytest.approx(12.6, abs=0.1) and d[2] == pytest.approx(11.8, abs=0.1)
    assert d[0] > d[1] > d[2]


def test_lambda_spacing_is_a_boresight_null():
    s = _solve(0.999, 90)
    assert s.directivity_boresight() < 0.05


@pytest.mark.parametrize("key,deg,spacing,radius", [
    ("corner_reflector_90", 90, 0.1875, 7e-4), ("corner_reflector_90", 90, 0.4375, 1.5e-4),
    ("corner_reflector_90", 90, 0.7125, 4e-3), ("corner_reflector_60", 60, 0.3125, 6e-4),
    ("corner_reflector_60", 60, 0.6625, 2.5e-3)])
def test_the_driving_point_impedance_is_the_exact_solutions(registry, key, deg, spacing, radius):
    """Off the fit's grid in both spacing and radius."""
    d = registry[key].synthesize(f0=F0, S_over_lambda=spacing, aw_over_lambda=radius)
    z = _solve(spacing, deg, radius).input_impedance
    got = complex(d.metrics["input_resistance_ohm"], d.metrics["input_reactance_ohm"])
    assert abs(got - z) < 3.0


def test_the_induced_emf_figure_was_not_what_the_feed_sees(registry):
    """126 ohm quoted at S = lambda/2; the real element presents about 151 + j49."""
    d = registry["corner_reflector_90"].synthesize(f0=F0, S_over_lambda=0.5)
    z = _solve(0.5, 90).input_impedance
    assert z.real / d.metrics["induced_emf_resistance_ohm"] > 1.15
    assert z.imag > 40.0
