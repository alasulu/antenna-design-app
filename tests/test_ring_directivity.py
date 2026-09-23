"""Annular ring patch directivity, against the integral it was fitted to.

The ring is the one patch whose directivity is genuinely two-dimensional: it
depends on both radii separately, not on their ratio. Two designs sharing a
k0*b but reached from different (eps_r, b/a) pairs differ by as much as 13%, so
a fit in the ratio cannot work and a fit in the two radii can. These tests pin
both facts.
"""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import jv, yv


def _dJ1(x):
    return 0.5 * (jv(0, x) - jv(2, x))


def _dY1(x):
    return 0.5 * (yv(0, x) - yv(2, x))


def _exact(k0, k, a, b):
    """4*pi*U(0)/P_rad for the two-edge magnetic ring currents."""
    def R(r):
        return jv(1, k * r) * _dY1(k * a) - yv(1, k * r) * _dJ1(k * a)

    wa, wb = -a * R(a), b * R(b)

    def shape(t):
        sa, sb = k0 * a * math.sin(t), k0 * b * math.sin(t)
        e_th = wb * (jv(0, sb) - jv(2, sb)) + wa * (jv(0, sa) - jv(2, sa))
        e_ph = math.cos(t) * (wb * (jv(0, sb) + jv(2, sb))
                              + wa * (jv(0, sa) + jv(2, sa)))
        return e_th, e_ph

    integ = quad(lambda t: (shape(t)[0] ** 2 + shape(t)[1] ** 2) * math.sin(t),
                 0, math.pi / 2, limit=400)[0]
    return 4.0 * (wb + wa) ** 2 / integ


@pytest.mark.parametrize("eps_r,ratio", [
    (2.2, 1.5), (2.2, 2.0), (2.2, 3.0), (4.4, 2.0), (10.2, 2.0), (13.0, 1.2),
])
def test_ring_directivity_matches_its_own_quadrature(eps_r, ratio, registry):
    design = registry["annular_ring_patch"].synthesize(
        f0=2e9, eps_r=eps_r, h=1.6e-3, ratio=ratio)
    k0 = 2 * math.pi * 2e9 / 2.99792458e8
    k = design.get("k_diel")
    want = _exact(k0, k, design.get("a_in"), design.get("b_out"))
    assert design.metrics["directivity_linear"] == pytest.approx(want, rel=2e-3)


def test_ring_directivity_falls_with_permittivity(registry):
    """A higher-permittivity ring is smaller, so it must be less directive."""
    previous = None
    for eps_r in (2.2, 4.4, 6.15, 10.2, 13.0):
        d = registry["annular_ring_patch"].synthesize(
            f0=2e9, eps_r=eps_r, h=1.6e-3, ratio=2.0)
        value = d.metrics["directivity_linear"]
        if previous is not None:
            assert value < previous, f"eps_r={eps_r} did not fall"
        previous = value


def test_ring_sits_just_below_the_circular_patch(registry):
    """Both are broadside resonant discs on the same board, and the ring has the
    smaller radiating aperture, so it should land slightly lower - close enough
    to be a real cross-check between two independently derived patterns."""
    for eps_r in (2.2, 4.4, 10.2):
        ring = registry["annular_ring_patch"].synthesize(
            f0=2e9, eps_r=eps_r, h=1.6e-3, ratio=2.0)
        disc = registry["circular_patch"].synthesize(
            f0=2e9, eps_r=eps_r, h=1.6e-3)
        r, c = ring.metrics["directivity_linear"], disc.metrics["directivity_linear"]
        assert r < c, f"eps_r={eps_r}: ring {r} should be below disc {c}"
        assert r / c > 0.90, f"eps_r={eps_r}: ring {r} implausibly far below {c}"


def test_ring_directivity_is_not_a_function_of_the_radius_ratio_alone(registry):
    """The reason the fit is in the two radii. Designs sharing a k0*b but
    reached from different (eps_r, b/a) pairs have genuinely different
    directivities; a one-variable fit would have to average them away."""
    k0 = 2 * math.pi * 2e9 / 2.99792458e8
    seen = []
    for eps_r, ratio in ((4.4, 1.2), (8.0, 2.5)):
        d = registry["annular_ring_patch"].synthesize(
            f0=2e9, eps_r=eps_r, h=1.6e-3, ratio=ratio)
        seen.append((k0 * d.get("b_out"), d.metrics["directivity_linear"]))
    (kb1, d1), (kb2, d2) = seen
    assert abs(kb1 / kb2 - 1) < 0.10, "these two should share a k0*b"
    assert abs(d1 / d2 - 1) > 0.05, (
        "and still differ in directivity, which is what makes it 2-D")
