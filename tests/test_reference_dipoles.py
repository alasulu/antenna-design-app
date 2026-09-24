"""The two reference dipoles: which numbers are definitions and which are the wire.

`half_wave_dipole` and `resonant_dipole` are the archetypes everything else is
measured against, and both were built entirely on the induced-EMF closed forms,
Balanis 4-70 and 4-79. Those are exact - for an ASSUMED sinusoidal current on a
vanishingly thin wire. A delta-gap-fed wire of finite radius presents something
else, and this file pins the distinction in both directions:

- the canonical 73.08 + j42.52 stays exactly what it is, because dBd and half
  the textbooks rest on it, and a test says so;
- the driving-point values sit beside it, and a test says they differ, by how
  much, and that the difference is the current's shape rather than the mesh.

The last point is the one that matters. Handing the MoM the sinusoidal current
reproduces the closed form to 0.3%, so the gap is physics, not numerics.
"""
from __future__ import annotations

import functools
import math

import numpy as np
import pytest

from otahub.num import mom

NSEG = 64
LAM = 2.99792458e8 / 3e8
EMF_R = 73.0796


def _l_spec(aw):
    t = math.log10(1.0 / aw)
    return -0.0040566 * t ** 2 + 0.0378414 * t + 0.3994411


@functools.lru_cache(maxsize=16)
def _at_spec_length(aw):
    return mom.solve(mom.dipole(_l_spec(aw), aw, NSEG))


def _res(registry, aw):
    return registry["resonant_dipole"].synthesize(f0=3e8, aw=aw * LAM)


def _hw(registry, aw):
    return registry["half_wave_dipole"].synthesize(f0=3e8, aw=aw * LAM)


# --------------------------------------------------- the gap is the current

@pytest.mark.parametrize("aw", [1e-4, 1e-3])
def test_the_closed_form_is_the_sinusoidal_current_exactly(aw, registry):
    """Hand the MoM the assumed sinusoid and it reproduces the spec's
    induced-EMF resistance. So the disagreement with the driving point is the
    real current's shape, not a numerical artefact."""
    L = _l_spec(aw)
    m = mom.dipole(L, aw, NSEG)
    z = np.array([m.node_of(n)[2] for n in range(m.n_basis)])
    current = np.sin(mom.K * (0.5 * L - np.abs(z))).astype(complex)
    sol = mom.MoMSolution(m, current, int(np.argmax(np.abs(current))))
    rr = 2.0 * mom.radiated_power(sol, 90, 90) / abs(current).max() ** 2
    assert rr == pytest.approx(
        _res(registry, aw).metrics["radiation_resistance_ohm"], rel=1e-2)


# ------------------------------------------------------- resonant dipole

@pytest.mark.parametrize("aw", [1e-5, 1e-4, 1e-3, 3e-3])
def test_resonant_dipole_driving_point_matches_the_solver(aw, registry):
    want = _at_spec_length(aw).input_impedance.real
    assert _res(registry, aw).metrics["input_resistance_ohm"] == pytest.approx(
        want, rel=1e-2)


def test_a_resonant_dipole_really_is_about_72_ohm(registry):
    """The folk figure was right; the closed form quoted in its place was not.
    Flat to about an ohm across two and a half decades of wire radius."""
    got = [_res(registry, a).metrics["input_resistance_ohm"]
           for a in (1e-5, 1e-4, 1e-3, 2.7e-3)]
    assert all(71.0 < r < 75.0 for r in got), got
    closed = [_res(registry, a).metrics["radiation_resistance_ohm"]
              for a in (1e-5, 1e-4, 1e-3, 2.7e-3)]
    assert all(c < r for c, r in zip(closed, got))
    assert closed[-1] < 62.0, "the closed form falls away on thick wire; the wire does not"


def test_vswr_now_uses_what_the_antenna_presents(registry):
    """In 50 ohm the old figure said about 1.26:1 for aw = 1e-3; a real one is
    nearer 1.46:1. Optimistic in exactly the number checked first."""
    d = _res(registry, 1e-3)
    assert d.metrics["vswr_in_50_ohm"] == pytest.approx(
        d.metrics["input_resistance_ohm"] / 50.0, rel=1e-9)
    assert d.metrics["vswr_in_50_ohm"] > 1.4


def test_residual_reactance_is_small_but_not_zero():
    """'Zero by construction' was true of the model, not the wire: the length
    was tuned to the closed form. A real wire there shows a few ohms."""
    xs = [_at_spec_length(a).input_impedance.imag for a in (1e-5, 1e-4, 1e-3, 5e-3)]
    assert all(abs(x) < 10.0 for x in xs), xs
    assert max(abs(x) for x in xs) > 1.0, xs


def test_bandwidth_triples_across_the_gauges(registry):
    """The quantitative form of 'fat dipoles are broadband'."""
    thin = _res(registry, 1e-5).metrics["fractional_bandwidth_vswr2"]
    fat = _res(registry, 5e-3).metrics["fractional_bandwidth_vswr2"]
    assert 0.05 < thin < 0.065
    assert fat / thin > 2.8


def test_resonant_directivity_is_just_below_the_half_wave_value(registry):
    """It used to borrow the half-wave 1.6409. Close, but a different antenna."""
    for aw in (1e-5, 1e-3, 5e-3):
        d = _res(registry, aw).metrics["directivity_linear"]
        assert 1.630 < d < 1.6409, d


# -------------------------------------------------------- half-wave dipole

def test_the_canonical_73_ohm_is_untouched(registry):
    """The reference dBd rests on. It must not move, whatever the wire."""
    for aw in (1e-5, 1e-3):
        d = _hw(registry, aw)
        assert d.metrics["radiation_resistance_ohm"] == pytest.approx(EMF_R, rel=1e-4)
        assert d.metrics["input_reactance_ohm"] == pytest.approx(42.5152, rel=1e-3)
        assert d.metrics["directivity_linear"] == pytest.approx(1.6409, rel=1e-4)


@pytest.mark.parametrize("aw", [1e-5, 1e-4, 1e-3, 2.7e-3])
def test_half_wave_driving_point_matches_the_solver(aw, registry):
    want = mom.input_impedance(mom.dipole(0.5, aw, NSEG))
    d = _hw(registry, aw)
    assert d.metrics["input_resistance_driving_point_ohm"] == pytest.approx(
        want.real, rel=1.5e-2)
    assert d.metrics["input_reactance_driving_point_ohm"] == pytest.approx(
        want.imag, rel=3e-2)


def test_the_driving_point_approaches_73_only_as_the_wire_vanishes(registry):
    """Monotone in the radius, and still 7% high on 1e-5 lambda wire."""
    got = [_hw(registry, a).metrics["input_resistance_driving_point_ohm"]
           for a in (2.7e-3, 1e-3, 1e-4, 1e-5)]
    assert got[0] > got[1] > got[2] > got[3] > EMF_R, got
    assert got[-1] / EMF_R - 1 > 0.05
