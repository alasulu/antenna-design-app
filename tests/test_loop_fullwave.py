"""The full-wave circular loop, against two solvers that share no machinery.

`one_wavelength_circular_loop` used to assert a fixed 1.09 lambda resonance,
100 ohm and 3.4 dBi, none of which survived contact with a solver. What is in
the spec now comes from a method-of-moments solve AND from an independent
Fourier-mode solution of the same loop - the modal one expands the current in
exp(j n phi) and never discretises the geometry at all. This file holds both
to account, and pins the one thing the old spec had backwards: a THICKER
conductor resonates at a LONGER circumference.
"""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import jv

from otahub.num import loop_modal, mom

K = 2.0 * math.pi
NOMINAL = 1.09          # the nominal circumference the thickness parameter uses


# ------------------------------------------------- the independent solver
#
# It lives in otahub/num/loop_modal.py, beside the MoM, because it is a
# reference model rather than a test fixture - the quad-loop tests anchor
# themselves to it too.

_modal_z = loop_modal.input_impedance
_modal_resonance = loop_modal.resonant_circumference


def _modal_directivity(circumference: float, a: float) -> float:
    m = mom.loop(circumference, a, 96)
    nodes = np.array([m.node_of(n) for n in range(m.n_basis)])
    phi = np.arctan2(nodes[:, 1], nodes[:, 0])
    cur = loop_modal.current_at(circumference, a, phi)
    return mom.directivity_towards(mom.MoMSolution(m, cur.astype(complex), 0),
                                   0.0, 0.0, 90, 90)


def _design(registry, b_over_lambda):
    f0 = 3e8
    return registry["one_wavelength_circular_loop"].synthesize(
        f0=f0, b=b_over_lambda * 2.99792458e8 / f0)


RADII = [1e-4, 1e-3, 3e-3]


# ------------------------------------------------------ the two agree

@pytest.mark.parametrize("a", RADII)
def test_modal_and_mom_agree_on_impedance(a):
    """A modal expansion and a discretised EFIE, at a circumference where the
    loop is emphatically not resonant, so both real and imaginary parts count."""
    modal = _modal_z(NOMINAL, a, 80)
    got = mom.solve(mom.loop(NOMINAL, a, 96), feed=0).input_impedance
    assert got.real == pytest.approx(modal.real, rel=0.01)
    assert got.imag == pytest.approx(modal.imag, rel=0.05, abs=3.0)


# -------------------------------------------------- the spec against them

@pytest.mark.parametrize("a", RADII)
def test_spec_resonance_matches_the_modal_solve(a, registry):
    d = _design(registry, a)
    want, _ = _modal_resonance(a)
    lam = 2.99792458e8 / 3e8
    assert d.get("resonant_circumference_m") / lam == pytest.approx(want, rel=3e-3)


@pytest.mark.parametrize("a", RADII)
def test_spec_resistance_matches_the_modal_solve(a, registry):
    d = _design(registry, a)
    _, want = _modal_resonance(a)
    assert d.metrics["input_resistance_ohm"] == pytest.approx(want, rel=5e-3)
    assert want > 130.0, "and it is nowhere near the 100 ohm this spec asserted"


@pytest.mark.parametrize("a", RADII)
def test_spec_directivity_matches_the_modal_current(a, registry):
    d = _design(registry, a)
    cres, _ = _modal_resonance(a)
    assert d.metrics["directivity_linear"] == pytest.approx(
        _modal_directivity(cres, a), rel=5e-3)


# --------------------------------------------------- what was backwards

def test_a_thicker_conductor_resonates_at_a_longer_circumference(registry):
    """The old spec's thickness correction ran the other way, and its validity
    block said so in words: 'thicker conductors resonate at a shorter
    circumference'. Both solvers say the opposite, so this is pinned."""
    got = [_design(registry, a).get("resonant_circumference_m") for a in RADII]
    assert got[0] < got[1] < got[2], got
    modal = [_modal_resonance(a)[0] for a in RADII]
    assert modal[0] < modal[1] < modal[2], modal


# ------------------------------------------------------- sanity in context

def test_the_loop_beats_a_half_wave_dipole_by_about_one_and_a_half_db(registry):
    for a in RADII:
        d = _design(registry, a)
        assert d.metrics["gain_over_dipole_db"] == pytest.approx(1.5, abs=0.15)
        assert d.metrics["directivity_linear"] > 1.64093


def test_broadside_is_the_main_beam():
    """Not obvious: with a UNIFORM current a loop this size peaks in its own
    plane instead. It is the cosinusoidal current that puts the beam
    broadside, so this checks the current and not just the geometry."""
    for a in RADII:
        cres, _ = _modal_resonance(a)
        sol = mom.solve(mom.loop(cres, a, 96), feed=0)
        broadside = mom.directivity_towards(sol, 0.0, 0.0, 90, 90)
        assert mom.directivity(sol, 90, 90) == pytest.approx(broadside, rel=3e-3)


def test_far_field_routine_on_a_circular_geometry():
    """The exact uniform-current loop directivity, 2*kb*J1(kb)^2 / int_0^2kb J2,
    checks the radiation integral on a ring rather than a straight wire."""
    for c in (0.3, 0.6, 1.0):
        m = mom.loop(c, 1e-3, 96)
        sol = mom.MoMSolution(m, np.ones(m.n_basis, dtype=complex), 0)
        want = 2 * c * jv(1, c) ** 2 / quad(lambda x: jv(2, x), 0, 2 * c)[0]
        assert mom.directivity(sol, 100, 100) == pytest.approx(want, rel=2e-3)


def test_a_circle_beats_a_square_of_the_same_perimeter():
    """More enclosed area for the same conductor. A cross-check between two
    geometries rather than against a reference."""
    per = 1.05
    circle = mom.directivity_towards(
        mom.solve(mom.loop(per, 1e-3, 96), feed=0), 0.0, 0.0, 80, 80)
    s = per / 4.0
    pts = []
    corners = [(-s / 2, -s / 2), (s / 2, -s / 2), (s / 2, s / 2), (-s / 2, s / 2)]
    for i in range(4):
        x0, y0 = corners[i]
        x1, y1 = corners[(i + 1) % 4]
        for j in range(16):
            t = j / 16
            pts.append((x0 + t * (x1 - x0), y0 + t * (y1 - y0), 0.0))
    msq = mom.WireModel([mom.Wire(np.array(pts), 1e-3, closed=True)])
    nodes = np.array([msq.node_of(n) for n in range(msq.n_basis)])
    feed = int(np.argmin(np.linalg.norm(nodes - np.array([0.0, -s / 2, 0.0]), axis=1)))
    square = mom.directivity_towards(mom.solve(msq, feed), 0.0, 0.0, 80, 80)
    assert 1.02 < circle / square < 1.20, (circle, square)


def test_small_loop_formulas_are_useless_here(registry):
    """A validity note the spec makes, checked rather than asserted: at one
    wavelength the uniform-current results are wrong by an order of magnitude
    and point the beam 90 degrees away."""
    d = _design(registry, 1e-3)
    c_over_lambda = d.get("resonant_circumference_m") * 3e8 / 2.99792458e8
    small_rr = 20.0 * math.pi ** 2 * c_over_lambda ** 4
    assert small_rr / d.metrics["input_resistance_ohm"] > 1.5
    assert d.metrics["directivity_linear"] > 1.4 * 1.5
