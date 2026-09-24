"""The square quad loop, anchored to the circle by a polygon sequence.

A square has no modal solution to check it against, so the verification here is
by continuity rather than by a second formula. A regular N-gon of a given
perimeter tends to the circle of the same perimeter, and the circle DOES have
an independent Fourier-mode solution. If the N-gon sequence lands on that
modal answer, then the corner handling, the non-collinear segments and the
closed-loop wrap are all anchored to something outside this solver, and the
square is simply N = 4 of the same validated sequence.

The three orderings between circle and square are worth as much as any single
number: the square encloses less area for the same perimeter, so it must
resonate LONGER, and show a LOWER resistance and a LOWER directivity. All three
come out that way without being arranged to.
"""
from __future__ import annotations

import functools
import math

import numpy as np
import pytest

from otahub.num import loop_modal, mom

pytestmark = pytest.mark.slow   # full-wave solves; skip with -m 'not slow'

NOMINAL = 1.0218      # the empirical perimeter the spec synthesises
LAMBDA = 2.99792458e8 / 3e8
RADII = [1e-4, 1e-3, 3e-3]


def _polygon(perimeter: float, nsides: int, a: float, per_side: int = 20,
             rotate: float = 0.0):
    """Regular N-gon of the given PERIMETER in the x-y plane."""
    side = perimeter / nsides
    r = side / (2.0 * math.sin(math.pi / nsides))
    ang = np.arange(nsides) * 2 * math.pi / nsides + rotate
    verts = np.stack([r * np.cos(ang), r * np.sin(ang), np.zeros(nsides)], axis=1)
    pts = []
    for i in range(nsides):
        p, q = verts[i], verts[(i + 1) % nsides]
        for j in range(per_side):
            pts.append(p + (j / per_side) * (q - p))
    model = mom.WireModel([mom.Wire(np.array(pts), a, closed=True)])
    return model, (per_side // 2) - 1      # feed at the midpoint of side 0


def _square(perimeter: float, a: float, per_side: int = 20):
    return _polygon(perimeter, 4, a, per_side, rotate=math.pi / 4)


@functools.lru_cache(maxsize=32)
def _square_resonance(a: float, per_side: int = 20, steps: int = 24):
    lo, hi = 0.95, 1.30
    for _ in range(steps):
        mid = 0.5 * (lo + hi)
        model, feed = _square(mid, a, per_side)
        if mom.solve(model, feed).input_impedance.imag < 0:
            lo = mid
        else:
            hi = mid
    model, feed = _square(lo, a, per_side)
    sol = mom.solve(model, feed)
    return lo, sol.input_impedance.real, mom.directivity_towards(
        sol, 0.0, 0.0, 90, 90)


def _design(registry, a):
    return registry["quad_loop_square"].synthesize(f0=3e8, b=a * LAMBDA)


# ------------------------------------------------ the anchor to the circle

def test_a_polygon_sequence_converges_on_the_modal_circle():
    """This is what licenses the square. The N-gon uses exactly the code path
    the square does - corners, non-collinear segments, a closed wrap - and its
    limit is a shape with an independent analytic solution."""
    per, a = 1.05, 1e-3
    got = []
    for n in (4, 8, 16, 48):
        model, feed = _polygon(per, n, a, max(96 // n, 2))
        got.append(mom.solve(model, feed).input_impedance)
    modal = loop_modal.input_impedance(per, a, 80)
    # monotone approach, and the last one is there
    errs = [abs(z - modal) for z in got]
    assert errs[0] > errs[1] > errs[2] > errs[3], errs
    assert got[-1].real == pytest.approx(modal.real, rel=0.01)
    assert got[-1].imag == pytest.approx(modal.imag, abs=3.0)


# ------------------------------------------------------- the spec's numbers

@pytest.mark.parametrize("a", RADII)
def test_spec_resonant_perimeter(a, registry):
    want, _, _ = _square_resonance(a)
    got = _design(registry, a).get("resonant_perimeter_m") / LAMBDA
    assert got == pytest.approx(want, rel=3e-3)


@pytest.mark.parametrize("a", RADII)
def test_spec_resistance_and_directivity(a, registry):
    _, want_r, want_d = _square_resonance(a)
    d = _design(registry, a)
    assert d.metrics["input_resistance_ohm"] == pytest.approx(want_r, rel=5e-3)
    assert d.metrics["directivity_linear"] == pytest.approx(want_d, rel=5e-3)


def test_the_nominal_perimeter_is_not_resonant():
    """1.0218 lambda is the ham-radio 1005/f rule - insulated wire, near ground,
    inside an array. An isolated bare-wire loop is strongly capacitive there,
    which is why the spec could not have had both that perimeter and a purely
    resistive feed."""
    for a in RADII:
        model, feed = _square(NOMINAL, a, 20)
        z = mom.solve(model, feed).input_impedance
        assert z.imag < -50.0, f"b={a}: {z}"


def test_a_thicker_conductor_resonates_at_a_longer_perimeter(registry):
    got = [_design(registry, a).get("resonant_perimeter_m") for a in RADII]
    assert got[0] < got[1] < got[2], got


# ------------------------------------------- the constant that held up

def test_three_point_three_dbi_was_right_after_all(registry):
    """Of the loop family's asserted constants this is the one that survived:
    3.3 dBi is essentially exact on thin wire. Worth a test of its own, so that
    a later refactor cannot quietly move it."""
    d = _design(registry, 1e-4)
    assert d.metrics["directivity_dbi"] == pytest.approx(3.3, abs=0.02)


# --------------------------------------------- circle against square

def test_the_square_loses_to_the_circle_on_all_three_counts(registry):
    """Same perimeter, less enclosed area: it must resonate longer, and show a
    lower resistance and a lower directivity. Three independent orderings, none
    of them arranged for."""
    for a in RADII:
        sq = _design(registry, a)
        ci = registry["one_wavelength_circular_loop"].synthesize(
            f0=3e8, b=a * LAMBDA)
        assert sq.get("resonant_perimeter_m") > ci.get("resonant_circumference_m")
        assert sq.metrics["input_resistance_ohm"] < ci.metrics["input_resistance_ohm"]
        assert sq.metrics["directivity_linear"] < ci.metrics["directivity_linear"]


def test_the_quarter_wave_transformer_still_wants_seventy_five_ohm_cable(registry):
    """sqrt(R*50) with the derived resistance, rather than with an assumed 120."""
    for a in RADII:
        z0 = _design(registry, a).metrics["matching_stub_z0_ohm"]
        assert 78.0 < z0 < 85.0, z0
