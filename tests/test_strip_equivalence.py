"""A flat strip against the tube of radius w/4 it is said to be equivalent to.

The slot family takes its length and resistance from the complementary strip
dipole through the round-wire laws at a = w/4. `otahub.num.strip` solves the
strip itself, twice over (a spectral Galerkin with edge-exact currents, and a
rooftop mixed-potential MoM that also takes a gap feed), and the open tube with
its rims resolved. What these tests pin: the equivalence holds for the body but
not the ends - the strip resonates about 0.1 w longer - and fed across the same
gap the two part further in resistance when the gap is small and the strip wide.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import dblquad

from otahub.num import strip as st
from otahub.num.patch_fdtd import extrapolate

DATA = json.loads((Path(__file__).parent / "data" / "strip_equivalence.json").read_text())


def _fed(r):
    ms = r["strip_meshes"]
    L = extrapolate([m["nx"] for m in ms], [m["L"] for m in ms])
    R = extrapolate([m["nx"] for m in ms], [m["R"] for m in ms])
    return L, R


@pytest.mark.parametrize("p", [(0.03, 0.001, -0.01, 0.02, -0.005, 0.003), (0.0, 0.003, -0.01, 0.0, -0.005, 0.003),
                               (-0.01, -0.005, -0.01, 0.02, -0.005, 0.003)])
def test_the_cell_integral_matches_quadrature(p):
    px, py, x1, x2, y1, y2 = p
    k = 2 * math.pi

    def part(fn):
        return dblquad(lambda y, x: fn(np.exp(-1j * k * math.hypot(x - px, y - py))) / (4 * math.pi * math.hypot(x - px, y - py)),
                       x1, x2, y1, y2, epsabs=1e-13)[0]

    want = part(np.real) + 1j * part(np.imag)
    assert st.rect_potential(px, py, x1, x2, y1, y2) == pytest.approx(want, rel=1e-6)


def test_the_two_strip_solvers_agree():
    """Natural resonance at w = 0.04: the rooftop MoM taken to zero cell size
    against the spectral one - nothing shared but the geometry."""
    c = DATA["rooftop_check"]
    rooftop = extrapolate([m["nx"] for m in c["meshes"]], [m["L"] for m in c["meshes"]])
    assert rooftop == pytest.approx(c["spectral"], rel=1e-3)
    assert c["spectral_larger_basis"] == pytest.approx(c["spectral"], rel=5e-4)


def test_the_strip_outlasts_its_tube_by_an_end_effect():
    """The strip resonates longer at every width, by a length proportional to w:
    the ends differ (a tube's rim is 1.57 w round), not the body."""
    nat = DATA["natural"]
    rel = [r["L_strip"] / r["L_tube"] - 1 for r in nat]
    assert all(x > 0 for x in rel) and rel == sorted(rel)
    for r in nat:
        if r["w"] >= 0.02:
            assert 0.08 < (r["L_strip"] - r["L_tube"]) / r["w"] < 0.11
    assert rel[-1] > 0.025                                          # 2.9% at w = 0.12


def test_an_ungraded_tube_converges_slowly_to_the_graded_one():
    rim = DATA["rim"]
    un = [r["L"] for r in rim["ungraded"]]
    graded = rim["graded"][-1]["L"]
    assert un == sorted(un, reverse=True) and un[-1] > graded
    assert rim["graded"][0]["L"] == pytest.approx(graded, abs=5e-5)
    # a first-order extrapolation from two ungraded meshes still lands close
    assert 2 * un[2] - un[1] == pytest.approx(graded, rel=1e-3)


def test_thin_wire_theory_sits_above_the_rim_resolved_tube():
    tw = DATA["thin_wire"]
    converged = tw["meshes"][-1]["L"]
    assert tw["meshes"][-2]["L"] == pytest.approx(converged, abs=5e-5)
    assert 0 < converged / max(r["L"] for r in tw["graded_tube_small_gap"]) - 1 < 0.002


@pytest.mark.parametrize("r", [r for r in DATA["fed"] if r["gap_over_a"] >= 2], ids=lambda r: f"w{r['w']}-g{r['gap_over_a']:g}a")
def test_fed_alike_at_moderate_gaps_the_strip_is_longer_and_as_resistive(r):
    L, R = _fed(r)
    assert 0.065 < (L - r["L_tube"]) / r["w"] < 0.10
    assert R / r["R_tube"] - 1 == pytest.approx(0, abs=0.025)


def test_fed_at_small_gaps_on_wide_strips_the_resistances_part():
    """The tube's ring gap holds more capacitance than the strip's straight one:
    on wide strips at gaps of a or less the strip's resistance is well below the
    tube's, and at w = 0.1, gap a/2, the tube has no resonance left at all while
    the strip still has one."""
    for r in DATA["fed"]:
        if r["w"] >= 0.08 and r["gap_over_a"] <= 1 and r["L_tube"] is not None:
            assert _fed(r)[1] / r["R_tube"] - 1 < -0.03
    last = next(r for r in DATA["fed"] if r["w"] == 0.1 and r["gap_over_a"] == 0.5)
    assert last["L_tube"] is None and 0.45 < _fed(last)[0] < 0.50


@pytest.mark.parametrize("w", sorted({r["w"] for r in DATA["fed"]}))
def test_the_slot_law_lies_inside_the_strips_own_feed_spread(w, registry):
    """What the slot family actually uses - the thin-wire length law at a = w/4 -
    against the strip it stands for, fed across every gap solved."""
    law = registry["resonant_dipole"].synthesize(f0=3e8, aw=w / 4 * 2.99792458e8 / 3e8).metrics["length_over_lambda"]
    Ls = [_fed(r)[0] for r in DATA["fed"] if r["w"] == w]
    assert min(Ls) <= law <= max(Ls)


@pytest.mark.slow
def test_the_rooftop_strip_reproduces_its_coarse_mesh():
    L = st.natural_length(lambda L: st.RooftopStrip(L, 0.04, 60, 6).matrix(), 0.45)
    assert L == pytest.approx(DATA["rooftop_check"]["meshes"][0]["L"], abs=2e-5)


@pytest.mark.slow
def test_the_graded_tube_reproduces_a_fed_row():
    r = next(r for r in DATA["fed"] if r["w"] == 0.04 and r["gap_over_a"] == 2)
    z = st.tube_impedance(r["L_tube"], 0.01, r["gap"])
    assert abs(z.imag) < 0.05 and z.real == pytest.approx(r["R_tube"], rel=1e-4)
