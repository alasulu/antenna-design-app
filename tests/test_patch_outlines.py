"""The circular and triangular patches solved in their own shapes (FDTD), and the size to build.

Both specs sized the patch by the cavity model with a fringing correction -
Balanis's for the disc, a + h/sqrt(eps_r) for the triangle - and called the
resulting resonance an upper bound on thick board. The FDTD with the outline
staircased onto the grid, across eps_r 2.2-10.2 and h sqrt(eps_r)/lambda0 up to
0.095, puts both patches 2-7% low on thick board, the triangle already 5% low
on thin board at eps_r 2.2. The specs now build the cavity size times a fitted
FDTD factor; the cavity physics keeps the cavity size.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from otahub.num import patch_fdtd as pf

C = 2.99792458e8
DATA = json.loads((Path(__file__).parent / "data" / "patch_outline_fdtd.json").read_text())
SPEC = {"circ": ("circular_patch", "a", "a_cavity", "radius_factor"),
        "tri": ("triangular_patch", "a_side", "a_side_cavity", "side_factor")}


def _extrap(pts):
    nh = np.array([p[0] for p in pts], float)
    r = np.array([p[1] for p in pts])
    return float(np.linalg.lstsq(np.stack([np.ones_like(nh), 1 / nh], 1), r, rcond=None)[0][0])


def test_a_rectangle_through_the_outline_engine_is_the_rectangle():
    a = pf.ringdown(2.2, 2, 20, 26, 0.02)
    b = pf.sheet_ringdown(2.2, 2, lambda x, y: (x <= 26) & (y <= 20) & (x >= 0) & (y >= 0), 26, 0.0, 20, 0.02, (0, 18))
    assert a == b


@pytest.mark.parametrize("kind", ["circ", "tri"])
@pytest.mark.parametrize("er,h", [(2.2, 0.0016), (4.4, 0.0032), (10.2, 0.00127)])
def test_the_spec_builds_the_cavity_size_times_the_fdtd_factor(registry, kind, er, h):
    key, built, cavity, factor = SPEC[kind]
    d = registry[key].synthesize(f0=3e9, eps_r=er, h=h)
    x = h * 3e9 * math.sqrt(er) / C
    c = DATA["fit"][kind]
    assert d.get(factor) == pytest.approx(math.exp(c[0] * x + c[1] * x / er + c[2] * x * x + c[3] * x * x / er), rel=1e-8)
    assert d.get(built) == pytest.approx(d.get(cavity) * d.get(factor), rel=1e-12)
    # the cavity model still closes on the cavity size: exactly for the triangle, and to Balanis's own
    # 0.3% for the disc, whose radius formula (14-64) only approximately inverts its resonance (14-66)
    assert d.metrics["resonant_frequency_check_hz"] == pytest.approx(3e9, rel=1e-6 if kind == "tri" else 0.005)


@pytest.mark.parametrize("kind", ["circ", "tri"])
def test_both_resonate_low_on_thick_board(kind):
    thick = [b["r_inf"] for b in DATA["boards"][kind] if len(b["nh"]) >= 3 and b["x"] >= 0.06]
    assert all(0.92 < r < 0.985 for r in thick), thick


def test_the_triangle_is_low_even_on_thin_soft_board():
    b = next(b for b in DATA["boards"]["tri"] if b["eps_r"] == 2.2 and abs(b["x"] - 0.035) < 1e-3)
    assert b["r_inf"] < 0.96


@pytest.mark.parametrize("kind,tol", [("circ", 0.011), ("tri", 0.011)])
def test_the_fit_follows_every_board_it_was_fitted_to(kind, tol):
    for b in DATA["boards"][kind]:
        if len(b["nh"]) >= 3:
            assert b["g_fit"] == pytest.approx(b["g"], rel=tol), b


@pytest.mark.parametrize("kind", ["circ", "tri"])
def test_patches_built_to_the_corrected_size_resonate_at_f0(kind):
    boards = {}
    for r in DATA["verify"]:
        if r["kind"] == kind:
            boards.setdefault(r["eps_r"], []).append((r["nh"], r["f_fdtd"]))
    assert len(boards) == 2
    for er, pts in boards.items():
        pts.sort()
        assert _extrap(pts) == pytest.approx(1.0, abs=0.011), (er, pts)
        assert _extrap([p for p in pts if p[0] >= 6]) == pytest.approx(1.0, abs=0.011), (er, pts)


def test_one_survey_run_live(registry):
    rec = next(r for r in DATA["survey"] if r["kind"] == "circ" and r["eps_r"] == 2.2 and r["nh"] == 4
               and abs(r["h_lam"] - 0.064049) < 1e-6)
    a = registry["circular_patch"].synthesize(f0=C, eps_r=2.2, h=rec["h_lam"]).get("a_cavity")
    cell = rec["h_lam"] / 4
    f, _ = pf.circular_ringdown(2.2, 4, a / cell, 0.95 * cell)
    assert f / cell == pytest.approx(rec["f_fdtd"], abs=1e-6)
