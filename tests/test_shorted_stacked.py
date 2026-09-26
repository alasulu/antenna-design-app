"""The shorted and stacked patches, full-wave - by FDTD, since the spectral
MoM has neither a shorting wall nor a second metal layer.

`otahub.num.patch_fdtd.shorted_ringdown` models a quarter-wave patch with a
full-width shorting wall in a half-space grid (no mirror plane along the
length); the stacked patch fits the quarter-space grid with a second sheet.
Both specs now build lengths corrected with the rectangular patch's full-wave
ratio; these tests pin what the FDTD says that correction does and does not do.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from otahub.num import patch_fdtd as pf

DATA = json.loads((Path(__file__).parent / "data" / "shorted_stacked_fdtd.json").read_text())


def _textbook_quarter_wave(er, hl, Ls, W):
    """Frequency (in f0) at which the textbook formula makes this length resonant."""
    ee = (er + 1) / 2 + ((er - 1) / 2) * (1 + 12 * hl / W) ** -0.5
    dL = 0.412 * hl * ((ee + 0.3) * (W / hl + 0.264)) / ((ee - 0.258) * (W / hl + 0.8))
    return 1 / (4 * math.sqrt(ee) * (Ls + dL))


def test_the_half_space_grid_is_the_quarter_space_grid_on_a_symmetric_patch():
    h = DATA["half_space_check"]
    assert h["half_space"] == pytest.approx(h["quarter_space"], rel=1e-6)


@pytest.mark.parametrize("board", DATA["shorted"], ids=lambda b: f"eps{b['eps_r']}")
def test_the_textbook_shorted_patch_resonates_lower_than_the_full_patch(board):
    runs = board["runs"]
    f = pf.extrapolate([r["nh"] for r in runs], [r["f_over_f0"] for r in runs])
    tb = _textbook_quarter_wave(board["eps_r"], board["h_lam"], board["L_lam"], board["W_lam"])
    shortfall = 1 - f / tb
    assert shortfall > 0.035                                  # 4-8% low
    assert shortfall > 1 - board["full_patch_ratio"] + 0.01    # and more than the full patch


def test_the_spec_corrections_are_the_full_patchs(registry):
    for f0, er, h in ((2.4e9, 4.4, 0.0016), (5.8e9, 10.2, 0.00127), (2.4e9, 2.2, 0.0016)):
        rect = registry["rectangular_patch"].synthesize(f0=f0, eps_r=er, h=h)
        ratio = rect.get("L") / rect.get("L_textbook")
        short = registry["quarter_wave_shorted_patch"].synthesize(f0=f0, eps_r=er, h=h)
        stack = registry["stacked_patch"].synthesize(f0=f0, eps_r=er, h=h)
        assert short.get("L") / short.get("L_textbook") == pytest.approx(ratio, rel=1e-12)
        assert stack.get("L") == pytest.approx(rect.get("L"), rel=1e-12)
        assert stack.get("L2") == pytest.approx(0.95 * stack.get("L"), rel=1e-12)


def test_the_default_stack_is_not_double_tuned():
    """Two coupled modes, but 25-38% apart at a 0.03-wavelength gap whatever the
    size ratio; at the default ratio the upper one is far above f0. Growing the
    parasitic lowers both - the opposite of the old note's 'make it smaller'."""
    scan = sorted(DATA["stacked"]["scan"], key=lambda s: s["size_ratio"])
    for s in scan:
        lo, hi = s["modes"][:2]
        assert hi - lo > 0.2
    assert scan[0]["modes"][1] > 1.3
    lows = [s["modes"][0] for s in scan]
    highs = [s["modes"][1] for s in scan]
    assert lows == sorted(lows, reverse=True) and highs == sorted(highs, reverse=True)


@pytest.mark.slow
def test_the_coarsest_shorted_run_reproduces():
    b = next(b for b in DATA["shorted"] if b["eps_r"] == 10.2)
    r = b["runs"][0]
    cell = b["h_lam"] / r["nh"]
    f, q = pf.shorted_ringdown(b["eps_r"], r["nh"], r["nl"], r["nw"], 0.9 * cell)
    assert f / cell == pytest.approx(r["f_over_f0"], rel=2e-3)
    assert q == pytest.approx(r["Q"], rel=0.02)
