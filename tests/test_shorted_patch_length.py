"""The quarter-wave shorted patch's length, solved with its shorting wall.

The spec built the textbook length times the full patch's full-wave ratio, which
corrects the open edge but has no term for the wall's own inductance, and quoted
the remainder as 1-3% from three thin boards. A half-space FDTD survey of the
spec's own patch over eleven boards (eps_r 2.2-10.2, h sqrt(eps_r)/lambda0
0.035-0.095), each extrapolated to zero cell size, found it resonating 1.2-11%
low - most on thick low-permittivity board, since the wall's inductance grows
with its height. The spec now multiplies in the wall's own factor, fitted to
vanish on thin board, and the patch built to it resonates at f0 on held-out boards.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import brentq

from otahub.num import patch_fdtd

C = 2.99792458e8
DATA = json.loads((Path(__file__).parent / "data" / "shorted_patch_survey.json").read_text())
HELD = {tuple(h) for h in DATA["held_out"]}


def _g(er, x):
    c = DATA["fit"]
    return math.exp(c[0] * x + c[1] * x / er + c[2] * x * x + c[3] * x * x / er)


def _extrapolate(pts):
    nh = np.array([p[0] for p in pts], float)
    r = np.array([p[1] for p in pts])
    A = np.stack([np.ones_like(nh), 1 / nh], 1)
    return float(np.linalg.lstsq(A, r, rcond=None)[0][0])


@pytest.mark.parametrize("f0,er,h", [(2.4e9, 4.4, 0.0016), (3e9, 2.2, 0.0064), (5.8e9, 10.2, 0.00127), (1e9, 3.0, 0.0066)])
def test_the_spec_builds_the_open_edge_length_times_the_wall_factor(registry, f0, er, h):
    d = registry["quarter_wave_shorted_patch"].synthesize(f0=f0, eps_r=er, h=h)
    x = h * f0 * math.sqrt(er) / C
    assert d.get("wall_length_factor") == pytest.approx(_g(er, x), rel=1e-9)
    assert d.get("L") == pytest.approx(d.get("L_open_edge_only") * d.get("wall_length_factor"), rel=1e-12)


def test_the_former_length_resonated_low_most_on_thick_soft_board():
    r = {(b["eps_r"], round(b["x"], 3)): b["r_inf"] for b in DATA["boards"] if not b["prev"]}
    assert all(v < 0.99 for v in r.values())
    assert r[(2.2, 0.095)] < 0.90                        # 11% low where the spec said 1-3%
    assert r[(2.2, 0.095)] < r[(4.4, 0.095)] < r[(10.2, 0.095)]
    for er in (2.2, 4.4, 10.2):
        assert r[(er, 0.035)] > r[(er, 0.065)] > r[(er, 0.095)]


def test_the_previous_spot_checks_come_out_as_recorded():
    """The earlier round's three boards, through the same inverse: 2.4, 1.4 and 2.7% low."""
    prev = {b["eps_r"]: b["r_inf"] for b in DATA["boards"] if b["prev"]}
    assert prev[2.2] == pytest.approx(0.9763, abs=5e-4)
    assert prev[4.4] == pytest.approx(0.9859, abs=5e-4)
    assert prev[10.2] == pytest.approx(0.9734, abs=5e-4)


def test_the_fit_follows_every_board():
    for b in DATA["boards"]:
        tol = 0.004 if (b["eps_r"], round(b["x"], 3)) not in HELD else 0.003
        assert _g(b["eps_r"], b["x"]) == pytest.approx(b["g"], rel=tol), b


def test_the_extrapolation_is_stable():
    """Least squares in 1/nh through every grid, and through 6 cells and finer only:
    within 0.6% on every board with three or more grids."""
    for b in DATA["boards"]:
        if len(b["nh"]) < 3:
            continue
        pts = [(n, r) for n, r in zip(b["nh"], b["ratios"]) if n >= 6]
        assert abs(_extrapolate(pts) - b["r_inf"]) < 0.006, b


def test_held_out_boards_resonate_at_f0():
    boards = {}
    for r in DATA["verify"]:
        boards.setdefault((r["eps_r"], round(r["h_lam"] * math.sqrt(r["eps_r"]), 3)), {})[r["nh"]] = r["ratio"]
    assert len(boards) == 3
    for key, grids in boards.items():
        r_inf = _extrapolate(sorted(grids.items()))
        assert r_inf == pytest.approx(1.0, abs=0.013), (key, r_inf)
    assert _extrapolate(sorted(boards[(3.0, 0.05)].items())) == pytest.approx(1.0, abs=0.003)
    assert _extrapolate(sorted(boards[(6.15, 0.08)].items())) == pytest.approx(1.0, abs=0.003)


def test_one_survey_run_live(registry):
    """The FDTD is deterministic: the thickest eps_r 2.2 board at 4 cells, rebuilt from
    the former length and solved again, lands on the recorded ratio."""
    rec = next(r for r in DATA["survey"] if r["eps_r"] == 2.2 and r["nh"] == 4 and abs(r["h_lam"] - 0.064049) < 1e-6)

    def former_L(f):
        return registry["quarter_wave_shorted_patch"].synthesize(f0=f, eps_r=2.2, h=rec["h_lam"]).get("L_open_edge_only")

    cell = rec["h_lam"] / 4
    nl, nw = int(round(former_L(C) / cell)), int(round(registry["quarter_wave_shorted_patch"]
                                                      .synthesize(f0=C, eps_r=2.2, h=rec["h_lam"]).get("W") / 2 / cell))
    assert (nl, nw) == (rec["nl"], rec["nw"])
    f_spec = brentq(lambda fr: former_L(fr * C) - nl * cell, 0.8, 0.0999 / (rec["h_lam"] * math.sqrt(2.2)), xtol=1e-10)
    f, _ = patch_fdtd.shorted_ringdown(2.2, 4, nl, nw, 0.9 * cell)
    assert (f / cell) / f_spec == pytest.approx(rec["ratio"], abs=1e-5)


def test_a_narrower_pifa_strip_makes_a_shorter_plate_and_none_outside_the_survey(registry):
    """Codex: a 1 mm strip got the full-width length (24.29 mm at 2.4 GHz). A narrow
    corner strip lengthens the current path, so the plate must shrink; where the
    shrunken plate would be wider than 1.6 times its length, or the strip narrower
    than W/30, the survey says nothing and neither does the spec."""
    pifa = registry["pifa"]
    full = pifa.synthesize(f0=2.4e9, h=0.006, W=0.02, Ws=0.02).get("L")
    assert full == pytest.approx(0.2425 * C / 2.4e9 - 0.006, rel=1e-12)
    lengths = [pifa.synthesize(f0=2.4e9, h=0.006, W=0.02, Ws=ws).get("L") for ws in (0.02, 0.012, 0.008, 0.005)]
    assert all(a > b for a, b in zip(lengths, lengths[1:]))
    assert math.isnan(pifa.synthesize(f0=2.4e9, h=0.006, W=0.02, Ws=0.001).get("L"))      # W/L would be 2.7
    assert math.isnan(pifa.synthesize(f0=1e9, h=0.01, W=0.04, Ws=0.001).get("L"))        # strip W/40


PIFA = json.loads((Path(__file__).parent / "data" / "pifa_strip_fdtd.json").read_text())


def _pifa_points():
    """Each strip run as a design in wavelengths, corrected by its plate's full-width run."""
    by = {}
    for r in PIFA["runs"]:
        by.setdefault(r["tag"], []).append(r)
    for tag, rs in by.items():
        full = next(r for r in rs if r["ns"] == r["nw"])
        eps = full["f"] * (full["nl"] + full["nh"]) / 0.2425 - 1
        if abs(eps) > 0.1:                      # W/L above 1.6: another mode, not used
            continue
        for r in rs:
            if r["ns"] < r["nw"]:
                lam = (1 + eps) / r["f"]
                yield r["nh"] / lam, r["nw"] / lam, r["nl"] / lam, r["ns"] / r["nw"]


def test_the_pifa_strip_law_is_the_whole_plate_survey(registry):
    """70 designs on 21 plates, normalised by each plate's full-width run. The
    synthesised length is checked, and the frequency at which the spec would build
    each surveyed plate is solved from the law: within 2.8% of where the FDTD put it
    (relative to that normalisation - the spec's notes give the absolute caveat)."""
    from scipy.optimize import brentq
    pifa = registry["pifa"]
    pts = list(_pifa_points())
    assert len(pts) == 70
    worst = 0.0
    for h, W, L, r in pts:                       # lengths in wavelengths at the corrected resonance
        def gap(lam):                            # the spec's L at wavelength lam, minus the plate's
            return pifa.synthesize(f0=C / lam, h=h, W=W, Ws=r * W).get("L") - L
        assert abs(gap(1.0)) / L < 0.06
        lams = np.linspace(0.93, 1.07, 29)       # a bracket inside the law's domain guard
        g = [gap(x) for x in lams]
        k = next(i for i in range(len(lams) - 1)
                 if math.isfinite(g[i]) and math.isfinite(g[i + 1]) and g[i] * g[i + 1] <= 0)
        worst = max(worst, abs(1 / brentq(gap, lams[k], lams[k + 1], xtol=1e-12) - 1))
    assert worst < 0.029


@pytest.mark.slow
def test_a_held_out_narrow_strip_pifa_resonates_where_it_was_designed(registry):
    """Not one of the survey's plates: h 0.04, W 0.12 and a W/6 strip, synthesised at
    f0, built in the whole-plate FDTD at four cells across h and corrected by its own
    full-width run. It lands within 0.1% (2% on a second design, half of that the
    plate length rounded to whole cells)."""
    lam = C / 1e9
    L = registry["pifa"].synthesize(f0=1e9, h=0.04 * lam, W=0.12 * lam, Ws=0.02 * lam).get("L") / lam
    nh, cell = 4, 0.01
    nl, nw, ns = round(L / cell), 12, 2
    f_full, _ = patch_fdtd.pifa_ringdown(nh, nl, nw, nw, 0.2425 / (nl + nh))
    f_strip, _ = patch_fdtd.pifa_ringdown(nh, nl, nw, ns, cell)
    eps = f_full * (nl + nh) / 0.2425 - 1
    assert cell * (1 + eps) / f_strip == pytest.approx(1.0, abs=0.02)


@pytest.mark.slow
def test_the_whole_plate_fdtd_is_the_mirrored_one_for_a_full_width_short():
    """pifa_ringdown models the whole plate with no magnetic wall, so it can take a
    strip at one corner; with the strip as wide as the plate it must reproduce the
    half model behind the 0.2425 calibration."""
    nh, nl, half = 3, 18, 5
    guess = 0.97 / (4 * (nl + nh))
    f_half, q_half = patch_fdtd.shorted_ringdown(1.0, nh, nl, half, guess)
    f_whole, q_whole = patch_fdtd.pifa_ringdown(nh, nl, 2 * half, 2 * half, guess)
    assert f_whole == pytest.approx(f_half, rel=1e-6) and q_whole == pytest.approx(q_half, rel=1e-4)
