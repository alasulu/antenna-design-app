"""The stacked patch's best air gap and parasitic length, as rules over the driven board.

The spec carried the best band a stack reaches on its board but not the design that
reaches it. tests/data/stacked_patch_rules.json searches the gap and the parasitic's
length on twenty driven boards - eps_r 2.2-10.2, t = h sqrt(eps_r)/lambda 0.0185-0.0385,
sixteen to fit and four held out - with the two-layer spectral MoM and a probe
(otahub.num.stacked). What it found:

- the best stack is one narrow family on every board: an air gap of 0.0675-0.085
  wavelengths under a parasitic 0.36-0.375 wavelengths long, whatever the board;
- the best sits at a cliff. Under it the gap is too tight, the two coupled modes part
  and the band about f0 breaks in two: 0.005 wavelengths less gap leaves under a
  third of the band on 18 of the 20 boards (next to nothing on 12), where 0.005 more
  keeps 83-97% of it;
- so the rule proposed is the fitted best moved back from the cliff - 0.005 more gap,
  0.005 shorter parasitic - and solved there it keeps most of the best band on every
  board, the held-out ones included;
- the best band is 8.2-21 times the textbook single patch's, not a constant 8.4: it
  grows with permittivity, and the four boards surveyed before came out 1-18% short
  of it on the coarser grid they were searched on.

Each expectation here is recomputed from the stored survey; the slow test solves a
held-out board at the rule's design again with the stock solver.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

C = 2.99792458e8
DATA = json.loads((Path(__file__).parent / "data" / "stacked_patch_rules.json").read_text())
BOARDS = DATA["boards"]
FIT = [b for b in BOARDS if b["role"] == "fit"]
HELD = [b for b in BOARDS if b["role"] == "held out"]

# The proposed rules, in t = h f0 sqrt(eps_r)/c and u = 1/eps_r (eps_r 2.2-10.2, t 0.0185-0.0385,
# an air gap): the fitted best design, and the rule - the best moved back from the cliff.
BEST_GAP = (0.060759, -0.15950, 0.17514, -0.29689)            # 1, t, u, u^2
BEST_L2 = (0.37954, -0.14383, 0.28505)                         # 1, u, u^2
RULE_GAP = (BEST_GAP[0] + 0.005,) + BEST_GAP[1:]
RULE_L2 = (BEST_L2[0] - 0.005,) + BEST_L2[1:]
# ln(best band) = c0 + c1 u + c2 u^2 + (c3 + c4 u) ln t + c5 (ln t)^2
BEST_BAND = (-3.7140, -0.74744, -6.6482, -2.0150, -1.9811, -0.51475)

# (fit boards, held out): the fits' largest errors are 0.0032 / 0.0024 in the gap,
# 0.0032 / 0.0037 in the length and 4.2% / 2.5% in the band
GAP_TOL = {False: 0.0035, True: 0.003}
L2_TOL = {False: 0.0035, True: 0.004}
BAND_TOL = {False: 0.045, True: 0.03}


def gap_of(c, t, u):
    return c[0] + c[1] * t + c[2] * u + c[3] * u * u


def l2_of(c, u):
    return c[0] + c[1] * u + c[2] * u * u


def band_of(t, u, c=BEST_BAND):
    lt = math.log(t)
    return math.exp(c[0] + c[1] * u + c[2] * u * u + (c[3] + c[4] * u) * lt + c[5] * lt * lt)


def _lsq(cols, y):
    return np.linalg.lstsq(np.stack(cols, 1), np.asarray(y), rcond=None)[0]


def test_the_best_stack_is_one_narrow_family():
    gaps = [b["best"]["gap"] for b in BOARDS]
    l2s = [b["best"]["l2"] for b in BOARDS]
    assert 0.0675 - 1e-9 <= min(gaps) and max(gaps) <= 0.085 + 1e-9
    assert 0.36 - 1e-9 <= min(l2s) and max(l2s) <= 0.375 + 1e-9
    # as a ratio to the driven patch it is anything but narrow: 1.15 on soft board, 2.4 on eps_r 10.2
    ratios = [b["best"]["ratio"] for b in BOARDS]
    assert min(ratios) < 1.16 and max(ratios) > 2.4


def test_the_rules_are_the_least_squares_fit_of_the_sixteen():
    t = np.array([b["t"] for b in FIT])
    u = 1 / np.array([b["eps_r"] for b in FIT])
    one = np.ones_like(t)
    cg = _lsq([one, t, u, u * u], [b["best"]["gap"] for b in FIT])
    cl = _lsq([one, u, u * u], [b["best"]["l2"] for b in FIT])
    lt = np.log(t)
    cb = _lsq([one, u, u * u, lt, u * lt, lt * lt], np.log([b["best"]["band"] for b in FIT]))
    assert cg == pytest.approx(BEST_GAP, rel=1e-4)
    assert cl == pytest.approx(BEST_L2, rel=1e-4)
    assert cb == pytest.approx(BEST_BAND, rel=1e-4)
    assert DATA["rules"]["best_gap"]["coef"] == pytest.approx(BEST_GAP, rel=1e-4)
    assert DATA["rules"]["best_l2"]["coef"] == pytest.approx(BEST_L2, rel=1e-4)


@pytest.mark.parametrize("b", BOARDS, ids=[b["name"] for b in BOARDS])
def test_the_fits_follow_every_board(b):
    """The search resolves the gap to 0.0025 and the length to 0.005 wavelengths."""
    t, u = b["t"], 1 / b["eps_r"]
    held = b["role"] == "held out"
    assert gap_of(BEST_GAP, t, u) == pytest.approx(b["best"]["gap"], abs=GAP_TOL[held])
    assert l2_of(BEST_L2, u) == pytest.approx(b["best"]["l2"], abs=L2_TOL[held])
    assert band_of(t, u) == pytest.approx(b["best"]["band"], rel=BAND_TOL[held])


def _band_at(curve_gaps, bands, gap):
    i = int(np.argmin(np.abs(np.asarray(curve_gaps) - gap)))
    assert abs(curve_gaps[i] - gap) < 1e-6
    return bands[i]


def test_the_best_sits_at_a_cliff():
    """At the best parasitic length: 0.005 wavelengths less gap, and 0.005 more."""
    gaps = DATA["settings"]["gaps"]
    under, over = [], []
    for b in BOARDS:
        curve = b["curves"][f"{b['best']['l2']:.4f}"]
        under.append(_band_at(gaps, curve, b["best"]["gap"] - 0.005) / b["best"]["band"])
        over.append(_band_at(gaps, curve, b["best"]["gap"] + 0.005) / b["best"]["band"])
    assert sum(r < 0.34 for r in under) >= 18
    assert sum(r < 0.02 for r in under) >= 12
    assert min(over) > 0.82 and max(over) < 0.98


def test_the_rule_keeps_the_band_on_every_board():
    """Both bands read off the solved tables - the best as the largest entry of every
    curve, the rule's as its own curve at its gap - and the summaries checked against
    them, so an edited summary cannot pass for retention."""
    kept = {}
    for b in BOARDS:
        best = max(max(curve) for curve in b["curves"].values())
        rule = _band_at(b["rule"]["gaps"], b["rule"]["bands"], b["rule"]["gap"])
        assert b["best"]["band"] == pytest.approx(best, abs=1e-9)
        assert b["rule"]["band"] == pytest.approx(rule, abs=1e-9)
        kept[b["name"]] = rule / best
    for b in BOARDS:
        t, u = b["t"], 1 / b["eps_r"]
        assert b["rule"]["gap"] == pytest.approx(gap_of(RULE_GAP, t, u), abs=1e-4)
        assert b["rule"]["l2"] == pytest.approx(l2_of(RULE_L2, u), abs=1e-5)
    fit = [kept[b["name"]] for b in FIT]
    held = [kept[b["name"]] for b in HELD]
    assert min(fit) > 0.90 and np.mean(fit) > 0.95                  # 0.91-0.99, the least on eps_r 10.2
    assert min(held) > 0.94 and np.mean(held) > 0.95                # 0.95-0.98


def test_the_rule_tolerates_its_own_error_on_the_held_out_boards():
    """0.005 wavelengths either way in the parasitic and in the gap about the rule's
    design - more than the fits' own error. Three held-out boards keep 80% of the best
    band over the whole box; on eps_r 6.15 the cliff is nearer, and the corner with both
    the gap and the parasitic short keeps 37%."""
    worst = {}
    for b in HELD:
        w = b["rule"]["band"]
        g0 = b["rule"]["gap"]
        for side in b["rule"]["l2_either_side"].values():
            for dg in (-0.005, 0.0, 0.005):
                w = min(w, _band_at(side["gaps"], side["bands"], g0 + dg))
        for dg in (-0.005, 0.005):
            w = min(w, _band_at(b["rule"]["gaps"], b["rule"]["bands"], g0 + dg))
        worst[b["eps_r"]] = w / b["best"]["band"]
    assert min(worst[e] for e in (2.5, 3.55, 4.0)) > 0.80
    assert 0.35 < worst[6.15] < 0.40
    # the gap alone 0.005 high, at the rule's parasitic: 87-91% on all four
    for b in HELD:
        hi = _band_at(b["rule"]["gaps"], b["rule"]["bands"], b["rule"]["gap"] + 0.005)
        assert 0.86 < hi / b["best"]["band"] < 0.92


def test_the_band_is_not_a_constant_multiple_of_the_single_patch():
    """The textbook single-patch band (the spec's single_patch_bandwidth_vswr2) against
    the best band solved: 8.2 times on thick eps_r 2.2, 21 on thick eps_r 10.2."""
    ratios = {}
    for b in BOARDS:
        single = 3.771 * (b["eps_r"] - 1) / b["eps_r"] ** 2 * b["h_lambda"] * b["W"] / b["L"]
        ratios[(b["eps_r"], b["t"])] = b["best"]["band"] / single
    assert min(ratios.values()) == pytest.approx(ratios[(2.2, 0.0385)], rel=1e-12)
    assert 8.0 < min(ratios.values()) < 8.4 and 21.0 < max(ratios.values()) < 21.5
    for t in (0.0185, 0.0252, 0.0318, 0.0385):
        r = [ratios[(e, t)] for e in (2.2, 3.0, 4.6, 10.2)]
        assert r[-1] > 1.5 * r[0]


def test_the_four_boards_surveyed_before_fell_short_of_the_best():
    """tests/data/stacked_patch_design.json stepped the gap by 0.02 and the ratio by about
    0.05-0.09: its best on each board is below the band rule's, by 1-18% - most on the
    default board, whose best it put at a 0.09 gap and ratio 1.1 (11.5%). On the board
    surveyed here next to it (eps_r 2.2, t 0.0185) the best is 13.7% at a 0.075 gap and
    a parasitic 0.37 long (ratio 1.15), and the rule's design keeps 13.0%."""
    old = json.loads((Path(__file__).parent / "data" / "stacked_patch_design.json").read_text())
    short = []
    for b in old["boards"]:
        t = b["h_lambda"] * math.sqrt(b["eps_r"])
        best = max(d["band"] for d in b["designs"])
        short.append(best / band_of(t, 1 / b["eps_r"]))
    assert max(short) < 1.0 and min(short) > 0.80


def test_the_spec_directivity_fit_holds_at_the_new_best_designs(registry):
    """The spec's directivity fit (made on the earlier survey's designs) against the
    directivity solved at each best design here, where its domain covers it (eps_r up to
    4.41); beyond, the solved values are 9.19-9.41 dBi."""
    lam = C / DATA["settings"]["f0_hz"]
    n = 0
    for b in BOARDS:
        d = registry["stacked_patch"].synthesize(f0=DATA["settings"]["f0_hz"], eps_r=b["eps_r"],
                                                 h=b["t"] * lam / math.sqrt(b["eps_r"]), eps_r2=1.0,
                                                 h2_over_lambda=b["best"]["gap"], size_ratio=b["best"]["ratio"])
        spec = d.metrics["directivity_dbi"]
        if math.isfinite(spec):
            n += 1
            assert spec == pytest.approx(b["directivity_dbi"], abs=0.05), b["name"]
        else:
            assert 9.15 < b["directivity_dbi"] < 9.45
    assert n >= 9


@pytest.mark.slow
def test_a_held_out_board_at_the_rules_design_solves_again():
    """eps_r 4.0, t 0.0285 at the rule's gap and parasitic, with the stock solver
    (patch_sdm.StackedPatch and probe_vector, splined by stacked.impedance_sweep)."""
    from otahub.core.registry import default_registry
    from otahub.num import patch_sdm as sdm
    from otahub.num import stacked as st
    b = next(x for x in HELD if x["eps_r"] == 4.0)
    s = DATA["settings"]
    f0 = s["f0_hz"]
    lam = C / f0
    h = b["t"] * lam / math.sqrt(b["eps_r"])
    d = default_registry()["stacked_patch"].synthesize(f0=f0, eps_r=b["eps_r"], h=h)
    L, W = d.get("L"), d.get("W")
    assert L / lam == pytest.approx(b["L"], rel=1e-9)
    L2 = b["rule"]["l2"] * lam
    p = sdm.StackedPatch(b["eps_r"], h, L, W, 1.0, b["rule"]["gap"] * lam, L2, W * L2 / L, **sdm.BASIS_FULL)
    rs = np.array(s["rs"])
    rf = np.arange(rs[0], rs[-1] + 1e-9, s["rf_step"])
    zi = st.impedance_sweep(p, f0, rs, np.array(s["feeds"]) * L / 2, s["probe_radius_m"], rf)
    band = max(st.vswr_band(rf, zi[:, q], 20.0)[0] for q in range(len(s["feeds"])))
    assert band == pytest.approx(b["rule"]["band"], abs=1e-3)
