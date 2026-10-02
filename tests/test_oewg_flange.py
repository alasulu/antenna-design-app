"""The flanged open-ended waveguide over its band: what the guide adds to the TE10 aperture.

The spec's directivity is the TE10 aperture - cos(pi x/a) across a, uniform across
b - in an infinite ground plane. The guide itself, in an infinite flange, is solved
here by FDTD (`horn_fdtd.sectoral_horn(..., flange=True)`: the far field of the
aperture plane's own E, M = -z x E imaged through the flange, which is exact for a
flange, so only the aperture field is approximated). It reads above the TE10
aperture everywhere in the single-mode band, from 0.06 dB near cutoff to 0.31 dB: the
r^-1/3 field at the broad walls' edges and the evanescent modes the flange excites,
which a pure TE10 aperture leaves out.

The survey (tests/data/oewg_flange_fdtd.json) runs a/lambda from cutoff to 1 at
b/a 0.3, 4/9 (the WR family), 1/2 and 0.6, every guide a whole number of cells. The
box, the absorber, the ring-down and the guide length move the answer by under
0.0002 dB; the grid does not - 36 cells across a reach 0.8 of the excess and 72
cells 0.9, converging at an order of about 1.2 (the edges' singular field) - so the
excess is fitted together with a grid term, and its continuum part F is what the
spec would carry. (The single WR-90 record in fdtd_farfield.json read 0.17 dB "on
three grids" because its b/a changed with the grid and cancelled the refinement; at
a fixed guide it is 0.21.)

The TE10 aperture the runs are compared with is integrated here independently of the
spec's fit (test_open_ended_waveguide._direct).
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from otahub.core.expr import evaluate
from otahub.num import horn_fdtd as hf
from tests.test_open_ended_waveguide import _direct

DATA = json.loads((Path(__file__).parent / "data" / "oewg_flange_fdtd.json").read_text())
OLD = json.loads((Path(__file__).parent / "data" / "fdtd_farfield.json").read_text())["oewg"]
RUNS = DATA["runs"]
FIT = DATA["fit"]
ASPECT = {"A": 4 / 9, "B": 1 / 2, "C": 0.3, "D": 0.6, "E": 0.4, "F": 0.55, "G": 0.45}
C0 = 2.99792458e8


def _te10_dbi(a, b):
    return 10 * math.log10(_direct(a, b))


def _delta(r):
    return r["fdtd_dbi"] - r["te10_dbi"]


def _design(a, b, N):
    """Columns of F = sum c_ij x^i y^j (x = a/lambda, y = b/lambda) and of the grid
    term h^order (g0 + g1 x + g2 y/x), h = 1/N; N = inf gives the continuum."""
    a, b, N = (np.atleast_1d(np.asarray(v, float)) for v in (a, b, N))
    h = (1.0 / N) ** FIT["order"]
    return np.stack([a ** i * b ** j for i, j in FIT["F_terms"]] + [h, h * a, h * b / a], 1)


def _columns(rows):
    return _design([r["a"] for r in rows], [r["b"] for r in rows], [r["N"] for r in rows])


COEF = np.array(FIT["F"] + FIT["G"])


def _fitted():
    return [r for r in RUNS if r["kind"] in ("survey", "grid")]


def _refit(rows):
    return np.linalg.lstsq(_columns(rows), np.array([_delta(r) for r in rows]), rcond=None)[0]


def excess_db(a, b):
    """The continuum excess F(a/lambda, b/lambda) over the TE10 aperture."""
    return float(_design(a, b, math.inf)[0] @ COEF)


# ------------------------------------------------------------------ the records

def test_every_record_is_the_guide_it_names():
    tags = [r["tag"] for r in RUNS]
    assert len(tags) == len(set(tags))
    for r in RUNS:
        # a whole number of cells across each half: a = 2p/N, b = 2q/N wavelengths
        assert all(isinstance(r[k], int) for k in ("p", "q", "N"))
        assert r["a"] == pytest.approx(2 * r["p"] / r["N"], rel=1e-12)
        assert r["b"] == pytest.approx(2 * r["q"] / r["N"], rel=1e-12)
        assert r["q"] / r["p"] == pytest.approx(ASPECT[r["family"]], rel=1e-12), r["tag"]
        assert 0.5 < r["a"] <= 1.0 + 1e-12                    # TE10 above cutoff, TE20 not
        assert r["kind"] in ("survey", "holdout", "grid", "coarse", "margin", "npml", "periods", "guide")
        assert r["tag"].startswith(f"{r['family']}-p{r['p']}q{r['q']}-N{r['N']}")
        if r["kind"] in ("survey", "holdout", "grid", "coarse") and not r["tag"].endswith("-cut"):
            assert all(r[k] == v for k, v in DATA["base_settings"].items()), r["tag"]
    assert all(r["p"] >= 18 for r in _fitted()) and all(r["p"] < 18 for r in RUNS if r["kind"] == "coarse")


def test_the_reference_is_the_integrated_te10_aperture():
    """Recomputed here for every guide, independently of the spec's fit."""
    seen: dict = {}
    for r in RUNS:
        key = (r["a"], r["b"])
        if key not in seen:
            seen[key] = _te10_dbi(*key)
        assert r["te10_dbi"] == pytest.approx(seen[key], abs=1e-6), r["tag"]


def test_the_old_wr90_record_is_one_of_the_survey_runs():
    """Same guide (0.75 x 1/3 wavelengths at 48 cells), same settings: the same number."""
    old = next(r for r in OLD if r["N"] == 48)
    new = next(r for r in RUNS if r["tag"] == "A-p18q8-N48")
    assert (new["a"], new["b"]) == pytest.approx((old["a"], old["b"]), rel=1e-12)
    assert new["fdtd_dbi"] == pytest.approx(old["fdtd_dbi"], abs=1e-9)


def test_box_absorber_ring_down_and_guide_length_change_nothing():
    """Twice the margin, 20 absorber cells, twice the periods, a 2.5-wavelength guide."""
    base = {(r["a"], r["b"], r["N"]): r for r in RUNS
            if r["kind"] in ("survey", "grid", "coarse") and not r["tag"].endswith("-cut")}
    variants = [r for r in RUNS if r["kind"] in ("margin", "npml", "periods", "guide")]
    assert {r["kind"] for r in variants} == {"margin", "npml", "periods", "guide"}
    for r in variants:
        assert r["fdtd_dbi"] == pytest.approx(base[(r["a"], r["b"], r["N"])]["fdtd_dbi"], abs=2e-4), r["tag"]


def test_the_grid_converges_at_an_order_above_one():
    """Five or six grids of one guide (18 to 144 cells across a): the excess grows with
    the cells, Delta(h) = Delta_inf + C h^q fits best with q between 1 and 4/3, and the
    finest ratio-2 triple (36-72-144 cells across a) gives about 1.2."""
    groups: dict = {}
    for r in RUNS:
        if r["kind"] in ("survey", "grid", "coarse") and not r["tag"].endswith("-cut"):
            groups.setdefault((r["a"], r["b"]), []).append(r)
    deep = [sorted(g, key=lambda r: r["N"]) for g in groups.values() if len(g) >= 5]
    assert len(deep) == 3
    qs = np.linspace(0.5, 2.0, 151)
    for g in deep:
        d = np.array([_delta(r) for r in g])
        h = np.array([1 / r["N"] for r in g])
        assert np.all(np.diff(d) > 0)
        rms = [np.linalg.lstsq(np.stack([np.ones_like(h), h ** q], 1), d, rcond=None)[1][0] for q in qs]
        assert 1.0 <= qs[int(np.argmin(rms))] <= 4 / 3, g[0]["tag"]
    top = {r["N"]: _delta(r) for r in groups[(1.0, 4 / 9)]}
    assert math.log2((top[72] - top[36]) / (top[144] - top[72])) == pytest.approx(FIT["order"], abs=0.05)


# ------------------------------------------------------------------ the fit

def test_the_fit_is_the_least_squares_fit_of_the_runs():
    rows = _fitted()
    c = _refit(rows)
    assert np.max(np.abs(_columns(rows) @ c - _columns(rows) @ COEF)) < 1e-7
    res = _columns(rows) @ COEF - np.array([_delta(r) for r in rows])
    assert np.max(np.abs(res)) == pytest.approx(FIT["train_max_db"], abs=1e-6)
    assert FIT["train_max_db"] < 0.002


def test_the_fit_holds_on_guides_it_never_saw():
    """b/a 0.4, 0.45 and 0.55 - between the surveyed aspects - and every third WR-family
    point left out of a refit."""
    hold = [r for r in RUNS if r["kind"] == "holdout"]
    assert {r["family"] for r in hold} == {"E", "F", "G"}
    err = np.max(np.abs(_columns(hold) @ COEF - [_delta(r) for r in hold]))
    assert err == pytest.approx(FIT["holdout_max_db"], abs=1e-6) and err < 0.0005
    left = [r for r in _fitted() if r["family"] in "AB" and r["kind"] == "survey" and r["N"] % 3 == 1]
    c2 = _refit([r for r in _fitted() if r not in left])
    assert np.max(np.abs(_columns(left) @ c2 - [_delta(r) for r in left])) < 0.001


def test_the_fit_reproduces_the_older_records_it_was_not_fitted_to():
    """fdtd_farfield.json's three WR-90 runs: b/a 0.467, 0.444 and 0.429 on 30, 36 and 42
    cells across a - read 0.17 dB above the TE10 aperture together only by cancellation."""
    for r in OLD:
        rr = dict(r, te10_dbi=_te10_dbi(r["a"], r["b"]))
        assert float(_columns([rr])[0] @ COEF) == pytest.approx(_delta(rr), abs=0.001), r["N"]
        assert excess_db(r["a"], r["b"]) - _delta(rr) > 0.025                 # the grid's share


def test_the_excess_is_positive_and_grows_from_cutoff():
    x = np.linspace(0.5, 1.0, 51)
    for ratio in (0.3, 4 / 9, 0.5, 0.6):
        f = np.array([excess_db(a, ratio * a) for a in x])
        assert np.all(f > 0.06) and np.all(f < 0.32)
        assert np.all(np.diff(f[:26]) > 0)                        # rising over the lower half of the band
    # a taller guide gains more: more of its aperture is near the broad walls' edges
    for a in x:
        assert excess_db(a, 0.3 * a) < excess_db(a, 0.45 * a) < excess_db(a, 0.6 * a)
    assert excess_db(0.7625, 0.3389) == pytest.approx(0.21, abs=0.005)          # WR-90 at 10 GHz


def test_mode_matching_converges_toward_the_same_excess():
    """Indicative - recorded from an independent mode-matching model of the same flanged
    aperture: it rises toward F with the modes, and at first order in 1/kc lands on it."""
    for m in FIT["modal_check"]["runs"]:
        f = excess_db(m["a"], m["b"])
        assert m["delta_kc16"] < m["delta_kc24"] < m["delta_kc32"] < f
        assert m["extrapolated_db"] == pytest.approx(f, abs=0.003)


# ------------------------------------------------------------------ the proposal

def _proposed(f0, a_wg, b_wg):
    return float(evaluate(DATA["proposal"]["directivity_expr"], {"f0": f0, "a_wg": a_wg, "b_wg": b_wg}))


def test_the_proposed_expression_is_the_te10_aperture_plus_the_excess():
    """Written in the spec's variables; its TE10 part is the spec's present fit (0.04%)."""
    lam = C0 / 10e9
    for a in np.linspace(0.5, 1.0, 11):
        for ratio in (0.3, 0.4, 4 / 9, 0.5, 0.6):
            b = ratio * a
            assert _proposed(10e9, a * lam, b * lam) == pytest.approx(_te10_dbi(a, b) + excess_db(a, b), abs=0.002), (a, ratio)
            ex = float(evaluate(DATA["proposal"]["excess_expr"], {"f0": 10e9, "a_wg": a * lam, "b_wg": b * lam}))
            assert ex == pytest.approx(excess_db(a, b), abs=1e-7)
    assert math.isnan(_proposed(15e9, 0.02286, 0.01016))         # TE20 propagates: outside, as now
    assert math.isnan(_proposed(10e9, 0.02286, 0.0160))          # b/a 0.7


def _interpolated_excess(a, ratio):
    """A known case's excess from the runs themselves: each survey run of the WR-family
    aspects (b/a 4/9 and 1/2) less its fitted grid term, interpolated linearly in a/lambda,
    then between the two aspects. It leans on the fit's grid term G, not on its polynomial F."""
    G, q = FIT["G"], FIT["order"]
    vals = []
    for fam in ("A", "B"):
        rs = sorted((r for r in RUNS if r.get("family") == fam and r.get("kind") == "survey"), key=lambda r: r["a"])
        xs = [r["a"] for r in rs]
        ys = [_delta(r) - (1 / r["N"]) ** q * (G[0] + G[1] * r["a"] + G[2] * r["b"] / r["a"]) for r in rs]
        vals.append(float(np.interp(a, xs, ys)))
    t = (ratio - 4 / 9) / (1 / 2 - 4 / 9)
    return vals[0] + (vals[1] - vals[0]) * t


@pytest.mark.parametrize("case", DATA["proposal"]["known_cases"], ids=lambda k: k["name"])
def test_the_proposed_known_cases(case):
    """Each value is the TE10 integral at the exact guide plus the runs' own excess -
    the bracketing WR-family runs, each less its grid term, interpolated in a/lambda (and
    between b/a 4/9 and 1/2) - not the fitted polynomial."""
    g = case["given"]
    a, b = g["a_wg"] * g["f0"] / C0, g["b_wg"] * g["f0"] / C0
    assert case["te10_dbi"] == pytest.approx(_te10_dbi(a, b), abs=1e-5)
    assert case["excess_db"] == pytest.approx(_interpolated_excess(a, b / a), abs=2e-6)   # recomputed, not trusted
    D = case["expect"]["directivity_dbi"]
    assert D == pytest.approx(case["te10_dbi"] + case["excess_db"], abs=1e-4)
    assert case["expect"]["effective_to_physical_area"] == pytest.approx(10 ** (D / 10) / (4 * math.pi * a * b), rel=5e-5)
    assert _proposed(g["f0"], g["a_wg"], g["b_wg"]) == pytest.approx(D, abs=0.002)
    assert case["tol_pct"] == 0.5


# ------------------------------------------------------------------ live

@pytest.mark.slow
def test_a_live_run_reproduces_its_record():
    r = next(r for r in RUNS if r["tag"] == DATA["live_tag"])
    got = hf.sectoral_horn(r["a"], r["b"], r["a"], r["b"], 0.0, cells_per_lambda=r["N"], guide_len=r["guide_len"],
                           margin=r["margin"], npml=r["npml"], periods=r["periods"], flange=True)
    assert 10 * math.log10(got["directivity"]) == pytest.approx(r["fdtd_dbi"], abs=1e-4)
