"""The stacked-patch design tool: `otahub.num.stacked` and `otahub stack`.

The band finder is checked on circuits whose VSWR-2 band is known exactly; the
survey behind the stacked patch's spec is in tests/data/stacked_patch_design.json,
spot-checked live against the spectral MoM, and the command is run on a narrow
sweep.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from otahub.num import stacked as st

RF = np.arange(0.8, 1.2 + 1e-9, 0.0005)


@pytest.mark.parametrize("q", [5.0, 10.0, 30.0])
def test_a_matched_series_resonator_has_the_textbook_band(q):
    """R = Z0 in series with L and C: VSWR <= 2 while |X| <= R/sqrt(2), a fractional
    band of 1/(sqrt(2) Q) to first order - exactly, in f/f0 - 1/(f/f0) = 1/(sqrt(2) Q)."""
    z = 50.0 + 1j * 50.0 * q * (RF - 1 / RF)
    u = 1 / (2 * math.sqrt(2) * q)                        # f/f0 - f0/f = +-2u at the band's edges
    hi, lo = u + math.sqrt(u * u + 1), -u + math.sqrt(u * u + 1)
    bw, centre, _ = st.vswr_band(RF, z, x_nonres=0.0, series=np.array([0.0]), about=None)
    assert bw == pytest.approx(2 * (hi - lo) / (hi + lo), abs=1e-3)      # the grid's own 0.0005
    assert centre == pytest.approx(0.5 * (hi + lo), abs=5e-4)
    # about f0 the lower edge sets it, 2 (1 - lo), unless a small series element recentres
    # the band - which can win back at most the rest of its full width
    about = st.vswr_band(RF, z, x_nonres=0.0)[0]
    assert 2 * (1 - lo) - 1e-3 <= about <= 2 * (hi - lo) / (hi + lo) + 1e-3


def test_the_series_element_takes_out_a_reactance():
    """The same resonator behind a fixed 30-ohm inductance at f0: the tuned series
    capacitor cancels it and restores the band."""
    q = 10.0
    z = 50.0 + 1j * 50.0 * q * (RF - 1 / RF)
    bare = st.vswr_band(RF, z, x_nonres=0.0)
    loaded = st.vswr_band(RF, z, x_nonres=30.0)
    assert loaded[0] == pytest.approx(bare[0], rel=0.06)               # the 30 ohm steepens it a little
    assert loaded[2] == pytest.approx(30.0, abs=1.5)                   # a capacitor of 30 ohm at f0
    assert st.vswr_band(RF, z, x_nonres=30.0, series=np.array([0.0]))[0] < 0.5 * bare[0]


def test_a_band_must_hold_f0():
    """A resonator tuned to 1.1 f0 has its band there, not at f0."""
    z = 50.0 + 1j * 500.0 * (RF / 1.1 - 1.1 / RF)
    assert st.vswr_band(RF, z, x_nonres=0.0, series=np.array([0.0]))[0] == 0.0
    assert st.vswr_band(RF, z, x_nonres=0.0, series=np.array([0.0]), about=None)[1] == pytest.approx(1.1, abs=2e-3)


def test_nothing_matches_a_short():
    bw, centre, xc = st.vswr_band(RF, np.full(len(RF), 1.0 + 0j), x_nonres=0.0)
    assert bw == 0.0 and math.isnan(centre) and math.isnan(xc)


@pytest.mark.slow
def test_the_command_on_a_narrow_sweep(capsys):
    from otahub.cli.main import main
    rc = main(["stack", "--f0", "2.4GHz", "--h", "1.6mm", "--gap", "0.09", "--ratio", "1.1",
               "--lo", "0.975", "--hi", "1.025", "--json"])
    assert rc == 0
    out = json.loads(capsys.readouterr().out)
    assert out["gap_m"] == pytest.approx(0.09 * 2.99792458e8 / 2.4e9, rel=1e-9)
    assert out["best"]["bandwidth"][1] > 0.04                  # the window itself is 5% wide
    assert 7.0 < 10 * math.log10(out["directivity"]) < 11.0


SURVEY = json.loads((Path(__file__).parent / "data" / "stacked_patch_design.json").read_text())
LAM = 2.99792458e8 / 2.4e9


def _synth(registry, board, design, **extra):
    return registry["stacked_patch"].synthesize(f0=2.4e9, eps_r=board["eps_r"], h=board["h_lambda"] * LAM, eps_r2=1.0,
                                                h2_over_lambda=design["gap"], size_ratio=design["L2"] / design["L"], **extra)


def test_the_spec_directivity_is_the_surveys(registry):
    """The fit against every solved design, the held-out board's included."""
    for board in SURVEY["boards"]:
        for design in board["designs"]:
            d = _synth(registry, board, design)
            assert d.get("L") / LAM == pytest.approx(design["L"], rel=1e-6)
            assert d.metrics["directivity_dbi"] == pytest.approx(design["directivity_dbi"], abs=0.1), (board["name"], design)


def test_the_best_band_is_about_eight_times_the_driven_patchs(registry):
    ratios = []
    for board in SURVEY["boards"]:
        best = max(d["band"] for d in board["designs"])
        d = _synth(registry, board, board["designs"][0])
        ratios.append(best / d.metrics["single_patch_bandwidth_vswr2"])
        assert d.metrics["best_stack_bandwidth_vswr2"] == pytest.approx(best, rel=0.08), board["name"]
    assert 7.5 < min(ratios) and max(ratios) < 9.2
    fit = [r for r, b in zip(ratios, SURVEY["boards"]) if b["role"] == "fit"]
    held = next(r for r, b in zip(ratios, SURVEY["boards"]) if b["role"] == "held out")
    assert sum(fit) / len(fit) == pytest.approx(held, rel=0.08)          # the rule without it predicts it


def test_the_default_is_the_default_boards_best_and_the_old_one_was_not(registry):
    board = next(b for b in SURVEY["boards"] if b["name"] == "default")
    best = max(board["designs"], key=lambda d: d["band"])
    spec = registry["stacked_patch"].synthesize(f0=2.4e9, eps_r=2.2, h=0.0016)
    assert best["gap"] == pytest.approx(spec.get("h2") / LAM, rel=1e-6)
    assert best["L2"] / best["L"] == pytest.approx(spec.get("size_ratio"), rel=1e-6)
    old = next(d for d in board["designs"] if abs(d["gap"] - 0.03) < 1e-9)
    assert old["band"] < 1.1 * spec.metrics["single_patch_bandwidth_vswr2"]   # a stack in name only
    assert best["band"] > 7 * old["band"]


def test_the_stored_survey_reproduces():
    """One design's patch-current impedance at f0, solved again."""
    from otahub.num import patch_sdm as sdm
    from otahub.core.registry import default_registry
    spot = SURVEY["spot"]
    d = default_registry()["stacked_patch"].synthesize(f0=2.4e9, eps_r=2.2, h=0.0016, eps_r2=1.0,
                                                       h2_over_lambda=spot["gap"], size_ratio=spot["ratio"])
    p = sdm.StackedPatch(2.2, 0.0016, d.get("L"), d.get("W"), 1.0, d.get("h2"), d.get("L2"), d.get("W2"), **sdm.BASIS_FULL)
    xs = np.array(spot["xp_fraction"][-2:]) * d.get("L") / 2
    v = sdm.probe_vector(p, 2.4e9, xs, spot["probe_radius_lambda"] * LAM)
    z = -np.einsum("ij,ij->i", v, np.linalg.solve(p.Z(2.4e9), v.T).T)
    for got, want in zip(z, spot["z_at_f0"][-2:]):
        assert got == pytest.approx(complex(*want), rel=5e-3)
