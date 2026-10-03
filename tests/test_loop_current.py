"""The small loops' efficiency and Q rest on the loop's own current, not its terminals.

A resistive wire loses r |I|^2 wherever the current is, so a loop's efficiency is
R_loop / (R_loop + R_loss) with R_loop the radiation resistance referred to the
mean-square current - and the loaded Q is the lossless radiation Q times that
efficiency. small_circular_loop and small_square_loop built both on the TERMINAL
resistance with the loss added at the terminals; once the current is not uniform
(the far side of a 0.3-wavelength loop carries 1.8 times the feed's current) that
counted the loss at the wrong current, and read efficiency and Q 8-12% high and
sized the wire 30% too thin. R_loop and Q_rad also do not depend on the feed,
where the terminal impedance of a thick wire does - so they hold on the thick
wire the terminal fits are NaN on. tests/data/small_loop_current.json holds the
surveys; this file refits them and checks the spec against every run.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from otahub.num import loop_modal

C0 = 2.99792458e8
DATA = json.loads((Path(__file__).parent / "data" / "small_loop_current.json").read_text())
CIRCLE, SQUARE = DATA["circle"]["runs"], DATA["square"]["runs"]


def _held(runs, key):
    sizes = sorted({r[key] for r in runs})
    return set(sizes[1::3])


def _circle(registry, C, ba):
    f0 = 3e8
    lam = C0 / f0
    return registry["small_circular_loop"].synthesize(f0=f0, C=C * lam, b=ba * C * lam / (2 * math.pi), N=1)


def _square(registry, P, bs):
    f0 = 3e8
    s = P * C0 / f0 / 4
    return registry["small_square_loop"].synthesize(f0=f0, s=s, b=bs * s, N=1)


def test_the_circle_survey_does_not_depend_on_the_feed():
    for r in CIRCLE:
        assert r["R_loop_60"] == pytest.approx(r["R_loop_40"], rel=2e-4)
        assert r["Q_rad_60"] == pytest.approx(r["Q_rad_40"], rel=2e-4)
    moved = max(abs(r["Z_60"][0] / r["Z_40"][0] - 1) for r in CIRCLE)
    assert moved > 0.3                                         # the terminal resistance: 35%


@pytest.mark.parametrize("which", ["fit", "held"])
def test_the_circle_spec_meets_every_run(registry, which):
    held = _held(CIRCLE, "C")
    runs = [r for r in CIRCLE if (r["C"] in held) == (which == "held")]
    for r in runs:
        d = _circle(registry, r["C"], r["ba"])
        assert d.metrics["radiation_resistance_loop_current_ohm"] == pytest.approx(r["R_loop_60"], rel=2.5e-3)
        assert d.metrics["radiation_quality_factor"] == pytest.approx(r["Q_rad_60"], rel=1.1e-3)


@pytest.mark.parametrize("which", ["fit", "held"])
def test_the_square_spec_meets_every_run(registry, which):
    held = _held(SQUARE, "P")
    runs = [r for r in SQUARE if (r["P"] in held) == (which == "held")]
    for r in runs:
        d = _square(registry, r["P"], r["bs"])
        assert d.metrics["radiation_resistance_loop_current_ohm"] == pytest.approx(r["R_loop"], rel=1.2e-3)
        assert d.metrics["radiation_quality_factor"] == pytest.approx(r["Q_rad"], rel=1.8e-3)


def test_the_circle_fit_is_the_surveys_own():
    """Refit ln(R_loop / 20 pi^2 C^4) in C, C^2, C^3 from the stored runs (every third
    size held out) and evaluate it against the spec's coefficients at every size."""
    held = _held(CIRCLE, "C")
    fit = [r for r in CIRCLE if r["C"] not in held]
    A = np.array([[r["C"] ** k for k in (1, 2, 3)] for r in fit])
    y = np.array([math.log(r["R_loop_60"] / (20 * math.pi ** 2 * r["C"] ** 4)) for r in fit])
    c = np.linalg.lstsq(A, y, rcond=None)[0]
    assert c == pytest.approx([-0.0344071148, 4.22693769, -1.76692031], rel=1e-6)


@pytest.mark.parametrize("C, ba", [(0.3, 0.05), (0.17, 0.005)])
def test_a_stored_circle_run_reproduces_live(C, ba):
    r = next(x for x in CIRCLE if abs(x["C"] - C) < 1e-9 and x["ba"] == ba)
    a = ba * C / (2 * math.pi)
    assert loop_modal.loop_current_resistance(C, a, 60) == pytest.approx(r["R_loop_60"], rel=1e-9)
    assert loop_modal.tuned_loop(C, a, 0.0, 60)["q"] == pytest.approx(r["Q_rad_60"], rel=1e-9)


@pytest.mark.parametrize("C, ba, f0", [(0.25, 0.01, 100e6), (0.3, 0.05, 100e6), (0.1, 0.02, 30e6)])
def test_efficiency_and_q_are_the_lossy_loops(registry, C, ba, f0):
    """The identities the spec uses, against the loop solved with the wire's
    resistance in every mode: efficiency P_rad/P_in and the loaded Q_Z."""
    lam = C0 / f0
    a = ba * C / (2 * math.pi)
    rs = math.sqrt(math.pi * f0 * 4e-7 * math.pi / 5.8e7)
    solved = loop_modal.tuned_loop(C, a, rs, 60)
    d = registry["small_circular_loop"].synthesize(f0=f0, C=C * lam, b=a * lam, N=1)
    assert d.metrics["radiation_efficiency"] == pytest.approx(solved["efficiency"], rel=3e-3)
    assert d.metrics["quality_factor"] == pytest.approx(solved["q"], rel=4e-3)


@pytest.mark.parametrize("shape, key, size", [("circle", "small_circular_loop", "C_over_lambda"),
                                              ("square", "small_square_loop", "P_over_lambda")])
def test_the_directivity_is_the_solved_currents(registry, shape, key, size):
    """Gain and aperture assumed the infinitesimal loop's 1.5; at 0.3 wavelengths the
    solved current gives 1.436 (0.19 dB less) and fills the axis to 0.47."""
    for r in DATA["directivity"][shape]:
        d = registry[key].synthesize(f0=3e8, **{size: r["x"]})
        assert 10 ** (d.metrics["directivity_dbi"] / 10) == pytest.approx(r["D"], rel=2e-4)
    if shape == "circle":
        assert DATA["directivity"]["circle"][-1]["D_axial"] > 0.45


def test_a_thick_loop_has_no_terminal_impedance_and_no_stand_in(registry):
    """Where the terminal impedance is NaN (the feed decides it), the textbook radiation
    resistance used to stand in for it - in terminal_impedance() and the GUI's tile."""
    from otahub.gui.models import key_figures
    lam = C0 / 1e8
    d = registry["small_circular_loop"].synthesize(f0=1e8, C=0.3 * lam, b=0.05 * 0.3 * lam / (2 * math.pi), N=1)
    assert d.terminal_impedance() is None
    assert not [f for f in key_figures(d.metrics, d.units) if f[0] == "Impedance"]
    assert d.metrics["quality_factor"] > 0                    # what does not depend on the feed is there


@pytest.mark.parametrize("given, expect_nan", [({"eta_target": 0.0, "N": 1}, True), ({"eta_target": 1.0, "N": 1}, True),
                                               ({"eta_target": 0.5, "N": 0}, True), ({"eta_target": 0.5, "N": 1}, False)])
def test_the_wire_sizing_refuses_what_it_cannot_mean(registry, given, expect_nan):
    d = registry["small_circular_loop"].synthesize(f0=1e8, C_over_lambda=0.2, **given)
    ratio = d.get("wire_ratio_eta")
    assert (ratio is None or math.isnan(ratio)) == expect_nan


def test_a_target_a_thinner_wire_beats_takes_the_thinnest_fitted_one(registry):
    """At 1 MHz a 0.3-wavelength loop meets 0.5 with b/a = 1.2e-4, under the Q fit's
    0.001: it takes 0.001, and its Q and bandwidth stay defined."""
    d = registry["small_circular_loop"].synthesize(f0=1e6, C_over_lambda=0.3, eta_target=0.5, N=1)
    assert d.get("wire_ratio_eta") == pytest.approx(0.001)
    assert d.metrics["radiation_efficiency"] > 0.5 and math.isfinite(d.metrics["quality_factor"])
