"""Designing a Potter horn: step and phasing length solved together, and which solution.

The whole-chain cascade (`waveguide_step`) showed the phasing guide is a TM11
resonator, so the share that arrives in phase belongs to step and phasing
jointly. `potter_design` solves them jointly for any horn - scanning the step,
following every in-phase phasing length as a branch, refining each crossing of
the target share - and `otahub potter` puts it in front of a designer.

What it found is that "the" joint design is several: at a share of 0.13 the
default horn has three, whose aperture phase moves -5.0, -8.1 and -11.0 degrees
per percent of frequency, so the cross-polar band of one is twice another's.
A guide one beat longer is in phase again and delivers another share - more
designs, all steeper. An FDTD of the flattest and steepest staircases, at three
frequencies, confirms the slopes to within 0.6 degree per percent.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from otahub.num import waveguide_step as ws

HERE = Path(__file__).parent / "data"
DATA = json.loads((HERE / "potter_design.json").read_text())
CHAIN = json.loads((HERE / "potter_chain.json").read_text())
DM = math.sqrt(30.0)


def test_a_free_phasing_length_is_the_whole_chain():
    ch = ws.PotterChain(1.1, 1.43, DM, 10.0)
    for ell in (0.2, 0.9, 2.4):
        share, psi = ch.at(ell)
        want = ws.potter_aperture(1.1, 1.43, ell, DM, 10.0)
        assert share == pytest.approx(want[0], abs=1e-10) and psi == pytest.approx(want[1], abs=1e-8)


def test_the_faster_cascade_is_the_same_cascade():
    """Mode profiles by the J1' recurrence, all modes at once, and a matrix product for
    the coupling integrals: the recorded joint design's aperture, unchanged."""
    d = CHAIN["design"]["stair_005"]
    share, psi = ws.potter_aperture(1.1, d["d_step"], d["ell"], DM, 10.0)
    assert share == pytest.approx(d["share"], abs=1e-9) and abs(psi - d["psi"]) < 1e-7


def test_the_search_refines_to_the_recorded_joint_design():
    got = ws.potter_design(1.1, DM, 10.0, 0.15, steps=np.arange(1.46, 1.505, 0.01))
    d = CHAIN["design"]["stair_005"]
    assert got[0]["d_step"] == pytest.approx(d["d_step"], abs=2e-4)
    assert got[0]["ell"] == pytest.approx(d["ell"], abs=1e-3)
    assert -6.0 < got[0]["phase_slope_deg_per_percent"] < -5.0


@pytest.mark.parametrize("i", [0, 1, 2])
def test_three_designs_deliver_the_same_share_in_phase(i):
    """Checked by the plain cascade, which knows nothing of the search."""
    r = DATA["share_013"][i]
    share, psi = ws.potter_aperture(1.1, r["d_step"], r["ell"], DM, 10.0)
    assert share == pytest.approx(0.13, abs=5e-4) and abs(psi) < 0.05


def test_the_designs_differ_twofold_in_bandwidth():
    slopes = [r["phase_slope_deg_per_percent"] for r in DATA["share_013"]]
    assert len(slopes) == 3 and all(s < 0 for s in slopes)
    assert slopes[-1] / slopes[0] > 2.0                       # listed flattest first
    steps = [r["d_step"] for r in DATA["share_013"]]
    assert max(steps) - min(steps) > 0.09                     # not one design found three times


def test_the_fdtd_confirms_the_slopes():
    runs = {r["tag"]: r for r in DATA["fdtd_slope"]}
    for r in runs.values():
        assert abs(r["fdtd_slope"] - r["cascade_slope"]) < 0.8
        for row in r["rows"]:
            assert abs(row["fdtd_psi"] - row["cascade_psi"]) < 2.5
            assert row["fdtd_share"] == pytest.approx(row["cascade_share"], abs=0.02)
    assert runs["steep"]["fdtd_slope"] / runs["flat"]["fdtd_slope"] > 1.9


def test_a_longer_guide_is_another_design():
    """A beat longer is in phase again, but the TM11 rattling in it has gone round a
    different phase: other steps, and every one of them steeper."""
    sols = DATA["share_015_two_periods"]
    beat = ws.PotterChain(1.1, 1.48, DM, 10.0).beat
    short = [r for r in sols if r["ell"] < beat]
    long_ = [r for r in sols if r["ell"] >= beat]
    assert len(short) == 1 and len(long_) >= 2
    assert all(abs(r["phase_slope_deg_per_percent"]) > 2 * abs(short[0]["phase_slope_deg_per_percent"]) for r in long_)


@pytest.mark.slow
def test_the_command_designs_the_spec_default(capsys):
    from otahub.cli.main import main
    assert main(["potter", "--f0", "10GHz", "--L", "0.3", "--json"]) == 0
    out = json.loads(capsys.readouterr().out)
    lam = out["wavelength_m"]
    fo = out["first_order"]
    assert fo["phase_deg"] == pytest.approx(26.5, abs=0.5) and fo["cross_pol_db"] > -23
    best = out["solutions"][0]
    assert best["d_step_m"] / lam == pytest.approx(1.4807, abs=2e-3)
    assert best["share"] == pytest.approx(0.15, abs=5e-4)
    assert out["target_cross_pol_db"] < -32
