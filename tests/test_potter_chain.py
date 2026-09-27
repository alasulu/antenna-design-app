"""The Potter horn as one chain: step, phasing guide and flare together.

The spec designs the step for a launched TM11 share and the phasing length from
the modes' local propagation constants. The whole chain (GSM cascade of mode
matching, `waveguide_step.potter_aperture`) shows the phasing guide is a TM11
resonator - TM11 is cut off in the input guide and partly reflected where the
flare begins - so the share and phase at the aperture belong to step and
phasing jointly. An FDTD of the same staircase (`bor_fdtd`) is the independent
check. These tests pin what that means for the spec's default horn.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from otahub.num import horn_pattern as hp
from otahub.num import waveguide_step as ws

DATA = json.loads((Path(__file__).parent / "data" / "potter_chain.json").read_text())
DM = math.sqrt(30.0)                                   # the 10-wavelength default's aperture


def _xpol(share, psi, fr=1.0):
    D, s = DM * fr, 0.375 * fr
    return hp.peak_cross(lambda t, f: hp.disc_pattern(hp.dual_mode(share, psi), D, s, t, f), D, math.pi / 4)


def test_the_phasing_guide_resonates():
    """Moving the phasing length through a beat period and more swings the
    aperture share from 0.04 to 0.27 about a launched 0.15."""
    s = [r["share"] for r in DATA["ell_scan"]]
    assert min(s) < 0.05 and max(s) > 0.25


def test_the_fdtd_sees_the_same_swing():
    rows = DATA["fdtd_ell"]
    for r in rows:
        assert abs((r["fdtd_psi"] - r["cascade_psi"] + 180) % 360 - 180) < 5
        assert r["fdtd_share"] == pytest.approx(r["cascade_share"], abs=0.045)
    fd = [r["fdtd_share"] for r in rows]
    assert max(fd) - min(fd) > 0.2                       # 0.13 to 0.34 at one step


@pytest.mark.parametrize("key", ["design_curve_10", "design_curve_5"])
def test_the_in_phase_share_is_not_the_launched_one(key):
    """For each step, the share at the first in-phase phasing length: it spans
    from nearly nothing to a quarter or more and does not rise steadily."""
    first = [r["in_phase"][0]["share"] for r in DATA[key] if r["in_phase"]]
    assert min(first) < 0.03 and max(first) > 0.23
    assert any(first[j] < max(first[:j]) - 0.02 for j in range(1, len(first)))   # it falls back somewhere


def test_the_joint_design_delivers_its_share_in_phase():
    d = DATA["design"]
    for k in ("stair_005", "stair_0025"):
        assert d[k]["share"] == pytest.approx(0.15, abs=5e-4) and abs(d[k]["psi"]) < 1e-3
    assert d["stair_0025"]["d_step"] == pytest.approx(d["stair_005"]["d_step"], rel=3e-3)
    assert d["stair_0025"]["ell"] == pytest.approx(d["stair_005"]["ell"], rel=0.015)
    conv = DATA["stair_convergence"]
    assert all(abs(r["psi"]) < 0.5 for r in conv) and conv[-1]["share"] - conv[0]["share"] < 0.01
    fo = DATA["first_order"]
    assert d["stair_0025"]["d_step"] > fo["d_step"] + 0.1 and d["stair_0025"]["ell"] > fo["ell"] + 0.25


def test_the_fdtd_follows_the_joint_design():
    for r in DATA["fdtd_design"]:
        assert abs(r["fdtd_psi"] - r["cascade_psi"]) < 1.0
        assert r["fdtd_share"] == pytest.approx(r["cascade_share"], abs=0.02)


def test_the_first_order_horn_misses_its_cross_polar_null():
    """The spec's own default, built first-order: -21 dB at f0 where the aperture
    model promises -32 for 0.15 in phase; its -30 dB window lies wholly above f0."""
    rows = {r["fr"]: r for r in DATA["band"]["first"]}
    assert _xpol(rows[1.0]["share"], rows[1.0]["psi"]) > -23
    assert _xpol(rows[1.02]["share"], rows[1.02]["psi"], 1.02) < -30
    assert _xpol(rows[0.99]["share"], rows[0.99]["psi"], 0.99) > -20


def test_the_joint_design_holds_its_null_across_a_band():
    rows = DATA["band"]["casc"]
    x = {r["fr"]: _xpol(r["share"], r["psi"], r["fr"]) for r in rows}
    assert x[1.0] < -32
    assert x[0.995] < -30 if 0.995 in x else x[0.99] < -29
    assert max(x[f] for f in (1.0, 1.01, 1.02)) < -30          # -30 dB from 0.992 to 1.028
    assert max(x[f] for f in (0.97, 0.98, 1.03, 1.04, 1.05)) < -25
    psi = {r["fr"]: r["psi"] for r in rows}
    assert -7 < (psi[1.01] - psi[0.99]) / 2 < -4                # deg per percent; the first-order horn's is -12.7


@pytest.mark.slow
def test_the_chain_reproduces_the_design():
    d = DATA["design"]["stair_005"]
    share, psi = ws.potter_aperture(1.1, d["d_step"], d["ell"], DM, 10.0)
    assert share == pytest.approx(d["share"], abs=1e-6) and abs(psi - d["psi"]) < 1e-4
