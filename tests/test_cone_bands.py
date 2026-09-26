"""Biconical and conical monopole bands, solved as bodies of revolution.

Both specs sized the cone with a quarter-wave slant and asserted an
'indicative' decade. `otahub.num.bor` solves the solid cones exactly. Matched to
its own characteristic impedance a bicone holds a decade from 10 to 65 degrees,
from a slant of 0.22-0.29 wavelengths; a conical monopole in 50 ohm holds it
only from 30 to 55 degrees, and below about 34 degrees its cutoff climbs
steeply because the cone's impedance leaves 50 ohm behind. Every check reads
the recorded solutions in tests/data or solves the spec's own dimensions live.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from otahub.num import bor

C0 = 2.99792458e8
DATA = json.loads((Path(__file__).parent / "data" / "cones_bor.json").read_text())


def _bic(registry, th, f0=1e9):
    return registry["biconical"].synthesize(f0=f0, theta_h=math.radians(th))


def _mono(registry, th, f_low=1e9):
    return registry["conical_monopole"].synthesize(f_low=f_low, cone_half_angle_deg=th)


def _vswr(z, z0):
    g = abs((z - z0) / (z + z0))
    return (1 + g) / (1 - g)


def test_bicone_tables_are_the_recorded_solutions(registry):
    for r in DATA["bic"]:
        d = _bic(registry, r["th"])
        assert d.get("k_slant") == pytest.approx(r["k"], rel=1e-4)
        assert d.metrics["bandwidth_ratio"] == pytest.approx(r["band"], rel=1e-3)
        assert d.metrics["worst_vswr_2_to_10_f0"] == pytest.approx(r["worst"], rel=1e-3)
        assert d.metrics["directivity_linear"] == pytest.approx(r["D"], rel=1e-4)


def test_monopole_tables_are_the_recorded_solutions(registry):
    for r in DATA["mono"]:
        d = _mono(registry, r["th"])
        assert d.get("k_slant") == pytest.approx(r["k"], rel=1e-4)
        assert d.metrics["bandwidth_ratio"] == pytest.approx(r["band"], rel=1e-3)
        assert d.metrics["directivity_linear"] == pytest.approx(r["D"], rel=1e-4)


def test_the_decade_is_the_bicones_and_only_sometimes_the_monopoles():
    bic = {r["th"]: r for r in DATA["bic"]}
    mono = {r["th"]: r for r in DATA["mono"]}
    assert all(bic[t]["band"] >= 9.99 for t in bic if t >= 10)
    assert all(mono[t]["band"] >= 9.99 for t in mono if 30 <= t <= 55)
    assert mono[20.0]["band"] < 2 and mono[65.0]["band"] < 2


def test_the_monopoles_cutoff_climbs_steeply_below_34_degrees():
    k = {r["th"]: r["k"] for r in DATA["mono"]}
    assert k[35.0] < 0.25 and k[34.0] < 0.25
    assert k[33.0] > 0.29 and k[30.0] > 0.39
    assert [k[t] for t in (30.0, 32.0, 32.5, 33.0, 34.0)] == sorted([k[t] for t in (30.0, 32.0, 32.5, 33.0, 34.0)], reverse=True)


def _bicone_profile(slant_lam, th, top_d_lam, gap_lam):
    t = math.radians(th)
    a = top_d_lam / 2
    g = gap_lam / 2
    zb = g + (slant_lam - a / math.sin(t)) * math.cos(t)
    rb = slant_lam * math.sin(t)
    return bor.profile([(rb, -zb), (a, -g), (a, g), (rb, zb)], slant_lam / 24, gap=(1, 2), gap_seg=4)


@pytest.mark.slow
def test_the_bicones_own_dimensions_reach_vswr_2_at_f0(registry):
    th, f0 = 30.0, 1e9
    d = _bic(registry, th, f0)
    lam = C0 / f0
    sol = bor.solve(_bicone_profile(d.get("Lc") / lam, th, d.get("cone_top_diameter") / lam, d.get("feed_gap") / lam))
    assert _vswr(sol.input_impedance, d.metrics["characteristic_impedance_ohm"]) == pytest.approx(2.0, abs=0.05)


@pytest.mark.slow
def test_the_monopoles_own_dimensions_reach_vswr_2_in_50_ohm_at_f_low(registry):
    """Solved as the cone plus its image: a bicone of twice the gap, half the impedance."""
    th, f_low = 40.0, 1e9
    d = _mono(registry, th, f_low)
    lam = C0 / f_low
    sol = bor.solve(_bicone_profile(d.get("slant") / lam, th, d.get("cone_top_diameter") / lam, 2 * d.get("feed_gap") / lam))
    assert _vswr(sol.input_impedance / 2, 50.0) == pytest.approx(2.0, abs=0.05)
    assert 2 * bor.directivity(sol)[0] == pytest.approx(d.metrics["directivity_linear"], rel=0.01)
