"""The discone, solved as a body of revolution rather than sized by rule of thumb.

It carried three conventions: a quarter-wave slant at f_low, a decade of VSWR
under 2, and 2.2 directivity. `otahub.num.bor` solves the solid cone and disc
exactly (axisymmetric surface currents, finite feed gap). The quarter wave is
short everywhere - 12% at the default, 58% on a 50 degree cone; the decade is
real for cones of 30 degrees or more with a coax-sized feed and not for narrower
ones; the directivity at f_low is 1.24-1.54, not 2.2. Every check below reads the
recorded survey in tests/data or rebuilds the antenna from the spec's own
dimensions and solves it live.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from otahub.num import bor

C0 = 2.99792458e8
DATA = json.loads((Path(__file__).parent / "data" / "discone_bor.json").read_text())


def _design(registry, th, ratio, f_low=100e6):
    return registry["discone"].synthesize(f_low=f_low, cone_half_angle_deg=th, disc_ratio=ratio)


def _profile_from_spec(d, freq):
    """The spec's own dimensions, in wavelengths at `freq`: cone truncated at
    cone_top_diameter with its top face at z = 0, disc feed_gap above, a feed
    tube of the top radius across the gap."""
    lam = C0 / freq
    rb = d.get("cone_base_diameter") / 2 / lam
    rt = d.get("cone_top_diameter") / 2 / lam
    h = d.get("cone_height") / lam
    s = d.get("feed_gap") / lam
    rd = d.get("disc_diameter") / 2 / lam
    slant = d.get("slant") / lam
    return bor.profile([(rb, -h), (rt, 0.0), (rt, s), (rd, s)], slant / 24, gap=(1, 2), gap_seg=3)


def _vswr(z, z0=50.0):
    g = abs((z - z0) / (z + z0))
    return (1 + g) / (1 - g)


@pytest.mark.parametrize("group", ["grid", "held"])
def test_low_cutoff_is_the_solved_one(registry, group):
    worst = 0.0
    for r in DATA[group]:
        d = _design(registry, r["th"], r["ratio"])
        worst = max(worst, abs(d.get("k_slant") / r["k"] - 1))
    assert worst < DATA["tolerance"]["k_" + group], worst


def test_tables_reproduce_the_survey_at_its_grid_points(registry):
    for r in DATA["grid"]:
        d = _design(registry, r["th"], r["ratio"])
        assert d.metrics["bandwidth_ratio"] == pytest.approx(r["band"], rel=2e-3)
        assert d.metrics["worst_vswr_2_to_10_f_low"] == pytest.approx(r["worst"], rel=2e-3)
        assert d.metrics["directivity_dbi"] == pytest.approx(10 * math.log10(r["D0"]), abs=0.01)


def test_the_decade_belongs_to_cones_of_30_degrees_and_more():
    """The claim this spec used to make as 'indicative', tested against every
    design surveyed: true from 30 degrees up at every disc ratio, false below."""
    for r in DATA["grid"]:
        if r["th"] >= 30:
            assert r["band"] >= 9.99 and r["worst"] <= 2.0 + 1e-3, r
        elif r["th"] <= 20:
            assert r["band"] < 5.0 and r["worst"] > 2.0, r


def test_a_quarter_wave_slant_is_always_short():
    assert min(r["k"] for r in DATA["grid"]) > 0.255
    assert max(r["k"] for r in DATA["grid"] if r["th"] >= 45) > 0.39


@pytest.mark.slow
@pytest.mark.parametrize("th,ratio", [(30.0, 0.7), (37.5, 0.65)])
def test_the_spec_dimensions_cross_vswr_2_at_f_low(registry, th, ratio):
    """Rebuilt from the spec's own outputs - slant, top diameter, spacing, disc -
    and solved at f_low: the VSWR must be 2 there, which ties the tables to the
    geometry the spec hands out."""
    f_low = 100e6
    d = _design(registry, th, ratio, f_low)
    sol = bor.solve(_profile_from_spec(d, f_low))
    assert _vswr(sol.input_impedance) == pytest.approx(2.0, abs=0.06)
    D, _ = bor.directivity(sol)
    assert 10 * math.log10(D) == pytest.approx(d.metrics["directivity_dbi"], abs=0.1)


def test_grading_the_free_edges_leaves_the_cutoffs_where_they_were():
    """The survey never graded the cone's base rim or the disc's rim. Rebuilt with
    both graded, the extrapolated low cutoff lands within 0.05% of the table at six
    designs; the ungraded rebuild reproduces the table's own extrapolation."""
    for r in DATA["graded_check"]["cutoff"]:
        assert r["ungraded_extrap"] == pytest.approx(r["k_table"], abs=3e-5)
        assert r["graded_extrap"] == pytest.approx(r["k_table"], rel=6e-4)
        assert abs(r["n48_lev8"] / r["k_table"] - 1) < abs(r["n48_lev0"] / r["k_table"] - 1)


def test_the_decade_is_marginal_where_the_ripple_touches_two():
    """Grading moves the worst in-band VSWR by under 0.01; at 30 degrees with a 0.9
    disc (and 25 with 0.6) it sits at 2.00 within the mesh's resolution."""
    w = {(r["th"], r["ratio"]): r for r in DATA["graded_check"]["worst"]}
    for r in w.values():
        assert abs(r["graded"] - r["ungraded"]) < 0.01
    for key in ((30.0, 0.9), (25.0, 0.6)):
        assert abs(w[key]["graded"] - 2.0) < 0.006
    assert w[(30.0, 0.8)]["graded"] < 1.98 and w[(25.0, 0.9)]["graded"] > 2.07
