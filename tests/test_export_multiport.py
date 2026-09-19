"""Builders that export more than one port, and why each extra one is there.

An extra port is never decoration. In every case here the model answers a
different question without it - a loaded whip becomes an unloaded one, a
travelling-wave wire becomes a standing-wave one - and the resulting file runs
cleanly either way. These tests pin the ports and the notes that explain them.
"""
from __future__ import annotations

import pytest

from otahub.export import build


def test_loaded_monopole_exports_a_port_for_its_coil(registry):
    """Port 2 is the loading coil's gap. Left open, the whip is just a short
    whip and the spec's resonance never happens."""
    model = build(registry["inductively_loaded_monopole"].synthesize(
        f0=10e6, h_over_lambda=0.05))
    assert len(model.ports) == 2
    joined = " ".join(model.notes).lower()
    assert "loading coil" in joined
    assert "loading_inductance_h" in joined, (
        "the note must name the spec value to put into that port")


def test_loaded_monopole_coil_gap_splits_the_rod(registry):
    design = registry["inductively_loaded_monopole"].synthesize(
        f0=10e6, h_over_lambda=0.05)
    model = build(design)
    lower = next(s for s in model.solids if s.name == "rod_lower")
    upper = next(s for s in model.solids if s.name == "rod_upper")
    assert upper.span[0] > lower.span[1], "the coil gap must be a real break"
    coil_port = model.ports[1]
    assert coil_port.start[2] == pytest.approx(lower.span[1], rel=1e-9)
    assert coil_port.end[2] == pytest.approx(upper.span[0], rel=1e-9)


def test_top_loaded_monopole_hat_sits_on_top_of_the_rod(registry):
    design = registry["top_loaded_monopole"].synthesize(
        f0=10e6, h_over_lambda=0.05, beta_top=0.6)
    model = build(design)
    rod = next(s for s in model.solids if s.name == "rod")
    hat = next(s for s in model.solids if s.name == "top_hat")
    assert hat.span[0] == pytest.approx(rod.span[1], rel=1e-9)
    assert hat.radius == pytest.approx(design.get("a_hat"), rel=1e-9)
    assert hat.radius > rod.radius


def test_top_loaded_monopole_says_the_hat_and_beta_are_not_linked(registry):
    """beta_top and the hat radius are both INPUTS to the spec, not derived
    from each other. A model that implied otherwise would mislead."""
    model = build(registry["top_loaded_monopole"].synthesize(
        f0=10e6, h_over_lambda=0.05, beta_top=0.6))
    joined = " ".join(model.notes).lower()
    assert "beta_top" in joined and "input" in joined


def test_multiturn_loop_builds_every_turn(registry):
    from otahub.export.base import Torus

    design = registry["multiturn_small_loop"].synthesize(
        f0=10e6, C_over_lambda=0.1, Rin_target=50.0)
    model = build(design)
    turns = [s for s in model.solids if isinstance(s, Torus)]
    assert len(turns) == int(round(design.get("N")))
    zs = sorted(t.centre[2] for t in turns)
    steps = [b - a for a, b in zip(zs, zs[1:])]
    assert max(steps) == pytest.approx(min(steps), rel=1e-9), "even pitch"


def test_multiturn_loop_warns_the_turns_are_not_connected(registry):
    """The rings are separate solids. Solved as drawn, that is one driven ring
    and N-1 parasitic ones, and the spec's N-squared resistance never appears."""
    model = build(registry["multiturn_small_loop"].synthesize(
        f0=10e6, C_over_lambda=0.1, Rin_target=50.0))
    joined = " ".join(model.notes).lower()
    assert "connect the turns" in joined
    assert "parasitic" in joined


def test_arbitrary_length_dipole_reuses_the_dipole_builder(registry):
    """It is the same two arms and a gap; a separate copy would be a second
    thing to keep in step."""
    a = build(registry["dipole_arbitrary_length"].synthesize(
        f0=300e6, L_over_lambda=0.5, aw=1e-4))
    b = build(registry["half_wave_dipole"].synthesize(f0=300e6, aw=1e-4))
    assert [s.name for s in a.solids] == [s.name for s in b.solids]
    # and at 0.5 lambda the two must produce the same arms
    for x, y in zip(a.solids, b.solids):
        assert x.span == pytest.approx(y.span, rel=1e-9)
