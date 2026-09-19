"""Lens exporter tests.

Both lens builders make their geometry out of concentric shells, which is the
one place where a boolean that quietly does nothing leaves a solid disc or a
solid ball behind - and a solid ball is not a Luneburg lens, it is a paperweight.
"""
from __future__ import annotations

import pytest

from otahub.export import build
from otahub.export.base import Subtract


def test_zone_plate_radii_agree_with_the_spec(registry):
    """The builder recomputes the zone radii from its own copy of
    r(m) = sqrt(m*lambda*F + (m*lambda/2)^2). If that ever drifts from the
    spec's, the exported plate focuses somewhere the design sheet does not say.
    """
    design = registry["fresnel_zone_plate"].synthesize(f0=30e9, F=0.15, M=4)
    model = build(design)
    radii = sorted(s.radius for s in model.solids)
    assert radii[0] == pytest.approx(design.get("r_first"), rel=1e-9)
    assert radii[-1] == pytest.approx(design.get("r_outer"), rel=1e-9)


def test_zone_plate_rings_are_annuli_not_discs(registry):
    model = build(registry["fresnel_zone_plate"].synthesize(f0=30e9, F=0.15, M=4))
    subs = [op for op in model.operations if isinstance(op, Subtract)]
    assert len(subs) == 2, "two metal zones out of four"
    for sub in subs:
        outer = next(s for s in model.solids if s.name == sub.target)
        inner = next(s for s in model.solids if s.name == sub.tools[0])
        assert inner.radius < outer.radius


@pytest.mark.parametrize("zones,rings", [(2, 1), (4, 2), (6, 3), (8, 4)])
def test_zone_plate_ring_count_follows_the_zone_count(zones, rings, registry):
    model = build(registry["fresnel_zone_plate"].synthesize(
        f0=30e9, F=0.15, M=zones))
    assert len(model.operations) == rings


def test_luneburg_shells_span_the_right_permittivity_range(registry):
    """n(r) = sqrt(2 - (r/R)^2), so permittivity runs from 2 at the core to 1 at
    the rim. A stepped lens must bracket that without ever leaving it."""
    model = build(registry["luneburg_lens"].synthesize(f0=30e9, D=0.3))
    eps = sorted(float(s.material.split("=")[1]) for s in model.solids
                 if s.material.startswith("eps_r="))
    assert len(eps) == 8
    assert 1.0 < eps[0] < 1.2, f"outermost shell {eps[0]} should sit near 1"
    assert 1.8 < eps[-1] < 2.0, f"core shell {eps[-1]} should sit near 2"
    assert all(a < b for a, b in zip(eps, eps[1:])), "must increase inward"


def test_luneburg_shells_are_hollowed_except_the_core(registry):
    design = registry["luneburg_lens"].synthesize(f0=30e9, D=0.3)
    model = build(design)
    subs = [op for op in model.operations if isinstance(op, Subtract)]
    assert len(subs) == 7, "seven shells hollowed, the innermost left solid"
    outer = max(s.radius for s in model.solids
                if s.material.startswith("eps_r="))
    assert outer == pytest.approx(design.get("R"), rel=1e-9)


def test_luneburg_shell_boundaries_are_equally_spaced_in_r_squared(registry):
    """Permittivity is exactly linear in r^2, so equal steps in r^2 give every
    shell the same permittivity span. Equal steps in r would not."""
    model = build(registry["luneburg_lens"].synthesize(f0=30e9, D=0.3))
    radii = sorted(s.radius for s in model.solids
                   if s.material.startswith("eps_r="))
    squares = [r * r for r in radii]
    steps = [b - a for a, b in zip(squares, squares[1:])]
    assert max(steps) / min(steps) == pytest.approx(1.0, rel=1e-9)


def test_lenses_carry_no_port_and_say_why(registry):
    """Neither lens is driven - both are illuminated by a separate feed. A model
    with no port and no explanation reads as unfinished."""
    for key, given in (("fresnel_zone_plate", {"f0": 30e9, "F": 0.15, "M": 4}),
                       ("luneburg_lens", {"f0": 30e9, "D": 0.3})):
        model = build(registry[key].synthesize(**given))
        assert not model.ports
        assert any("feed" in n.lower() for n in model.notes), key
