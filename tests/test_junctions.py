"""Wire junctions: where the ends of several wires meet at one node.

Before this, wires that merely touched were not connected at all, which put a
whole class of real antennas out of reach - a top hat of radial wires, a
ground plane of radials, a discone built as a wire cage. Each basis function now
stores the orientation of its two halves, so a wire can END or START at a node
and still carry current through it, and K coincident wire ends get K - 1
junction functions carrying current from the first wire into each of the rest.
Kirchhoff's current law then holds by construction.

The checks, in increasing order of what they cover:

- refactoring the layout changed nothing for any existing geometry, bit for bit;
- a dipole split at its centre is the same dipole, whichever way each arm runs;
- a bent wire in three pieces, one reversed, is the same bent wire;
- a symmetric T splits its current equally between its arms, balances power,
  and keeps the impedance matrix exactly symmetric;
- a ground-plane vertical on four radials reproduces the two familiar rules of
  thumb once it is actually resonated - about 22 ohm with the radials flat,
  about 50 with them drooped 45 degrees.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from otahub.num import mom

L, A, N = 0.47, 1e-3, 40


def _split_dipole(reverse_lower=False, reverse_upper=False):
    z = np.linspace(-L / 2, L / 2, N + 1)
    pts = np.stack([np.zeros_like(z), np.zeros_like(z), z], axis=1)
    lower, upper = pts[: N // 2 + 1], pts[N // 2:]
    if reverse_lower:
        lower = lower[::-1].copy()
    if reverse_upper:
        upper = upper[::-1].copy()
    return mom.WireModel([mom.Wire(lower, A), mom.Wire(upper, A)])


@pytest.mark.parametrize("rl,ru", [(False, False), (False, True), (True, False), (True, True)])
def test_a_split_dipole_is_the_same_dipole_whichever_way_its_arms_run(rl, ru):
    """The orientation logic, isolated: a junction of two wires is just a bend."""
    want = mom.input_impedance(mom.dipole(L, A, N))
    m = _split_dipole(rl, ru)
    got = mom.solve(m, m.junction_at([0, 0, 0])).input_impedance
    assert abs(got - want) / abs(want) < 1e-12


def test_a_bent_wire_in_three_pieces_is_the_same_bent_wire():
    arc = mom.arc(0.6, 3.0, A, 36).wires[0].nodes
    one = mom.WireModel([mom.Wire(arc, A)])
    three = mom.WireModel([mom.Wire(arc[:13], A),
                           mom.Wire(arc[12:25][::-1].copy(), A),
                           mom.Wire(arc[24:], A)])
    feed_pt = one.node_of(one.n_basis // 2)
    k = int(np.argmin(np.linalg.norm(
        np.array([three.node_of(i) for i in range(three.n_basis)]) - feed_pt, axis=1)))
    z1 = mom.solve(one, one.n_basis // 2).input_impedance
    z3 = mom.solve(three, k).input_impedance
    assert abs(z3 - z1) / abs(z1) < 1e-12


def _tee():
    vert = np.stack([np.zeros(21), np.zeros(21), np.linspace(-0.3, 0.0, 21)], axis=1)
    right = np.stack([np.linspace(0, 0.2, 11), np.zeros(11), np.zeros(11)], axis=1)
    left = np.stack([np.linspace(0, -0.2, 11), np.zeros(11), np.zeros(11)], axis=1)
    return mom.WireModel([mom.Wire(vert, A), mom.Wire(right, A), mom.Wire(left, A)])


def test_a_symmetric_tee_splits_its_current_equally():
    m = _tee()
    j = m.junction_at([0, 0, 0])
    assert len(j) == 2, "three wires meet, so two junction functions"
    sol = mom.solve(m, 9)
    assert abs(sol.currents[j[0]] - sol.currents[j[1]]) / abs(sol.currents[j[0]]) < 1e-12


def test_a_tee_balances_power_and_stays_reciprocal():
    m = _tee()
    sol = mom.solve(m, 9)
    assert mom.radiated_power(sol, 50, 50) == pytest.approx(sol.circuit_power, rel=1e-4)
    Z = mom.impedance_matrix(m)
    assert np.abs(Z - Z.T).max() == 0.0


def test_wires_that_do_not_touch_are_not_joined():
    """Junctions come only from coincident END points."""
    d1 = mom.dipole(0.47, A, 20).wires[0]
    d2 = mom.Wire(d1.nodes + np.array([0.1, 0.0, 0.0]), A)
    assert mom.WireModel([d1, d2]).junctions == []
    with pytest.raises(ValueError):
        mom.WireModel([d1, d2]).junction_at([0, 0, 0])


def _ground_plane(droop_deg, scale):
    q = 0.25 * scale
    n = 13
    wires = [mom.Wire(np.stack([np.zeros(n), np.zeros(n), np.linspace(0, q, n)], axis=1), A)]
    dr = math.radians(droop_deg)
    for k in range(4):
        phi = k * math.pi / 2
        r = np.linspace(0, q, n)
        wires.append(mom.Wire(np.stack([r * math.cos(dr) * math.cos(phi),
                                        r * math.cos(dr) * math.sin(phi),
                                        -r * math.sin(dr)], axis=1), A))
    m = mom.WireModel(wires)
    return mom.solve(m, m.junction_at([0, 0, 0]))


@pytest.mark.parametrize("droop,want", [(0.0, 22.0), (45.0, 50.0)])
def test_the_ground_plane_vertical_reproduces_its_rules_of_thumb(droop, want):
    """A five-wire junction, fed there. At RESONANCE - the unresonated numbers
    are not comparable, carrying tens of ohms of reactance - the familiar
    figures come out: about 22 ohm flat, about 50 drooped 45 degrees."""
    lo, hi = 0.85, 1.15
    for _ in range(24):
        mid = 0.5 * (lo + hi)
        if _ground_plane(droop, mid).input_impedance.imag < 0:
            lo = mid
        else:
            hi = mid
    sol = _ground_plane(droop, lo)
    assert sol.input_impedance.real == pytest.approx(want, rel=0.08)
    assert mom.radiated_power(sol, 50, 50) == pytest.approx(sol.circuit_power, rel=1e-4)
