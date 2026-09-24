"""The rhombic, solved with its termination in circuit.

The spec's directivity came from an idealised four-leg model: an unattenuated
travelling wave on every leg and a perfect termination. A real rhombic's
current decays as it radiates, and its termination is only as good as the
match to a characteristic impedance set by the wire. Solved here as it is built
- one closed wire, the feed at one acute vertex and the resistor at the other -
with the power that disappears into the resistor accounted for exactly.

The wire radius turned out to be the parameter the spec had left out. It sets
the termination, the share of the power that termination eats, and therefore
the gain, while barely touching the directivity.
"""
from __future__ import annotations

import functools
import math

import numpy as np
import pytest

from otahub.num import mom

pytestmark = pytest.mark.slow   # full-wave solves; skip with -m 'not slow'

MESH = 16          # directivity reads ~0.15 dB low here; fine for the physics


@functools.lru_cache(maxsize=32)
def _solve(leg, rt, radius=1e-4, half=None, seg=MESH):
    m, feed, term = mom.rhombic_model(leg, half, radius, seg)
    return m, term, mom.solve(m, feed, loads={term: rt})


def _front_back(sol):
    def u(ph):
        a, b = mom.far_field(sol, np.array([math.pi / 2]), np.array([ph]))
        return abs(a[0]) ** 2 + abs(b[0]) ** 2
    return 10 * math.log10(u(0.0) / u(math.pi))


def test_power_into_the_resistor_is_accounted_for():
    """Input power = radiated + dissipated in the termination, exactly. This
    is what makes the efficiency figure below trustworthy."""
    m, term, sol = _solve(4.0, 600.0)
    p_term = 0.5 * abs(sol.currents[term]) ** 2 * 600.0
    assert mom.radiated_power(sol, 50, 100) + p_term == pytest.approx(
        sol.circuit_power, rel=1e-4)


def test_it_fires_towards_the_terminated_apex():
    _, _, sol = _solve(4.0, 600.0)
    assert _front_back(sol) > 15.0


def test_without_a_termination_the_front_to_back_collapses():
    """An open rhombic reflects its wave and radiates both ways."""
    _, _, sol = _solve(4.0, 1e9)
    assert _front_back(sol) < 8.0


@pytest.mark.parametrize("radius", [3e-5, 1e-3])
def test_the_best_termination_depends_on_the_wire(radius):
    """Thinner wire, higher characteristic impedance, higher optimum: the
    front-to-back peaks at a higher resistance on the thinner wire."""
    fbs = {rt: _front_back(_solve(4.0, rt, radius)[2]) for rt in (300.0, 450.0, 700.0, 1000.0)}
    best = max(fbs, key=fbs.get)
    if radius < 1e-4:
        assert best >= 700.0, fbs
    else:
        assert best <= 450.0, fbs


def test_a_thinner_wire_sends_more_power_to_the_termination():
    """So 'about half the power is dissipated' is a statement about thin wire."""
    effs = []
    for radius, rt in ((3e-5, 800.0), (1e-3, 400.0)):
        m, term, sol = _solve(4.0, rt, radius)
        effs.append(mom.radiated_power(sol, 50, 100) / sol.circuit_power)
    assert effs[0] < effs[1]
    assert 0.45 < effs[0] < 0.6


def test_directivity_sits_below_the_idealised_model():
    """The idealised model has no attenuation along the legs; the real current
    does, so the real directivity is lower - but not by much at this length."""
    _, _, sol = _solve(4.0, 650.0, seg=28)
    d = 10 * math.log10(4 * math.pi * _u_axis(sol) / mom.pattern_power(sol, 50, 100)[1])
    ideal = 10 * math.log10(9.2038 * 4.0 ** 1.0404)
    assert 0.0 < ideal - d < 0.6


def _u_axis(sol):
    a, b = mom.far_field(sol, np.array([math.pi / 2]), np.array([0.0]))
    return abs(a[0]) ** 2 + abs(b[0]) ** 2


# ------------------------------------------------------------ the spec itself

LAM_10MHZ = 2.99792458e8 / 1e7


def _spec(registry, L=4.0, a_l=1e-4):
    return registry["rhombic"].synthesize(f0=1e7, L_over_lambda=L, aw=a_l * LAM_10MHZ)


def test_spec_termination_follows_the_wire(registry):
    got = [_spec(registry, a_l=a).metrics["termination_resistance_ohm"]
           for a in (1e-5, 1e-4, 1e-3, 3e-3)]
    assert got[0] > got[1] > got[2] > got[3], got
    assert 600 * 0.9 < got[1] < 600 * 1.15, "600 ohm is right for typical HF wire"
    assert got[2] < 450, "and well off for fat wire"


def test_spec_efficiency_and_gain(registry):
    thin, fat = _spec(registry, a_l=1e-5).metrics, _spec(registry, a_l=3e-3).metrics
    assert thin["radiation_efficiency"] < 0.5 < fat["radiation_efficiency"]
    for m in (thin, fat):
        assert m["gain_dbi"] == pytest.approx(
            m["directivity_dbi"] + 10 * math.log10(m["radiation_efficiency"]), abs=1e-9)


def test_spec_directivity_sits_a_few_tenths_below_the_idealised_model(registry):
    for L in (2.0, 4.0, 8.0, 12.0):
        m = _spec(registry, L=L).metrics
        gap = 10 * math.log10(m["directivity_travelling_wave_model"] / m["directivity_linear"])
        assert 0.15 < gap < 0.6, (L, gap)


@pytest.mark.parametrize("leg", [4.0, 8.0])
def test_spec_directivity_matches_a_fresh_solve(leg, registry):
    """At the fine mesh the fit was made on, with the spec's own termination."""
    m = _spec(registry, L=leg).metrics
    model, feed, term = mom.rhombic_model(leg, radius=1e-4, seg_per_lambda=28)
    sol = mom.solve(model, feed, loads={term: m["termination_resistance_ohm"]})
    d = 4 * math.pi * _u_axis(sol) / mom.pattern_power(sol, 50, 100)[1]
    assert m["directivity_linear"] == pytest.approx(d, rel=0.02)


def test_the_best_angle_is_wider_than_the_alignment_angle():
    """The alignment angle lines each leg's cone up with the axis for an
    UNATTENUATED wave. With the real, decaying current the axial directivity
    peaks a couple of degrees wider."""
    def axial(half):
        _, _, sol = _solve(4.0, 640.0, 1e-4, half, 28)
        return 4 * math.pi * _u_axis(sol) / mom.pattern_power(sol, 50, 100)[1]
    align = axial(None)
    wider = axial(27.0)
    assert 10 * math.log10(wider / align) > 0.2
