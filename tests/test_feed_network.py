"""Feed networks coupled to the wire solver, and the LPDA they made reachable.

A non-radiating network - a feeder, a phasing line, a stub - joins the wires in
admittance form, as NEC's TL cards do. Three checks license it, each against
something that shares no code with the coupling:

1. a dipole behind a length of line must reproduce the textbook impedance
   transformation exactly;
2. two dipoles tied in parallel through a stiff network must match the plain
   solver driving both gaps at once;
3. a lossless feeder must deliver exactly the power the pattern integral says
   is radiated.

The LPDA then gets two checks of its own that need no reference at all: it must
fire towards its apex, and must stop doing so when the feeder is NOT transposed.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from otahub.num import mom

pytestmark = pytest.mark.slow   # full-wave solves; skip with -m 'not slow'


# ---------------------------------------------------------- the coupling

@pytest.mark.parametrize("z0,bl", [(50.0, 0.3), (100.0, 1.0), (300.0, 2.2)])
def test_a_line_transforms_a_dipole_exactly(z0, bl):
    m = mom.dipole(0.47, 1e-3, 40)
    zl = mom.input_impedance(m)
    ns = mom.solve_network(m, [m.n_basis // 2], mom.tl_admittance(z0, bl), 1)
    t = math.tan(bl)
    want = z0 * (zl + 1j * z0 * t) / (z0 + 1j * zl * t)
    assert abs(ns.input_impedance - want) / abs(want) < 1e-10


def test_a_stiff_parallel_tie_matches_driving_both_gaps():
    d1 = mom.dipole(0.47, 1e-3, 30).wires[0]
    m = mom.WireModel([d1, mom.Wire(d1.nodes + np.array([0.15, 0.0, 0.0]), 1e-3)])
    p1, p2 = 14, m.n_basis // 2 + 14
    V = np.zeros(m.n_basis, dtype=complex)
    V[p1] = V[p2] = 1.0
    I = np.linalg.solve(mom.impedance_matrix(m), V)
    direct = 1.0 / (I[p1] + I[p2])
    g = 1e9
    ns = mom.solve_network(m, [p1, p2], np.array([[g, -g], [-g, g]], dtype=complex), 0)
    assert ns.input_impedance == pytest.approx(direct, rel=1e-5)


def test_a_lossless_feeder_conserves_power():
    m = mom.dipole(0.47, 1e-3, 40)
    ns = mom.solve_network(m, [m.n_basis // 2], mom.tl_admittance(75.0, 1.3), 1)
    assert mom.radiated_power(ns.solution, 70, 70) == pytest.approx(
        ns.generator_power, rel=1e-4)


def test_the_network_must_cover_every_port():
    m = mom.dipole(0.47, 1e-3, 20)
    with pytest.raises(ValueError):
        mom.solve_network(m, [3, 5], np.zeros((1, 1), dtype=complex), 0)


# ------------------------------------------------------------- the LPDA

def _lpda(tau=0.9, transposed=True, z0=100.0, back=5, front=9, scale=1.0):
    sigma = 0.243 * tau - 0.051
    return mom.solve_network(*mom.lpda_model(
        tau, sigma, 0.5 * tau ** (-back) * scale, back + front + 1,
        feeder_z0=z0, transposed=transposed))


def _fwd_back(ns):
    sol = ns.solution
    _, total, _ = mom.pattern_power(sol, 60, 80)

    def u(ph):
        a, b = mom.far_field(sol, np.array([math.pi / 2]), np.array([ph]))
        return abs(a[0]) ** 2 + abs(b[0]) ** 2
    return 4 * math.pi * u(math.pi) / total, u(math.pi) / u(0.0)


def test_an_lpda_conserves_power_through_its_feeder():
    ns = _lpda()
    assert mom.radiated_power(ns.solution, 60, 80) == pytest.approx(
        ns.generator_power, rel=1e-3)


def test_it_fires_towards_the_apex_and_needs_the_transposition():
    """What the transposition actually buys, measured over one log-period.

    The spec used to say that without the 180-degree reversal the array "fires
    backwards". It does not, reliably: it loses its front-to-back ratio
    altogether - roughly bidirectional, within a few dB either way - and its
    input impedance swings tenfold inside a single period. It stops being a
    frequency-independent antenna, which is the whole point of one.
    """
    def over_period(transposed):
        fbs, zs = [], []
        for i in range(4):
            ns = _lpda(transposed=transposed, scale=0.9 ** (-i / 4))
            fbs.append(10 * math.log10(_fwd_back(ns)[1]))
            zs.append(abs(ns.input_impedance))
        return fbs, zs

    fb_t, z_t = over_period(True)
    fb_n, z_n = over_period(False)
    assert min(fb_t) > 10.0, fb_t
    assert max(z_t) / min(z_t) < 1.4, z_t
    assert max(abs(f) for f in fb_n) < 6.0, fb_n
    assert max(z_n) / min(z_n) > 3.0, z_n


def test_the_feeder_sets_the_input_resistance(registry):
    """The design relation the spec now carries, against the solver."""
    for z0 in (50.0, 100.0, 200.0):
        got = np.mean([_lpda(z0=z0, scale=0.9 ** (-i / 4)).input_impedance.real
                       for i in range(4)])
        d = registry["lpda"].synthesize(f_low=1e8, f_high=3e8, tau=0.9,
                                        feeder_z0=z0, length_over_diameter=125.0)
        assert d.metrics["input_resistance_ohm"] == pytest.approx(got, rel=0.09)
