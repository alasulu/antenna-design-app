"""The multi-turn small loop: its inductance, and where its small-loop laws stop.

Two findings pinned here. The spec tuned with Wheeler's current-sheet formula,
good for coils longer than about 0.8 of their radius, while its own close
winding makes them a few tenths of that: against N coaxial rings (Maxwell's
mutual inductance, carried here) Wheeler reads up to 52% low there. And every
electrical law it carries assumes a uniform current, which holds only while the
whole winding is electrically small - the method of moments on closed windings
puts the limit near 0.05 wavelengths of wire, where the spec's own 50 ohm design
had asked for five wavelengths.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.special import ellipe, ellipk

from otahub.num import mom

MU0 = 4e-7 * math.pi
DATA = json.loads((Path(__file__).parent / "data" / "multiturn_loop_mom.json").read_text())


def _rings(a, b, N, pitch):
    def mutual(d):
        k2 = 4 * a * a / (4 * a * a + d * d)
        k = math.sqrt(k2)
        return MU0 * a * ((2 / k - k) * ellipk(k2) - (2 / k) * ellipe(k2))

    L = N * MU0 * a * (math.log(8 * a / b) - 2.0)
    for i in range(N):
        for j in range(i + 1, N):
            L += 2 * mutual(pitch * (j - i))
    return L


def _wheeler(a, N, l):
    return MU0 * N ** 2 * math.pi * a ** 2 / (l + 0.9 * a)


@pytest.mark.parametrize("a,b,N,l", [(0.05, 0.0005, 40, 0.1), (0.02, 0.0002, 50, 0.05)])
def test_rings_meet_wheeler_on_long_coils(a, b, N, l):
    assert _wheeler(a, N, l) == pytest.approx(_rings(a, b, N, l / N), rel=0.01)


@pytest.mark.parametrize("a,b,N", [(0.3, 0.001, 10), (0.12, 0.0015, 6), (0.5, 0.0005, 4)])
def test_wheeler_reads_low_on_the_specs_close_winding(a, b, N):
    l = 2.2 * N * b
    assert _wheeler(a, N, l) / _rings(a, b, N, 2.2 * b) < 0.92


@pytest.mark.parametrize("a,b,N,pf", [(0.07, 0.0004, 7, 2.2), (0.4, 0.0009, 18, 3.5), (0.02, 0.0001, 33, 5.0)])
def test_the_spec_inductance_is_the_ring_models(registry, a, b, N, pf):
    d = registry["multiturn_small_loop"].synthesize(f0=1e6, C=2 * math.pi * a, b=b, N=N, l_coil=N * pf * b,
                                                    C_over_lambda=0.01, Rin_target=1.0)
    assert d.metrics["inductance_H"] == pytest.approx(_rings(a, b, N, pf * b), rel=0.03)
    assert d.get("C_tune") == pytest.approx(1 / ((2 * math.pi * 1e6) ** 2 * d.metrics["inductance_H"]), rel=1e-9)


@pytest.mark.slow
def test_the_mom_meets_the_rings_well_below_self_resonance():
    """A closed 4-turn loop, 0.04 wavelengths of wire: X/omega against the rings."""
    lam, a, b, N, pitch = 2 * math.pi * 0.05 / 0.01, 0.05, 0.0005, 4, 0.003
    n = N * 32
    t = np.linspace(0, 2 * math.pi * N, n + 1)
    helix = np.stack([a * np.cos(t), a * np.sin(t), pitch * t / (2 * math.pi)], axis=1) / lam
    s = np.linspace(0, 1, 9)[:, None]
    lead = helix[-1][None, :] + (helix[0] - helix[-1])[None, :] * s
    sol = mom.solve(mom.WireModel([mom.Wire(helix, b / lam), mom.Wire(lead, b / lam)]), (n - 1) + 3)
    L = sol.input_impedance.imag / (2 * math.pi * 2.99792458e8 / lam)
    lead_L = MU0 / (2 * math.pi) * N * pitch * (math.log(2 * N * pitch / b) - 1)
    assert L == pytest.approx(_rings(a, b, N, pitch) + lead_L, rel=0.01)


def test_one_turn_of_the_mom_model_is_the_loop_solution():
    for s in DATA["single"]:
        assert s["mom"][0] == pytest.approx(s["loop_modal"][0], rel=0.03)
        assert s["mom"][1] == pytest.approx(s["loop_modal"][1], rel=0.01)


def test_the_uniform_current_laws_hold_to_about_005_wavelengths_of_wire():
    for r in DATA["sweep"]:
        dev = r["R"] / r["Rr"] - 1
        if r["wire"] <= 0.05 + 1e-9:
            assert abs(dev) < 0.16, r
        if r["wire"] >= 0.1 - 1e-9:
            assert dev > 0.3, r


def test_past_the_limit_the_spec_refuses(registry):
    """The 50 ohm target at C = 0.1 lambda asks for 51 turns and 5.1 wavelengths of
    wire: the uniform-current law still says 50 ohm, the design metrics say NaN."""
    d = registry["multiturn_small_loop"].synthesize(f0=1e7, C_over_lambda=0.1, Rin_target=50.0, b=0.0005)
    m = d.metrics
    assert d.get("N") == 51 and m["winding_wire_over_lambda"] == pytest.approx(5.1)
    assert m["radiation_resistance_uniform_current_ohm"] == pytest.approx(51 ** 2 * 20 * math.pi ** 2 * 1e-4, rel=1e-9)
    for k in ("radiation_resistance_ohm", "radiation_efficiency", "quality_factor", "fractional_bandwidth_vswr2"):
        assert math.isnan(m[k]), k
    assert m["max_radiation_resistance_small_ohm"] < 0.01
