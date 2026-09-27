"""Two wire archetypes no round had examined, against the wire solvers.

quarter_wave_monopole carried only the induced-EMF pair, 36.54 + j21.26, which
assumes a sinusoidal current on a vanishing wire. half_wave_dipole already
distinguishes that from what a real wire presents at its terminals; by image
theory the monopole over perfect ground presents exactly half its dipole's
driving point, checked here against Hallen's equation solved live.

normal_mode_helix's axial ratio is Kraus's 2 S lambda / C^2, derived for a
uniform current. The spec's own default carries 2.1 wavelengths of wire, so
the current is anything but uniform - and the MoM shows the ratio holds anyway,
since every turn's loop and dipole parts scale with the same local current.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from otahub.num import mom
from otahub.num.hallen import hallen_dipole

C = 2.99792458e8


@pytest.mark.parametrize("a", [1e-5, 3e-4, 2e-3])
def test_the_monopole_is_half_its_dipole(registry, a):
    aw = a * C / 3e8
    mono = registry["quarter_wave_monopole"].synthesize(f0=3e8, aw=aw).metrics
    dip = registry["half_wave_dipole"].synthesize(f0=3e8, aw=aw).metrics
    for k in ("input_resistance_driving_point_ohm", "input_reactance_driving_point_ohm"):
        assert mono[k] == pytest.approx(0.5 * dip[k], rel=1e-12)
    assert mono["radiation_resistance_ohm"] == pytest.approx(0.5 * dip["radiation_resistance_ohm"], rel=1e-9)


def test_the_monopole_driving_point_against_hallen(registry):
    a = 3e-4
    z = hallen_dipole(0.5, a, 120) / 2
    m = registry["quarter_wave_monopole"].synthesize(f0=3e8, aw=a * C / 3e8).metrics
    assert m["input_resistance_driving_point_ohm"] == pytest.approx(z.real, rel=0.01)
    assert m["input_reactance_driving_point_ohm"] == pytest.approx(z.imag, rel=0.015)
    assert m["input_resistance_driving_point_ohm"] > m["radiation_resistance_ohm"] * 1.08


def _helix_axial_ratio(D, S, N, a=1e-3, seg_per_turn=24):
    n = int(N * seg_per_turn)
    n += n % 2
    t = np.linspace(0, 2 * math.pi * N, n + 1)
    nodes = np.stack([D / 2 * np.cos(t), D / 2 * np.sin(t), S * t / (2 * math.pi) - S * N / 2], axis=1)
    sol = mom.solve(mom.WireModel([mom.Wire(nodes, a)]), n // 2 - 1)
    ph = np.linspace(0, 2 * math.pi, 37)[:-1]
    et, ep = mom.far_field(sol, np.full_like(ph, math.pi / 2), ph)
    return float(np.median(np.abs(et) / np.abs(ep)))


@pytest.mark.parametrize("D,S,N", [(0.01, 0.002, 20), (0.006, 0.001, 4), (0.01, 0.0016, 3)])
def test_the_helix_axial_ratio_holds_even_on_a_long_wire(registry, D, S, N):
    d = registry["normal_mode_helix"].synthesize(f0=1e9, D_helix=D, S=S, N=N)
    lam = C / 1e9
    got = _helix_axial_ratio(D / lam, S / lam, N)
    assert got == pytest.approx(d.metrics["axial_ratio"], rel=0.015)
