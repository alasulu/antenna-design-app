"""Single-turn small loops against their real current.

The circular and square small-loop specs used the uniform-current laws out to
a third of a wavelength round, claiming the quartic resistance law was within
2.3% there - a comparison with the uniform-current integral, not with a real
loop. Solved with its actual current (the Fourier-series loop solution for the
circle, the method of moments for the square, which meet each other on the
circle), the loop's resistance is already 12-17% above the law at 0.1
wavelengths and several times it at a third: the current's cos(phi) part
radiates like an electric dipole, growing as the size squared.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from otahub.num import loop_modal, mom

C0 = 2.99792458e8


def _mom_circle(Cl, b, nseg=64):
    a = Cl / (2 * math.pi)
    t = 2 * math.pi * np.arange(nseg) / nseg
    m = mom.WireModel([mom.Wire(np.stack([a * np.cos(t), a * np.sin(t), 0 * t], axis=1), b, closed=True)])
    return mom.solve(m, 0).input_impedance


def _mom_square(P, b):
    s = P / 4
    n = max(8, min(40, int(s / (4 * b))))
    n += n % 2
    c = np.array([[-s / 2, -s / 2, 0], [s / 2, -s / 2, 0], [s / 2, s / 2, 0], [-s / 2, s / 2, 0], [-s / 2, -s / 2, 0]])
    pts = [c[0]] + [c[k] + (c[k + 1] - c[k]) * j / n for k in range(4) for j in range(1, n + 1)]
    m = mom.WireModel([mom.Wire(np.array(pts)[:-1], b, closed=True)])
    return mom.solve(m, n // 2 - 1).input_impedance


@pytest.mark.parametrize("Cl", [0.05, 0.1, 0.2])
def test_the_two_loop_solvers_agree(Cl):
    b = Cl / (2 * math.pi) / 100
    z1, z2 = _mom_circle(Cl, b), loop_modal.input_impedance(Cl, b)
    assert z1.real == pytest.approx(z2.real, rel=0.01) and z1.imag == pytest.approx(z2.imag, rel=0.006)


@pytest.mark.parametrize("Cl,ba", [(0.07, 0.004), (0.13, 0.02), (0.24, 0.0015)])
def test_the_circular_loop_driving_point(registry, Cl, ba):
    lam = C0 / 300e6
    d = registry["small_circular_loop"].synthesize(f0=300e6, C=Cl * lam, b=ba * Cl * lam / (2 * math.pi), N=1)
    z = loop_modal.input_impedance(Cl, ba * Cl / (2 * math.pi))
    assert d.metrics["input_resistance_driving_point_ohm"] == pytest.approx(z.real, rel=0.025)
    assert d.metrics["input_reactance_driving_point_ohm"] == pytest.approx(z.imag, rel=0.012)
    assert d.get("C_tune") == pytest.approx(1 / (2 * math.pi * 300e6 * z.imag), rel=0.015)


@pytest.mark.parametrize("P,bs", [(0.08, 0.005), (0.15, 0.02), (0.26, 0.002)])
def test_the_square_loop_driving_point(registry, P, bs):
    lam = C0 / 300e6
    d = registry["small_square_loop"].synthesize(f0=300e6, s=P * lam / 4, b=bs * P * lam / 4, N=1)
    z = _mom_square(P, bs * P / 4)
    assert d.metrics["input_resistance_driving_point_ohm"] == pytest.approx(z.real, rel=0.03)
    assert d.metrics["input_reactance_driving_point_ohm"] == pytest.approx(z.imag, rel=0.015)


def test_the_uniform_current_law_understates_a_tenth_wave_loop_and_fails_at_a_third():
    for ba in (0.001, 0.01, 0.05):
        z = loop_modal.input_impedance(0.1, ba * 0.1 / (2 * math.pi))
        assert 0.11 < z.real / (20 * math.pi ** 2 * 0.1 ** 4) - 1 < 0.18
        z3 = loop_modal.input_impedance(1 / 3, ba / (6 * math.pi))
        assert z3.real / (20 * math.pi ** 2 / 81) > 4.5


def test_past_030_the_driving_point_is_nan(registry):
    lam = C0 / 300e6
    d = registry["small_circular_loop"].synthesize(f0=300e6, C=0.33 * lam, b=0.001, N=1)
    assert math.isnan(d.metrics["input_resistance_driving_point_ohm"]) and math.isnan(d.metrics["radiation_efficiency"])
    assert d.metrics["radiation_resistance_ohm"] == pytest.approx(20 * math.pi ** 2 * 0.33 ** 4, rel=1e-6)
