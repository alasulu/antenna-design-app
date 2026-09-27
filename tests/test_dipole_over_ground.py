"""The horizontal dipole over ground: its feed impedance against the real wire.

The spec gave only the induced-EMF resistance, R11 - R12 at twice the height,
which assumes a sinusoidal current on a vanishing wire - so it ignored the wire
radius it asked for, and had no reactance at all. The dipole and its reversed
image, solved together by the method of moments, present 5-50% more resistance
(14% more a quarter wave up on 1e-4 wavelength wire), and Hallen's equation with
the image in its kernel agrees once the feed gap is the same width. The obvious
shortcut - the free dipole's driving point minus the induced-EMF mutual
impedance - works well up high and fails low, where the image reshapes the
current. The zenith directivity, a pattern property, survives to 0.04 dB.
"""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.special import sici

from otahub.core.constants import ETA0
from otahub.num import mom

C = 2.99792458e8
K = 2 * math.pi


def _pair(h, a, n=64):
    """The half-wave dipole at height h and its image, reversed, fed together (wavelengths)."""
    z = np.linspace(-0.25, 0.25, n + 1)
    real = np.stack([z, 0 * z, 0 * z + h], axis=1)
    image = np.stack([z, 0 * z, 0 * z - h], axis=1)[::-1]
    return mom.solve(mom.WireModel([mom.Wire(real, a), mom.Wire(image, a)]), (n // 2 - 1, (n - 1) + n // 2 - 1))


def _mom(h, a):
    s = _pair(h, a)
    return complex(1 / s.currents[31])


def _hallen_over_ground(a, h, n=64, nq=48):
    """Hallen's equation for the half-wave dipole; the image enters the kernel as -G(R')."""
    hh = 0.25
    z = np.linspace(-hh, hh, n + 1)
    d = z[1] - z[0]
    nb = n - 1
    zm = np.concatenate([z[1:n], [hh]])
    x, w = np.polynomial.legendre.leggauss(nq)
    A = np.zeros((nb + 1, nb + 1), dtype=complex)
    for j in range(nb):
        c = z[j + 1]
        for lo, hi, rise in ((c - d, c, True), (c, c + d, False)):
            s = 0.5 * (hi - lo) * (x + 1.0) + lo
            f = (s - (c - d)) / d if rise else ((c + d) - s) / d
            dz = zm[:, None] - s[None, :]
            R, Ri = np.sqrt(dz ** 2 + a * a), np.sqrt(dz ** 2 + 4 * h * h)
            g = np.exp(-1j * K * R) / (4 * math.pi * R) - np.exp(-1j * K * Ri) / (4 * math.pi * Ri)
            A[:, j] += 0.5 * (hi - lo) * (g * f[None, :] * w[None, :]).sum(axis=1)
    A[:, nb] = (1j / ETA0) * np.cos(K * zm)
    sol = np.linalg.solve(A, -(1j / ETA0) * 0.5 * np.sin(K * np.abs(zm)))
    return complex(1.0 / sol[n // 2 - 1])


def _z12(d):
    """Induced-EMF mutual impedance of side-by-side half-wave dipoles (Balanis 8-72)."""
    u0, u1, u2 = K * d, K * (math.hypot(d, 0.5) + 0.5), K * (math.hypot(d, 0.5) - 0.5)
    si, ci = zip(*(sici(u) for u in (u0, u1, u2)))
    return ETA0 / (4 * math.pi) * complex(2 * ci[0] - ci[1] - ci[2], -(2 * si[0] - si[1] - si[2]))


def _spec(registry, h, a, L=0.5):
    return registry["dipole_over_ground"].synthesize(f0=C, h_over_lambda=h, L_over_lambda=L, aw=a).metrics


@pytest.mark.parametrize("h,a", [(0.06, 3e-5), (0.3, 5e-4), (0.9, 2e-3)])
def test_hallen_with_the_image_meets_the_mom(h, a):
    zm, zh = _mom(h, a), _hallen_over_ground(a, h)
    assert zh.real == pytest.approx(zm.real, rel=2e-3) and zh.imag == pytest.approx(zm.imag, rel=3e-3)


@pytest.mark.parametrize("h,a", [(0.052, 1.5e-5), (0.11, 7e-4), (0.18, 2.4e-3), (0.33, 6e-5), (0.72, 1.1e-3), (1.7, 3e-4)])
def test_the_spec_is_the_dipole_and_its_image_solved_together(registry, h, a):
    z, m = _mom(h, a), _spec(registry, h, a)
    assert m["input_resistance_driving_point_ohm"] == pytest.approx(z.real, rel=6e-3)
    assert m["input_reactance_driving_point_ohm"] == pytest.approx(z.imag, rel=4e-3)


@pytest.mark.parametrize("a", [1e-5, 4e-4, 2.7e-3])
def test_the_self_term_is_the_free_dipole(registry, a):
    m = _spec(registry, 0.25, a)
    z = mom.input_impedance(mom.dipole(0.5, a, 64))
    assert abs(m["self_impedance_driving_point_ohm"] - z) < 0.005
    hw = registry["half_wave_dipole"].synthesize(f0=C, aw=a).metrics
    assert m["self_impedance_driving_point_ohm"].real == pytest.approx(hw["input_resistance_driving_point_ohm"], rel=0.01)


def test_the_induced_emf_resistance_reads_low(registry):
    m = _spec(registry, 0.25, 1e-4)
    assert m["input_resistance_ohm"] / _mom(0.25, 1e-4).real == pytest.approx(0.878, abs=0.003)
    for h in (0.07, 0.4, 1.3):
        for a in (1e-5, 2.7e-3):
            z = _mom(h, a)
            assert _spec(registry, h, a)["input_resistance_ohm"] < 0.96 * z.real, (h, a)
    fat = _mom(0.25, 2.7e-3).real
    assert _spec(registry, 0.25, 2.7e-3)["input_resistance_ohm"] < 0.72 * fat


@pytest.mark.parametrize("a,low,high", [(1e-5, 1.7, 0.035), (2.7e-3, 3.7, 0.135)])
def test_the_free_dipole_minus_the_mutual_impedance_fails_low(registry, a, low, high):
    hw = registry["half_wave_dipole"].synthesize(f0=C, aw=a).metrics
    z11 = hw["input_resistance_driving_point_ohm"] + 1j * hw["input_reactance_driving_point_ohm"]
    assert (z11 - _z12(0.1)).real / _mom(0.05, a).real > low
    errs = [(z11 - _z12(2 * h)).real / _mom(h, a).real - 1 for h in (0.2, 0.3, 0.5, 0.75, 1.2)]
    assert max(abs(e) for e in errs) < high


@pytest.mark.parametrize("h", [0.05, 0.25])
@pytest.mark.parametrize("a", [1e-5, 2.7e-3])
def test_the_zenith_directivity_holds_for_the_real_current(registry, h, a):
    """Twice the pair's free-space directivity overhead: the upper half carries half its power."""
    got = 10 * math.log10(2 * mom.directivity_towards(_pair(h, a), 0.0, 0.0, 64, 64))
    assert 0.0 < got - _spec(registry, h, a)["zenith_directivity_dbi"] < 0.05


@pytest.mark.parametrize("h,a,L", [(0.25, 1e-4, 0.48), (0.04, 1e-4, 0.5), (2.1, 1e-4, 0.5), (0.25, 5e-3, 0.5), (0.25, 5e-6, 0.5)])
def test_outside_its_range_the_driving_point_is_nan(registry, h, a, L):
    m = _spec(registry, h, a, L)
    assert math.isnan(m["input_resistance_driving_point_ohm"]) and math.isnan(m["input_reactance_driving_point_ohm"])
    assert math.isfinite(m["input_resistance_ohm"])
