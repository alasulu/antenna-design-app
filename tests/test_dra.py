"""Dielectric resonator antennas, checked against their complex natural frequencies.

Two findings shaped this file. The hemispherical DRA took its resonance and Q
from the peak and half-power width of the Mie coefficient |b1| on the real
frequency axis; at Q near 9 the scattering background drags that peak off the
pole, so the resonance sat 1.2% high and Q 12% low. And the rectangular DRA
solved the dielectric-waveguide transcendental along its HEIGHT - the mode of an
isolated resonator whose magnetic dipole is vertical, which a ground plane
shorts out. It matched the physical mode only where the imaged brick happens to
be a cube, which the default is.

Every reference here comes from `otahub.num.dra`: the exact sphere pole, and a
3-D FDTD ringdown that is itself checked against that pole before it is trusted
with a brick.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import brentq

from otahub.num import dra

F0 = 1e10
K0 = 2 * math.pi * F0 / 2.99792458e8


# ------------------------------------------------------------ the exact pole

@pytest.mark.parametrize("eps_r", [6.0, 10.0, 25.0, 50.0, 100.0])
def test_sphere_pole_solves_the_characteristic_equation(eps_r):
    x = dra.sphere_pole(eps_r)
    m = math.sqrt(eps_r)
    scale = abs(m * dra._dpsi(m * x) * dra._xi(x))
    assert abs(dra._te1(x, m)) < 1e-10 * scale
    assert x.imag > 0          # decays for exp(+j w t)


def test_sphere_pole_tends_to_one_internal_half_wavelength():
    """sqrt(eps_r)*k0*a -> pi from below: the high-permittivity limit of
    m psi'(m x)/psi(m x) = xi'(x)/xi(x) is d/dz[z psi1(z)] = z sin z = 0."""
    ratios = [math.sqrt(e) * dra.sphere_pole(e).real / math.pi for e in (10.0, 100.0, 1000.0, 10000.0)]
    assert all(b > a for a, b in zip(ratios, ratios[1:]))
    assert ratios[-1] == pytest.approx(1.0, abs=2e-3)


def _b1(x, m):
    """Mie magnetic-dipole scattering coefficient on the real axis (Bohren & Huffman 4.53)."""
    psi, dpsi, xi, dxi = dra._psi, dra._dpsi, dra._xi, dra._dxi
    return (psi(m * x) * dpsi(x) - m * psi(x) * dpsi(m * x)) / (psi(m * x) * dxi(x) - m * xi(x) * dpsi(m * x))


def test_the_scattering_peak_is_not_the_pole_at_low_q():
    """What the hemispherical spec used to quote, reproduced and shown to be off.
    The |b1| peak and its half-power width recover the pole only as Q grows."""
    m = math.sqrt(10.0)
    x = np.linspace(0.85, 1.05, 200001)
    p = np.abs(_b1(x, m)) ** 2
    peak = x[np.argmax(p)]
    band = x[p >= p.max() / 2]
    q_hp = peak / (band[-1] - band[0])
    k0a, q = dra.sphere_resonance(10.0)
    assert peak == pytest.approx(0.95112, abs=2e-5)
    assert q_hp == pytest.approx(8.039, abs=2e-3)
    assert peak / k0a - 1 > 0.01           # 1.2% high
    assert q_hp / q - 1 < -0.10            # 12% low


@pytest.mark.parametrize("eps_r", [6.0, 8.0, 10.0, 14.0, 20.0, 33.0, 50.0, 75.0, 100.0])
def test_hemisphere_follows_the_exact_pole(registry, eps_r):
    d = registry["hemispherical_dra"].synthesize(f0=F0, eps_r=eps_r)
    k0a, q = dra.sphere_resonance(eps_r)
    assert d.get("k0a") == pytest.approx(k0a, rel=2.5e-3)
    assert d.metrics["radiation_q"] == pytest.approx(q, rel=2.5e-3)


# ------------------------------------------------ the FDTD, checked on the pole

@pytest.mark.slow
@pytest.mark.parametrize("eps_r", [10.0, 50.0])
def test_fdtd_ringdown_reproduces_the_sphere_pole(eps_r):
    """A staircased hemisphere 10 cells in radius, no Mie theory anywhere in
    the solver. This is what licenses the brick results below."""
    k0a, q = dra.sphere_resonance(eps_r)
    k_fd, q_fd = dra.hemisphere(eps_r, 10.0)
    assert k_fd == pytest.approx(k0a, rel=3e-3)
    assert q_fd == pytest.approx(q, rel=0.02)


# ------------------------------------------------------------- the brick

DATA = json.loads((Path(__file__).parent / "data" / "rectangular_dra_fdtd.json").read_text())


def _brick(registry, eps_r, aw, aL):
    return registry["rectangular_dra"].synthesize(f0=F0, eps_r=eps_r, aspect_wd=aw, aspect_Ld=aL)


def _dwm_k0d(eps_r, aw, aL):
    """The dielectric waveguide model solved by bisection, transcendental along
    the dipole (L), magnetic walls across w and on the imaged top face."""
    P = (math.pi / 2) ** 2 * aL ** 2 * (1 / aw ** 2 + 0.25)
    g = lambda v: v * math.tan(v) - math.sqrt(max(((eps_r - 1) * P - v * v) / eps_r, 0.0))
    v = brentq(g, 1e-12, math.pi / 2 - 1e-12)
    return 2 * math.sqrt((P + v * v) / eps_r) / aL


@pytest.mark.parametrize("group", ["train", "test"])
def test_brick_resonance_and_q_are_the_recorded_fdtd_ringdowns(registry, group):
    """Every ringdown behind the fit, and every held-out one, reproduced by the
    spec. A ringdown is (eps_r, A along the dipole, B across it, k0*d, Q)."""
    worst_k = worst_q = 0.0
    for r in DATA[group]:
        d = _brick(registry, r["eps"], r["B"], r["A"])
        worst_k = max(worst_k, abs(d.get("k0d") / r["k0d"] - 1))
        worst_q = max(worst_q, abs(d.metrics["radiation_q"] / r["Q"] - 1))
    assert worst_k < DATA["tolerance"]["k0d"], worst_k
    assert worst_q < DATA["tolerance"]["Q"], worst_q


def test_the_other_mode_is_the_same_brick_turned_round(registry):
    """The mode along w of an (aw, aL) brick is the designed mode of the
    (aL, aw) brick: same resonator, so the sizes must scale by the frequency."""
    for eps_r, aw, aL in ((10.0, 3.0, 1.5), (25.0, 1.2, 2.6), (6.0, 2.0, 1.0)):
        d = _brick(registry, eps_r, aw, aL)
        turned = _brick(registry, eps_r, aL, aw)
        ratio = d.metrics["resonant_frequency_other_mode_hz"] / F0
        assert turned.get("d") / d.get("d") == pytest.approx(ratio, rel=1e-12)
    square = _brick(registry, 10.0, 2.0, 2.0)
    assert square.metrics["resonant_frequency_other_mode_hz"] == pytest.approx(F0, rel=1e-12)


@pytest.mark.parametrize("eps_r,aw,aL", [(6.0, 1.0, 3.0), (10.0, 2.0, 2.0), (10.0, 3.0, 1.5),
                                         (22.0, 1.4, 2.3), (50.0, 3.0, 1.0), (50.0, 1.0, 1.0)])
def test_dwm_comparison_is_the_bisected_transcendental(registry, eps_r, aw, aL):
    d = _brick(registry, eps_r, aw, aL)
    assert d.metrics["resonance_residual"] == pytest.approx(0.0, abs=2e-3)
    f_dwm = F0 * _dwm_k0d(eps_r, aw, aL) / d.get("k0d")
    assert d.metrics["resonant_frequency_dwm_hz"] == pytest.approx(f_dwm, rel=2e-4)


def test_the_old_orientation_was_a_mode_the_ground_plane_forbids(registry):
    """The previous spec ran the transcendental along the HEIGHT: the even mode
    of the imaged 2d-tall resonator, whose tangential E is largest exactly where
    the ground plane sits. It agreed with the physical mode only when the imaged
    brick is a cube, w = L = 2d - the default. On a tall brick it was 14% out."""
    def old_k0d(eps_r, aw, aL):
        P = math.pi ** 2 * (1 / aw ** 2 + 1 / aL ** 2)
        u = brentq(lambda u: u * math.tan(u) - math.sqrt(max(((eps_r - 1) * P - u * u) / eps_r, 0.0)),
                   1e-12, math.pi / 2 - 1e-12)
        return math.sqrt((P + u * u) / eps_r)
    assert old_k0d(10.0, 2.0, 2.0) == pytest.approx(_dwm_k0d(10.0, 2.0, 2.0), rel=1e-9)
    tall = _brick(registry, 8.0, 1.0, 1.0)
    assert old_k0d(8.0, 1.0, 1.0) / tall.get("k0d") - 1 > 0.12


@pytest.mark.slow
@pytest.mark.parametrize("eps_r,A,B,height", [(10.0, 2.0, 2.0, 12), (10.0, 1.5, 3.0, 12),
                                              (10.0, 3.0, 1.5, 12), (40.0, 2.0, 2.0, 8)])
def test_brick_fdtd_reproduces_live(registry, eps_r, A, B, height):
    """The arbiter run live on four bricks, dipole along A: the default, both
    modes of an asymmetric brick, and a high-permittivity one on a coarse grid."""
    d = _brick(registry, eps_r, B, A)
    k0d, q = dra.brick(eps_r, int(A * height / 2), int(B * height / 2), height, d.get("k0d"))
    assert d.get("k0d") == pytest.approx(k0d, rel=4e-3)
    assert d.metrics["radiation_q"] == pytest.approx(q, rel=0.02)
