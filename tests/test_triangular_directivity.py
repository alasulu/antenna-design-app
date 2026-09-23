"""Equilateral triangular patch directivity, against the integral it was fitted to.

The triangle's dominant cavity mode is not quoted from a table here: the six
plane waves of magnitude k = 4*pi/(3a) that sit on the hexagonal star span a
two-dimensional null space of the magnetic-wall condition, and the member
symmetric about a median is the one a feed on that median excites. All three
walls radiate. Ez restricted to a straight side is a sum of exponentials, so
each wall's far-field integral is closed form and the only numerical step is
the hemisphere integral.

The test that matters most is the last one: as the patch shrinks it becomes a
horizontal magnetic dipole over a ground plane, so D must go to exactly 3.
Nothing in the derivation was arranged to make that happen.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

A = 1.0
V1 = np.array([0.0, 0.0])
V2 = np.array([A, 0.0])
V3 = np.array([A * 0.5, A * math.sqrt(3) / 2])
K = 4 * math.pi / (3 * A)
KV = K * np.column_stack([np.cos(np.deg2rad(np.arange(0, 360, 60))),
                          np.sin(np.deg2rad(np.arange(0, 360, 60)))])
EDGES = ((V1, V2), (V2, V3), (V3, V1))
CEN = (V1 + V2 + V3) / 3


def _outward(p, q):
    d = q - p
    n = np.array([d[1], -d[0]]) / np.linalg.norm(d)
    return -n if np.dot(n, CEN - (p + q) / 2) > 0 else n


def _walls():
    """The dominant symmetric mode, as twelve exponentials per wall."""
    nout = [_outward(p, q) for p, q in EDGES]
    rows = []
    for (p, q), n in zip(EDGES, nout):
        for s in np.linspace(0.005, 0.995, 60):
            rows.append(1j * (KV @ n) * np.exp(1j * (KV @ (p + s * (q - p)))))
    null = np.linalg.svd(np.array(rows))[2][-2:].conj()
    perm, t = [3, 2, 1, 0, 5, 4], np.array([A, 0.0])
    sym = np.zeros_like(null)
    for i in range(2):
        for j in range(6):
            sym[i, perm[j]] = null[i, j] * np.exp(1j * KV[j] @ t)
    w, v = np.linalg.eig(np.linalg.solve(null.conj() @ null.T, null.conj() @ sym.T))
    c = v[:, np.argmin(np.abs(w - 1.0))] @ null

    def ez(pts):
        return (np.exp(1j * (np.atleast_2d(pts) @ KV.T)) @ c).real

    c = c / max(np.abs(ez(p + np.linspace(0, 1, 201)[:, None] * (q - p))).max()
                for p, q in EDGES)
    out = []
    for (p, q), n in zip(EDGES, nout):
        d = q - p
        base = c * np.exp(1j * (KV @ p))
        out.append((p, d, float(np.linalg.norm(d)), np.array([-n[1], n[0]]),
                    np.concatenate([0.5 * base, 0.5 * np.conj(base)]),
                    np.concatenate([KV @ d, -(KV @ d)])))
    return out


WALLS = _walls()


def _radiation_vector(k0, rx, ry):
    lx = np.zeros(np.shape(rx), complex)
    ly = np.zeros(np.shape(rx), complex)
    for p, d, length, tan, amp, rate in WALLS:
        dr = k0 * (d[0] * rx + d[1] * ry)
        acc = sum(a * np.exp(0.5j * (r + dr)) * np.sinc((r + dr) / (2 * np.pi))
                  for a, r in zip(amp, rate))
        s = length * np.exp(1j * k0 * (p[0] * rx + p[1] * ry)) * acc
        lx += tan[0] * s
        ly += tan[1] * s
    return lx, ly


def _exact(k0a, nth=120, nph=240):
    ct, wct = np.polynomial.legendre.leggauss(nth)
    ct, wct = 0.5 * (ct + 1), 0.5 * wct
    st = np.sqrt(1 - ct ** 2)
    ph = np.arange(nph) * 2 * np.pi / nph
    cp, sp = np.cos(ph)[None, :], np.sin(ph)[None, :]
    lx, ly = _radiation_vector(k0a, st[:, None] * cp, st[:, None] * sp)
    lth = (lx * cp + ly * sp) * ct[:, None]
    lph = -lx * sp + ly * cp
    total = float(wct @ (np.abs(lth) ** 2 + np.abs(lph) ** 2).sum(1)) * (2 * np.pi / nph)
    lx0, ly0 = _radiation_vector(k0a, np.array(0.0), np.array(0.0))
    return 4 * np.pi * float(abs(lx0) ** 2 + abs(ly0) ** 2) / total


@pytest.mark.parametrize("eps_r", [1.0, 2.2, 3.0, 4.4, 6.15, 10.2, 13.0])
def test_triangle_directivity_matches_its_own_integral(eps_r, registry):
    design = registry["triangular_patch"].synthesize(f0=2e9, eps_r=eps_r, h=1.6e-3)
    want = _exact(4 * math.pi / (3 * math.sqrt(eps_r)))
    assert design.metrics["directivity_linear"] == pytest.approx(want, rel=1e-3)


@pytest.mark.parametrize("h", [0.5e-3, 1.6e-3, 3.0e-3])
def test_k0_a_eff_is_fixed_by_permittivity_alone(h, registry):
    """The fringing correction cancels exactly out of k0*a_eff, which is why the
    directivity is a function of one variable and not three."""
    for eps_r in (2.2, 4.4, 10.2):
        design = registry["triangular_patch"].synthesize(f0=3e9, eps_r=eps_r, h=h)
        assert design.metrics["k0_a_eff"] == pytest.approx(
            4 * math.pi / (3 * math.sqrt(eps_r)), rel=1e-9)


def test_directivity_falls_with_permittivity(registry):
    previous = None
    for eps_r in (1.0, 2.2, 4.4, 6.15, 10.2, 13.0):
        value = registry["triangular_patch"].synthesize(
            f0=2e9, eps_r=eps_r, h=1.6e-3).metrics["directivity_linear"]
        if previous is not None:
            assert value < previous, f"eps_r={eps_r} did not fall"
        previous = value


def test_small_patch_limit_is_the_magnetic_dipole(registry):
    """A patch much smaller than a wavelength is a horizontal magnetic dipole
    over a ground plane: D = 3 exactly. The shipped form D = 3 + x^2*g(x) keeps
    that limit outside the fitted range, where a bare polynomial would not."""
    assert _exact(0.00419) == pytest.approx(3.0, abs=1e-5)
    previous = 4.0
    for eps_r in (50.0, 200.0, 1000.0, 1e4):
        value = registry["triangular_patch"].synthesize(
            f0=2e9, eps_r=eps_r, h=1.6e-3).metrics["directivity_linear"]
        assert 3.0 <= value < previous
        previous = value
    assert value == pytest.approx(3.0, abs=5e-3)


def test_triangle_matches_the_rectangle_it_replaces(registry):
    """The spec used to say the triangle runs 'roughly 1 dB below' a rectangular
    patch. It does not: the two agree within 0.1 dB across the whole substrate
    range, which is the actual case for the shape - a third less board area at
    no cost in directivity. Two independently derived patterns agreeing this
    closely is also a real cross-check on both."""
    for eps_r in (1.0, 2.2, 4.4, 10.2, 13.0):
        tri = registry["triangular_patch"].synthesize(
            f0=2e9, eps_r=eps_r, h=1.6e-3).metrics["directivity_linear"]
        rect = registry["rectangular_patch"].synthesize(
            f0=2e9, eps_r=eps_r, h=1.6e-3).metrics["directivity_linear"]
        assert abs(10 * math.log10(tri / rect)) < 0.15, (
            f"eps_r={eps_r}: triangle {tri:.4f} vs rectangle {rect:.4f}")


def test_area_saving_is_real(registry):
    """The reason to accept the awkward feed geometry in the first place."""
    design = registry["triangular_patch"].synthesize(f0=1e9, eps_r=2.2, h=1.6e-3)
    assert 0.55 < design.metrics["area_ratio_vs_rectangular"] < 0.70
