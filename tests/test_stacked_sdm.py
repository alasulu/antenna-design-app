"""The spectral MoM with a second layer and a probe feed - the stacked patch's solver.

`patch_sdm.StackedPatch` puts a parasitic patch on a second layer above the driven
one; the Green's function couples them through the two-layer transmission-line
network `layered_z`. `probe_vector` couples the patch currents to a probe from the
ground to the driven patch. These are the checks that need no FDTD: the network
against the single slab it must reduce to, the stacked matrix against the single
patch, the far field against the spectral power, and the probe coupling against its
own truncation and the cavity model's position law.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from otahub.num import patch_sdm as sdm

C = 2.99792458e8
F0 = 1e9
LAM = C / F0
K0 = 2 * math.pi / LAM
KR = np.array([0.3, 1.2, 1.6 + 0.1j, 5.0, 50.0]) * K0      # through, near and past the poles, on and off the axis
LOWER = (2.2, 0.0128 * LAM, 0.32 * LAM, 0.3968 * LAM)
UPPER = (1.0, 0.032 * LAM, 0.3072 * LAM, 0.3712 * LAM)


def test_air_above_the_driven_layer_is_the_single_slab():
    z = sdm.layered_z(KR, K0, 2.2, 0.0016, 1.0, 0.004)
    tm, te = sdm.RectPatch(2.2, 0.0016, 0.04, 0.05).znode(KR, K0)
    assert np.max(np.abs(z[0][0, 0] / tm - 1)) < 1e-12
    assert np.max(np.abs(z[1][0, 0] / te - 1)) < 1e-12


def test_one_dielectric_throughout_is_the_thicker_slab():
    z = sdm.layered_z(KR, K0, 2.2, 0.0016, 2.2, 0.004)
    tm, te = sdm.RectPatch(2.2, 0.0056, 0.04, 0.05).znode(KR, K0)
    assert np.max(np.abs(z[0][1, 1] / tm - 1)) < 1e-12
    assert np.max(np.abs(z[1][1, 1] / te - 1)) < 1e-12


def test_the_coupling_is_reciprocal():
    """Z21 is found going up from node 1; the same term found coming down from node 2,
    through the shorted layer below, must agree."""
    from otahub.core.constants import EPS0, MU0
    er1, h1, er2, h2 = 2.2, 0.0016, 1.3, 0.004
    z = sdm.layered_z(KR, K0, er1, h1, er2, h2)
    w = K0 * C
    kz = [-1j * np.sqrt(KR * KR - e * K0 * K0 + 0j) for e in (1.0, er1, er2)]
    for t, (Y0, Y1, Y2) in enumerate(((w * EPS0 / kz[0], w * EPS0 * er1 / kz[1], w * EPS0 * er2 / kz[2]),
                                      (kz[0] / (w * MU0), kz[1] / (w * MU0), kz[2] / (w * MU0)))):
        yd1 = -1j * Y1 / np.tan(kz[1] * h1)
        c2, s2 = np.cos(kz[2] * h2), np.sin(kz[2] * h2)
        assert np.max(np.abs(z[t][1, 1] / (c2 + 1j * (yd1 / Y2) * s2) / z[t][1, 0] - 1)) < 1e-12


def test_the_stack_holds_the_single_patch_under_air():
    """With air in the second layer the driven patch's block is the lone patch's matrix."""
    s = sdm.StackedPatch(*LOWER, *UPPER, **sdm.BASIS_FULL)
    lone = sdm.RectPatch(*LOWER, **sdm.BASIS_FULL).Z(F0)
    assert np.max(np.abs(s.Z(F0)[:s.n1, :s.n1] - lone)) < 1e-12 * np.max(np.abs(lone))


def test_the_far_field_carries_the_spectral_power():
    """The two-layer far-field factors, integrated over the sphere, against the space-wave
    part of the reaction - for an arbitrary current on both patches."""
    s = sdm.StackedPatch(*LOWER, *UPPER, **sdm.BASIS_FULL)
    v = np.random.default_rng(3).normal(size=len(s.basis))
    for f in (0.99 * F0, 1.39 * F0):
        p_far, _ = s.far_field(v, f)
        assert p_far == pytest.approx(0.5 * v @ s.R_space(f) @ v, rel=1e-9)


def test_the_probe_coupling_converges():
    s = sdm.StackedPatch(*LOWER, *UPPER, **sdm.BASIS_FULL)
    a, b = (sdm.probe_vector(s, F0, 0.096 * LAM, 0.0005 * LAM, kmax=k) for k in (200.0, 400.0))
    assert np.max(np.abs(a - b)) < 3e-3 * np.max(np.abs(b))
    assert np.max(np.abs(a[s.n1:] - b[s.n1:])) < 1e-6 * np.max(np.abs(b))   # the parasitic's converge at once


def test_a_probe_on_the_centre_line_couples_to_nothing():
    p = sdm.RectPatch(*LOWER, **sdm.BASIS_FULL)
    assert np.all(sdm.probe_vector(p, F0, 0.0, 0.0005 * LAM) == 0)


def test_the_resonant_resistance_follows_the_cavity_law():
    """R at resonance against the probe's distance from the centre: sin^2(pi x/L), the
    cavity model's law, and an edge resistance near its 1/(2G) of about 210 ohm."""
    p = sdm.RectPatch(*LOWER, **sdm.BASIS_FULL)
    f = 1.0040795 * F0
    xs = np.array([0.04, 0.08, 0.12]) * LAM
    v = sdm.probe_vector(p, f, xs, 0.000432 * LAM)
    r = -np.einsum("ij,ij->i", v, np.linalg.solve(p.Z(f), v.T).T).real
    edge = r / np.sin(math.pi * xs / (0.32 * LAM)) ** 2
    assert np.ptp(edge) < 0.04 * edge.mean()
    assert 190 < edge.mean() < 240
