"""Subarray-level steering (`otahub.arrays.subarrays`).

With one phase per subarray the pattern factorises exactly into the subarray's
own (broadside) pattern times an array factor at the subarray period; the
quantisation lobes sit at the period's grating-lobe positions, weighted by the
subarray pattern, and the beam pays the subarray pattern's value as scan loss.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from otahub.arrays import elements as E
from otahub.arrays import subarrays as S
from otahub.arrays.planar import planar_array_factor, planar_directivity, rectangular_lattice

NX = NY = 16
D = 0.5
SX = SY = 4


def _setup(scan):
    pos = rectangular_lattice(NX, NY, D)
    g = S.rect_subarray_groups(NX, NY, SX, SY)
    return pos, g, S.subarray_steering(pos, g, None, scan, 0.0, "subarray")


def _sub(u):
    x = math.pi * D * u
    return 1.0 if abs(x) < 1e-12 else abs(math.sin(SX * x) / (SX * math.sin(x)))


def test_groups_and_centres():
    g = S.rect_subarray_groups(NX, NY, SX, SY)
    assert np.bincount(g).tolist() == [SX * SY] * (NX * NY // (SX * SY))
    c = S.subarray_centres(rectangular_lattice(NX, NY, D), g)
    assert np.allclose(np.sort(np.unique(np.round(c[:, 0], 9))), (np.arange(4) - 1.5) * SX * D)
    with pytest.raises(ValueError):
        S.rect_subarray_groups(10, 10, 4, 4)


@pytest.mark.parametrize("scan", [5.0, 20.0])
def test_the_pattern_is_subarray_pattern_times_the_array_of_centres(scan):
    pos, g, a = _setup(scan)
    u = np.linspace(-1, 1, 801)
    th, ph = np.arcsin(np.abs(u)), np.where(u < 0, math.pi, 0.0)
    full = planar_array_factor(pos, a, th, ph)
    sub = planar_array_factor(pos[g == 0] - pos[g == 0].mean(axis=0), np.ones(SX * SY), th, ph)
    cen = S.subarray_centres(pos, g)
    u0 = math.sin(math.radians(scan))
    arr = np.sum(np.exp(2j * math.pi * np.outer(u - u0, cen[:, 0])), axis=1)
    assert np.max(np.abs(full - sub * arr)) < 1e-11 * np.max(np.abs(full))


@pytest.mark.parametrize("scan", [3.0, 8.0, 20.0])
def test_the_quantisation_lobe_is_weighed_by_the_subarray_pattern(scan):
    """At the subarray period's grating-lobe position the array factor of the
    centres is as strong as at the beam, so the lobe-to-beam ratio is the
    subarray pattern's there over its value at the beam: -17.9 dB at 3 degrees,
    -7.9 at 8, and at 20 degrees on 2-wavelength subarrays the lobe is 6.4 dB ABOVE the beam."""
    pos, g, a = _setup(scan)
    u0 = math.sin(math.radians(scan))
    uq = u0 - 1.0 / (SX * D)
    val = lambda u: abs(planar_array_factor(pos, a, math.asin(abs(u)), math.pi if u < 0 else 0.0)) ** 2
    assert 10 * math.log10(val(uq) / val(u0)) == pytest.approx(20 * math.log10(_sub(uq) / _sub(u0)), abs=1e-9)
    if scan == 20.0:
        assert val(uq) > val(u0)


@pytest.mark.parametrize("scan", [3.0, 8.0, 20.0])
def test_subarray_steering_pays_the_subarray_pattern_in_directivity(scan):
    """The beam's intensity falls by S(u0)^2 exactly; the directivity by that to
    within a quarter decibel, the quantisation lobes shifting the radiated power."""
    pos, g, a = _setup(scan)
    a_el = S.subarray_steering(pos, g, None, scan, 0.0, "element")
    d_el = S.directivity_toward(pos, a_el, scan)
    assert d_el == pytest.approx(planar_directivity(pos, np.ones(len(pos)), scan), rel=1e-9)
    loss = 10 * math.log10(S.directivity_toward(pos, a, scan) / d_el)
    assert loss == pytest.approx(20 * math.log10(_sub(math.sin(math.radians(scan)))), abs=0.25)


def test_directivity_toward_takes_an_element():
    pos, g, a = _setup(0.0)
    assert S.directivity_toward(pos, a, 0.0, element=E.cosine(1.0)) == pytest.approx(
        planar_directivity(pos, np.ones(len(pos)), element=E.cosine(1.0)), rel=1e-9)
