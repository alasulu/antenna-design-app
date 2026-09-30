"""Circular apertures and thinned arrays (`otahub.arrays.layouts`).

Taylor's circular distribution is checked by transforming it back: its Hankel
transform must be the closed-form pattern it was built from. Sampled on a
lattice clipped to a circle it must hold the design level, which depends on
taking the aperture radius right. A thinned array's mean pattern and mean
directivity are exact expectations, checked by Monte Carlo over seeded
thinnings.
"""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import j0

from otahub.arrays import elements as E
from otahub.arrays import layouts as Ly
from otahub.arrays.planar import (planar_array_factor, planar_directivity, planar_pattern_cut,
                                  rectangular_lattice, triangular_lattice)


def _peak_sidelobe_db(pos, w, phi_deg):
    U = planar_pattern_cut(pos, w, phi_deg, n_theta=8001).U[:, 0]
    i0 = int(np.argmax(U))
    k = i0
    while k + 1 < len(U) and U[k + 1] < U[k]:
        k += 1
    j = i0
    while j - 1 >= 0 and U[j - 1] < U[j]:
        j -= 1
    return 10 * math.log10(max(U[k:].max(), U[:j + 1].max()) / U[i0])


def test_ring_and_rings():
    r = Ly.ring(12, 2.0)
    assert np.allclose(np.hypot(r[:, 0], r[:, 1]), 2.0)
    rings = Ly.concentric_rings([0.5, 1.0, 1.5], 0.5)
    assert len(rings) == 1 + 6 + 12 + 18
    clipped = Ly.clip_to_circle(rectangular_lattice(9, 9, 0.5), 1.0)
    assert len(clipped) == 13 and np.all(np.hypot(clipped[:, 0], clipped[:, 1]) <= 1.0)


@pytest.mark.parametrize("sll,nbar", [(-25.0, 4), (-30.0, 5), (-40.0, 8)])
def test_the_taylor_distribution_transforms_back_to_its_pattern(sll, nbar):
    g = lambda p: float(Ly.taylor_circular_distribution(np.array([p]), sll, nbar)[0])
    norm = quad(lambda p: g(p) * p, 0, 1, limit=200)[0]
    for u in (0.3, 0.9, 1.5, 2.2, 3.1, 4.4, 6.0):
        got = quad(lambda p: g(p) * j0(math.pi * u * p) * p, 0, 1, limit=400)[0] / norm
        assert got == pytest.approx(float(Ly.taylor_circular_pattern(u, sll, nbar)), abs=1e-10)


@pytest.mark.parametrize("sll,nbar", [(-30.0, 5), (-40.0, 8)])
def test_the_taylor_pattern_holds_its_design_level(sll, nbar):
    u = np.linspace(0.0, 12.0, 24001)
    F = np.abs(Ly.taylor_circular_pattern(u, sll, nbar))
    i = int(np.argmax(np.diff(np.sign(np.diff(F))) > 0)) + 1
    assert 20 * math.log10(F[i:].max()) == pytest.approx(sll, abs=0.6)
    assert Ly.taylor_circular_distribution(np.array([0.0]), sll, nbar)[0] == pytest.approx(1.0)


@pytest.mark.parametrize("sll,nbar", [(-30, 5), (-40, 8)])
def test_the_taylor_pattern_is_finite_at_the_displaced_nulls(sll, nbar):
    """At u = mu_n the uniform-disc factor and its divisor both vanish; the
    pattern used to return -inf or nan there. Against the Hankel transform."""
    from scipy.special import jn_zeros
    g = lambda p: float(Ly.taylor_circular_distribution(np.array([p]), sll, nbar)[0])
    norm = quad(lambda p: g(p) * p, 0, 1, epsabs=0, epsrel=1e-12, limit=200)[0]
    for mu in jn_zeros(1, nbar - 1) / math.pi:
        hank = quad(lambda p: g(p) * j0(math.pi * mu * p) * p, 0, 1,
                    epsabs=0, epsrel=1e-12, limit=200)[0] / norm
        assert float(Ly.taylor_circular_pattern(np.array([mu]), sll, nbar)[0]) == pytest.approx(
            hank, rel=1e-9, abs=1e-14)


@pytest.mark.parametrize("lattice", ["rect", "tri"])
def test_a_sampled_circular_taylor_array_holds_the_design_with_the_equal_area_radius(lattice):
    d = 0.5
    if lattice == "rect":
        pos, cell = Ly.clip_to_circle(rectangular_lattice(41, 41, d), 8.0), d * d
    else:
        pos, cell = Ly.clip_to_circle(triangular_lattice(41, 48, d), 8.0), math.sqrt(3) / 2 * d * d
    w = Ly.circular_taylor(pos, sidelobe_db=-40.0, nbar=8, cell_area=cell)
    worst = max(_peak_sidelobe_db(pos, w, phi) for phi in (0.0, 30.0, 45.0))
    assert worst < -38.5
    # the outermost element's radius plus half a spacing tapers too little and costs decibels
    r_edge = float(np.hypot(pos[:, 0], pos[:, 1]).max()) + d / 2
    loose = Ly.circular_taylor(pos, r_edge, -40.0, 8)
    assert max(_peak_sidelobe_db(pos, loose, phi) for phi in (0.0, 30.0, 45.0)) > worst + 2.0


def _thinned_case():
    pos = Ly.clip_to_circle(rectangular_lattice(41, 41, 0.5), 8.0)
    return pos, Ly.circular_taylor(pos, sidelobe_db=-35.0, nbar=6, cell_area=0.25)


def test_the_mean_thinned_pattern_is_the_taper_on_a_floor():
    pos, dens = _thinned_case()
    psi = np.linspace(-math.pi / 2, math.pi / 2, 361)
    th, ph = np.abs(psi), np.where(psi < 0, math.pi, 0.0)
    want = Ly.thinned_expected_power(pos, dens, th, ph)
    acc = np.zeros_like(want)
    M = 200
    for s in range(M):
        keep = Ly.thin(dens, seed=s)
        acc += np.abs(planar_array_factor(pos[keep], np.ones(keep.sum()), th, ph)) ** 2
    side = np.abs(psi) > 0.2
    rel = acc[side] / M / want[side] - 1
    assert abs(rel.mean()) < 0.03 and np.sqrt((rel ** 2).mean()) < 0.12       # 200 samples: about 7% noise
    floor = float(np.sum(dens * (1 - dens)))
    assert 10 * math.log10(floor / want.max()) == pytest.approx(-29.3, abs=0.3)


@pytest.mark.parametrize("element,scan", [(None, 0.0), (E.cosine(1.0), 30.0)], ids=["isotropic", "cos-30deg"])
def test_the_expected_thinned_directivity_is_the_mean_of_the_realisations(element, scan):
    pos, dens = _thinned_case()
    want = Ly.thinned_expected_directivity(pos, dens, scan, 0.0, element)
    Ds = np.array([planar_directivity(pos[k], np.ones(k.sum()), scan, 0.0, element=element)
                   for k in (Ly.thin(dens, seed=s) for s in range(120))])
    assert want == pytest.approx(Ds.mean(), abs=3.5 * Ds.std() / math.sqrt(len(Ds)))
    # thinning to about 45% of the elements roughly halves the filled taper's directivity
    assert 0.4 < want / planar_directivity(pos, dens, scan, 0.0, element=element) < 0.7


def test_thin_is_reproducible_and_honours_the_density():
    p = np.linspace(0.0, 1.0, 5001)
    k1, k2 = Ly.thin(p, seed=3), Ly.thin(p, seed=3)
    assert np.array_equal(k1, k2) and not k1[0] and k1[-1]
    assert k1.sum() == pytest.approx(p.sum(), rel=0.03)


def test_the_planar_command_builds_circular_and_thinned_arrays(capsys):
    from otahub.cli.main import main
    assert main(["planar", "--nx", "33", "--ny", "33", "--circle", "8", "--taper", "taylor-circular",
                 "--sll", "35", "--nbar", "6", "--element", "cos"]) == 0
    out = capsys.readouterr().out
    assert out.startswith("797-element rectangular planar array, clipped to a circle of radius 8 lambda")
    sl = float(out.split("first sidelobe")[1].split("dB")[0])
    assert -36.0 < sl < -33.5
    assert main(["planar", "--nx", "33", "--ny", "33", "--circle", "8", "--taper", "taylor-circular",
                 "--sll", "35", "--nbar", "6", "--thin", "7"]) == 0
    out = capsys.readouterr().out
    assert "thinned (seed 7):" in out and "mean sidelobe floor -29.3 dB" in out
    assert main(["planar", "--taper", "bogus"]) == 1
