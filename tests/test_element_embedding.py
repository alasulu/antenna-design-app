"""Planar arrays with an element pattern folded in.

`otahub.arrays.elements` turns the element's power pattern into azimuthal
harmonics, which makes the array's power integral a one-dimensional Bessel
integral per separation. Checked against: brute-force integration over the
sphere, the isotropic sinc kernel it must reduce to, closed-form element
directivities (cos^q, a short dipole, a horizontal dipole over ground from its
image), and the aperture limit 4 pi A cos(theta)/lambda^2 that a large array of
ideal (cos theta) elements must reach.
"""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import simpson

from otahub.arrays import elements as E
from otahub.arrays.planar import (planar_array_factor, planar_beam_cut, planar_directivity,
                                  rectangular_lattice, triangular_lattice)


def _brute_power(el, pos, w, scan_t, scan_p, nt=1601, nph=1600):
    th = np.linspace(0, math.pi if not el.half_space else math.pi / 2, nt)
    ph = 2 * math.pi * np.arange(nph) / nph
    T, P = np.meshgrid(th, ph, indexing="ij")
    af = planar_array_factor(pos, w, T, P, scan_t, scan_p)
    f = el(T, P) * np.abs(af) ** 2 * np.sin(T)
    return simpson(f.sum(axis=1) * (2 * math.pi / nph), x=th)


@pytest.mark.parametrize("el", [E.short_dipole("x", 0.25), E.cosine(1.5), E.short_dipole("y")], ids=lambda e: e.name)
def test_the_kernel_matches_brute_force_integration(el):
    for d in ([0.5, 0.0], [0.35, 0.6], [2.0, -1.0]):
        th = np.linspace(0, math.pi / 2 if el.half_space else math.pi, 2001)
        ph = 2 * math.pi * np.arange(2000) / 2000
        T, P = np.meshgrid(th, ph, indexing="ij")
        f = el(T, P) * np.exp(2j * math.pi * np.sin(T) * (d[0] * np.cos(P) + d[1] * np.sin(P))) * np.sin(T)
        want = simpson(f.sum(axis=1) * (2 * math.pi / 2000), x=th)
        assert E.power_kernel(el, np.array([d]))[0] == pytest.approx(want, rel=1e-6, abs=1e-9)


def test_an_isotropic_element_reproduces_the_sinc_formula():
    rng = np.random.default_rng(7)
    pos, w = rng.uniform(-2, 2, (14, 2)), rng.uniform(0.3, 1.0, 14)
    for t, p in ((0, 0), (25, 40), (60, 110)):
        assert planar_directivity(pos, w, t, p, element=E.isotropic()) == pytest.approx(
            planar_directivity(pos, w, t, p), rel=1e-9)


@pytest.mark.parametrize("q", [0.0, 1.0, 1.5, 2.0])
def test_a_cosine_element_on_its_own(q):
    assert E.element_directivity(E.cosine(q)) == pytest.approx(2 * (q + 1), rel=1e-6)


def test_a_short_dipole_on_its_own():
    assert E.element_directivity(E.short_dipole("x")) == pytest.approx(1.5, rel=1e-12)


@pytest.mark.parametrize("h", [0.1, 0.25, 0.4])
def test_a_horizontal_dipole_over_ground_matches_its_image_closed_form(h):
    """Balanis's horizontal dipole above a PEC plane: radiated power in the
    bracket 2/3 - sin(2kh)/(2kh) - cos(2kh)/(2kh)^2 + sin(2kh)/(2kh)^3, so the
    broadside directivity is 4 sin^2(kh) over it."""
    x = 4 * math.pi * h
    bracket = 2 / 3 - math.sin(x) / x - math.cos(x) / x ** 2 + math.sin(x) / x ** 3
    want = 4 * math.sin(2 * math.pi * h) ** 2 / bracket
    assert E.element_directivity(E.short_dipole("x", h)) == pytest.approx(want, rel=1e-9)


@pytest.mark.parametrize("el,scan", [(E.short_dipole("x", 0.25), (25.0, 30.0)), (E.cosine(1.0), (40.0, 0.0)),
                                     (E.short_dipole("y"), (15.0, 70.0))], ids=["dipole-over-ground", "cos", "free-dipole"])
def test_array_directivity_matches_brute_force(el, scan):
    pos = triangular_lattice(3, 3, 0.6)
    w = np.array([1.0, 0.8, 0.6, 0.9, 1.0, 0.7, 0.5, 0.8, 1.0])
    t, p = scan
    af0 = abs(np.sum(w)) ** 2
    want = 4 * math.pi * float(el(math.radians(t), math.radians(p))) * af0 / _brute_power(el, pos, w, t, p)
    assert planar_directivity(pos, w, t, p, element=el) == pytest.approx(want, rel=1e-5)


@pytest.mark.parametrize("scan", [0.0, 30.0, 50.0])
def test_a_large_array_of_ideal_elements_reaches_the_aperture_limit(scan):
    """4 pi A cos(theta)/lambda^2 for a uniform 24 x 24 array at half a
    wavelength; the isotropic-doubled 'ground plane' figure falls 2-4% short."""
    pos = rectangular_lattice(24, 24, 0.5)
    w = np.ones(len(pos))
    limit = 4 * math.pi * 12.0 ** 2 * math.cos(math.radians(scan))
    ideal = planar_directivity(pos, w, scan, 0.0, element=E.cosine(1.0))
    assert ideal == pytest.approx(limit, rel=0.01)
    doubled = planar_directivity(pos, w, scan, 0.0, half_space=True)
    assert 0.02 < 1 - doubled / ideal < 0.045


def test_the_element_pattern_suppresses_a_grating_lobe_at_the_horizon():
    """One-wavelength spacing puts grating lobes on the horizon at broadside;
    ideal elements null them there, isotropic ones spend power on them."""
    pos = rectangular_lattice(8, 8, 1.0)
    w = np.ones(len(pos))
    ideal = planar_directivity(pos, w, element=E.cosine(1.0))
    doubled = planar_directivity(pos, w, half_space=True)
    assert 3.5 < ideal / doubled < 4.3                       # 3.9: 5.9 dB


def test_a_beam_cut_carries_the_element_pattern():
    pos = rectangular_lattice(6, 6, 0.5)
    w = np.ones(len(pos))
    el = E.cosine(2.0)
    bare = planar_beam_cut(pos, w, "scan", 20.0, 0.0)
    shaped = planar_beam_cut(pos, w, "scan", 20.0, 0.0, element=el)
    psi = bare.theta - math.pi / 2
    ratio = shaped.U[:, 0] / np.where(bare.U[:, 0] > 0, bare.U[:, 0], 1.0)
    ok = bare.U[:, 0] > 1e-9 * bare.U.max()
    assert np.allclose(ratio[ok], np.cos(psi[ok]) ** 2, rtol=1e-9)


def test_an_element_and_a_ground_plane_flag_do_not_mix():
    with pytest.raises(ValueError):
        planar_directivity(rectangular_lattice(2, 2, 0.5), np.ones(4), half_space=True, element=E.cosine(1.0))


def test_the_planar_command_takes_an_element(capsys):
    from otahub.cli.main import main
    assert main(["planar", "--nx", "4", "--ny", "4", "--element", "dipole-x:0.25"]) == 0
    out = capsys.readouterr().out
    assert "short dipole along x, 0.25 lambda over ground, 7.17 dBi" in out
    assert main(["planar", "--element", "bogus"]) == 1
