"""The open-ended waveguide's directivity, integrated rather than read off the
large-aperture formula.

A TE10 aperture (E across b, cos(pi x/a) across a) in an infinite ground plane
radiates into the half space as the magnetic current 2E. Its pattern is
integrated over the hemisphere here, independently of the spec's fit; a second
route sums the same aperture as a grid of magnetic dipoles through the exact
planar power kernel. The spec's old (8/pi^2) 4 pi a b/lambda^2 is only the
large-aperture limit - it vanishes for a small aperture, whose directivity is
the magnetic dipole's 3 - and it read WR-90 1.0-3.3 dB low.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

K = 2 * math.pi
C = 2.99792458e8


def _direct(a, b, nt=300, nph=600):
    x, w = np.polynomial.legendre.leggauss(nt)
    th = 0.25 * math.pi * (x + 1)
    wt = 0.25 * math.pi * w
    ph = 2 * math.pi * (np.arange(nph) + 0.5) / nph
    T, P = np.meshgrid(th, ph, indexing="ij")
    X = K * a / 2 * np.sin(T) * np.cos(P)
    Y = K * b / 2 * np.sin(T) * np.sin(P)
    with np.errstate(divide="ignore", invalid="ignore"):
        Fx = np.where(np.abs(np.abs(X) - math.pi / 2) < 1e-9, 1 / math.pi, np.cos(X) / ((math.pi / 2) ** 2 - X ** 2))
        Fy = np.where(np.abs(Y) < 1e-12, 1.0, np.sin(Y) / np.where(np.abs(Y) < 1e-12, 1.0, Y))
    U = (np.sin(P) ** 2 + np.cos(T) ** 2 * np.cos(P) ** 2) * (Fx * Fy) ** 2
    Prad = np.sum(U * np.sin(T) * wt[:, None]) * 2 * math.pi / nph
    return 4 * math.pi * (2 / math.pi) ** 4 / Prad


def _kernel(a, b, n=20):
    from otahub.arrays import elements as E
    from otahub.arrays.planar import planar_directivity
    el = E.custom(lambda t, f: 1 - (np.sin(t) * np.cos(f)) ** 2, half_space=True, mmax=2)
    ny = max(4, int(round(n * b / a)))
    xs = (np.arange(n) + 0.5) / n * a - a / 2
    ys = (np.arange(ny) + 0.5) / ny * b - b / 2
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    pos = np.column_stack([X.ravel(), Y.ravel()])
    return planar_directivity(pos, np.cos(math.pi * pos[:, 0] / a), element=el)


def test_a_vanishing_aperture_is_a_magnetic_dipole_on_the_plane():
    assert _direct(0.01, 0.005) == pytest.approx(3.0, rel=1e-4)


def test_a_large_aperture_approaches_the_taper_efficiency():
    a, b = 6.0, 3.0
    assert _direct(a, b, 600, 1200) / (4 * math.pi * a * b) == pytest.approx(8 / math.pi ** 2, rel=0.02)


@pytest.mark.parametrize("a,b", [(0.62, 0.28), (0.85, 0.42)])
def test_two_routes_agree(a, b):
    assert 10 * math.log10(_kernel(a, b)) == pytest.approx(10 * math.log10(_direct(a, b)), abs=0.01)


@pytest.mark.parametrize("f0,A,B", [(9.0e9, 0.02286, 0.01016), (11.5e9, 0.02286, 0.01016),
                                    (5.5e9, 0.034849, 0.015799), (15e9, 0.015799, 0.007899)])
def test_the_spec_is_the_integrated_aperture_plus_the_flanges_excess(registry, f0, A, B):
    """The TE10 aperture integrated here, plus the flanged guide's excess as the FDTD
    survey fitted it (tests/data/oewg_flange_fdtd.json, its own coefficients - not the
    spec's string)."""
    d = registry["open_ended_waveguide"].synthesize(f0=f0, a_wg=A, b_wg=B)
    a, b = A * f0 / C, B * f0 / C
    D = _direct(a, b) * 10 ** (_flange_excess_db(a, b) / 10)
    assert d.metrics["directivity_dbi"] == pytest.approx(10 * math.log10(D), abs=0.01)
    assert d.metrics["effective_to_physical_area"] == pytest.approx(D / (4 * math.pi * a * b), rel=0.003)


def _flange_excess_db(a, b):
    import json
    from pathlib import Path
    fit = json.loads((Path(__file__).parent / "data" / "oewg_flange_fdtd.json").read_text())["fit"]
    return sum(c * a ** i * b ** j for c, (i, j) in zip(fit["F"], fit["F_terms"]))


def test_the_old_formula_was_one_to_three_and_a_half_decibels_low(registry):
    for f0, low in ((8.2e9, 3.49), (10e9, 2.32), (12.4e9, 1.26)):
        d = registry["open_ended_waveguide"].synthesize(f0=f0, a_wg=0.02286, b_wg=0.01016)
        old = 10 * math.log10((8 / math.pi ** 2) * 4 * math.pi * 0.02286 * 0.01016 * (f0 / C) ** 2)
        assert d.metrics["directivity_dbi"] - old == pytest.approx(low, abs=0.06)
        assert d.metrics["effective_to_physical_area"] > 1.0


def test_outside_the_single_mode_range_it_is_nan(registry):
    d = registry["open_ended_waveguide"].synthesize(f0=15e9, a_wg=0.02286, b_wg=0.01016)   # a = 1.14 lambda: TE20 propagates
    assert math.isnan(d.metrics["directivity_dbi"])
    lam = C / 10e9                                   # a = 0.9, b = 0.54 lambda: TE01 propagates, though b/a is 0.6
    assert math.isnan(registry["open_ended_waveguide"].synthesize(f0=10e9, a_wg=0.9 * lam, b_wg=0.54 * lam)
                      .metrics["directivity_dbi"])
