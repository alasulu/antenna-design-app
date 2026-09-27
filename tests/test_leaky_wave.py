"""The leaky-wave line source against its own pattern, integrated directly.

The source: an x-directed magnetic line current exp(-(alpha + j beta) x) over
0..L in a ground plane, radiating into the half space above it. With theta'
measured from the line's axis its element pattern is sin^2(theta') and the half
space is half the azimuth, so directivity is a one-dimensional integral - done
here with quad, independently of the spec's fitted closed form. The findings
pinned: directivity does not fall with scan (the spec had a cos(theta) that a
line source does not have), load loss belongs in the gain, and the exponential
aperture's beam is wider than the uniform one's.
"""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq

K = 2 * math.pi
C = 2.99792458e8


def _af2(bn, an, L, s):
    g = K * (an + 1j * bn) - 1j * K * s
    return abs((1 - np.exp(-g * L)) / g) ** 2


def _directivity(bn, an, L, ground=True):
    el = (lambda t: math.sin(t) ** 2) if ground else (lambda t: 1.0)
    U = lambda t: el(t) * _af2(bn, an, L, math.cos(t))
    tb = math.acos(bn)
    pts = [max(1e-9, tb - 0.2), tb, min(math.pi - 1e-9, tb + 0.2)]
    P = quad(lambda t: U(t) * math.sin(t), 0, math.pi, points=pts, limit=2000, epsrel=1e-10)[0]
    P *= math.pi if ground else 2 * math.pi
    ts = np.linspace(max(1e-6, tb - 0.3), min(math.pi - 1e-6, tb + 0.3), 20001)
    return 4 * math.pi * max(U(t) for t in ts) / P


def _beta(f0, a=0.02286):
    return math.sqrt(1 - (C / (2 * a * f0)) ** 2)


def _design(registry, f0, L, an):
    return registry["leaky_wave_line_source"].synthesize(f0=f0, a=0.02286, L_over_lambda=L, alpha_norm=an)


def test_the_integral_is_right_on_a_uniform_broadside_line():
    assert 10 * math.log10(_directivity(0.0, 1e-9, 10.0, ground=False)) == pytest.approx(10 * math.log10(20.0), abs=0.06)


@pytest.mark.parametrize("f0,L,an", [(8.6e9, 12.0, 0.015), (11.0e9, 25.0, 0.008), (9.5e9, 30.0, 0.004)])
def test_the_spec_directivity_is_the_integrated_one(registry, f0, L, an):
    want = 10 * math.log10(_directivity(_beta(f0), an, L))
    d = _design(registry, f0, L, an)
    assert d.metrics["directivity_dbi"] == pytest.approx(want, abs=0.08)
    eff = 1 - math.exp(-2 * an * K * L)
    assert d.metrics["gain_dbi"] == pytest.approx(want + 10 * math.log10(eff), abs=0.08)


def test_directivity_holds_as_the_beam_scans(registry):
    """Across X band the beam swings from 37 to 58 degrees; the line source's
    directivity moves by under 0.05 dB. The old cos(theta) would have cost 1.8 dB."""
    vals = [10 * math.log10(_directivity(_beta(f), 0.018, 10.0)) for f in (8.2e9, 10e9, 12.4e9)]
    assert max(vals) - min(vals) < 0.05
    spec = [_design(registry, f, 10.0, 0.018).metrics["directivity_dbi"] for f in (8.2e9, 10e9, 12.4e9)]
    assert max(spec) - min(spec) < 0.1


def test_the_ground_plane_is_worth_three_decibels():
    bn = _beta(10e9)
    diff = 10 * math.log10(_directivity(bn, 0.018, 10.0, True) / _directivity(bn, 0.018, 10.0, False))
    assert diff == pytest.approx(3.0, abs=0.1)


def _k_inf(aL):
    F = lambda x: abs((1 - np.exp(-(aL + 1j * x))) / (aL + 1j * x)) ** 2
    F0 = F(1e-12)
    return 2 * brentq(lambda x: F(x) - 0.5 * F0, 1e-9, 2 * math.pi) / (2 * math.pi)


@pytest.mark.parametrize("an,L", [(0.001, 10.0), (0.018, 10.0), (0.03, 16.0)])
def test_the_beamwidth_carries_the_leak_taper(registry, an, L):
    aL = an * K * L
    d = _design(registry, 9e9, L, an)
    bn = _beta(9e9)
    want = math.degrees(_k_inf(aL) / (L * math.cos(math.asin(bn))))
    assert d.metrics["beamwidth_deg"] == pytest.approx(want, rel=2e-3)          # the table interpolates the exact factor
    assert _k_inf(1e-6) == pytest.approx(0.8859, abs=1e-4)                      # the uniform aperture's


def test_the_old_directivity_was_five_decibels_low_at_the_default(registry):
    """2(L/lambda) cos(theta) (1 - exp(-2 alpha L)) at the default gives 10.7 dBi; the
    source integrated gives 15.7 - 3 dB for the half space, 1.8 for the cos(theta),
    0.5 for the load loss that belongs in the gain."""
    bn = _beta(10e9)
    old = 10 * math.log10(2 * 10 * math.cos(math.asin(bn)) * (1 - math.exp(-2 * 0.018 * K * 10)))
    new = _design(registry, 10e9, 10.0, 0.018).metrics["directivity_dbi"]
    assert old == pytest.approx(10.70, abs=0.01) and new - old == pytest.approx(5.0, abs=0.1)


# ------------------------------------------------------------------ the travelling-wave slot array
# Same mistake, discrete form: N broad-wall slots on the guide wall. Checked against
# the array computed exactly by the planar power kernel with the slot's element
# pattern (a magnetic dipole along the guide, half space) - otahub.arrays.

def _slot_array(N, d, th):
    from otahub.arrays import elements as E
    from otahub.arrays.planar import planar_beam_cut, planar_directivity
    from otahub.core.pattern import hpbw_deg
    slot = E.custom(lambda t, f: 1 - (np.sin(t) * np.cos(f)) ** 2, half_space=True, mmax=2)
    pos = np.column_stack([(np.arange(N) - (N - 1) / 2) * d, np.zeros(N)])
    phi = 0.0 if th >= 0 else 180.0
    D = planar_directivity(pos, np.ones(N), abs(th), phi, element=slot)
    bw = hpbw_deg(planar_beam_cut(pos, np.ones(N), "scan", abs(th), phi, n_theta=20001, element=slot), 0.0)
    return D, bw


@pytest.mark.parametrize("f0,N,dlg", [(8.8e9, 16, 0.58), (10.5e9, 24, 0.44), (12e9, 30, 0.61)])
def test_the_slot_array_is_4_n_d_over_lambda_and_its_beam_n_d_wide(registry, f0, N, dlg):
    d = registry["waveguide_slot_array_travelling_wave"].synthesize(
        f0=f0, a_wg=0.02286, b_wg=0.01016, N=N, d_over_lambdag=dlg, load_fraction=0.05)
    lam = C / f0
    th = d.metrics["beam_from_broadside_deg"]
    D, bw = _slot_array(N, d.get("spacing") / lam, th)
    assert d.metrics["directivity_dbi"] == pytest.approx(10 * math.log10(D), abs=0.09)
    assert d.metrics["gain_dbi"] == pytest.approx(10 * math.log10(D * 0.95), abs=0.09)
    assert d.metrics["beamwidth_deg"] == pytest.approx(bw, rel=3e-3)


def test_the_slot_array_directivity_is_nan_once_a_grating_lobe_is_real(registry):
    """0.75 guide wavelengths is 0.993 lambda0 - legal by the old 'below 1' rule -
    but it steers the beam 14.6 degrees off broadside, which puts a grating lobe
    in real space, and the array keeps only 0.69 of 4 N d/lambda."""
    f0 = 10e9
    d = registry["waveguide_slot_array_travelling_wave"].synthesize(
        f0=f0, a_wg=0.02286, b_wg=0.01016, N=20, d_over_lambdag=0.75, load_fraction=0.05)
    assert d.get("spacing") > d.metrics["grating_lobe_free_spacing_m"]
    assert math.isnan(d.metrics["directivity_dbi"])
    D, _ = _slot_array(20, d.get("spacing") / (C / f0), d.metrics["beam_from_broadside_deg"])
    assert D / (4 * 20 * d.get("spacing") / (C / f0)) == pytest.approx(0.687, abs=0.01)
