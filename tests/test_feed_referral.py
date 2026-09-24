"""Balanis 4-79 referred to the feed, and what forgetting it did.

Both induced-EMF closed forms are referred to the current MAXIMUM of the
assumed sinusoid. For any dipole length other than lambda/2 that point is not
the feed, and the feed-point value is the current-maximum one divided by
sin^2(kL/2). The catalogue applied that to the resistance of
`dipole_arbitrary_length` and not to the reactance on the next line, so the two
described different points on the same antenna. `inductively_loaded_monopole`
copied the reactance without the referral, and for a short whip that is a
factor of ten: the coil it prescribed would not resonate the antenna, and the
efficiency it reported was four times too high.

The independent check is the method of moments by image theory - a monopole on
an infinite PEC plane has exactly half a dipole's impedance - with Hallen's
equation as a second opinion.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from otahub.num import mom
from otahub.num.hallen import hallen_dipole as _hallen_dipole

LAM = 2.99792458e8 / 3e8


def _arb(registry, L, aw=1e-4):
    return registry["dipole_arbitrary_length"].synthesize(
        f0=3e8, L_over_lambda=L, aw=aw * LAM)


# ------------------------------------------------- dipole_arbitrary_length

@pytest.mark.parametrize("L", [0.1, 0.2, 0.3, 0.4, 0.6])
def test_input_reactance_matches_the_solver_at_the_feed(L, registry):
    got = _arb(registry, L).metrics["input_reactance_ohm"]
    want = mom.input_impedance(mom.dipole(L, 1e-4, 80)).imag
    assert got == pytest.approx(want, rel=0.06)


@pytest.mark.parametrize("L", [0.1, 0.25, 0.4])
def test_resistance_and_reactance_now_refer_to_the_same_point(L, registry):
    """The inconsistency itself. Both ratios must be 1/sin^2(kL/2)."""
    m = _arb(registry, L).metrics
    r_ratio = m["input_resistance_ohm"] / m["radiation_resistance_current_max_ohm"]
    x_ratio = m["input_reactance_ohm"] / m["reactance_current_max_ohm"]
    assert r_ratio == pytest.approx(x_ratio, rel=1e-12)
    assert r_ratio == pytest.approx(1.0 / math.sin(math.pi * L) ** 2, rel=1e-12)


def test_the_unreferred_value_is_nowhere_near_the_wire(registry):
    """So that nobody 'simplifies' the referral away again."""
    m = _arb(registry, 0.1).metrics
    want = mom.input_impedance(mom.dipole(0.1, 1e-4, 80)).imag
    assert abs(m["reactance_current_max_ohm"] / want) < 0.15


def test_at_half_a_wavelength_the_referral_is_the_identity(registry):
    """Why the half-wave dipole and quarter-wave monopole were never wrong."""
    m = _arb(registry, 0.5).metrics
    assert m["input_reactance_ohm"] == pytest.approx(m["reactance_current_max_ohm"], rel=1e-12)
    assert m["input_reactance_ohm"] == pytest.approx(
        registry["half_wave_dipole"].synthesize(f0=3e8, aw=1e-4 * LAM)
        .metrics["input_reactance_ohm"], rel=1e-9)


# ----------------------------------------------- inductively_loaded_monopole

def _whip(registry, h, aw_l, **kw):
    lam = 2.99792458e8 / 1e7
    return registry["inductively_loaded_monopole"].synthesize(
        f0=1e7, h_over_lambda=h, aw=aw_l * lam, coil_q=200, r_ground=0.0, **kw)


@pytest.mark.parametrize("h", [0.03, 0.05, 0.1, 0.15])
def test_the_prescribed_coil_actually_resonates_the_whip(h, registry):
    """The builder's test: wind the coil the spec gives, then measure. The
    solver's whip plus that coil must be within a few percent of resonance."""
    d = _whip(registry, h, 3.3e-5)
    x_coil = 2 * math.pi * 1e7 * d.metrics["loading_inductance_H"]
    x_whip = 0.5 * mom.input_impedance(mom.dipole(2 * h, 3.3e-5, 80)).imag
    assert abs(x_whip + x_coil) / abs(x_whip) < 0.06


def test_hallen_agrees_on_the_whips_reactance(registry):
    d = _whip(registry, 0.05, 1e-4)
    hallen = 0.5 * _hallen_dipole(0.1, 1e-4, 160).imag
    assert d.metrics["input_reactance_unloaded_ohm"] == pytest.approx(hallen, rel=0.06)


def test_coil_loss_now_dwarfs_the_radiation_resistance(registry):
    """'Over half the radiation resistance' was the old lesson. It is nearer
    six times, which is the real reason short loaded whips are inefficient."""
    d = _whip(registry, 0.05, 3.3e-5)
    ratio = d.metrics["coil_loss_resistance_ohm"] / d.metrics["radiation_resistance_ohm"]
    assert 4.0 < ratio < 8.0, ratio


def test_efficiency_is_a_quarter_of_what_it_used_to_claim(registry):
    d = _whip(registry, 0.05, 3.3e-5, sigma=5.8e7)
    assert 0.10 < d.metrics["radiation_efficiency"] < 0.20
    assert d.metrics["efficiency_db"] < -7.0


# --------------------------------------------------- the solver's own limit

def test_the_delta_gap_limit_on_short_wires_is_real_and_recorded():
    """Why the loaded whip's RESISTANCE was left to the closed form. On a short
    wire the feed node carries a current excess over the smooth distribution
    that grows under mesh refinement; the reactance hardly moves, the
    resistance does. If this ever stops being true, the note in mom.py and the
    decision not to fit R should be revisited."""
    excess, xs = [], []
    for n in (20, 80):
        m = mom.dipole(0.1, 1e-4, n)
        s = mom.solve(m)
        z = np.array([m.node_of(k)[2] for k in range(m.n_basis)])
        sel = (np.abs(z) > 0.01) & (np.abs(z) < 0.035)
        c = np.polyfit(np.abs(z[sel]), np.abs(s.currents[sel]), 1)
        excess.append(abs(s.feed_current) / np.polyval(c, 0.0) - 1)
        xs.append(s.input_impedance.imag)
    assert excess[0] > 0.05 and excess[1] > excess[0]
    assert abs(xs[1] / xs[0] - 1) < 0.03
