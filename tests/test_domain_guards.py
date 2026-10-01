"""Fitted formulas return NaN outside the domain they were fitted over.

An external review found fits evaluated far outside their data with no warning,
returning impossible values: a negative horn efficiency, a negative corner-
reflector resistance, a sidelobe of -112 dB, a directivity below isotropic, a
bandwidth copied from the nearest table end. Each case here is one of those
designs, beside one inside the domain that must still give a number.
"""
from __future__ import annotations

import math

import pytest

from otahub.num import waveguide_step as ws

C = 2.99792458e8


def _nan(x) -> bool:
    return x is None or (isinstance(x, float) and math.isnan(x)) or (hasattr(x, "real") and math.isnan(abs(x)))


@pytest.mark.parametrize("key,metric,outside,inside", [
    ("conical_horn", "aperture_efficiency", dict(f0=1e10, L=3, flare=2), dict(f0=1e10, L=0.3, flare=1)),
    ("corrugated_conical_horn", "aperture_efficiency", dict(f0=1e10, L=3, flare=2), dict(f0=1e10, L=0.3, flare=1)),
    ("prime_focus_parabolic", "peak_sidelobe_db", dict(f0=1e10, D=1, f_over_D=0.4, edge_taper_db=-30),
     dict(f0=1e10, D=1, f_over_D=0.4, edge_taper_db=-11)),
    ("prime_focus_parabolic", "hpbw_deg", dict(f0=1e10, D=1, f_over_D=0.4, d_blockage=0.3),
     dict(f0=1e10, D=1, f_over_D=0.4, d_blockage=0.2)),
    ("corner_reflector_90", "input_resistance_ohm", dict(f0=3e8, S_over_lambda=0.05), dict(f0=3e8, S_over_lambda=0.5)),
    ("axial_mode_helix", "gain_dbi", dict(f0=1e9, C_over_lambda=1.2, N=20, pitch_deg=12),
     dict(f0=1e9, C_over_lambda=1.0, N=10, pitch_deg=13)),
    ("half_wave_dipole", "input_resistance_driving_point_ohm", dict(f0=100e3), dict(f0=300e6, aw=1e-4)),
    ("resonant_dipole", "input_resistance_ohm", dict(f0=100e3), dict(f0=300e6, aw=1e-4)),
    ("v_antenna_travelling", "directivity_linear", dict(f0=1e8, L_over_lambda=0.1), dict(f0=1e8, L_over_lambda=4)),
    ("rhombic", "directivity_linear", dict(f0=1e8, L_over_lambda=0.1), dict(f0=1e8, L_over_lambda=4)),
    ("conical_monopole", "bandwidth_ratio", dict(f_low=1e9, cone_half_angle_deg=5), dict(f_low=1e9, cone_half_angle_deg=45)),
    ("annular_ring_patch", "radiation_q", dict(f0=1e9, eps_r=2.2, h=0.003, ratio=1.06),
     dict(f0=1e9, eps_r=2.2, h=0.003, ratio=2.0)),
    ("cylindrical_dra", "radiation_q", dict(f0=1e10, eps_r=50, aspect=6), dict(f0=1e10, eps_r=50, aspect=2)),
])
def test_a_fit_answers_only_inside_its_domain(registry, key, metric, outside, inside):
    out = registry[key].synthesize(**outside).get(metric)
    ins = registry[key].synthesize(**inside).get(metric)
    assert _nan(out), (key, metric, out)
    assert not _nan(ins), (key, metric, ins)


def test_the_cylindrical_dra_has_no_geometry_off_its_fdtd_domain(registry):
    """k0a itself is the fit, so a flat puck gets no radius rather than a guess."""
    assert _nan(registry["cylindrical_dra"].synthesize(f0=1e10, eps_r=50, aspect=6).get("a"))


@pytest.mark.parametrize("d_in,share", [(0.85, 0.30), (0.9, 0.20), (0.8625, 0.26), (1.13, 0.05), (0.97, 0.41)])
def test_the_potter_step_launches_the_share_asked_for(registry, d_in, share):
    """Rows every 0.025 lambda read 0.011 off near the launchable limit (0.30 asked,
    0.289 launched) and 0.005 between rows; now within 0.002, checked by mode
    matching at the synthesised step, on and off the rows."""
    d = registry["conical_horn_dual_mode"].synthesize(f0=1e10, L=0.3, tm11_fraction=share, d_in_over_lambda=d_in)
    step = d.get("d_step") / (C / 1e10)
    assert ws.tm11_launch(d_in, step)[0] == pytest.approx(share, abs=0.002)


@pytest.mark.parametrize("d_in,share", [(1.2, 0.005), (0.85, 0.47), (0.85, 0.25)])
def test_the_potter_step_refuses_what_it_cannot_launch(registry, d_in, share):
    """Shares outside the table used to be clamped to its ends (0.005 asked, 0.030
    launched); below the limit the step would sit within 0.005 lambda of TM11
    cutoff."""
    d = registry["conical_horn_dual_mode"].synthesize(f0=1e10, L=0.3, tm11_fraction=share, d_in_over_lambda=d_in)
    assert _nan(d.get("d_step"))
