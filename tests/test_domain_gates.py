"""Where a fitted or calibrated rule stops, it says NaN - the cases the project
agents' audit found answering with a number instead, each just inside and just
outside its edge."""
from __future__ import annotations

import math

import pytest

from otahub.num import loop_modal

C0 = 2.99792458e8


def _nan(x) -> bool:
    return x is None or (isinstance(x, float) and math.isnan(x))


def test_the_slot_laws_stop_where_their_checks_stop(registry):
    """slot.json copied resonant_dipole's resistance cubic without its gate: w/L = 1e-6
    read 38.8 ohm (a 915 ohm slot), and the length law ran on past w/L = 0.22."""
    hws = registry["half_wave_slot"]
    assert _nan(hws.synthesize(f0=3e8, w_over_L=1e-6).get("Rd_res"))
    assert _nan(hws.synthesize(f0=3e8, w_over_L=0.3).get("L"))
    d = hws.synthesize(f0=3e8, w_over_L=0.2)
    assert not _nan(d.get("L")) and not _nan(d.get("Rd_res"))


def test_a_thick_wire_on_a_larger_loop_has_no_driving_point_figure(registry):
    """At C = 0.3 lambda and b/a = 0.05 the delta-gap impedance follows the feed model
    (the Fourier series 8.1, 10.9, then 18.8 - j1021 ohm at 40, 60, 120 modes), and the
    spec's fit had frozen the 60-mode answer. Inside (C/lambda)^2 (b/a) <= 0.0009 the
    feed moves it under 2.5%."""
    lam = C0 / 1e9
    loop = registry["small_circular_loop"]
    for C, ba, inside in ((0.3, 0.05, False), (0.2, 0.05, False), (0.3, 0.01, True), (0.1, 0.05, True)):
        a = ba * C / (2 * math.pi)
        z40, z60 = loop_modal.input_impedance(C, a, 40), loop_modal.input_impedance(C, a, 60)
        moved = abs(z60.real / z40.real - 1)
        d = loop.synthesize(f0=1e9, C=C * lam, b=a * lam, N=1)
        r = d.get("input_resistance_driving_point_ohm")
        if inside:
            assert moved < 0.025 and not _nan(r)
        else:
            assert moved > 0.08 and _nan(r)


def test_a_single_turn_wire_is_never_the_impossible_fallback(registry):
    """The uniform-current fallback built a 36 m wire for a 0.16 m loop at eta 0.5."""
    d = registry["small_circular_loop"].synthesize(f0=3e6, C_over_lambda=0.01, eta_target=0.5, N=1)
    assert _nan(d.get("b"))


def test_the_stacked_patch_takes_its_own_edge_board(registry):
    """h = t lambda / sqrt(eps_r) at t = 0.0185 came back as 0.018499999999999996 and
    a surveyed board read NaN."""
    h = 0.0185 * C0 / 2.4e9 / math.sqrt(2.2)
    d = registry["stacked_patch"].synthesize(f0=2.4e9, eps_r=2.2, h=h, eps_r2=1.0)
    assert d.get("best_stack_bandwidth_vswr2") == pytest.approx(0.137, rel=0.05)


@pytest.mark.parametrize("M, band, beam", [(1, False, False), (2, False, True), (3, True, True),
                                           (50, True, True), (60, False, False)])
def test_the_zone_plate_fits_hold_where_the_solver_agrees(registry, M, band, beam):
    """At M = 1 the bandwidth fit read 45% high (opaque) and 65% low (phase reversal)."""
    d = registry["fresnel_zone_plate"].synthesize(f0=30e9, F=0.1, M=M, phase_levels=2)
    assert _nan(d.get("gain_bandwidth_1db")) != band
    assert _nan(d.get("hpbw_deg")) != beam


def test_a_blocker_as_wide_as_the_dish_has_no_efficiency(registry):
    """Past D the blockage integral counted blocked field again: 28 dBi at 1.5 D."""
    pf = registry["prime_focus_parabolic"]
    assert not _nan(pf.synthesize(f0=10e9, D=1.0, d_blockage=0.99).get("eta_blockage"))
    for db in (1.0, 1.5):
        d = pf.synthesize(f0=10e9, D=1.0, d_blockage=db)
        assert _nan(d.get("eta_blockage")) and _nan(d.get("gain_dbi"))


def test_the_open_guide_is_nan_at_cutoff(registry):
    lam = C0 / 10e9
    oe = registry["open_ended_waveguide"]
    assert _nan(oe.synthesize(f0=10e9, a_wg=0.5 * lam, b_wg=0.2 * lam).get("directivity_dbi"))
    assert not _nan(oe.synthesize(f0=10e9, a_wg=0.5001 * lam, b_wg=0.2 * lam).get("directivity_dbi"))


def test_the_pifa_full_width_short_keeps_to_the_surveyed_plates(registry):
    lam = C0 / 1e9
    p = registry["pifa"]
    assert not _nan(p.synthesize(f0=1e9, h=0.03 * lam, W=0.1 * lam, Ws=0.1 * lam).get("L"))
    assert _nan(p.synthesize(f0=1e9, h=0.08 * lam, W=0.1 * lam, Ws=0.1 * lam).get("L"))     # h past 0.065
    assert _nan(p.synthesize(f0=1e9, h=0.03 * lam, W=0.1 * lam, Ws=0.12 * lam).get("L"))    # strip wider than plate


@pytest.mark.parametrize("key, metric, aperture, k", [("e_plane_sectoral_horn", "hpbw_e_deg", "b1", 2),
                                                      ("h_plane_sectoral_horn", "hpbw_h_deg", "a1", 3)])
def test_the_sectoral_beamwidth_stops_at_the_fits_top(registry, key, metric, aperture, k):
    """The beamwidth fits span apertures of 1.5-40 wavelengths and ran on past 40.
    At the optimum flare the aperture is sqrt(k lambda rho): k = 2 in the E plane, 3 in H."""
    a = registry[key]
    lam = C0 / 10e9
    inside = a.synthesize(f0=10e9, a_wg=0.02286, b_wg=0.01016, rho=(30 * lam) ** 2 / (k * lam))
    outside = a.synthesize(f0=10e9, a_wg=0.02286, b_wg=0.01016, rho=(50 * lam) ** 2 / (k * lam))
    assert inside.get(aperture) / lam == pytest.approx(30, rel=0.02) and not _nan(inside.get(metric))
    assert outside.get(aperture) / lam == pytest.approx(50, rel=0.02) and _nan(outside.get(metric))
