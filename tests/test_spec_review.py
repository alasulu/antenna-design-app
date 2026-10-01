"""Spec faults an external review found, each pinned to an independent number.

The reference values here were computed outside the specs - by the review's own
scripts or from the closed forms quoted beside them - so a spec that drifts back
to the old formula fails.
"""
from __future__ import annotations

import math

import pytest

C = 2.99792458e8
ETA0 = 376.730313412


def test_the_half_wave_slot_reports_the_slot_it_designs(registry):
    """It synthesised the resonant length but reported the impedance of an
    exactly-half-wave slot (363 - j211 ohm); its own resonant law says 468."""
    d = registry["half_wave_slot"].synthesize(f0=3e8)
    z = complex(d.metrics["input_impedance_ohm"])
    assert z.real == pytest.approx(d.metrics["resonant_resistance_ohm"], rel=1e-12) and z.imag == 0
    assert d.terminal_impedance()[0] == z
    assert complex(d.metrics["complement_impedance_ohm"]) == pytest.approx(362.8 - 210.9j, abs=0.2)


def test_the_cavity_is_a_quarter_guide_wavelength_deep(registry):
    """Its TE10 mode, not free space: 135.586 mm at 1 GHz for a 0.6 lambda cavity
    (74.95 mm was lambda0/4); at 0.5 lambda the mode is at cutoff."""
    d = registry["cavity_backed_slot"].synthesize(f0=1e9)
    assert d.get("cavity_depth") == pytest.approx(0.135586, rel=1e-5)
    assert d.metrics["input_resistance_ohm"] == pytest.approx(ETA0 ** 2 / (2 * d.get("Rd_res")), rel=1e-6)
    assert d.metrics["input_resistance_ohm"] > 900                      # not 2 x 363
    assert math.isnan(registry["cavity_backed_slot"].synthesize(f0=1e9, cavity_width_over_lambda=0.5).get("cavity_depth"))


def test_the_lpda_follows_carrel_and_the_spacing_it_is_given(registry):
    d = registry["lpda"].synthesize(f_low=1e8, f_high=1e9, tau=0.9)
    assert d.get("N_elements") == 28
    assert d.get("B_ar") == pytest.approx(1.1 + 7.7 * 0.01 / math.tan(math.radians(d.get("alpha_deg"))), rel=1e-12)
    loose = registry["lpda"].synthesize(f_low=1e8, f_high=1e9, tau=0.9, sigma=0.1)
    assert loose.get("alpha_deg") == pytest.approx(14.03624, rel=1e-6)
    assert math.isnan(loose.metrics["directivity_dbi"])                # Carrel's figure is for the optimum


def test_the_long_wire_terminates_in_its_characteristic_impedance(registry):
    """It was set to the radiation resistance, 210 ohm, VSWR 2.4 on a 500 ohm line."""
    d = registry["long_wire_travelling"].synthesize(f0=3e7, height=10.0, aw=0.001)
    assert d.metrics["termination_resistance_ohm"] == pytest.approx(60 * math.log(2 * 10.0 / 0.001), rel=1e-12)
    assert "termination_resistance_ohm" not in registry["long_wire_travelling"].synthesize(f0=3e7).metrics


def test_the_chu_sphere_encloses_the_hat(registry):
    """Radius h left the hat outside: 35.43 where sqrt(h^2 + a_hat^2) gives 17.84."""
    lam = C / 1e7
    d = registry["top_loaded_monopole"].synthesize(f0=1e7, h_over_lambda=0.05, hat_radius_over_lambda=0.04)
    assert d.get("a_hat") / lam == pytest.approx(0.04)
    assert d.metrics["chu_limit_q"] == pytest.approx(17.84181, rel=1e-5)


def test_the_circular_patch_probe_stays_on_the_copper(registry):
    """rho_frac is on the fringing-corrected radius: at 0.9 the probe stood 0.31 mm
    off a 10 GHz Balanis patch."""
    off = registry["circular_patch"].synthesize(f0=1e10, eps_r=2.2, h=0.001588, rho_frac=0.9)
    assert math.isnan(off.get("probe_radius_m"))
    on = registry["circular_patch"].synthesize(f0=1e10, eps_r=2.2, h=0.001588, rho_frac=0.7)
    assert on.get("probe_radius_m") < on.get("a")


def test_the_triangle_reports_the_area_it_builds(registry):
    d = registry["triangular_patch"].synthesize(f0=1e9, eps_r=2.2, h=0.0192014)
    assert d.metrics["area_m2"] == pytest.approx(math.sqrt(3) / 4 * d.get("a_side") ** 2, rel=1e-12)
    assert d.metrics["area_m2"] == pytest.approx(0.005498495, rel=1e-5)     # not the cavity side's 0.006424


def test_a_four_level_zone_plates_finest_feature_is_a_quarter_wave_step(registry):
    """r(4) - r(3.5) = 5.45755 mm at 30 GHz, F 0.15 m; the whole zone is 11.26."""
    d = registry["fresnel_zone_plate"].synthesize(f0=30e9, F=0.15, M=4, phase_levels=4)
    assert d.metrics["outer_zone_width_m"] == pytest.approx(0.00545755, rel=1e-5)
    amp = registry["fresnel_zone_plate"].synthesize(f0=30e9, F=0.15, M=4, phase_levels=1)
    assert amp.metrics["outer_zone_width_m"] == pytest.approx(0.0112569, rel=1e-4)


def test_a_cavity_backed_spiral_gains_directivity_not_gain(registry):
    d = registry["archimedean_spiral"].synthesize(**registry["archimedean_spiral"].spec.known_cases[0].given)
    free = d.metrics["directivity_dbi"]
    assert d.metrics["cavity_backed_directivity_dbi"] == pytest.approx(free + 10 * math.log10(2))
    assert d.metrics["cavity_backed_gain_dbi"] == pytest.approx(free)


def test_a_slot_array_near_cutoff_counts_its_grating_lobes(registry):
    """At 7 GHz in WR-90 the half-guide-wavelength spacing is 1.43 lambda0; 2Nd/lambda
    read 15.35 dBi where the twelve elements' array factor gives 10.59."""
    d = registry["waveguide_slot_array_resonant"].synthesize(f0=7e9, a_wg=0.02286, b_wg=0.01016, N=12)
    assert d.metrics["spacing_over_lambda0"] == pytest.approx(1.428404, rel=1e-5)
    assert d.metrics["array_factor_directivity_dbi"] == pytest.approx(10.5895, abs=1e-3)


def test_the_diagonal_horn_claims_no_gain_over_a_square_pyramidal_one(registry):
    """Their aperture efficiencies are equal at every phase error."""
    from otahub.num import horn_pattern as hp
    d = registry["diagonal_horn"].synthesize(f0=1e10, a_ap=0.06, R_axial=0.2)
    assert "gain_advantage_over_pyramidal_db" not in d.metrics
    assert d.metrics["aperture_efficiency_with_phase_error"] == pytest.approx(
        hp.diagonal_efficiency(d.metrics["max_phase_error_wavelengths"]), rel=1e-5)


@pytest.mark.parametrize("b_over_lambda", [1e-4, 1e-3])
def test_a_loop_built_to_its_circumference_presents_the_impedance_reported(registry, b_over_lambda):
    """C was the nominal 1.09 lambda while the impedance was the resonant loop's:
    built to C it presented 162.5 + j130.5 ohm against a reported 138.3. C is the
    resonant circumference now; the independent modal solver agrees."""
    from otahub.num import loop_modal
    lam = C / 3e8
    d = registry["one_wavelength_circular_loop"].synthesize(f0=3e8, b=b_over_lambda * lam)
    z = loop_modal.input_impedance(d.get("C") / lam, b_over_lambda)
    assert z.real == pytest.approx(complex(d.metrics["input_impedance_ohm"]).real, rel=0.005)
    assert abs(z.imag) < 0.01 * z.real
    assert d.get("C_textbook") / lam == pytest.approx(1.09)


def test_the_halo_closes_its_own_ring(registry):
    """pi Dm was 1.0267 m where Lc + g was 0.9878 at 146 MHz."""
    d = registry["halo_loop"].synthesize(f0=146e6, b=0.0041067)
    assert math.pi * d.get("Dm") == pytest.approx(d.get("Lc") + d.get("g"), rel=1e-12)
