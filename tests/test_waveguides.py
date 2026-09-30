"""Waveguide and transmission-line checks against datasheets and limits."""
import math

import pytest

from otahub.waveguides import lines as L
from otahub.waveguides.circular import CircularWaveguide, te_root, tm_root
from otahub.waveguides.rectangular import (WR_SERIES, RectangularWaveguide,
                                           recommended_band, standard)


# ------------------------------------------------------------- rectangular

def test_wr90_te10_cutoff_matches_datasheet():
    assert standard("WR-90").dominant_cutoff_hz == pytest.approx(6.557e9, rel=1e-3)


@pytest.mark.parametrize("name,fc_ghz", [
    ("WR-284", 2.078), ("WR-137", 4.301), ("WR-90", 6.557),
    ("WR-62", 9.488), ("WR-28", 21.077), ("WR-10", 59.015),
])
def test_wr_series_cutoffs(name, fc_ghz):
    assert standard(name).dominant_cutoff_hz == pytest.approx(fc_ghz * 1e9, rel=2e-3)


def test_recommended_band_reproduces_published_wr90_range():
    """The 1.25x-1.90x heuristic should land on the published 8.2-12.4 GHz."""
    lo, hi = recommended_band("WR-90")
    assert lo == pytest.approx(8.2e9, rel=0.01)
    assert hi == pytest.approx(12.4e9, rel=0.01)


def test_standard_guide_aspect_ratios_are_near_but_not_always_exactly_two():
    """WR-90 is 0.9 x 0.4 inch, so 2.25 - the 2:1 rule of thumb is not universal.

    It matters: at 2:1 the TE20 and TE01 cutoffs coincide, while at 2.25 TE20
    arrives first and sets the top of the single-mode band on its own.
    """
    assert standard("WR-62").aspect_ratio == pytest.approx(2.0, rel=1e-6)
    assert standard("WR-28").aspect_ratio == pytest.approx(2.0, rel=1e-6)
    assert standard("WR-90").aspect_ratio == pytest.approx(2.25, rel=1e-6)
    assert standard("WR-42").aspect_ratio == pytest.approx(0.420 / 0.170, rel=1e-6)
    for name in WR_SERIES:
        assert 1.9 <= standard(name).aspect_ratio <= 2.5


def test_te01_sets_the_band_edge_only_for_squat_guides():
    wr90 = standard("WR-90")                      # aspect 2.25
    assert wr90.single_mode_band_hz[1] == pytest.approx(wr90.cutoff(2, 0))
    assert wr90.cutoff(0, 1) > wr90.cutoff(2, 0)
    wr62 = standard("WR-62")                      # aspect 2.00
    assert wr62.cutoff(0, 1) == pytest.approx(wr62.cutoff(2, 0), rel=1e-3)


def test_te20_and_te01_coincide_for_two_to_one_guides():
    g = RectangularWaveguide(0.02, 0.01)
    assert g.cutoff(2, 0) == pytest.approx(g.cutoff(0, 1), rel=1e-12)
    assert g.cutoff(2, 0) == pytest.approx(2 * g.cutoff(1, 0), rel=1e-12)


def test_the_walls_keep_their_own_permeability_when_the_filling_is_magnetic():
    """A ferrite filling (mu_r = 4) changes the wave, not the copper: the wall
    surface resistance stays sqrt(pi f mu0 / sigma). It used to take the
    filling's mu_r and doubled the loss."""
    f, sigma = 10e9, 5.8e7
    g = RectangularWaveguide(0.02286, 0.01016, mu_r=4.0, sigma=sigma)
    rs = math.sqrt(math.pi * f * 4e-7 * math.pi * 1.00000000055 / sigma)
    eta = 376.730313668 * 2.0
    fc = 2.99792458e8 / (2 * 0.02286 * 2.0)
    r = (fc / f) ** 2
    alpha = rs / (0.01016 * eta * math.sqrt(1 - r)) * (1 + 2 * 0.01016 / 0.02286 * r)
    assert g.conductor_attenuation(f) == pytest.approx(alpha, rel=1e-6)
    ferrite_walls = RectangularWaveguide(0.02286, 0.01016, mu_r=4.0, wall_mu_r=4.0)
    assert ferrite_walls.conductor_attenuation(f) == pytest.approx(2 * alpha, rel=1e-6)
    c = CircularWaveguide(0.02, mu_r=4.0)
    assert c.te01_attenuation(f) * 2 == pytest.approx(
        CircularWaveguide(0.02, mu_r=4.0, wall_mu_r=4.0).te01_attenuation(f), rel=1e-12)


def test_guide_wavelength_exceeds_free_space_and_diverges_at_cutoff():
    g = standard("WR-90")
    lam0 = 2.99792458e8 / 10e9
    assert g.guide_wavelength(10e9) > lam0
    assert g.guide_wavelength(10e9) == pytest.approx(0.039707, rel=1e-3)
    assert math.isinf(g.guide_wavelength(g.dominant_cutoff_hz))


def test_te_wave_impedance_exceeds_eta0_and_tm_falls_below():
    g = standard("WR-90")
    assert g.wave_impedance(10e9, family="TE") > 376.7
    assert g.wave_impedance(10e9, 1, 1, family="TM") < 376.7


def test_wr90_attenuation_is_about_a_tenth_of_a_db_per_metre():
    assert standard("WR-90").attenuation_db_per_m(10e9) == pytest.approx(0.11, abs=0.02)


def test_attenuation_is_infinite_below_cutoff():
    assert math.isinf(standard("WR-90").attenuation_db_per_m(5e9))


def test_perfect_conductor_has_no_conductor_loss():
    g = RectangularWaveguide(0.02286, 0.01016, sigma=float("inf"))
    assert g.conductor_attenuation(10e9) == 0.0


def test_dielectric_filling_lowers_the_cutoff():
    air = RectangularWaveguide(0.02286, 0.01016)
    filled = RectangularWaveguide(0.02286, 0.01016, eps_r=4.0)
    assert filled.dominant_cutoff_hz == pytest.approx(air.dominant_cutoff_hz / 2)


def test_mode_ordering_puts_te10_first():
    modes = standard("WR-90").modes(20e9)
    assert modes[0].label == "TE10"
    assert [m.f_cutoff_hz for m in modes] == sorted(m.f_cutoff_hz for m in modes)


def test_tm_modes_require_both_indices_nonzero():
    labels = {m.label for m in standard("WR-90").modes(30e9) if m.family == "TM"}
    assert "TM10" not in labels and "TM01" not in labels
    assert "TM11" in labels


def test_aspect_convention_is_enforced():
    with pytest.raises(ValueError, match="a >= b"):
        RectangularWaveguide(0.01, 0.02)


def test_unknown_standard_guide_is_reported():
    with pytest.raises(KeyError, match="unknown guide"):
        standard("WR-999")


# ---------------------------------------------------------------- circular

def test_bessel_roots_match_textbook():
    assert te_root(1, 1) == pytest.approx(1.8412, rel=1e-4)
    assert tm_root(0, 1) == pytest.approx(2.4049, rel=1e-4)
    assert te_root(0, 1) == pytest.approx(3.8317, rel=1e-4)


def test_te11_is_dominant_in_circular_guide():
    g = CircularWaveguide(0.01)
    assert g.modes(30e9)[0].label == "TE11"


def test_circular_single_mode_bandwidth_ratio_is_1_306():
    """TM01/TE11 = 2.4049/1.8412, markedly worse than rectangular's 2:1."""
    g = CircularWaveguide(0.01)
    lo, hi = g.single_mode_band_hz
    assert hi / lo == pytest.approx(2.4048 / 1.8412, rel=1e-3)


def test_te01_and_tm11_are_degenerate():
    """A real degeneracy: both come from the 3.8317 root, which is why TE01
    is so hard to keep pure in practice."""
    g = CircularWaveguide(0.01)
    assert g.cutoff("TE", 0, 1) == pytest.approx(g.cutoff("TM", 1, 1), rel=1e-9)


def test_te01_attenuation_falls_with_frequency():
    """The distinguishing property of TE01: loss decreases without limit."""
    g = CircularWaveguide(0.01)
    fc = g.cutoff("TE", 0, 1)
    losses = [g.te01_attenuation(m * fc) for m in (1.5, 3.0, 6.0, 12.0)]
    assert losses == sorted(losses, reverse=True)


# ------------------------------------------------------------------- lines

def test_coax_impedance_matches_the_standard_50_ohm_geometry():
    assert L.coax_impedance(1.0, 2.3023) == pytest.approx(50.0, rel=1e-3)
    assert L.coax_impedance(1.0, 3.3483, eps_r=2.1) == pytest.approx(50.0, rel=1e-3)


def test_coax_optima_are_the_textbook_values():
    o = L.coax_optimum_ratios()
    assert o["min_attenuation_z0_air"] == pytest.approx(76.7, abs=0.2)
    assert o["max_power_z0_air"] == pytest.approx(30.0, abs=0.2)
    assert o["max_power_ratio"] == pytest.approx(math.sqrt(math.e))


def test_coax_rejects_inverted_geometry():
    with pytest.raises(ValueError):
        L.coax_impedance(2.0, 1.0)


def test_microstrip_50_ohm_on_fr4():
    """Standard result: ~3.0 mm wide on 1.6 mm FR-4."""
    w = L.microstrip_width_for(50.0, 1.6e-3, 4.4)
    assert w == pytest.approx(3.0e-3, rel=0.05)


def test_microstrip_synthesis_round_trips_exactly():
    for z0 in (25.0, 50.0, 75.0, 100.0):
        w = L.microstrip_width_for(z0, 1.6e-3, 4.4)
        assert L.microstrip_impedance(w, 1.6e-3, 4.4) == pytest.approx(z0, rel=1e-8)


def test_microstrip_is_continuous_across_w_equals_h():
    """The two-branch Hammerstad forms jumped 0.4% at w/h = 1, leaving 71 ohm
    on 1.6 mm FR-4 with no width at all; Hammerstad-Jensen is continuous."""
    w = L.microstrip_width_for(71.0, 1.6e-3, 4.4)
    assert L.microstrip_impedance(w, 1.6e-3, 4.4) == pytest.approx(71.0, rel=1e-8)
    below = L.microstrip_impedance(1.6e-3 * (1 - 1e-9), 1.6e-3, 4.4)
    above = L.microstrip_impedance(1.6e-3 * (1 + 1e-9), 1.6e-3, 4.4)
    assert below == pytest.approx(above, rel=1e-8)


@pytest.mark.parametrize("er", [1.0, 2.2, 4.4, 10.2])
def test_microstrip_agrees_with_wheeler_1977(er):
    """Wheeler's independent zero-thickness formula is good to about 1%."""
    import numpy as np
    for u in np.geomspace(0.05, 20, 15):
        x = 4 / u
        a = (14 + 8 / er) / 11 * x
        wheeler = 376.730313668 / (2 * math.pi * math.sqrt(2 * (1 + er))) * math.log(
            1 + x * (a + math.sqrt(a * a + math.pi ** 2 * (1 + 1 / er) / 2)))
        assert L.microstrip_impedance(u, 1.0, er) == pytest.approx(wheeler, rel=0.012)


def test_microstrip_eps_eff_lies_between_air_and_substrate():
    ee = L.microstrip_eps_eff(3.0e-3, 1.6e-3, 4.4)
    assert 1.0 < ee < 4.4


def test_wider_strip_gives_lower_impedance():
    z = [L.microstrip_impedance(w, 1.6e-3, 4.4) for w in (0.5e-3, 1e-3, 3e-3, 10e-3)]
    assert z == sorted(z, reverse=True)


def test_stripline_matches_the_familiar_approximation():
    """Elliptic form vs 30*pi/sqrt(er) * b/(w + 0.441b) at w/b = 0.5."""
    exact = L.stripline_impedance(0.5, 1.0, 2.2)
    approx = (30 * math.pi / math.sqrt(2.2)) / (0.5 + 0.441)
    assert exact == pytest.approx(approx, rel=0.01)


def test_very_wide_and_very_narrow_lines_keep_their_digits():
    """k' = tanh(pi w / 2b) for stripline and 1 - k = 2s/(w + 2s) for CPW are
    formed directly; 1 - k^2 in floating point rounded to 1 and returned a
    zero-ohm stripline at w/b = 15."""
    wide = L.stripline_impedance(15.0, 1.0, 2.2)
    assert wide == pytest.approx(30 * math.pi / math.sqrt(2.2) / (15 + 0.441), rel=2e-3)
    assert 0 < L.stripline_impedance(40.0, 1.0, 2.2) < wide
    assert math.isfinite(L.stripline_impedance(1e-3, 1.0, 2.2))
    assert 0 < L.cpw_impedance(1e4, 1.0, 12.9) < L.cpw_impedance(1e2, 1.0, 12.9)
    assert math.isfinite(L.cpw_impedance(1e-6, 1.0, 12.9))


def test_stripline_and_cpw_ratios_run_the_right_way():
    """Both use K/K' but inverted relative to each other; limits prove which."""
    assert L.stripline_impedance(10.0, 1.0, 2.2) < L.stripline_impedance(0.1, 1.0, 2.2)
    assert L.cpw_impedance(10.0, 1.0, 12.9) < L.cpw_impedance(0.1, 1.0, 12.9)


def test_cpw_eps_eff_is_the_arithmetic_mean_with_air():
    assert L.cpw_eps_eff(12.9) == pytest.approx(6.95)


def test_quarter_wave_transformer_is_the_geometric_mean():
    assert L.quarter_wave_transformer(50, 100) == pytest.approx(70.7107, rel=1e-5)
    assert L.quarter_wave_transformer(50, 50) == pytest.approx(50.0)


def test_wr_dimensions_are_exact_inch_conversions():
    """The WR number is the broad wall in HUNDREDTHS of an inch - WR-90 is
    0.90 in, WR-2300 is 23.0 in. Storing rounded millimetres perturbs the
    aspect ratio (WR-28 came out 1.997 rather than 2.000)."""
    a, b = WR_SERIES["WR-90"]
    assert a == pytest.approx(0.900 * 0.0254, rel=1e-12)
    assert b == pytest.approx(0.400 * 0.0254, rel=1e-12)
    for name, (a, _) in WR_SERIES.items():
        # The designation is a ROUNDED nominal: WR-22 is really 0.224 in and
        # WR-62 is 0.622 in, so the label can be off by a couple of percent.
        hundredths = float(name.split("-")[1])
        assert a == pytest.approx(hundredths / 100.0 * 0.0254, rel=2e-2)
