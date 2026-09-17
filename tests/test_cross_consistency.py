"""Cross-consistency between independently written archetypes.

The `known_cases` harness checks each archetype against numbers a human chose.
The dimensional audit checks every formula against the fact that physics has no
preferred length. This file checks archetypes against EACH OTHER, at the places
where the same physics reaches two specs by different routes.

It exists because of two finds. `corner_reflector_90` shipped for four sessions
with its null where its optimum is, passing its own cited cases the whole time.
`rectangular_patch_inset` agreed with Balanis to 1.7% while both of its
conductance terms were wrong, in cancelling directions. Neither a citation nor a
unit check would catch those; a second, independent route to the same number
does.
"""
from __future__ import annotations

import math

import pytest

ETA0 = 376.730313412
F = 300e6


@pytest.fixture(scope="module")
def syn(registry):
    def _syn(key, **kw):
        return registry[key].synthesize(**kw)
    return _syn


# ------------------------------------------------------------- the wire family
# Five archetypes are, in one limit or another, the same half-wave dipole.

def test_arbitrary_length_dipole_reduces_to_the_half_wave_one(syn):
    hw = syn("half_wave_dipole", f0=F, aw=1e-4)
    arb = syn("dipole_arbitrary_length", f0=F, L_over_lambda=0.5, aw=1e-4)
    assert arb.metrics["radiation_resistance_current_max_ohm"] == pytest.approx(
        hw.metrics["radiation_resistance_ohm"], rel=1e-9)
    assert arb.metrics["directivity_linear"] == pytest.approx(
        hw.metrics["directivity_linear"], rel=1e-3)


def test_arbitrary_length_dipole_reduces_to_the_short_one(syn):
    """As L -> 0 the exact Si/Ci solution must land on 1.5 and on the
    20*pi^2*(L/lambda)^2 triangular-current resistance."""
    short = syn("short_dipole", f0=F, L_over_lambda=0.01)
    arb = syn("dipole_arbitrary_length", f0=F, L_over_lambda=0.01, aw=1e-4)
    assert arb.metrics["directivity_linear"] == pytest.approx(
        short.metrics["directivity_linear"], rel=5e-3)
    assert arb.metrics["input_resistance_ohm"] == pytest.approx(
        short.metrics["radiation_resistance_ohm"], rel=2e-2)


def test_monopole_is_exactly_half_a_dipole(syn):
    hw = syn("half_wave_dipole", f0=F, aw=1e-4)
    mono = syn("quarter_wave_monopole", f0=F, aw=1e-4)
    assert mono.metrics["radiation_resistance_ohm"] == pytest.approx(
        hw.metrics["radiation_resistance_ohm"] / 2, rel=1e-9)
    assert mono.metrics["directivity_linear"] == pytest.approx(
        hw.metrics["directivity_linear"] * 2, rel=1e-9)


def test_folded_dipole_is_four_times_a_plain_one(syn):
    hw = syn("half_wave_dipole", f0=F, aw=1e-4)
    fold = syn("folded_dipole", f0=F, N=2, aw=1e-4)
    assert fold.metrics["input_resistance_ohm"] == pytest.approx(
        4 * hw.metrics["radiation_resistance_ohm"], rel=5e-3)


def test_turnstile_on_axis_matches_a_single_dipole(syn):
    """Splitting the power two ways is exactly offset by the two orthogonal
    polarisations adding on axis. If this drifts, one of the two has changed."""
    hw = syn("half_wave_dipole", f0=F, aw=1e-4)
    turn = syn("turnstile_dipole", f0=F, L_over_lambda=0.5, aw=1e-4)
    assert turn.metrics["axial_directivity_linear"] == pytest.approx(
        hw.metrics["directivity_linear"], rel=1e-3)


def test_unloaded_top_loaded_monopole_is_the_textbook_short_monopole(syn):
    d = syn("top_loaded_monopole", f0=10e6, h_over_lambda=0.05, beta_top=0.0)
    assert d.metrics["radiation_resistance_ohm"] == pytest.approx(
        40 * math.pi ** 2 * 0.05 ** 2, rel=1e-6)


# ------------------------------------------------------------------- Babinet
# A slot and its complementary dipole are two specs written from one principle.

def test_slot_times_its_complementary_dipole_is_eta_squared_over_four(syn):
    hw = syn("half_wave_dipole", f0=F, aw=1e-4)
    slot = syn("half_wave_slot", f0=F)
    z_dipole = abs(complex(hw.metrics["radiation_resistance_ohm"],
                           hw.metrics["input_reactance_ohm"]))
    product = slot.metrics["input_impedance_magnitude_ohm"] * z_dipole
    assert product == pytest.approx(ETA0 ** 2 / 4, rel=2e-3)


def test_folded_slot_reduces_to_the_plain_slot(syn):
    slot = syn("half_wave_slot", f0=F)
    fs = syn("folded_slot", f0=F, N=1)
    assert fs.metrics["input_resistance_ohm"] == pytest.approx(
        slot.metrics["resonant_resistance_ohm"], rel=1e-9)


def test_folding_a_slot_divides_where_folding_a_dipole_multiplies(syn):
    """The dual relationship, checked against the folded DIPOLE rather than
    asserted: one steps up by N^2, the other down by the same factor."""
    hw = syn("half_wave_dipole", f0=F, aw=1e-4)
    fold = syn("folded_dipole", f0=F, N=3, aw=1e-4)
    fs1 = syn("folded_slot", f0=F, N=1)
    fs3 = syn("folded_slot", f0=F, N=3)
    dipole_step = fold.metrics["input_resistance_ohm"] / \
        hw.metrics["radiation_resistance_ohm"]
    slot_step = fs1.metrics["input_resistance_ohm"] / \
        fs3.metrics["input_resistance_ohm"]
    assert dipole_step == pytest.approx(9.0, rel=5e-3)
    assert slot_step == pytest.approx(dipole_step, rel=5e-3)


# ----------------------------------------------------------------------- cones

def test_conical_monopole_is_half_a_biconical(syn):
    bic = syn("biconical", f0=1e9, theta_h=math.radians(10.0))
    con = syn("conical_monopole", f_low=1e9, cone_half_angle_deg=5.0)
    assert con.metrics["characteristic_impedance_ohm"] == pytest.approx(
        bic.metrics["characteristic_impedance_ohm"] / 2, rel=1e-9)


def test_discone_quotes_the_same_biconical_impedance_as_the_biconical_spec(syn):
    dis = syn("discone", f_low=100e6, cone_half_angle_deg=30.0)
    bic = syn("biconical", f0=1e9, theta_h=math.radians(60.0))
    assert dis.metrics["biconical_impedance_ohm"] == pytest.approx(
        bic.metrics["characteristic_impedance_ohm"], rel=1e-9)


# --------------------------------------- archetypes against the library modules

def test_waveguide_archetypes_agree_with_the_waveguide_module(syn):
    from otahub.waveguides.rectangular import standard

    wr90 = standard("WR-90")
    guide = dict(a_wg=0.02286, b_wg=0.01016)
    oewg = syn("open_ended_waveguide", f0=10e9, **guide)
    assert oewg.metrics["cutoff_te10_hz"] == pytest.approx(
        wr90.dominant_cutoff_hz, rel=1e-9)
    assert oewg.metrics["guide_wavelength_m"] == pytest.approx(
        wr90.guide_wavelength(10e9), rel=1e-9)
    leaky = syn("leaky_wave_line_source", f0=10e9, a=0.02286)
    assert leaky.metrics["guide_wavelength_m"] == pytest.approx(
        wr90.guide_wavelength(10e9), rel=1e-9)


def test_slot_array_reuses_the_single_slot_conductance(syn):
    """Written months apart from the same Stevenson result; a drift in either
    would be invisible inside its own known cases."""
    guide = dict(a_wg=0.02286, b_wg=0.01016)
    array = syn("waveguide_slot_array_resonant", f0=10e9, N=12, **guide)
    single = syn("waveguide_longitudinal_slot", f0=10e9, x1=0.003, **guide)
    assert array.metrics["g1_normalised"] == pytest.approx(
        single.metrics["g1_normalised"], rel=1e-12)


def test_travelling_wave_slot_array_degenerates_to_the_resonant_one(syn):
    guide = dict(a_wg=0.02286, b_wg=0.01016)
    tw = syn("waveguide_slot_array_travelling_wave", f0=10e9,
             d_over_lambdag=0.5, **guide)
    res = syn("waveguide_slot_array_resonant", f0=10e9, N=12, **guide)
    assert tw.metrics["beam_from_broadside_deg"] == pytest.approx(0.0, abs=1e-9)
    assert tw.metrics["broadside_spacing_m"] == pytest.approx(
        res.parameters["spacing"], rel=1e-12)


def test_potter_horn_cutoffs_are_the_circular_guide_bessel_zeros(syn):
    pot = syn("conical_horn_dual_mode", f0=10e9, L=0.3)
    k0 = 2 * math.pi * 10e9 / 2.99792458e8
    assert pot.metrics["te11_cutoff_diameter_m"] == pytest.approx(
        2 * 1.8412 / k0, rel=1e-3)
    assert pot.metrics["tm11_cutoff_diameter_m"] == pytest.approx(
        2 * 3.8317 / k0, rel=1e-3)


def test_horns_of_equal_length_share_an_aperture(syn):
    smooth = syn("conical_horn", f0=10e9, L=0.3)
    corrugated = syn("corrugated_conical_horn", f0=10e9, L=0.3)
    potter = syn("conical_horn_dual_mode", f0=10e9, L=0.3)
    assert corrugated.parameters["dm"] == pytest.approx(
        smooth.parameters["dm"], rel=1e-12)
    assert potter.parameters["dm"] == pytest.approx(
        smooth.parameters["dm"], rel=1e-12)


def test_slot_array_factor_matches_the_arrays_module(syn):
    array = syn("waveguide_slot_array_resonant", f0=10e9, N=12,
                a_wg=0.02286, b_wg=0.01016)
    d_over_lambda = array.parameters["spacing"] * 10e9 / 2.99792458e8
    expected = 2 * 12 * d_over_lambda
    assert 10 ** (array.metrics["array_factor_directivity_dbi"] / 10) == \
        pytest.approx(expected, rel=1e-9)


# ------------------------------------------------------------ corner reflectors

def test_narrower_corner_gives_more_gain(syn):
    """Physically required, and the check that would have caught the shipped
    90-degree array factor having its null where its optimum is."""
    c90 = syn("corner_reflector_90", f0=F, S_over_lambda=0.5)
    c60 = syn("corner_reflector_60", f0=F, S_over_lambda=0.5)
    assert c60.metrics["directivity_dbi"] > c90.metrics["directivity_dbi"] + 1.0


def test_corner_reflectors_peak_where_they_should_and_null_at_a_wavelength(syn):
    for key in ("corner_reflector_90", "corner_reflector_60"):
        null = syn(key, f0=F, S_over_lambda=0.9999)
        assert abs(null.metrics["array_factor_field"]) < 0.02, (
            f"{key} should null at S = lambda")
    half = syn("corner_reflector_90", f0=F, S_over_lambda=0.5)
    assert half.metrics["array_factor_field"] == pytest.approx(4.0, rel=1e-6), \
        "the 90-degree corner peaks at S = lambda/2, it does not null there"


def test_corner_reflector_self_resistance_is_the_isolated_dipole(syn):
    hw = syn("half_wave_dipole", f0=F, aw=1e-4)
    for key in ("corner_reflector_90", "corner_reflector_60"):
        d = syn(key, f0=F, S_over_lambda=0.5)
        assert d.metrics["self_resistance_ohm"] == pytest.approx(
            hw.metrics["radiation_resistance_ohm"], rel=1e-9)


# ----------------------------------------------------------------------- DRAs

def test_three_dra_shapes_agree_on_size_and_directivity(syn):
    """Three unrelated models - exact Mie for the hemisphere, published curve
    fits for the cylinder, magnetic-wall algebra for the brick. At the same
    permittivity and frequency they must land on comparable volumes."""
    volumes = {
        "hemispherical": syn("hemispherical_dra", f0=10e9, eps_r=10.0),
        "cylindrical": syn("cylindrical_dra", f0=10e9, eps_r=10.0, aspect=1.0),
        "rectangular": syn("rectangular_dra", f0=10e9, eps_r=10.0,
                           aspect_wd=2.0, aspect_Ld=2.0),
    }
    vols = {k: d.metrics["volume_mm3"] for k, d in volumes.items()}
    assert max(vols.values()) / min(vols.values()) < 2.5, vols
    # all three radiate as a horizontal magnetic dipole over a ground plane
    dbi = [d.metrics["directivity_dbi"] for d in volumes.values()]
    assert max(dbi) - min(dbi) < 1e-9


def test_dra_q_rises_steeply_with_permittivity_in_every_shape(syn):
    """The eps_r^1.3 law, reached independently by Mie theory and by published
    fits. A shape whose Q did not rise would be a broken model."""
    for key, extra, metric in (
            ("hemispherical_dra", {}, "radiation_q"),
            ("cylindrical_dra", {"aspect": 1.0}, "radiation_q"),
            ("rectangular_dra", {"aspect_wd": 2.0, "aspect_Ld": 2.0},
             "radiation_q_indicative")):
        low = syn(key, f0=10e9, eps_r=10.0, **extra).metrics[metric]
        high = syn(key, f0=10e9, eps_r=20.0, **extra).metrics[metric]
        assert high / low == pytest.approx(2 ** 1.3, rel=0.15), (
            f"{key}: Q went {low:.3f} -> {high:.3f} when eps_r doubled")


# ------------------------------------------------------------- inset patch

def test_inset_patch_conductances_match_direct_quadrature(syn):
    """Both terms checked separately against the integrals they approximate.
    The spec once agreed with Balanis to 1.7% with BOTH terms wrong, in
    cancelling directions, so the sum alone proves nothing."""
    from scipy.integrate import quad
    from scipy.special import jv, sici

    d = syn("rectangular_patch_inset", f0=10e9, eps_r=2.2, h=0.001588,
            Z_target=50.0)
    k0 = 2 * math.pi * 10e9 / 2.99792458e8
    W, L = d.get("W"), d.get("L")

    X = k0 * W
    i1 = -2 + math.cos(X) + X * sici(X)[0] + math.sin(X) / X

    def integrand(th, bessel):
        c = math.cos(th)
        inner = X / 2 if abs(c) < 1e-12 else math.sin(X / 2 * c) / c
        weight = jv(0, k0 * L * math.sin(th)) if bessel else 1.0
        return inner ** 2 * math.sin(th) ** 3 * weight

    i12 = quad(integrand, 0, math.pi, args=(True,), limit=400)[0]
    g1_exact = i1 / (120 * math.pi ** 2)
    g12_exact = i12 / (120 * math.pi ** 2)

    assert d.metrics["self_conductance_S"] == pytest.approx(g1_exact, rel=1e-9)
    assert d.metrics["mutual_conductance_S"] == pytest.approx(g12_exact, rel=5e-3)
    assert d.metrics["edge_resistance_ohm"] == pytest.approx(
        1 / (2 * (g1_exact + g12_exact)), rel=5e-3)
