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
    """Four times, but four times WHAT.

    This used to demand exactly 4 x 73.079, and the folded dipole's resistance
    is now solved rather than defined as that product - so it no longer lands
    there, and should not. 73.079 is the induced-EMF figure for an assumed
    sinusoidal current; a real dipole's RESONANT resistance is about 71.9, and
    the folded dipole is four times that. Against the induced-EMF constant the
    ratio is therefore about 3.91 rather than 4.00. The tight 4:1 between
    resonant resistances is checked directly in test_folded_dipole.py.
    """
    hw = syn("half_wave_dipole", f0=F, aw=1e-4)
    fold = syn("folded_dipole", f0=F, N=2, aw=1e-4)
    ratio = (fold.metrics["input_resistance_ohm"]
             / hw.metrics["radiation_resistance_ohm"])
    assert ratio == pytest.approx(4.0, rel=0.03)
    assert ratio < 4.0, "the induced-EMF 73.079 overstates a resonant dipole"


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

def test_the_slot_borrows_the_dipole_constants_the_dipole_spec_computes(syn):
    """What this can check, and what it cannot.

    It used to be called a Babinet test - "slot times its complementary dipole
    is eta squared over four" - but the slot spec DEFINES its impedance as
    eta0^2 / (4 Z_dipole) from dipole constants it carries itself, so that
    product is eta0^2/4 by construction and asserting it checked nothing about
    Babinet. Babinet cannot be checked here at all: the solver models wires,
    not a slot in a conducting screen.

    What CAN drift is the borrowed constants. The slot carries the dipole's
    impedance as fixed numbers; the dipole spec computes it. They must agree.
    """
    hw = syn("half_wave_dipole", f0=F, aw=1e-4)
    slot = syn("half_wave_slot", f0=F)
    assert slot.parameters["Rd"] == pytest.approx(hw.metrics["radiation_resistance_ohm"], rel=2e-3)
    assert slot.parameters["Xd"] == pytest.approx(hw.metrics["input_reactance_ohm"], rel=2e-3)


def test_folded_slot_reduces_to_the_plain_slot(syn):
    slot = syn("half_wave_slot", f0=F)
    fs = syn("folded_slot", f0=F, N=1)
    assert fs.metrics["input_resistance_ohm"] == pytest.approx(
        slot.metrics["resonant_resistance_ohm"], rel=1e-9)


def test_folding_a_slot_divides_where_folding_a_dipole_multiplies(syn):
    """The dual relationship, checked against the folded DIPOLE rather than
    asserted: one steps up by N^2, the other down by the same factor.

    Compared within each archetype rather than between them. The folded
    dipole's absolute resistance is now solved, the folded slot's is still
    asserted, and dividing one by the other would test the difference between
    a derived number and a stated one instead of the duality.
    """
    fd1 = syn("folded_dipole", f0=F, N=1, aw=1e-4)
    fd3 = syn("folded_dipole", f0=F, N=3, aw=1e-4)
    fs1 = syn("folded_slot", f0=F, N=1)
    fs3 = syn("folded_slot", f0=F, N=3)
    # Each archetype's own N-scaling, so a derived absolute value on one side
    # is not being compared against an asserted one on the other.
    dipole_step = (fd3.metrics["input_resistance_ohm"]
                   / fd1.metrics["input_resistance_ohm"])
    slot_step = (fs1.metrics["input_resistance_ohm"]
                 / fs3.metrics["input_resistance_ohm"])
    assert dipole_step == pytest.approx(9.0, rel=5e-3)
    assert slot_step == pytest.approx(dipole_step, rel=5e-3)


# ----------------------------------------------------------------------- cones

@pytest.mark.parametrize("half", [5.0, 30.0, 47.0])
def test_conical_monopole_is_half_a_biconical(syn, half):
    """Same physical cone, same HALF angle.

    This used to compare a biconical at theta_h = 10 degrees with a conical
    monopole at 5 - doubled on one side, because the biconical spec was using
    the full-angle form cot(theta/4) against a half-angle input. The test had
    been written to make the two specs agree rather than to catch that they
    disagreed. Both now use the half angle, and the impedance halves by the image.
    """
    bic = syn("biconical", f0=1e9, theta_h=math.radians(half))
    con = syn("conical_monopole", f_low=1e9, cone_half_angle_deg=half)
    assert con.metrics["characteristic_impedance_ohm"] == pytest.approx(
        bic.metrics["characteristic_impedance_ohm"] / 2, rel=1e-9)


def test_at_47_degrees_the_image_relation_holds_at_the_cutoff_too(syn):
    """The two specs quote directivity at their own cutoffs, and those differ:
    the monopole's is where VSWR in 50 ohm reaches 2, the bicone's where VSWR
    against its own Zc does. At 47 degrees the monopole's Zc IS 50 ohm, so the
    two cutoffs are the same frequency, and image theory must then hold
    exactly: the same slant, twice the directivity."""
    bic = syn("biconical", f0=1e9, theta_h=math.radians(47.0))
    con = syn("conical_monopole", f_low=1e9, cone_half_angle_deg=47.0)
    assert con.get("slant") == pytest.approx(bic.get("Lc"), rel=2e-4)
    assert con.metrics["directivity_linear"] == pytest.approx(2 * bic.metrics["directivity_linear"], rel=5e-4)


def test_discone_quotes_the_same_biconical_impedance_as_the_biconical_spec(syn):
    """Same half angle on both sides now; it used to pair a 30 degree discone
    with a 60 degree biconical to cover the biconical spec's angle error."""
    dis = syn("discone", f_low=100e6, cone_half_angle_deg=30.0)
    bic = syn("biconical", f0=1e9, theta_h=math.radians(30.0))
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
    """Three unrelated models - the exact Mie pole for the hemisphere, published
    curve fits for the cylinder, FDTD ringdowns for the brick. At the same
    permittivity and frequency they must land on comparable volumes.

    The tolerance is deliberately tight. While the brick used an
    all-magnetic-wall resonance it was 1.41x the other two and this test had to
    allow 2.5x to pass, which made it nearly useless as a guard. With the
    transcendental solved the spread is 1.06-1.14x across eps_r from 8 to 40,
    so 1.25x is now a real constraint.
    """
    volumes = {
        "hemispherical": syn("hemispherical_dra", f0=10e9, eps_r=10.0),
        "cylindrical": syn("cylindrical_dra", f0=10e9, eps_r=10.0, aspect=1.0),
        "rectangular": syn("rectangular_dra", f0=10e9, eps_r=10.0,
                           aspect_wd=2.0, aspect_Ld=2.0),
    }
    vols = {k: d.metrics["volume_mm3"] for k, d in volumes.items()}
    assert max(vols.values()) / min(vols.values()) < 1.25, vols
    # all three radiate as a horizontal magnetic dipole over a ground plane
    dbi = [d.metrics["directivity_dbi"] for d in volumes.values()]
    assert max(dbi) - min(dbi) < 1e-9


def test_dra_q_rises_steeply_with_permittivity_in_every_shape(syn):
    """Q rises a little faster than eps_r^1.1 between 10 and 20 in every shape:
    the hemisphere's exact pole and the brick's and the cylinder's ringdowns
    all go up 2.23-2.25 times - not the 2^1.3 the published cylinder fit
    assumed. A shape whose Q did not rise would be a broken model."""
    for key, extra, metric in (
            ("hemispherical_dra", {}, "radiation_q"),
            ("cylindrical_dra", {"aspect": 1.0}, "radiation_q"),
            ("rectangular_dra", {"aspect_wd": 2.0, "aspect_Ld": 2.0}, "radiation_q")):
        low = syn(key, f0=10e9, eps_r=10.0, **extra).metrics[metric]
        high = syn(key, f0=10e9, eps_r=20.0, **extra).metrics[metric]
        # three unrelated routes (Mie pole, brick ringdowns, puck ringdowns)
        # agree on the ratio to 1%; the tolerance leaves room for the fits
        assert high / low == pytest.approx(2.24, rel=0.02), (
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


def test_dra_shapes_stay_close_across_the_permittivity_range(syn):
    """Not one lucky point: the three models must track each other as eps_r
    moves, which is what would expose a wrong exponent in any of them."""
    for eps_r in (8.0, 10.0, 20.0, 40.0):
        vols = [
            syn("hemispherical_dra", f0=10e9, eps_r=eps_r).metrics["volume_mm3"],
            syn("cylindrical_dra", f0=10e9, eps_r=eps_r,
                aspect=1.0).metrics["volume_mm3"],
            syn("rectangular_dra", f0=10e9, eps_r=eps_r, aspect_wd=2.0,
                aspect_Ld=2.0).metrics["volume_mm3"],
        ]
        assert max(vols) / min(vols) < 1.25, f"eps_r={eps_r}: {vols}"


def test_rectangular_dra_beats_the_magnetic_wall_model_it_replaced(syn):
    """The oversize factor is reported so the improvement is visible, not just
    claimed. It must be real and it must be bounded."""
    for eps_r in (8.0, 20.0):
        d = syn("rectangular_dra", f0=10e9, eps_r=eps_r, aspect_wd=2.0,
                aspect_Ld=2.0)
        assert 1.02 < d.metrics["magnetic_wall_oversize"] < 1.30
        assert d.metrics["height_magnetic_wall_m"] > d.get("d")


# ---------------------------------------------------------------------- horns

def test_horn_gain_agrees_with_the_published_directivity_expression(syn):
    """Each sectoral horn carries its efficiency as a closed form AND Balanis's
    directivity expression as published. They are algebraically the same thing,
    so any drift between them is a transcription error."""
    for key, guide, flares in (
            ("e_plane_sectoral_horn", {"a_wg": 0.02286}, (0.6, 0.8, 1.0, 1.2, 1.4)),
            ("h_plane_sectoral_horn", {"b_wg": 0.01016}, (0.6, 1.0, 1.4))):
        for flare in flares:
            d = syn(key, f0=10e9, rho=0.3, flare=flare, **guide)
            from_aperture = 10 ** (d.metrics["gain_dbi"] / 10)
            assert d.metrics["directivity_balanis_form"] == pytest.approx(
                from_aperture, rel=1e-12), f"{key} at flare {flare}"


def test_the_optimum_flare_is_actually_a_maximum(syn):
    """The point of an optimum-gain horn. A spec with a PINNED aperture
    efficiency cannot express this at all - its gain rises without limit as the
    flare grows, which would recommend over-flaring indefinitely."""
    for key, guide in (("e_plane_sectoral_horn", {"a_wg": 0.02286}),
                       ("h_plane_sectoral_horn", {"b_wg": 0.01016})):
        gains = {f: syn(key, f0=10e9, rho=0.3, flare=f, **guide).metrics["gain_dbi"]
                 for f in (0.7, 0.85, 1.0, 1.15, 1.3)}
        assert gains[1.0] == max(gains.values()), f"{key}: {gains}"
        assert gains[1.3] < gains[1.0] - 0.5, "over-flaring must cost real gain"


def test_optimum_horn_phase_errors_are_the_canonical_quarter_and_three_eighths(syn):
    e = syn("e_plane_sectoral_horn", f0=10e9, a_wg=0.02286, rho=0.3, flare=1.0)
    h = syn("h_plane_sectoral_horn", f0=10e9, b_wg=0.01016, rho=0.3, flare=1.0)
    assert e.metrics["max_phase_error_wavelengths"] == pytest.approx(0.25, rel=1e-9)
    assert h.metrics["max_phase_error_wavelengths"] == pytest.approx(0.375, rel=1e-9)


def test_pyramidal_efficiency_is_the_two_sectoral_ones_combined(syn):
    """eta_P = (pi^2/8)*eta_E*eta_H. The pi^2/8 removes one copy of the TE10
    cosine taper, which both sectoral efficiencies contain."""
    p = syn("pyramidal_horn", f0=10e9, G_target=20.0)
    assert p.metrics["aperture_efficiency"] == pytest.approx(
        math.pi ** 2 / 8 * p.metrics["eta_e_plane"] * p.metrics["eta_h_plane"],
        rel=1e-12)
    # and its two halves must equal the standalone sectoral horns at the optimum
    e = syn("e_plane_sectoral_horn", f0=10e9, a_wg=0.02286, rho=0.3, flare=1.0)
    h = syn("h_plane_sectoral_horn", f0=10e9, b_wg=0.01016, rho=0.3, flare=1.0)
    assert p.metrics["eta_e_plane"] == pytest.approx(
        e.metrics["aperture_efficiency"], rel=1e-9)
    assert p.metrics["eta_h_plane"] == pytest.approx(
        h.metrics["aperture_efficiency"], rel=1e-9)


def test_horn_efficiency_depends_only_on_the_flare_ratio(syn):
    """A dimensionless quadratic-phase problem cannot care about frequency or
    absolute size, only about how far the flare is from optimum."""
    for rho in (0.1, 0.3, 1.0):
        for f0 in (6e9, 10e9, 18e9):
            d = syn("e_plane_sectoral_horn", f0=f0, a_wg=0.02286, rho=rho, flare=1.0)
            assert d.metrics["aperture_efficiency"] == pytest.approx(
                0.64870263, rel=1e-6), f"rho={rho}, f0={f0}"


def test_circular_patch_directivity_matches_its_own_quadrature(syn):
    """The thin-substrate fit must track the integral it was fitted to, across
    substrates. The value it replaced was a hard-coded 6.3 that was 84% high on
    eps_r 10.2."""
    from scipy.integrate import quad
    from scipy.special import jv

    def exact(k0ae):
        integ = quad(lambda t: ((jv(0, k0ae*math.sin(t)) - jv(2, k0ae*math.sin(t)))**2
                                + math.cos(t)**2
                                * (jv(0, k0ae*math.sin(t)) + jv(2, k0ae*math.sin(t)))**2)
                     * math.sin(t), 0, math.pi/2, limit=400)[0]
        return 4.0/integ

    for eps_r, h, f0 in ((2.2, 1.588e-3, 10e9), (4.4, 1.6e-3, 2.4e9),
                         (10.2, 1.27e-3, 5.8e9)):
        d = syn("circular_patch", f0=f0, eps_r=eps_r, h=h)
        # the free-space edge-current integral is the THIN-substrate directivity
        assert d.metrics["directivity_thin_substrate_linear"] == pytest.approx(
            exact(d.metrics["k0_ae"]), rel=3e-3), f"eps_r={eps_r}"


def test_circular_patch_feed_law_is_exactly_one_at_the_edge(syn):
    """R(rho0)/R_edge must be 1 when the probe IS at the edge, and it must not
    depend on permittivity - k*a_e is 1.8412 at resonance for every substrate."""
    for eps_r, h, f0 in ((2.2, 1.588e-3, 10e9), (10.2, 1.27e-3, 5.8e9)):
        edge = syn("circular_patch", f0=f0, eps_r=eps_r, h=h, rho_frac=1.0)
        assert edge.metrics["feed_resistance_ratio"] == pytest.approx(1.0, rel=1e-9)
    ratios = [syn("circular_patch", f0=f, eps_r=e, h=h,
                  rho_frac=0.3).metrics["feed_resistance_ratio"]
              for e, h, f in ((2.2, 1.588e-3, 10e9), (4.4, 1.6e-3, 2.4e9),
                              (10.2, 1.27e-3, 5.8e9))]
    assert max(ratios) - min(ratios) < 1e-12, "permittivity must cancel out"


def test_circular_and_rectangular_patches_track_each_other(syn):
    """Two resonant patches of similar size on the same board must land within a
    few percent, and both must fall with permittivity.

    This test previously asserted that the circular patch is the LESS directive
    of the two. That was not physics - it was an artifact of comparing two
    invented constants, 6.3 against 6.6. With both integrated from their own
    cavity models they agree to within 2% and cross over around eps_r 6, which
    is well inside either model's accuracy. Asserting an ordering here would be
    asserting numerical noise.
    """
    previous = None
    for eps_r, h, f0 in ((2.2, 1.588e-3, 10e9), (4.4, 1.6e-3, 2.4e9),
                         (10.2, 1.27e-3, 5.8e9)):
        circ = syn("circular_patch", f0=f0, eps_r=eps_r, h=h)
        rect = syn("rectangular_patch", f0=f0, eps_r=eps_r, h=h)
        dc = circ.metrics["directivity_linear"]
        dr = rect.metrics["directivity_linear"]
        assert abs(dc/dr - 1) < 0.05, f"eps_r={eps_r}: circular {dc}, rectangular {dr}"
        if previous is not None:
            assert dc < previous[0] and dr < previous[1], (
                "both must fall as permittivity rises and the patch shrinks")
        previous = (dc, dr)


def test_the_two_slot_model_leaves_out_the_side_walls(syn):
    """The closed form 2*D1/(1 + G12/G1) this spec carried, as a direct 2-D
    integration of the two-slot pattern. The complete cavity model - every
    wall, or the patch current, the same field on a thin board - is above it,
    by most on low-permittivity board: the non-radiating edges cancel in the
    principal planes, not in the total power."""
    import numpy as np

    def sinc(x):
        return np.where(np.abs(x) < 1e-12, 1.0,
                        np.sin(x)/np.where(np.abs(x) < 1e-12, 1, x))

    def grid(k0W, k0Le, n=900):
        th = np.linspace(0, np.pi/2, n)
        ph = np.linspace(0, 2*np.pi, 2*n)
        T, P = np.meshgrid(th, ph, indexing="ij")
        cy = np.sin(T)*np.sin(P)
        slot = np.sqrt(np.maximum(1 - cy**2, 0.0))*sinc(k0W/2*cy)
        u = (slot*2*np.cos(k0Le/2*np.sin(T)*np.cos(P)))**2
        prad = np.trapezoid(np.trapezoid(u*np.sin(T), ph, axis=1), th, axis=0)
        return 4*np.pi*u.max()/prad

    ratios = []
    for eps_r, h, f0 in ((2.2, 1.588e-3, 10e9), (4.4, 1.6e-3, 2.4e9),
                         (10.2, 1.27e-3, 5.8e9)):
        d = syn("rectangular_patch", f0=f0, eps_r=eps_r, h=h)
        k0 = 2*math.pi*f0/2.99792458e8
        ratios.append(d.metrics["directivity_thin_substrate_linear"]/grid(k0*d.get("W"), k0*d.get("L_eff")))
    assert 1.11 > ratios[0] > ratios[1] > ratios[2] > 1.0, ratios


def test_a_narrow_patch_slot_approaches_a_magnetic_dipole(syn):
    """One slot's half-space directivity is (k0 W)^2/I1, which must tend to 3.0
    as the slot narrows - a magnetic dipole's 1.5 doubled by the ground plane.
    That limit is what says the two-slot normalisation is right."""
    from scipy.special import sici

    for k0W in (0.02, 0.1, 0.3):
        i1 = -2 + math.cos(k0W) + k0W*sici(k0W)[0] + math.sin(k0W)/k0W
        assert k0W**2/i1 == pytest.approx(3.0, rel=0.02), f"k0W={k0W}"


def test_conical_horn_uniform_phase_limit_is_the_te11_aperture_value(syn):
    """0.836 for a TE11 circular aperture is in every textbook. Reaching it as
    the flare vanishes is what says the aperture integral is set up right."""
    d = syn("conical_horn", f0=10e9, L=0.3, flare=0.1)
    assert d.metrics["aperture_efficiency"] == pytest.approx(0.8368, abs=0.002)


def test_corrugated_horn_uniform_phase_limit_is_the_published_069(syn):
    """The 0.69 quoted for corrugated horns is the uniform-phase efficiency of
    a J0(2.405 rho/a) taper. Derived, not assumed."""
    d = syn("corrugated_conical_horn", f0=10e9, L=0.3, flare=0.1)
    assert d.metrics["aperture_efficiency"] == pytest.approx(0.6916, abs=0.002)


def test_corrugating_a_horn_does_not_buy_gain_at_the_same_flare(syn):
    """The shipped spec claimed db10(0.69/0.51) = +1.31 dB, which compared the
    corrugated horn at ZERO phase error against the smooth horn at its optimum.
    Like for like, the heavier J0 taper is slightly BEHIND on raw gain near the
    optimum flare - corrugation buys pattern symmetry and cross-polarisation."""
    at_optimum = syn("corrugated_conical_horn", f0=10e9, L=0.3, flare=1.0)
    assert -1.0 < at_optimum.metrics["gain_advantage_over_smooth_db"] < 0.0
    # but the J0 taper tolerates phase error better, so it wins when over-flared
    over = syn("corrugated_conical_horn", f0=10e9, L=0.3, flare=1.3)
    assert over.metrics["gain_advantage_over_smooth_db"] > \
        at_optimum.metrics["gain_advantage_over_smooth_db"]


def test_every_horn_efficiency_falls_away_from_its_optimum_flare(syn):
    """True of all four flared horns now, and expressible by none of them while
    their efficiencies were pinned."""
    for key, extra in (("conical_horn", {"L": 0.3}),
                       ("corrugated_conical_horn", {"L": 0.3})):
        etas = {f: syn(key, f0=10e9, flare=f, **extra).metrics["aperture_efficiency"]
                for f in (0.4, 0.7, 1.0, 1.3)}
        assert etas[0.4] > etas[0.7] > etas[1.0] > etas[1.3], f"{key}: {etas}"


def test_shorted_patch_is_more_than_the_single_slot_of_the_full_patch(syn):
    """A quarter-wave shorted patch has one radiating edge where the full patch
    has two - but its side edges carry in-phase fields and radiate too, so its
    thin-substrate directivity sits BELOW the single slot the spec used to
    carry, by less as the patch shrinks."""
    from scipy.special import sici

    shortfall = []
    for eps_r, h, f0 in ((2.2, 1.588e-3, 10e9), (4.4, 1.6e-3, 2.4e9),
                         (10.2, 1.27e-3, 5.8e9)):
        d = syn("quarter_wave_shorted_patch", f0=f0, eps_r=eps_r, h=h)
        k0w = 2*math.pi*f0*d.get("W")/2.99792458e8
        i1 = -2 + math.cos(k0w) + k0w*sici(k0w)[0] + math.sin(k0w)/k0w
        shortfall.append(1 - d.metrics["directivity_thin_substrate_linear"]/(k0w**2/i1))
        # and the full-patch comparison must match the standalone full patch
        full = syn("rectangular_patch", f0=f0, eps_r=eps_r, h=h)
        assert d.metrics["full_patch_directivity_linear"] == pytest.approx(
            full.metrics["directivity_linear"], rel=1e-9)
    assert 0.03 < shortfall[2] < shortfall[1] < shortfall[0] < 0.2, shortfall


def test_shorting_a_patch_costs_less_than_three_db_and_varies_with_substrate(syn):
    """The spec assumed a flat 3 dB. Two slots less than a half wavelength apart
    never arrayed perfectly, so removing one always costs LESS - and much less
    on high-permittivity board, where they were nearly coincident."""
    losses = {}
    for eps_r, h, f0 in ((2.2, 1.588e-3, 10e9), (4.4, 1.6e-3, 2.4e9),
                         (10.2, 1.27e-3, 5.8e9)):
        d = syn("quarter_wave_shorted_patch", f0=f0, eps_r=eps_r, h=h)
        loss = d.metrics["directivity_lost_to_shorting_db"]
        assert 0.0 < loss < 3.0, f"eps_r={eps_r} lost {loss} dB"
        losses[eps_r] = loss
    assert losses[2.2] > losses[4.4] > losses[10.2], (
        f"the cost must fall as the slots crowd together: {losses}")


@pytest.mark.parametrize("w_over_L", [0.02, 0.05])
def test_a_slot_resonates_where_its_complementary_dipole_does(syn, w_over_L):
    """Babinet, on the one thing it fixes exactly: the resonance.

    Held open as a strict expected failure for a round, because the slot used a
    thin dipole's 0.4785 lambda whatever its width. The complementary dipole has
    radius w/4, and the slot now takes its length AND its resonant resistance
    from resonant_dipole at that radius, so both sides must agree.
    """
    slot = syn("half_wave_slot", f0=F, w_over_L=w_over_L)
    lam = 2.99792458e8 / F
    dip = syn("resonant_dipole", f0=F, aw=slot.parameters["w"] / 4)
    assert slot.parameters["L"] / lam == pytest.approx(
        dip.metrics["length_over_lambda"], rel=0.005)
    assert slot.metrics["resonant_resistance_ohm"] == pytest.approx(
        376.730313412 ** 2 / (4 * dip.metrics["input_resistance_ohm"]), rel=0.005)


def test_a_wider_slot_is_shorter_and_lower_in_resistance(syn):
    """The direction the old flat constants could not show: a wider slot has a
    fatter complement, which resonates shorter and at MORE resistance - so the
    slot, its inverse, sits at LESS."""
    narrow = syn("half_wave_slot", f0=F, w_over_L=0.02)
    wide = syn("half_wave_slot", f0=F, w_over_L=0.05)
    assert wide.parameters["L"] < narrow.parameters["L"]
    assert (wide.metrics["resonant_resistance_ohm"]
            < narrow.metrics["resonant_resistance_ohm"])
