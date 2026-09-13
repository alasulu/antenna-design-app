"""First-principles pattern maths, checked against closed-form textbook values.

These tests are the independent yardstick for the spec-driven analysis
formulas: if a spec claims a directivity that this module contradicts, the
disagreement is real and one of the two is wrong.
"""
import math

import numpy as np
import pytest

from otahub.core import pattern as p


def test_short_dipole_directivity_is_exactly_three_halves():
    assert p.short_dipole().directivity() == pytest.approx(1.5, rel=1e-4)


def test_half_wave_dipole_matches_balanis():
    # Balanis: D0 = 1.6409 (2.15 dBi) for the L = lambda/2 thin dipole.
    d = p.finite_dipole(0.5).directivity()
    assert d == pytest.approx(1.6409, rel=1e-3)
    assert 10 * math.log10(d) == pytest.approx(2.15, abs=0.01)


def test_full_wave_dipole_directivity():
    assert p.finite_dipole(1.0).directivity() == pytest.approx(2.41, rel=1e-2)


def test_half_wave_dipole_hpbw_is_78_degrees():
    assert p.hpbw_deg(p.finite_dipole(0.5)) == pytest.approx(78.0, abs=0.5)


def test_short_dipole_hpbw_is_90_degrees():
    assert p.hpbw_deg(p.short_dipole()) == pytest.approx(90.0, abs=0.5)


@pytest.mark.parametrize("q,expected", [(0.5, 4.0), (1.0, 6.0), (2.0, 10.0), (3.0, 14.0)])
def test_cosine_q_feed_directivity_follows_2_2q_plus_1(q, expected):
    """D = 2(2q+1) for a cos^q field taper over the forward hemisphere."""
    assert p.cosine_q(q).directivity() == pytest.approx(expected, rel=2e-3)


def test_uniform_line_source_sidelobe_is_minus_13_26_db():
    """The canonical uniform-illumination sidelobe level."""
    sll = p.first_sidelobe_db(p.uniform_line_source(10))
    assert sll == pytest.approx(-13.26, abs=0.1)


def test_uniform_line_source_sidelobe_is_independent_of_length():
    """-13.26 dB is a property of uniform illumination, not of aperture size."""
    for L in (5, 10, 20):
        assert p.first_sidelobe_db(p.uniform_line_source(L)) == pytest.approx(-13.26, abs=0.15)


def test_uniform_line_source_directivity_approaches_2L_over_lambda():
    for L in (10, 20):
        assert p.uniform_line_source(L).directivity() == pytest.approx(2 * L, rel=0.02)


def test_uniform_line_source_hpbw_follows_50_8_lambda_over_L():
    assert p.hpbw_deg(p.uniform_line_source(10)) == pytest.approx(5.08, abs=0.05)


def test_radiated_power_of_isotropic_pattern_is_four_pi():
    # Trapezoidal integration of sin(theta) on the default grid converges to
    # 4*pi from below with O(h^2) error (~6e-6 here), so the tolerance tracks
    # the quadrature scheme rather than claiming exactness.
    theta, phi = p.make_grid()
    iso = p.Pattern(theta, phi, np.ones((theta.size, phi.size)))
    assert iso.radiated_power() == pytest.approx(4 * math.pi, rel=1e-4)
    assert iso.directivity() == pytest.approx(1.0, rel=1e-4)


def test_isotropic_quadrature_error_shrinks_as_the_grid_refines():
    """Guards the quadrature itself: refining the grid must improve accuracy."""
    errors = []
    for n in (91, 361, 1441):
        theta, phi = p.make_grid(n, n // 2 + 1)
        iso = p.Pattern(theta, phi, np.ones((theta.size, phi.size)))
        errors.append(abs(iso.radiated_power() - 4 * math.pi))
    assert errors[0] > errors[1] > errors[2], f"not converging: {errors}"


def test_mismatched_grid_is_rejected():
    theta, phi = p.make_grid(11, 7)
    with pytest.raises(ValueError, match="does not match grid"):
        p.Pattern(theta, phi, np.ones((5, 5)))


def test_non_radiating_pattern_has_no_directivity():
    theta, phi = p.make_grid(21, 11)
    with pytest.raises(ValueError, match="radiates no power"):
        p.Pattern(theta, phi, np.zeros((21, 11))).directivity()


# --- engine behaviour that the loop-spec repair depended on -----------------

def test_skin_effect_surface_resistance_two_ways_agree():
    """Rs = sqrt(pi f mu / sigma) must equal 1/(sigma * delta).

    The loop spec shipped two known cases whose loss resistances were wrong by
    30.9x and 3.3x; this identity is what settled which side was right.
    """
    from otahub.core import constants as k
    for f, sigma in ((14.2e6, 5.8e7), (146e6, 3.54e7), (2.4e9, 5.8e7)):
        rs = k.surface_resistance(f, sigma)
        assert rs == pytest.approx(1.0 / (sigma * k.skin_depth(f, sigma)), rel=1e-12)
