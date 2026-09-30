"""The body-of-revolution solver, checked before it is trusted with a discone.

`otahub.num.bor` solves solid rotationally symmetric surfaces - cones, discs,
tubes - with the exact ring kernel and a finite feed gap. Three checks that
share nothing with each other: a thin tube against the wire code's exact
kernel, which models the same tube; power in against power radiated through an
independent far-field integral; and a long wide bicone against Schelkunoff's
characteristic impedance.
"""
from __future__ import annotations

import math

import pytest

from otahub.num import bor, mom


def _tube(length, a, gap, seg):
    return bor.profile([(a, -length / 2), (a, -gap / 2), (a, gap / 2), (a, length / 2)],
                       seg, gap=(1, 2), gap_seg=2)


def _bicone(slant, half_angle_deg, af=0.005, gap=0.01, seg=0.02):
    th = math.radians(half_angle_deg)
    st, ct = math.sin(th), math.cos(th)
    return bor.profile([(af + slant * st, -gap / 2 - slant * ct), (af, -gap / 2),
                        (af, gap / 2), (af + slant * st, gap / 2 + slant * ct)],
                       seg, gap=(1, 2), gap_seg=2)


@pytest.mark.slow
def test_thin_tube_matches_the_wire_codes_exact_kernel():
    """A tube of radius a carrying a uniform ring current is exactly what the
    wire code's exact kernel models; the two agree to the meshes' own spread."""
    zb = bor.solve(_tube(0.47, 0.002, 0.004, 0.005)).input_impedance
    zw = mom.input_impedance(mom.dipole(0.47, 0.002, 192), exact=True)
    assert zb.real == pytest.approx(zw.real, rel=0.005)
    assert abs(zb.imag - zw.imag) < 1.0


@pytest.mark.parametrize("make", [lambda: _tube(0.47, 0.002, 0.01, 0.02),
                                  lambda: _bicone(0.6, 30.0),
                                  lambda: bor.profile([(0.3, -0.5), (0.01, -0.02), (0.01, 0.02), (0.25, 0.02)],
                                                      0.03, gap=(1, 2), gap_seg=3)])
def test_power_in_is_power_radiated(make):
    sol = bor.solve(make())
    assert bor.radiated_power(sol) == pytest.approx(sol.circuit_power, rel=1e-4)


@pytest.mark.slow
def test_a_long_wide_bicone_presents_its_characteristic_impedance():
    """Schelkunoff: an infinite bicone is a transmission line of impedance
    (eta0/pi) ln cot(theta/2). A wide finite one reflects little at its ends,
    so its resistance must sit on that value within a few percent."""
    zc = 120.0 * math.log(1.0 / math.tan(math.radians(47.0) / 2))
    for slant in (1.25, 2.0):
        z = bor.solve(_bicone(slant, 47.0)).input_impedance
        assert z.real == pytest.approx(zc, rel=0.03)
        assert abs(z.imag) < 0.15 * zc


def test_profile_rejects_the_axis():
    with pytest.raises(ValueError):
        bor.profile([(0.0, -0.2), (0.01, 0.0), (0.01, 0.01), (0.2, 0.01)], 0.02, gap=(1, 2))


def test_the_far_field_takes_a_single_angle():
    sol = bor.solve(_tube(0.47, 0.001, 0.02, 0.02))
    many = bor.far_field(sol, [0.3, math.pi / 2])
    assert bor.far_field(sol, math.pi / 2) == pytest.approx(many[1], rel=1e-14)
    assert bor.far_field(sol, math.pi / 2).shape == ()
