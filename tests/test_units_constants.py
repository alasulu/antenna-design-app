import math

import pytest

from otahub.core import constants as k
from otahub.core import units as u


def test_free_space_impedance_is_consistent_with_eps0_mu0():
    assert k.ETA0 == pytest.approx(math.sqrt(k.MU0 / k.EPS0))
    assert k.ETA0 == pytest.approx(376.730313412, rel=1e-9)


def test_eps0_matches_codata():
    assert k.EPS0 == pytest.approx(8.8541878128e-12, rel=1e-8)


def test_wavelength_and_wavenumber_are_inverse():
    lam = k.wavelength(2.4e9)
    assert lam == pytest.approx(0.12491, rel=1e-4)
    assert k.wavenumber(2.4e9) == pytest.approx(2 * math.pi / lam)


def test_wavelength_shortens_in_dielectric():
    assert k.wavelength(1e9, eps_r=4.0) == pytest.approx(k.wavelength(1e9) / 2.0)


def test_nonpositive_frequency_rejected():
    with pytest.raises(ValueError):
        k.wavelength(0)


def test_copper_skin_depth_at_1ghz():
    # Textbook value: ~2.1 um in copper at 1 GHz.
    d = k.skin_depth(1e9, k.CONDUCTIVITY["copper"])
    assert d == pytest.approx(2.09e-6, rel=0.05)


def test_surface_resistance_relates_to_skin_depth():
    sigma = k.CONDUCTIVITY["copper"]
    f = 1e9
    assert k.surface_resistance(f, sigma) == pytest.approx(
        1.0 / (sigma * k.skin_depth(f, sigma)), rel=1e-9)


@pytest.mark.parametrize("value,unit,si", [
    (1.0, "mm", 1e-3), (1.0, "in", 0.0254), (2.4, "GHz", 2.4e9), (180.0, "deg", math.pi),
])
def test_to_si_round_trip(value, unit, si):
    assert u.to_si(value, unit) == pytest.approx(si)
    assert u.from_si(si, unit) == pytest.approx(value)


def test_db_helpers_round_trip():
    assert u.db10(u.undb10(6.0)) == pytest.approx(6.0)
    assert u.db20(100.0) == pytest.approx(40.0)
    assert u.db10(0) == float("-inf")


def test_cross_kind_conversion_is_refused():
    with pytest.raises(KeyError):
        u.convert(1.0, "mm", "GHz")
