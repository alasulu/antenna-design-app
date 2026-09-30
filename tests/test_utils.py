"""Network and matching maths, checked against Pozar's worked examples."""
import cmath
import math

import numpy as np
import pytest

from otahub.utils import matching as M
from otahub.utils import network as N

POZAR = "Pozar, Microwave Engineering 4th ed."


# ----------------------------------------------------------------- network

def test_matched_load_has_zero_reflection():
    assert N.reflection_coefficient(50, 50) == 0
    assert N.vswr(0) == 1.0
    assert math.isinf(N.return_loss_db(0))


def test_open_and_short_reflect_completely():
    assert abs(N.reflection_coefficient(1e18, 50)) == pytest.approx(1.0, abs=1e-9)
    assert N.reflection_coefficient(0, 50) == -1


@pytest.mark.parametrize("s,gamma,rl_db,ml_db", [
    (1.5, 0.2, 13.9794, 0.1773),
    (2.0, 1 / 3, 9.5424, 0.5115),
    (3.0, 0.5, 6.0206, 1.2494),
])
def test_the_standard_mismatch_table(s, gamma, rl_db, ml_db):
    assert N.gamma_from_vswr(s) == pytest.approx(gamma, rel=1e-6)
    assert N.return_loss_db(gamma) == pytest.approx(rl_db, abs=1e-3)
    assert N.mismatch_loss_db(gamma) == pytest.approx(ml_db, abs=1e-3)


def test_gamma_and_impedance_round_trip():
    for z in (25 + 10j, 100 - 50j, 75, 200 + 300j):
        g = N.reflection_coefficient(z, 50)
        assert N.impedance_from_gamma(g, 50) == pytest.approx(z, rel=1e-12)


def test_vswr_below_one_is_rejected():
    with pytest.raises(ValueError, match="cannot be below 1"):
        N.gamma_from_vswr(0.5)


def test_return_loss_is_positive_and_s11_is_negative():
    """The sign convention that trips people: they are the same number, mirrored."""
    z = 100 + 0j
    assert N.return_loss_db(N.reflection_coefficient(z, 50)) > 0
    assert N.s11_db(z, 50) == pytest.approx(-N.return_loss_db(N.reflection_coefficient(z, 50)))


def test_two_port_conversions_round_trip():
    z = np.array([[60 + 30j, 20], [20, 45 - 10j]])
    assert N.s_to_z(N.z_to_s(z)) == pytest.approx(z, abs=1e-12)
    assert N.abcd_to_z(N.z_to_abcd(z)) == pytest.approx(z, abs=1e-12)
    assert N.y_to_z(N.z_to_y(z)) == pytest.approx(z, abs=1e-12)
    assert N.abcd_to_s(N.z_to_abcd(z)) == pytest.approx(N.z_to_s(z), abs=1e-12)


def test_reciprocal_network_has_unit_abcd_determinant():
    z = np.array([[60 + 30j, 20], [20, 45 - 10j]])       # symmetric => reciprocal
    assert abs(np.linalg.det(N.z_to_abcd(z))) == pytest.approx(1.0, rel=1e-12)


def test_cascading_series_then_shunt_matches_hand_calculation():
    abcd = N.cascade(N.series_impedance(25j), N.shunt_admittance(-0.01j))
    z_in = N.abcd_to_z(abcd)[0, 0]
    # Series 25j into shunt -0.01j terminated open: Zin = 25j + 1/(-0.01j)
    assert z_in == pytest.approx(25j + 100j, rel=1e-12)


def test_quarter_wave_line_inverts_the_load():
    assert N.input_impedance(100, 70.710678, math.pi / 2) == pytest.approx(50, rel=1e-6)
    assert N.input_impedance(25, 50, math.pi / 2) == pytest.approx(100, rel=1e-12)


def test_half_wave_line_repeats_the_load():
    for z in (30 + 40j, 12, 200 - 75j):
        assert N.input_impedance(z, 50, math.pi) == pytest.approx(z, rel=1e-9)


def test_an_open_circuit_load_is_handled_as_a_limit():
    """An open circuit reflects +1, and a line in front of it is -j Z0 cot(bl);
    a load at the line's pole gives an open circuit, not a ZeroDivisionError."""
    open_load = N.impedance_from_gamma(1)
    assert N.reflection_coefficient(open_load) == 1
    assert N.input_impedance(open_load, 50, math.pi / 4) == pytest.approx(-50j)
    assert N.input_impedance(open_load, 50, math.pi / 2) == 0
    pole = 1j * 50 / math.tan(0.2)
    assert math.isinf(N.input_impedance(pole, 50, 0.2).real)


def test_quarter_wave_of_a_short_is_an_open():
    assert math.isinf(N.input_impedance(0, 50, math.pi / 2).real)


def test_non_square_network_is_rejected():
    with pytest.raises(ValueError, match="2x2"):
        N.z_to_s(np.eye(3))


# ---------------------------------------------------------------- matching

def test_l_section_reproduces_pozar_example_5_1():
    """200 - j100 to 100 ohm at 500 MHz: C=0.92pF/L=38.8nH and C=2.61pF/L=46.1nH."""
    sections = M.l_section(200 - 100j, 100.0, 500e6)
    assert len(sections) == 2
    values = sorted((s.shunt.kind, round(s.shunt.value * 1e12, 2) if s.shunt.kind == "C"
                     else round(s.shunt.value * 1e9, 1)) for s in sections)
    assert ("C", 0.92) in values
    assert ("L", 46.1) in values
    series = sorted((s.series.kind, round(s.series.value * 1e9, 1) if s.series.kind == "L"
                     else round(s.series.value * 1e12, 2)) for s in sections)
    assert ("L", 39.0) in series
    assert ("C", 2.6) in series


@pytest.mark.parametrize("z_load", [
    200 - 100j, 25 + 15j, 10 - 5j, 150 + 200j, 75 + 0j, 5 + 1j, 300 - 400j, 20 - 60j])
@pytest.mark.parametrize("z0", [50.0, 75.0, 100.0])
def test_every_l_section_solution_actually_matches(z_load, z0):
    """Forward-evaluate each network: it must land on Z0 exactly.

    This is the test that caught two sign errors - an inverted sqrt(RL/Z0) in
    the RL > Z0 branch, and a wrongly-flipped susceptance sign in the other.
    """
    for section in M.l_section(z_load, z0, 1e9):
        assert section.achieved == pytest.approx(complex(z0, 0.0), abs=1e-9 * z0)


def test_l_section_covers_both_branches():
    assert M.l_section(200 + 0j, 50.0, 1e9)[0].topology.startswith("shunt")
    assert M.l_section(20 + 0j, 50.0, 1e9)[0].topology.startswith("series")


def test_a_low_resistance_load_can_take_both_topologies():
    """25 + j50 to 50 ohm: RL < Z0, so the textbook picks series-first, but
    its conductance (0.01 S) is below 1/Z0 and shunt-first matches too - four
    solutions, found here by root-finding on each topology independently."""
    from scipy.optimize import brentq
    z = 25 + 50j
    sections = M.l_section(z, 50.0)
    assert [s.topology.split(",")[0] for s in sections] == [
        "series at load", "series at load", "shunt at load", "shunt at load"]
    b = [brentq(lambda b: (1 / (1 / z + 1j * b)).real - 50, lo, hi)
         for lo, hi in ((0.001, 0.015), (0.015, 0.05))]
    x = sorted(-(1 / (1 / z + 1j * bb)).imag for bb in b)
    shunt_x = sorted(s.series.reactance for s in sections[2:])
    assert shunt_x == pytest.approx(x, rel=1e-9)
    assert shunt_x == pytest.approx([-61.23724357, 61.23724357], rel=1e-8)


def test_a_one_element_match_is_kept():
    """50 + j50 to 50 ohm needs only a -j50 series capacitor; the shunt
    susceptance of that solution is zero and it used to be dropped."""
    sections = M.l_section(50 + 50j, 50.0)
    one = [s for s in sections if s.shunt.value == 0.0]
    assert len(one) == 1 and one[0].series.reactance == pytest.approx(-50.0)
    assert abs(one[0].achieved - 50) < 1e-9
    assert len(M.l_section(50 + 0j, 50.0)) == 1          # nothing to do, once


def test_negative_load_resistance_is_rejected():
    with pytest.raises(ValueError, match="positive"):
        M.l_section(-10 + 5j, 50.0, 1e9)


def test_single_stub_reproduces_pozar_example_5_2():
    """60 - j80 to 50 ohm: d=0.110λ with l=0.095λ, and d=0.260λ with l=0.405λ."""
    solutions = M.single_stub(60 - 80j, 50.0, "short")
    pairs = [(round(s.line_length_lambda, 3), round(s.stub_length_lambda, 3))
             for s in solutions]
    assert (0.110, 0.095) in pairs
    assert (0.259, 0.405) in pairs or (0.260, 0.405) in pairs


@pytest.mark.parametrize("kind", ["short", "open"])
@pytest.mark.parametrize("z_load", [60 - 80j, 25 + 40j, 120 - 30j, 15 + 5j])
def test_every_stub_solution_actually_matches(kind, z_load):
    solutions = M.single_stub(z_load, 50.0, kind)
    assert solutions, f"no {kind} stub solution for {z_load}"
    for s in solutions:
        assert abs(N.reflection_coefficient(s.achieved, 50.0)) < 1e-3


def test_open_and_short_stubs_differ_by_a_quarter_wavelength():
    shorts = M.single_stub(60 - 80j, 50.0, "short")
    opens = M.single_stub(60 - 80j, 50.0, "open")
    for a, b in zip(shorts, opens):
        assert a.line_length_lambda == pytest.approx(b.line_length_lambda, abs=1e-3)
        delta = abs(a.stub_length_lambda - b.stub_length_lambda)
        assert delta == pytest.approx(0.25, abs=1e-3)


def test_quarter_wave_transformer_is_the_geometric_mean():
    assert M.quarter_wave(100, 50) == pytest.approx(math.sqrt(5000))


def test_quarter_wave_refuses_a_reactive_load():
    """Applying it to a complex load is a quiet, common mistake."""
    with pytest.raises(ValueError, match="real load"):
        M.quarter_wave(100 + 30j, 50)


def test_quarter_wave_bandwidth_narrows_as_the_ratio_grows():
    widths = [M.quarter_wave_bandwidth(z, 50.0) for z in (75, 150, 300, 600)]
    assert widths == sorted(widths, reverse=True)


def test_element_from_reactance_picks_the_right_component():
    assert M.element_from_reactance(100.0, 1e9).kind == "L"
    assert M.element_from_reactance(-100.0, 1e9).kind == "C"
    assert M.element_from_susceptance(0.01, 1e9).kind == "C"
    assert M.element_from_susceptance(-0.01, 1e9).kind == "L"


def test_inductor_value_matches_its_reactance():
    el = M.element_from_reactance(122.474, 500e6)
    assert el.value * 2 * math.pi * 500e6 == pytest.approx(122.474, rel=1e-9)
