"""Touchstone reading and writing.

The tests that matter here are the format's two traps: two-port files store
their matrix column-major while every other size is row-major, and a frequency
point may wrap across any number of lines. Both are silent failures - a
transposed two-port and a mis-chunked multi-port both parse without complaint
and give wrong answers.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from otahub.utils.touchstone import (Network, compare_to_prediction,
                                     read_touchstone, write_touchstone)


# ------------------------------------------------------- the two-port ordering

TWO_PORT_RI = """! deliberately asymmetric, so a transpose cannot hide
# GHz S RI R 50
1.0   0.10 0.11   0.20 0.21   0.30 0.31   0.40 0.41
2.0   0.50 0.51   0.60 0.61   0.70 0.71   0.80 0.81
"""


def test_two_port_files_are_column_major():
    """`freq S11 S21 S12 S22`. Reading it row-major transposes the network,
    swapping forward gain with reverse isolation - which for an amplifier is
    the difference between 20 dB of gain and 20 dB of isolation."""
    n = read_touchstone(TWO_PORT_RI, n_ports=2)
    assert n.s[0, 0, 0] == pytest.approx(0.10 + 0.11j)   # S11: first pair
    assert n.s[0, 1, 0] == pytest.approx(0.20 + 0.21j)   # S21: SECOND pair
    assert n.s[0, 0, 1] == pytest.approx(0.30 + 0.31j)   # S12: THIRD pair
    assert n.s[0, 1, 1] == pytest.approx(0.40 + 0.41j)   # S22: fourth


def test_version_1_1_can_override_the_two_port_order():
    text = TWO_PORT_RI.replace("# GHz", "[Two-Port Data Order] 21_12\n# GHz")
    n = read_touchstone(text, n_ports=2)
    assert n.s[0, 1, 0] == pytest.approx(0.30 + 0.31j), "21_12 means row-major"
    assert n.s[0, 0, 1] == pytest.approx(0.20 + 0.21j)


def test_three_port_is_row_major_not_column_major():
    """The exception applies to two ports ONLY."""
    rows = []
    for f in (1.0, 2.0):
        cells = [f"{i/100:.2f} 0.0" for i in range(9)]
        rows.append(f"{f} " + " ".join(cells[:3]))
        rows.append("  " + " ".join(cells[3:6]))
        rows.append("  " + " ".join(cells[6:9]))
    n = read_touchstone("# GHz S RI R 50\n" + "\n".join(rows), n_ports=3)
    assert n.n_ports == 3
    assert n.s[0, 0, 1] == pytest.approx(0.01)   # second value is row 0, col 1
    assert n.s[0, 1, 0] == pytest.approx(0.03)   # fourth value is row 1, col 0


# --------------------------------------------------------------- line wrapping

def test_a_frequency_point_may_wrap_across_lines():
    """Real instruments wrap. Parsing line-by-line works until it doesn't."""
    tight = read_touchstone(TWO_PORT_RI, n_ports=2)
    wrapped = read_touchstone("""# GHz S RI R 50
1.0   0.10 0.11
      0.20 0.21   0.30 0.31
      0.40 0.41
2.0
  0.50 0.51 0.60 0.61
  0.70 0.71
  0.80 0.81
""", n_ports=2)
    assert np.allclose(wrapped.s, tight.s)
    assert np.allclose(wrapped.frequency_hz, tight.frequency_hz)


# -------------------------------------------------------------- option parsing

@pytest.mark.parametrize("unit,scale", [("Hz", 1.0), ("kHz", 1e3),
                                        ("MHz", 1e6), ("GHz", 1e9)])
def test_every_frequency_unit_scales(unit, scale):
    n = read_touchstone(f"# {unit} S RI R 50\n1.0 0.1 0.0\n2.0 0.2 0.0\n",
                        n_ports=1)
    assert n.frequency_hz[0] == pytest.approx(scale)


def test_magnitude_angle_and_decibel_formats_agree():
    ri = read_touchstone("# GHz S RI R 50\n1.0 0.3 0.4\n", n_ports=1)
    mag, ang = 0.5, math.degrees(math.atan2(0.4, 0.3))
    ma = read_touchstone(f"# GHz S MA R 50\n1.0 {mag} {ang}\n", n_ports=1)
    db = read_touchstone(
        f"# GHz S DB R 50\n1.0 {20*math.log10(mag)} {ang}\n", n_ports=1)
    assert ma.s[0, 0, 0] == pytest.approx(ri.s[0, 0, 0], abs=1e-12)
    assert db.s[0, 0, 0] == pytest.approx(ri.s[0, 0, 0], abs=1e-9)


def test_reference_impedance_is_honoured():
    for text, want in (("# GHz S RI R 50\n1.0 0.0 0.0\n", 50.0),
                       ("# GHz S RI R 75\n1.0 0.0 0.0\n", 75.0)):
        n = read_touchstone(text, n_ports=1)
        assert n.z0 == want
        # S11 = 0 means the port is matched to its OWN reference
        assert n.impedance_at_port(0)[0] == pytest.approx(want)


def test_missing_option_line_uses_the_specified_defaults_and_says_so():
    n = read_touchstone("1.0 0.5 0.0\n2.0 0.4 0.0\n", n_ports=1)
    assert n.frequency_hz[0] == pytest.approx(1e9), "default unit is GHz"
    assert n.z0 == 50.0
    assert any("assumed" in c for c in n.comments), (
        "a guessed option line must be recorded, not silently applied")


def test_comments_are_kept():
    n = read_touchstone("! measured 2026-09-17\n# GHz S RI R 50\n1.0 0.1 0.0\n",
                        n_ports=1)
    assert any("measured 2026-09-17" in c for c in n.comments)


def test_port_count_comes_from_the_filename_suffix(tmp_path):
    f = tmp_path / "device.s2p"
    f.write_text(TWO_PORT_RI)
    assert read_touchstone(f).n_ports == 2


def test_port_count_is_inferred_when_nothing_declares_it():
    assert read_touchstone(TWO_PORT_RI).n_ports == 2
    assert read_touchstone("# GHz S RI R 50\n1.0 0.1 0.2\n").n_ports == 1


# ------------------------------------------------------------------ round trip

@pytest.mark.parametrize("fmt", ["ri", "ma", "db"])
@pytest.mark.parametrize("ports", [1, 2, 3])
def test_write_then_read_returns_what_went_in(fmt, ports):
    rng = np.random.default_rng(7)
    freq = np.array([1e9, 2e9, 3.5e9])
    s = (rng.uniform(-0.8, 0.8, (3, ports, ports))
         + 1j * rng.uniform(-0.8, 0.8, (3, ports, ports)))
    original = Network(freq, s, z0=50.0)
    back = read_touchstone(write_touchstone(original, fmt=fmt), n_ports=ports)
    assert np.allclose(back.frequency_hz, freq)
    assert np.allclose(back.s, s, atol=1e-7), f"{fmt} round trip at {ports} ports"


def test_round_trip_survives_a_nonstandard_reference_impedance():
    net = Network(np.array([1e9]), np.array([[[0.2 + 0.3j]]]), z0=75.0)
    back = read_touchstone(write_touchstone(net), n_ports=1)
    assert back.z0 == 75.0
    assert back.s[0, 0, 0] == pytest.approx(0.2 + 0.3j)


# --------------------------------------------------------------- derived views

def test_impedance_reflection_and_vswr_are_mutually_consistent():
    z_true = 73.1 + 42.5j
    z0 = 50.0
    gamma = (z_true - z0) / (z_true + z0)
    n = read_touchstone(
        f"# GHz S RI R 50\n1.0 {gamma.real} {gamma.imag}\n", n_ports=1)
    assert n.impedance_at_port(0)[0] == pytest.approx(z_true)
    assert n.s_db(0)[0] == pytest.approx(20 * math.log10(abs(gamma)))
    assert n.vswr(0)[0] == pytest.approx((1 + abs(gamma)) / (1 - abs(gamma)))


def test_one_port_z_matches_the_network_module():
    z_true = 30 - 20j
    gamma = (z_true - 50) / (z_true + 50)
    n = read_touchstone(f"# GHz S RI R 50\n1.0 {gamma.real} {gamma.imag}\n",
                        n_ports=1)
    assert n.z[0, 0, 0] == pytest.approx(z_true, rel=1e-9)


def test_resonances_finds_the_dips_and_ignores_shallow_ones():
    freq = np.linspace(1e9, 3e9, 201)
    # a deep dip at 2 GHz and a shallow one at 2.6
    mag = 1 - 0.95 * np.exp(-((freq - 2e9) / 4e7) ** 2) \
            - 0.3 * np.exp(-((freq - 2.6e9) / 4e7) ** 2)
    s = mag.reshape(-1, 1, 1).astype(complex)
    net = Network(freq, s)
    found = net.resonances(threshold_db=-10.0)
    assert len(found) == 1
    assert found[0] == pytest.approx(2e9, rel=1e-3)


def test_interpolation_refuses_to_extrapolate():
    n = read_touchstone(TWO_PORT_RI, n_ports=2)
    n.at(1.5e9)                       # inside, fine
    with pytest.raises(ValueError, match="outside"):
        n.at(5e9)


# -------------------------------------------------------------- failure modes

def test_a_truncated_file_is_rejected_not_silently_padded():
    with pytest.raises(ValueError, match="whole number"):
        read_touchstone("# GHz S RI R 50\n1.0 0.1 0.2 0.3\n", n_ports=2)


def test_the_wrong_port_count_is_caught():
    with pytest.raises(ValueError, match="whole number"):
        read_touchstone(TWO_PORT_RI, n_ports=3)


def test_an_empty_file_is_rejected():
    with pytest.raises(ValueError, match="no data"):
        read_touchstone("! nothing but a comment\n# GHz S RI R 50\n", n_ports=1)


def test_unsorted_frequencies_are_rejected():
    with pytest.raises(ValueError, match="increase"):
        read_touchstone("# GHz S RI R 50\n2.0 0.1 0\n1.0 0.2 0\n", n_ports=1)


def test_an_unrecognised_option_token_is_reported():
    with pytest.raises(ValueError, match="unrecognised"):
        read_touchstone("# GHz S RI R 50 WAT\n1.0 0.1 0\n", n_ports=1)


def test_g_and_h_parameter_files_say_they_are_unsupported():
    with pytest.raises(ValueError, match="not implement"):
        read_touchstone("# GHz H RI R 50\n1.0 0.1 0 0.2 0 0.3 0 0.4 0\n",
                        n_ports=2)


def test_z_parameter_files_convert():
    """Z-parameter files are rarer but legal, and the conversion must land
    back on the same impedance."""
    n = read_touchstone("# GHz Z RI R 50\n1.0 73.1 42.5\n", n_ports=1)
    assert n.impedance_at_port(0)[0] == pytest.approx(73.1 + 42.5j, rel=1e-9)


# ------------------------------------------------------- against an archetype

def test_comparison_puts_a_measurement_beside_a_prediction(registry):
    """The point of the module: an archetype predicts, a solver or VNA
    measures, and the two get stated in the same terms."""
    design = registry["half_wave_dipole"].synthesize(f0=300e6, aw=1e-4)
    z_pred = complex(design.metrics["radiation_resistance_ohm"],
                     design.metrics["input_reactance_ohm"])
    gamma = (z_pred - 50) / (z_pred + 50)
    # a "measurement" 5% high in resistance
    z_meas = complex(z_pred.real * 1.05, z_pred.imag)
    g_meas = (z_meas - 50) / (z_meas + 50)
    net = read_touchstone(
        f"# MHz S RI R 50\n250 {gamma.real} {gamma.imag}\n"
        f"300 {g_meas.real} {g_meas.imag}\n"
        f"350 {gamma.real} {gamma.imag}\n", n_ports=1)
    out = compare_to_prediction(net, z_pred, 300e6)
    assert out["resistance_error_pct"] == pytest.approx(5.0, abs=0.01)
    assert out["reactance_error_ohm"] == pytest.approx(0.0, abs=1e-6)
    assert out["measured_vswr"] > 1.0
    assert out["predicted_s11_db"] == pytest.approx(
        20 * math.log10(abs(gamma)), abs=1e-9)


def test_a_genuinely_ambiguous_file_says_so_rather_than_guessing():
    """99 values with a nine-wide first line is 11 two-port points or 3
    four-port points. Both divide evenly, so there is nothing to choose
    between them and the reader must say so."""
    lines = ["# GHz S RI R 50"]
    for i in range(11):
        lines.append(f"{i+1}.0 " + " ".join(["0.1 0.0"] * 4))
    with pytest.raises(ValueError, match="2 and 4 ports"):
        read_touchstone("\n".join(lines))
    # naming it resolves the ambiguity, which is what the suffix is for
    assert read_touchstone("\n".join(lines), n_ports=2).n_ports == 2
