"""Touchstone reading and writing.

The tests that matter here are the format's traps: two-port files store
their matrix column-major while every other size is row-major, a frequency
point may wrap across any number of lines, Version 1 Z and Y data are
normalised while Version 2's are not, and noise data may follow a two-port's
network data. Both are silent failures - a
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


def _v2(order: str) -> str:
    return ("[Version] 2.0\n# GHz S RI R 50\n[Number of Ports] 2\n"
            f"[Two-Port Data Order] {order}\n[Number of Frequencies] 2\n"
            "[Network Data]\n" + TWO_PORT_RI.split("R 50\n", 1)[1] + "[End]\n")


def test_version_2_names_the_two_port_order():
    """Touchstone 2.1: `21_12` is Version 1's own order, N11 N21 N12 N22;
    `12_21` is the row-major alternative. Reading either backwards swaps S21
    with S12 without a complaint."""
    same = read_touchstone(_v2("21_12"))
    assert np.allclose(same.s, read_touchstone(TWO_PORT_RI, n_ports=2).s)
    rows = read_touchstone(_v2("12_21"))
    assert rows.s[0, 0, 1] == pytest.approx(0.20 + 0.21j), "12_21 is row-major"
    assert rows.s[0, 1, 0] == pytest.approx(0.30 + 0.31j)


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


def test_version_1_z_and_y_data_are_normalised():
    """In a Version 1 file Z and Y are normalised to the option line's R, so
    `1 0` is a matched load, not a one-ohm short."""
    z = read_touchstone("# GHz Z RI R 50\n1.0 1.462 0.85\n", n_ports=1)
    assert z.impedance_at_port(0)[0] == pytest.approx(73.1 + 42.5j, rel=1e-9)
    for param in ("Z", "Y"):
        matched = read_touchstone(f"# GHz {param} RI R 50\n1 1 0\n", n_ports=1)
        assert abs(matched.s[0, 0, 0]) < 1e-12, param


def test_the_specifications_own_z_example_reads_the_same_in_both_versions():
    """Examples 10 and 11 of Touchstone 2.1 are one network twice: Version 1
    normalised to R 75, Version 2 in ohms (where [Reference] 20 has no say
    over Z data)."""
    v1 = read_touchstone("""! 1-port Z-parameter file, multiple frequency points
# MHz Z MA R 75
100    0.99   -4
200    0.80   -22
300    0.707  -45
400    0.40   -62
500    0.01   -89
""", n_ports=1)
    v2 = read_touchstone("""[Version] 2.1
# MHz Z MA
[Number of Ports] 1
[Number of Frequencies] 5
[Reference] 20.0
[Network Data]
100    74.25    -4
200    60      -22
300    53.025  -45
400    30      -62
500     0.75   -89
[End]
""")
    assert v2.z0 == 20.0
    assert np.allclose(v1.impedance_at_port(0), v2.impedance_at_port(0), rtol=1e-12)
    assert v1.impedance_at_port(0)[0] == pytest.approx(
        74.25 * np.exp(-1j * np.radians(4)), rel=1e-12)


def test_version_2_references_span_lines_and_renormalise():
    """[Reference] may continue on the next line. Unequal per-port references
    are renormalised to one, and the renormalised S must describe the same
    network: checked through Z = sqrt(R)(I+S)(I-S)^-1 sqrt(R), written out."""
    s = np.array([[0.2 + 0.1j, 0.5 - 0.2j], [0.5 - 0.2j, -0.1 + 0.3j]])
    cells = " ".join(f"{float(v.real)!r} {float(v.imag)!r}"
                     for v in (s[0, 0], s[1, 0], s[0, 1], s[1, 1]))
    net = read_touchstone(f"""[Version] 2.0
# GHz S RI R 50
[Number of Ports] 2
[Two-Port Data Order] 21_12
[Reference]
50
25
[Network Data]
1.0 {cells}
""")
    assert net.z0 == 50.0
    r = np.diag(np.sqrt([50.0, 25.0]))
    z = r @ (np.eye(2) + s) @ np.linalg.inv(np.eye(2) - s) @ r
    assert np.allclose(net.z[0], z, rtol=1e-12)
    assert any("renormalised" in c for c in net.comments)


def test_version_2_lower_matrix_format_expands_symmetrically():
    rows = ["1.0 0.11 0", "0.21 0 0.22 0", "0.31 0 0.32 0 0.33 0"]
    net = read_touchstone("[Version] 2.0\n# GHz S RI R 50\n[Number of Ports] 3\n"
                          "[Matrix Format] Lower\n[Network Data]\n"
                          + "\n".join(rows) + "\n[End]\n")
    assert net.s[0, 0, 2] == pytest.approx(0.31) and net.s[0, 2, 0] == pytest.approx(0.31)
    assert net.s[0, 1, 2] == pytest.approx(0.32) and net.s[0, 2, 2] == pytest.approx(0.33)


NOISE_V1 = """! 2-port network, S-parameter and noise data
! Default MA format, GHz frequencies, 50-ohm reference, S-parameters
#
! NETWORK PARAMETERS
2  0.95  -26  3.57 157 0.04 76 0.66 -14
22 0.60 -144  1.30  40 0.14 40 0.56 -85
! NOISE PARAMETERS
4  0.7 0.64  69 0.38
18 2.7 0.46 -33 0.40
"""


def test_version_1_noise_data_follows_the_network_data():
    """Example 19 of Touchstone 2.1: the noise block starts where frequency
    stops increasing, five values a line, Rn normalised to R."""
    for net in (read_touchstone(NOISE_V1, n_ports=2), read_touchstone(NOISE_V1)):
        assert net.n_ports == 2 and len(net.frequency_hz) == 2
        assert net.s[1, 1, 0] == pytest.approx(1.30 * np.exp(1j * np.radians(40)))
        assert net.noise.shape == (2, 4)
        assert net.noise[0, 0].real == pytest.approx(4e9)
        assert net.noise[1, 1].real == pytest.approx(2.7)
        assert net.noise[0, 2] == pytest.approx(0.64 * np.exp(1j * np.radians(69)))
        assert net.noise[0, 3].real == pytest.approx(0.38 * 50)


def test_version_2_noise_data_is_marked_and_not_normalised():
    """Example 18 of Touchstone 2.1, with its [Reference] 50 25.0."""
    net = read_touchstone("""[Version] 2.1
#
[Number of Ports] 2
[Two-Port Data Order] 21_12
[Number of Frequencies] 2
[Number of Noise Frequencies] 2
[Reference] 50 25.0
[Network Data]
2  0.95  -26 3.57 157 0.04 76 0.66 -14
22 0.60 -144 1.30  40 0.14 40 0.56 -85
[Noise Data]
4  0.7 0.64  69 19
18 2.7 0.46 -33 20
[End]
""")
    assert net.noise[:, 3].real == pytest.approx([19.0, 20.0])
    assert net.z0 == 50.0


def test_close_frequencies_stay_distinct():
    """1 GHz and 1 GHz + 1 Hz written to nine figures in GHz were the same
    line twice, and the file then failed to read back."""
    net = Network([1e9, 1e9 + 1], np.zeros((2, 1, 1), complex))
    back = read_touchstone(write_touchstone(net), n_ports=1)
    assert np.allclose(back.frequency_hz, net.frequency_hz, rtol=1e-15, atol=0)


def test_five_ports_and_up_write_at_most_four_pairs_a_line():
    """Version 1 allows four pairs a line; a longer matrix row continues on
    the next line, and each row starts a new one."""
    rng = np.random.default_rng(3)
    s = (rng.normal(size=(2, 6, 6)) + 1j * rng.normal(size=(2, 6, 6))) * 0.1
    net = Network([1e9, 2e9], s)
    text = write_touchstone(net)
    data = [ln for ln in text.splitlines() if ln and ln[0] not in "!#"]
    widths = [len(ln.split()) for ln in data]
    assert max(widths) <= 9                           # frequency + four pairs
    assert widths[:12] == [9, 4] + [8, 4] * 5         # 4 + 2 pairs per row
    back = read_touchstone(text, n_ports=6)
    assert np.allclose(back.s, s, rtol=1e-8, atol=1e-12)
    assert read_touchstone(text).n_ports == 6, "the continuation line's width"


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
