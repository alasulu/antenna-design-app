"""The command line: arguments read in the units the specs declare, and the
figures it prints being the ones it claims to print."""
from __future__ import annotations

import math

import numpy as np
import pytest

from otahub.cli.main import _kv, main


def test_an_angle_with_a_unit_stays_in_the_specs_degrees():
    """`flare_deg=90deg` was parsed to 1.5708 (radians) and handed to a spec
    that reads degrees - a different bowtie, silently."""
    units = {"flare_deg": "deg", "L": "m"}
    assert _kv(["flare_deg=90deg"], units) == {"flare_deg": pytest.approx(90.0)}
    assert _kv(["flare_deg=1.5707963267948966rad"], units)["flare_deg"] == pytest.approx(90.0)
    assert _kv(["flare_deg=90"], units)["flare_deg"] == 90.0
    assert _kv(["L=12mm"], units)["L"] == pytest.approx(0.012)
    with pytest.raises(Exception, match="angle"):
        _kv(["L=12deg"], units)


def test_bowtie_flare_in_degrees_either_way(capsys):
    assert main(["synth", "bowtie", "--set", "f_low=1GHz", "--set", "flare_deg=90", "--json"]) == 0
    plain = capsys.readouterr().out
    assert main(["synth", "bowtie", "--set", "f_low=1GHz", "--set", "flare_deg=90deg", "--json"]) == 0
    assert capsys.readouterr().out == plain


def test_conductivity_takes_its_unit_and_bad_values_do_not_traceback(capsys):
    assert main(["synth", "small_circular_loop", "--f0", "14.2MHz",
                 "--set", "sigma=58e6S/m", "--set", "C_over_lambda=0.1"]) == 0
    capsys.readouterr()
    assert main(["synth", "small_circular_loop", "--f0", "14.2MHz", "--set", "sigma=lots"]) == 2
    assert "bad --set" in capsys.readouterr().err


def test_a_filtered_list_counts_its_own_families(capsys):
    assert main(["list", "--family", "wire"]) == 0
    assert "1 family(ies)" in capsys.readouterr().out


def test_check_shows_the_absolute_tolerance_that_decided(capsys):
    assert main(["check", "--key", "resonant_dipole", "-v"]) == 0
    out = capsys.readouterr().out
    line = next(ln for ln in out.splitlines() if "input_reactance_ohm" in ln)
    assert "tol 6 absolute" in line and "off by" in line
    assert "known-case expectations pass" in out


def test_touchstone_comparison_uses_the_complex_driving_point_impedance(tmp_path, capsys, registry):
    """It used to take the radiation resistance and zero reactance - 0.02
    ohm for a small loop whose terminals show 0.19 + j204."""
    d = registry["small_circular_loop"].synthesize(f0=14.2e6, C_over_lambda=0.1, N=1, b=0.002)
    z = complex(d.metrics["input_impedance_ohm"])
    g = (z - 50) / (z + 50)
    f = tmp_path / "loop.s1p"
    f.write_text("# MHz S RI R 50\n" + "".join(
        f"{mhz} {g.real!r} {g.imag!r}\n" for mhz in (14.0, 14.2, 14.4)))
    assert main(["touchstone", str(f), "--compare", "small_circular_loop", "--at", "14.2e6",
                 "--set", "C_over_lambda=0.1", "--set", "N=1", "--set", "b=2mm"]) == 0
    out = capsys.readouterr().out
    assert "input_impedance_ohm" in out
    assert "resistance off by +0.00%" in out


def test_planar_uses_dy_and_the_real_steering_direction(capsys):
    """d = 0.5, dy = 1.1 at broadside has a full grating lobe along y; only d
    was checked. And 0.6 lambda scanned to 45 deg at phi 45 is clean, though
    past the all-azimuth limit."""
    assert main(["planar", "--nx", "4", "--ny", "4", "--d", "0.5", "--dy", "1.1"]) == 0
    assert "A GRATING LOBE IS IN REAL SPACE" in capsys.readouterr().out
    assert main(["planar", "--nx", "4", "--ny", "4", "--d", "0.6", "--scan", "45",
                 "--scan-phi", "45"]) == 0
    out = capsys.readouterr().out
    assert "A GRATING LOBE" not in out and "there is none" in out


def test_thinning_a_uniform_array_does_not_crash(capsys):
    assert main(["planar", "--nx", "4", "--ny", "4", "--thin", "7"]) == 0
    assert "nothing was thinned" in capsys.readouterr().out
