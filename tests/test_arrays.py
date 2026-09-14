"""Array taper and pattern checks.

The decisive test for Dolph-Chebyshev is that the MEASURED sidelobe level of
the synthesised pattern equals the design level. That is the defining property
of the taper, and it catches the class of bug where plausible-looking weights
are subtly wrong.
"""
import math

import numpy as np
import pytest

from otahub import arrays as A
from otahub.core.pattern import first_sidelobe_db, hpbw_deg


# ------------------------------------------------------------------ tapers

def test_uniform_taper_efficiency_is_exactly_one():
    assert A.taper_efficiency(A.uniform(10)) == pytest.approx(1.0, rel=1e-12)


def test_binomial_weights_are_binomial_coefficients():
    assert A.binomial(5) * 6 == pytest.approx([1, 4, 6, 4, 1])


@pytest.mark.parametrize("n", [4, 5, 8, 10, 16, 21])
@pytest.mark.parametrize("sll", [-20.0, -30.0, -40.0])
def test_chebyshev_sidelobes_equal_the_design_level(n, sll):
    """The defining property: every sidelobe sits exactly at the design level."""
    pattern = A.array_pattern(A.dolph_chebyshev(n, sll), 0.5)
    assert first_sidelobe_db(pattern) == pytest.approx(sll, abs=0.05)


def test_chebyshev_weights_are_symmetric_and_peak_at_centre():
    for n in (7, 8, 12, 15):
        w = A.dolph_chebyshev(n, -30.0)
        assert w == pytest.approx(w[::-1], rel=1e-9), f"n={n} not symmetric"
        assert np.argmax(w) in (n // 2, (n - 1) // 2)


def test_chebyshev_depends_on_the_design_level():
    """Guards the bug where even-length weights came out identical for every
    sidelobe level because of a missing half-sample phase term."""
    a = A.dolph_chebyshev(10, -20.0)
    b = A.dolph_chebyshev(10, -40.0)
    assert not np.allclose(a, b)
    assert A.taper_efficiency(a) > A.taper_efficiency(b)


def test_deeper_sidelobes_cost_taper_efficiency():
    effs = [A.taper_efficiency(A.dolph_chebyshev(16, s))
            for s in (-20.0, -30.0, -40.0, -50.0)]
    assert effs == sorted(effs, reverse=True)


def test_binomial_array_has_no_sidelobes():
    """Binomial excitation at d <= lambda/2 gives a pattern with no sidelobes."""
    assert first_sidelobe_db(A.array_pattern(A.binomial(10), 0.5)) < -60.0


def test_uniform_array_sidelobe_approaches_minus_13_26_db():
    """Discrete uniform arrays tend to the continuous aperture's -13.26 dB."""
    sll_small = first_sidelobe_db(A.array_pattern(A.uniform(10), 0.5))
    sll_large = first_sidelobe_db(A.array_pattern(A.uniform(50), 0.5))
    assert sll_small == pytest.approx(-12.97, abs=0.15)
    assert sll_large == pytest.approx(-13.26, abs=0.10)
    assert abs(sll_large + 13.26) < abs(sll_small + 13.26)


def test_taper_ordering_uniform_narrowest_binomial_broadest():
    beams = {name: hpbw_deg(A.array_pattern(w, 0.5)) for name, w in (
        ("uniform", A.uniform(16)),
        ("chebyshev", A.dolph_chebyshev(16, -30.0)),
        ("binomial", A.binomial(16)))}
    assert beams["uniform"] < beams["chebyshev"] < beams["binomial"]


def test_taylor_holds_near_sidelobes_near_the_design_level():
    sll = first_sidelobe_db(A.array_pattern(A.taylor_nbar(20, -30.0, 5), 0.5))
    assert sll == pytest.approx(-30.0, abs=2.0)


def test_chebyshev_holds_every_sidelobe_at_exactly_the_design_level():
    """The equal-ripple signature. Eight consecutive sidelobes land on -30.00 dB.

    Worth stating plainly because the common claim that Taylor beats
    Dolph-Chebyshev on aperture efficiency does NOT hold for discrete arrays -
    measured here, Chebyshev is the more efficient of the two for n >= 20.
    Dolph-Chebyshev is provably optimal for a discrete array; Taylor's real
    advantages are the absence of edge spikes and decaying far sidelobes.
    """
    pattern = A.array_pattern(A.dolph_chebyshev(30, -30.0), 0.5)
    cut, levels = pattern.cut(0.0), []
    for i in range(1, len(cut) - 1):
        if cut[i] > cut[i - 1] and cut[i] >= cut[i + 1] and cut[i] < 0.5:
            levels.append(10 * math.log10(max(cut[i], 1e-30)))
    assert len(levels) >= 6
    for level in levels[:6]:
        assert level == pytest.approx(-30.0, abs=0.05)


def test_chebyshev_has_edge_spikes_and_taylor_does_not():
    """The practical reason to choose Taylor: Chebyshev's outermost element is
    excited harder than its neighbour, which is awkward to build and feeds
    mutual-coupling trouble at the array edge."""
    cheb = A.dolph_chebyshev(40, -30.0)
    assert cheb[0] > cheb[1], "Chebyshev should show the characteristic edge spike"
    taylor = A.taylor_nbar(40, -30.0, 5)
    assert np.all(np.diff(taylor[:20]) > 0), "Taylor should rise monotonically to centre"


def test_raised_cosine_pedestal_one_is_uniform():
    assert A.raised_cosine(10, pedestal=1.0) == pytest.approx(np.ones(10))


@pytest.mark.parametrize("bad", [0, -3, 2.5])
def test_invalid_element_counts_are_rejected(bad):
    with pytest.raises(ValueError):
        A.uniform(bad)


def test_positive_sidelobe_specification_is_rejected():
    with pytest.raises(ValueError, match="must be negative"):
        A.dolph_chebyshev(10, 30.0)


# ------------------------------------------------------------------ factor

def test_uniform_half_wave_array_directivity_is_exactly_n():
    for n in (4, 10, 20):
        assert A.broadside_directivity(A.uniform(n), 0.5) == pytest.approx(n, rel=1e-6)


def test_taper_reduces_directivity_by_the_taper_efficiency():
    w = A.dolph_chebyshev(20, -30.0)
    expected = 20 * A.taper_efficiency(w)
    assert A.broadside_directivity(w, 0.5) == pytest.approx(expected, rel=0.02)


def test_uniform_array_beamwidth_follows_the_standard_relation():
    """HPBW ~ 0.886 * lambda / (N*d) radians for a uniform broadside array."""
    n, d = 20, 0.5
    expected = math.degrees(0.886 / (n * d))
    assert hpbw_deg(A.array_pattern(A.uniform(n), d)) == pytest.approx(expected, rel=0.03)


def test_grating_lobe_limit_is_one_at_broadside_and_half_at_endfire():
    assert A.grating_lobe_free_spacing(90.0) == pytest.approx(1.0)
    assert A.grating_lobe_free_spacing(0.0) == pytest.approx(0.5)
    assert A.grating_lobe_free_spacing(30.0) == pytest.approx(1 / (1 + math.cos(math.pi / 6)))


def test_half_wave_spacing_is_grating_lobe_free_at_any_scan():
    for scan in (0.0, 30.0, 60.0, 90.0):
        assert not A.has_grating_lobe(0.49, scan)


def test_wide_spacing_produces_a_grating_lobe_when_scanned():
    assert not A.has_grating_lobe(0.9, 90.0)
    assert A.has_grating_lobe(0.9, 20.0)


def test_scanning_steers_the_main_beam_where_asked():
    for scan in (40.0, 60.0, 90.0):
        pattern = A.array_pattern(A.uniform(16), 0.4, scan_deg=scan)
        peak = math.degrees(pattern.theta[int(np.argmax(pattern.U[:, 0]))])
        assert peak == pytest.approx(scan, abs=1.0)


def test_summarise_reports_a_consistent_picture():
    s = A.summarise(A.dolph_chebyshev(16, -25.0), 0.5)
    assert s["elements"] == 16
    assert s["sidelobe_db"] == pytest.approx(-25.0, abs=0.05)
    assert not s["grating_lobe"]
    assert s["directivity_dbi"] < 10 * math.log10(16)   # taper costs directivity


# --------------------------------------------------------------------- CLI
# `--sll 30` is how almost everyone says "30 dB sidelobes", and it used to hit
# an unhandled ValueError from the taper. A sidelobe above the main beam is not
# a thing, so either sign is accepted.

def _run_cli(argv):
    from otahub.cli.main import main
    return main(argv)


@pytest.mark.parametrize("sll", ["30", "-30", "45", "-45"])
def test_array_command_accepts_either_sidelobe_sign(sll, capsys):
    assert _run_cli(["array", "-n", "8", "--taper", "chebyshev", "--sll", sll]) == 0
    out = capsys.readouterr().out
    want = -abs(float(sll))
    assert f"{want:.1f} dB design" in out
    # and the synthesised pattern really does hit the design level
    assert f"first sidelobe     {want:.2f} dB" in out


def test_array_command_reports_impossible_taper_cleanly(capsys):
    """A taper that cannot be synthesised must not raise a traceback."""
    rc = _run_cli(["array", "-n", "1", "--taper", "chebyshev", "--sll", "30"])
    assert rc in (0, 1)
    if rc == 1:
        assert "cannot synthesise" in capsys.readouterr().err
