"""Computed bandwidth for the loop family, and what could not be computed.

Four archetypes carried `1/(Q_meas*sqrt(2))` with Q simply asserted. Three of
them are solvable and now are. The fourth, `alford_loop`, is not, for a reason
worth stating rather than papering over: its bandwidth is set by the corner
loading network that forces the current uniform, and that network is not part
of the geometry the spec describes. It keeps its assumed Q and says so.

The circular loop is checked three ways - the MoM's direct VSWR walk, the MoM's
Q-derivative, and the independent Fourier-mode solver's own walk. Where a modal
solution exists there is no reason not to use it, and the three agreeing to
under 1% is a stronger statement than any one of them.
"""
from __future__ import annotations

import functools
import math

import numpy as np
import pytest

from otahub.num import loop_modal, mom

NSEG = 48
LAM = 1.0


def _bw_q(zfn):
    s0 = mom.resonant_scale(zfn, 0.9, 1.1)
    z = zfn(s0)
    bw, _ = mom.vswr_bandwidth(zfn, s0, z.real, step=8e-3)
    return bw, mom.antenna_q(zfn, s0)


@functools.lru_cache(maxsize=16)
def _circle(a):
    cres, _ = loop_modal.resonant_circumference(a, 60, 26)
    mom_bw, q = _bw_q(lambda s: mom.solve(
        mom.loop(cres * s, a * s, NSEG), 0).input_impedance)
    modal_bw, _ = _bw_q(lambda s: loop_modal.input_impedance(cres * s, a * s, 60))
    return mom_bw, modal_bw, q


def _square(perimeter, a, per_side=14):
    side = perimeter / 4.0
    r = side / (2.0 * math.sin(math.pi / 4))
    ang = np.arange(4) * math.pi / 2 + math.pi / 4
    verts = np.stack([r * np.cos(ang), r * np.sin(ang), np.zeros(4)], axis=1)
    pts = []
    for i in range(4):
        p, q = verts[i], verts[(i + 1) % 4]
        for j in range(per_side):
            pts.append(p + (j / per_side) * (q - p))
    return (mom.WireModel([mom.Wire(np.array(pts), a, closed=True)]),
            (per_side // 2) - 1)


def _design(registry, key, **kw):
    return registry[key].synthesize(f0=3e8, **kw)


LAMBDA_M = 2.99792458e8 / 3e8


# ------------------------------------------------------- three routes agree

@pytest.mark.parametrize("a", [1e-4, 1e-3, 3e-3])
def test_the_circular_loops_bandwidth_by_three_routes(a):
    """Two independent SOLVERS and one independent formula."""
    mom_bw, modal_bw, q = _circle(a)
    assert modal_bw == pytest.approx(mom_bw, rel=0.01)
    assert 0.7071 / q == pytest.approx(mom_bw, rel=0.02)


# ------------------------------------------------------ the spec's numbers

@pytest.mark.parametrize("a", [1e-4, 1e-3, 3e-3])
def test_spec_bandwidth_and_q_for_the_circular_loop(a, registry):
    mom_bw, _, q = _circle(a)
    d = _design(registry, "one_wavelength_circular_loop", b=a * LAMBDA_M)
    assert d.metrics["fractional_bandwidth_vswr2"] == pytest.approx(mom_bw, rel=0.02)
    assert d.metrics["quality_factor"] == pytest.approx(q, rel=0.02)


def test_the_assumed_q_was_badly_wrong_for_both_full_wave_loops(registry):
    """Q = 12 and Q = 14 were asserted. The truth runs 4.3 to 10.2, so the old
    bandwidths were up to two thirds too narrow."""
    for key, old_bw in (("one_wavelength_circular_loop", 0.05893),
                        ("quad_loop_square", 0.05051)):
        for a in (1e-4, 3e-3):
            d = _design(registry, key, b=a * LAMBDA_M)
            assert d.metrics["fractional_bandwidth_vswr2"] > old_bw
            assert 4.0 < d.metrics["quality_factor"] < 11.0


def test_the_halos_assumed_q_erred_the_other_way(registry):
    """Q = 25 flattered it. The halo is genuinely higher-Q than that over most
    of its range, so the old figure promised bandwidth it does not have."""
    narrow = registry["halo_loop"].synthesize(
        f0=3e8, g=0.005 * LAMBDA_M, b=0.001 * LAMBDA_M)
    assert narrow.metrics["fractional_bandwidth_vswr2"] < 0.0283
    assert narrow.metrics["quality_factor"] > 25.0


def test_a_halo_is_far_higher_q_than_a_full_wave_loop(registry):
    """Roughly five times, which is the price of the size."""
    halo = registry["halo_loop"].synthesize(
        f0=3e8, g=0.015 * LAMBDA_M, b=0.002 * LAMBDA_M)
    loop = _design(registry, "one_wavelength_circular_loop", b=0.002 * LAMBDA_M)
    assert halo.metrics["quality_factor"] / loop.metrics["quality_factor"] > 4.0


def test_bandwidth_widens_with_the_conductor(registry):
    """The oldest rule there is, now asserted for three archetypes at once."""
    for key, kw in (("one_wavelength_circular_loop", {}),
                    ("quad_loop_square", {}),
                    ("halo_loop", {"g": 0.02 * LAMBDA_M})):
        got = [registry[key].synthesize(f0=3e8, b=a * LAMBDA_M, **kw)
               .metrics["fractional_bandwidth_vswr2"]
               for a in (3e-4, 1e-3, 3e-3)]
        assert got[0] < got[1] < got[2], (key, got)


# --------------------------------------------------- the Alford equivalence

@pytest.mark.parametrize("perimeter", [0.5, 0.75, 1.0])
def test_alfords_equal_area_circle_is_better_than_a_few_percent(perimeter):
    """The spec's own validity note called this 'good to a few percent'. Against
    a uniform-current solve of the actual square it is inside 0.75% on
    resistance and 0.1% on directivity."""
    from scipy.integrate import quad
    from scipy.special import jv
    ka = 2.0 * math.sqrt(math.pi) * (perimeter / 4.0)
    i2 = quad(lambda x: jv(2, x), 0, 2 * ka)[0]
    rr_circle = 376.730313412 * (math.pi / 2.0) * ka * i2
    d_circle = 2.0 * ka * jv(1, ka) ** 2 / i2

    m, feed = _square(perimeter, 2e-3, 20)
    sol = mom.MoMSolution(m, np.ones(m.n_basis, dtype=complex), feed)
    rr_square = 2.0 * mom.radiated_power(sol, 100, 100)
    assert rr_square == pytest.approx(rr_circle, rel=0.01)
    assert mom.directivity(sol, 100, 100) == pytest.approx(d_circle, rel=2e-3)


def test_a_uniform_current_square_has_essentially_no_azimuth_ripple():
    """So the Alford loop's asserted 0.5 dB is ENTIRELY an allowance for
    imperfect current forcing, not a property of the shape. Worth pinning,
    because it tells a builder where the whole ripple budget goes."""
    m, feed = _square(1.0, 2e-3, 20)
    sol = mom.MoMSolution(m, np.ones(m.n_basis, dtype=complex), feed)
    phi = np.linspace(0, 2 * math.pi, 361)
    e_th, e_ph = mom.far_field(sol, np.full(361, math.pi / 2), phi)
    u = np.abs(e_th) ** 2 + np.abs(e_ph) ** 2
    assert 10 * math.log10(u.max() / u.min()) < 0.05


def test_alford_still_says_its_bandwidth_is_assumed(registry):
    """It is the one of the four that cannot be derived, and the spec should go
    on saying so rather than quietly acquiring a fitted number."""
    spec = registry["alford_loop"].spec
    rule = next(r for r in spec.analysis
                if r.metric == "fractional_bandwidth_vswr2")
    assert "Q_meas" in rule.expr
    assert "ASSUMED" in (rule.notes or "")
