"""Single-turn small loops against their real current.

The circular and square small-loop specs used the uniform-current laws out to
a third of a wavelength round, claiming the quartic resistance law was within
2.3% there - a comparison with the uniform-current integral, not with a real
loop. Solved with its actual current (the Fourier-series loop solution for the
circle, the method of moments for the square, which meet each other on the
circle), the loop's resistance is already 12-17% above the law at 0.1
wavelengths and several times it at a third: the current's cos(phi) part
radiates like an electric dipole, growing as the size squared.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from otahub.num import loop_modal, mom

C0 = 2.99792458e8


def _mom_circle(Cl, b, nseg=64):
    a = Cl / (2 * math.pi)
    t = 2 * math.pi * np.arange(nseg) / nseg
    m = mom.WireModel([mom.Wire(np.stack([a * np.cos(t), a * np.sin(t), 0 * t], axis=1), b, closed=True)])
    return mom.solve(m, 0).input_impedance


def _mom_square(P, b, n=None):
    s = P / 4
    if n is None:
        n = max(8, min(40, int(s / (4 * b))))
        n += n % 2
    c = np.array([[-s / 2, -s / 2, 0], [s / 2, -s / 2, 0], [s / 2, s / 2, 0], [-s / 2, s / 2, 0], [-s / 2, -s / 2, 0]])
    pts = [c[0]] + [c[k] + (c[k + 1] - c[k]) * j / n for k in range(4) for j in range(1, n + 1)]
    m = mom.WireModel([mom.Wire(np.array(pts)[:-1], b, closed=True)])
    return mom.solve(m, n // 2 - 1).input_impedance


@pytest.mark.parametrize("Cl", [0.05, 0.1, 0.2])
def test_the_two_loop_solvers_agree(Cl):
    b = Cl / (2 * math.pi) / 100
    z1, z2 = _mom_circle(Cl, b), loop_modal.input_impedance(Cl, b)
    assert z1.real == pytest.approx(z2.real, rel=0.01) and z1.imag == pytest.approx(z2.imag, rel=0.006)


@pytest.mark.parametrize("Cl,ba", [(0.07, 0.004), (0.13, 0.02), (0.24, 0.0015)])
def test_the_circular_loop_driving_point(registry, Cl, ba):
    lam = C0 / 300e6
    d = registry["small_circular_loop"].synthesize(f0=300e6, C=Cl * lam, b=ba * Cl * lam / (2 * math.pi), N=1)
    z = loop_modal.input_impedance(Cl, ba * Cl / (2 * math.pi))
    assert d.metrics["input_resistance_driving_point_ohm"] == pytest.approx(z.real, rel=0.025)
    assert d.metrics["input_reactance_driving_point_ohm"] == pytest.approx(z.imag, rel=0.012)
    assert d.get("C_tune") == pytest.approx(1 / (2 * math.pi * 300e6 * z.imag), rel=0.015)


@pytest.mark.parametrize("P,bs", [(0.08, 0.005), (0.15, 0.02), (0.26, 0.002)])
def test_the_square_loop_driving_point(registry, P, bs):
    lam = C0 / 300e6
    d = registry["small_square_loop"].synthesize(f0=300e6, s=P * lam / 4, b=bs * P * lam / 4, N=1)
    z = _mom_square(P, bs * P / 4)
    assert d.metrics["input_resistance_driving_point_ohm"] == pytest.approx(z.real, rel=0.03)
    assert d.metrics["input_reactance_driving_point_ohm"] == pytest.approx(z.imag, rel=0.015)


def test_the_uniform_current_law_understates_a_tenth_wave_loop_and_fails_at_a_third():
    for ba in (0.001, 0.01, 0.05):
        z = loop_modal.input_impedance(0.1, ba * 0.1 / (2 * math.pi))
        assert 0.11 < z.real / (20 * math.pi ** 2 * 0.1 ** 4) - 1 < 0.18
        z3 = loop_modal.input_impedance(1 / 3, ba / (6 * math.pi))
        assert z3.real / (20 * math.pi ** 2 / 81) > 4.5


def test_past_030_the_driving_point_is_nan(registry):
    lam = C0 / 300e6
    d = registry["small_circular_loop"].synthesize(f0=300e6, C=0.33 * lam, b=0.001, N=1)
    assert math.isnan(d.metrics["input_resistance_driving_point_ohm"]) and math.isnan(d.metrics["radiation_efficiency"])
    assert d.metrics["radiation_resistance_ohm"] == pytest.approx(20 * math.pi ** 2 * 0.33 ** 4, rel=1e-6)


# --------------------------------------------------------------- tuned loop

MU0 = 4e-7 * math.pi
RS = lambda f0: math.sqrt(math.pi * f0 * MU0 / 5.8e7)       # copper


def _tuned_and_swept(zr, rl0):
    """Q_Z and the VSWR-2 band of a loop series-tuned at f0, by sweeping its
    impedance at fixed geometry: zr(k) is the lossless loop at k = f/f0, the skin
    loss grows as sqrt(k) and the capacitor's reactance falls as 1/k."""
    from scipy.optimize import brentq
    x1 = zr(1.0).imag
    z = lambda k: zr(k) + rl0 * math.sqrt(k) - 1j * x1 / k
    r = z(1.0).real
    q = abs(z(1 + 1e-4) - z(1 - 1e-4)) / 2e-4 / (2 * r)
    gam = lambda k: abs((z(k) - r) / (z(k) + r))
    v = lambda k: (1 + gam(k)) / (1 - gam(k)) - 2.0
    return q, brentq(v, 1.0, 1 + 3 / q, xtol=1e-12) - brentq(v, 1 - 3 / q, 1.0, xtol=1e-12)


def _lossy_circle(Cl, ba, f0):
    """The loop solved with the wire's resistance in every mode (R_s/(2 pi b) per
    unit length, R_s as sqrt(f)), at k = f/f0 with the geometry fixed."""
    return lambda k: loop_modal.input_impedance(Cl * k, ba * Cl * k / (2 * math.pi), 60,
                                                RS(f0 * k) / (2 * math.pi * ba * Cl * k / (2 * math.pi)))


@pytest.mark.parametrize("f0,Cl,ba", [(300e6, 0.3, 0.01), (3e6, 0.17, 0.004), (300e6, 0.1, 0.05), (300e6, 0.25, 0.001),
                                      (100e6, 0.3, 0.05)])
def test_the_tuned_circular_loop_q_is_its_swept_impedance(registry, f0, Cl, ba):
    """X/R is the Q of an ideal inductor. The real loop's reactance rises faster
    than omega L toward its first resonance, so the Q of the tuned loop - the
    slope of its impedance - is higher, and its band narrower: at a third of
    the way to resonance X/R overstated the VSWR-2 band by 65%. The loss is in
    every mode of the solved loop, as a resistive wire's is: added at the
    terminals instead, it was counted at the feed's current and the Q read 8-12%
    high where the current is not uniform. The thick wire (0.3, 0.05) is in too:
    its terminal impedance depends on the feed, its Q does not."""
    lam = C0 / f0
    b = ba * Cl * lam / (2 * math.pi)
    d = registry["small_circular_loop"].synthesize(f0=f0, C_over_lambda=Cl, b=b, N=1)
    q, bw = _tuned_and_swept(_lossy_circle(Cl, ba, f0), 0.0)
    assert d.metrics["quality_factor"] == pytest.approx(q, rel=6e-3)
    assert d.metrics["fractional_bandwidth_vswr2"] == pytest.approx(bw, rel=6e-3)
    z = d.metrics["input_impedance_ohm"]                     # the terminals, loss included
    if Cl >= 0.25 and not math.isnan(z.real):
        assert z.imag / z.real < 0.8 * q                     # X/R: 0.61 and 0.75 of the tuned Q


@pytest.mark.slow
@pytest.mark.parametrize("f0,P,bs", [(300e6, 0.3, 0.01), (30e6, 0.12, 0.03)])
def test_the_tuned_square_loop_q_is_its_swept_impedance(registry, f0, P, bs):
    """The square's loaded Q against its own method of moments with the wire's loss
    INSIDE it (r times the rooftop functions' Gram matrix), swept at a fixed mesh:
    the loss added at the terminals instead counted it at the feed's current, and the
    Q read 11% high at P = 0.3 wavelengths."""
    lam = C0 / f0
    s = P * lam / 4
    d = registry["small_square_loop"].synthesize(f0=f0, P_over_lambda=P, b=bs * s, N=1)
    n = max(8, min(40, round(1 / (4 * bs))))                # fixed: the sweep scales the loop, not the mesh
    n += n % 2
    q, bw = _tuned_and_swept(lambda k: _mom_square_lossy(P * k, bs * P * k / 4, n, RS(f0 * k)), 0.0)
    assert d.metrics["quality_factor"] == pytest.approx(q, rel=8e-3)
    assert d.metrics["fractional_bandwidth_vswr2"] == pytest.approx(bw, rel=8e-3)


def _mom_square_lossy(P, b, n, rs):
    """_mom_square with a resistive wire: r = R_s / (2 pi b) per wavelength of wire, in
    the Galerkin matrix as r <f_m, f_n> over every segment the two rooftops share."""
    s = P / 4
    c = np.array([[-s / 2, -s / 2, 0], [s / 2, -s / 2, 0], [s / 2, s / 2, 0], [-s / 2, s / 2, 0], [-s / 2, -s / 2, 0]])
    pts = [c[0]] + [c[k] + (c[k + 1] - c[k]) * j / n for k in range(4) for j in range(1, n + 1)]
    m = mom.WireModel([mom.Wire(np.array(pts)[:-1], b, closed=True)])
    Z = mom.impedance_matrix(m)
    halves = [mom._half_weights(m, i) for i in range(m.n_basis)]
    r = rs / (2 * math.pi * b)
    for i in range(m.n_basis):
        for j in range(m.n_basis):
            for p, (a1, b1), _ in halves[i]:
                for q, (a2, b2), _ in halves[j]:
                    if p == q:
                        L = float(m.seg_len[p])
                        Z[i, j] += r * (a1 * a2 * L + (a1 * b2 + b1 * a2) * L ** 2 / 2 + b1 * b2 * L ** 3 / 3)
    V = np.zeros(m.n_basis, dtype=complex)
    V[n // 2 - 1] = 1.0
    return 1.0 / np.linalg.solve(Z, V)[n // 2 - 1]


@pytest.mark.parametrize("f0,Cl,eta", [(300e6, 0.2, 0.5), (30e6, 0.15, 0.7), (3e6, 0.12, 0.25), (100e6, 0.25, 0.8),
                                      (100e6, 0.3, 0.97)])
def test_the_synthesised_wire_meets_the_efficiency_target_on_the_real_loop(registry, f0, Cl, eta):
    """The wire was sized first by the uniform-current law (too fat), then through
    the terminal resistance (too thin: b/a 0.0087 for 0.5 at 0.2 wavelengths, where
    0.0123 is needed). Checked forward on the loop solved with its loss in every
    mode: P_rad / P_in at the wire chosen."""
    d = registry["small_circular_loop"].synthesize(f0=f0, C_over_lambda=Cl, eta_target=eta, N=1)
    b_over_lam = d.get("b") / (C0 / f0)
    solved = loop_modal.tuned_loop(Cl, b_over_lam, RS(f0), 60)["efficiency"]
    assert solved == pytest.approx(eta, abs=0.003)


def test_the_efficiency_counts_the_loss_where_the_current_is(registry):
    """At 0.3 wavelengths the far side carries 1.8 times the feed's current; the
    radiation resistance referred to the loop's mean-square current is 2.2 ohm,
    the terminal resistance 5.8. Efficiency takes the first."""
    a = 0.01 * 0.3 / (2 * math.pi)
    modes = loop_modal.mode_currents(0.3, a, 60)
    far = sum(v * (-1) ** n for n, v in modes.items()) / sum(modes.values())
    assert 1.7 < abs(far) < 1.9
    r_loop, r_term = loop_modal.loop_current_resistance(0.3, a), loop_modal.input_impedance(0.3, a).real
    assert 2.15 < r_loop < 2.25 and 2.3 < r_term / r_loop < 2.8
    lam = C0 / 300e6
    d = registry["small_circular_loop"].synthesize(f0=300e6, C=0.3 * lam, b=a * lam, N=1)
    assert d.metrics["radiation_resistance_loop_current_ohm"] == pytest.approx(r_loop, rel=3e-3)


def test_a_target_past_the_thin_wire_is_nan(registry):
    d = registry["small_circular_loop"].synthesize(f0=3e6, C_over_lambda=0.01, eta_target=0.5, N=1)
    assert math.isnan(d.get("wire_ratio_eta"))
