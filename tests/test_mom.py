"""Validation of the thin-wire method of moments.

A solver is only worth what it has been checked against, so this file is the
solver's licence. Six independent checks, in increasing order of how much of
the code they exercise:

1. the impedance matrix alone, driven with a prescribed sinusoid, reproduces
   the induced-EMF 73.0796 + j42.5152 ohm - no solve involved;
2. the far-field code, given that same prescribed current, returns the same
   73.08 ohm and D = 1.64093 - no matrix involved;
3. the pattern integral and the circuit power 0.5*Re(V I*) agree, which ties
   the solved current to the radiated field;
4. limits that are exact: D = 1.5 for a short dipole and for a small loop, and
   Rr = 20*pi^2*(C/lambda)^4 for a small loop;
5. a thin dipole resonates just below half a wavelength at about 72 ohm;
6. a folded dipole shows the 4:1 transformation, which also exercises a closed
   multi-conductor path with corners.

A seventh test records something that cost real time to establish: the classic
73.08 + j42.52 is the induced-EMF value for an ASSUMED sinusoidal current, and
it is NOT the driving-point impedance of a delta-gap-fed wire of finite radius.
The two differ by 18% at a/lambda = 0.001. An independent Hallen solve - a
different integral equation, no divergence term, point matching instead of
Galerkin - agrees with the EFIE, not with the textbook constant.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from otahub.num import mom
from otahub.num.hallen import hallen_dipole as _hallen_dipole

EMF_R, EMF_X = 73.0796, 42.5152      # induced-EMF half-wave dipole [ohm]


def _sinusoid(model, length):
    z = np.array([model.node_of(n)[2] for n in range(model.n_basis)])
    return np.sin(mom.K * (0.5 * length - np.abs(z))).astype(complex)


# ------------------------------------------------------------------ 1. matrix

def test_matrix_reproduces_the_induced_emf_impedance():
    """I^T Z I with the sinusoidal current is the induced-EMF integral, which
    has a published closed form. This touches the matrix and nothing else."""
    m = mom.dipole(0.5, 1e-5, 120)
    I = _sinusoid(m, 0.5)
    z = complex((I @ mom.impedance_matrix(m) @ I) / I.max() ** 2)
    assert z.real == pytest.approx(EMF_R, rel=2e-3)
    assert z.imag == pytest.approx(EMF_X, rel=2e-2)


def test_matrix_is_symmetric():
    """Reciprocity. Galerkin testing makes it exact, so any asymmetry is a bug
    in the assembly rather than a tolerance question."""
    Z = mom.impedance_matrix(mom.loop(0.5, 1e-3, 24))
    assert np.abs(Z - Z.T).max() == 0.0


# ------------------------------------------------------------- 2. far field

def test_far_field_of_a_prescribed_sinusoid():
    """Same reference number by a completely different route: no matrix, just
    the radiation integral of a current that was handed to it."""
    m = mom.dipole(0.5, 1e-3, 80)
    I = _sinusoid(m, 0.5)
    sol = mom.MoMSolution(m, I, int(np.argmax(np.abs(I))))
    rr = 2.0 * mom.radiated_power(sol, 100, 100) / abs(I).max() ** 2
    assert rr == pytest.approx(EMF_R, rel=2e-3)
    assert mom.directivity(sol, 100, 100) == pytest.approx(1.64093, rel=2e-3)


# --------------------------------------------------------- 3. power balance

@pytest.mark.parametrize("length", [0.1, 0.5, 1.0])
def test_pattern_power_matches_circuit_power(length):
    """0.5*Re(V I*) at the terminals against the integral of the far field over
    the sphere. They share no algebra, so this ties the whole chain together."""
    sol = mom.solve(mom.dipole(length, 1e-3, 60))
    assert mom.radiated_power(sol, 70, 70) == pytest.approx(
        sol.circuit_power, rel=1e-4)


# ----------------------------------------------------------------- 4. limits

@pytest.mark.parametrize("length", [0.02, 0.1])
def test_short_dipole_directivity_is_three_halves(length):
    sol = mom.solve(mom.dipole(length, 1e-4, 30))
    assert mom.directivity(sol, 70, 70) == pytest.approx(1.5, rel=5e-3)


def test_small_loop_directivity_is_also_three_halves():
    """A small loop is a magnetic dipole; the pattern is the same doughnut."""
    sol = mom.solve(mom.loop(0.05, 1e-4, 36))
    assert mom.directivity(sol, 70, 70) == pytest.approx(1.5, rel=5e-3)


def test_small_loop_radiation_resistance_tends_to_the_closed_form():
    """Rr = 20*pi^2*(C/lambda)^4 is exact for a UNIFORM current, so the ratio
    must approach 1 as the loop shrinks and the current flattens - and must not
    be assumed at sizes where it does not. At C = 0.1 lambda the delta-gap
    current already varies by 6% and Rr is 11% above the uniform value."""
    ratios = {}
    for c in (0.1, 0.05, 0.025):
        sol = mom.solve(mom.loop(c, c * 1e-3, 36))
        rr = 2.0 * sol.circuit_power / abs(sol.feed_current) ** 2
        ratios[c] = rr / (20.0 * math.pi ** 2 * c ** 4)
    assert ratios[0.025] == pytest.approx(1.0, abs=0.02)
    assert ratios[0.05] < ratios[0.1]
    assert ratios[0.1] > 1.05


# -------------------------------------------------------------- 5. resonance

def test_thin_dipole_resonates_just_below_half_a_wavelength():
    """The textbook shortening: a wire dipole is resonant near 0.475 lambda at
    a/lambda = 0.001, with a resistance a little under the 73 ohm figure."""
    lo, hi = 0.42, 0.52
    for _ in range(22):
        mid = 0.5 * (lo + hi)
        if mom.input_impedance(mom.dipole(mid, 1e-3, 40)).imag < 0:
            lo = mid
        else:
            hi = mid
    z = mom.input_impedance(mom.dipole(lo, 1e-3, 40))
    assert 0.465 < lo < 0.485, f"resonant length {lo}"
    assert 65.0 < z.real < 80.0, f"resonant resistance {z.real}"


# ----------------------------------------------------------- 6. folded dipole

def test_folded_dipole_steps_the_impedance_up_four_times():
    """The classic 4:1, and the only geometry here with corners and two
    parallel conductors, so it exercises the non-collinear path as well."""
    length = 0.46
    m = mom.folded_dipole_wire(length, 0.01, 1e-3, 80)
    nodes = np.array([m.node_of(n) for n in range(m.n_basis)])
    feed = int(np.argmin(np.abs(nodes[:, 0]) + np.abs(nodes[:, 1] - 0.005)))
    folded = mom.solve(m, feed).input_impedance
    plain = mom.input_impedance(mom.dipole(length, 1e-3, 80))
    assert abs(folded / plain) == pytest.approx(4.0, rel=0.05)


# ------------------------------------------ 7. what 73 ohm actually refers to

def test_delta_gap_impedance_is_not_the_induced_emf_value():
    """The trap this solver walked into, recorded so nobody walks into it twice.

    73.08 + j42.52 is the induced-EMF result for an ASSUMED sinusoidal current
    on a vanishingly thin wire. The driving-point impedance of a delta-gap-fed
    wire of a/lambda = 0.001 is about 86 ohm - nearly 20% higher - because the
    real current is visibly fatter than a sinusoid near the ends. An entirely
    independent Hallen solve agrees with the EFIE, not with the constant.
    """
    efie = mom.input_impedance(mom.dipole(0.5, 1e-3, 80))
    hallen = _hallen_dipole(0.5, 1e-3, 120)
    assert efie.real == pytest.approx(hallen.real, rel=0.03)
    assert efie.imag == pytest.approx(hallen.imag, rel=0.05)
    assert efie.real > 1.15 * EMF_R, (
        "the delta-gap value should sit well above the induced-EMF one")


def test_feed_lands_on_a_node_whatever_segment_count_is_asked_for():
    """An odd segment count has no node at the centre, so the feed would sit
    half a segment off and quietly break the symmetry. dipole() rounds up."""
    for asked in (39, 40, 41):
        m = mom.dipole(0.5, 1e-3, asked)
        sol = mom.solve(m)
        assert m.node_of(sol.feed)[2] == pytest.approx(0.0, abs=1e-12)


def test_segments_shorter_than_a_few_radii_leave_the_approximation():
    """A fine mesh on a fat wire is not more accurate - it is outside the
    reduced kernel. Resonant resistance on a 0.006-wavelength wire is stable
    while segments are several radii long and runs away below one radius."""
    def res_r(n):
        lo, hi = 0.40, 0.53
        for _ in range(26):
            mid = 0.5 * (lo + hi)
            if mom.input_impedance(mom.dipole(mid, 0.006, n)).imag < 0:
                lo = mid
            else:
                hi = mid
        return mom.input_impedance(mom.dipole(lo, 0.006, n)).real
    ok = [res_r(n) for n in (10, 16)]              # 7.8 and 4.9 radii
    assert ok[1] == pytest.approx(ok[0], rel=0.03)
    with pytest.warns(UserWarning, match="shorter than their wire radius"):
        broken = res_r(100)                        # 0.8 radii
    assert broken > 1.3 * ok[1]


# ------------------------------------------------------------ exact kernel

def test_the_exact_kernel_is_off_by_default():
    """Every result derived with the reduced kernel must reproduce bit for bit."""
    m = mom.dipole(0.47, 1e-3, 40)
    assert mom.input_impedance(m) == mom.input_impedance(m, exact=False)
    assert mom.input_impedance(m) == complex(69.06656211407594, -10.49855834396867)


def test_the_exact_kernel_changes_little_on_thin_wire():
    m = mom.dipole(0.47, 1e-3, 41)
    reduced, exact = mom.input_impedance(m), mom.input_impedance(m, exact=True)
    assert abs(exact - reduced) / abs(reduced) < 0.005


def test_the_exact_kernel_holds_a_fat_wire_the_reduced_kernel_loses():
    """The same 0.006-wavelength wire, meshed from 3.7 radii per segment down to
    0.8. With the exact kernel the resonance stays put and the resistance creeps
    about 3% (the delta gap, not the kernel); with the reduced kernel the finest
    mesh knocks the reactance 16 ohm off resonance."""
    exact = []
    for n in (21, 41, 61, 97):
        if n < 97:
            m = mom.dipole(0.465, 0.006, n)
        else:
            with pytest.warns(UserWarning, match="exact=True"):
                m = mom.dipole(0.465, 0.006, n)
        exact.append(mom.input_impedance(m, exact=True))
    assert all(abs(z.imag) < 1.5 for z in exact), exact
    rs = [z.real for z in exact]
    assert max(rs) / min(rs) < 1.045, rs
    assert mom.input_impedance(m).imag < -10.0


# ------------------------------------------------ found by review, fixed

def test_wires_of_different_radii_give_the_same_answer_in_either_order():
    """The reduced kernel took the source segment's radius and only the upper
    triangle was assembled, so two parallel wires of 0.003 and 0.001 lambda
    changed resistance by 2.4% when listed the other way round. The pair now
    uses (a_p^2 + a_q^2)/2, and equal radii are untouched to the bit."""
    from otahub.num import mom

    def model(order, exact):
        w1 = mom.Wire(np.linspace((0, 0, -0.235), (0, 0, 0.235), 21), 0.003)
        w2 = mom.Wire(np.linspace((0.01, 0, -0.235), (0.01, 0, 0.235), 21), 0.001)
        m = mom.WireModel([w1, w2] if order == 0 else [w2, w1])
        feed = int(np.argmin([np.linalg.norm(m.node_of(n)) for n in range(m.n_basis)]))
        return mom.input_impedance(m, feed, exact=exact)

    for exact in (False, True):
        a, b = model(0, exact), model(1, exact)
        assert abs(a - b) < 1e-12 * abs(a), exact


def test_the_far_field_is_exact_beside_broadside():
    """The closed-form segment integral cancelled catastrophically a hair off
    broadside (60% wrong at 1e-8 rad). A lone rooftop on nodes -0.1, 0, 0.2
    against its current integrated by quadrature, at every offset."""
    from scipy.integrate import quad
    from otahub.num import mom
    m = mom.WireModel([mom.Wire(np.array([(0, 0, -0.1), (0, 0, 0), (0, 0, 0.2)]), 1e-4)])
    cur = np.zeros(m.n_basis, complex)
    cur[0] = 1.0
    sol = mom.MoMSolution(m, cur, 0)
    tri = lambda z: 1 + z / 0.1 if z < 0 else 1 - z / 0.2
    for d in (0.0, 1e-8, 1e-7, 1e-4, 0.05, 0.3):
        th = math.pi / 2 - d
        e_th, _ = mom.far_field(sol, np.array([th]), np.array([0.0]))
        ct = math.cos(th)
        re = quad(lambda z: tri(z) * math.cos(2 * math.pi * z * ct), -0.1, 0.2, points=[0], epsabs=0, epsrel=1e-13)[0]
        im = quad(lambda z: tri(z) * math.sin(2 * math.pi * z * ct), -0.1, 0.2, points=[0], epsabs=0, epsrel=1e-13)[0]
        assert abs(e_th[0]) == pytest.approx(abs(complex(re, im)) * math.sin(th), rel=1e-12), d


def test_a_downward_reactance_crossing_is_found_not_its_bracket_end():
    """Between 0.85 and 1.15 wavelengths a dipole's reactance falls through zero
    (the antiresonance); the bisection assumed it rose and returned 1.15."""
    from scipy.optimize import brentq
    from otahub.num import mom
    z = lambda s: mom.solve(mom.dipole(s, .001, 40)).input_impedance
    root = mom.resonant_scale(z)
    assert root == pytest.approx(brentq(lambda s: z(s).imag, 0.85, 1.15, xtol=1e-10), abs=1e-7)



def test_collinear_wires_of_different_radii_are_order_independent_with_the_exact_kernel():
    """The exact kernel's cutoff used the source radius alone; two collinear wires
    of 0.003 and 0.001 lambda gave 12.028 - j325.589 or 12.023 - j325.608 ohm
    with the list reversed."""
    from otahub.num import mom

    def z(order):
        w1 = mom.Wire(np.linspace((0, 0, -0.235), (0, 0, 0.0), 11), 0.003)
        w2 = mom.Wire(np.linspace((0, 0, 0.008), (0, 0, 0.243), 11), 0.001)
        m = mom.WireModel([w1, w2] if order == 0 else [w2, w1])
        feed = int(np.argmin([np.linalg.norm(m.node_of(n) - np.array([0, 0, -0.1175]))
                              for n in range(m.n_basis)]))
        return mom.input_impedance(m, feed, exact=True)

    assert abs(z(0) - z(1)) < 1e-9 * abs(z(0))


@pytest.mark.parametrize("zero_at", [0.85, 1.15])
def test_a_resonance_on_the_bracket_end_is_found(zero_at):
    from otahub.num import mom
    for sign in (1, -1):
        assert mom.resonant_scale(lambda s: 50 + 1j * sign * (s - zero_at)) == zero_at
