"""A finite feed gap for the wire method of moments.

The delta gap - 1 V across one node - has a capacitance that grows as the mesh
refines, so on fat wire the driving point never converges: the resonant
resistance of a 0.005-wavelength dipole climbs about 1% per mesh doubling. The
finite gap is the model `bor` uses for tubes: V/g impressed along a gap of
physical length g, the input current averaged over it. It converges, and it
meets the body-of-revolution tube solver - a different formulation, the current
on the tube's surface, rims resolved - across the tube table of resonances.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from otahub.num import mom

HERE = Path(__file__).parent / "data"
TABLE = json.loads((HERE / "wire_gap_bor.json").read_text())


def _res(f, L0=0.46, L1=0.44):
    x0, x1 = f(L0).imag, f(L1).imag
    for _ in range(40):
        L2 = L1 - x1 * (L1 - L0) / (x1 - x0)
        L0, x0, L1 = L1, x1, L2
        z = f(L1)
        x1 = z.imag
        if abs(L1 - L0) < 1e-7:
            return L1, z.real
    raise RuntimeError("no resonance")


def _gap(L, a, g, seg):
    m, gs = mom.gapped_dipole(L, a, g, seg)
    return mom.solve(m, exact=True, gap=gs).input_impedance


def test_the_gap_weights():
    m = mom.dipole(0.5, 1e-3, 40)
    w = mom.gap_vector(m, (19, 20))
    assert np.flatnonzero(w).tolist() == [18, 19, 20]
    assert w[[18, 19, 20]] == pytest.approx([0.25, 0.5, 0.25]) and w.sum() == pytest.approx(1.0)


def test_power_balances():
    """The gap-averaged current makes V conj(I)/2 the power delivered. On thin wire it
    balances the far field to 1e-7; on fat wire both feeds share the thin-wire model's
    own (ka)^2 mismatch - the kernel couples surface currents, the far field radiates a
    line current - and the finite gap adds nothing to it."""
    def mismatch(a, finite):
        if finite:
            m, g = mom.gapped_dipole(0.47, a, 0.01, 0.005)
            s = mom.solve(m, exact=True, gap=g)
        else:
            s = mom.solve(mom.dipole(0.47, a, 94), exact=True)
        return s.circuit_power / mom.radiated_power(s, 90, 90) - 1
    assert abs(mismatch(1e-4, True)) < 1e-6
    assert mismatch(3e-3, True) == pytest.approx(mismatch(3e-3, False), abs=1e-6)


def test_on_thin_wire_a_small_gap_is_the_delta_gap():
    zd = mom.input_impedance(mom.dipole(0.5, 1e-4, 80))
    m, g = mom.gapped_dipole(0.5, 1e-4, 0.0125, 0.00625)
    zg = mom.solve(m, gap=g).input_impedance
    assert zg.real == pytest.approx(zd.real, rel=0.01) and zg.imag == pytest.approx(zd.imag, rel=0.01)


def test_a_finite_gap_converges_where_the_delta_gap_drifts():
    a = 0.005
    R_gap = [_res(lambda L: _gap(L, a, 2 * a, seg))[1] for seg in (0.01, 0.005)]
    assert abs(R_gap[1] - R_gap[0]) / R_gap[0] < 0.002

    def delta(L, n):
        z = np.linspace(-L / 2, L / 2, n + 1)
        return mom.solve(mom.WireModel([mom.Wire(np.stack([0 * z, 0 * z, z], 1), a)]), n // 2 - 1,
                         exact=True).input_impedance
    R_delta = [_res(lambda L: delta(L, n))[1] for n in (50, 100)]
    assert (R_delta[1] - R_delta[0]) / R_delta[0] > 0.008


def test_the_tube_solver_agrees_across_its_table():
    """Resonant length and resistance of the finite-gap wire, extrapolated in the
    mesh, against the body-of-revolution tube with its rims resolved."""
    for r in TABLE["rows"]:
        if r["bor_L"] is None or r["a"] > TABLE["agree_to_a"]:
            continue
        assert r["L_extrap"] == pytest.approx(r["bor_L"], rel=TABLE["L_tol"]), r
        assert r["R_extrap"] == pytest.approx(r["bor_R"], rel=TABLE["R_tol"]), r


def test_one_row_live():
    r = next(x for x in TABLE["rows"] if x["a"] == 0.005 and x["gap_over_a"] == 2.0)
    L, R = _res(lambda L: _gap(L, 0.005, 0.01, 0.0025))
    assert L == pytest.approx(r["L_0.0025"], abs=1e-6) and R == pytest.approx(r["R_0.0025"], abs=1e-4)
