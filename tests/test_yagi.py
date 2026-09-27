"""The Yagi-Uda's gain, solved rather than read off a table in the wrong unit.

The spec fitted NBS Technical Note 688's six optimised designs, whose gains are
tabulated over a half-wave dipole, and called the result dBi. Built and solved
in free space - by the thin-wire method of moments and, independently, by
coupled Hallen equations carried here - the designs sit about 2 dB above that
fit, and NBS's own figures read as dBd sit within a third of a decibel of them.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from otahub.core.constants import ETA0
from otahub.num import mom

DATA = json.loads((Path(__file__).parent / "data" / "yagi_nbs.json").read_text())
DESIGNS = {d["boom"]: d for d in DATA["designs"]}
SOLVED = {s["boom"]: s for s in DATA["solved"]}
A = 0.0085 / 2
K = 2 * math.pi


def _geometry(boom, driven=0.47):
    d = DESIGNS[boom]
    xs = [0.0, 0.2] + [0.2 + d["spacing"] * (i + 1) for i in range(len(d["directors"]))]
    return xs, [d["reflector"], driven] + d["directors"]


def _mom(boom, seg):
    xs, Ls = _geometry(boom)
    wires, feed, base = [], None, 0
    for k, (x, L) in enumerate(zip(xs, Ls)):
        n = max(8, int(math.ceil(L / seg)))
        n += n % 2
        z = np.linspace(-L / 2, L / 2, n + 1)
        wires.append(mom.Wire(np.stack([np.full_like(z, x), np.zeros_like(z), z], axis=1), A))
        if k == 1:
            feed = base + n // 2 - 1
        base += n - 1
    sol = mom.solve(mom.WireModel(wires), feed, exact=True)
    return 10 * math.log10(mom.directivity_towards(sol, math.pi / 2, 0.0, 120, 120))


def _hallen(boom, n=40, nq=24):
    """Coupled Hallen equations for parallel z-directed wires, point matching:
    sum_j int I_j G = -(j/eta)[B_i cos kz + (V_i/2) sin k|z|] on every wire."""
    xs, Ls = _geometry(boom)
    W = len(xs)
    zs = [np.linspace(-L / 2, L / 2, n + 1) for L in Ls]
    zms = [np.concatenate([z[1:n], [z[-1]]]) for z in zs]
    nb = n - 1
    NU = W * nb + W
    M = np.zeros((NU, NU), dtype=complex)
    rhs = np.zeros(NU, dtype=complex)
    x, w = np.polynomial.legendre.leggauss(nq)
    for i in range(W):
        rows = slice(i * n, (i + 1) * n)
        for j in range(W):
            z = zs[j]
            dz = z[1] - z[0]
            rho2 = A * A if i == j else (xs[i] - xs[j]) ** 2
            for b in range(nb):
                c = z[b + 1]
                acc = 0
                for lo, hi, rise in ((c - dz, c, True), (c, c + dz, False)):
                    s = 0.5 * (hi - lo) * (x + 1.0) + lo
                    f = (s - (c - dz)) / dz if rise else ((c + dz) - s) / dz
                    R = np.sqrt((zms[i][:, None] - s[None, :]) ** 2 + rho2)
                    acc = acc + 0.5 * (hi - lo) * ((np.exp(-1j * K * R) / (4 * math.pi * R)) * f * w).sum(axis=1)
                M[rows, j * nb + b] = acc
        M[rows, W * nb + i] = (1j / ETA0) * np.cos(K * zms[i])
        if i == 1:
            rhs[rows] = -(1j / ETA0) * 0.5 * np.sin(K * np.abs(zms[i]))
    sol = np.linalg.solve(M, rhs)
    cur = [np.concatenate([[0.0], sol[j * nb:(j + 1) * nb], [0.0]]) for j in range(W)]

    def far(th, ph):
        out = 0
        for x0, z, I in zip(xs, zs, cur):
            e = np.exp(1j * K * (np.multiply.outer(np.cos(th), z) + x0 * (np.sin(th) * np.cos(ph))[..., None]))
            out = out + np.trapezoid(I * e, z, axis=-1)
        return np.sin(th) * out

    th = np.linspace(0, math.pi, 121)
    ph = 2 * math.pi * np.arange(240) / 240
    T, P = np.meshgrid(th, ph, indexing="ij")
    U = np.abs(far(T, P)) ** 2
    total = np.trapezoid((U * np.sin(T)).sum(axis=1) * 2 * math.pi / 240, th)
    return 10 * math.log10(4 * math.pi * abs(far(np.array(math.pi / 2), np.array(0.0))) ** 2 / total)


@pytest.mark.parametrize("boom", [0.4, 0.8, 1.2])
def test_the_recorded_solutions_reproduce_live(boom):
    s = SOLVED[boom]
    assert _mom(boom, 0.02) == pytest.approx(s["mom_002"], abs=1e-6)
    assert _hallen(boom) == pytest.approx(s["hallen_40"], abs=0.02)


def test_the_two_solvers_agree():
    for s in DATA["solved"]:
        assert s["hallen_80"] == pytest.approx(s["mom_extrap"], abs=0.1)


def test_nbs_tabulated_gain_over_a_dipole():
    """Read as dBi, NBS's figures sit 1.8-2.2 dB under the solved designs; read
    as dB over a half-wave dipole they sit within a third of a decibel."""
    for s in DATA["solved"]:
        q = DESIGNS[s["boom"]]["quoted_dbd"]
        assert 1.8 < s["mom_extrap"] - q < 2.2
        assert -0.1 < q + 2.15 - s["mom_extrap"] < 0.35


@pytest.mark.parametrize("boom", [0.4, 0.8, 1.2, 2.2, 3.2, 4.2])
def test_the_spec_gain_and_beams_are_the_solved_designs(registry, boom):
    d = registry["yagi_uda"].synthesize(f0=3e8, boom_over_lambda=boom)
    s = SOLVED[boom]
    assert d.metrics["gain_dbi"] == pytest.approx(s["mom_extrap"], abs=0.09)
    assert d.metrics["gain_dbd"] == pytest.approx(s["mom_extrap"] - 2.15, abs=0.09)
    assert d.metrics["hpbw_e_deg"] == pytest.approx(s["hpbw_e"], rel=0.02)
    assert d.metrics["hpbw_h_deg"] == pytest.approx(s["hpbw_h"], rel=0.03)
    assert d.metrics["director_count_nbs"] == len(DESIGNS[boom]["directors"])
    assert d.metrics["front_to_back_db"] > 0


def test_the_old_beamwidth_was_far_too_wide_on_short_booms():
    s = SOLVED[0.4]
    assert 55.0 / math.sqrt(0.4) / s["hpbw_e"] - 1 > 0.4
