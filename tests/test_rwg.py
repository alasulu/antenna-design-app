"""The planar RWG solver, and what it says about the three archetypes no solver had examined.

`otahub.num.rwg` is a method of moments on triangles for flat metal in free space.
It is checked against `strip.RooftopStrip` - rooftops on rectangles, a different
basis and different singular integrals - on the same strip with the same finite
gap, where the two meet in the mesh limit to 0.04% (narrow) and 0.2% (wide) in
resistance (reactance within 0.6% and 1.3% of |Z|), and against the wire MoM at the strip's
equivalent radius w/4.
Then it solves the spec bowtie and the two spirals.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from otahub.num import rwg

DATA = json.loads((Path(__file__).parent / "data" / "rwg_planar.json").read_text())
ETA0 = 376.730313668


def _extrap(a, b):
    """First order in the cell size from meshes of 60 and 80 cells: x80 + 3 (x80 - x60)."""
    return b + 3 * (b - a)


def test_the_closed_form_potentials():
    """Against a fine midpoint rule, at a point inside the triangle and two outside."""
    v = np.array([[[0.0, 0.0], [1.0, 0.1], [0.3, 0.8]]])
    pts = np.array([[[0.4, 0.3], [0.9, 0.9], [-0.3, 0.2]]])
    S0, S1 = rwg.tri_potentials(v, pts)
    n = 1200
    u, w = np.meshgrid((np.arange(n) + 0.5) / n, (np.arange(n) + 0.5) / n, indexing="ij")
    m = u + w < 1
    q = v[0, 0] + np.outer(u[m], v[0, 1] - v[0, 0]) + np.outer(w[m], v[0, 2] - v[0, 0])
    e1, e2 = v[0, 1] - v[0, 0], v[0, 2] - v[0, 0]
    dA = abs(e1[0] * e2[1] - e1[1] * e2[0]) / n ** 2            # each sample carries 2A/n^2
    for k, r in enumerate(pts[0]):
        d = q - r
        R = np.hypot(d[:, 0], d[:, 1])
        assert S0[0, k] == pytest.approx(np.sum(1 / R) * dA, rel=3e-3)
        assert S1[0, k] == pytest.approx((d / R[:, None]).sum(0) * dA, rel=3e-3, abs=1e-3)


def test_the_streamed_matrix_is_the_assembled_one():
    m = rwg.strip_mesh(0.3, 0.02, 12, 2)
    C = m.spread()
    M = rwg.corner_matrix(m)
    Z = rwg.impedance_matrix(m)
    assert np.max(np.abs(Z - np.asarray(C.T @ (C.T @ M.T).T))) < 1e-12 * np.max(np.abs(Z))


def test_power_balances():
    m = rwg.strip_mesh(0.47, 0.02, 40, 4)
    s = rwg.solve(m, rwg.gap_vector(m, 0.0, 0.47 / 20))
    assert s.circuit_power == pytest.approx(rwg.radiated_power(s), rel=1e-6)


@pytest.mark.parametrize("W", [0.02, 0.06])
def test_the_rooftop_strip_agrees_in_the_mesh_limit(W):
    rows = [r for r in DATA["strip"] if r["W"] == W and "nx" in r]
    by = {r["nx"]: r for r in rows}
    for part, tol in ((0, 0.003), (1, None)):
        rt = _extrap(by[60]["rooftop"][part], by[80]["rooftop"][part])
        rw = _extrap(by[60]["rwg"][part], by[80]["rwg"][part])
        if tol is not None:
            assert rw == pytest.approx(rt, rel=tol)
        else:
            # the reactance is not monotone in the mesh on the wide strip, so its extrapolation is
            # rougher: within 1.5% of |Z| (0.6% on the narrow strip)
            assert abs(rw - rt) < 0.015 * abs(complex(*by[80]["rwg"]))
    changes = [abs(by[80][k][0] - by[40][k][0]) for k in ("rooftop", "rwg")]
    assert changes[1] < changes[0]                              # the triangles converge faster


def test_one_strip_live():
    r = next(x for x in DATA["strip"] if x.get("W") == 0.02 and x.get("nx") == 40)
    m = rwg.strip_mesh(0.47, 0.02, 40, 4)
    z = rwg.solve(m, rwg.gap_vector(m, 0.0, 0.47 / 20)).input_impedance
    assert z.real == pytest.approx(r["rwg"][0], rel=1e-9) and z.imag == pytest.approx(r["rwg"][1], rel=1e-9)


def test_the_equivalent_radius_holds_on_a_narrow_strip_not_a_wide_one():
    for W, tol in ((0.02, 0.015), (0.06, None)):
        by = {r["nx"]: r for r in DATA["strip"] if r["W"] == W and "nx" in r}
        strip_R = _extrap(by[60]["rwg"][0], by[80]["rwg"][0])
        wire_R = next(r for r in DATA["strip"] if r["W"] == W and "wire_equivalent" in r)["wire_equivalent"][0]
        if tol:
            assert wire_R == pytest.approx(strip_R, rel=tol)
        else:
            assert wire_R / strip_R - 1 > 0.03                  # w/4 overstates a wide strip, as strip.py found


def test_the_bowtie_is_converged():
    coarse = next(r for r in DATA["bowtie_sweep"] if r["fr"] == 4.0)
    for fine in DATA["bowtie_fine_4"]:
        assert fine["Z"][0] == pytest.approx(coarse["Z"][0], rel=0.01)
        assert abs(fine["Z"][1] - coarse["Z"][1]) < 5
    fl = {}
    for r in DATA["bowtie_flare_f_low"]:
        fl.setdefault(r["flare"], []).append(r)
    for flare, (a, b) in fl.items():
        assert b["Z"][0] == pytest.approx(a["Z"][0], rel=0.015) and b["D"] == pytest.approx(a["D"], rel=0.01)


def test_a_large_self_complementary_bowtie_approaches_mushiake():
    """Mushiake's eta0/2 is for an infinite sheet; the finite 90-degree bowtie, once its
    wings are several half-wavelengths long, ripples about it."""
    high = [r["Z"][0] for r in DATA["bowtie_sweep"] if r["fr"] >= 2.5]
    assert np.mean(high) == pytest.approx(ETA0 / 2, rel=0.12)
    assert all(abs(r["Z"][1]) < 45 for r in DATA["bowtie_sweep"] if r["fr"] >= 2.5)


def test_the_bowtie_beam_splits_above_twice_f_low():
    sweep = {r["fr"]: r for r in DATA["bowtie_sweep"]}
    assert all(1.5 < sweep[f]["D_broadside"] < 2.3 for f in sweep if f <= 2.0)
    assert sweep[3.0]["D_broadside"] < 0.5 and sweep[3.5]["D_broadside"] < 0.05
    assert all(sweep[f]["D_max"] > 2.0 for f in sweep if f >= 2.5)


def test_power_balances_on_every_solved_antenna():
    for key in ("bowtie_sweep", "equi_sweep", "arch_sweep"):
        for r in DATA[key]:
            assert abs(r["power_mismatch"]) < 1e-3, (key, r["fr"])        # 7e-4 on the coarsest spiral mesh


def _sweep(kind):
    return {r["fr"]: r for r in DATA[f"{kind}_sweep"]}


@pytest.mark.parametrize("kind", ["equi", "arch"])
def test_the_spirals_are_converged(kind):
    rows = [r for r in DATA["spiral_convergence_2"] if r["kind"] == kind]
    Rs, ARs, Ds = [r["Z"][0] for r in rows], [r["AR_db"] for r in rows], [r["D"] for r in rows]
    assert max(Rs) / min(Rs) - 1 < 0.035
    assert max(ARs) - min(ARs) < 0.1 and max(Ds) / min(Ds) - 1 < 0.005


@pytest.mark.parametrize("kind", ["equi", "arch"])
def test_at_f_low_a_spiral_is_not_yet_circularly_polarised(kind):
    """The specs put f_low where the outer circumference is one wavelength. There the
    axial ratio is 16-21 dB; it falls through 3 dB only well above."""
    s = _sweep(kind)
    assert s[1.0]["AR_db"] > 15
    frs = sorted(s)
    assert all(s[a]["AR_db"] >= s[b]["AR_db"] - 0.2 for a, b in zip(frs, frs[1:]))     # falls with frequency
    assert all(s[f]["AR_db"] < 1.0 for f in frs if f >= 3.0)
    edge = DATA[f"{kind}_cp_edge"]["edge"]
    assert (2.0 < edge < 2.5) if kind == "equi" else (1.25 < edge < 1.5)


@pytest.mark.parametrize("kind", ["equi", "arch"])
def test_spiral_directivity_is_well_above_the_old_figure(kind):
    """The specs said about 1.5 (1.8 dBi) per side; solved, the broadside directivity
    runs 2.7-3.5 dBi at f_low rising to about 6 dBi at four times it."""
    s = _sweep(kind)
    assert all(s[f]["D_broadside"] > 1.8 for f in s)
    assert s[3.0]["D_broadside"] > 3.2 and s[4.0]["D_broadside"] > 3.9
