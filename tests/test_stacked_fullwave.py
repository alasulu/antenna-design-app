"""The stacked patch, full-wave: the two-layer spectral MoM against the FDTD.

The FDTD rings the stack down (`patch_fdtd.stacked_ringdown`) and feeds it through a
probe (`patch_fdtd.probe_impedance`); the spectral MoM (`patch_sdm.StackedPatch`,
`probe_vector`) solves the same geometry. They share nothing but the geometry. The
runs are in tests/data/stacked_patch_fullwave.json; the SDM values there are checked
live at a spot, the FDTD's coarsest ringdown is reproduced, and the probe port is
checked live on a monopole against the wire MoM.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import least_squares

from otahub.num import patch_fdtd as pf
from otahub.num import patch_sdm as sdm

DATA = json.loads((Path(__file__).parent / "data" / "stacked_patch_fullwave.json").read_text())
LAM = 2.99792458e8 / 1e9
STACK = (2.2, 0.0128 * LAM, 0.32 * LAM, 0.3968 * LAM, 1.0, 0.032 * LAM, 0.3072 * LAM, 0.3712 * LAM)


def _pole(r, z, f_guess, q_guess):
    """One conjugate pole pair plus a quadratic background, fitted to Z over r (f/f0)."""
    s = 1j * np.asarray(r)

    def model(x):
        p, res = complex(x[0], x[1]), complex(x[2], x[3])
        back = sum(complex(x[4 + 2 * k], x[5 + 2 * k]) * (s - 1j * f_guess) ** k for k in range(3))
        return res / (s - p) + np.conj(res) / (s - np.conj(p)) + back

    x0 = [-f_guess / (2 * q_guess), f_guess, 0.1, 0.0] + [0.0] * 6
    x = least_squares(lambda x: np.concatenate([(model(x) - z).real, (model(x) - z).imag]), x0, max_nfev=20000).x
    return x[1], x[1] / (-2 * x[0])


def _series(rows, key):
    r = np.array([d["r"] for d in rows])
    return r, np.array([complex(*d[key]) for d in rows])


def _peak(r, R):
    k = int(np.argmax(R))
    a, b, c = R[k - 1], R[k], R[k + 1]
    d = 0.5 * (a - c) / (a - 2 * b + c)
    return r[k] + d * (r[1] - r[0]), b - 0.25 * (a - c) * d


def _extrapolated(mode):
    rows = DATA["ringdowns"]
    return pf.extrapolate([d["nh"] for d in rows], [d["modes"][mode] for d in rows])


def test_the_ringdowns_converge_at_first_order():
    rows = DATA["ringdowns"][:3]
    for mode in (0, 1):
        f = [d["modes"][mode] for d in rows]
        assert (f[1] - f[0]) / (f[2] - f[1]) == pytest.approx(3.0, rel=0.1)   # 1/nh: (1/2 - 1/4)/(1/4 - 1/6)


def test_the_lower_mode_is_the_sdms():
    r, z = _series(DATA["sdm"]["stack"], "zA")
    m = (r >= 0.94) & (r <= 1.06)
    f_sdm, q_sdm = _pole(r[m], z[m], 0.993, 34)
    assert _extrapolated(0) == pytest.approx(f_sdm, rel=5e-4)
    assert DATA["ringdowns"][-1]["q"][0] == pytest.approx(q_sdm, rel=0.05)


def test_the_upper_mode_is_the_sdms_pole_not_its_characteristic_zero():
    """For a mode of Q 13 the characteristic eigenvalue crosses zero 1.2% above the
    natural frequency the ringdown measures; the SDM's own pole agrees with the FDTD."""
    r, z = _series(DATA["sdm"]["stack"], "zA")
    m = (r >= 1.28) & (r <= 1.48)
    f_sdm, _ = _pole(r[m], z[m], 1.375, 13)
    assert _extrapolated(1) == pytest.approx(f_sdm, rel=4e-3)
    # the characteristic zero, 1.3932 (test_the_characteristic_zeros_of_the_stack), sits above
    assert 1.3932 / f_sdm > 1.01


@pytest.mark.slow
def test_the_characteristic_zeros_of_the_stack():
    p = sdm.StackedPatch(*STACK, **sdm.BASIS_FULL)
    f, v = sdm.stacked_resonance(p, 0.99e9)
    assert f / 1e9 == pytest.approx(0.99356, rel=2e-4)
    assert np.linalg.norm(v[:p.n1]) / np.linalg.norm(v) > 0.99            # the driven patch's mode
    f, v = sdm.stacked_resonance(p, 1.37e9)
    assert f / 1e9 == pytest.approx(1.39316, rel=2e-4)
    assert v[0] * v[p.n1] < 0 and np.linalg.norm(v[:p.n1]) / np.linalg.norm(v) < 0.8   # both patches, in antiphase


@pytest.mark.parametrize("structure", ["lone", "stack"])
def test_a_single_probe_meets_the_sdm(structure):
    run = next(d for d in DATA["probe"]["runs"] if d["kind"] == "single" and d["structure"] == structure)
    r = np.array(run["r"])
    zf = np.array([complex(*v) for v in run["z"]])
    rs, zs = _series(DATA["sdm"][structure], "zA")
    if structure == "lone":
        zs = zs + _series(DATA["sdm"]["lone"], "zB")[1]
    ff, rf_ = _peak(r[r < 1.15], zf.real[r < 1.15])
    fs, rs_ = _peak(rs[rs < 1.15], zs.real[rs < 1.15])
    assert rf_ == pytest.approx(rs_, rel=0.05)
    assert ff == pytest.approx(fs, rel=0.008)                 # the 4-cell grid reads 0.5-0.7% low


def test_the_quarter_space_probe_is_a_pair():
    """The electric wall at y = 0 mirrors the probe into a second, reversed one: the TM10
    mode's share of the impedance doubles."""
    run = next(d for d in DATA["probe"]["runs"] if d["kind"] == "pair")
    r = np.array(run["r"])
    zf = np.array([complex(*v) for v in run["z"]])
    rs, zs = _series(DATA["sdm"]["stack"], "zA")
    assert _peak(r[r < 1.2], zf.real[r < 1.2])[1] == pytest.approx(2 * _peak(rs[rs < 1.2], zs.real[rs < 1.2])[1], rel=0.05)


def test_the_probe_reactance_is_about_half_the_parallel_plate_formula():
    """Above resonance the FDTD's reactance exceeds the SDM's patch currents (both classes)
    by 13-14 ohm - the probe's own; the parallel-plate formula gives 27. Grid-dependent (a
    one-cell wire) and indicative, which is why the stack's bandwidth is quoted tuned."""
    run = next(d for d in DATA["probe"]["runs"] if d["kind"] == "single" and d["structure"] == "lone")
    fd = {round(r, 2): complex(*z) for r, z in zip(run["r"], run["z"])}
    rest = []
    for d in DATA["sdm"]["lone"]:
        r = round(d["r"], 2)
        if 1.07 <= r <= 1.15:
            rest.append(fd[r].imag - complex(*d["zA"]).imag - complex(*d["zB"]).imag)
    assert 12.0 < min(rest) and max(rest) < 15.5
    assert sdm.probe_reactance(1.1e9, 2.2, 0.0128 * LAM, 0.000432 * LAM) > 1.8 * max(rest)


def test_the_stored_sdm_sweep_reproduces():
    p = sdm.StackedPatch(*STACK, **sdm.BASIS_FULL)
    row = next(d for d in DATA["sdm"]["stack"] if abs(d["r"] - 1.0) < 1e-9)
    v = sdm.probe_vector(p, 1e9, 0.096 * LAM, 0.000432 * LAM)
    z = -v @ np.linalg.solve(p.Z(1e9), v)
    assert z == pytest.approx(complex(*row["zA"]), rel=5e-3)


@pytest.mark.slow
def test_the_coarsest_stack_ringdown_reproduces():
    d = DATA["ringdowns"][0]
    c = d["cells"]
    cell = 0.0128 / c[0]
    got = pf.stacked_ringdown(2.2, c[0], c[1], c[2], 1.0, c[3], c[4], c[5], 0.97 * cell)
    assert [f / cell for f, _ in got] == pytest.approx(d["modes"], rel=1e-6)


@pytest.mark.parametrize("half", [False, True])
def test_the_probe_port_on_a_monopole(half):
    """A 10-cell monopole at 40 cells a wavelength against the wire MoM with the same
    one-cell finite gap: the resistance to 20%, the grid's; the quarter-space grid's wall
    at y = 0 images it into an antiphase pair, which the MoM models too."""
    from otahub.num import mom
    N, hm, jp = 40, 10, 10
    fr = np.array([0.95, 1.05])
    z = pf.probe_impedance(1.0, hm, 0, 0, jp, fr / N, 1.0 / N, air=N // 2, periods=40, half=half)
    for f, zz in zip(fr, z):
        s = f / N
        zs = np.linspace(-hm * s, hm * s, 2 * hm + 1)
        wires = [mom.Wire(np.stack([0 * zs, 0 * zs + y, zs], 1), 0.135 * s) for y in ((0.0,) if half else (jp * s, -jp * s))]
        m = mom.WireModel(wires)
        g1 = [hm - 1, hm]
        v = mom.gap_vector(m, g1)
        if not half:
            v = v - mom.gap_vector(m, [q + 2 * hm for q in g1])
        cur = np.linalg.solve(mom.impedance_matrix(m), v.astype(complex))
        ref = 0.5 / (mom.gap_vector(m, g1) @ cur)
        assert zz.real == pytest.approx(ref.real, rel=0.2)


def test_a_double_tuned_stack_meets_the_fdtd():
    """The design the survey found (gap 0.09, parasitic 1.1 times the driven patch, probe
    near the radiating edge): the FDTD's own probe and the SDM with a tuned series element
    give the same double-tuned band, about 12.4%."""
    from scipy.interpolate import CubicSpline
    from otahub.num import stacked as st
    dt = DATA["double_tuned"]
    r = np.array(dt["fdtd"]["r"])
    zf = np.array([complex(*v) for v in dt["fdtd"]["z"]])
    zs = np.array([complex(*v) for v in dt["sdm"]["zA"]]) + np.array([complex(*v) for v in dt["sdm"]["zB"]])
    peaks = lambda z: [j for j in range(1, len(r) - 1) if z[j].real > z[j - 1].real and z[j].real > z[j + 1].real]
    pf_, ps_ = peaks(zf), peaks(zs)
    assert len(pf_) == 2 and len(ps_) == 2                                # two resistance peaks: double-tuned
    for a, b in zip(pf_, ps_):
        assert zf[a].real == pytest.approx(zs[b].real, rel=0.05)
    fine = np.arange(r[0], r[-1] + 1e-9, 0.0025)
    band_f = st.vswr_band(fine, CubicSpline(r, zf)(fine), 0.0, about=None)
    band_s = st.vswr_band(fine, CubicSpline(r, zs)(fine), 13.5, about=None)
    assert band_f[0] == pytest.approx(band_s[0], abs=0.01)
    assert 0.11 < band_f[0] < 0.14
