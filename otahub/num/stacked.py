"""A probe-fed stacked patch, designed with the two-layer spectral MoM.

The stack's bandwidth is a property of its input impedance, which the spectral MoM
gives through `patch_sdm.StackedPatch` and `probe_vector`: the patch currents'
reaction -v^T Z^-1 v, checked against the FDTD with a single probe (resonant
resistance within 4%). The matrix and the probe couplings are smooth in frequency -
the resonances come from the solve - so they are computed on a coarse grid and
splined (0.05 ohm through resonance at a grid of 0.025 f0).

What the MoM does not give is the probe's own reactance: its self-term diverges
without an attachment mode, and a disc attachment leaves the probe's charge where
the patch's edge-conforming currents cannot move it. The FDTD puts it (with the
x-even currents) near 20 ohm on 0.0128 wavelengths of eps_r 2.2 where the
parallel-plate formula says 27. So the band is taken with a series capacitor tuned
for the widest band, the way a thick probe-fed stack is matched; it moves by under a
point for a non-resonant reactance anywhere from 10 to 35 ohm.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.interpolate import CubicSpline

from . import patch_sdm as sdm

__all__ = ["impedance_sweep", "vswr_band", "design"]


def impedance_sweep(p, f0: float, rs, xps, a: float, rf, progress=None):
    """Patch-current input impedance (len(rf), len(xps)) at frequencies rf * f0 for
    probes at xps (metres from the centre along the resonant direction), radius a, from
    the matrix and couplings at rs * f0, splined."""
    Zs, Vs = [], []
    for r in rs:
        Zs.append(p.Z(r * f0))
        Vs.append(sdm.probe_vector(p, r * f0, xps, a))
        if progress:
            progress(f"solved {r:.3f} f0")
    Zi = CubicSpline(rs, np.array(Zs), axis=0)(rf)
    Vi = CubicSpline(rs, np.array(Vs), axis=0)(rf)
    out = np.empty((len(rf), len(xps)), complex)
    for k in range(len(rf)):
        out[k] = -np.einsum("qi,iq->q", Vi[k], np.linalg.solve(Zi[k], Vi[k].T))
    return out


def vswr_band(rf, z, x_nonres: float = 20.0, z0: float = 50.0, vswr: float = 2.0,
              series=np.linspace(-60.0, 120.0, 181), about: float | None = 1.0):
    """The widest VSWR <= vswr band for the impedance z(rf) plus a non-resonant reactance
    x_nonres * rf and one series element chosen for it: a capacitor of reactance s at f0
    (-j s / rf) for s > 0, an inductor (+j |s| rf) for s < 0. With `about` (f0, by
    default) the band is the one a design there can use - symmetric about it, 2 min(about
    - lo, hi - about) / about; with None, the widest band anywhere, about its own centre.
    (fractional width, centre in rf units, s), or (0, nan, nan) when nothing matches."""
    lim = (vswr - 1) / (vswr + 1)
    best = (0.0, math.nan, math.nan)
    base = z + 1j * x_nonres * rf
    for s_ in series:
        zz = base - 1j * s_ / rf if s_ >= 0 else base - 1j * s_ * rf
        g = np.abs((zz - z0) / (zz + z0))
        ok = np.concatenate([[False], g <= lim, [False]]).astype(int)
        edges = np.flatnonzero(np.diff(ok))
        for lo, hi in zip(edges[::2], edges[1::2] - 1):
            if hi <= lo:
                continue
            if about is None:
                bw, c = 2 * (rf[hi] - rf[lo]) / (rf[hi] + rf[lo]), 0.5 * (rf[hi] + rf[lo])
            elif rf[lo] <= about <= rf[hi]:
                bw, c = 2 * min(about - rf[lo], rf[hi] - about) / about, about
            else:
                continue
            if bw > best[0]:
                best = (bw, c, float(s_))
    return best


def design(eps_r: float, h: float, L: float, W: float, eps_r2: float, h2: float, L2: float, W2: float,
           f0: float, a: float, xp_fractions=(0.45, 0.6, 0.72, 0.84, 0.95), r_lo: float = 0.8,
           r_hi: float = 1.25, x_nonres=(10.0, 20.0, 35.0), progress=None) -> dict:
    """Solve a stack: the widest VSWR-2 band over the probe positions (fractions of the
    driven patch's half length), for each assumed non-resonant reactance, and the
    broadside directivity of the probe-driven currents at the best band's centre."""
    p = sdm.StackedPatch(eps_r, h, L, W, eps_r2, h2, L2, W2, **sdm.BASIS_FULL)
    rs = np.arange(r_lo, r_hi + 1e-9, 0.025)
    rf = np.arange(r_lo, r_hi + 1e-9, 0.0025)
    xps = np.asarray(xp_fractions) * L / 2
    zi = impedance_sweep(p, f0, rs, xps, a, rf, progress)
    rows = []
    for q, frac in enumerate(xp_fractions):
        bands = [vswr_band(rf, zi[:, q], xn) for xn in x_nonres]
        rows.append(dict(xp_fraction=float(frac), bandwidth=[b[0] for b in bands], centre=[b[1] for b in bands],
                         series_reactance_ohm=[b[2] for b in bands]))
    k = int(np.argmax([r["bandwidth"][len(x_nonres) // 2] for r in rows]))
    best = rows[k]
    peaks = [float(rf[j]) for j in range(1, len(rf) - 1)
             if zi[j, k].real > zi[j - 1, k].real and zi[j, k].real > zi[j + 1, k].real]
    out = dict(feeds=rows, best=best, resistance_peaks=peaks, x_nonres=list(x_nonres))
    fc = best["centre"][len(x_nonres) // 2]
    if best["bandwidth"][len(x_nonres) // 2] > 0:
        f = fc * f0
        v = sdm.probe_vector(p, f, xps[k], a)
        cur = -np.linalg.solve(p.Z(f), v)
        p_far, u0 = p.far_field(cur, f)
        out["directivity"] = 4 * math.pi * u0 / p_far
    return out
