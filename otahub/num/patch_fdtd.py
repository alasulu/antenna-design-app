"""A rectangular patch on an infinite grounded slab by FDTD ringdown - an
arbiter for `patch_sdm` that shares nothing with it but the geometry.

The Yee grid, CPML and matrix-pencil fit are `dra`'s, already checked against
the exact sphere pole. One quarter of the space is modelled: a magnetic wall at
x = 0 and an electric wall at y = 0, so with the patch's resonant length along
y only modes whose field is odd in y and even in x exist, of which TM10 is the
lowest. The slab runs into the CPML, so surface waves leave the grid and the Q
it returns is the total - space wave and surface wave together. The patch is a
PEC sheet on a node plane: its tangential E is zeroed after every update.

The sheet's edge sits on the grid, and the resonance converges at first order
in the cell size; `extrapolate` takes runs at several resolutions of the same
geometry to the zero-cell limit. Grid units: cell size 1, c = 1.
"""
from __future__ import annotations

import math

import numpy as np

from .dra import _Yee, matrix_pencil

__all__ = ["ringdown", "extrapolate"]


def ringdown(eps_r: float, nh: int, nl: int, nw: int, f_guess: float, air: int | None = None,
             npml: int = 12, periods: float = 30.0, dtype=np.float32) -> tuple[float, float]:
    """(frequency in cycles per cell-time, total Q) of the patch's TM10 mode.

    nh is the slab height, nl the half length (along y, the resonant direction)
    and nw the half width, all in cells. f_guess need only be near the mode."""
    lam = 1.0 / f_guess
    air = air or int(0.15 * lam) + 4
    shape = (nw + air + npml, nl + air + npml, nh + air + npml)
    g = _Yee(shape, lambda x, y, z: np.where(z < nh, eps_r, 1.0) + 0 * x + 0 * y, npml, dtype)
    tau = 0.5 / f_guess
    t0 = 3.0 * tau
    steps = int((t0 + 3.0 * tau + periods / f_guess) / g.dt)
    start = int((t0 + 3.0 * tau) / g.dt)
    js, ks = nl - 2, nh // 2                 # source and probe under the radiating edge
    probe = np.empty(steps)
    for n in range(steps):
        g.step()
        g.Ex[:nw, :nl + 1, nh] = 0.0         # the sheet: Ex at x = i + 1/2 < nw, y = j <= nl
        g.Ey[:nw + 1, :nl, nh] = 0.0         #            Ey at x = i <= nw, y = j + 1/2 < nl
        t = (n + 1) * g.dt
        g.Ez[0, js, ks] += math.exp(-((t - t0) / tau) ** 2) * math.sin(2 * math.pi * f_guess * (t - t0))
        probe[n] = g.Ez[1, js - 1, ks]
    dec = max(1, int(1.0 / (20.0 * f_guess * g.dt)))
    s, a = matrix_pencil(probe[start::dec], dec * g.dt, modes=20)
    f = s.imag / (2 * math.pi)
    q = s.imag / (-2.0 * s.real)
    ok = (f > 0.5 * f_guess) & (f < 1.6 * f_guess) & (q > 2)
    k = int(np.flatnonzero(ok)[np.argmax(np.abs(a)[ok])])
    return float(f[k]), float(q[k])


def extrapolate(cells, values) -> float:
    """Zero-cell limit of values measured at cell sizes proportional to 1/cells,
    assuming the first-order error the sheet edge gives (a least-squares line in
    1/cells when there are more than two)."""
    x = 1.0 / np.asarray(cells, dtype=float)
    A = np.column_stack([np.ones_like(x), x])
    return float(np.linalg.lstsq(A, np.asarray(values, dtype=float), rcond=None)[0][0])
