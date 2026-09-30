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

from ..core.constants import ETA0
from .dra import _Yee, matrix_pencil

__all__ = ["ringdown", "stacked_ringdown", "probe_impedance", "extrapolate", "shorted_ringdown", "sheet_ringdown", "circular_ringdown", "triangular_ringdown",
           "sheet_directivity"]


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


def stacked_ringdown(eps_r: float, nh: int, nl: int, nw: int, eps_r2: float, nh2: int, nl2: int, nw2: int,
                     f_guess: float, air: int | None = None, npml: int = 12, periods: float = 40.0,
                     n_modes: int = 2, dtype=np.float32) -> list[tuple[float, float]]:
    """[(frequency, total Q)] of a stacked patch's coupled TM10 modes, lowest first: the
    driven sheet (half length nl, half width nw) on nh cells of eps_r, the parasitic
    (nl2, nw2) on nh2 cells of eps_r2 above it, quarter space as `ringdown`. The source and
    probe sit under the driven sheet's radiating edge; the n_modes strongest are kept."""
    lam = 1.0 / f_guess
    air = air or int(0.15 * lam) + 4
    top = nh + nh2
    shape = (max(nw, nw2) + air + npml, max(nl, nl2) + air + npml, top + air + npml)
    g = _Yee(shape, lambda x, y, z: np.where(z < nh, eps_r, np.where(z < top, eps_r2, 1.0)) + 0 * x + 0 * y,
             npml, dtype)
    tau = 0.5 / f_guess
    t0 = 3.0 * tau
    steps = int((t0 + 3.0 * tau + periods / f_guess) / g.dt)
    start = int((t0 + 3.0 * tau) / g.dt)
    js, ks = nl - 2, nh // 2
    probe = np.empty(steps)
    for n in range(steps):
        g.step()
        g.Ex[:nw, :nl + 1, nh] = 0.0
        g.Ey[:nw + 1, :nl, nh] = 0.0
        g.Ex[:nw2, :nl2 + 1, top] = 0.0
        g.Ey[:nw2 + 1, :nl2, top] = 0.0
        t = (n + 1) * g.dt
        g.Ez[0, js, ks] += math.exp(-((t - t0) / tau) ** 2) * math.sin(2 * math.pi * f_guess * (t - t0))
        probe[n] = g.Ez[1, js - 1, ks]
    dec = max(1, int(1.0 / (20.0 * f_guess * g.dt)))
    s, a = matrix_pencil(probe[start::dec], dec * g.dt, modes=20)
    f = s.imag / (2 * math.pi)
    q = s.imag / (-2.0 * s.real)
    # the window matters: a mode of the parasitic near 1.9-1.95 f0 rings harder at the probe
    # than the weak upper coupled mode, and within 2 f_guess displaces it
    ok = np.flatnonzero((f > 0.5 * f_guess) & (f < 1.6 * f_guess) & (q > 2))
    keep = ok[np.argsort(-np.abs(a)[ok])][:n_modes]
    return sorted((float(f[k]), float(q[k])) for k in keep)


def probe_impedance(eps_r: float, nh: int, nl: int, nw: int, jp: int, freqs, f_guess: float,
                    upper: tuple | None = None, R: float = 50.0, half: bool = False, air: int | None = None,
                    npml: int = 12, periods: float = 60.0, dtype=np.float32) -> np.ndarray:
    """Input impedance (ohm) at `freqs` (cycles per cell-time) of a probe feeding the sheet
    of `ringdown` (half length nl along y, half width nw, on nh cells of eps_r), jp cells
    from the centre along the resonant direction, on the centre line x = 0; `upper` =
    (eps_r2, nh2, nl2, nw2) adds the parasitic sheet of `stacked_ringdown`.

    The probe is a wire of zero thickness (Ez held at zero from z = 1 to the sheet), fed in
    the bottom cell by a Norton source: a current with R ohm across it, the resistor
    stepped semi-implicitly so it damps without a stability limit. V is read across that
    cell and I round the wire one cell up, so the source and the gap's own capacitance
    stay outside the impedance. A zero-thickness wire has an equivalent radius set by the
    cell, about a seventh of it.

    The quarter-space grid's electric wall at y = 0 mirrors the probe into a second one,
    reversed, at -jp: two probes in antiphase, which doubles the TM10 mode's share of the
    impedance. `half` models y in full (CPML on both ends, as `shorted_ringdown`), the
    sheets centred on j0: a single probe, as built."""
    lam = 1.0 / f_guess
    air = air or int(0.15 * lam) + 4
    er2, nh2, nl2, nw2 = upper if upper else (1.0, 0, 0, 0)
    top = nh + nh2
    j0 = (npml + air + max(nl, nl2)) if half else 0              # the patch centre on the grid
    ny = (2 * j0) if half else (max(nl, nl2) + air + npml)
    shape = (max(nw, nw2) + air + npml, ny, top + air + npml)
    g = (_YeeHalf if half else _Yee)(shape, lambda x, y, z: np.where(z < nh, eps_r, np.where(z < top, er2, 1.0))
                                     + 0 * x + 0 * y, npml, dtype)
    lo, lo2 = (j0 - nl, j0 - nl2) if half else (0, 0)
    jp = j0 + jp
    c = float(g.cEz[0, jp - 1, 0])
    beta = c * ETA0 / (2.0 * R)
    tau = 0.25 / f_guess
    t0 = 4.0 * tau
    steps = int((t0 + 3.0 * tau + periods / f_guess) / g.dt)
    freqs = np.atleast_1d(np.asarray(freqs, float))
    V = np.zeros(len(freqs), complex)
    I = np.zeros(len(freqs), complex)
    w = -2j * math.pi * freqs * g.dt
    for n in range(steps):
        e_old = float(g.Ez[0, jp, 0])
        g.step()
        curl = (float(g.Ez[0, jp, 0]) - e_old) / c
        t = (n + 1) * g.dt
        js = math.exp(-((t - t0) / tau) ** 2) * math.sin(2 * math.pi * f_guess * (t - t0))
        g.Ez[0, jp, 0] = (e_old * (1.0 - beta) + c * (curl - js)) / (1.0 + beta)
        g.Ez[0, jp, 1:nh] = 0.0
        g.Ex[:nw, lo:j0 + nl + 1, nh] = 0.0
        g.Ey[:nw + 1, lo:j0 + nl, nh] = 0.0
        if upper:
            g.Ex[:nw2, lo2:j0 + nl2 + 1, top] = 0.0
            g.Ey[:nw2 + 1, lo2:j0 + nl2, top] = 0.0
        V += -float(g.Ez[0, jp, 0]) * np.exp(w * (n + 1))
        I += (2.0 * float(g.Hy[0, jp, 1]) - float(g.Hx[0, jp, 1]) + float(g.Hx[0, jp - 1, 1])) * np.exp(w * (n + 0.5))
    return ETA0 * V / I


def extrapolate(cells, values) -> float:
    """Zero-cell limit of values measured at cell sizes proportional to 1/cells,
    assuming the first-order error the sheet edge gives (a least-squares line in
    1/cells when there are more than two)."""
    x = 1.0 / np.asarray(cells, dtype=float)
    A = np.column_stack([np.ones_like(x), x])
    return float(np.linalg.lstsq(A, np.asarray(values, dtype=float), rcond=None)[0][0])


# ---------------------------------------------------------------- a shorted (quarter-wave) patch
#
# A shorted patch has no mirror plane along its length, so y is modelled in full:
# the grid's y = 0 electric wall sits behind a CPML layer on that side too, and
# the short is a finite PEC plane - under the patch only, ground to patch.

def _cpml_both(n: int, npml: int, half: bool, dt: float, m: int = 3, kappa_max: float = 5.0):
    pos = np.arange(n) + (0.5 if half else 0.0)
    last = n if half else n - 1
    rho = np.maximum(np.clip((pos - (last - npml)) / npml, 0.0, 1.0), np.clip((npml - pos) / npml, 0.0, 1.0))
    sig = 0.8 * (m + 1) * rho ** m
    kap = 1.0 + (kappa_max - 1.0) * rho ** m
    b = np.exp(-sig / kap * dt)
    c = np.where(sig > 0, (b - 1.0) / kap, 0.0)
    return kap, b, c, 0


class _YeeHalf(_Yee):
    """Half-space grid: PMC at x = 0, ground at z = 0, CPML on both y ends and at the far x and z ends."""

    def __init__(self, shape, eps_cell, npml: int, dtype=np.float32):
        super().__init__(shape, eps_cell, npml, dtype)
        ny = shape[1]
        self.p_int[1] = _cpml_both(ny + 1, npml, False, self.dt)
        self.p_half[1] = _cpml_both(ny, npml, True, self.dt)


def shorted_ringdown(eps_r: float, nh: int, nl: int, nw: int, f_guess: float, air: int | None = None,
                     npml: int = 12, periods: float = 30.0, dtype=np.float32) -> tuple[float, float]:
    """(frequency in cycles per cell-time, total Q) of a quarter-wave patch shorted along one
    end: nh slab height, nl length from the short to the open edge, nw half width, in cells."""
    lam = 1.0 / f_guess
    air = air or int(0.15 * lam) + 4
    j0 = npml + air                              # the short's plane
    shape = (nw + air + npml, j0 + nl + air + npml, nh + air + npml)
    g = _YeeHalf(shape, lambda x, y, z: np.where(z < nh, eps_r, 1.0) + 0 * x + 0 * y, npml, dtype)
    tau = 0.5 / f_guess
    t0 = 3.0 * tau
    steps = int((t0 + 3.0 * tau + periods / f_guess) / g.dt)
    start = int((t0 + 3.0 * tau) / g.dt)
    js, ks = j0 + nl - 2, nh // 2                # under the open edge
    probe = np.empty(steps)
    for n in range(steps):
        g.step()
        g.Ex[:nw, j0:j0 + nl + 1, nh] = 0.0      # the patch sheet
        g.Ey[:nw + 1, j0:j0 + nl, nh] = 0.0
        g.Ex[:nw, j0, :nh + 1] = 0.0             # the shorting wall: tangential E on the plane y = j0
        g.Ez[:nw + 1, j0, :nh] = 0.0
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


# ---------------------------------------------------------------- any outline: the circle and the triangle
#
# The sheet is a mask over the tangential-E sites of its plane: an Ex edge (x = i + 1/2,
# y = j) or Ey edge (x = i, y = j + 1/2) is zeroed where it lies inside the outline, so a
# curved or slanted edge becomes a staircase that converges as the cells shrink. The
# x = 0 plane is a magnetic wall either way; `half` models y in full (CPML on both
# ends, like `shorted_ringdown`) for outlines with no mirror plane along y.

def _sheet_run(eps_r: float, nh: int, inside, x_max: float, y_min: float, y_max: float, f_guess: float,
               src: tuple[int, int], half: bool = False, short_y: float | None = None, box_gap: int | None = None,
               air: int | None = None, npml: int = 12, periods: float = 30.0, dtype=np.float32):
    """The engine behind `sheet_ringdown` and `sheet_directivity`: returns (f, Q, box), box
    being a `horn_fdtd.BoxTransform` round the sheet (box_gap cells clear of it) whose
    frequency is set after the run - it records the whole run at a band of trial
    frequencies and keeps the one nearest the resonance found."""
    lam = 1.0 / f_guess
    air = air or int(0.15 * lam) + 4
    j0 = (npml + air - int(math.floor(y_min))) if half else 0          # the outline's y = 0 on the grid
    nx = int(math.ceil(x_max)) + air + npml
    ny = j0 + int(math.ceil(y_max)) + air + npml
    shape = (nx, ny, nh + air + npml)
    grid = _YeeHalf if half else _Yee
    g = grid(shape, lambda x, y, z: np.where(z < nh, eps_r, 1.0) + 0 * x + 0 * y, npml, dtype)
    i = np.arange(nx + 1)
    j = np.arange(ny + 1)
    ex = inside((i[:nx] + 0.5)[:, None], (j[:, None] - j0).T)           # (nx, ny + 1)
    ey = inside(i[:, None] + 0 * j[None, :ny], (j[None, :ny] + 0.5 - j0) + 0 * i[:, None])   # (nx + 1, ny)
    ex, ey = np.broadcast_to(ex, (nx, ny + 1)), np.broadcast_to(ey, (nx + 1, ny))
    if short_y is not None:                      # a full-height wall on the grid plane y = short_y, where the sheet meets it
        js = int(round(short_y)) + j0
        wx = ex[:, js].copy()                    # Ex edges of the wall, ground to sheet
        wz = np.zeros(nx + 1, bool)              # Ez sites beside them
        wz[:nx] |= wx
        wz[1:] |= wx
    tau = 0.5 / f_guess
    t0 = 3.0 * tau
    steps = int((t0 + 3.0 * tau + periods / f_guess) / g.dt)
    start = int((t0 + 3.0 * tau) / g.dt)
    si, sj = src[0], src[1] + j0
    ks = nh // 2
    dec = max(1, int(1.0 / (20.0 * f_guess * g.dt)))
    probe = np.empty(steps)
    boxes = []
    if box_gap is not None:
        from .horn_fdtd import BoxTransform
        bi0 = int(math.ceil(x_max)) + box_gap
        bj0 = j0 + int(math.ceil(y_max)) + box_gap
        bjl = (j0 + int(math.floor(y_min)) - box_gap) if half else None
        bk0 = nh + box_gap
        boxes = [BoxTransform(bi0, bj0, bk0, f_guess * r, g.dt, j_lo=bjl) for r in np.linspace(0.8, 1.2, 41)]
    for n in range(steps):
        g.step()
        g.Ex[:, :, nh][ex] = 0.0
        g.Ey[:, :, nh][ey] = 0.0
        if short_y is not None:
            g.Ex[:, js, :nh + 1][wx] = 0.0
            g.Ez[:, js, :nh][wz] = 0.0
        t = (n + 1) * g.dt
        g.Ez[si, sj, ks] += math.exp(-((t - t0) / tau) ** 2) * math.sin(2 * math.pi * f_guess * (t - t0))
        probe[n] = g.Ez[si + 1, sj - 1, ks]
        for b in boxes:
            b.record(g, n)
    s, a = matrix_pencil(probe[start::dec], dec * g.dt, modes=20)
    f = s.imag / (2 * math.pi)
    q = s.imag / (-2.0 * s.real)
    ok = (f > 0.5 * f_guess) & (f < 1.6 * f_guess) & (q > 2)
    k = int(np.flatnonzero(ok)[np.argmax(np.abs(a)[ok])])
    box = min(boxes, key=lambda b: abs(b.f - f[k])) if boxes else None
    return float(f[k]), float(q[k]), box


def sheet_ringdown(eps_r: float, nh: int, inside, x_max: float, y_min: float, y_max: float, f_guess: float,
                   src: tuple[int, int], half: bool = False, air: int | None = None, npml: int = 12,
                   periods: float = 30.0, dtype=np.float32) -> tuple[float, float]:
    """(frequency in cycles per cell-time, total Q) of the dominant mode of a sheet on the
    slab whose outline is inside(x, y) -> bool, in cells, x >= 0 (the x = 0 magnetic wall
    halves it), y from y_min to y_max. src = (i, j): the source and probe column, in the
    outline's own coordinates, near an edge where the mode is strong."""
    f, q, _ = _sheet_run(eps_r, nh, inside, x_max, y_min, y_max, f_guess, src, half=half, air=air, npml=npml,
                         periods=periods, dtype=dtype)
    return f, q


def circular_ringdown(eps_r: float, nh: int, radius: float, f_guess: float, **kw) -> tuple[float, float]:
    """The circular patch's TM11 mode: a disc of `radius` cells (not necessarily whole),
    centred on the origin - odd in y and even in x, so the quarter-space grid holds it."""
    return sheet_ringdown(eps_r, nh, lambda x, y: x * x + y * y <= radius * radius, radius, 0.0, radius,
                          f_guess, (0, int(radius) - 2), **kw)


def triangular_ringdown(eps_r: float, nh: int, side: float, f_guess: float, **kw) -> tuple[float, float]:
    """The equilateral triangle's TM10 mode, the member of the degenerate pair even about
    the median through its apex: apex on +y, centroid at the origin, side `side` cells."""
    ht = side * math.sqrt(3.0) / 2.0
    ya, yb = 2.0 * ht / 3.0, -ht / 3.0

    def inside(x, y):
        return (y >= yb) & (y <= ya) & (np.abs(x) <= 0.5 * side * (ya - y) / ht)
    return sheet_ringdown(eps_r, nh, inside, 0.5 * side, yb, ya, f_guess, (0, int(ya) - 2), half=True, **kw)


def sheet_directivity(eps_r: float, nh: int, inside, x_max: float, y_min: float, y_max: float, f_guess: float,
                      src: tuple[int, int], half: bool = False, short_y: float | None = None, box_gap: int = 6,
                      **kw) -> dict:
    """Ring-down plus the closed box transform (horn_fdtd.BoxTransform) round the sheet: the
    resonance, total Q and the directivity at broadside, over the half space above the
    ground. The box's side faces cross the slab, where free-space equivalent currents are
    not exact - use it on air plates (a PIFA), where they are. The DFT frequency is the
    nearest of 41 trial frequencies (0.8-1.2 f_guess) to the resonance, within 0.5% of it."""
    f, q, box = _sheet_run(eps_r, nh, inside, x_max, y_min, y_max, f_guess, src, half=half, short_y=short_y,
                           box_gap=box_gap, **kw)
    d_peak, th, ph = box.peak()
    return dict(f=f, q=q, directivity=box.directivity(), peak=d_peak, peak_theta=th, peak_phi=ph, f_dft=box.f)
