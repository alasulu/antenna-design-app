"""A flat strip in free space - the Babinet complement of a slot - by two methods
of moments that share nothing but the geometry, and the open tube it is said to
be equivalent to.

Why this exists. The slot family takes its resonant length and resistance from
the complementary dipole, a strip of width w, through the round-wire laws at
the equivalent radius a = w/4. That radius is exact for the strip's body in the
thin limit (the static kernel); nothing had checked it on WIDE strips, where the
ends are no longer negligible.

`SpectralStrip` is `patch_sdm`'s Galerkin with the ground and slab removed: the
sheet sees free space on both sides, so the TM and TE impedances at its plane are
kz/(2 omega eps0) and omega mu0/(2 kz). The edge-exact entire-domain basis, the
contour past the branch point and the 1/k_max extrapolation are `patch_sdm`'s.
It has no feed: its answer is the natural (characteristic-mode) resonance.

`RooftopStrip` is a space-domain mixed-potential MoM (Glisson and Wilton):
rooftops on a rectangular mesh graded toward the edges, pulse charges on the
cells, the static part of each cell integral in closed form and the rest by
Gauss-Legendre. It converges at first order in the cell size, so answers are
taken from two or more meshes to the zero-cell limit (`patch_fdtd.extrapolate`).
It carries a finite gap feed - a uniform impressed field V/g across the full
width over |x| < g/2, the input current the gap-averaged current - as `bor` does
for a tube, so strip and tube can be compared fed identically.

`tube_profile` is the open tube for `bor`, its segments halved toward both rims:
the current on a thin tube's open rim is edge-singular, and an ungraded mesh
converges slowly there (the natural resonance at a = 0.01 moved 0.4486, 0.4467,
0.4458, 0.4449 as the segment halved from 0.02 to 0.0025, where a graded one
sits at 0.4446 from 0.01 down).

Lengths in wavelengths, k = 2 pi.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq

from ..core.constants import C0, EPS0, ETA0, MU0
from . import bor
from .patch_sdm import BASIS_FULL, RectPatch

__all__ = ["SpectralStrip", "RooftopStrip", "rect_potential", "dominant", "natural_length",
           "tube_profile", "tube_impedance"]

K = 2.0 * math.pi
_GL = np.polynomial.legendre.leggauss(4)


class SpectralStrip(RectPatch):
    """Length L along x, width W along y, in wavelengths; no ground, no slab."""

    def __init__(self, L: float, W: float, bx=BASIS_FULL["bx"], by=BASIS_FULL["by"]):
        super().__init__(1.0, 1.0, L, W, bx=bx, by=by)      # eps_r = 1 keeps the contour past k0

    def znode(self, kr, k0: float):
        w = k0 * C0
        kz = -1j * np.sqrt(kr * kr - k0 * k0 + 0j)
        return kz / (2 * w * EPS0), w * MU0 / (2 * kz)

    def matrix(self, kmax: float = 1600.0):
        """Galerkin matrix at f = c, so that metres are wavelengths."""
        return self.Z(C0, kmax=kmax)

    def far_field(self, *args, **kwargs):
        raise NotImplementedError("patch_sdm's far field is the slab's; a free strip has none here")


def _F(x, y):
    ax, ay = np.abs(x), np.abs(y)
    t1 = np.where(ax > 0, x * np.arcsinh(y / np.where(ax > 0, ax, 1.0)), 0.0)
    t2 = np.where(ay > 0, y * np.arcsinh(x / np.where(ay > 0, ay, 1.0)), 0.0)
    return t1 + t2


def rect_potential(px, py, x1, x2, y1, y2):
    """Integral of exp(-jkR)/(4 pi R) over the rectangle [x1,x2]x[y1,y2] at in-plane
    points (px, py), broadcasting: 1/R in closed form, the smooth rest by 4x4 Gauss."""
    a1, a2, b1, b2 = x1 - px, x2 - px, y1 - py, y2 - py
    static = (_F(a2, b2) - _F(a1, b2) - _F(a2, b1) + _F(a1, b1)) / (4 * math.pi)
    t, w = _GL
    acc = 0.0
    for ti, wi in zip(t, w):
        xx = 0.5 * (a1 + a2) + 0.5 * (a2 - a1) * ti
        for tj, wj in zip(t, w):
            yy = 0.5 * (b1 + b2) + 0.5 * (b2 - b1) * tj
            R = np.sqrt(xx * xx + yy * yy)
            Rs = np.where(R > 1e-12, R, 1.0)
            acc = acc + wi * wj * np.where(R > 1e-12, (np.exp(-1j * K * Rs) - 1) / Rs, -1j * K)
    return static + acc * 0.25 * (a2 - a1) * (b2 - b1) / (4 * math.pi)


class RooftopStrip:
    """Length L along x, width W along y, in wavelengths, on nx x ny cells: cosine
    grading across the width, and along the length a blend (weight `xgrade`) of
    uniform and cosine spacing."""

    def __init__(self, L: float, W: float, nx: int, ny: int, xgrade: float = 0.5, gap: tuple[float, int] | None = None):
        if gap is None:
            j = np.arange(nx + 1)
            self.xs = (1 - xgrade) * np.linspace(-0.5 * L, 0.5 * L, nx + 1) - xgrade * 0.5 * L * np.cos(math.pi * j / nx)
        else:
            # gap = (g, m): |x| < g/2 in 2m equal cells; each side beyond it nx/2 - m cells, graded toward the end
            g, m = gap
            t = np.linspace(0.0, 1.0, nx // 2 - m + 1)
            out = 0.5 * g + (0.5 * (L - g)) * ((1 - xgrade) * t + xgrade * np.sin(0.5 * math.pi * t))
            half = np.concatenate([np.linspace(0.0, 0.5 * g, m + 1)[:-1], out])
            self.xs = np.concatenate([-half[::-1], half[1:]])
            nx = len(self.xs) - 1
        self.ys = -0.5 * W * np.cos(math.pi * np.arange(ny + 1) / ny)
        self.L, self.W, self.nx, self.ny = L, W, nx, ny
        self.xc, self.yc = 0.5 * (self.xs[1:] + self.xs[:-1]), 0.5 * (self.ys[1:] + self.ys[:-1])
        self.dx, self.dy = np.diff(self.xs), np.diff(self.ys)
        cid = lambda i, jj: i * ny + jj
        self.bx = np.array([(i, jj) for i in range(1, nx) for jj in range(ny)])   # x-edge at xs[i], row jj
        self.by = np.array([(i, jj) for i in range(nx) for jj in range(1, ny)])   # y-edge at ys[jj], column i
        # each rooftop carries unit current out of its `minus` cell into its `plus` cell
        self.minus = np.array([cid(i - 1, jj) for i, jj in self.bx] + [cid(i, jj - 1) for i, jj in self.by])
        self.plus = np.array([cid(i, jj) for i, jj in self.bx] + [cid(i, jj) for i, jj in self.by])

    def matrix(self) -> np.ndarray:
        xs, ys, xc, yc, dx, dy = self.xs, self.ys, self.xc, self.yc, self.dx, self.dy
        ci, cj = np.divmod(np.arange(self.nx * self.ny), self.ny)
        # cell-averaged potential kernel, observation at cell centres
        P = rect_potential(xc[ci][:, None], yc[cj][:, None], xs[ci][None, :], xs[ci + 1][None, :],
                           ys[cj][None, :], ys[cj + 1][None, :]) / (dx[ci] * dy[cj])[None, :]
        mi, pl = self.minus, self.plus
        Zphi = P[np.ix_(pl, pl)] - P[np.ix_(pl, mi)] - P[np.ix_(mi, pl)] + P[np.ix_(mi, mi)]
        nbx, nb = len(self.bx), len(mi)
        Za = np.zeros((nb, nb), complex)
        i, jj = self.bx[:, 0], self.bx[:, 1]
        Za[:nbx, :nbx] = (xc[i] - xc[i - 1])[:, None] * rect_potential(
            xs[i][:, None], yc[jj][:, None], xc[i - 1][None, :], xc[i][None, :],
            ys[jj][None, :], ys[jj + 1][None, :]) / dy[jj][None, :]
        if len(self.by):
            i, jj = self.by[:, 0], self.by[:, 1]
            Za[nbx:, nbx:] = (yc[jj] - yc[jj - 1])[:, None] * rect_potential(
                xc[i][:, None], ys[jj][:, None], xs[i][None, :], xs[i + 1][None, :],
                yc[jj - 1][None, :], yc[jj][None, :]) / dx[i][None, :]
        Z = 1j * ETA0 * K * Za + ETA0 * Zphi / (1j * K)      # j omega mu = j eta k, 1/(j omega eps) = eta/(j k)
        return 0.5 * (Z + Z.T)                                # reciprocity; the razor testing is not quite symmetric

    def gap(self, cells: int) -> float:
        """Length of the gap spanning `cells` cells either side of the centre (nx even)."""
        c = self.nx // 2
        return float(self.xs[c + cells] - self.xs[c - cells])

    def gap_vector(self, g: float) -> np.ndarray:
        """(1/g) times the overlap of each x-rooftop with |x| < g/2: the excitation
        for 1 V across the gap, and the weights of the gap-averaged current."""
        e = np.zeros(len(self.minus))
        xs = self.xs

        def part(a, b, rising):
            lo, hi = max(a, -g / 2), min(b, g / 2)
            if hi <= lo:
                return 0.0
            if rising:
                return ((hi - a) ** 2 - (lo - a) ** 2) / (2 * (b - a))
            return ((b - lo) ** 2 - (b - hi) ** 2) / (2 * (b - a))

        for m, (i, _) in enumerate(self.bx):
            e[m] = (part(xs[i - 1], xs[i], True) + part(xs[i], xs[i + 1], False)) / g
        return e

    def input_impedance(self, gap_cells: int) -> complex:
        e = self.gap_vector(self.gap(gap_cells))
        return complex(1.0 / (e @ np.linalg.solve(self.matrix(), e.astype(complex))))


def dominant(Z: np.ndarray):
    """(eigenvalue, current) of the characteristic mode that radiates most. R is
    lifted by its own round-off so that the generalised problem stays definite."""
    R = Z.real
    ev = np.linalg.eigvalsh(R)
    R = R + 2 * (max(0.0, -ev[0]) + 1e-9 * ev[-1]) * np.eye(len(R))
    lam, V = eigh(Z.imag, R)
    j = int(np.argmin(np.linalg.norm(V, axis=0)))
    return lam[j], V[:, j]


def natural_length(matrix_of_L, L0: float, span: float = 0.02) -> float:
    """Length at which the dominant characteristic eigenvalue crosses zero, for a
    function returning the MoM matrix at a given length (it rises with L)."""
    f = lambda L: dominant(matrix_of_L(L))[0]
    lo, hi = L0 - span, L0 + span
    while f(lo) > 0:
        lo, hi = lo - span, lo
    while f(hi) < 0:
        lo, hi = hi, hi + span
    return brentq(f, lo, hi, xtol=2e-6)


def tube_profile(L: float, a: float, seg: float, gap: float, levels: int = 8, gap_seg: int = 4) -> bor.Profile:
    """Open tube of radius a and length L for `bor`, fed across a gap at the centre,
    its segments halved `levels` times toward both rims."""
    ends = np.cumsum([seg * 2.0 ** -k for k in range(levels, 0, -1)])
    v = ([(a, -L / 2)] + [(a, -L / 2 + z) for z in ends] + [(a, -gap / 2), (a, gap / 2)]
         + [(a, L / 2 - z) for z in ends[::-1]] + [(a, L / 2)])
    return bor.profile(v, seg, gap=(levels + 1, levels + 2), gap_seg=gap_seg)


def tube_impedance(L: float, a: float, gap: float, seg: float = 0.005) -> complex:
    return bor.solve(tube_profile(L, a, min(seg, 2 * a), gap)).input_impedance
