"""Cavity model of a corner-truncated square patch, by finite elements.

An arbiter for the single-feed circularly polarised patch. Under a thin patch
the field is E_z = psi(x, y) with magnetic side walls, so the modes are the
Neumann eigenfunctions of the patch outline:

    -laplacian(psi) = k^2 psi   inside,    d psi / dn = 0   on the edge.

Solved here with linear triangles on the truncated square itself, so nothing
about the cut is assumed: not the first-order perturbation that gives the
classic dS/S = 1/(2Q) rule, not which modes the cut splits, not where the
circular polarisation lands. A probe at a node excites mode m with amplitude
psi_m(feed) / (k_m^2 - k^2 (1 - j/Q)), one loss Q shared by every mode as in
the usual cavity model. Each mode radiates at broadside through its edge
magnetic current, whose moment is the closed integral of psi_m along the edge;
the broadside polarisation is the sum of those moments weighted by the mode
amplitudes, and the axial ratio follows from it.

Geometry is in units of the side, L = 1; the uncut square's first resonance is
k = pi, and `TruncatedSquare.k_square` is the mesh's own value of it, which is
what frequency offsets should be measured against. The corners at (0, 0) and
(1, 1) are cut by right isosceles triangles of leg c = m/n, removing c^2 of
area in all. The mesh follows the cut exactly and splits every cell along the
anti-diagonal. That keeps the x <-> y reflection but not the 90-degree
rotation, and the reflection alone does not hold cos(pi x) +- cos(pi y)
together (they are even and odd under it), so the uncut square's pair splits
on the mesh too - by 8e-5 at 4 cells a side, 6e-7 at 12, falling as h^4,
against the 2.7e-2 of a 1/12 cut. The cut's split is the one measured.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from scipy.optimize import minimize_scalar

__all__ = ["TruncatedSquare", "truncated_square", "mesh"]


def mesh(n: int, m: int):
    """Nodes, CCW triangles and a node index grid for the truncated square."""
    if not 0 <= m < n:
        raise ValueError("need 0 <= m < n")
    ij = np.array([(i, j) for i in range(n + 1) for j in range(n + 1)
                   if m <= i + j <= 2 * n - m])
    index = -np.ones((n + 1, n + 1), dtype=int)
    index[ij[:, 0], ij[:, 1]] = np.arange(len(ij))
    tris = []
    for i in range(n):
        for j in range(n):
            s = i + j
            if s >= m and s + 1 <= 2 * n - m:
                tris.append((index[i, j], index[i + 1, j], index[i, j + 1]))
            if s + 1 >= m and s + 2 <= 2 * n - m:
                tris.append((index[i + 1, j], index[i + 1, j + 1], index[i, j + 1]))
    return ij / n, np.array(tris), index


def _assemble(p: np.ndarray, t: np.ndarray):
    x, y = p[t, 0], p[t, 1]
    b = np.stack([y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]], 1)
    c = np.stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]], 1)
    area = 0.5 * (b[:, 0] * c[:, 1] - b[:, 1] * c[:, 0])
    ke = (b[:, :, None] * b[:, None, :] + c[:, :, None] * c[:, None, :]) / (4 * area[:, None, None])
    me = area[:, None, None] / 12.0 * (np.ones((3, 3)) + np.eye(3))[None]
    rows = np.repeat(t, 3, axis=1).ravel()
    cols = np.tile(t, (1, 3)).ravel()
    size = (len(p), len(p))
    return (sp.coo_matrix((ke.ravel(), (rows, cols)), size).tocsc(),
            sp.coo_matrix((me.ravel(), (rows, cols)), size).tocsc())


def _modes(p, t, count):
    stiff, mass = _assemble(p, t)
    lam, vec = sla.eigsh(stiff, k=count, M=mass, sigma=-1.0)
    order = np.argsort(lam)
    return lam[order], vec[:, order]


def _edge_moments(p, t, vec):
    """Closed integral of psi_m along the outline, CCW: each mode's radiating moment."""
    e = np.concatenate([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
    key = np.sort(e, axis=1)
    _, inv, cnt = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    edge = e[cnt[inv.ravel()] == 1]
    step = p[edge[:, 1]] - p[edge[:, 0]]
    return (0.5 * (vec[edge[:, 0]] + vec[edge[:, 1]])).T @ step


@dataclass
class TruncatedSquare:
    n: int
    m: int
    points: np.ndarray
    index: np.ndarray
    lam: np.ndarray          # eigenvalues k_m^2, ascending; lam[0] = 0 is the static mode
    modes: np.ndarray        # mass-orthonormal eigenvectors, one column per mode
    moments: np.ndarray      # (modes, 2) broadside radiating moments
    k_square: float          # the UNCUT square's first resonance on the same mesh

    @property
    def cut_over_side(self) -> float:
        return self.m / self.n

    @property
    def split(self) -> float:
        """Fractional separation of the two modes the cut splits apart, relative
        to their geometric mean - the frequency CP forms near, so that perfect
        CP needs Q close to 1/split."""
        kl, kh = math.sqrt(self.lam[1]), math.sqrt(self.lam[2])
        return (kh - kl) / math.sqrt(kl * kh)

    def _node(self, feed):
        i, j = (int(round(v * self.n)) for v in feed)
        node = self.index[i, j]
        if node < 0:
            raise ValueError(f"feed {feed} is outside the patch")
        return node

    def broadside_field(self, k: float, q: float, feed) -> np.ndarray:
        """Complex broadside polarisation vector (up to a common factor)."""
        psi = self.modes[self._node(feed)]
        return (psi / (self.lam - k * k * (1.0 - 1j / q))) @ self.moments

    def axial_ratio_db(self, k: float, q: float, feed) -> float:
        ex, ey = self.broadside_field(k, q, feed)
        total = abs(ex) ** 2 + abs(ey) ** 2
        circ = 2.0 * np.imag(ex * np.conj(ey))
        lin = math.sqrt(max(total * total - circ * circ, 0.0))
        return 20.0 * math.log10(math.sqrt((total + lin) / max(total - lin, 1e-300)))

    def pair_impedance(self, k: float, q: float, feed) -> complex:
        """The two split modes' impedance at the probe, relative scale, with the
        probe's own reactance (every other mode) taken as tuned out."""
        psi = self.modes[self._node(feed)][1:3]
        return complex(1j * k * np.sum(psi ** 2 / (self.lam[1:3] - k * k * (1.0 - 1j / q))))

    def best_axial_ratio(self, q: float, feed) -> tuple[float, float]:
        """(lowest AR in dB, the k where it occurs) between the two split modes."""
        lo, hi = math.sqrt(self.lam[1]), math.sqrt(self.lam[2])
        r = minimize_scalar(lambda k: self.axial_ratio_db(k, q, feed), bounds=(lo, hi),
                            method="bounded", options={"xatol": 1e-12})
        return float(r.fun), float(r.x)

    def circular_q(self, feed=(0.2, 0.5)) -> tuple[float, float]:
        """The loss Q at which this cut gives perfect circular polarisation, and
        the k of that CP centre."""
        q0 = 1.0 / self.split
        r = minimize_scalar(lambda q: self.best_axial_ratio(q, feed)[0],
                            bounds=(0.9 * q0, 1.1 * q0), method="bounded",
                            options={"xatol": 1e-6})
        return float(r.x), self.best_axial_ratio(r.x, feed)[1]


def truncated_square(n: int, m: int, count: int = 30) -> TruncatedSquare:
    """Solve the truncated square with n cells per side and a cut of m cells."""
    p, t, index = mesh(n, m)
    lam, vec = _modes(p, t, count)
    p0, t0, _ = mesh(n, 0)
    k_square = math.sqrt(_modes(p0, t0, 3)[0][1])
    return TruncatedSquare(n, m, p, index, lam, vec, _edge_moments(p, t, vec), k_square)
