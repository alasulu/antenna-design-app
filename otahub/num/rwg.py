"""A flat metal sheet in free space by the method of moments on triangles - for the
antennas no other solver here can hold: the bowtie and the two spirals.

Rao-Wilton-Glisson basis functions on the interior edges of a triangle mesh in
the plane z = 0, Galerkin testing of the mixed-potential EFIE,

    Z_mn = j k eta int int f_m . f_n G  -  j (eta / k) int int (div f_m)(div f_n) G,

with G = exp(-jkR) / (4 pi R). Every pair of triangles is integrated by a 7-point
rule on each; where a source triangle is near the observation point the static
1/R part is taken out and integrated in closed form (Wilton et al. 1984, in the
plane of the sheet), the smooth remainder by the rule. The matrix is assembled
triangle-corner by triangle-corner (M, 3T x 3T) and mapped onto the edge
functions (Z = C^T M C), so the fill is a handful of array operations.

Feeds: a delta gap across a line of edges, or a finite gap - a uniform impressed
field V/g across a band of triangles, the input current the matching weighted
current, so V conj(I)/2 is the delivered power - as `bor`, `strip` and `mom` do.

Lengths in wavelengths, k = 2 pi. It is checked against `strip.RooftopStrip`, a
rectangular-rooftop MoM on the same strip with the same gap, and against the wire
MoM through the strip's equivalent radius w/4.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy import sparse

from ..core.constants import ETA0

__all__ = ["PlanarMesh", "structured", "merge", "strip_mesh", "gap_vector", "edge_gap_vector", "solve", "RwgSolution", "far_field",
           "radiated_power", "tri_potentials"]

K = 2.0 * math.pi

# Dunavant's 7-point rule, degree 5: barycentric coordinates and weights (summing to 1)
_A1, _B1, _W1 = 0.0597158717, 0.4701420641, 0.1323941527
_A2, _B2, _W2 = 0.7974269853, 0.1012865073, 0.1259391805
_BARY = np.array([[1 / 3, 1 / 3, 1 / 3],
                  [_A1, _B1, _B1], [_B1, _A1, _B1], [_B1, _B1, _A1],
                  [_A2, _B2, _B2], [_B2, _A2, _B2], [_B2, _B2, _A2]])
_WQ = np.array([0.225, _W1, _W1, _W1, _W2, _W2, _W2])


@dataclass
class PlanarMesh:
    """Triangles in the plane: nodes (P, 2), tris (T, 3). Triangles are turned
    counter-clockwise; every edge shared by exactly two triangles carries one RWG
    function, from its plus triangle (current leaving the free vertex) to its minus."""
    nodes: np.ndarray
    tris: np.ndarray
    area: np.ndarray = field(init=False)
    qp: np.ndarray = field(init=False)          # quadrature points (T, 7, 2)
    plus: np.ndarray = field(init=False)        # (N, 2): triangle, local index of its free vertex
    minus: np.ndarray = field(init=False)
    length: np.ndarray = field(init=False)
    edge_nodes: np.ndarray = field(init=False)  # (N, 2) the edge's two node indices

    def __post_init__(self) -> None:
        self.nodes = np.asarray(self.nodes, dtype=float)
        tris = np.array(self.tris, dtype=int)
        v = self.nodes[tris]
        cross = (v[:, 1, 0] - v[:, 0, 0]) * (v[:, 2, 1] - v[:, 0, 1]) - (v[:, 1, 1] - v[:, 0, 1]) * (v[:, 2, 0] - v[:, 0, 0])
        if np.any(np.abs(cross) < 1e-15):
            raise ValueError("degenerate triangle")
        flip = cross < 0
        tris[flip] = tris[flip][:, [0, 2, 1]]
        self.tris = tris
        self.area = 0.5 * np.abs(cross)
        v = self.nodes[tris]
        self.qp = np.einsum("qk,tkd->tqd", _BARY, v)
        owners: dict[tuple[int, int], list[tuple[int, int]]] = {}
        for t, (a, b, c) in enumerate(tris):
            for local, (i, j) in ((2, (a, b)), (0, (b, c)), (1, (c, a))):      # local index of the opposite vertex
                owners.setdefault((min(i, j), max(i, j)), []).append((t, local))
        plus, minus, length, en = [], [], [], []
        for (i, j), own in owners.items():
            if len(own) > 2:
                raise ValueError("an edge is shared by more than two triangles")
            if len(own) == 2:
                plus.append(own[0]); minus.append(own[1])
                length.append(float(np.hypot(*(self.nodes[i] - self.nodes[j]))))
                en.append((i, j))
        self.plus, self.minus = np.array(plus), np.array(minus)
        self.length, self.edge_nodes = np.array(length), np.array(en)

    @property
    def n_tri(self) -> int:
        return len(self.tris)

    @property
    def n_basis(self) -> int:
        return len(self.length)

    def spread(self) -> sparse.csr_matrix:
        """C (3T x N): basis n = +l g(T+, free+) - l g(T-, free-), g = (r - v)/(2A)."""
        n = self.n_basis
        rows = np.concatenate([3 * self.plus[:, 0] + self.plus[:, 1], 3 * self.minus[:, 0] + self.minus[:, 1]])
        cols = np.concatenate([np.arange(n), np.arange(n)])
        vals = np.concatenate([self.length, -self.length])
        return sparse.csr_matrix((vals, (rows, cols)), shape=(3 * self.n_tri, n))

    def h(self) -> float:
        """Longest edge."""
        v = self.nodes[self.tris]
        return float(max(np.hypot(*(v[:, i] - v[:, (i + 1) % 3]).T).max() for i in range(3)))


def structured(xy, nu: int, nv: int) -> tuple[np.ndarray, np.ndarray]:
    """Nodes and triangles of a parametric quadrilateral patch: xy(u, v) -> (x, y)
    arrays on [0, 1]^2, cut into nu x nv cells, each split along alternating diagonals."""
    u, v = np.meshgrid(np.linspace(0, 1, nu + 1), np.linspace(0, 1, nv + 1), indexing="ij")
    x, y = xy(u, v)
    nodes = np.stack([x.ravel(), y.ravel()], 1)
    idx = np.arange((nu + 1) * (nv + 1)).reshape(nu + 1, nv + 1)
    tris = []
    for i in range(nu):
        for j in range(nv):
            a, b, c, d = idx[i, j], idx[i + 1, j], idx[i + 1, j + 1], idx[i, j + 1]
            tris += [(a, b, c), (a, c, d)] if (i + j) % 2 == 0 else [(a, b, d), (b, c, d)]
    return nodes, np.array(tris)


def merge(*parts, tol: float = 1e-9) -> PlanarMesh:
    """Join (nodes, tris) patches, fusing coincident nodes so the patches connect."""
    nodes, tris, keymap = [], [], {}
    for pn, pt in parts:
        local = []
        for p in pn:
            key = (round(p[0] / tol), round(p[1] / tol))
            if key not in keymap:
                keymap[key] = len(nodes)
                nodes.append(p)
            local.append(keymap[key])
        tris.append(np.array(local)[pt])
    return PlanarMesh(np.array(nodes), np.concatenate(tris))


def strip_mesh(L: float, W: float, nx: int, ny: int, graded: bool = True) -> PlanarMesh:
    """A strip along x, |x| < L/2, |y| < W/2; nx even puts a column of nodes at x = 0.
    `graded` spaces the nodes across the width as -W/2 cos(pi v), toward the edges where
    the current is singular."""
    ymap = (lambda v: -0.5 * W * np.cos(math.pi * v)) if graded else (lambda v: (v - 0.5) * W)
    return merge(structured(lambda u, v: ((u - 0.5) * L, ymap(v)), nx, ny))


def tri_potentials(v: np.ndarray, rho: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """In the plane of the triangles (d = 0), for triangles v (Q, 3, 2) counter-clockwise
    and points rho (Q, P, 2): S0 = int 1/R dS' and S1 = int (rho' - rho)/R dS'
    (Wilton, Rao, Glisson et al. 1984 with the observation point in the plane)."""
    S0 = np.zeros(rho.shape[:2])
    S1 = np.zeros(rho.shape)
    for i in range(3):
        a, b = v[:, i][:, None, :], v[:, (i + 1) % 3][:, None, :]
        e = b - a
        L = np.linalg.norm(e, axis=-1, keepdims=True)
        lhat = e / L
        uhat = np.stack([lhat[..., 1], -lhat[..., 0]], -1)
        lp = np.sum((b - rho) * lhat, -1)
        lm = np.sum((a - rho) * lhat, -1)
        p = np.sum((a - rho) * uhat, -1)
        Rp = np.linalg.norm(b - rho, axis=-1)
        Rm = np.linalg.norm(a - rho, axis=-1)
        # R + l cancels when the point lies beyond an end of the edge's line; (R + l)(R - l) = p^2
        # gives it without the cancellation, which otherwise rounds to zero near that line
        sp = np.where(lp < 0, p * p / np.maximum(Rp - lp, 1e-300), Rp + lp)
        sm = np.where(lm < 0, p * p / np.maximum(Rm - lm, 1e-300), Rm + lm)
        with np.errstate(divide="ignore", invalid="ignore"):
            lg = np.log(sp / sm)
        lg = np.where((np.abs(p) > 1e-14 * L[..., 0]) & np.isfinite(lg), lg, 0.0)
        S0 += p * lg
        S1 += (0.5 * (p * p * lg + lp * Rp - lm * Rm))[..., None] * uhat
    return S0, S1


def _corner_blocks(mesh: PlanarMesh, near: float = 2.5, chunk: int | None = None):
    """Rows of M (3T x 3T) over triangle corners, a chunk of test triangles at a time:
    M[(p,a),(q,b)] with g = (r - v)/(2A) on each. Yields (first test triangle, rows)."""
    T = mesh.n_tri
    chunk = chunk or max(2, min(24, int(2.5e6 / (49 * T))))    # about a gigabyte of work arrays at most
    V = mesh.nodes[mesh.tris]                                    # (T, 3, 2)
    A = mesh.area
    cen = V.mean(1)
    qp, wq = mesh.qp, _WQ
    ws = wq[None, :] * A[:, None]                                # (T, 7) source weights
    reach = near * mesh.h()
    for p0 in range(0, T, chunk):
        ps = np.arange(p0, min(T, p0 + chunk))
        r = qp[ps]                                               # (c, 7, 2)
        d = r[:, :, None, None, :] - qp[None, None, :, :, :]     # (c, 7, T, 7, 2)
        R = np.linalg.norm(d, axis=-1)
        near_pair = np.linalg.norm(cen[ps][:, None, :] - cen[None, :, :], axis=-1) < reach   # (c, T)
        with np.errstate(divide="ignore", invalid="ignore"):
            G = np.exp(-1j * K * R) / (4 * math.pi * R)
            Gs = np.where(R > 1e-14, (np.exp(-1j * K * R) - 1.0) / (4 * math.pi * R), -1j * K / (4 * math.pi))
        Gk = np.where(near_pair[:, None, :, None], Gs, G)
        I0 = np.einsum("citj,tj->cit", Gk, ws)                                   # (c, 7, T)
        I1 = np.einsum("citj,citjd,tj->citd", Gk, -d, ws)                        # int (r' - r) G
        cp, tq = np.nonzero(near_pair)
        if len(cp):
            S0, S1 = tri_potentials(V[tq], r[cp])                 # (n, 7), (n, 7, 2)
            I0[cp, :, tq] += S0 / (4 * math.pi)
            I1[cp, :, tq, :] += S1 / (4 * math.pi)
        w = wq[None, :] * A[ps][:, None]                          # (c, 7)
        a0 = np.einsum("ci,cit->ct", w, I0)
        a1 = np.einsum("ci,cid,cit->ctd", w, r, I0)
        a2 = np.einsum("ci,ci,cit->ct", w, np.sum(r * r, -1), I0)
        b0 = np.einsum("ci,cid,citd->ct", w, r, I1)
        b1 = np.einsum("ci,citd->ctd", w, I1)
        vp = V[ps]                                                # (c, 3, 2) test corners
        vq = V                                                    # (T, 3, 2) source corners
        vpvq = np.einsum("cad,tbd->ctab", vp, vq)
        Aterm = (b0[:, :, None, None] - np.einsum("cad,ctd->cta", vp, b1)[:, :, :, None]
                 + a2[:, :, None, None] - np.einsum("ctd,tbd->ctb", a1, vq)[:, :, None, :]
                 - np.einsum("cad,ctd->cta", vp, a1)[:, :, :, None] + vpvq * a0[:, :, None, None])
        Aterm = Aterm / (4 * A[ps][:, None, None, None] * A[None, :, None, None])
        Phi = (a0 / (A[ps][:, None] * A[None, :]))[:, :, None, None]
        blk = 1j * K * ETA0 * Aterm - 1j * (ETA0 / K) * Phi       # (c, T, 3, 3)
        yield p0, blk.transpose(0, 2, 1, 3).reshape(3 * len(ps), 3 * T)


def corner_matrix(mesh: PlanarMesh, **kw) -> np.ndarray:
    """M in full - 3T x 3T, so only for small meshes; `impedance_matrix` streams it."""
    T = mesh.n_tri
    M = np.zeros((3 * T, 3 * T), dtype=complex)
    for p0, rows in _corner_blocks(mesh, **kw):
        M[3 * p0:3 * p0 + rows.shape[0]] = rows
    return M


def impedance_matrix(mesh: PlanarMesh, **kw) -> np.ndarray:
    """Z = C^T M C, assembled a chunk of M's rows at a time: each chunk is mapped onto
    the edge functions it touches and added in, so M (9T^2 entries) is never held."""
    C = mesh.spread().tocsc()
    Cr = mesh.spread()
    Z = np.zeros((mesh.n_basis, mesh.n_basis), dtype=complex)
    for p0, rows in _corner_blocks(mesh, **kw):
        sub = Cr[3 * p0:3 * p0 + rows.shape[0]]                   # (3c, N) sparse
        touched = np.unique(sub.indices)
        X = np.asarray((C.T @ rows.T).T)                            # rows @ C: (3c, N)
        Z[touched] += sub[:, touched].toarray().T @ X
    return Z


def gap_vector(mesh: PlanarMesh, x0: float, g: float, direction=(1.0, 0.0), among=None) -> np.ndarray:
    """A finite gap: 1 V impressed uniformly, V/g along `direction`, over the triangles
    whose centroids lie within |(r - x0 direction) . direction| < g/2 - only among the
    triangle indices `among` if given (a spiral's arms cross any straight band through
    its centre). With g -> the width of one cell either side of a line of nodes it
    tends to a delta gap."""
    dvec = np.asarray(direction, float) / np.hypot(*direction)
    V = mesh.nodes[mesh.tris]
    cen = V.mean(1) @ dvec
    inside = np.abs(cen - x0) < g / 2
    if among is not None:
        keep = np.zeros(mesh.n_tri, bool)
        keep[np.asarray(among)] = True
        inside &= keep
    corner = np.zeros(3 * mesh.n_tri)
    for t in np.flatnonzero(inside):
        r = mesh.qp[t]
        for a in range(3):
            gfun = (r - V[t, a]) / (2 * mesh.area[t])            # (7, 2)
            corner[3 * t + a] = mesh.area[t] * np.sum(_WQ * (gfun @ dvec)) / g
    return mesh.spread().T @ corner


def edge_gap_vector(mesh: PlanarMesh, on_line, direction) -> np.ndarray:
    """A delta gap: 1 V across the edges whose both nodes satisfy `on_line(xy) -> bool`,
    positive where the edge function's current (plus triangle to minus) runs along
    `direction`. V_m = +-l_m, and the input current is the same weighting of the edge
    currents (Rao, Wilton and Glisson's feed)."""
    d = np.asarray(direction, float)
    V = mesh.nodes[mesh.tris]
    cen = V.mean(1)
    ok = on_line(mesh.nodes)
    w = np.zeros(mesh.n_basis)
    for n, (i, j) in enumerate(mesh.edge_nodes):
        if ok[i] and ok[j]:
            s = np.sign((cen[mesh.minus[n, 0]] - cen[mesh.plus[n, 0]]) @ d)
            w[n] = s * mesh.length[n]
    if not np.any(w):
        raise ValueError("no edge lies on the feed line")
    return w


@dataclass
class RwgSolution:
    mesh: PlanarMesh
    currents: np.ndarray
    excitation: np.ndarray

    @property
    def input_current(self) -> complex:
        return complex(self.excitation @ self.currents)

    @property
    def input_impedance(self) -> complex:
        return 1.0 / self.input_current

    @property
    def circuit_power(self) -> float:
        return 0.5 * float(np.real(np.conj(self.input_current)))


def solve(mesh: PlanarMesh, excitation: np.ndarray, Z: np.ndarray | None = None) -> RwgSolution:
    from scipy.linalg import solve as _solve
    own = Z is None
    Z = impedance_matrix(mesh) if own else Z
    return RwgSolution(mesh, _solve(Z, excitation.astype(complex), overwrite_a=own, check_finite=False), excitation)


def far_field(sol: RwgSolution, theta, phi):
    """(E_theta, E_phi) up to the common factor -j eta k exp(-jkr)/(4 pi r), per volt."""
    th, ph = np.broadcast_arrays(np.asarray(theta, float), np.asarray(phi, float))
    mesh = sol.mesh
    h = mesh.spread() @ sol.currents                             # corner amplitudes (3T,)
    V = mesh.nodes[mesh.tris]
    J = np.einsum("ta,tqad->tqd", h.reshape(-1, 3), (mesh.qp[:, :, None, :] - V[:, None, :, :])) / (2 * mesh.area[:, None, None])
    wJ = J * (_WQ[None, :, None] * mesh.area[:, None, None])     # (T, 7, 2)
    pts = mesh.qp.reshape(-1, 2)
    wJ = wJ.reshape(-1, 2)
    st, ct, sp, cp = np.sin(th).ravel(), np.cos(th).ravel(), np.sin(ph).ravel(), np.cos(ph).ravel()
    phase = np.exp(1j * K * (np.outer(st * cp, pts[:, 0]) + np.outer(st * sp, pts[:, 1])))
    N = phase @ wJ                                               # (P, 2)
    e_th = (N[:, 0] * cp + N[:, 1] * sp) * ct
    e_ph = -N[:, 0] * sp + N[:, 1] * cp
    return e_th.reshape(th.shape), e_ph.reshape(th.shape)


def radiated_power(sol: RwgSolution, nth: int = 60, nph: int = 72) -> float:
    """Integrated over the sphere: U = r^2 |E|^2 / (2 eta) = eta k^2 |N_perp|^2 / (32 pi^2), per volt."""
    x, w = np.polynomial.legendre.leggauss(nth)
    th = np.arccos(x)
    ph = 2 * math.pi * np.arange(nph) / nph
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    et, ep = far_field(sol, TH, PH)
    U = (np.abs(et) ** 2 + np.abs(ep) ** 2) * ETA0 * K ** 2 / (32 * math.pi ** 2)
    return float(np.sum(U * w[:, None]) * 2 * math.pi / nph)


def directivity(sol: RwgSolution, theta: float, phi: float, nth: int = 60, nph: int = 72) -> float:
    et, ep = far_field(sol, np.array([theta]), np.array([phi]))
    U = (abs(et[0]) ** 2 + abs(ep[0]) ** 2) * ETA0 * K ** 2 / (32 * math.pi ** 2)
    return 4 * math.pi * U / radiated_power(sol, nth, nph)
