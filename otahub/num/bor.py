"""Bodies of revolution by the method of moments - solid cones, discs and tubes.

Why this exists. A discone, a biconical antenna or a conical monopole is a
SURFACE, and the wire-cage stand-ins `mom` can build do not converge in
impedance: a 47 degree cage still moved 3-5% between 24 and 32 wires. A body of
revolution driven symmetrically carries only an axisymmetric current along its
generating curve, so the surface problem is one-dimensional and can be solved
exactly - no cage, no thin-wire kernel.

Formulation. EFIE, mixed-potential form, rotationally symmetric (m = 0) current
flowing along the generating curve. The unknown is the TOTAL ring current I(t)
at arc length t, rooftop basis, Galerkin testing, exactly as `mom` does for a
wire, with the wire kernel replaced by its ring averages:

    K0(t, t') = (1/2pi) int exp(-jkR)/R dphi'
    K1(t, t') = (1/2pi) int cos(phi') exp(-jkR)/R dphi'

K1 carries the radial parts of the two tangents, which turn by phi' round the
ring; the axial parts do not. For each phi' a straight generating segment is a
straight line in space, so the static part 1/R is integrated along it in
closed form (asinh and sqrt, as in `mom`), and the logarithmic singularity that
leaves in phi' is removed by phi' = pi u^2. Only the smooth remainder
(exp(-jkR) - 1)/R goes to Gauss-Legendre.

Feed. A FINITE gap: a uniform impressed field V/g over the stretch of the
generating curve that forms the gap. The input current is the gap-averaged
current, which makes V*conj(I)/2 the power the source delivers. A delta gap is
the limit g -> 0 and inherits its trouble on a fat body - the gap capacitance
diverges - so a gap of physical length is the honest model of a real feed.

Geometry in wavelengths, k = 2*pi, as in `mom`. A disc or a cone reaching the
axis must stop at a small radius: the ring there has no circumference.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from ..core.constants import ETA0

__all__ = ["Profile", "profile", "solve", "BorSolution", "far_field", "radiated_power", "directivity"]

K = 2.0 * math.pi
_GL_OUT = np.polynomial.legendre.leggauss(8)
_GL_IN = np.polynomial.legendre.leggauss(12)
_GL_PHI = np.polynomial.legendre.leggauss(40)


@dataclass
class Profile:
    """An open generating curve (rho, z), rho >= 0, as a polyline of segments.

    The current vanishes at both ends. `gap` is a pair of arc lengths (t0, t1)
    along the curve that carry the impressed field; it must fall on segment
    boundaries, which `profile` arranges.
    """
    nodes: np.ndarray
    gap: tuple[float, float]
    seg_len: np.ndarray = field(init=False)
    seg_t: np.ndarray = field(init=False)
    t_node: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        d = np.diff(self.nodes, axis=0)
        self.seg_len = np.linalg.norm(d, axis=1)
        if np.any(self.seg_len <= 0):
            raise ValueError("zero-length segment in profile")
        if np.any(self.nodes[:, 0] <= 0):
            raise ValueError("a profile must stay off the axis (rho > 0)")
        self.seg_t = d / self.seg_len[:, None]
        self.t_node = np.concatenate([[0.0], np.cumsum(self.seg_len)])

    @property
    def n_seg(self) -> int:
        return len(self.seg_len)

    @property
    def n_basis(self) -> int:
        return self.n_seg - 1


def profile(vertices, seg: float, gap: tuple[int, int] | None = None, gap_seg: int = 2) -> Profile:
    """Polyline through `vertices` [(rho, z), ...], each straight run cut into
    segments no longer than `seg`. `gap` = (i, i+1) marks the run between two
    vertices as the feed gap; it is cut into `gap_seg` segments whatever its
    length."""
    v = np.asarray(vertices, dtype=float)
    pts, t0, t1 = [v[0]], None, None
    for i in range(len(v) - 1):
        L = float(np.linalg.norm(v[i + 1] - v[i]))
        is_gap = gap is not None and i == gap[0]
        n = gap_seg if is_gap else max(1, int(math.ceil(L / seg)))
        if is_gap:
            t0 = sum(np.linalg.norm(np.diff(np.array(pts), axis=0), axis=1)) if len(pts) > 1 else 0.0
        for j in range(1, n + 1):
            pts.append(v[i] + (v[i + 1] - v[i]) * j / n)
        if is_gap:
            t1 = t0 + L
    if gap is None:
        raise ValueError("a profile needs a feed gap")
    return Profile(np.array(pts), (float(t0), float(t1)))


def _moments(pr: Profile):
    """M[w][a][b][p, q] = int_p int_q s^a s'^b K_w ds' ds, w = 0 (K0) or 1 (K1),
    s and s' measured from each segment's start."""
    ns = pr.n_seg
    xo, wo = _GL_OUT
    xi, wi = _GL_IN
    xp, wp = _GL_PHI
    u = 0.5 * (xp + 1.0)
    phi = math.pi * u ** 2
    wphi = wp * u                       # (1/2pi) int_0^2pi dphi = int_0^1 2u du; GL on [0,1] halves the weights
    so = 0.5 * pr.seg_len[:, None] * (xo[None, :] + 1.0)                  # (ns, no)
    obs = pr.nodes[:-1, None, :] + so[:, :, None] * pr.seg_t[:, None, :]  # (ns, no, 2)
    orho = obs[..., 0].reshape(-1)
    oz = obs[..., 1].reshape(-1)
    P = orho.size
    ob3 = np.stack([orho, np.zeros(P), oz], axis=1)                       # observer at phi = 0
    cph, sph = np.cos(phi), np.sin(phi)
    M = np.zeros((2, 2, 2, ns, ns), dtype=complex)
    for q in range(ns):
        (ra, za), L = pr.nodes[q], pr.seg_len[q]
        dr, dz = pr.seg_t[q]
        # source line for each phi: start A(phi), unit direction e(phi)
        A = np.stack([ra * cph, ra * sph, np.full_like(phi, za)], axis=1)          # (F,3)
        e = np.stack([dr * cph, dr * sph, np.full_like(phi, dz)], axis=1)          # (F,3)
        d = ob3[:, None, :] - A[None, :, :]                                        # (P,F,3)
        u0 = np.einsum('pfk,fk->pf', d, e)
        h2 = np.maximum(np.einsum('pfk,pfk->pf', d, d) - u0 ** 2, 1e-30)
        h = np.sqrt(h2)
        v1, v2 = -u0, L - u0
        asinh = np.arcsinh(v2 / h) - np.arcsinh(v1 / h)
        root = np.sqrt(v2 ** 2 + h2) - np.sqrt(v1 ** 2 + h2)
        si = 0.5 * L * (xi + 1.0)
        vv = si[None, None, :] - u0[:, :, None]
        R = np.sqrt(vv ** 2 + h2[:, :, None])
        rem = ((np.exp(-1j * K * R) - 1.0) / R) * wi * (0.5 * L)
        in0 = asinh + rem.sum(axis=2)                                   # (P,F)  weight 1
        in1 = (u0 * asinh + root) + (rem * si).sum(axis=2)               # weight s'
        for w, ker in ((0, wphi), (1, wphi * cph)):
            for b, inner in ((0, in0), (1, in1)):
                g = (inner @ ker).reshape(ns, -1)                        # (ns, no)
                M[w, 0, b, :, q] = (0.5 * pr.seg_len) * (g * wo).sum(axis=1)
                M[w, 1, b, :, q] = (0.5 * pr.seg_len) * (g * so * wo).sum(axis=1)
    return M


def _halves(pr: Profile, n: int):
    """Basis n sits on node n+1: rising over segment n, falling over n+1.
    [(segment, alpha, beta, sigma)] with current alpha + beta*s on the segment."""
    p, q = n, n + 1
    Lp, Lq = pr.seg_len[p], pr.seg_len[q]
    return [(p, 0.0, 1.0 / Lp, 1.0 / Lp), (q, 1.0, -1.0 / Lq, -1.0 / Lq)]


def impedance_matrix(pr: Profile) -> np.ndarray:
    M = _moments(pr)
    nb = pr.n_basis
    halves = [_halves(pr, n) for n in range(nb)]
    tr, tz = pr.seg_t[:, 0], pr.seg_t[:, 1]
    Z = np.zeros((nb, nb), dtype=complex)
    pref = 1j * ETA0 / (4.0 * math.pi * K)
    for m in range(nb):
        for n in range(m, nb):
            acc = 0.0 + 0.0j
            for p, ap, bp, sp in halves[m]:
                for q, aq, bq, sq in halves[n]:
                    def ff(w):
                        return (ap * aq * M[w, 0, 0, p, q] + ap * bq * M[w, 0, 1, p, q]
                                + bp * aq * M[w, 1, 0, p, q] + bp * bq * M[w, 1, 1, p, q])
                    acc += K ** 2 * (tr[p] * tr[q] * ff(1) + tz[p] * tz[q] * ff(0)) - sp * sq * M[0, 0, 0, p, q]
            Z[m, n] = Z[n, m] = pref * acc
    return Z


def _gap_vector(pr: Profile) -> np.ndarray:
    """V_m = (1/g) int_gap f_m dt for a 1 V gap: the basis functions' overlap
    with the impressed field, and also the weights of the gap-averaged current."""
    t0, t1 = pr.gap
    g = t1 - t0
    V = np.zeros(pr.n_basis)
    for n in range(pr.n_basis):
        for p, a, b, _ in _halves(pr, n):
            s0, s1 = pr.t_node[p], pr.t_node[p + 1]
            lo, hi = max(s0, t0), min(s1, t1)
            if hi > lo:
                x0, x1 = lo - s0, hi - s0
                V[n] += (a * (x1 - x0) + 0.5 * b * (x1 ** 2 - x0 ** 2)) / g
    return V


@dataclass(frozen=True)
class BorSolution:
    profile: Profile
    currents: np.ndarray        # node currents I_n at the basis nodes [A for 1 V]
    weights: np.ndarray

    @property
    def input_impedance(self) -> complex:
        return complex(1.0 / (self.weights @ self.currents))

    @property
    def circuit_power(self) -> float:
        return 0.5 * float(np.real(np.conj(self.weights @ self.currents)))


def solve(pr: Profile) -> BorSolution:
    V = _gap_vector(pr)
    return BorSolution(pr, np.linalg.solve(impedance_matrix(pr), V.astype(complex)), V)


def far_field(sol: BorSolution, theta: np.ndarray) -> np.ndarray:
    """E_theta * r * exp(jkr) [V] against theta, from the ring currents:
    A_theta picks up j*J1(k rho sin theta) from the radial part of the current
    and -sin(theta)*J0(...) from the axial part."""
    from scipy.special import j0, j1
    pr = sol.profile
    x, w = np.polynomial.legendre.leggauss(12)
    I = np.concatenate([[0.0], sol.currents, [0.0]])
    th = np.asarray(theta, dtype=float)[:, None]
    out = np.zeros(th.shape[0], dtype=complex)
    for p in range(pr.n_seg):
        s = 0.5 * pr.seg_len[p] * (x + 1.0)
        cur = I[p] + (I[p + 1] - I[p]) * s / pr.seg_len[p]
        rho = pr.nodes[p, 0] + pr.seg_t[p, 0] * s
        z = pr.nodes[p, 1] + pr.seg_t[p, 1] * s
        arg = K * rho[None, :] * np.sin(th)
        comp = (pr.seg_t[p, 0] * np.cos(th) * 1j * j1(arg) - pr.seg_t[p, 1] * np.sin(th) * j0(arg))
        out += (cur * comp * np.exp(1j * K * z[None, :] * np.cos(th)) * w).sum(axis=1) * 0.5 * pr.seg_len[p]
    return -1j * ETA0 * K / (4.0 * math.pi) * out


def radiated_power(sol: BorSolution, n: int = 721) -> float:
    th = np.linspace(0.0, math.pi, n)
    e = np.abs(far_field(sol, th)) ** 2
    return float(np.trapezoid(e * np.sin(th), th) * 2.0 * math.pi / (2.0 * ETA0))


def directivity(sol: BorSolution, n: int = 721) -> tuple[float, float]:
    """(peak directivity, elevation of the peak in degrees above the plane z = 0).
    The pattern is omnidirectional in azimuth, so the peak is a cone of
    directions; positive elevation is towards +z."""
    th = np.linspace(0.0, math.pi, n)
    u = np.abs(far_field(sol, th)) ** 2 / (2.0 * ETA0)
    p = float(np.trapezoid(u * np.sin(th), th) * 2.0 * math.pi)
    k = int(np.argmax(u))
    return 4.0 * math.pi * float(u[k]) / p, 90.0 - math.degrees(th[k])
