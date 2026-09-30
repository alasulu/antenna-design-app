"""Thin-wire method of moments: an independent full-wave reference.

Why this exists. Every archetype in the catalogue is a closed form, and a
closed form can only be checked against another closed form or against
published data. Several archetypes have neither - a one-wavelength loop, a quad
loop, a folded dipole's transformation ratio - and have been carrying
engineering constants labelled "indicative" because there was no way to do
better without a commercial solver. This module is that way.

Formulation. Electric-field integral equation on thin wires, mixed-potential
form, triangular (rooftop) basis with Galerkin testing:

    Z_mn = (j*eta/(4*pi*k)) * [ k^2 * sum_pq (t_p . t_q) <f_m, f_n>_pq
                                - sum_pq sigma_p sigma_q <1, 1>_pq ]

where <w, w'>_pq is the double line integral of w(s) w'(s') exp(-jkR)/R over
segments p and q, and sigma is the basis function's divergence, constant on
each half. The reduced (thin-wire) kernel R = sqrt(|r - r'|^2 + a^2) keeps the
integrand finite; the singular part is then removed analytically rather than
hoped away by quadrature:

    integral w(s')/R ds'  is exact for linear w, in asinh and sqrt,

and only the smooth remainder (exp(-jkR) - 1)/R goes to Gauss-Legendre. That
matters: with a = 0.001 lambda and segments of 0.02 lambda the kernel peak is
twenty times narrower than a segment, and plain quadrature silently loses the
self term.

Units. Geometry is in WAVELENGTHS and k = 2*pi throughout, so a solution is
frequency-independent; impedances come out in ohms. A structure given in metres
is scaled by the caller.

What it does not do. No ground plane (image theory stands in, exactly, for an
infinite PEC plane), no dielectric, no distributed loss, and the reduced kernel
by default (the exact one is opt-in, for straight runs - see below). Junctions
of any number of wires ARE supported, but
only at wire end points. Each of
those is a real limit, not an approximation that washes out, so the validation
suite pins the cases where they do not bite.

A limit of the delta gap specifically, measured rather than assumed: on an
electrically SHORT wire the feed node carries a local current excess over the
smooth distribution - about 8% at 0.1 lambda and 17% at 0.04 lambda - and it
grows as the mesh resolves it. Reactance barely notices (about 0.6% per mesh
doubling), but the resistance referred to that node reads roughly 10% low and
keeps drifting. For short-antenna RESISTANCE, trust the closed form over this
solver; for reactance, impedance of anything near resonance, patterns and
power, the delta gap is fine.

A limit of the reduced kernel, also measured: keep every segment at least about
three wire radii long. On a 0.006-wavelength wire the resonant resistance is
72-74 ohm while segments are 3 to 8 radii long, drifts to 79 ohm at 1.3 radii,
and reads 115 ohm - with the resonance moved 6% - at 0.8 radii. A fine mesh on a
fat wire is not more accurate; it is outside the approximation. WireModel warns
when any segment is shorter than its own radius, which is the clearly broken
end of that range.

`exact=True` lifts that limit on straight runs. It puts the exact kernel -
current and observer both spread round the circumference - on every self and
near collinear pair, and the same 0.006-wavelength wire then holds its
resonance from 3.7 radii per segment down to 0.8 (reactance within 0.6 ohm at a
fixed length, where the reduced kernel falls 16 ohm), its resistance creeping 3%
as the delta gap sharpens. Bends and junctions keep the reduced kernel. It is
off by default so every result derived with the reduced kernel reproduces bit
for bit; on thin wire the two differ by about 0.35%. What it does NOT cure is
the delta gap itself on a very fat wire: at a = 0.015 lambda neither kernel
converges, and that needs a finite-gap feed model.

`solve(..., gap=segments)` is that model, the one `bor` uses for tubes: 1 V
impressed uniformly along a run of segments of physical length g, and the input
current taken as the gap-averaged current, so that V conj(I)/2 is the power the
source delivers. `gapped_dipole` builds a straight dipole whose gap is exactly
two segments. A delta gap is its g -> 0 limit; a gap of fixed width converges as
the mesh is refined, where the delta gap's own capacitance keeps growing.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from ..core.constants import ETA0

__all__ = [
    "Wire", "WireModel", "MoMSolution", "dipole", "loop", "arc", "halo",
    "folded_dipole_wire", "helix_over_ground", "lpda_model", "rhombic_model",
    "top_hat_monopole", "bicone_cage",
    "tl_admittance", "NetworkSolution", "solve_network",
    "solve", "input_impedance", "far_field", "directivity", "radiated_power",
    "gap_vector", "gapped_dipole",
]

K = 2.0 * math.pi          # wavenumber, geometry being in wavelengths
_GL_OUTER = np.polynomial.legendre.leggauss(12)
_GL_INNER = np.polynomial.legendre.leggauss(24)


@dataclass(frozen=True)
class Wire:
    """A polyline of PEC wire. `nodes` is (M, 3) in wavelengths."""

    nodes: np.ndarray
    radius: float
    closed: bool = False

    def __post_init__(self) -> None:
        arr = np.asarray(self.nodes, dtype=float)
        if arr.ndim != 2 or arr.shape[1] != 3 or arr.shape[0] < 2:
            raise ValueError("nodes must be (M, 3) with M >= 2")
        if self.radius <= 0:
            raise ValueError("radius must be positive")
        object.__setattr__(self, "nodes", arr)


@dataclass
class WireModel:
    """Wires discretised into segments and rooftop basis functions.

    Basis function n is centred on a node and spans two segments, rising to 1
    at the node. Each of its two halves is stored with its orientation, because
    at a junction a wire may END at the node or START there, and the current
    has to flow through regardless:

      - half 0 carries current TOWARDS the node, half 1 AWAY from it;
      - `basis_at_b[n, h]` says whether the node is the segment's end point
        (its `b`) or its start point (its `a`).

    An open wire of N segments carries N-1 interior functions, which forces the
    current to zero at a FREE end; a closed wire carries N.

    Junctions. Wherever the ends of K >= 2 open wires coincide, K-1 junction
    functions are added, each carrying current from the first of those wires
    (the lowest wire index) through the node into one of the others. That is
    Kirchhoff's current law by construction: the reference wire carries the sum
    of the rest. Two wires meeting end to end is simply a bend, and reproduces
    one polyline exactly - which is how this was checked. Junctions are found
    only at wire END points; a wire that should branch part-way along must be
    split there.
    """

    wires: list[Wire]
    seg_a: np.ndarray = field(init=False)      # segment start points
    seg_b: np.ndarray = field(init=False)      # segment end points
    seg_t: np.ndarray = field(init=False)      # unit tangents
    seg_len: np.ndarray = field(init=False)
    seg_rad: np.ndarray = field(init=False)
    basis: np.ndarray = field(init=False)      # (Nb, 2) segment indices (towards, away)
    basis_at_b: np.ndarray = field(init=False)  # (Nb, 2) node is the segment's end?
    basis_node: np.ndarray = field(init=False)  # (Nb, 3) node positions
    junctions: list = field(init=False)        # [(point, (basis indices,...)), ...]

    def __post_init__(self) -> None:
        a, b, rad = [], [], []
        basis, at_b, node = [], [], []
        ends = []            # (wire index, point, segment index, node at b?)
        for wi, w in enumerate(self.wires):
            base = len(a)
            nodes = w.nodes
            pairs = list(zip(nodes[:-1], nodes[1:]))
            if w.closed:
                pairs.append((nodes[-1], nodes[0]))
            for p, q in pairs:
                a.append(p)
                b.append(q)
                rad.append(w.radius)
            n = len(pairs)
            if w.closed:
                for i in range(n):
                    basis.append((base + i, base + (i + 1) % n))
                    at_b.append((True, False))
                    node.append(pairs[i][1])
            else:
                for i in range(n - 1):
                    basis.append((base + i, base + i + 1))
                    at_b.append((True, False))
                    node.append(pairs[i][1])
                ends.append((wi, np.asarray(pairs[0][0], float), base, False))
                ends.append((wi, np.asarray(pairs[-1][1], float), base + n - 1, True))
        self.seg_a = np.array(a, dtype=float)
        self.seg_b = np.array(b, dtype=float)
        d = self.seg_b - self.seg_a
        self.seg_len = np.linalg.norm(d, axis=1)
        if np.any(self.seg_len <= 0):
            raise ValueError("a segment has zero length")
        self.seg_t = d / self.seg_len[:, None]
        self.seg_rad = np.array(rad, dtype=float)
        short = self.seg_len < self.seg_rad
        if np.any(short):
            import warnings
            warnings.warn(
                f"{int(short.sum())} segment(s) shorter than their wire radius: the "
                f"reduced thin-wire kernel is not valid there and results can be "
                f"badly wrong. Keep segments at least ~3 radii long, or solve a "
                f"straight wire with exact=True.",
                UserWarning, stacklevel=2)

        # group coincident open-wire ends into junctions
        tol = 1e-6 * float(self.seg_len.min())
        used = [False] * len(ends)
        self.junctions = []
        for i, (wi, pi, si, bi) in enumerate(ends):
            if used[i]:
                continue
            group = [i]
            for j in range(i + 1, len(ends)):
                if not used[j] and np.linalg.norm(ends[j][1] - pi) <= tol:
                    group.append(j)
            if len(group) < 2:
                continue
            for j in group:
                used[j] = True
            group.sort(key=lambda k: ends[k][0])
            ref = ends[group[0]]
            made = []
            for k in group[1:]:
                other = ends[k]
                basis.append((ref[2], other[2]))
                at_b.append((ref[3], other[3]))
                node.append(pi)
                made.append(len(basis) - 1)
            self.junctions.append((pi.copy(), tuple(made)))

        self.basis = np.array(basis, dtype=int).reshape(-1, 2)
        self.basis_at_b = np.array(at_b, dtype=bool).reshape(-1, 2)
        self.basis_node = np.array(node, dtype=float).reshape(-1, 3)
        if len(self.basis) == 0:
            raise ValueError("model has no basis functions; wires are too coarse")

    @property
    def n_basis(self) -> int:
        return len(self.basis)

    def node_of(self, n: int) -> np.ndarray:
        """Position of the node basis function `n` is centred on."""
        return self.basis_node[n]

    def junction_at(self, point) -> tuple[int, ...]:
        """The junction functions at `point`, for feeding there: a delta gap on
        the reference wire (the lowest-indexed wire meeting at the junction)
        drives all of them together. Pass the tuple to `solve` as the feed."""
        p = np.asarray(point, float)
        for q, idx in self.junctions:
            if np.linalg.norm(q - p) <= 1e-6 * float(self.seg_len.min()) + 1e-12:
                return idx
        raise ValueError(f"no junction at {tuple(p)}")


# ---------------------------------------------------------------- kernels

_GL_PHI = np.polynomial.legendre.leggauss(24)


def _inner_exact(u0: np.ndarray, L: float, a: float):
    """Inner integrals over a straight source segment with the EXACT kernel.

    The reduced kernel puts the source current on the axis and the observer
    on the surface, R = sqrt(v^2 + a^2). The exact kernel spreads both round
    the circumference, and after one integration by symmetry that is
    R(phi) = sqrt(v^2 + 4 a^2 sin^2(phi/2)) averaged over phi. It only differs
    from the reduced one where the axial separation is a few radii or less -
    the self and near terms of a fat wire - and there it is the difference
    between a converged answer and one that runs away as the mesh is refined.

    The phi average has a log singularity at phi = 0 when observer and source
    coincide; phi = pi t^2 turns it into t log t, which Gauss-Legendre handles.
    `u0` are the observers' axial positions measured from the segment's start.
    Returns (integral of G, integral of s' G), each averaged over phi.
    """
    xt, wt = _GL_PHI
    tt = 0.5 * (xt + 1.0)
    wphi = wt * tt                   # (1/pi) int_0^pi d phi -> int_0^1 2t dt, GL on [0,1]
    phi = math.pi * tt ** 2
    rho2 = (2.0 * a * np.sin(0.5 * phi)) ** 2                         # (F,)
    rho = np.sqrt(rho2)
    xi, wi = _GL_INNER
    si = 0.5 * L * (xi + 1.0)
    v1 = -u0[:, None]                                                 # (P,1)
    v2 = L - u0[:, None]
    asinh = np.arcsinh(v2 / rho[None, :]) - np.arcsinh(v1 / rho[None, :])      # (P,F)
    root = np.sqrt(v2 ** 2 + rho2[None, :]) - np.sqrt(v1 ** 2 + rho2[None, :])
    v = si[None, None, :] - u0[:, None, None]                         # (P,1,G)
    R = np.sqrt(v ** 2 + rho2[None, :, None])                         # (P,F,G)
    rem = ((np.exp(-1j * K * R) - 1.0) / R) * wi[None, None, :] * (0.5 * L)
    in0 = asinh + rem.sum(axis=2)
    in1 = (u0[:, None] * asinh + root) + (rem * si[None, None, :]).sum(axis=2)
    return in0 @ wphi, in1 @ wphi


def _segment_moments(model: WireModel, exact: bool = False):
    """M[a][b][p, q] = double integral of s^a s'^b exp(-jkR)/R over segments p, q.

    Every basis function's weight on a segment is linear in the arc length, so
    these four arrays are all the geometry any matrix element needs, and each
    segment pair is touched once instead of up to four times. The inner
    integral keeps the exact static part; only the smooth remainder is
    quadratured.
    """
    ns = len(model.seg_len)
    xo, wo = _GL_OUTER
    xi, wi = _GL_INNER

    # observation points: (ns, no, 3), one set per outer quadrature node
    so = 0.5 * model.seg_len[:, None] * (xo[None, :] + 1.0)              # (ns, no)
    obs = (model.seg_a[:, None, :] + so[:, :, None] * model.seg_t[:, None, :])
    flat = obs.reshape(-1, 3)

    M = np.zeros((2, 2, ns, ns), dtype=complex)
    # The reduced kernel's radius for a PAIR of segments is the rms of the
    # two, (a_p^2 + a_q^2)/2: symmetric, so wires of different radii give the
    # same impedance in either order (only the upper triangle is assembled),
    # and exactly a^2 - to the bit - when the radii are equal.
    rad2_obs = np.repeat(model.seg_rad ** 2, len(xo))
    for q in range(ns):
        A, u, L = model.seg_a[q], model.seg_t[q], model.seg_len[q]
        a2 = 0.5 * (rad2_obs + model.seg_rad[q] ** 2)
        d = flat - A
        u0 = d @ u
        rho2 = np.maximum(np.einsum('ij,ij->i', d, d) - u0 ** 2 + a2, a2)
        rho = np.sqrt(rho2)
        v1, v2 = -u0, L - u0
        asinh = np.arcsinh(v2 / rho) - np.arcsinh(v1 / rho)
        root = np.sqrt(v2 ** 2 + rho2) - np.sqrt(v1 ** 2 + rho2)

        si = 0.5 * L * (xi + 1.0)                                        # (ni,)
        v = si[None, :] - u0[:, None]
        R = np.sqrt(v ** 2 + rho2[:, None])
        rem = ((np.exp(-1j * K * R) - 1.0) / R) * wi[None, :] * (0.5 * L)

        inner0 = asinh + rem.sum(axis=1)                        # weight 1
        inner1 = (u0 * asinh + root) + (rem * si[None, :]).sum(axis=1)   # weight s'
        if exact:
            aq = model.seg_rad[q]
            mid_q = A + 0.5 * L * u
            for p in range(ns):
                if abs(float(model.seg_t[p] @ u)) < 1.0 - 1e-9:
                    continue
                mid_p = 0.5 * (model.seg_a[p] + model.seg_b[p])
                off = mid_p - mid_q
                if np.linalg.norm(off - (off @ u) * u) > 1e-9 * L:
                    continue                      # parallel but not on the same axis
                if abs(off @ u) > 0.5 * (L + model.seg_len[p]) + 6.0 * aq:
                    continue                      # far enough for the reduced kernel
                rows = slice(p * len(xo), (p + 1) * len(xo))
                e0, e1 = _inner_exact(u0[rows], L,
                                      math.sqrt(0.5 * (aq ** 2 + model.seg_rad[p] ** 2)))
                inner0 = inner0.copy()
                inner1 = inner1.copy()
                inner0[rows], inner1[rows] = e0, e1
        for b, inner in ((0, inner0), (1, inner1)):
            g = inner.reshape(ns, -1)                                    # (ns, no)
            M[0, b, :, q] = (0.5 * model.seg_len) * (g * wo[None, :]).sum(axis=1)
            M[1, b, :, q] = (0.5 * model.seg_len) * (g * so * wo[None, :]).sum(axis=1)
    return M


def _half_weights(model: WireModel, n: int):
    """For basis n: [(segment, (alpha, beta), sigma)] for its two halves.

    The current on a half is (alpha + beta*s) along the segment's own tangent,
    s measured from the segment's start - so a half whose current runs against
    the tangent simply carries negative weights - and sigma is its divergence.
    Towards-the-node halves always have sigma = +1/L and away halves -1/L,
    whichever way the segment happens to be parameterised: divergence is a
    scalar and does not care.
    """
    out = []
    for h in (0, 1):
        p = int(model.basis[n, h])
        L = float(model.seg_len[p])
        at_b = bool(model.basis_at_b[n, h])
        if h == 0:                       # towards the node
            w = (0.0, 1.0 / L) if at_b else (-1.0, 1.0 / L)
            sigma = 1.0 / L
        else:                            # away from the node
            w = (1.0, -1.0 / L) if not at_b else (0.0, -1.0 / L)
            sigma = -1.0 / L
        out.append((p, w, sigma))
    return out


def impedance_matrix(model: WireModel, exact: bool = False) -> np.ndarray:
    """The Galerkin EFIE matrix [ohm].

    `exact` switches the self and near collinear terms to the exact kernel,
    which is what a FAT wire needs; for thin wire it changes nothing that
    matters, and it is off by default so that every result already derived
    with the reduced kernel reproduces exactly.
    """
    M = _segment_moments(model, exact)
    halves = [_half_weights(model, n) for n in range(model.n_basis)]
    dot = model.seg_t @ model.seg_t.T
    nb = model.n_basis
    Z = np.zeros((nb, nb), dtype=complex)
    pref = 1j * ETA0 / (4.0 * math.pi * K)
    for m in range(nb):
        for n in range(m, nb):
            acc = 0.0 + 0.0j
            for p, (ap, bp), sp in halves[m]:
                for q, (aq, bq), sq in halves[n]:
                    ff = (ap * aq * M[0, 0, p, q] + ap * bq * M[0, 1, p, q]
                          + bp * aq * M[1, 0, p, q] + bp * bq * M[1, 1, p, q])
                    acc += K ** 2 * dot[p, q] * ff - sp * sq * M[0, 0, p, q]
            Z[m, n] = Z[n, m] = pref * acc
    return Z


# ---------------------------------------------------------------- solution

@dataclass(frozen=True)
class MoMSolution:
    model: WireModel
    currents: np.ndarray
    feed: int | tuple[int, ...]
    gap_weights: np.ndarray | None = None

    @property
    def feed_current(self) -> complex:
        """Terminal current. A junction feed drives several functions that all
        share the reference wire, and its terminal current is their sum. A finite
        gap's is the current averaged over the gap."""
        if self.gap_weights is not None:
            return complex(self.gap_weights @ self.currents)
        if isinstance(self.feed, tuple):
            return complex(sum(self.currents[k] for k in self.feed))
        return complex(self.currents[self.feed])

    @property
    def input_impedance(self) -> complex:
        return 1.0 / self.feed_current

    @property
    def circuit_power(self) -> float:
        """P = 0.5 Re(V I*) for a 1 V gap, delta or finite [W]."""
        return 0.5 * float(np.real(np.conj(self.feed_current)))


def solve(model: WireModel, feed: int | tuple[int, ...] | None = None,
          loads: dict[int, complex] | None = None,
          exact: bool = False, gap=None) -> MoMSolution:
    """Delta-gap excitation of one basis function; 1 V across the gap.

    `gap`, a run of segment indices, replaces the delta gap with a FINITE one:
    see `gap_vector`. `feed` is then ignored.

    `loads` puts a SERIES lumped impedance at a basis function, the way a real
    wire code does it: the load's voltage drop is Z_L * I_m, so Z_L adds
    straight onto the diagonal. A capacitor tuning a halo's gap, a loading coil
    partway up a whip, a resistive termination on a travelling-wave wire are
    all this one line. A very large load is the same as cutting the wire there,
    which is how it gets checked.
    """
    nb = model.n_basis
    if feed is None:
        feed = nb // 2
    feeds = feed if isinstance(feed, tuple) else (feed,)
    for f in feeds:
        if not 0 <= f < nb:
            raise ValueError(f"feed index {f} outside 0..{nb - 1}")
    Z = impedance_matrix(model, exact)
    for m, zl in (loads or {}).items():
        if not 0 <= m < nb:
            raise ValueError(f"load index {m} outside 0..{nb - 1}")
        Z[m, m] += zl
    if gap is not None:
        w = gap_vector(model, gap)
        return MoMSolution(model, np.linalg.solve(Z, w.astype(complex)), int(np.argmax(w)), w)
    V = np.zeros(nb, dtype=complex)
    for f in feeds:
        V[f] = 1.0
    return MoMSolution(model, np.linalg.solve(Z, V), feed)


def gap_vector(model: WireModel, segments) -> np.ndarray:
    """A finite gap: 1 V impressed uniformly, V/g, along the tangents of `segments`
    (consecutive segments of total length g, running the same way).

    V_m = (1/g) int_gap f_m ds - each basis function's overlap with the impressed
    field, and also its weight in the gap-averaged input current, so the power
    V conj(I_in)/2 is what the source delivers. With the gap one segment either
    side of a node, the node's function gets 1/2 and its neighbours 1/4 each.
    """
    segs = {int(p) for p in segments}
    g = float(sum(model.seg_len[p] for p in segs))
    w = np.zeros(model.n_basis)
    for n in range(model.n_basis):
        for p, (alpha, beta), _ in _half_weights(model, n):
            if p in segs:
                L = float(model.seg_len[p])
                w[n] += alpha * L + beta * L * L / 2.0
    return w / g


def input_impedance(model: WireModel, feed: int | None = None,
                    loads: dict[int, complex] | None = None,
                    exact: bool = False) -> complex:
    return solve(model, feed, loads, exact).input_impedance


# ---------------------------------------------------------------- radiation

def _segment_phasor(A, u, L, alpha, beta, rhat):
    """integral_0^L (alpha + beta*s) exp(j k rhat.r(s)) ds, in closed form.

    `rhat` is (P, 3); the result is (P,). The g -> 0 branch is the broadside
    limit of a segment, which is not a rare case - it is every segment
    perpendicular to the observation direction, i.e. the whole of a dipole's
    main beam - so it is taken by np.where rather than left to divide by zero.
    """
    g = K * (rhat @ u)
    phase = np.exp(1j * K * (rhat @ A))
    # Near broadside the closed form subtracts two terms of size L/(g L)^2
    # and loses every digit (60% wrong a 1e-8 rad from broadside); below
    # |g L| = 0.1 take the series, sum_n (jg)^n L^(n+1+b) / (n! (n+1+b)),
    # which ten terms carry to 1e-16.
    small = np.abs(g * L) < 0.1
    gs = np.where(small, 1.0, g)
    e = np.exp(1j * gs * L)
    i0 = (e - 1.0) / (1j * gs)
    i1 = (L * e) / (1j * gs) - (e - 1.0) / (1j * gs) ** 2
    if np.any(small):
        x = 1j * np.where(small, g, 0.0) * L
        term = np.ones_like(x)
        s0 = np.zeros_like(x)
        s1 = np.zeros_like(x)
        for n in range(10):
            s0 = s0 + term / (n + 1)
            s1 = s1 + term / (n + 2)
            term = term * x / (n + 1)
        i0 = np.where(small, L * s0, i0)
        i1 = np.where(small, L * L * s1, i1)
    return phase * (alpha * i0 + beta * i1)


def far_field(sol: MoMSolution, theta: np.ndarray, phi: np.ndarray):
    """(E_theta, E_phi) up to a common constant, on the given angle arrays."""
    th = np.asarray(theta, dtype=float)
    ph = np.asarray(phi, dtype=float)
    st, ct, sp, cp = np.sin(th), np.cos(th), np.sin(ph), np.cos(ph)
    rhat = np.stack([st * cp, st * sp, ct], axis=-1)
    flat = rhat.reshape(-1, 3)
    model = sol.model

    # accumulate per segment: each one carries a linear current weight summed
    # over the one or two basis functions it belongs to
    weights = np.zeros((len(model.seg_len), 2), dtype=complex)
    for n in range(model.n_basis):
        In = sol.currents[n]
        if In == 0:
            continue
        for p, (alpha, beta), _ in _half_weights(model, n):
            weights[p, 0] += In * alpha
            weights[p, 1] += In * beta

    N = np.zeros(flat.shape, dtype=complex)
    for p, (alpha, beta) in enumerate(weights):
        if alpha == 0 and beta == 0:
            continue
        vals = _segment_phasor(model.seg_a[p], model.seg_t[p], model.seg_len[p],
                               alpha, beta, flat)
        N += vals[:, None] * model.seg_t[p][None, :]
    N = N.reshape(rhat.shape)
    e_th = (N[..., 0] * cp + N[..., 1] * sp) * ct - N[..., 2] * st
    e_ph = -N[..., 0] * sp + N[..., 1] * cp
    return e_th, e_ph


def _hemisphere_grid(nth: int, nph: int):
    ct, wct = np.polynomial.legendre.leggauss(nth)
    th = np.arccos(ct)
    ph = np.arange(nph) * 2.0 * math.pi / nph
    return th[:, None] + 0.0 * ph[None, :], 0.0 * th[:, None] + ph[None, :], wct


def pattern_power(sol: MoMSolution, nth: int = 60, nph: int = 60):
    """(U on a Gauss-Legendre sphere, total power, quadrature weights)."""
    TH, PH, wct = _hemisphere_grid(nth, nph)
    e_th, e_ph = far_field(sol, TH, PH)
    U = np.abs(e_th) ** 2 + np.abs(e_ph) ** 2
    total = float(wct @ U.sum(axis=1)) * (2.0 * math.pi / nph)
    return U, total, wct


def radiated_power(sol: MoMSolution, nth: int = 60, nph: int = 60) -> float:
    """Radiated power [W] from the far field, scaled to match the 1 V source.

    The scale factor is fixed by |E|^2 = (eta k^2 / (32 pi^2 r^2)) |N|^2.
    """
    _, total, _ = pattern_power(sol, nth, nph)
    return total * ETA0 * K ** 2 / (32.0 * math.pi ** 2)


def directivity(sol: MoMSolution, nth: int = 60, nph: int = 60) -> float:
    U, total, _ = pattern_power(sol, nth, nph)
    return 4.0 * math.pi * float(U.max()) / total


def directivity_towards(sol: MoMSolution, theta: float, phi: float,
                        nth: int = 60, nph: int = 60) -> float:
    _, total, _ = pattern_power(sol, nth, nph)
    e_th, e_ph = far_field(sol, np.array(theta), np.array(phi))
    u = float(abs(e_th) ** 2 + abs(e_ph) ** 2)
    return 4.0 * math.pi * u / total


# ---------------------------------------------------------------- geometries

def dipole(length: float, radius: float = 1e-3, segments: int = 40) -> WireModel:
    """Straight z-directed dipole centred on the origin.

    `segments` is rounded UP to an even number. With N segments there are N+1
    nodes, so a node sits at z = 0 only when N is even - and the feed has to be
    at a node, because that is where the basis functions are centred. An odd
    count puts the source half a segment off centre and quietly breaks the
    symmetry of the problem.
    """
    if segments % 2:
        segments += 1
    z = np.linspace(-0.5 * length, 0.5 * length, segments + 1)
    nodes = np.stack([np.zeros_like(z), np.zeros_like(z), z], axis=1)
    return WireModel([Wire(nodes, radius)])


def gapped_dipole(length: float, radius: float, gap: float, seg: float) -> tuple[WireModel, tuple[int, int]]:
    """Straight z-directed dipole with a centred feed gap of physical length `gap`,
    cut into exactly two segments; each arm cut into segments no longer than `seg`.
    Returns (model, gap segments) for `solve(model, gap=...)`. Segments shorter than
    the radius need `exact=True`, which straight wires take."""
    arm = 0.5 * (length - gap)
    n = max(1, int(math.ceil(arm / seg)))
    z = np.concatenate([np.linspace(-0.5 * length, -0.5 * gap, n + 1), [0.0],
                        np.linspace(0.5 * gap, 0.5 * length, n + 1)])
    nodes = np.stack([np.zeros_like(z), np.zeros_like(z), z], axis=1)
    return WireModel([Wire(nodes, radius)]), (n, n + 1)


def loop(circumference: float, radius: float = 1e-3, segments: int = 48,
         plane: str = "xy") -> WireModel:
    """Closed circular loop of the given circumference, centred on the origin."""
    b = circumference / (2.0 * math.pi)
    t = np.arange(segments) * 2.0 * math.pi / segments
    x, y = b * np.cos(t), b * np.sin(t)
    zero = np.zeros_like(x)
    nodes = {"xy": np.stack([x, y, zero], axis=1),
             "xz": np.stack([x, zero, y], axis=1),
             "yz": np.stack([zero, x, y], axis=1)}[plane]
    return WireModel([Wire(nodes, radius, closed=True)])


def folded_dipole_wire(length: float, spacing: float, radius: float = 1e-3,
                       segments: int = 40) -> WireModel:
    """Two parallel conductors shorted at both ends: one closed rectangle."""
    n = max(segments // 2, 4)
    xs = np.linspace(-0.5 * length, 0.5 * length, n + 1)
    top = np.stack([xs, np.full_like(xs, 0.5 * spacing), np.zeros_like(xs)], axis=1)
    bot = np.stack([xs[::-1], np.full_like(xs, -0.5 * spacing),
                    np.zeros_like(xs)], axis=1)
    nodes = np.vstack([top, bot])
    return WireModel([Wire(nodes, radius, closed=True)])


def arc(length: float, subtend: float, radius: float = 1e-3,
        segments: int = 60) -> WireModel:
    """An open wire of the given LENGTH bent into a circular arc, in the x-y
    plane, symmetric about +x, with its midpoint towards +x.

    `subtend` is the angle the arc covers, so 0 is a straight wire along y and
    just under 2*pi is a ring with a small gap. Bending a half-wave dipole round
    until its ends nearly meet is exactly what a halo is, and running the bend
    continuously from zero is how the bent-wire path gets checked: at subtend 0
    it has to reproduce `dipole`, which is already validated two ways.
    """
    if subtend < 1e-6:
        y = np.linspace(-0.5 * length, 0.5 * length, segments + 1)
        nodes = np.stack([np.zeros_like(y), y, np.zeros_like(y)], axis=1)
        return WireModel([Wire(nodes, radius)])
    curve = length / subtend
    t = np.linspace(-0.5 * subtend, 0.5 * subtend, segments + 1)
    nodes = np.stack([curve * np.cos(t), curve * np.sin(t),
                      np.zeros_like(t)], axis=1)
    return WireModel([Wire(nodes, radius)])


def halo(circumference: float, gap: float, radius: float = 1e-3,
         segments: int = 60) -> WireModel:
    """A ring of the given circumference broken by a gap of arc length `gap`.

    The conductor is circumference - gap long, which is the identity the
    `halo_loop` spec got wrong: it subtracted the gap twice.
    """
    if not 0 < gap < circumference:
        raise ValueError("gap must be a fraction of the circumference")
    return arc(circumference - gap,
               2.0 * math.pi * (1.0 - gap / circumference), radius, segments)


# ------------------------------------------------------- frequency behaviour
#
# Geometry here is in wavelengths, so scaling every dimension by `s` IS moving
# the frequency by `s`. A frequency sweep therefore costs nothing beyond
# rebuilding the model, and the functions below all take a callable
# `z_of_scale(s) -> complex` rather than a model, so they work for any
# structure the caller can parameterise.


def vswr(z: complex, z0: float) -> float:
    """Standing-wave ratio of `z` against a real reference `z0`."""
    g = abs((z - z0) / (z + z0))
    return float("inf") if g >= 1.0 else (1.0 + g) / (1.0 - g)


def resonant_scale(z_of_scale, lo: float = 0.85, hi: float = 1.15,
                   steps: int = 26) -> float | None:
    """Frequency scale where the reactance crosses zero, or None if it does not
    inside the bracket. Returning None rather than a bracket end matters: a
    structure with no resonance in range should say so, not report an edge."""
    below = z_of_scale(lo).imag < 0
    if (z_of_scale(hi).imag < 0) == below:
        return None
    # keep the side whose sign `lo` started with: a reactance that crosses
    # DOWNWARD (an antiresonance) used to walk to the bracket's end
    for _ in range(steps):
        mid = 0.5 * (lo + hi)
        if (z_of_scale(mid).imag < 0) == below:
            lo = mid
        else:
            hi = mid
    return lo


def antenna_q(z_of_scale, scale0: float, h: float = 2e-3) -> float:
    """Q = (w0 / 2 R0) |dZ/dw|, after Yaghjian and Best.

    Scale is frequency, so d/d(scale) at scale0 = 1 already carries the w0.
    This is the narrowband route to bandwidth; `vswr_bandwidth` is the direct
    one, and the two agreeing is the check that either is worth anything.
    """
    z0 = z_of_scale(scale0)
    dz = (z_of_scale(scale0 + h) - z_of_scale(scale0 - h)) / (2.0 * h)
    return float(abs(dz) * scale0 / (2.0 * z0.real))


def vswr_bandwidth(z_of_scale, scale0: float, z0: float, target: float = 2.0,
                   step: float = 5e-3, span: float = 0.6):
    """(fractional bandwidth, (lower edge, upper edge)) about `scale0`.

    Walks outwards in `step` until VSWR exceeds `target`, then bisects. The
    walk matters: bisecting a bracket chosen by guesswork finds whichever
    crossing happens to be inside it, and a wideband antenna has several.
    """
    edges = []
    for direction in (-1.0, +1.0):
        lo = hi = scale0
        while abs(hi - scale0) < span:
            hi += direction * step
            if vswr(z_of_scale(hi), z0) > target:
                break
        else:
            return float("inf"), (None, None)
        for _ in range(30):
            mid = 0.5 * (lo + hi)
            if vswr(z_of_scale(mid), z0) < target:
                lo = mid
            else:
                hi = mid
        edges.append(hi)
    return (edges[1] - edges[0]) / scale0, (edges[0], edges[1])


def helix_over_ground(circumference: float = 1.0, pitch_deg: float = 13.0,
                      turns: int = 10, radius: float = 0.005,
                      feed_height: float = 0.02, seg_per_turn: int = 16):
    """Axial-mode helix on an infinite PEC ground plane, by image theory.

    Returns (model, feed index). The model is ONE continuous open wire: the
    image helix below the plane, a short vertical feed wire straddling z = 0,
    and the real helix above. Current continuity through the bends produces
    exactly the image-current rules - horizontal components reversed, vertical
    ones not - so nothing special is needed, and the far field of the result is
    mirror-symmetric about z = 0 to about 1e-7, which is how that was checked.

    Two consequences for the caller, both factors of two:
      - the physical antenna sees HALF the full structure's input impedance;
      - it radiates into the upper half-space only, so its directivity is
        TWICE the full structure's free-space directivity.

    The ground plane is infinite; a real helix sits on a finite disc or cup,
    which mostly changes the back lobe rather than the forward gain.
    """
    R = circumference / (2.0 * math.pi)
    S = circumference * math.tan(math.radians(pitch_deg))
    t = np.linspace(0.0, 2.0 * math.pi * turns, turns * seg_per_turn + 1)
    up = np.stack([R * np.cos(t), R * np.sin(t),
                   feed_height + S * t / (2.0 * math.pi)], axis=1)
    down = up[::-1].copy()
    down[:, 2] *= -1.0
    zf = np.linspace(-feed_height, feed_height, 5)[1:-1]
    feed = np.stack([np.full_like(zf, R), np.zeros_like(zf), zf], axis=1)
    model = WireModel([Wire(np.vstack([down, feed, up]), radius)])
    nodes = np.array([model.node_of(n) for n in range(model.n_basis)])
    idx = int(np.argmin(np.abs(nodes[:, 2]) + 10.0 * np.abs(nodes[:, 0] - R)
                        + 10.0 * np.abs(nodes[:, 1])))
    return model, idx


# ------------------------------------------------------- feed networks
#
# A non-radiating network - a feeder line, a phasing section, a matching
# stub - coupled to the wires in admittance form, which is how NEC's TL cards
# work. With delta-gap voltages V at the port basis functions the antenna's
# port currents are Y_ant V, Y_ant = (Z^-1)[ports, ports]; the network adds its
# own admittance over the same nodes plus any internal ones (a generator node
# that is not on a wire, say); and Kirchhoff's current law at every node gives
# (Y_net + Y_ant) V = J for injected source currents J.
#
# Checked three ways before use: a dipole behind a length of line reproduces
# the textbook impedance transformation to machine precision; two dipoles tied
# in parallel through a stiff network match the plain solver driving both gaps
# at once; and a lossless feeder delivers exactly the power the pattern
# integral says is radiated.


def tl_admittance(z0: float, electrical_length: float,
                  transposed: bool = False) -> np.ndarray:
    """2x2 admittance of a lossless line, electrical length beta*l in radians.

    `transposed` reverses the second port - the crisscross an LPDA's feeder
    makes between neighbouring elements, without which the array fires
    backwards.
    """
    s = 1.0 / (1j * z0 * math.sin(electrical_length))
    y12 = -s if not transposed else s
    return np.array([[math.cos(electrical_length) * s, y12],
                     [y12, math.cos(electrical_length) * s]], dtype=complex)


@dataclass(frozen=True)
class NetworkSolution:
    solution: MoMSolution
    node_voltages: np.ndarray
    source: int

    @property
    def input_impedance(self) -> complex:
        """Voltage at the source node for the unit injected current."""
        return complex(self.node_voltages[self.source])

    @property
    def generator_power(self) -> float:
        return 0.5 * self.input_impedance.real


def solve_network(model: WireModel, ports: list[int], y_net: np.ndarray,
                  source: int) -> NetworkSolution:
    """Drive the network with 1 A injected at node `source`.

    Network nodes 0 .. len(ports)-1 are the antenna ports, in the order given;
    any further rows of `y_net` are internal nodes. The returned solution's
    wire currents give patterns and radiated power exactly as a plain solve's
    do.
    """
    P = len(ports)
    n = y_net.shape[0]
    if n < P or y_net.shape != (n, n):
        raise ValueError("y_net must be square and cover every port")
    Z = impedance_matrix(model)
    E = np.zeros((model.n_basis, P), dtype=complex)
    for k, m in enumerate(ports):
        E[m, k] = 1.0
    X = np.linalg.solve(Z, E)
    Y = np.array(y_net, dtype=complex)
    Y[:P, :P] += X[ports, :]
    J = np.zeros(n, dtype=complex)
    J[source] = 1.0
    V = np.linalg.solve(Y, J)
    return NetworkSolution(MoMSolution(model, X @ V[:P], ports[0]), V, source)


def lpda_model(tau: float, sigma: float, longest: float, elements: int,
               length_over_diameter: float = 125.0, feeder_z0: float = 100.0,
               transposed: bool = True, seg_per_lambda: int = 24):
    """Log-periodic dipole array on Carrel's apex geometry, in wavelengths.

    Returns (model, ports, y_net, source). Elements are z-directed, element 0 is
    the LONGEST, the apex is at the origin and the array fires towards it
    (-x). Each element's radius scales with its length, so the structure is
    truly log-periodic; the feeder is a line of the physical spacing between
    neighbours, transposed; the generator is at the shortest element.
    """
    tan_a = (1.0 - tau) / (4.0 * sigma)
    lengths = [longest * tau ** k for k in range(elements)]
    wires, ports, base = [], [], 0
    for L in lengths:
        nseg = max(8, int(math.ceil(seg_per_lambda * L)))
        nseg += nseg % 2
        z = np.linspace(-0.5 * L, 0.5 * L, nseg + 1)
        nodes = np.stack([np.full_like(z, L / (2.0 * tan_a)), np.zeros_like(z), z],
                         axis=1)
        wires.append(Wire(nodes, L / (2.0 * length_over_diameter)))
        ports.append(base + nseg // 2 - 1)
        base += nseg - 1
    y = np.zeros((elements, elements), dtype=complex)
    for k in range(elements - 1):
        d = lengths[k] * (1.0 - tau) / (2.0 * tan_a)
        idx = [k, k + 1]
        y[np.ix_(idx, idx)] += tl_admittance(feeder_z0, 2.0 * math.pi * d, transposed)
    return WireModel(wires), ports, y, elements - 1


def rhombic_model(leg: float, half_angle_deg: float | None = None,
                  radius: float = 1e-4, seg_per_lambda: int = 28):
    """Free-space planar rhombic, in wavelengths. Returns (model, feed, term).

    One closed wire of four legs in the x-y plane: the feed at the acute vertex
    on the origin, the terminating resistor's basis function at the far acute
    vertex on +x, towards which the antenna fires. Pass the resistor to
    `solve` as `loads={term: R}`. The default half-angle is the alignment
    angle, acos(1 - 0.371/leg), where each leg's cone lines up with the axis.

    Mesh matters here more than on a dipole: at 12 segments per wavelength the
    axial directivity reads 0.35 dB low on a 4-wavelength leg, and 28 is
    within about 0.05 dB of converged.
    """
    if half_angle_deg is None:
        half_angle_deg = math.degrees(math.acos(max(0.0, 1.0 - 0.371 / leg)))
    th = math.radians(half_angle_deg)
    p0 = np.zeros(3)
    p1 = np.array([leg * math.cos(th), leg * math.sin(th), 0.0])
    p2 = np.array([2.0 * leg * math.cos(th), 0.0, 0.0])
    p3 = np.array([leg * math.cos(th), -leg * math.sin(th), 0.0])
    n = max(8, int(round(seg_per_lambda * leg)))
    pts = [a + (j / n) * (b - a)
           for a, b in ((p0, p1), (p1, p2), (p2, p3), (p3, p0)) for j in range(n)]
    model = WireModel([Wire(np.array(pts), radius, closed=True)])
    nodes = np.array([model.node_of(k) for k in range(model.n_basis)])
    return (model, int(np.argmin(np.linalg.norm(nodes - p0, axis=1))),
            int(np.argmin(np.linalg.norm(nodes - p2, axis=1))))


def top_hat_monopole(height: float, hat_radius: float, radials: int,
                     radius: float = 1e-4, segments: int = 40):
    """Monopole with a hat of radial wires, on an infinite PEC plane by images.

    Returns (model, feed, top_junction). The full structure is one vertical
    wire from -height to +height fed at its centre, `radials` wires at the top
    and their images at the bottom, joined by junctions; the physical antenna
    sees half its impedance and radiates into the upper half-space only.
    `radials` = 0 gives a bare whip. Lengths in wavelengths.

    Read the base current from the vertical AWAY from the feed: on a wire this
    short the delta gap puts a local excess on the feed node itself (see the
    note at the top of this module).
    """
    z = np.linspace(-height, height, segments + 1)
    wires = [Wire(np.stack([np.zeros_like(z), np.zeros_like(z), z], axis=1), radius)]
    if radials > 0 and hat_radius > 0:
        step = 2.0 * height / segments
        # cap the radial mesh: tying it to the vertical's step gives a wide hat
        # on a short whip thousands of segments; 12 keeps the length ratio
        # across the junction at 5 or less, which thin-wire codes tolerate
        nr = max(3, min(int(math.ceil(hat_radius / step)), 12))
        for zz in (height, -height):
            for k in range(radials):
                phi = 2.0 * math.pi * k / radials
                r = np.linspace(0.0, hat_radius, nr + 1)
                wires.append(Wire(np.stack([r * math.cos(phi), r * math.sin(phi),
                                            np.full_like(r, zz)], axis=1), radius))
    model = WireModel(wires)
    top = model.junction_at([0.0, 0.0, height]) if model.junctions else ()
    return model, segments // 2 - 1, top


def bicone_cage(slant: float, half_angle_deg: float, wires: int = 16,
                radius: float = 0.002, gap: float = 0.01, seg: float = 0.03):
    """Biconical antenna as a wire cage, in wavelengths. Returns (model, feed).

    Two cones of `wires` wires each, half angle measured from the axis, joined
    at their apexes to a short feed wire of length 2*gap that carries the delta
    gap. A conical monopole on an infinite PEC plane is half of it: half the
    impedance, twice the directivity.

    Directivity converges quickly in the wire count - 0.2% between 16 and 24 -
    but the IMPEDANCE of a cage does not: it was still moving 3-5% between 24
    and 32 wires. Use a cage for patterns, not for the input impedance of a
    solid cone.
    """
    th = math.radians(half_angle_deg)
    parts = [Wire(np.array([[0.0, 0.0, -gap], [0.0, 0.0, 0.0], [0.0, 0.0, gap]]), radius)]
    ns = max(4, int(math.ceil(slant / seg)))
    t = np.linspace(0.0, slant, ns + 1)
    for sgn in (1.0, -1.0):
        for k in range(wires):
            ph = 2.0 * math.pi * k / wires
            parts.append(Wire(np.stack([t * math.sin(th) * math.cos(ph),
                                        t * math.sin(th) * math.sin(ph),
                                        sgn * (gap + t * math.cos(th))], axis=1), radius))
    return WireModel(parts), 0
