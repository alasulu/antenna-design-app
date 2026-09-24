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

What it does not do. No ground plane, no dielectric, no loss, no junctions of
more than two wires, and the reduced kernel rather than the exact one. Each of
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
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from ..core.constants import ETA0

__all__ = [
    "Wire", "WireModel", "MoMSolution", "dipole", "loop", "arc", "halo",
    "folded_dipole_wire", "helix_over_ground", "lpda_model",
    "tl_admittance", "NetworkSolution", "solve_network",
    "solve", "input_impedance", "far_field", "directivity", "radiated_power",
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

    Basis function n is centred on an interior node and spans the two segments
    meeting there, rising to 1 at the node. An open wire of N segments carries
    N-1 of them, which is what forces the current to zero at the free ends; a
    closed wire carries N, the node indices wrapping.
    """

    wires: list[Wire]
    seg_a: np.ndarray = field(init=False)      # segment start points
    seg_b: np.ndarray = field(init=False)      # segment end points
    seg_t: np.ndarray = field(init=False)      # unit tangents
    seg_len: np.ndarray = field(init=False)
    seg_rad: np.ndarray = field(init=False)
    basis: np.ndarray = field(init=False)      # (Nb, 2) segment indices (minus, plus)

    def __post_init__(self) -> None:
        a, b, rad, basis = [], [], [], []
        for w in self.wires:
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
            else:
                for i in range(n - 1):
                    basis.append((base + i, base + i + 1))
        self.seg_a = np.array(a, dtype=float)
        self.seg_b = np.array(b, dtype=float)
        d = self.seg_b - self.seg_a
        self.seg_len = np.linalg.norm(d, axis=1)
        if np.any(self.seg_len <= 0):
            raise ValueError("a segment has zero length")
        self.seg_t = d / self.seg_len[:, None]
        self.seg_rad = np.array(rad, dtype=float)
        self.basis = np.array(basis, dtype=int).reshape(-1, 2)
        if len(self.basis) == 0:
            raise ValueError("model has no basis functions; wires are too coarse")

    @property
    def n_basis(self) -> int:
        return len(self.basis)

    def node_of(self, n: int) -> np.ndarray:
        """Position of the node basis function `n` is centred on."""
        return self.seg_b[self.basis[n, 0]]


# ---------------------------------------------------------------- kernels

def _segment_moments(model: WireModel):
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
    for q in range(ns):
        A, u, L = model.seg_a[q], model.seg_t[q], model.seg_len[q]
        a2 = model.seg_rad[q] ** 2
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
        for b, inner in ((0, inner0), (1, inner1)):
            g = inner.reshape(ns, -1)                                    # (ns, no)
            M[0, b, :, q] = (0.5 * model.seg_len) * (g * wo[None, :]).sum(axis=1)
            M[1, b, :, q] = (0.5 * model.seg_len) * (g * so * wo[None, :]).sum(axis=1)
    return M


def _half_weights(model: WireModel, n: int):
    """For basis n: (segment, (alpha, beta) of f, sigma).

    On the minus segment f rises 0 -> 1 towards the node and the divergence is
    +1/L; on the plus segment f falls 1 -> 0 away from it and the divergence is
    -1/L. Both segments already point the way the current flows, which is the
    reason for orienting them so - there is no tangent sign to track.
    """
    m, pl = model.basis[n]
    Lm, Lp = model.seg_len[m], model.seg_len[pl]
    return [(m, (0.0, 1.0 / Lm), 1.0 / Lm),
            (pl, (1.0, -1.0 / Lp), -1.0 / Lp)]


def impedance_matrix(model: WireModel) -> np.ndarray:
    """The Galerkin EFIE matrix [ohm]."""
    M = _segment_moments(model)
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
    feed: int

    @property
    def feed_current(self) -> complex:
        return complex(self.currents[self.feed])

    @property
    def input_impedance(self) -> complex:
        return 1.0 / self.feed_current

    @property
    def circuit_power(self) -> float:
        """P = 0.5 Re(V I*) for a 1 V delta gap [W]."""
        return 0.5 * float(np.real(np.conj(self.feed_current)))


def solve(model: WireModel, feed: int | None = None,
          loads: dict[int, complex] | None = None) -> MoMSolution:
    """Delta-gap excitation of one basis function; 1 V across the gap.

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
    if not 0 <= feed < nb:
        raise ValueError(f"feed index {feed} outside 0..{nb - 1}")
    Z = impedance_matrix(model)
    for m, zl in (loads or {}).items():
        if not 0 <= m < nb:
            raise ValueError(f"load index {m} outside 0..{nb - 1}")
        Z[m, m] += zl
    V = np.zeros(nb, dtype=complex)
    V[feed] = 1.0
    return MoMSolution(model, np.linalg.solve(Z, V), feed)


def input_impedance(model: WireModel, feed: int | None = None,
                    loads: dict[int, complex] | None = None) -> complex:
    return solve(model, feed, loads).input_impedance


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
    small = np.abs(g) < 1e-9
    gs = np.where(small, 1.0, g)
    e = np.exp(1j * gs * L)
    i0 = np.where(small, L, (e - 1.0) / (1j * gs))
    i1 = np.where(small, 0.5 * L * L,
                  (L * e) / (1j * gs) - (e - 1.0) / (1j * gs) ** 2)
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
    if z_of_scale(lo).imag * z_of_scale(hi).imag > 0:
        return None
    for _ in range(steps):
        mid = 0.5 * (lo + hi)
        if z_of_scale(mid).imag < 0:
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
