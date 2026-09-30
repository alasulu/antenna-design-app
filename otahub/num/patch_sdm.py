"""A rectangular patch on a grounded slab by the spectral-domain method of
moments - a full-wave arbiter for the patch family.

`patch_q` radiates the cavity model's current through the slab; this module
does not assume the current. The patch (L along x, W along y, centred on the
origin) carries x- and y-directed currents expanded in entire-domain basis
functions that meet the edge conditions exactly,

    J_x = U_n(2x/L) sqrt(1 - (2x/L)^2) * T_j(2y/W) / sqrt(1 - (2y/W)^2)
    J_y = U_m(2y/W) sqrt(1 - (2y/W)^2) * T_j(2x/L) / sqrt(1 - (2x/L)^2)

- a current normal to an edge vanishes as the square root of the distance, one
parallel to it is singular as its inverse - with closed-form spectra (Bessel
functions). Galerkin's method is applied with the exact spectral Green's
function of the grounded slab from its transverse equivalent network. The
spectral integrals run over a contour lifted above the real axis past the
surface-wave poles, then along it; the slowly converging reactance is
extrapolated from two truncations (Richardson, error ~ 1/k_max), which share
their common range. For the dominant (TM10) symmetry class the integrand is
even in both spectral axes and one quadrant is integrated.

A patch's natural mode is the characteristic mode (X v = lambda R v) that
radiates most, and its resonance is where lambda = 0. From its current come the
broadside directivity, the radiation efficiency (space wave against space plus
surface wave) and the Q, from d lambda / d omega.

Checked against independent answers: a vanishing current element against
Jackson and Alexopoulos's thin-slab space and surface waves; an infinite
microstrip line's effective permittivity against Kirschning and Jansen; the
radiated power computed from the spectral reaction and from the far field; and
the resonance and Q of a thick eps_r 10.2 patch against an FDTD ringdown
(`patch_fdtd`) extrapolated in cell size. A basis of sines along the current
(linear, not square-root, at the open ends) was tried first and wandered by
+-1.5% in resonance; this one moves by 1e-4 from its default set to the next.
SI units, exp(+j omega t).
"""
from __future__ import annotations

import math

import numpy as np
from scipy.linalg import eigh
from scipy.special import jv, yv

from ..core.constants import C0, EPS0, ETA0, MU0

__all__ = ["RectPatch", "BASIS_SMALL", "BASIS_FULL", "BASIS_EVEN", "resonance", "resonant_length", "mode_metrics", "slab_gf",
           "layered_z", "StackedPatch", "tracked_mode", "stacked_resonance", "probe_vector", "probe_reactance",
           "input_impedance"]

# (n, j) for J_x, (m, j) for J_y; the TM10 class: J_x even-even, J_y odd-odd
BASIS_SMALL = dict(bx=((0, 0), (0, 2), (0, 4)), by=())
BASIS_FULL = dict(bx=((0, 0), (0, 2), (0, 4), (2, 0), (2, 2), (4, 0)), by=((1, 1), (1, 3), (3, 1)))
# the partner class, E_z even in x (J_x odd in x, J_y even in x): no resonance near TM10's
BASIS_EVEN = dict(bx=((1, 0), (1, 2), (3, 0)), by=((1, 0), (1, 2), (3, 0)))


def _usqrt_ft(k, a_len: float, n: int):
    """FT of U_n(2x/A) sqrt(1 - (2x/A)^2) on [-A/2, A/2]: (A/2) pi (n+1) j^n J_{n+1}(z)/z, z = kA/2."""
    z = np.asarray(k, dtype=complex) * a_len / 2
    small = np.abs(z) < 1e-8
    zz = np.where(small, 1.0, z)
    lim = 0.5 if n == 0 else 0.0
    return 0.5 * a_len * math.pi * (n + 1) * (1j ** n) * np.where(small, lim, jv(n + 1, zz) / zz)


def _cheb_ft(k, a_len: float, j: int):
    """FT of T_j(2x/A) / sqrt(1 - (2x/A)^2) on [-A/2, A/2]: (pi A/2) j^j J_j(k A/2)."""
    return 0.5 * math.pi * a_len * (1j ** j) * jv(j, np.asarray(k, dtype=complex) * a_len / 2)


def _sine_ft(k, a_len: float, n: int):
    """FT of sin(n pi (x + A/2)/A) on [-A/2, A/2] (kept for line and dipole checks)."""
    a = n * math.pi / a_len
    k = np.asarray(k, dtype=complex)
    k = np.where(np.abs(a * a - k * k) < 1e-9 * a * a, k + 1e-6 * a, k)
    return np.exp(-1j * k * a_len / 2) * a * (1 - (-1) ** n * np.exp(1j * k * a_len)) / (a * a - k * k)


def slab_gf(theta: float, eps_r: float, k0h: float):
    """Far-field factors of a horizontal current on the slab: TM (G) and TE (F)."""
    st, ct = math.sin(theta), math.cos(theta)
    N = np.sqrt(eps_r - st ** 2 + 0j)
    t = np.tan(k0h * N)
    return 2j * N * t * ct / (1j * N * t + eps_r * ct), 2j * t * ct / (1j * t * ct + N)


class RectPatch:
    def __init__(self, eps_r: float, h: float, L: float, W: float, bx=((0, 0),), by=()):
        self.er, self.h, self.L, self.W = eps_r, h, L, W
        self.basis = [("x", n, j) for n, j in bx] + [("y", m, j) for m, j in by]
        # one quadrant suffices when every function is in the TM10 class, or every one in
        # its partner with E_z even in x (J_x odd in x) - the class a probe off the centre
        # line also excites
        self.quadrant = (all((n % 2 == 0 and j % 2 == 0) if d == "x" else (n % 2 == 1 and j % 2 == 1)
                             for d, n, j in self.basis)
                         or all((n % 2 == 1 and j % 2 == 0) if d == "x" else (n % 2 == 1 and j % 2 == 0)
                                for d, n, j in self.basis))
        self.x_even_ez = any(d == "x" and n % 2 == 1 for d, n, j in self.basis)

    def ft(self, kx, ky):
        """(n_basis, 2, ...) spectra of the basis functions."""
        out = []
        for d, n, j in self.basis:
            if d == "x":
                f = _usqrt_ft(kx, self.L, n) * _cheb_ft(ky, self.W, j)
                out.append(np.stack([f, np.zeros_like(f)]))
            else:
                f = _usqrt_ft(ky, self.W, n) * _cheb_ft(kx, self.L, j)
                out.append(np.stack([np.zeros_like(f), f]))
        return np.array(out)

    def znode(self, kr, k0: float):
        """TM and TE impedances at the slab surface, air above, shorted slab below."""
        w = k0 * C0
        kz0 = -1j * np.sqrt(kr * kr - k0 * k0 + 0j)
        kz1 = -1j * np.sqrt(kr * kr - self.er * k0 * k0 + 0j)
        x = kz1 * self.h
        # cos/sin overflows for |Im x| in the hundreds (a thick slab far out in the
        # evanescent spectrum); the stable form takes over there, leaving every
        # ordinary value exactly as it was
        with np.errstate(over="ignore", invalid="ignore"):
            cot = np.where(np.abs(np.imag(x)) > 20, _cot(np.asarray(x, dtype=complex)), np.cos(x) / np.sin(x))
        ytm = w * EPS0 / kz0 - 1j * (w * EPS0 * self.er / kz1) * cot
        yte = kz0 / (w * MU0) - 1j * (kz1 / (w * MU0)) * cot
        return 1 / ytm, 1 / yte

    def _alphas(self, kr):
        # a product of two basis spectra has azimuthal harmonics up to ~ |kr| max(L, W)
        n = int(4 * math.ceil((1.25 * abs(kr) * max(self.L, self.W) + 32) / 4))
        if self.quadrant:
            m = n // 4
            al = 0.5 * math.pi * np.arange(m + 1) / m
            w = np.full(m + 1, 0.5 * math.pi / m)
            w[0] = w[-1] = 0.25 * math.pi / m
            return al, 4 * w
        return 2 * math.pi * np.arange(n) / n, np.full(n, 2 * math.pi / n)

    def _accumulate(self, nodes, weights, k0: float):
        nb = len(self.basis)
        Z = np.zeros((nb, nb), complex)
        for kr, wk in zip(nodes, weights):
            al, wa = self._alphas(kr)
            ca, sa = np.cos(al), np.sin(al)
            Fp = self.ft(kr * ca, kr * sa)
            Fm = self.ft(-kr * ca, -kr * sa)
            ztm, zte = self.znode(kr, k0)
            gxx = ztm * ca ** 2 + zte * sa ** 2
            gyy = ztm * sa ** 2 + zte * ca ** 2
            gxy = (ztm - zte) * ca * sa
            Ex = (gxx * Fp[:, 0] + gxy * Fp[:, 1]) * wa
            Ey = (gxy * Fp[:, 0] + gyy * Fp[:, 1]) * wa
            Z += wk * kr * (Fm[:, 0] @ Ex.T + Fm[:, 1] @ Ey.T)
        return Z / (4 * math.pi ** 2)

    def _nodes(self, k0: float, a: float, b: float, panels: int):
        x, w = np.polynomial.legendre.leggauss(32)
        nodes, weights = [], []
        edges = np.geomspace(a, b, panels + 1)
        for lo, hi in zip(edges[:-1], edges[1:]):
            nodes += list(0.5 * (hi - lo) * (x + 1) + lo)
            weights += list(0.5 * (hi - lo) * w)
        return np.array(nodes), np.array(weights)

    def Z(self, f: float, kmax: float = 400.0, n_ell: int = 400, panels: int = 48):
        """Galerkin matrix: the spectral integral truncated at kmax/2 and kmax * k0,
        extrapolated in 1/k_max (the two truncations share their common range)."""
        k0 = 2 * math.pi * f / C0
        K1 = k0 * (math.sqrt(self.er) + 1.0)          # past every surface-wave pole
        t, wt = np.polynomial.legendre.leggauss(n_ell)
        t = 0.5 * math.pi * (t + 1)
        wt = 0.5 * math.pi * wt
        b = 0.2 * k0
        ell = (0.5 * K1 * (1 - np.cos(t)) + 1j * b * np.sin(t), wt * (0.5 * K1 * np.sin(t) + 1j * b * np.cos(t)))
        near = self._nodes(k0, K1, 0.5 * kmax * k0, panels // 2)
        far = self._nodes(k0, 0.5 * kmax * k0, kmax * k0, panels // 2)
        A = self._accumulate(ell[0], ell[1], k0) + self._accumulate(*near, k0)
        return A + 2 * self._accumulate(*far, k0)

    def Z_trunc(self, f: float, kmax: float, n_ell: int = 400, panels: int = 48):
        """The spectral integral truncated at kmax * k0, not extrapolated."""
        k0 = 2 * math.pi * f / C0
        K1 = k0 * (math.sqrt(self.er) + 1.0)
        t, wt = np.polynomial.legendre.leggauss(n_ell)
        t = 0.5 * math.pi * (t + 1)
        wt = 0.5 * math.pi * wt
        b = 0.2 * k0
        Z = self._accumulate(0.5 * K1 * (1 - np.cos(t)) + 1j * b * np.sin(t),
                             wt * (0.5 * K1 * np.sin(t) + 1j * b * np.cos(t)), k0)
        return Z + self._accumulate(*self._nodes(k0, K1, kmax * k0, panels), k0)

    def R_space(self, f: float, n: int = 200):
        """The space-wave (radiated) part of Re Z: 0 < kr < k0, as kr = k0 sin(tau)."""
        k0 = 2 * math.pi * f / C0
        t, wt = np.polynomial.legendre.leggauss(n)
        t = 0.25 * math.pi * (t + 1)
        wt = 0.25 * math.pi * wt
        return self._accumulate(k0 * np.sin(t) + 0j, k0 * np.cos(t) * wt, k0).real

    def far_field(self, v, f: float, nt: int = 96, nph: int = 192):
        """(space-wave power by far-field integration, broadside intensity)."""
        k0 = 2 * math.pi * f / C0
        xt, wt = np.polynomial.legendre.leggauss(nt)
        th = 0.25 * math.pi * (xt + 1)
        wth = 0.25 * math.pi * wt
        ph = 2 * math.pi * np.arange(nph) / nph
        total = 0.0
        for t_, w_ in zip(th, wth):
            G, F = slab_gf(t_, self.er, k0 * self.h)
            kx, ky = k0 * math.sin(t_) * np.cos(ph), k0 * math.sin(t_) * np.sin(ph)
            J = np.tensordot(v, self.ft(kx, ky), axes=(0, 0))
            Jr = J[0] * np.cos(ph) + J[1] * np.sin(ph)
            Jp = -J[0] * np.sin(ph) + J[1] * np.cos(ph)
            total += w_ * math.sin(t_) * float(np.sum(abs(G) ** 2 * abs(Jr) ** 2 + abs(F) ** 2 * abs(Jp) ** 2)) * 2 * math.pi / nph
        c = ETA0 * k0 ** 2 / (32 * math.pi ** 2)
        G0, _ = slab_gf(0.0, self.er, k0 * self.h)
        J0 = np.tensordot(v, self.ft(np.array([0.0]), np.array([0.0])), axes=(0, 0))[:, 0]
        return c * total, c * abs(G0) ** 2 * (abs(J0[0]) ** 2 + abs(J0[1]) ** 2)


def _dominant(Z):
    R = Z.real + 1e-9 * np.max(np.diag(Z.real)) * np.eye(len(Z))
    lam, V = eigh(Z.imag, R)
    j = int(np.argmin(np.linalg.norm(V, axis=0)))     # R-normalised: the strongest radiator
    return lam[j], V[:, j]


def _secant(g, x0: float, x1: float, tol: float = 2e-6, it: int = 20) -> float:
    g0, g1 = g(x0), g(x1)
    for _ in range(it):
        x2 = x1 - g1 * (x1 - x0) / (g1 - g0)
        if abs(x2 / x1 - 1) < tol:
            return x2
        x0, g0, x1, g1 = x1, g1, x2, g(x2)
    raise ValueError("resonance search did not converge")


def resonance(p: RectPatch, f_guess: float) -> float:
    """Frequency at which the dominant characteristic mode's eigenvalue is zero."""
    return _secant(lambda f: _dominant(p.Z(f))[0], f_guess, 1.01 * f_guess)


def resonant_length(eps_r: float, h: float, W: float, f0: float, L_guess: float, basis=None) -> float:
    """The patch length whose dominant mode resonates at f0, width held."""
    basis = basis or BASIS_FULL
    return _secant(lambda L: _dominant(RectPatch(eps_r, h, L, W, **basis).Z(f0))[0], L_guess, 0.99 * L_guess)


def mode_metrics(p: RectPatch, f: float) -> dict:
    """Directivity, radiation efficiency and Q of the dominant mode at f."""
    Z = p.Z(f)
    _, v = _dominant(Z)
    lp = v @ p.Z(f * 1.002).imag @ v
    lm = v @ p.Z(f * 0.998).imag @ v
    p_total = 0.5 * v @ Z.real @ v
    p_space = 0.5 * v @ p.R_space(f) @ v
    q_total = (lp - lm) / (0.008 * 2 * p_total)        # (omega/2) d lambda / d omega
    p_far, u0 = p.far_field(v, f)
    return dict(directivity=4 * math.pi * u0 / p_far, efficiency=p_space / p_total, q_total=q_total,
                q_radiation=q_total * p_total / p_space, far_over_spectral=p_far / p_space, current=v / v[0])


# ------------------------------------------------------------ two layers and a probe

def _cot(x):
    """cot, stable for any imaginary part (the evanescent spectrum's large |kz h|)."""
    e = np.exp(np.where(x.imag < 0, -2j * x, 2j * x))
    return np.where(x.imag < 0, 1j * (1 + e) / (1 - e), -1j * (1 + e) / (1 - e))


def _tan(x):
    e = np.exp(np.where(x.imag < 0, -2j * x, 2j * x))
    return np.where(x.imag < 0, -1j * (1 - e) / (1 + e), 1j * (1 - e) / (1 + e))


def _sec(x):
    e = np.exp(np.where(x.imag < 0, -1j * x, 1j * x))
    return 2 * e / (1 + e * e)


def layered_z(kr, k0: float, er1: float, h1: float, er2: float, h2: float):
    """TM and TE transfer impedances [[Z11, Z12], [Z21, Z22]] at the two interfaces of a
    grounded two-layer slab: short at z = 0, eps_r1 for h1 (node 1), eps_r2 for h2 (node 2),
    air above. Z_ij is the voltage at node i per unit current injected at node j. Written
    in tan, cot and sec so it holds however evanescent the spectrum."""
    w = k0 * C0
    kr = np.asarray(kr, complex)
    kz0 = -1j * np.sqrt(kr * kr - k0 * k0 + 0j)
    kz1 = -1j * np.sqrt(kr * kr - er1 * k0 * k0 + 0j)
    kz2 = -1j * np.sqrt(kr * kr - er2 * k0 * k0 + 0j)
    out = []
    for Y0, Y1, Y2 in ((w * EPS0 / kz0, w * EPS0 * er1 / kz1, w * EPS0 * er2 / kz2),
                       (kz0 / (w * MU0), kz1 / (w * MU0), kz2 / (w * MU0))):
        yd1 = -1j * Y1 * _cot(kz1 * h1)                                  # the shorted layer below node 1
        t2 = _tan(kz2 * h2)
        yu1 = Y2 * (Y0 + 1j * Y2 * t2) / (Y2 + 1j * Y0 * t2)              # layer 2 loaded by air, from node 1
        yd2 = Y2 * (yd1 + 1j * Y2 * t2) / (Y2 + 1j * yd1 * t2)            # layer 2 loaded by layer 1, from node 2
        z11 = 1 / (yd1 + yu1)
        z21 = z11 * _sec(kz2 * h2) / (1 + 1j * (Y0 / Y2) * t2)
        out.append(np.array([[z11, z21], [z21, 1 / (yd2 + Y0)]]))
    return out


class StackedPatch(RectPatch):
    """A driven patch L x W on (eps_r, h) and a parasitic L2 x W2 on a second layer
    (eps_r2, h2) above it, both centred on the origin and resonant along x. The basis is
    each patch's own, the lower's first; the Green's function couples them through
    `layered_z`."""

    def __init__(self, eps_r: float, h: float, L: float, W: float, eps_r2: float, h2: float, L2: float,
                 W2: float, bx=((0, 0),), by=()):
        self.lower = RectPatch(eps_r, h, L, W, bx, by)
        self.upper = RectPatch(eps_r2, h2, L2, W2, bx, by)
        self.er1, self.h1, self.er2, self.h2 = eps_r, h, eps_r2, h2
        self.er, self.h = max(eps_r, eps_r2), h           # the contour clears the densest layer's poles
        self.L, self.W = max(L, L2), max(W, W2)           # the angular sampling follows the larger patch
        self.n1 = len(self.lower.basis)
        self.basis = self.lower.basis + self.upper.basis
        self.quadrant = self.lower.quadrant and self.upper.quadrant
        self.x_even_ez = self.lower.x_even_ez

    def ft(self, kx, ky):
        return np.concatenate([self.lower.ft(kx, ky), self.upper.ft(kx, ky)])

    def _accumulate(self, nodes, weights, k0: float):
        nb, n1 = len(self.basis), self.n1
        sl = (slice(0, n1), slice(n1, nb))
        Z = np.zeros((nb, nb), complex)
        for kr, wk in zip(nodes, weights):
            al, wa = self._alphas(kr)
            ca, sa = np.cos(al), np.sin(al)
            Fp = (self.lower.ft(kr * ca, kr * sa), self.upper.ft(kr * ca, kr * sa))
            Fm = (self.lower.ft(-kr * ca, -kr * sa), self.upper.ft(-kr * ca, -kr * sa))
            ztm, zte = layered_z(kr, k0, self.er1, self.h1, self.er2, self.h2)
            for i in range(2):
                for j in range(2):
                    gxx = ztm[i, j] * ca ** 2 + zte[i, j] * sa ** 2
                    gyy = ztm[i, j] * sa ** 2 + zte[i, j] * ca ** 2
                    gxy = (ztm[i, j] - zte[i, j]) * ca * sa
                    Ex = (gxx * Fp[j][:, 0] + gxy * Fp[j][:, 1]) * wa
                    Ey = (gxy * Fp[j][:, 0] + gyy * Fp[j][:, 1]) * wa
                    Z[sl[i], sl[j]] += wk * kr * (Fm[i][:, 0] @ Ex.T + Fm[i][:, 1] @ Ey.T)
        return Z / (4 * math.pi ** 2)

    def _factors(self, theta: float, k0: float):
        """Far-field factors (TM, TE) of a current on each layer: the top node's voltage
        radiates, so they are 2 cos(theta) Z_2i Y0 (TM) and 2 Z_2i Y0 (TE)."""
        kr = k0 * math.sin(theta)
        ztm, zte = layered_z(kr, k0, self.er1, self.h1, self.er2, self.h2)
        w = k0 * C0
        y0tm, y0te = w * EPS0 / (k0 * math.cos(theta)), k0 * math.cos(theta) / (w * MU0)
        return ([2 * math.cos(theta) * ztm[1, i] * y0tm for i in range(2)], [2 * zte[1, i] * y0te for i in range(2)])

    def far_field(self, v, f: float, nt: int = 96, nph: int = 192):
        """(space-wave power by far-field integration, broadside intensity)."""
        k0 = 2 * math.pi * f / C0
        n1 = self.n1
        xt, wt = np.polynomial.legendre.leggauss(nt)
        th = 0.25 * math.pi * (xt + 1)
        wth = 0.25 * math.pi * wt
        ph = 2 * math.pi * np.arange(nph) / nph
        total = 0.0
        for t_, w_ in zip(th, wth):
            G, F = self._factors(t_, k0)
            kx, ky = k0 * math.sin(t_) * np.cos(ph), k0 * math.sin(t_) * np.sin(ph)
            J = (np.tensordot(v[:n1], self.lower.ft(kx, ky), axes=(0, 0)),
                 np.tensordot(v[n1:], self.upper.ft(kx, ky), axes=(0, 0)))
            Er = sum(G[i] * (J[i][0] * np.cos(ph) + J[i][1] * np.sin(ph)) for i in range(2))
            Ep = sum(F[i] * (-J[i][0] * np.sin(ph) + J[i][1] * np.cos(ph)) for i in range(2))
            total += w_ * math.sin(t_) * float(np.sum(abs(Er) ** 2 + abs(Ep) ** 2)) * 2 * math.pi / nph
        c = ETA0 * k0 ** 2 / (32 * math.pi ** 2)
        G0, _ = self._factors(0.0, k0)
        z = np.array([0.0])
        J0 = (np.tensordot(v[:n1], self.lower.ft(z, z), axes=(0, 0))[:, 0],
              np.tensordot(v[n1:], self.upper.ft(z, z), axes=(0, 0))[:, 0])
        E0 = G0[0] * J0[0] + G0[1] * J0[1]
        return c * total, c * (abs(E0[0]) ** 2 + abs(E0[1]) ** 2)


def _modes(Z):
    R = Z.real + 1e-9 * np.max(np.diag(Z.real)) * np.eye(len(Z))
    lam, V = eigh(Z.imag, R)
    return lam, V, R


def tracked_mode(p: RectPatch, f: float, v_ref):
    """(eigenvalue, current) at f of the characteristic mode nearest v_ref in the R inner product."""
    lam, V, R = _modes(p.Z(f))
    j = int(np.argmax(np.abs(v_ref @ R @ V) / np.linalg.norm(V, axis=0)))
    return lam[j], V[:, j]


def stacked_resonance(p: RectPatch, f_guess: float, tol: float = 2e-6, it: int = 25):
    """(frequency, current) of the characteristic mode that, of the three strongest
    radiators at f_guess, has the eigenvalue nearest zero - followed by its current as
    the frequency moves, so a stack's two coupled modes are found one at a time."""
    lam, V, _ = _modes(p.Z(f_guess))
    j = min(np.argsort(np.linalg.norm(V, axis=0))[:3], key=lambda i: abs(lam[i]))
    f1, (g1, v1) = f_guess, (lam[j], V[:, j])
    f2 = 1.01 * f_guess
    g2, v2 = tracked_mode(p, f2, v1)
    for _ in range(it):
        f3 = f2 - g2 * (f2 - f1) / (g2 - g1)
        if abs(f3 / f2 - 1) < tol:
            return f3, v2
        f1, g1, v1, f2 = f2, g2, v2, f3
        g2, v2 = tracked_mode(p, f2, v1)
    raise ValueError("resonance search did not converge")


def _probe_terms(p: RectPatch, kr, k0: float):
    """[(slice of basis, TM transfer impedance from its layer to the probe's)], eps_r, h of the probe's layer."""
    if isinstance(p, StackedPatch):
        ztm, _ = layered_z(kr, k0, p.er1, p.h1, p.er2, p.h2)
        n1, nb = p.n1, len(p.basis)
        return [(slice(0, n1), ztm[0, 0]), (slice(n1, nb), ztm[0, 1])], p.er1
    return [(slice(0, len(p.basis)), p.znode(kr, k0)[0])], p.er


def _odd_in_kx(p: RectPatch) -> np.ndarray:
    """Per basis function: is its radial spectrum odd in kx? J_x even in x (n even) or
    J_y odd in x (j odd) makes it so - the TM10 class, coupled through sin(kx xp); the
    others (E_z even in x) through cos(kx xp)."""
    return np.array([(n % 2 == 0) if d == "x" else (j % 2 == 1) for d, n, j in p.basis])


def _probe_acc(p: RectPatch, nodes, weights, k0: float, xp, a: float):
    """Per-node contributions (len(nodes), len(xp), n_basis) to the probe coupling."""
    out = np.zeros((len(nodes), len(xp), len(p.basis)), complex)
    odd = _odd_in_kx(p)
    for m, (kr, wk) in enumerate(zip(nodes, weights)):
        al, wa = p._alphas(kr)
        ca, sa = np.cos(al), np.sin(al)
        terms, er1 = _probe_terms(p, kr, k0)
        F = p.ft(kr * ca, kr * sa)
        ju = (F[:, 0] * ca + F[:, 1] * sa) * wa                          # (n_basis, n_alpha)
        # e^{-j kx xp}: its odd part for the TM10 class, its even part for the other,
        # chosen per function (one flag for the whole basis zeroed TM10's coupling
        # as soon as a single even function joined it)
        arg = kr * np.multiply.outer(xp, ca)
        proj = np.empty((len(xp), len(p.basis)), complex)
        if odd.any():
            proj[:, odd] = np.sin(arg) @ ju[odd].T
        if (~odd).any():
            proj[:, ~odd] = 1j * np.cos(arg) @ ju[~odd].T
        rad = jv(0, kr * a) if a else 1.0
        for sl, zt in terms:
            out[m, :, sl] = wk * kr * proj[:, sl] * (-kr * zt * rad / (er1 * k0 * k0 - kr * kr))
    return out / (4 * math.pi ** 2)


def probe_vector(p: RectPatch, f: float, xp, a: float = 0.0, kmax: float = 200.0, n_ell: int = 400):
    """Coupling of each basis function to a unit probe current from the ground to the
    driven patch at (xp, 0), radius a: <J_n, -E_probe>, the Galerkin matrix's own sign;
    (len(xp), n_basis) for several positions at once. Per spectral component, a patch
    current's field integrated up the probe is j kr Z_TM / kz^2 times its component along
    kr (the probe layer's transmission-line current). The tail falls as 1/k^2 but
    oscillates with k xp: the partial integrals are averaged over the last period of
    k xp, which cancels the oscillation, at kmax/2 and kmax, and those two extrapolated
    in 1/k^2."""
    xs = np.atleast_1d(np.asarray(xp, float))
    all_odd = bool(_odd_in_kx(p).all())
    if not xs.any() and all_odd:                          # the centre line: odd symmetry, no coupling
        out = np.zeros((len(xs), len(p.basis)), complex)
        return out[0] if np.ndim(xp) == 0 else out
    k0 = 2 * math.pi * f / C0
    K1 = k0 * (math.sqrt(p.er) + 1.0)
    t, wt = np.polynomial.legendre.leggauss(n_ell)
    t = 0.5 * math.pi * (t + 1)
    wt = 0.5 * math.pi * wt
    b = 0.2 * k0
    base = _probe_acc(p, 0.5 * K1 * (1 - np.cos(t)) + 1j * b * np.sin(t),
                      wt * (0.5 * K1 * np.sin(t) + 1j * b * np.cos(t)), k0, xs, a).sum(axis=0)
    width = 2 * math.pi / max(np.abs(xs).max(), 0.05 * C0 / f) / 8   # eight panels to the shortest period of k |xp|
    x, w = np.polynomial.legendre.leggauss(16)
    edges = np.arange(K1, kmax * k0 + 0.5 * width, width)
    nodes = (0.5 * (edges[1:] - edges[:-1])[:, None] * (x + 1) + edges[:-1, None]).ravel()
    weights = (0.5 * (edges[1:] - edges[:-1])[:, None] * w).ravel()
    acc = _probe_acc(p, nodes, weights, k0, xs, a)
    partial = np.cumsum(acc.reshape(len(edges) - 1, 16, len(xs), -1).sum(axis=1), axis=0)
    out = np.empty((len(xs), len(p.basis)), complex)
    for q, x_ in enumerate(xs):
        if x_ == 0.0 and all_odd:                         # the centre line: odd symmetry, no coupling
            out[q] = 0.0
            continue
        # panels in one period of k |x_| (the sign of x_ does not change the period)
        n = max(1, int(round(2 * math.pi / abs(x_) / width))) if x_ else 1
        part, ends = partial[:, q], edges[1:]
        # the half-way window needs a whole period past K1: near the centre line the
        # period outgrows half the range (the window then wrapped round, and returned
        # nan or a wrong mean), so the tail is carried further for that probe alone
        need = 2.0 * (K1 + (n + 2) * width)
        if need > ends[-1]:
            if need > 50.0 * kmax * k0:
                raise ValueError(f"a probe {abs(x_):.3g} m from the centre line needs the spectral "
                                 f"tail beyond {need / k0:.0f} k0; move it further out")
            more = np.arange(edges[-1], need + width, width)
            if len(more) > 1:
                mn = (0.5 * (more[1:] - more[:-1])[:, None] * (x + 1) + more[:-1, None]).ravel()
                mw = (0.5 * (more[1:] - more[:-1])[:, None] * w).ravel()
                extra = _probe_acc(p, mn, mw, k0, xs[q:q + 1], a)[:, 0]
                part = np.concatenate([part, part[-1] + np.cumsum(
                    extra.reshape(len(more) - 1, 16, -1).sum(axis=1), axis=0)])
                ends = np.concatenate([ends, more[1:]])

        def mean_at(K):
            k = int(np.searchsorted(ends, K - 1e-9 * K))
            return part[k - n + 1:k + 1].mean(axis=0), ends[k] - 0.5 * n * width

        full, kf = mean_at(ends[-1])
        half, kh = mean_at(0.5 * ends[-1])
        out[q] = base[q] + full + (full - half) / ((kf / kh) ** 2 - 1)
    return out[0] if np.ndim(xp) == 0 else out


def probe_reactance(f: float, eps_r: float, h: float, a: float) -> float:
    """The probe's own reactance, from the parallel plate: -(eta0 k0 h / 4) Y0(k a), k in the
    dielectric - the probe inductance of the cavity model. The spectral self-term of a
    probe ending on a patch diverges without an attachment mode; this stands in for it."""
    k0 = 2 * math.pi * f / C0
    return -(ETA0 * k0 * h / 4) * float(yv(0, k0 * math.sqrt(eps_r) * a))


def input_impedance(p: RectPatch, f: float, xp: float, a: float) -> complex:
    """Input impedance of a probe of radius a at (xp, 0) under the driven patch: the
    probe reactance plus the patch currents' reaction, -v^T Z^-1 v."""
    v = probe_vector(p, f, xp, a)
    er, h = (p.er1, p.h1) if isinstance(p, StackedPatch) else (p.er, p.h)
    zm = -np.einsum("...i,...i->...", v, np.linalg.solve(p.Z(f), v.T).T)
    return 1j * probe_reactance(f, er, h, a) + zm
