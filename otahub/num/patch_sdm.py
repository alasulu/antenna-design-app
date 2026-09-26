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
from scipy.special import jv

from ..core.constants import C0, EPS0, ETA0, MU0

__all__ = ["RectPatch", "BASIS_SMALL", "BASIS_FULL", "resonance", "resonant_length", "mode_metrics", "slab_gf"]

# (n, j) for J_x, (m, j) for J_y; the TM10 class: J_x even-even, J_y odd-odd
BASIS_SMALL = dict(bx=((0, 0), (0, 2), (0, 4)), by=())
BASIS_FULL = dict(bx=((0, 0), (0, 2), (0, 4), (2, 0), (2, 2), (4, 0)), by=((1, 1), (1, 3), (3, 1)))


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
        # one quadrant suffices when every function is in the TM10 class
        self.quadrant = all((n % 2 == 0 and j % 2 == 0) if d == "x" else (n % 2 == 1 and j % 2 == 1)
                            for d, n, j in self.basis)

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
        cot = np.cos(x) / np.sin(x)
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
