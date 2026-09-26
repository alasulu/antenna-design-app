"""A rectangular patch on a grounded slab by the spectral-domain method of
moments - a full-wave arbiter for the patch family.

`patch_q` radiates the cavity model's current through the slab; this module
does not assume the current. The patch (L along x, W along y, centred on the
origin) carries x- and y-directed currents expanded in entire-domain basis
functions,

    J_x = sin(n pi (x + L/2)/L) T_j(2y/W) / sqrt(1 - (2y/W)^2)
    J_y = sin(m pi (y + W/2)/W) T_j(2x/L) / sqrt(1 - (2x/L)^2)

(normal components vanish at the edges, parallel ones carry the edge
singularity), and Galerkin's method is applied with the exact spectral Green's
function of the grounded slab from its transverse equivalent network. The
spectral integrals run over a contour lifted above the real axis past the
surface-wave poles, then along it; the slowly converging reactance is
extrapolated from two truncations (Richardson, error ~ 1/k_max).

A patch's natural mode is taken as the characteristic mode (X v = lambda R v)
that radiates most, and its resonance as lambda = 0. From its current come the
broadside directivity (U(0) against the space-wave power), the radiation
efficiency (space wave against space plus surface wave) and the Q, from
d lambda / d omega.

It is checked on problems with independent answers: a vanishing current
element against Jackson and Alexopoulos's closed forms for the space and
surface waves of a dipole on a thin slab; an infinite microstrip line's
effective permittivity against Kirschning and Jansen; and the space-wave power
computed twice, from the spectral reaction and from the far field. SI units,
exp(+j omega t).
"""
from __future__ import annotations

import math

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
from scipy.special import jv

from ..core.constants import C0, EPS0, ETA0, MU0

__all__ = ["RectPatch", "BASIS_SMALL", "BASIS_FULL", "resonance", "mode_metrics", "slab_gf"]

BASIS_SMALL = dict(bx=((1, 0), (1, 2), (1, 4)), by=())
BASIS_FULL = dict(bx=((1, 0), (1, 2), (1, 4), (3, 0), (3, 2), (5, 0), (1, 6), (3, 4), (5, 2)),
                  by=((2, 1), (2, 3), (2, 5), (4, 1), (4, 3), (6, 1), (6, 3)))


def _sine_ft(k, a_len: float, n: int):
    """FT of sin(n pi (x + A/2)/A) on [-A/2, A/2], kernel exp(j k x)."""
    a = n * math.pi / a_len
    k = np.asarray(k, dtype=complex)
    k = np.where(np.abs(a * a - k * k) < 1e-9 * a * a, k + 1e-6 * a, k)
    return np.exp(-1j * k * a_len / 2) * a * (1 - (-1) ** n * np.exp(1j * k * a_len)) / (a * a - k * k)


def _cheb_ft(k, a_len: float, j: int):
    """FT of T_j(2x/A) / sqrt(1 - (2x/A)^2) on [-A/2, A/2]: (pi A/2) j^j J_j(k A/2)."""
    return 0.5 * math.pi * a_len * (1j ** j) * jv(j, np.asarray(k, dtype=complex) * a_len / 2)


def slab_gf(theta: float, eps_r: float, k0h: float):
    """Far-field factors of a horizontal current on the slab: TM (G) and TE (F)."""
    st, ct = math.sin(theta), math.cos(theta)
    N = np.sqrt(eps_r - st ** 2 + 0j)
    t = np.tan(k0h * N)
    return 2j * N * t * ct / (1j * N * t + eps_r * ct), 2j * t * ct / (1j * t * ct + N)


class RectPatch:
    def __init__(self, eps_r: float, h: float, L: float, W: float, bx=((1, 0),), by=()):
        self.er, self.h, self.L, self.W = eps_r, h, L, W
        self.basis = [("x", n, j) for n, j in bx] + [("y", m, j) for m, j in by]

    def ft(self, kx, ky):
        """(n_basis, 2, ...) spectra of the basis functions."""
        out = []
        for d, n, j in self.basis:
            if d == "x":
                f = _sine_ft(kx, self.L, n) * _cheb_ft(ky, self.W, j)
                out.append(np.stack([f, np.zeros_like(f)]))
            else:
                f = _sine_ft(ky, self.W, n) * _cheb_ft(kx, self.L, j)
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

    def _n_alpha(self, kr) -> int:
        # a product of two basis spectra has azimuthal harmonics up to ~ |kr| max(L, W)
        return int(4 * math.ceil((1.25 * abs(kr) * max(self.L, self.W) + 32) / 4))

    def _accumulate(self, nodes, weights, k0: float):
        nb = len(self.basis)
        Z = np.zeros((nb, nb), complex)
        for kr, wk in zip(nodes, weights):
            na = self._n_alpha(kr)
            al = 2 * math.pi * np.arange(na) / na
            ca, sa = np.cos(al), np.sin(al)
            Fp = self.ft(kr * ca, kr * sa)
            Fm = self.ft(-kr * ca, -kr * sa)
            ztm, zte = self.znode(kr, k0)
            gxx = ztm * ca ** 2 + zte * sa ** 2
            gyy = ztm * sa ** 2 + zte * ca ** 2
            gxy = (ztm - zte) * ca * sa
            Ex = gxx * Fp[:, 0] + gxy * Fp[:, 1]
            Ey = gxy * Fp[:, 0] + gyy * Fp[:, 1]
            Z += wk * kr * (Fm[:, 0] @ Ex.T + Fm[:, 1] @ Ey.T) * (2 * math.pi / na)
        return Z / (4 * math.pi ** 2)

    def Z_trunc(self, f: float, kmax: float, n_ell: int = 400, panels: int = 48):
        """Galerkin matrix with the spectral integral truncated at kmax * k0."""
        k0 = 2 * math.pi * f / C0
        K1 = k0 * (math.sqrt(self.er) + 1.0)          # past every surface-wave pole
        t, wt = np.polynomial.legendre.leggauss(n_ell)
        t = 0.5 * math.pi * (t + 1)
        wt = 0.5 * math.pi * wt
        b = 0.2 * k0
        nodes = list(0.5 * K1 * (1 - np.cos(t)) + 1j * b * np.sin(t))
        weights = list(wt * (0.5 * K1 * np.sin(t) + 1j * b * np.cos(t)))
        edges = K1 * np.geomspace(1, kmax * k0 / K1, panels + 1)
        x, w = np.polynomial.legendre.leggauss(32)
        for a, c in zip(edges[:-1], edges[1:]):
            nodes += list(0.5 * (c - a) * (x + 1) + a)
            weights += list(0.5 * (c - a) * w)
        return self._accumulate(np.array(nodes), np.array(weights), k0)

    def Z(self, f: float, kmax: float = 400.0):
        """Galerkin matrix, the truncation extrapolated from kmax and kmax/2."""
        return 2 * self.Z_trunc(f, kmax) - self.Z_trunc(f, kmax / 2)

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


def resonance(p: RectPatch, f_guess: float, span=(0.93, 1.07), step: float = 0.02) -> float:
    """Frequency at which the dominant characteristic mode's eigenvalue crosses zero."""
    g = lambda r: _dominant(p.Z(r * f_guess))[0]
    grid = np.arange(span[0], span[1] + 1e-9, step)
    vals = [g(r) for r in grid]
    for i in range(1, len(grid)):
        if vals[i - 1] < 0 <= vals[i]:
            return brentq(g, grid[i - 1], grid[i], xtol=2e-5) * f_guess
    raise ValueError("no resonance of the dominant mode in the span")


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
