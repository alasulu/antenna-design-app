"""Radiation Q of a microstrip patch on its grounded substrate - an arbiter for
the patch family's bandwidths.

Q = omega * W / P. The stored energy W comes from the cavity mode under the
patch, E_z = V0 * psi(x, y) / h. The radiated power P does NOT come from the
cavity model's edge slots radiating as if the substrate were not there; it is
the space wave of the patch's own surface current J = grad(psi) * V0 / (j omega
mu0 h) lying on a grounded slab of thickness h. Each plane-wave direction sees
the slab through the transverse-equivalent-network factors (reciprocity):

    TM:  G(theta) = 2j N t cos(theta) / (j N t + eps_r cos(theta))
    TE:  F(theta) = 2j t cos(theta) / (j t cos(theta) + N)

with N = sqrt(eps_r - sin^2 theta) and t = tan(k0 h N); in free space G would
be cos(theta) and F would be 1. This is the physics Jackson and Alexopoulos's
closed-form patch bandwidth was derived from, before its thin-substrate
expansion - so for a rectangular patch that formula is the check, and for other
shapes this is the derivation. Surface waves are not included: this is the
radiation Q, as the bandwidth formulas it is compared with are.

A shape is described by quadrature points inside it with psi and grad(psi) at
each; `rectangle`, `disc` and `triangle` build the dominant modes. SI units.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.special import j1, jvp

from ..core.constants import C0, EPS0, ETA0, MU0

__all__ = ["radiation_q", "rectangle", "disc", "triangle", "jackson_q"]


def radiation_q(pts, w, psi, grad, eps_r: float, h: float, f: float,
                n_theta: int = 96, n_phi: int = 192) -> float:
    """Radiation Q of the mode sampled at `pts` (N, 2) with weights `w`."""
    omega = 2 * math.pi * f
    k0 = omega / C0
    energy = 0.5 * EPS0 * eps_r / h * float(np.sum(psi ** 2 * w))        # 2 * W_e, V0 = 1
    J = grad / (1j * omega * MU0 * h)                                    # (N, 2)
    xt, wt = np.polynomial.legendre.leggauss(n_theta)
    th = 0.25 * math.pi * (xt + 1.0)                                     # [0, pi/2]
    wth = 0.25 * math.pi * wt
    ph = 2 * math.pi * np.arange(n_phi) / n_phi
    total = 0.0
    for i, t_ in enumerate(th):
        st, ct = math.sin(t_), math.cos(t_)
        N = np.sqrt(eps_r - st ** 2 + 0j)
        tt = np.tan(k0 * h * N)
        G = 2j * N * tt * ct / (1j * N * tt + eps_r * ct)
        F = 2j * tt * ct / (1j * tt * ct + N)
        kx, ky = k0 * st * np.cos(ph), k0 * st * np.sin(ph)
        phase = np.exp(1j * (np.outer(kx, pts[:, 0]) + np.outer(ky, pts[:, 1])))   # (P, N)
        Jx = phase @ (J[:, 0] * w)
        Jy = phase @ (J[:, 1] * w)
        Jr = Jx * np.cos(ph) + Jy * np.sin(ph)
        Jp = -Jx * np.sin(ph) + Jy * np.cos(ph)
        dens = abs(G) ** 2 * np.abs(Jr) ** 2 + abs(F) ** 2 * np.abs(Jp) ** 2
        total += wth[i] * st * float(np.sum(dens)) * (2 * math.pi / n_phi)
    power = ETA0 * k0 ** 2 / (32 * math.pi ** 2) * total
    return omega * energy / power


def _gauss_box(a: float, b: float, n: int, m: int):
    xa, wa = np.polynomial.legendre.leggauss(n)
    xb, wb = np.polynomial.legendre.leggauss(m)
    x = 0.5 * a * (xa + 1)
    y = 0.5 * b * (xb + 1)
    X, Y = np.meshgrid(x, y, indexing="ij")
    W = np.outer(0.5 * a * wa, 0.5 * b * wb)
    return np.column_stack([X.ravel(), Y.ravel()]), W.ravel()


def rectangle(L: float, W: float, n: int = 40):
    """TM10: psi = cos(pi x / L) on 0 < x < L, 0 < y < W."""
    pts, w = _gauss_box(L, W, n, 16)
    kx = math.pi / L
    psi = np.cos(kx * pts[:, 0])
    grad = np.column_stack([-kx * np.sin(kx * pts[:, 0]), np.zeros(len(pts))])
    return pts, w, psi, grad


def disc(a: float, n: int = 48):
    """TM11: psi = J1(k rho) cos(phi), k a = 1.8412 (J1' = 0 at the rim)."""
    k = 1.8411837813 / a
    xr, wr = np.polynomial.legendre.leggauss(n)
    rho = 0.5 * a * (xr + 1)
    phi = 2 * math.pi * np.arange(2 * n) / (2 * n)
    R, P = np.meshgrid(rho, phi, indexing="ij")
    Wt = np.outer(0.5 * a * wr * rho, np.full(2 * n, 2 * math.pi / (2 * n)))
    pts = np.column_stack([(R * np.cos(P)).ravel(), (R * np.sin(P)).ravel()])
    kr = k * R
    psi = (j1(kr) * np.cos(P)).ravel()
    dr = (k * jvp(1, kr) * np.cos(P)).ravel()
    with np.errstate(invalid="ignore", divide="ignore"):
        dp = np.where(R > 0, -j1(kr) * np.sin(P) / R, 0.0).ravel()
    c, s = np.cos(P).ravel(), np.sin(P).ravel()
    grad = np.column_stack([dr * c - dp * s, dr * s + dp * c])
    return pts, Wt.ravel(), psi, grad


def triangle(a: float, n: int = 60):
    """Dominant mode of an equilateral triangle of side a, vertices (0,0), (a,0)
    and (a/2, a*sqrt(3)/2). It is not quoted: six plane waves of magnitude
    k = 4 pi / (3a) on the hexagonal star are combined through the null space of
    the magnetic-wall condition d psi / dn = 0 sampled along all three walls.
    The mode is a degenerate pair; by the triangle's symmetry every member has
    the same stored energy and radiated power, so any one gives the Q."""
    k = 4 * math.pi / (3 * a)
    kv = k * np.column_stack([np.cos(np.deg2rad(np.arange(0, 360, 60))),
                              np.sin(np.deg2rad(np.arange(0, 360, 60)))])
    v = [np.array([0.0, 0.0]), np.array([a, 0.0]), np.array([a / 2, a * math.sqrt(3) / 2])]
    cen = sum(v) / 3
    rows = []
    for p0, q0 in ((v[0], v[1]), (v[1], v[2]), (v[2], v[0])):
        d = q0 - p0
        nrm = np.array([d[1], -d[0]]) / np.linalg.norm(d)
        if nrm @ (cen - (p0 + q0) / 2) > 0:
            nrm = -nrm
        for t in np.linspace(0.01, 0.99, 50):
            rows.append(1j * (kv @ nrm) * np.exp(1j * (kv @ (p0 + t * d))))
    c = np.linalg.svd(np.array(rows))[2][-1].conj()
    # Duffy map of the unit square onto the triangle
    xg, wg = np.polynomial.legendre.leggauss(n)
    u = 0.5 * (xg + 1)
    U, V = np.meshgrid(u, u, indexing="ij")
    e1, e2 = v[1] - v[0], v[2] - v[1]
    X = U * e1[0] + U * V * e2[0]
    Y = U * e1[1] + U * V * e2[1]
    jac = abs(e1[0] * e2[1] - e1[1] * e2[0])
    Wt = (np.outer(0.5 * wg, 0.5 * wg) * U * jac).ravel()
    pts = np.column_stack([X.ravel(), Y.ravel()])
    ph = np.exp(1j * (pts @ kv.T))
    psi = (ph @ c).real
    grad = np.column_stack([(ph @ (1j * kv[:, 0] * c)).real, (ph @ (1j * kv[:, 1] * c)).real])
    scale = np.abs(psi).max()
    return pts, Wt, psi / scale, grad / scale


def jackson_q(L: float, W: float, eps_r: float, h: float, f: float) -> float:
    """Jackson & Alexopoulos (1991) radiation Q of a rectangular patch, full
    form: Q = (3/16) (eps_r / (p c1)) (lambda0 / h) (L / W)."""
    k0 = 2 * math.pi * f / C0
    c1 = 1 - 1 / eps_r + 2 / (5 * eps_r ** 2)
    p = 1 - 0.16605 / 20 * (k0 * W) ** 2 + 0.00761 * 3 / 560 * (k0 * W) ** 4 - 0.0914153 / 10 * (k0 * L) ** 2
    return 3 / 16 * eps_r / (p * c1) * (C0 / f / h) * (L / W)


def directivity(pts, w, psi, grad, eps_r: float, h: float, f: float,
                n_theta: int = 96, n_phi: int = 192) -> float:
    """Broadside directivity of the same space wave: 4 pi U(0) / P, a check on
    the radiation integral against the cavity model's edge-slot patterns."""
    omega = 2 * math.pi * f
    k0 = omega / C0
    J = grad / (1j * omega * MU0 * h)
    N0 = math.sqrt(eps_r)
    t0 = math.tan(k0 * h * N0)
    G0 = 2j * N0 * t0 / (1j * N0 * t0 + eps_r)
    F0 = 2j * t0 / (1j * t0 + N0)
    Jx0, Jy0 = np.sum(J[:, 0] * w), np.sum(J[:, 1] * w)
    u0 = ETA0 * k0 ** 2 / (32 * math.pi ** 2) * (abs(G0) ** 2 * abs(Jx0) ** 2 + abs(F0) ** 2 * abs(Jy0) ** 2)
    q = radiation_q(pts, w, psi, grad, eps_r, h, f, n_theta, n_phi)
    energy = 0.5 * EPS0 * eps_r / h * float(np.sum(psi ** 2 * w))
    power = omega * energy / q
    return 4 * math.pi * u0 / power
