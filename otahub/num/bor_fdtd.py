"""Body-of-revolution FDTD at azimuthal order m = 1 - an arbiter for
`waveguide_step` that shares nothing with it.

A field with cos(phi)/sin(phi) dependence reduces Maxwell's equations to six
components on a (rho, z) grid:

    E_rho = er cos(phi), E_phi = ep sin(phi), E_z = ez cos(phi)
    H_rho = hr sin(phi), H_phi = hp cos(phi), H_z = hz sin(phi)

staggered Yee-fashion with rho = i (cell 1, c = 1). On the axis ez and hz are
odd in rho, which fixes the two axis updates; walls of constant radius and the
step face sit on grid lines exactly, so there is no staircase. Both ends are
CPML. A TE11 pulse is launched in the smaller guide; transmitted TE11 and TM11
are read by projection onto the continuum mode profiles one wavelength past the
step, where the evanescent modes have died, and referred back to the step plane
with their own propagation constants. The incident wave comes from a run with no
step. At the ~1% level it agrees with mode matching on a 60-cell wavelength and
converges on it as the cell shrinks (order about 1.4 in the power share, set by
the field singularity at the step's corner).
"""
from __future__ import annotations

import math

import numpy as np
from scipy.special import jn_zeros, jnp_zeros, jv, jvp

__all__ = ["run", "project", "step_launch"]


def _cpml(n: int, npml: int, half: bool, dt: float, m: int = 3, alpha_max: float = 0.05):
    pos = np.arange(n) + (0.5 if half else 0.0)
    last = n if half else n - 1
    rho = np.maximum(np.clip((npml - pos) / npml, 0, 1), np.clip((pos - (last - npml)) / npml, 0, 1))
    sig = 0.8 * (m + 1) * rho ** m
    alp = alpha_max * (1 - rho)
    b = np.exp(-(sig + alp) * dt)
    c = np.where(sig > 0, sig * (b - 1) / (sig + alp), 0.0)
    return b, c


def _te11(radius: float, r):
    kc = jnp_zeros(1, 1)[0] / radius
    with np.errstate(invalid="ignore", divide="ignore"):
        jr = np.where(r > 0, jv(1, kc * r) / np.where(r > 0, r, 1), kc / 2)
    return -jr, kc * jvp(1, kc * r)


def run(nr_a: int, nr_b: int, zstep: int, nz: int, f: float, z_src: int, z_probe: int,
        step: bool = True, npml: int = 40, dt: float = 0.3, periods: float = 90.0, wall=None):
    """DFT at f (cycles per cell-time) of (er, ep) across the row z_probe.
    The guide has radius nr_a for k < zstep and nr_b after (nr_a throughout if not step),
    or follows `wall`, an integer radius for every plane k = 0..nz, if one is given: wherever
    it widens, the annulus between the two radii on that plane is a PEC face."""
    if wall is None:
        wall = np.where(np.arange(nz + 1) < zstep, nr_a, nr_b) if step else np.full(nz + 1, nr_a)
    wall = np.asarray(wall, dtype=int)
    nr_b = int(wall.max())
    NR = nr_b + 1
    er = np.zeros((NR, nz + 1)); ep = np.zeros((NR + 1, nz + 1)); ez = np.zeros((NR + 1, nz))
    hr = np.zeros((NR + 1, nz)); hp = np.zeros((NR, nz)); hz = np.zeros((NR, nz + 1))
    ri = np.arange(NR + 1, dtype=float)
    rh = np.arange(NR) + 0.5
    bi, ci = _cpml(nz + 1, npml, False, dt)
    bh, ch = _cpml(nz, npml, True, dt)
    psi_hr = np.zeros((NR + 1, nz)); psi_hp = np.zeros((NR, nz))
    psi_er = np.zeros((NR, nz + 1)); psi_ep = np.zeros((NR + 1, nz + 1))
    ep_mask = (ri[:, None] < wall[None, :]).astype(float)
    ez_mask = (ri[:, None] < np.minimum(wall[:-1], wall[1:])[None, :]).astype(float)
    er_mask = (rh[:, None] < wall[None, :]).astype(float)
    for k in np.flatnonzero(wall[1:] > wall[:-1]) + 1:   # every plane where the guide widens
        er_mask[:, k] = (rh < wall[k - 1]).astype(float)
        ep_mask[:, k] = (ri < wall[k - 1]).astype(float)
    s_er = _te11(nr_a, rh)[0] * (rh < nr_a)
    s_ep = _te11(nr_a, ri)[1] * (ri < nr_a)
    inv_ri = np.where(ri > 0, 1.0 / np.where(ri > 0, ri, 1), 0.0)[:, None]
    tau = 6.0 / f
    t0 = 3 * tau
    w = 2 * math.pi * f
    acc_r = np.zeros(NR, complex)
    acc_p = np.zeros(NR + 1, complex)
    for n in range(int((t0 + periods / f) / dt)):
        t = n * dt
        d = ep[:, 1:] - ep[:, :-1]
        psi_hr = bh * psi_hr + ch * d
        ezr = ez * inv_ri
        ezr[0, :] = ez[1, :]                          # ez odd in rho: ez/rho -> ez(1)/1 on the axis
        hr += dt * (ezr + d + psi_hr)
        d = er[:, 1:] - er[:, :-1]
        psi_hp = bh * psi_hp + ch * d
        hp -= dt * (d + psi_hp - (ez[1:, :] - ez[:-1, :]))
        rep = ri[:, None] * ep
        hz -= dt * ((rep[1:, :] - rep[:-1, :]) / rh[:, None] + er / rh[:, None])
        d = np.zeros((NR, nz + 1)); d[:, 1:-1] = hp[:, 1:] - hp[:, :-1]
        psi_er = bi * psi_er + ci * d
        er += dt * (hz / rh[:, None] - (d + psi_er))
        d = np.zeros((NR + 1, nz + 1)); d[:, 1:-1] = hr[:, 1:] - hr[:, :-1]
        psi_ep = bi * psi_ep + ci * d
        dhz = np.zeros((NR + 1, nz + 1)); dhz[1:NR, :] = hz[1:, :] - hz[:-1, :]
        dhz[0, :] = 2 * hz[0, :]                      # hz odd in rho
        ep += dt * (d + psi_ep - dhz)
        rhp = rh[:, None] * hp
        dr = np.zeros((NR + 1, nz)); dr[1:NR, :] = rhp[1:, :] - rhp[:-1, :]
        ez += dt * (dr - hr) * inv_ri
        ez[0, :] = 0.0
        g = math.exp(-((t - t0) / tau) ** 2) * math.sin(w * (t - t0))
        er[:, z_src] += g * s_er
        ep[:, z_src] += g * s_ep
        er *= er_mask; ep *= ep_mask; ez *= ez_mask
        ph = complex(math.cos(w * t), -math.sin(w * t))
        acc_r += er[:, z_probe] * ph
        acc_p += ep[:, z_probe] * ph
    return acc_r, acc_p


def project(er_row, ep_row, radius: float, kind: str, kc: float) -> complex:
    """Amplitude of a mode from a sampled cross-section (er at half radii, ep at integer)."""
    rh = np.arange(len(er_row)) + 0.5
    ri = np.arange(len(ep_row), dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        jr_h = jv(1, kc * rh) / rh
        jr_i = np.where(ri > 0, jv(1, kc * ri) / np.where(ri > 0, ri, 1), kc / 2)
    if kind == "TE":
        mr, mp = -jr_h, kc * jvp(1, kc * ri)
    else:
        mr, mp = -kc * jvp(1, kc * rh), jr_i
    wr = rh * (rh < radius)
    wp = (ri * (ri < radius)).astype(float)
    wp[np.isclose(ri, radius)] = 0.5 * radius
    return complex((np.sum(er_row * mr * wr) + np.sum(ep_row * mp * wp)) / (np.sum(mr * mr * wr) + np.sum(mp * mp * wp)))


def step_launch(nr_a: int, nr_b: int, cells_per_wavelength: int) -> tuple[float, float]:
    """(TM11 share of the transmitted power, TM11 phase relative to TE11 at the
    step, degrees) for guide radii in cells."""
    N = cells_per_wavelength
    f = 1.0 / N
    nz, zstep, zsrc = int(8.5 * N), int(3.7 * N), int(1.3 * N)
    npml = int(40 * N / 60)
    ref = run(nr_a, nr_b, zstep, nz, f, zsrc, zstep, step=False, npml=npml)
    a_inc = project(*ref, nr_a, "TE", jnp_zeros(1, 1)[0] / nr_a)
    out = run(nr_a, nr_b, zstep, nz, f, zsrc, zstep + N, step=True, npml=npml)
    k0 = 2 * math.pi * f
    kte, ktm = jnp_zeros(1, 1)[0] / nr_b, jn_zeros(1, 1)[0] / nr_b
    bte, btm = math.sqrt(k0 ** 2 - kte ** 2), math.sqrt(k0 ** 2 - ktm ** 2)
    t_te = project(*out, nr_b, "TE", kte) * np.exp(1j * bte * N) / a_inc
    t_tm = project(*out, nr_b, "TM", ktm) * np.exp(1j * btm * N) / a_inc
    # power share with the continuum mode norms and wave admittances (eta0 = 1 units)
    x = np.polynomial.legendre.leggauss(400)
    r = 0.5 * nr_b * (x[0] + 1)
    wq = 0.5 * nr_b * x[1]
    def norm(kind, kc):
        with np.errstate(invalid="ignore", divide="ignore"):
            jr = np.where(r > 0, jv(1, kc * r) / r, kc / 2)
        er_, ep_ = (-jr, kc * jvp(1, kc * r)) if kind == "TE" else (-kc * jvp(1, kc * r), jr)
        return float(np.sum((er_ ** 2 + ep_ ** 2) * r * wq))
    p_te = abs(t_te) ** 2 * norm("TE", kte) * bte / k0
    p_tm = abs(t_tm) ** 2 * norm("TM", ktm) * k0 / btm
    return p_tm / (p_te + p_tm), math.degrees(np.angle(t_tm / t_te))
