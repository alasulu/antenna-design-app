"""A step between two coaxial circular waveguides by mode matching - an arbiter
for the dual-mode (Potter) horn.

TE11 arrives from the smaller guide (radius a) at the plane where it widens to
radius b. Both sides are expanded in the m = 1 modes that the symmetry allows,
TE1n and TM1n, all x-polarised at the axis:

    TE: e_t = z x grad(J1(kc r) sin(phi)),  kc = j'_1n / radius
    TM: e_t = -grad(J1(kc r) cos(phi)),     kc = j_1n / radius

Transverse E is matched over the larger cross-section (zero on the step face)
and transverse H over the aperture of the smaller guide, which gives the
reflection R and transmission T of every mode. The number of modes on each side
is kept in proportion to the radii (relative convergence). The coupling
integrals are done by quadrature. Power balances to 1e-8 by construction of the
method, which is a check on the integrals, not on the method; the method is
checked against `bor_fdtd`, an FDTD of the same step.

What the horn needs from it: the share of the transmitted power carried by TM11
and the phase of TM11 relative to TE11 at the step. SI units, exp(+j omega t).
"""
from __future__ import annotations

import math

import numpy as np
from scipy.special import jn_zeros, jnp_zeros, jv, jvp

from ..core.constants import C0, EPS0, MU0

__all__ = ["modes", "step", "tm11_launch", "cascade", "aperture_share"]


def modes(radius: float, n_te: int, n_tm: int):
    """[(kind, kc)] for TE1n then TM1n; TE11 is first."""
    return [("TE", x / radius) for x in jnp_zeros(1, n_te)] + [("TM", x / radius) for x in jn_zeros(1, n_tm)]


def _radial(kind: str, kc: float, r):
    """(e_rho, e_phi) radial profiles; e_rho multiplies cos(phi), e_phi sin(phi)."""
    x = kc * r
    with np.errstate(invalid="ignore", divide="ignore"):
        jr = np.where(r > 0, jv(1, x) / np.where(r > 0, r, 1), kc / 2)
    if kind == "TE":
        return -jr, kc * jvp(1, x)
    return -kc * jvp(1, x), jr


def _profiles(mode_list, r):
    """(n_modes, 2, len(r)) radial profiles."""
    return np.array([np.stack(_radial(k, kc, r)) for k, kc in mode_list])


def _gram(A, B, a: float, n: int = 400):
    """pi * int_0^a (e_rho e_rho' + e_phi e_phi') r dr for every pair of A and B."""
    x, w = np.polynomial.legendre.leggauss(n)
    r = 0.5 * a * (x + 1)
    wr = 0.5 * a * w * r
    PA, PB = _profiles(A, r), _profiles(B, r)
    return math.pi * (np.einsum("ict,jct,t->ij", PA, PB, wr))


def _admittance(kind: str, kc: float, k0: float) -> complex:
    beta = complex(np.sqrt(k0 * k0 - kc * kc + 0j))
    if beta.imag > 0:
        beta = -beta                                   # evanescent: decaying for exp(-j beta z)
    w = k0 * C0
    return beta / (w * MU0) if kind == "TE" else w * EPS0 / beta


def step(a: float, b: float, f: float, n_modes: int = 12):
    """Scattering of TE11 incident from radius a into radius b > a at frequency f.

    Returns (R, T, modes_a, modes_b, powers) with powers = dict(incident, reflected,
    transmitted) per mode, in the method's common units."""
    k0 = 2 * math.pi * f / C0
    s = b / a
    M1 = modes(a, n_modes, n_modes)
    M2 = modes(b, max(1, round(n_modes * s)), max(1, round(n_modes * s)))
    X = _gram(M1, M2, a)
    N1 = np.diag(_gram(M1, M1, a)).copy()
    N2 = np.diag(_gram(M2, M2, b)).copy()
    Y1 = np.array([_admittance(*p, k0) for p in M1])
    Y2 = np.array([_admittance(*q, k0) for q in M2])
    e = np.zeros(len(M1))
    e[0] = 1.0
    D1 = np.diag(N1 * Y1)
    M = X @ np.diag(Y2 / N2) @ X.T
    R = np.linalg.solve(D1 + M, (D1 - M) @ e)
    T = (X.T @ (e + R)) / N2
    powers = dict(incident=N1[0] * Y1[0].real, reflected=np.abs(R) ** 2 * N1 * Y1.real,
                  transmitted=np.abs(T) ** 2 * N2 * Y2.real)
    return R, T, M1, M2, powers


def tm11_launch(d_in: float, d_step: float, n_modes: int = 12) -> tuple[float, float, float]:
    """(TM11 share of the transmitted power, TM11 phase relative to TE11 at the
    step in degrees, TE11 power reflected) for guide diameters in wavelengths."""
    R, T, M1, M2, P = step(0.5 * d_in, 0.5 * d_step, C0, n_modes)
    i = next(k for k, m in enumerate(M2) if m[0] == "TM")
    pt = P["transmitted"]
    return (float(pt[i] / (pt[0] + pt[i])), float(math.degrees(np.angle(T[i] / T[0]))),
            float(P["reflected"][0] / P["incident"]))


# ---------------------------------------------------------------- cascades: a stepped or flared horn

def _section(radius: float, per: float, k0: float):
    n = max(3, int(round(per * radius)))
    M = modes(radius, n, n)
    Y = np.array([_admittance(*m, k0) for m in M])
    N = np.diag(_gram(M, M, radius)).copy()
    beta = np.array([complex(np.sqrt(k0 * k0 - kc * kc + 0j)) for _, kc in M])
    beta = np.where(beta.imag > 0, -beta, beta)
    return M, Y, N, beta


def _gsm_step(a, b, s1, s2):
    """Generalised scattering matrix of a widening from radius a to b; port 1 is the narrow side."""
    M1, Y1, N1, _ = s1
    M2, Y2, N2, _ = s2
    X = _gram(M1, M2, a)
    P = (X / N2[None, :]).T
    Q = (X * Y2[None, :]) / (N1 * Y1)[:, None]
    I1, I2 = np.eye(len(M1)), np.eye(len(M2))
    inv = np.linalg.inv(I1 + Q @ P)
    S11 = inv @ (I1 - Q @ P)
    S12 = 2 * inv @ Q
    S21 = P @ (I1 + S11)
    return S11, S12, S21, P @ S12 - I2


def _star(A, B):
    A11, A12, A21, A22 = A
    B11, B12, B21, B22 = B
    I = np.eye(A22.shape[0])
    X1 = np.linalg.inv(I - B11 @ A22)
    X2 = np.linalg.inv(I - A22 @ B11)
    return A11 + A12 @ X1 @ B11 @ A21, A12 @ X1 @ B12, B21 @ X2 @ A21, B22 + B21 @ X2 @ A22 @ B12


def cascade(profile, f: float, per: float = 18.0):
    """A horn as a staircase: profile = [(radius, length), ...] with non-decreasing radius
    (lengths and radii in metres). Returns (S, first section, last section) where S is the
    generalised scattering matrix (S11, S12, S21, S22) between the first and last sections'
    modes, each section = (modes, admittances, norms, betas). The number of modes grows in
    proportion to the radius (per modes of each family per metre of radius at 1 m wavelength)."""
    k0 = 2 * math.pi * f / C0
    per_m = per * k0 / (2 * math.pi)
    S = first = prev = None
    prev_a = None
    for a, length in profile:
        sec = _section(a, per_m, k0)
        d = np.diag(np.exp(-1j * sec[3] * length))
        z = np.zeros_like(d)
        if prev is None:
            first = sec
            S = (z, d, d, z)
        else:
            if a > prev_a * (1 + 1e-12):
                S = _star(S, _gsm_step(prev_a, a, prev, sec))
            S = _star(S, (z, d, d, z))
        prev, prev_a = sec, a
    return S, first, prev


def aperture_share(b, section, radius: float, apex_distance: float, f: float, n: int = 800):
    """TM11 share and TM11-to-TE11 phase (degrees) of the aperture field built from the mode
    amplitudes b of the last section, after removing a spherical phase front centred
    apex_distance behind the aperture - the form the aperture model (horn_pattern) takes."""
    k0 = 2 * math.pi * f / C0
    M = section[0]
    x, w = np.polynomial.legendre.leggauss(n)
    r = 0.5 * radius * (x + 1)
    w = 0.5 * radius * w
    field = np.tensordot(b, _profiles(M, r), axes=(0, 0))
    field = field * np.exp(1j * k0 * (np.sqrt(apex_distance ** 2 + r ** 2) - apex_distance))
    out = {}
    for kind, chi in (("TE", 1.8411837813), ("TM", 3.8317059702)):
        m = _profiles([(kind, chi / radius)], r)[0]
        den = np.sum((m[0] ** 2 + m[1] ** 2) * r * w)
        amp = np.sum((field[0] * m[0] + field[1] * m[1]) * r * w) / den
        out[kind] = (amp, abs(amp) ** 2 * math.pi * den * _admittance(kind, chi / radius, k0).real)
    return (float(out["TM"][1] / (out["TE"][1] + out["TM"][1])),
            float(math.degrees(np.angle(out["TM"][0] / out["TE"][0]))))
