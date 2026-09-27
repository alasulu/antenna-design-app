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

The launched share is not the aperture's. TM11 is cut off in the input guide, so
the step reflects it totally, and the start of the flare - where TM11, barely
above cutoff, sees its impedance change fast - reflects part of it: the phasing
guide is a TM11 resonator. `potter_aperture` runs the whole chain; on a
10-wavelength horn from a 1.1-wavelength guide the share arriving in phase at the
aperture swings from 0.001 to 0.24 as the step goes from 1.25 to 1.55
wavelengths, not monotonically, and an FDTD of the same staircase follows the
cascade to a degree in phase.
"""
from __future__ import annotations

import math
from functools import lru_cache

import numpy as np
from scipy.special import j0, j1, jn_zeros, jnp_zeros

from ..core.constants import C0, EPS0, MU0

__all__ = ["modes", "step", "tm11_launch", "cascade", "aperture_share", "potter_profile", "potter_aperture",
           "PotterChain", "potter_design"]


def modes(radius: float, n_te: int, n_tm: int):
    """[(kind, kc)] for TE1n then TM1n; TE11 is first."""
    return [("TE", x / radius) for x in jnp_zeros(1, n_te)] + [("TM", x / radius) for x in jn_zeros(1, n_tm)]


def _radial(kind: str, kc: float, r):
    """(e_rho, e_phi) radial profiles; e_rho multiplies cos(phi), e_phi sin(phi)."""
    p = _profiles([(kind, kc)], np.asarray(r, dtype=float))[0]
    return p[0], p[1]


def _profiles(mode_list, r):
    """(n_modes, 2, len(r)) radial profiles, all modes at once.

    J1(x)/r and kc J1'(x), with J1' = J0 - J1/x (the recurrence, exact, and far
    cheaper than a general derivative formula); both tend to kc/2 on the axis."""
    r = np.asarray(r, dtype=float)
    kc = np.array([k for _, k in mode_list], dtype=float)[:, None]
    te = np.array([kind == "TE" for kind, _ in mode_list])[:, None]
    x = kc * r[None, :]
    safe = np.where(x > 0, x, 1.0)
    jr = np.where(x > 0, kc * j1(safe) / safe, kc / 2)          # J1(kc r)/r
    dj = np.where(x > 0, kc * (j0(safe) - j1(safe) / safe), kc / 2)   # kc J1'(kc r)
    e_rho = np.where(te, -jr, -dj)
    e_phi = np.where(te, dj, jr)
    return np.stack([e_rho, e_phi], axis=1)


@lru_cache(maxsize=None)
def _gauss(n: int):
    return np.polynomial.legendre.leggauss(n)


def _gram(A, B, a: float, n: int = 400):
    """pi * int_0^a (e_rho e_rho' + e_phi e_phi') r dr for every pair of A and B."""
    x, w = _gauss(n)
    r = 0.5 * a * (x + 1)
    wr = 0.5 * a * w * r
    PA, PB = _profiles(A, r), _profiles(B, r)
    nt = 2 * len(r)
    return math.pi * ((PA * wr).reshape(len(A), nt) @ PB.reshape(len(B), nt).T)


def _norms(A, a: float, n: int = 400):
    """The diagonal of _gram(A, A, a), without the rest of it."""
    x, w = _gauss(n)
    r = 0.5 * a * (x + 1)
    PA = _profiles(A, r)
    return math.pi * np.sum(PA * PA * (0.5 * a * w * r), axis=(1, 2))


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
    N = _norms(M, radius)
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
    x, w = _gauss(n)
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


def potter_profile(d_in: float, d_step: float, ell: float, d_ap: float, slant: float, stair: float = 0.05,
                   lead: float = 0.3):
    """A Potter horn as a staircase for `cascade`, in wavelengths (so metres at f = c):
    `lead` of the input guide, the step to d_step, `ell` of phasing guide, then a cone
    from d_step to the aperture d_ap whose slant length from its apex is `slant`, cut into
    stairs of axial length at most `stair` (each at the cone's radius mid-stair).
    Returns (profile, axial distance from the apex to the aperture)."""
    a1, a2 = d_step / 2, d_ap / 2
    tan_t = (a2 / slant) / math.sqrt(1 - (a2 / slant) ** 2)
    lf = (a2 - a1) / tan_t
    n = max(4, int(math.ceil(lf / stair)))
    cone = [((a1 + (i - 0.5) * lf / n * tan_t) if i < n else a2, lf / n) for i in range(1, n + 1)]
    return [(d_in / 2, lead), (a1, ell)] + cone, math.sqrt(slant ** 2 - a2 ** 2)


def potter_aperture(d_in: float, d_step: float, ell: float, d_ap: float, slant: float, f_ratio: float = 1.0,
                    stair: float = 0.05, per: float = 18.0) -> tuple[float, float]:
    """(TM11 share, TM11-to-TE11 phase in degrees) at the aperture of that horn, TE11 in,
    at f_ratio times the frequency at which the dimensions are in wavelengths - the
    whole chain at once: step, phasing guide and flare, and the TM11 that rattles between
    the step (where it is cut off) and the start of the flare."""
    prof, apex = potter_profile(d_in, d_step, ell, d_ap, slant, stair)
    S, _, last = cascade(prof, C0 * f_ratio, per)
    return aperture_share(S[2][:, 0], last, d_ap / 2, apex, C0 * f_ratio)


# ---------------------------------------------------------------- designing a Potter horn jointly

class PotterChain:
    """The whole chain at one step diameter, with the phasing length left free.

    Everything before the phasing guide (input guide and step) and everything after
    it (the start of the flare and the cone) are cascaded once; the phasing guide
    between them is only a phase per mode, so a new phasing length costs two small
    star products instead of a whole cascade. Dimensions in wavelengths at f_ratio = 1.
    """

    def __init__(self, d_in: float, d_step: float, d_ap: float, slant: float, f_ratio: float = 1.0,
                 stair: float = 0.05, per: float = 18.0, lead: float = 0.3):
        prof, self.apex = potter_profile(d_in, d_step, 0.0, d_ap, slant, stair, lead)
        self.f = C0 * f_ratio
        self.radius = d_ap / 2
        self.A = cascade(prof[:2], self.f, per)[0]
        self.B, first, self.last = cascade([(d_step / 2, 0.0)] + prof[2:], self.f, per)
        self.beta = first[3]
        M = first[0]
        i_tm = next(k for k, m in enumerate(M) if m[0] == "TM")
        self.beat = 2 * math.pi / (self.beta[0].real - self.beta[i_tm].real)   # TE11-TM11 beat length [m]

    def at(self, ell: float) -> tuple[float, float]:
        """(TM11 share, TM11-to-TE11 phase in degrees) at the aperture for phasing length ell."""
        d = np.diag(np.exp(-1j * self.beta * ell))
        z = np.zeros_like(d)
        S = _star(_star(self.A, (z, d, d, z)), self.B)
        return aperture_share(S[2][:, 0], self.last, self.radius, self.apex, self.f)

    def in_phase(self, phase_deg: float = 0.0, n: int = 48, periods: int = 1) -> list[tuple[float, float]]:
        """Every phasing length in `periods` beat periods [0, periods * beat) that puts
        TM11 at phase_deg, with the share it delivers there: [(ell, share), ...]. A guide
        one beat longer is in phase again but delivers another share - the TM11 that
        rattles in it has gone round a different phase - so each period is its own design."""
        from scipy.optimize import brentq

        def err(e):
            return (self.at(e)[1] - phase_deg + 180.0) % 360.0 - 180.0

        es = np.linspace(0.0, periods * self.beat, periods * n + 1)
        v = [err(e) for e in es]
        out = []
        for i in range(periods * n):
            if v[i] <= 0.0 < v[i + 1] or v[i] >= 0.0 > v[i + 1]:
                if abs(v[i + 1] - v[i]) < 180.0:                        # a crossing, not the wrap
                    e = brentq(err, es[i], es[i + 1], xtol=1e-9)
                    out.append((e, self.at(e)[0]))
        return out


def potter_design(d_in: float, d_ap: float, slant: float, share: float = 0.15, phase_deg: float = 0.0,
                  steps=None, stair: float = 0.05, per: float = 18.0, slope: bool = True,
                  periods: int = 1, progress=None) -> list[dict]:
    """Step diameter and phasing length that deliver `share` of TM11 at `phase_deg` at the
    aperture, the whole chain solved together. Dimensions in wavelengths.

    The share that arrives in phase is not monotonic in the step (the phasing guide is a
    TM11 resonator), so there is not one answer: the step is scanned over `steps`
    (default 1.23 to 1.70 wavelengths, 0.01 apart), every in-phase phasing length within a
    beat period (or `periods` of them: a guide a beat longer is another design, usually a
    narrower-band one) is followed from one step to the next as a branch, and each crossing
    of the target share is refined. Returns every solution found, each with the aperture
    phase's slope against frequency (degrees per percent, the number that sets the
    cross-polar bandwidth) when `slope`, flattest first.
    """
    if steps is None:
        steps = np.round(np.arange(max(1.23, d_in + 0.02), min(1.70, d_ap - 0.05) + 1e-9, 0.01), 6)
    note = progress or (lambda *_: None)
    rows = []
    for ds in steps:
        ch = PotterChain(d_in, ds, d_ap, slant, 1.0, stair, per)
        rows.append((float(ds), ch.beat, ch.in_phase(phase_deg, periods=periods)))
        note(f"step {ds:.3f}: " + ", ".join(f"ell {e:.3f} share {p:.3f}" for e, p in rows[-1][2]))

    # link roots into branches by continuity of the phasing length (modulo the beat)
    branches: list[list[tuple[float, float, float]]] = []
    last_index: list[int] = []
    for i, (ds, beat, roots) in enumerate(rows):
        used = set()
        for b, br in enumerate(branches):
            if last_index[b] != i - 1:
                continue
            _, pe, _ = br[-1]
            best = None
            for k, (e, p) in enumerate(roots):
                if k in used:
                    continue
                gap = abs(e - pe)
                if gap < 0.3 and (best is None or gap < best[0]):
                    best = (gap, k)
            if best is not None:
                used.add(best[1])
                br.append((ds, roots[best[1]][0], roots[best[1]][1]))
                last_index[b] = i
        for k, (e, p) in enumerate(roots):
            if k not in used:
                branches.append([(ds, e, p)])
                last_index.append(i)

    from scipy.optimize import brentq
    out = []
    for br in branches:
        for (d0, e0, p0), (d1, e1, p1) in zip(br, br[1:]):
            if (p0 - share) * (p1 - share) > 0:
                continue
            cache = {}

            def share_err(ds, e_guess=e0):
                ch = PotterChain(d_in, ds, d_ap, slant, 1.0, stair, per)
                roots = ch.in_phase(phase_deg, periods=periods)
                if not roots:
                    raise ValueError("branch lost")
                e, p = min(roots, key=lambda r: abs(r[0] - e_guess))
                cache[ds] = (e, p)
                return p - share

            try:
                ds = brentq(share_err, d0, d1, xtol=1e-5)
            except ValueError:
                continue
            e, p = cache[ds]
            sol = dict(d_step=ds, ell=e, share=p, phase_deg=phase_deg, launched_share=tm11_launch(d_in, ds)[0])
            if slope:
                ph = [PotterChain(d_in, ds, d_ap, slant, fr, stair, per).at(e)[1] for fr in (0.99, 1.01)]
                sol["phase_slope_deg_per_percent"] = ((ph[1] - ph[0] + 180.0) % 360.0 - 180.0) / 2
            note(f"solution: step {ds:.4f} ell {e:.4f} share {p:.4f}"
                 + (f" slope {sol['phase_slope_deg_per_percent']:+.2f} deg/%" if slope else ""))
            out.append(sol)
    if slope:
        out.sort(key=lambda r: abs(r["phase_slope_deg_per_percent"]))
    return out

