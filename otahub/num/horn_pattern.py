"""Horn beamwidths by aperture integration - an arbiter for the horn family.

Every horn here is an aperture field with a quadratic phase error across it,
exp(-j 2 pi s (r/r_edge)^2), s being the edge phase error in wavelengths. The
far field is the aperture's Fourier transform times the Huygens obliquity
factor (1 + cos theta)/2; the half-power beamwidth is found by root-finding on
that pattern. Four apertures, in units of the free-space wavelength:

``slit_uniform``   E plane of an E-plane sectoral horn (uniform across b1)
``slit_cosine``    H plane of an H-plane sectoral horn (cos(pi x/a1))
``te11``           a smooth conical horn, E or H plane, from the circular
                   guide's TE11 field (J1(v)/v and J1'(v) mixed by azimuth)
``he11``           a corrugated horn's balanced hybrid mode, J0(2.405 r/a),
                   the same in every plane

Two-mode apertures transform both field components, which gives the Ludwig-3
co- and cross-polar patterns from the same integral:

``dual_mode``      Potter's horn: TE11 plus TM11 in phase at the aperture, a
                   given fraction of the power in TM11 (``disc_pattern``,
                   ``disc_efficiency``)
``diagonal``       the diagonal horn, y cos(pi x/a) + x cos(pi y/a) on a square
                   (``diagonal_efficiency``)
``peak_cross``     the highest cross-polar level in a plane

The integrals are done by Gauss-Legendre quadrature, not in the Fresnel
integrals the sectoral patterns (and the diagonal horn's efficiency) also have
closed forms in - which is what makes those closed forms a check on this module.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq
from scipy.special import j0, j1, jvp

__all__ = ["slit_uniform", "slit_cosine", "te11", "he11", "hpbw", "boresight_is_peak",
           "dual_mode", "disc_pattern", "disc_efficiency", "diagonal", "diagonal_efficiency", "peak_cross"]

_X, _W = np.polynomial.legendre.leggauss(400)
_XR, _WR = np.polynomial.legendre.leggauss(160)
_NPH = 96


def _slit(D: float, s: float, amp, sin_t: float) -> complex:
    x = 0.5 * D * _X                                   # aperture coordinate, wavelengths
    phase = -2 * math.pi * s * (2 * x / D) ** 2 + 2 * math.pi * x * sin_t
    return complex(np.sum(amp(x, D) * np.exp(1j * phase) * _W) * 0.5 * D)


def slit_uniform(D: float, s: float, theta: float) -> complex:
    """E-plane field of a uniform line aperture D wavelengths wide."""
    return (1 + math.cos(theta)) / 2 * _slit(D, s, lambda x, D: np.ones_like(x), math.sin(theta))


def slit_cosine(D: float, s: float, theta: float) -> complex:
    """H-plane field of a cosine-tapered line aperture D wavelengths wide."""
    return (1 + math.cos(theta)) / 2 * _slit(D, s, lambda x, D: np.cos(math.pi * x / D), math.sin(theta))


def _disc(D: float, s: float, ex, kx: float, ky: float) -> complex:
    a = D / 2
    r = 0.5 * a * (_XR + 1)
    wr = 0.5 * a * _WR * r
    ph = 2 * math.pi * np.arange(_NPH) / _NPH
    R, P = np.meshgrid(r, ph, indexing="ij")
    field = ex(R / a, P) * np.exp(-2j * math.pi * s * (R / a) ** 2)
    kern = np.exp(2j * math.pi * (kx * R * np.cos(P) + ky * R * np.sin(P)))
    return complex(np.sum(field * kern * wr[:, None]) * (2 * math.pi / _NPH))


def _te11_x(u, phi):
    v = 1.8411837813 * u
    with np.errstate(invalid="ignore", divide="ignore"):
        jv = np.where(v > 1e-12, j1(v) / v, 0.5)
    return jv * np.cos(phi) ** 2 + jvp(1, v) * np.sin(phi) ** 2


def te11(D: float, s: float, theta: float, plane: str) -> complex:
    """Co-polar field of a TE11 aperture D wavelengths across; 'E' is the
    plane containing the aperture E field at the centre, 'H' the one across it."""
    st = math.sin(theta)
    kx, ky = (st, 0.0) if plane == "E" else (0.0, st)
    return (1 + math.cos(theta)) / 2 * _disc(D, s, _te11_x, kx, ky)


def he11(D: float, s: float, theta: float) -> complex:
    st = math.sin(theta)
    return (1 + math.cos(theta)) / 2 * _disc(D, s, lambda u, phi: j0(2.404825557695773 * u), st, 0.0)


def boresight_is_peak(pattern, D: float, n: int = 200) -> bool:
    """True if nothing within the first null region outshines broadside."""
    p0 = abs(pattern(0.0)) ** 2
    th = np.linspace(0.0, min(math.pi / 2, math.asin(min(1.0, 2.0 / D))), n)
    return all(abs(pattern(t)) ** 2 <= p0 * (1 + 1e-9) for t in th)


def hpbw(pattern, D: float) -> float:
    """Full half-power beamwidth in degrees; the pattern is a function of theta."""
    p0 = abs(pattern(0.0)) ** 2
    f = lambda t: abs(pattern(t)) ** 2 / p0 - 0.5
    hi = min(math.pi / 2, math.asin(min(1.0, 4.0 / D)))
    ts = np.linspace(1e-6, hi, 240)
    v = [f(t) for t in ts]
    i = next((k for k in range(1, len(ts)) if v[k] < 0), None)
    if i is None:
        return float("nan")            # no half-power point: the beam has broken up
    return 2.0 * math.degrees(brentq(f, ts[i - 1], ts[i], xtol=1e-10))


# ------------------------------------------------ vector apertures (co / cross)
#
# With the Huygens obliquity factor the Ludwig-3 co- and cross-polar far fields
# of an aperture are just the Fourier transforms of its two field components,
# times (1 + cos theta)/2. So a two-component aperture needs nothing more than
# the scalar one: transform E_x and E_y separately.

def _disc_xy(D, s, fx, fy, kx, ky):
    a = D / 2
    r = 0.5 * a * (_XR + 1)
    wr = 0.5 * a * _WR * r
    ph = 2 * math.pi * np.arange(_NPH) / _NPH
    R, P = np.meshgrid(r, ph, indexing="ij")
    q = np.exp(-2j * math.pi * s * (R / a) ** 2) * wr[:, None] * (2 * math.pi / _NPH)
    kern = np.exp(2j * math.pi * (kx * R * np.cos(P) + ky * R * np.sin(P)))
    return (complex(np.sum(fx(R / a, P) * q * kern)), complex(np.sum(fy(R / a, P) * q * kern)))


_CHI_TE, _CHI_TM = 1.8411837813, 3.8317059702


def _te11_fields(u, phi):
    v = _CHI_TE * u
    with np.errstate(invalid="ignore", divide="ignore"):
        jv = np.where(v > 1e-12, j1(v) / v, 0.5)
    d = jvp(1, v)
    return jv * np.cos(phi) ** 2 + d * np.sin(phi) ** 2, (jv - d) * np.sin(phi) * np.cos(phi)


def _tm11_fields(u, phi):
    w = _CHI_TM * u
    with np.errstate(invalid="ignore", divide="ignore"):
        jw = np.where(w > 1e-12, j1(w) / w, 0.5)
    d = jvp(1, w)
    return d * np.cos(phi) ** 2 + jw * np.sin(phi) ** 2, (d - jw) * np.sin(phi) * np.cos(phi)


def _mode_power(fields):
    """Integral of |E_t|^2 over the unit-radius aperture (wave impedance ~ eta0
    for both modes on an aperture many wavelengths across)."""
    r = 0.5 * (_XR + 1)
    ph = 2 * math.pi * np.arange(_NPH) / _NPH
    R, P = np.meshgrid(r, ph, indexing="ij")
    ex, ey = fields(R, P)
    return float(np.sum((np.abs(ex) ** 2 + np.abs(ey) ** 2) * (0.5 * _WR * r)[:, None]) * 2 * math.pi / _NPH)


def dual_mode(p_tm: float, psi_deg: float = 0.0):
    """(fx, fy) of TE11 plus TM11, with a fraction p_tm of the aperture power in
    TM11, in the sense that tapers the E plane (Potter), TM11 leading TE11 by
    psi_deg - zero when the horn is on its design frequency."""
    if not 0.0 <= p_tm <= 1.0:
        raise ValueError(f"p_tm is a power fraction in [0, 1], got {p_tm}")
    pte, ptm = _mode_power(_te11_fields), _mode_power(_tm11_fields)
    # TE11 at unit amplitude and TM11 scaled to carry p_tm of the power; at
    # p_tm = 1 there is no TE11 left, so it is TM11 alone at TE11's power
    # (this used to return pure TE11 there)
    a = 1.0 if p_tm < 1 else 0.0
    b = math.sqrt(p_tm / (1 - p_tm) * pte / ptm) if p_tm < 1 else math.sqrt(pte / ptm)
    b = b * complex(math.cos(math.radians(psi_deg)), math.sin(math.radians(psi_deg)))
    # At the E-plane rim TE11's x-field is J1(1.841)/1.841 = +0.316 and TM11's is
    # J1'(3.832) = -0.403, so adding TM11 with the SAME sign at the centre is
    # what tapers the E plane towards the H plane's cosine-like edge
    fx = lambda u, phi: a * _te11_fields(u, phi)[0] + b * _tm11_fields(u, phi)[0]
    fy = lambda u, phi: a * _te11_fields(u, phi)[1] + b * _tm11_fields(u, phi)[1]
    return fx, fy


def disc_pattern(fields, D: float, s: float, theta: float, phi: float) -> tuple[complex, complex]:
    """(co, cross) Ludwig-3 far field, x-referenced, of a circular aperture."""
    st = math.sin(theta)
    cx, cy = _disc_xy(D, s, fields[0], fields[1], st * math.cos(phi), st * math.sin(phi))
    ob = (1 + math.cos(theta)) / 2
    return ob * cx, ob * cy


def disc_efficiency(fields, s: float) -> float:
    """Co-polar aperture efficiency: |int E_x|^2 / (A int |E|^2), on a unit-radius aperture."""
    cx, _ = _disc_xy(2.0, s, fields[0], fields[1], 0.0, 0.0)
    power = _mode_power(lambda u, phi: (fields[0](u, phi), fields[1](u, phi)))
    return abs(cx) ** 2 / (math.pi * power)


_XS, _WS = np.polynomial.legendre.leggauss(160)


def _square_xy(D, s, kx, ky):
    """Diagonal horn: E = y cos(pi x / a) + x cos(pi y / a), quadratic phase
    exp(-j 2 pi s ((2x/a)^2 + (2y/a)^2)) - s is the error at an edge's centre,
    2s at a corner. Returns the transforms of (E_x, E_y)."""
    x = 0.5 * D * _XS
    w = 0.5 * D * _WS
    X, Y = np.meshgrid(x, x, indexing="ij")
    Wt = np.outer(w, w)
    q = np.exp(-2j * math.pi * s * ((2 * X / D) ** 2 + (2 * Y / D) ** 2)) * Wt
    kern = np.exp(2j * math.pi * (kx * X + ky * Y))
    ex = np.cos(math.pi * Y / D)
    ey = np.cos(math.pi * X / D)
    return complex(np.sum(ex * q * kern)), complex(np.sum(ey * q * kern))


def diagonal(D: float, s: float, theta: float, phi: float) -> tuple[complex, complex]:
    """(co, cross) far field of the diagonal horn, Ludwig-3 referenced to the
    aperture diagonal (x + y)/sqrt(2) along which it is polarised."""
    st = math.sin(theta)
    cx, cy = _square_xy(D, s, st * math.cos(phi), st * math.sin(phi))
    ob = (1 + math.cos(theta)) / 2
    return ob * (cx + cy) / math.sqrt(2), ob * (cx - cy) / math.sqrt(2)


def diagonal_efficiency(s: float) -> float:
    D = 2.0
    cx, cy = _square_xy(D, s, 0.0, 0.0)
    # int |E|^2 = int cos^2(pi y/D) + cos^2(pi x/D) over the square = D^2
    return abs((cx + cy) / math.sqrt(2)) ** 2 / (D * D * D * D)


def peak_cross(pattern, D: float, phi: float, n: int = 200) -> float:
    """Highest cross-polar level in the plane phi, dB below the co-polar peak
    in that plane - boresight for a well-behaved horn, but a large phase error
    moves the co-polar maximum off axis, and normalising to boresight then
    overstated the cross-polar level (1.7 dB for a diagonal horn at s = 0.8)."""
    th = np.linspace(0.0, math.pi / 2 * 0.999, n)
    vals = [pattern(t, phi) for t in th]
    co = max(abs(v[0]) for v in vals)
    x = max(abs(v[1]) for v in vals)
    return 20 * math.log10(x / co) if x > 0 else float("-inf")
