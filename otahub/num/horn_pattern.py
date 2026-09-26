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

The integrals are done by Gauss-Legendre quadrature, not in the Fresnel
integrals the sectoral patterns also have closed forms in - which is what makes
those closed forms a check on this module.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq
from scipy.special import j0, j1, jvp

__all__ = ["slit_uniform", "slit_cosine", "te11", "he11", "hpbw", "boresight_is_peak"]

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
