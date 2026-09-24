"""Prime-focus paraboloid under a cos^n feed, by direct aperture integration.

An arbiter for `prime_focus_parabolic`. The feed's power pattern is
2(n+1) cos^n(theta); a ray leaving it at theta lands on the aperture at
rho = 2 f tan(theta/2), and the 1/r spreading of the spherical wave over the
path 2f/(1 + cos theta) makes the aperture field

    A(rho) = cos^(n/2)(theta) * (1 + cos theta)/2.

Everything here works from that aperture field directly - the taper efficiency
as |int A dS|^2 / (S int |A|^2 dS), the pattern as its Hankel transform - so it
shares no algebra with the feed-angle integral the spec evaluates. Lengths are
in units of the dish radius.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import j0

__all__ = ["half_angle", "feed_exponent", "efficiencies", "far_field", "beam"]


def half_angle(f_over_d: float) -> float:
    """Angle the dish rim subtends at the focus [rad]."""
    return 2.0 * math.atan(1.0 / (4.0 * f_over_d))


def feed_exponent(edge_taper_db: float, f_over_d: float) -> float:
    """n such that the rim sits edge_taper_db below the centre, spreading included."""
    t0 = half_angle(f_over_d)
    if t0 >= math.pi / 2:
        raise ValueError("f/D must exceed 0.25 for a finite edge taper")
    return 2.0 * (edge_taper_db / 20.0 - math.log10((1 + math.cos(t0)) / 2)) / math.log10(math.cos(t0))


def _aperture(n: float, f_over_d: float, m: int = 1500):
    x, w = np.polynomial.legendre.leggauss(m)
    rho, wr = 0.5 * (x + 1.0), 0.5 * w
    t = 2.0 * np.arctan(rho / (4.0 * f_over_d))          # f = 2 f_over_d radii
    return rho, wr, np.cos(t) ** (n / 2.0) * (1.0 + np.cos(t)) / 2.0


def efficiencies(n: float, f_over_d: float) -> tuple[float, float]:
    """(taper, spillover). Taper from the aperture field; spillover by quadrature
    of the feed pattern over the rim angle, not by its closed form."""
    rho, wr, A = _aperture(n, f_over_d)
    taper = (2 * math.pi * (A * rho * wr).sum()) ** 2 / (math.pi * 2 * math.pi * (A * A * rho * wr).sum())
    t0 = half_angle(f_over_d)
    x, w = np.polynomial.legendre.leggauss(400)
    th = 0.5 * t0 * (x + 1)
    inside = ((2 * (n + 1) * np.cos(th) ** n * np.sin(th)) * w).sum() * 0.5 * t0
    return float(taper), float(inside) / 2.0               # the whole pattern integrates to 2


def far_field(n: float, f_over_d: float):
    """Normalised pattern F(u), u = k a sin(theta) with a the dish radius."""
    rho, wr, A = _aperture(n, f_over_d)
    f0 = (A * rho * wr).sum()
    return lambda u: float((A * j0(u * rho) * rho * wr).sum() / f0)


def beam(n: float, f_over_d: float) -> tuple[float, float]:
    """(HPBW in units of lambda/D degrees, first sidelobe in dB) - the first sidelobe
    being the lobe between the first and second nulls."""
    F = far_field(n, f_over_d)
    u3 = brentq(lambda u: F(u) ** 2 - 0.5, 0.1, 5.0)
    us = np.arange(u3, 30.0, 0.02)
    vals = np.array([F(u) for u in us])
    flips = np.nonzero(np.sign(vals[1:]) != np.sign(vals[:-1]))[0]
    n1 = brentq(F, us[flips[0]], us[flips[0] + 1])
    n2 = brentq(F, us[flips[1]], us[flips[1] + 1])
    r = minimize_scalar(lambda u: -abs(F(u)), bounds=(n1, n2), method="bounded",
                        options={"xatol": 1e-9})
    return 2.0 * math.degrees(u3) / math.pi, 20.0 * math.log10(abs(F(r.x)))
