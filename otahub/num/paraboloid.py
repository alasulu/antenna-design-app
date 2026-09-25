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

__all__ = ["half_angle", "feed_exponent", "efficiencies", "far_field", "beam",
           "offset_geometry", "offset_rim_angles", "offset_taper", "offset_beamwidths",
           "cylinder_feed_exponent", "cylinder_efficiencies", "cylinder_beam",
           "dual_trace", "dual_efficiencies", "blocked_beam"]


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


# --------------------------------------------------------------- offset dish

def offset_geometry(f_over_d_parent: float, h0_over_d: float) -> tuple[float, float]:
    """(feed tilt from the parent axis, rim half-angle theta*) [rad]. The rim of an
    offset paraboloid is a circular cone as seen from the focus, about the
    bisector of the upper and lower rim angles - checked in the tests, not assumed."""
    up = 2.0 * math.atan((h0_over_d + 0.5) / (2.0 * f_over_d_parent))
    lo = 2.0 * math.atan((h0_over_d - 0.5) / (2.0 * f_over_d_parent))
    return 0.5 * (up + lo), 0.5 * (up - lo)


def _offset_aperture(f_over_d_parent, h0_over_d, n, m=80):
    """Aperture field A = cos^(n/2)(psi)/r on the projected circle (D = 1), where psi is
    the angle from the tilted feed axis and r the focus-to-dish distance."""
    tilt, _ = offset_geometry(f_over_d_parent, h0_over_d)
    F = f_over_d_parent
    xr, wr = np.polynomial.legendre.leggauss(m)
    xp, wp = np.polynomial.legendre.leggauss(2 * m)
    rr, wrr = 0.25 * (xr + 1.0), 0.25 * wr
    pp, wpp = math.pi * (xp + 1.0), math.pi * wp
    R, P = np.meshgrid(rr, pp, indexing="ij")
    W = np.outer(wrr, wpp) * R
    xl, yl = R * np.cos(P), R * np.sin(P)
    x = h0_over_d + xl
    rho2 = x * x + yl * yl
    zc = F - rho2 / (4.0 * F)
    r = np.sqrt(rho2 + zc * zc)
    cospsi = (x * math.sin(tilt) + zc * math.cos(tilt)) / r
    return xl, yl, W, np.clip(cospsi, 0.0, None) ** (n / 2.0) / r


def offset_rim_angles(f_over_d_parent: float, h0_over_d: float, points: int = 73) -> np.ndarray:
    """Angle of each rim ray from the tilted feed axis [rad], for the circularity check."""
    tilt, _ = offset_geometry(f_over_d_parent, h0_over_d)
    F = f_over_d_parent
    axis = np.array([math.sin(tilt), 0.0, math.cos(tilt)])
    out = []
    for ph in np.linspace(0.0, 2.0 * math.pi, points):
        x, y = h0_over_d + 0.5 * math.cos(ph), 0.5 * math.sin(ph)
        v = np.array([x, y, F - (x * x + y * y) / (4.0 * F)])
        out.append(math.acos(float(np.clip(v @ axis / np.linalg.norm(v), -1.0, 1.0))))
    return np.array(out)


def offset_taper(f_over_d_parent: float, h0_over_d: float, n: float) -> float:
    """Taper efficiency of the offset dish, by 2-D integration over its aperture."""
    _, _, W, A = _offset_aperture(f_over_d_parent, h0_over_d, n)
    return float((A * W).sum() ** 2 / ((math.pi / 4.0) * (A * A * W).sum()))


def offset_beamwidths(f_over_d_parent: float, h0_over_d: float, n: float) -> tuple[float, float]:
    """HPBW in the plane of the offset and across it, in units of lambda/D degrees."""
    xl, yl, W, A = _offset_aperture(f_over_d_parent, h0_over_d, n)
    aw, total = (A * W).ravel(), float((A * W).sum())
    widths = []
    for c in (xl.ravel(), yl.ravel()):
        p = lambda u: abs((aw * np.exp(1j * u * c)).sum() / total) ** 2 - 0.5
        widths.append(math.degrees(brentq(p, 0.1, 20.0) - brentq(p, -20.0, -0.1)) / (2.0 * math.pi))
    return widths[0], widths[1]


# ------------------------------------------------------- parabolic cylinder

def cylinder_feed_exponent(edge_taper_db: float, f_over_w: float) -> float:
    """n for a line feed whose transverse power pattern is cos^n: the rim sits
    edge_taper_db below the centre, cylindrical spreading sqrt((1+cos)/2) included."""
    t0 = half_angle(f_over_w)
    if t0 >= math.pi / 2:
        raise ValueError("f/W must exceed 0.25 for a finite edge taper")
    return 2.0 * (edge_taper_db / 20.0 - 0.5 * math.log10((1 + math.cos(t0)) / 2)) / math.log10(math.cos(t0))


def _cylinder_aperture(n: float, f_over_w: float, m: int = 1500):
    """A(y) = sqrt(cos^n(theta)/r) over y in [-1/2, 1/2] (W = 1): a cylindrical wave."""
    x, w = np.polynomial.legendre.leggauss(m)
    y, wy = 0.5 * x, 0.5 * w
    t = 2.0 * np.arctan(y / (2.0 * f_over_w))
    r = 2.0 * f_over_w / (1.0 + np.cos(t))
    return y, wy, np.sqrt(np.cos(t) ** n / r)


def cylinder_efficiencies(n: float, f_over_w: float) -> tuple[float, float]:
    """(taper, spillover): taper from the aperture field along y, spillover by
    quadrature of the feed's transverse pattern."""
    y, wy, A = _cylinder_aperture(n, f_over_w)
    taper = float((A * wy).sum() ** 2 / (A * A * wy).sum())
    t0 = half_angle(f_over_w)
    x, w = np.polynomial.legendre.leggauss(400)
    inside = (np.cos(0.5 * t0 * (x + 1)) ** n * w).sum() * 0.5 * t0
    whole = (np.cos(0.25 * math.pi * (x + 1)) ** n * w).sum() * 0.25 * math.pi
    return taper, float(inside / whole)


def cylinder_beam(n: float, f_over_w: float) -> tuple[float, float]:
    """(HPBW in units of lambda/W degrees, PEAK sidelobe in dB) in the focusing plane.
    The peak, not the first: under heavy taper the first sidelobe collapses as its
    two nulls merge, and a later lobe carries the peak."""
    y, wy, A = _cylinder_aperture(n, f_over_w)
    aw = A * wy
    f0 = aw.sum()
    us = np.arange(0.0, 60.0, 0.005)
    v = np.cos(np.outer(us, y)) @ aw / f0
    F = lambda u: float(np.cos(u * y) @ aw / f0)
    i3 = int(np.argmax(v ** 2 < 0.5))
    u3 = brentq(lambda u: F(u) ** 2 - 0.5, us[i3 - 1], us[i3])
    flips = np.nonzero(np.sign(v[1:]) != np.sign(v[:-1]))[0]
    n1 = brentq(F, us[flips[0]], us[flips[0] + 1])
    k = int(np.argmax(np.abs(v) * (us > n1)))
    r = minimize_scalar(lambda u: -abs(F(u)), bounds=(us[k] - 0.01, us[k] + 0.01),
                        method="bounded", options={"xatol": 1e-10})
    return 2.0 * math.degrees(u3) / (2.0 * math.pi), 20.0 * math.log10(abs(F(r.x)))


# ----------------------------------------------------------- dual reflectors

def dual_trace(kind: str, magnification: float, theta: float) -> float:
    """Trace one ray from the feed off the subreflector [rad in, rad out].

    Prime focus at the origin, main dish toward -z, feed at z = -2c (c = 1).
    Cassegrain: hyperboloid branch nearer the prime focus, |P-F2| - |P-F1| = 2a,
    the ray leaving as if from the prime focus. Gregorian: ellipsoid beyond it,
    |P-F1| + |P-F2| = 2a, the ray passing through it. Returns psi, the angle
    from the prime focus toward the dish, measured from the -z axis.
    """
    M = magnification
    e = (M + 1) / (M - 1) if kind == "cassegrain" else (M - 1) / (M + 1)
    c = 1.0
    a = c / e
    f2 = np.array([0.0, 0.0, -2.0 * c])
    d = np.array([math.sin(theta), 0.0, math.cos(theta)])
    if kind == "cassegrain":
        t = brentq(lambda s: s - np.linalg.norm(f2 + s * d) - 2 * a, 1e-12, 100.0 * c)
        v = f2 + t * d
    else:
        t = brentq(lambda s: s + np.linalg.norm(f2 + s * d) - 2 * a, 1e-12, 2 * a)
        v = -(f2 + t * d)
    v = v / np.linalg.norm(v)
    return math.acos(-v[2])


def dual_efficiencies(kind: str, magnification: float, f_over_d: float, ds_over_d: float,
                      n: float, m: int = 4000) -> tuple[float, float, float]:
    """(taper, spillover past the subreflector, blockage) from TRACED rays: each
    feed angle is mapped to an aperture radius through the real subreflector and
    main dish, and power is conserved along ray tubes. No equivalent paraboloid."""
    tf = 2.0 * math.atan(1.0 / (4.0 * magnification * f_over_d))    # the feed's rim angle
    th = np.linspace(1e-7, tf, m)
    psi = np.array([dual_trace(kind, magnification, t) for t in th])
    rho = 2.0 * f_over_d * np.tan(psi / 2.0)                         # D = 1
    drho = np.gradient(rho, th)
    G = 2.0 * (n + 1.0) * np.cos(th) ** n
    A = np.sqrt(G * np.sin(th) / (rho * drho))
    dS = 2.0 * math.pi * rho * drho
    full = np.trapezoid(A * dS, th)
    taper = full ** 2 / ((math.pi / 4.0) * np.trapezoid(A * A * dS, th))
    spill = np.trapezoid(G * np.sin(th), th) / 2.0
    unblocked = np.trapezoid(np.where(rho >= ds_over_d / 2.0, A, 0.0) * dS, th)
    return float(taper), float(spill), float((unblocked / full) ** 2)


def blocked_beam(edge_taper_db: float, f_over_d: float, ds_over_d: float,
                 m: int = 1500) -> tuple[float, float]:
    """(HPBW in lambda/D degrees, PEAK sidelobe dB) of a paraboloid aperture with a
    central blocked disc - the equivalent paraboloid of a dual reflector. The peak,
    because blockage raises the sidelobes and can move the highest off the first."""
    t0 = half_angle(f_over_d)
    n = 2.0 * (edge_taper_db / 20.0 - math.log10((1 + math.cos(t0)) / 2)) / math.log10(math.cos(t0))
    rb = ds_over_d / 2.0
    x, w = np.polynomial.legendre.leggauss(m)
    rho = rb + (0.5 - rb) * 0.5 * (x + 1.0)
    wr = (0.5 - rb) * 0.5 * w
    t = 2.0 * np.arctan(rho / (2.0 * f_over_d))
    aw = np.cos(t) ** (n / 2.0) * (1.0 + np.cos(t)) / 2.0 * rho * wr
    f0 = aw.sum()
    us = np.arange(0.0, 40.0, 0.005)
    v = j0(np.outer(us, rho)) @ aw / f0
    F = lambda u: float(j0(u * rho) @ aw / f0)
    i3 = int(np.argmax(v ** 2 < 0.5))
    u3 = brentq(lambda u: F(u) ** 2 - 0.5, us[i3 - 1], us[i3])
    flips = np.nonzero(np.sign(v[1:]) != np.sign(v[:-1]))[0]
    n1 = brentq(F, us[flips[0]], us[flips[0] + 1])
    k = int(np.argmax(np.abs(v) * (us > n1)))
    r = minimize_scalar(lambda u: -abs(F(u)), bounds=(us[k] - 0.01, us[k] + 0.01),
                        method="bounded", options={"xatol": 1e-10})
    return 2.0 * math.degrees(u3) / (2.0 * math.pi), 20.0 * math.log10(abs(F(r.x)))
