"""Single-surface collimating lenses by ray tracing - an arbiter for the lens specs.

Feed at the focus, lens axis +z, a curved first face and a flat back. For index
n the equal-path face is r(t) = |n - 1| F / |n cos t - 1|: a hyperbola for a
dielectric (n > 1), an ellipse for a metal-plate lens (n < 1). Each ray is
refracted by Snell's law at that face, with the normal taken from the curve
itself, and followed to the flat face; nothing assumes the closed-form mapping
the specs use. Power conserved along ray tubes then gives the aperture field
for a feed of power pattern 2(n+1)cos^n, and from it the taper efficiency and
the beam. Lengths in units of F.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq
from scipy.special import j0

__all__ = ["rim_limit", "trace", "aperture", "efficiencies", "hpbw"]


def rim_limit(n: float) -> float:
    """Largest usable rim angle [rad]: the hyperbola runs to infinity at cos t = 1/n;
    the ellipse's aperture radius peaks, and folds back, at cos t = n."""
    return math.acos(1.0 / n) if n > 1 else math.acos(n)


def _surface(n: float, t: float) -> float:
    return (n - 1.0) / (n * math.cos(t) - 1.0)


def trace(n: float, t: float, h: float = 1e-7) -> tuple[np.ndarray, float]:
    """Refract the feed ray at angle t: (direction inside the lens, aperture height)."""
    r = _surface(n, t)
    p = np.array([r * math.sin(t), r * math.cos(t)])
    r2, r1 = _surface(n, t + h), _surface(n, t - h)
    tang = np.array([r2 * math.sin(t + h) - r1 * math.sin(t - h), r2 * math.cos(t + h) - r1 * math.cos(t - h)])
    tang /= np.linalg.norm(tang)
    normal = np.array([-tang[1], tang[0]])
    d = p / np.linalg.norm(p)
    if normal @ d < 0:
        normal = -normal
    s_t = (d - (normal @ d) * normal) / n
    out = s_t + math.sqrt(max(1.0 - s_t @ s_t, 0.0)) * normal
    return out / np.linalg.norm(out), float(p[0])


def aperture(n: float, trim: float, nf: float, m: int = 3000):
    """(theta, rho, drho/dtheta, aperture field) from traced rays, rho in units of F."""
    t = np.linspace(1e-6, trim, m)
    rho = np.array([trace(n, x)[1] for x in t])
    drho = np.gradient(rho, t)
    return t, rho, drho, np.sqrt(np.cos(t) ** nf * np.sin(t) / (rho * drho))


def efficiencies(n: float, trim: float, nf: float) -> tuple[float, float]:
    """(taper, spillover) for a 2(n+1)cos^n feed and rim angle trim."""
    t, rho, drho, A = aperture(n, trim, nf)
    dS = 2.0 * math.pi * rho * drho
    R = rho[-1]
    taper = np.trapezoid(A * dS, t) ** 2 / (math.pi * R * R * np.trapezoid(A * A * dS, t))
    spill = np.trapezoid(2 * (nf + 1) * np.cos(t) ** nf * np.sin(t), t) / 2.0
    return float(taper), float(spill)


def hpbw(n: float, trim: float, nf: float) -> float:
    """Half-power beamwidth in units of lambda/D degrees, from the Hankel transform."""
    t, rho, drho, A = aperture(n, trim, nf)
    R = rho[-1]
    rn = rho / R
    aw = A * rho * drho
    wt = np.gradient(t)
    aw = aw * wt
    f0 = aw.sum()
    F = lambda u: float(j0(u * rn) @ aw / f0)
    us = np.arange(0.0, 12.0, 0.002)
    v = j0(np.outer(us, rn)) @ aw / f0
    i3 = int(np.argmax(v ** 2 < 0.5))
    return 2.0 * math.degrees(brentq(lambda u: F(u) ** 2 - 0.5, us[i3 - 1], us[i3])) / math.pi
