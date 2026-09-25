"""Fresnel zone plate antennas by Kirchhoff diffraction - an arbiter for the spec.

A feed at the focus with power pattern 2(n+1)cos^n(theta) illuminates a flat
plate a focal length F away. The plate multiplies the incident spherical wave
by its transmission t: 0 or 1 by zone for an opaque (Soret) plate, +-1 for a
phase-reversal plate, and a staircase of phase steps for a multilevel one. The
field on the far side radiates by the scalar Kirchhoff integral with obliquity
(cos(incidence) + cos(observation))/2, and the gain is referred to the feed's
TOTAL power, so spillover past the rim is simply power that never reaches the
integral.

Everything is done on a Cartesian grid over the plate: zone membership comes
from each point's own path length to the focus, and the pattern from the full
2-D radiation integral. Nothing here uses the feed-angle substitution or the
first-order grating result the spec relies on. Lengths in wavelengths.
"""
from __future__ import annotations

import math

import numpy as np

__all__ = ["zone_radius", "transmission", "plate_field", "gain", "gain_bandwidth", "pattern", "hpbw"]


def zone_radius(m: float, F: float) -> float:
    """Radius at which the path to the focus exceeds F by m/2 wavelengths."""
    return math.sqrt(m * F + (m / 2.0) ** 2)


def transmission(excess, levels: int):
    """Plate transmission for a path excess (wavelengths) over the focal length.

    levels = 1: opaque rings on the odd half-wave zones, the central zone open.
    levels >= 2: each 1/levels-wavelength band of path is advanced by the same
    step, so the phase left behind is a sawtooth of height 2*pi/levels;
    levels = 2 is the phase-reversal plate.
    """
    if levels < 2:
        return (np.floor(2.0 * excess) % 2 < 0.5).astype(float)
    return np.exp(2j * math.pi / levels * np.floor(levels * excess))


def plate_field(F: float, M: int, levels: int, n_feed: float, cells: int = 1600,
                freq: float = 1.0):
    """(x, y, dA, field after the plate, cos(incidence)) on the plate's quarter
    x, y >= 0 - the pattern cuts used here are symmetric in both - with the
    feed's total power normalised to 1.

    freq is the operating frequency over the design frequency: the zones stay
    where the design wavelength put them, the phase is that of the new one,
    and lengths stay in DESIGN wavelengths (the gain functions rescale)."""
    r_out = zone_radius(M, F)
    h = r_out / cells
    c = (np.arange(cells) + 0.5) * h
    x, y = np.meshgrid(c, c, indexing="ij")
    rho2 = x ** 2 + y ** 2
    inside = rho2 < r_out ** 2
    r = np.sqrt(rho2 + F * F)
    cos_i = F / r
    # feed amplitude sqrt(U) with U = 2(n+1)cos^n / (4 pi), spherical spreading 1/r
    amp = np.sqrt(2.0 * (n_feed + 1.0) * cos_i ** n_feed / (4.0 * math.pi)) / r
    field = amp * np.exp(-2j * math.pi * freq * r) * transmission(r - F, levels) * inside
    return x, y, h * h, field, cos_i


def gain(F: float, M: int, levels: int, n_feed: float, cells: int = 1600,
         freq: float = 1.0) -> float:
    """Boresight gain over the feed's total power, linear. At freq != 1 the
    Kirchhoff integral's 1/wavelength scaling makes it freq^2 times the sum."""
    x, y, dA, field, cos_i = plate_field(F, M, levels, n_feed, cells, freq)
    s = 4.0 * np.sum(field * 0.5 * (cos_i + 1.0)) * dA          # four quadrants
    return float(4.0 * math.pi * freq ** 2 * abs(s) ** 2)


def gain_bandwidth(F: float, M: int, levels: int, n_feed: float, drop_db: float = 1.0,
                   cells: int = 800) -> float:
    """Fractional bandwidth over which the gain stays within drop_db of its
    peak, the plate fixed and the feed pattern held constant."""
    from scipy.optimize import brentq, minimize_scalar
    g = lambda f: gain(F, M, levels, n_feed, cells, f)
    w = min(0.85, 3.0 / M)
    peak = minimize_scalar(lambda f: -g(f), bounds=(1 - 0.3 / M, 1 + 0.3 / M), method="bounded").x
    thr = g(peak) * 10 ** (-drop_db / 10)
    lo = brentq(lambda f: g(f) - thr, 1 - w, peak)
    hi = brentq(lambda f: g(f) - thr, peak, 1 + w)
    return hi - lo


def pattern(F: float, M: int, levels: int, n_feed: float, angles_deg, cells: int = 1200):
    """Gain against angle in the x-z plane, linear, from the 2-D integral."""
    x, y, dA, field, cos_i = plate_field(F, M, levels, n_feed, cells)
    out = []
    for a in np.radians(np.asarray(angles_deg, dtype=float)):
        k = 2.0 * math.pi * math.sin(a)
        ob = 0.5 * (cos_i + math.cos(a))
        # the quarter plate mirrored in x: exp(jkx) + exp(-jkx) = 2cos(kx); in y, twice
        s = 4.0 * np.sum(field * ob * np.cos(k * x)) * dA
        out.append(4.0 * math.pi * abs(s) ** 2)
    return np.array(out)


def hpbw(F: float, M: int, levels: int, n_feed: float, cells: int = 1200) -> float:
    """Full half-power beamwidth in degrees (scalar model, so one plane serves)."""
    D = 2.0 * zone_radius(M, F)
    guess = 58.0 / D
    a = np.linspace(0.0, 2.0 * guess, 161)
    g = pattern(F, M, levels, n_feed, a, cells)
    g = g / g[0]
    i = int(np.argmax(g < 0.5))
    return 2.0 * float(np.interp(0.5, [g[i], g[i - 1]], [a[i], a[i - 1]]))
