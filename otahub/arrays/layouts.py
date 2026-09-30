"""Circular apertures and thinned arrays - the layouts a rectangular grid cannot
describe.

Circular layouts: a single ring, concentric rings, and any lattice clipped to a
circle. A circular aperture wants a circularly symmetric taper, and the one
that plays the part Dolph and Taylor play for a line is Taylor's circular
distribution (T. T. Taylor, "Design of circular apertures for narrow beamwidth
and low sidelobes", IRE Trans. AP-8, 1960). Its pattern, in u = (2a/lambda)
sin(theta) for an aperture of radius a, is

    F(u) = [2 J1(pi u)/(pi u)] * prod_{n<nbar} (1 - u^2/u_n^2) / (1 - u^2/mu_n^2)

- the uniform disc's pattern with its first nbar - 1 nulls mu_n (zeros of
J1(pi mu)) moved to u_n = mu_nbar sqrt((A^2 + (n - 1/2)^2)/(A^2 + (nbar -
1/2)^2)), A = arccosh(R)/pi, R the main-beam to sidelobe voltage ratio - and
the aperture distribution that radiates it is a finite Fourier-Bessel series,

    g(p) = sum_{m<nbar} F(mu_m)/J0(pi mu_m)^2 J0(pi mu_m p),    p = r/a,

normalised here to 1 at the centre. `circular_taylor` samples it at the
elements.

Thinned arrays: switching elements off at random, with the probability of
keeping each one following a taper (a statistical density taper), gives a
sparse array whose average pattern is exact:

    E|AF(u)|^2 = |sum p_i exp(j k r_i.u)|^2 + sum p_i (1 - p_i)

- the density-tapered array's pattern on a flat floor. `thinned_expected_power`
evaluates it; a thinned array's sidelobes scatter about that floor, not about
the taper's design level.

Positions in wavelengths, theta from the array normal, as in `planar`.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.special import j0, j1, jn_zeros

__all__ = ["ring", "concentric_rings", "clip_to_circle", "taylor_circular_pattern",
           "taylor_circular_distribution", "equal_area_radius", "circular_taylor", "thin", "thinned_expected_power", "thinned_expected_directivity"]


def ring(n: int, radius: float, phase0_deg: float = 0.0) -> np.ndarray:
    """`n` elements equally spaced on a circle of `radius` wavelengths."""
    if n < 1 or radius < 0:
        raise ValueError("need n >= 1 and a non-negative radius")
    a = math.radians(phase0_deg) + 2 * math.pi * np.arange(n) / n
    return np.column_stack([radius * np.cos(a), radius * np.sin(a)])


def concentric_rings(radii, spacing: float, centre: bool = True) -> np.ndarray:
    """Rings at the given radii, each holding as many elements as keep the arc
    spacing at or just above `spacing` (wavelengths); a centre element if asked."""
    pts = [np.zeros((1, 2))] if centre else []
    for r in radii:
        if r <= 0:
            continue
        pts.append(ring(max(1, int(math.floor(2 * math.pi * r / spacing))), r))
    return np.vstack(pts)


def clip_to_circle(positions, radius: float) -> np.ndarray:
    """The elements of a lattice that lie within `radius` of its centre."""
    pos = np.asarray(positions, dtype=float)
    return pos[np.hypot(pos[:, 0], pos[:, 1]) <= radius * (1 + 1e-12)]


def _taylor_nulls(sidelobe_db: float, nbar: int):
    R = 10 ** (abs(sidelobe_db) / 20)
    A = math.acosh(R) / math.pi
    mu = jn_zeros(1, nbar) / math.pi                     # zeros of J1(pi mu)
    sigma = mu[nbar - 1] / math.sqrt(A * A + (nbar - 0.5) ** 2)
    un = sigma * np.sqrt(A * A + (np.arange(1, nbar) - 0.5) ** 2)
    return mu, un


def taylor_circular_pattern(u, sidelobe_db: float = -30.0, nbar: int = 5):
    """Taylor's circular pattern F(u), F(0) = 1.

    At a displaced null u = mu_m the factor 2 J1(pi u)/(pi u) and 1 - u^2/mu_m^2
    both vanish; the limit (J2 = -J0 at a zero of J1) is
    -J0(pi mu_m) prod_n (1 - mu_m^2/u_n^2) / prod_{n != m} (1 - mu_m^2/mu_n^2),
    taken within a relative 1e-8 of mu_m where the quotient loses its digits.
    """
    mu, un = _taylor_nulls(sidelobe_db, nbar)
    u = np.asarray(u, dtype=float)
    x = math.pi * u
    base = np.where(np.abs(x) < 1e-12, 1.0, 2 * j1(x) / np.where(np.abs(x) < 1e-12, 1.0, x))
    with np.errstate(divide="ignore", invalid="ignore"):
        for n in range(nbar - 1):
            base = base * (1 - (u / un[n]) ** 2) / (1 - (u / mu[n]) ** 2)
    for m in range(nbar - 1):
        at = np.abs(np.abs(u) - mu[m]) <= 1e-8 * mu[m]
        if np.any(at):
            limit = (-j0(math.pi * mu[m]) * np.prod(1 - (mu[m] / un) ** 2)
                     / np.prod([1 - (mu[m] / mu[n]) ** 2
                                for n in range(nbar - 1) if n != m]))
            base = np.where(at, limit, base)
    return base


def taylor_circular_distribution(p, sidelobe_db: float = -30.0, nbar: int = 5):
    """The aperture distribution g(p), p = r/a in [0, 1], normalised to g(0) = 1."""
    mu, un = _taylor_nulls(sidelobe_db, nbar)
    p = np.asarray(p, dtype=float)
    mus = np.concatenate([[0.0], mu[:nbar - 1]])
    g = np.zeros_like(p)
    for m in mus:
        if m == 0.0:
            Fm = 1.0
        else:
            num = np.prod(1 - (m / un) ** 2)
            den = np.prod([1 - (m / mu[n]) ** 2 for n in range(nbar - 1) if not math.isclose(mu[n], m)])
            Fm = -j0(math.pi * m) * num / den
        g = g + Fm / j0(math.pi * m) ** 2 * j0(math.pi * m * p)
    g0 = sum((1.0 if m == 0.0 else
              -j0(math.pi * m) * np.prod(1 - (m / un) ** 2)
              / np.prod([1 - (m / mu[n]) ** 2 for n in range(nbar - 1) if not math.isclose(mu[n], m)]))
             / j0(math.pi * m) ** 2 for m in mus)
    return g / g0


def equal_area_radius(n_elements: int, cell_area: float) -> float:
    """Radius of the disc with the elements' total cell area: sqrt(N A_cell / pi)."""
    return math.sqrt(n_elements * cell_area / math.pi)


def circular_taylor(positions, radius: float | None = None, sidelobe_db: float = -30.0, nbar: int = 5,
                    cell_area: float | None = None) -> np.ndarray:
    """Taylor circular weights sampled at each element for an aperture of
    `radius` wavelengths - by default the equal-area radius from `cell_area`
    (d^2 for a square lattice, (sqrt(3)/2) s^2 for a triangular one). That choice
    matters: on a lattice clipped to a circle it keeps a -40 dB design within
    about a decibel, where taking the outermost element's radius plus half a
    spacing loses 2-5 dB."""
    pos = np.asarray(positions, dtype=float)
    if radius is None:
        if cell_area is None:
            raise ValueError("give the aperture radius or the lattice cell area")
        radius = equal_area_radius(len(pos), cell_area)
    p = np.clip(np.hypot(pos[:, 0], pos[:, 1]) / radius, 0.0, 1.0)
    return taylor_circular_distribution(p, sidelobe_db, nbar)


def thin(density, seed: int | None = None) -> np.ndarray:
    """Keep each element with probability density_i (clipped to [0, 1]); a boolean mask."""
    p = np.clip(np.asarray(density, dtype=float), 0.0, 1.0)
    return np.random.default_rng(seed).random(p.shape) < p


def thinned_expected_power(positions, density, theta, phi) -> np.ndarray:
    """E|AF|^2 over thinnings with keep probabilities `density`, unit weights on
    the elements kept: the density-tapered pattern plus sum p (1 - p)."""
    pos = np.asarray(positions, dtype=float)
    p = np.clip(np.asarray(density, dtype=float), 0.0, 1.0)
    theta, phi = np.asarray(theta, dtype=float), np.asarray(phi, dtype=float)
    u, v = np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi)
    af = np.zeros(np.broadcast(u, v).shape, dtype=complex)
    for (x, y), pi in zip(pos, p):
        af = af + pi * np.exp(2j * math.pi * (x * u + y * v))
    return np.abs(af) ** 2 + float(np.sum(p * (1 - p)))


def thinned_expected_directivity(positions, density, scan_theta_deg: float = 0.0, scan_phi_deg: float = 0.0,
                                 element=None) -> float:
    """Directivity toward the scan direction from the expected numerator and the
    expected radiated power over thinnings (unit weights on the elements kept):

        [(sum p)^2 + sum p(1-p)] / [p K p^T + sum p(1-p) K(0)]   (times 4 pi p_el/...)

    with K the planar power kernel (sinc, or an element's). A ratio of
    expectations: the mean directivity of the realisations agrees with it to a
    fraction of a percent once a few hundred elements are kept."""
    from .elements import isotropic, power_kernel
    from .planar import steering_phase
    el = element or isotropic()
    pos = np.asarray(positions, dtype=float)
    p = np.clip(np.asarray(density, dtype=float), 0.0, 1.0)
    diff = np.round((pos[:, None, :] - pos[None, :, :]).reshape(-1, 2), 9)
    uniq, inv = np.unique(diff, axis=0, return_inverse=True)
    Kd = power_kernel(el, uniq)[inv.ravel()].reshape(len(pos), len(pos))
    q = p * steering_phase(pos, scan_theta_deg, scan_phi_deg)
    den = float(np.real(q @ Kd @ np.conj(q))) + float(np.sum(p * (1 - p))) * float(np.real(Kd[0, 0]))
    num = float(p.sum()) ** 2 + float(np.sum(p * (1 - p)))
    pk = float(el(math.radians(scan_theta_deg), math.radians(scan_phi_deg)))
    return 4 * math.pi * pk * num / den
