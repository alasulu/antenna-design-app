"""Element patterns for planar arrays, and the exact power integral with one.

The planar module's directivity assumes isotropic elements, where the radiated
power is a double sum over sinc(k |r_m - r_n|). With an element power pattern
p(theta, phi) the kernel becomes

    K(d) = integral of p(u) exp(j k d.u) dOmega

and for a separation d in the array plane (azimuth phi_d) the azimuthal
integral is exact once p is written as azimuthal harmonics c_m(theta):

    integral of exp(j m phi) exp(j x cos(phi - phi_d)) dphi = 2 pi j^m J_m(x) exp(j m phi_d),

x = k d sin(theta). What is left is one integral in u = cos(theta), taken by
Gauss-Legendre over each hemisphere separately, so a pattern cut off at the
horizon (a ground-backed element) has its kink at a node boundary. The
harmonics come from an FFT of p over azimuth, exact for patterns that are
trigonometric polynomials in phi of degree up to `mmax`.

What an element pattern here is: the power pattern of ONE element, theta from
the array normal. Mutual coupling enters only if the pattern supplied is an
embedded (active) element pattern, measured or simulated with its neighbours
present; the built-in ones are isolated elements.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.special import jv

__all__ = ["ElementPattern", "isotropic", "cosine", "short_dipole", "custom", "power_kernel",
           "element_directivity"]

K = 2.0 * math.pi


@dataclass(frozen=True)
class ElementPattern:
    """Power pattern p(theta, phi) >= 0 (theta from the array normal), whether it
    is zero below the array plane, and its highest azimuthal harmonic."""
    name: str
    power: Callable
    half_space: bool
    mmax: int = 0

    def __call__(self, theta, phi):
        theta = np.asarray(theta, dtype=float)
        p = np.asarray(self.power(theta, np.asarray(phi, dtype=float)), dtype=float)
        if self.half_space:
            p = np.where(np.cos(theta) >= 0, p, 0.0)
        return p


def isotropic() -> ElementPattern:
    return ElementPattern("isotropic", lambda t, f: np.ones(np.broadcast(t, f).shape), False, 0)


def cosine(q: float = 1.0) -> ElementPattern:
    """cos^q(theta) in power over the upper hemisphere, nothing below: a
    ground-backed element. q = 1 is the ideal element of a large matched array,
    whose gain follows the projected cell area."""
    if q < 0:
        raise ValueError("q must be non-negative")
    return ElementPattern(f"cos^{q:g}", lambda t, f: np.clip(np.cos(t), 0.0, None) ** q
                          * np.ones(np.broadcast(t, f).shape), True, 0)


def short_dipole(axis: str = "x", height: float | None = None) -> ElementPattern:
    """A short dipole in the array plane along x or y; with `height` (in
    wavelengths) it sits that far above an infinite ground plane, whose image
    multiplies the power by 4 sin^2(k h cos theta) over the upper hemisphere."""
    ax = axis.strip().lower()
    if ax not in ("x", "y"):
        raise ValueError("axis must be 'x' or 'y'")

    def p(t, f):
        s = np.sin(t) * (np.cos(f) if ax == "x" else np.sin(f))
        out = 1.0 - s * s
        if height is not None:
            out = out * 4.0 * np.sin(K * height * np.cos(t)) ** 2
        return out

    name = f"short dipole along {ax}" + ("" if height is None else f", {height:g} lambda over ground")
    return ElementPattern(name, p, height is not None, 2)


def custom(power: Callable, half_space: bool = False, mmax: int = 16, name: str = "custom") -> ElementPattern:
    """Any power pattern p(theta, phi), vectorised; `mmax` bounds its azimuthal
    harmonics (raise it for a pattern with fine azimuthal detail)."""
    return ElementPattern(name, power, half_space, mmax)


def _harmonics(el: ElementPattern, theta: np.ndarray):
    """c_m(theta) for m = -mmax..mmax, shape (len(theta), 2*mmax + 1)."""
    M = 4 * el.mmax + 8
    phi = 2 * math.pi * np.arange(M) / M
    c = np.fft.fft(el(theta[:, None], phi[None, :]), axis=1) / M
    m = np.arange(-el.mmax, el.mmax + 1)
    return c[:, m % M], m


def _nodes(el: ElementPattern, n: int):
    """(theta, weight in u = cos theta) for each hemisphere the element radiates into."""
    x, w = np.polynomial.legendre.leggauss(n)
    us, ws = [0.5 * (x + 1)], [0.5 * w]
    if not el.half_space:
        us.append(-0.5 * (x + 1))
        ws.append(0.5 * w)
    u = np.concatenate(us)
    return np.arccos(np.clip(u, -1.0, 1.0)), np.concatenate(ws)


def power_kernel(el: ElementPattern, sep, n: int | None = None) -> np.ndarray:
    """K(d) for in-plane separations `sep` (..., 2) in wavelengths."""
    sep = np.asarray(sep, dtype=float)
    d = np.hypot(sep[..., 0], sep[..., 1])
    phid = np.arctan2(sep[..., 1], sep[..., 0])
    if n is None:
        n = int(K * float(d.max(initial=0.0)) / 2) + 48
    theta, wu = _nodes(el, n)
    c, m = _harmonics(el, theta)
    x = K * d.ravel()[:, None] * np.sin(theta)[None, :]            # (P, T)
    out = np.zeros(d.size, dtype=complex)
    for k, mm in enumerate(m):
        if not np.any(c[:, k]):
            continue
        mm = int(mm)
        Jm = jv(abs(mm), x) * (-1.0) ** (mm if mm < 0 else 0)        # J_{-m} = (-1)^m J_m
        out += (2 * math.pi * 1j ** mm) * np.exp(1j * mm * phid.ravel()) * (Jm @ (wu * c[:, k]))
    return out.reshape(d.shape)


def element_directivity(el: ElementPattern, theta: float = 0.0, phi: float = 0.0) -> float:
    """4 pi p / (integral of p): one element on its own, toward (theta, phi) in radians."""
    return float(4 * math.pi * el(theta, phi) / power_kernel(el, np.zeros((1, 2)))[0].real)
