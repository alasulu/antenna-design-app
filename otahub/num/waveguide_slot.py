"""Shunt conductance of a longitudinal broad-wall slot, from first principles.

An arbiter for Stevenson's formula. A slot of length L and width w sits in the
broad wall of an a x b guide, offset x1 from the centreline, with the aperture
field E_x = (V/w) cos(pi z/L) that a resonant slot carries. Three steps, each
done by quadrature rather than by the closed form it checks:

1. Reciprocity. The slot's equivalent magnetic current M_z excites TE10 in both
   directions with amplitude c = (1/P_n) * integral of H_z(-+) M_z over the
   slot, P_n = 2 * integral of (E_t x H_t).z over the cross-section. Because it
   couples through H_z, which is even in z for both directions, the two
   scattered waves are equal: the slot is a SHUNT element.
2. Radiation. Outside, the same current - doubled by the ground plane -
   radiates into a half space. Its one-sided conductance G_ext = 2 P_ext / V^2
   is integrated over the hemisphere.
3. Power balance for a shunt conductance g on the line: the scattered wave is
   (g/2) times the total wave at the slot, and the power g absorbs is what the
   slot radiates. Hence g = 2 E_s^2 a b / (Z_TE V^2 G_ext).

Nothing assumes Stevenson's algebra, which is what makes it a check on it: with
L = lambda/2 it reproduces g = 2.09 (lambda_g/lambda)(a/b)
cos^2(pi lambda/(2 lambda_g)) sin^2(pi x1/a) to the rounding of the 2.09, and
not the printing with the wavelength ratio inverted. SI units throughout.
"""
from __future__ import annotations

import math

import numpy as np

from ..core.constants import C0, ETA0, MU0

__all__ = ["guide_wavelength", "one_sided_conductance", "shunt_conductance"]

_GL = np.polynomial.legendre.leggauss


def guide_wavelength(f: float, a: float) -> float:
    lam = C0 / f
    return lam / math.sqrt(1.0 - (lam / (2.0 * a)) ** 2)


def _slot_axis(length: float, n: int):
    z, w = _GL(n)
    return 0.5 * length * z, 0.5 * length * w


def one_sided_conductance(f: float, length: float | None = None, n: int = 400) -> float:
    """G_ext [S] of a slot with a cosine aperture field, radiating into one half space."""
    lam = C0 / f
    k = 2.0 * math.pi / lam
    length = lam / 2.0 if length is None else length
    zs, wz = _slot_axis(length, n)
    im = 2.0 * np.cos(math.pi * zs / length)            # V = 1, image-doubled
    th, wt = _GL(n)
    th = 0.5 * math.pi * (th + 1.0)
    F = (im[None, :] * np.exp(1j * k * np.outer(np.cos(th), zs)) * wz[None, :]).sum(axis=1)
    p_full = 2.0 * math.pi * ((k ** 2 * np.sin(th) ** 3 * np.abs(F) ** 2
                               / (32.0 * math.pi ** 2 * ETA0)) * wt).sum() * 0.5 * math.pi
    return 2.0 * (0.5 * p_full)


def shunt_conductance(f: float, a: float, b: float, x1: float, w: float | None = None,
                      length: float | None = None, n: int = 400) -> float:
    """Normalised shunt conductance g of the slot, by reciprocity and power balance."""
    lam = C0 / f
    k = 2.0 * math.pi / lam
    kc = math.pi / a
    if k <= kc:
        raise ValueError("TE10 is cut off at this frequency")
    beta = math.sqrt(k * k - kc * kc)
    omega = 2.0 * math.pi * f
    zte = omega * MU0 / beta
    w = 1e-6 * a if w is None else w
    length = lam / 2.0 if length is None else length
    # mode normalisation for unit E0: E_y = sin(pi x/a), H_x = -sin(pi x/a)/Z_TE
    xg, wg = _GL(200)
    xs = 0.5 * a * (xg + 1.0)
    pn = 2.0 * ((np.sin(math.pi * xs / a) ** 2 / zte) * wg).sum() * 0.5 * a * b
    # overlap of |H_z| = (pi/(omega mu a)) cos(pi x/a) with M_z = (V/w) cos(pi z/L)
    xq, wq = _GL(40)
    xsl = a / 2.0 + x1 + 0.5 * w * xq
    hz = (math.pi / (omega * MU0 * a)) * np.cos(math.pi * xsl / a)
    zs, wz = _slot_axis(length, n)
    mz = np.cos(math.pi * zs / length) / w
    overlap = abs((hz * wq).sum() * 0.5 * w * ((mz * np.exp(1j * beta * zs)) * wz).sum())
    e_s = overlap / pn
    return 2.0 * e_s ** 2 * a * b / (zte * one_sided_conductance(f, length, n))
