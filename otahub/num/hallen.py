"""Hallen's integral equation for a straight centre-fed dipole - a second opinion.

    integral I(z') exp(-jkR)/(4 pi R) dz' = -(j/eta) [ B cos(kz) + (V/2) sin(k|z|) ]

Deliberately shares nothing with `mom.py`: a different integral equation, an
extra scalar unknown B in place of a divergence term, and point matching in
place of Galerkin testing. It exists to arbitrate when the method of moments
disagrees with a closed form, which it has done twice - once to show that the
textbook 73.08 ohm is an induced-EMF figure rather than a driving-point one,
and once to confirm a short whip's reactance was ten times what a spec said.

It has a well-known weakness worth knowing: with the reduced kernel it becomes
ill-conditioned on fat wires under fine meshes (at a = 0.005 lambda it collapses
by N = 400). Keep it to thin wire, where it is reliable.

Geometry in wavelengths, k = 2*pi, as in `mom`.
"""
from __future__ import annotations

import math

import numpy as np

from ..core.constants import ETA0

__all__ = ["hallen_dipole"]

K = 2.0 * math.pi


def hallen_dipole(length: float, a: float, n: int = 120, nq: int = 64) -> complex:
    """Driving-point impedance of a delta-gap-fed dipole [ohm].

    Triangular basis on the interior nodes, current forced to zero at the tips,
    matched at the interior nodes plus one tip.
    """
    h = 0.5 * length
    z = np.linspace(-h, h, n + 1)
    d = z[1] - z[0]
    nb = n - 1
    zm = np.concatenate([z[1:n], [h]])
    x, w = np.polynomial.legendre.leggauss(nq)
    A = np.zeros((nb + 1, nb + 1), dtype=complex)
    for j in range(nb):
        c = z[j + 1]
        for lo, hi, rise in ((c - d, c, True), (c, c + d, False)):
            s = 0.5 * (hi - lo) * (x + 1.0) + lo
            f = (s - (c - d)) / d if rise else ((c + d) - s) / d
            R = np.sqrt((zm[:, None] - s[None, :]) ** 2 + a * a)
            A[:, j] += (0.5 * (hi - lo)) * (
                (np.exp(-1j * K * R) / (4 * math.pi * R))
                * f[None, :] * w[None, :]).sum(axis=1)
    A[:, nb] = (1j / ETA0) * np.cos(K * zm)
    rhs = -(1j / ETA0) * 0.5 * np.sin(K * np.abs(zm))
    sol = np.linalg.solve(A, rhs)
    return complex(1.0 / np.concatenate([[0.0], sol[:nb], [0.0]])[n // 2])
