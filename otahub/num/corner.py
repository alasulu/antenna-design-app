"""A dipole in an infinite PEC corner, solved exactly by images and the MoM.

A corner of angle 180/n degrees turns a dipole on its bisector into 2n sources
on a circle of radius S with alternating signs; the image set leaves no
tangential E on either plate. Solving that set as one method-of-moments
problem - every dipole driven, with its image sign - gives the real current on
the real dipole, not an assumed sinusoid, so the driving-point impedance and
the pattern both come out of Maxwell rather than induced EMF.

The image system radiates the same power into each of its 2n wedges, so the
real reflector's directivity is 2n times the image system's full-space value,
and the power the image system radiates is 2n times what the real feed
delivers - which `CornerSolution.power_balance` checks.

Geometry in wavelengths, as in `mom`.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from . import mom

__all__ = ["CornerSolution", "solve_corner", "image_set"]


def image_set(spacing: float, corner_deg: float) -> list[tuple[float, float, int]]:
    """(x, y, sign) of the real dipole (first) and its images."""
    n = 360.0 / corner_deg
    if abs(n - round(n)) > 1e-9 or round(n) % 2:
        raise ValueError("the corner angle must divide 180 degrees for images to close")
    n = int(round(n))
    return [(spacing * math.cos(2 * math.pi * i / n), spacing * math.sin(2 * math.pi * i / n),
             (-1) ** i) for i in range(n)]


@dataclass
class CornerSolution:
    solution: mom.MoMSolution
    wedges: int

    @property
    def input_impedance(self) -> complex:
        return self.solution.input_impedance

    def directivity_boresight(self, nth: int = 90, nph: int = 180) -> float:
        _, total, _ = mom.pattern_power(self.solution, nth, nph)
        e_th, e_ph = mom.far_field(self.solution, np.array(math.pi / 2), np.array(0.0))
        u = float(abs(e_th) ** 2 + abs(e_ph) ** 2)
        return self.wedges * 4.0 * math.pi * u / total

    def power_balance(self, nth: int = 90, nph: int = 180) -> float:
        """Image-system radiated power over (wedges x real feed power), normalised by
        the same ratio for an isolated dipole - 1 when the images are right."""
        _, total, _ = mom.pattern_power(self.solution, nth, nph)
        iso = mom.solve(mom.dipole(0.5, 1e-3, 40))
        _, t_iso, _ = mom.pattern_power(iso, nth, nph)
        cal = t_iso / iso.circuit_power
        return total / (self.wedges * self.solution.circuit_power) / cal


def solve_corner(spacing: float, corner_deg: float = 90.0, radius: float = 1e-3,
                 length: float = 0.5, segments: int = 40, exact: bool = True) -> CornerSolution:
    """Vertical dipole at `spacing` from the vertex on the bisector, lambda units."""
    images = image_set(spacing, corner_deg)
    wires = []
    for x, y, _ in images:
        z = np.linspace(-length / 2, length / 2, segments + 1)
        wires.append(mom.Wire(np.stack([np.full_like(z, x), np.full_like(z, y), z], 1), radius))
    model = mom.WireModel(wires)
    nodes = np.array([model.node_of(b) for b in range(model.n_basis)])
    feeds = [int(np.argmin(np.linalg.norm(nodes - np.array([x, y, 0.0]), axis=1)))
             for x, y, _ in images]
    Z = mom.impedance_matrix(model, exact)
    V = np.zeros(model.n_basis, dtype=complex)
    for k, (_, _, sign) in zip(feeds, images):
        V[k] = sign
    currents = np.linalg.solve(Z, V)
    return CornerSolution(mom.MoMSolution(model, currents, feeds[0]), len(images))
