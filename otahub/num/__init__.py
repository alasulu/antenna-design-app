"""Numerical reference models.

Nothing in here is an archetype. These are independent solvers used to CHECK
the closed forms the specs ship, so that a number written into a spec has been
confirmed by a route that shares no algebra with it.
"""
from __future__ import annotations

from . import loop_modal
from .mom import (
    Wire, WireModel, MoMSolution, dipole, loop, folded_dipole_wire,
    impedance_matrix, solve, input_impedance, far_field, pattern_power,
    radiated_power, directivity, directivity_towards,
)

__all__ = [
    "loop_modal",
    "Wire", "WireModel", "MoMSolution", "dipole", "loop", "folded_dipole_wire",
    "impedance_matrix", "solve", "input_impedance", "far_field",
    "pattern_power", "radiated_power", "directivity", "directivity_towards",
]
