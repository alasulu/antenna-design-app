"""Numerical reference models.

Nothing in here is an archetype. These are independent solvers used to CHECK
the closed forms the specs ship, so that a number written into a spec has been
confirmed by a route that shares no algebra with it.
"""
from __future__ import annotations

from . import hallen, loop_modal
from .mom import (
    Wire, WireModel, MoMSolution, dipole, loop, arc, halo, folded_dipole_wire,
    helix_over_ground, lpda_model, rhombic_model,
    top_hat_monopole, tl_admittance, NetworkSolution,
    solve_network,
    impedance_matrix, solve, input_impedance, far_field, pattern_power,
    radiated_power, directivity, directivity_towards,
    vswr, resonant_scale, antenna_q, vswr_bandwidth,
)

__all__ = [
    "hallen", "loop_modal",
    "Wire", "WireModel", "MoMSolution", "dipole", "loop", "arc", "halo",
    "folded_dipole_wire", "helix_over_ground", "lpda_model",
    "rhombic_model", "top_hat_monopole",
    "tl_admittance", "NetworkSolution", "solve_network",
    "impedance_matrix", "solve", "input_impedance", "far_field",
    "pattern_power", "radiated_power", "directivity", "directivity_towards",
    "vswr", "resonant_scale", "antenna_q", "vswr_bandwidth",
]
