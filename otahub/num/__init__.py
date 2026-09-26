"""Numerical reference models.

Nothing in here is an archetype. These are independent solvers used to CHECK
the closed forms the specs ship, so that a number written into a spec has been
confirmed by a route that shares no algebra with it.
"""
from __future__ import annotations

from . import bor, corner, dra, hallen, horn_pattern, lens, loop_modal, paraboloid, patch_cavity, patch_q, waveguide_slot, zone_plate
from .mom import (
    Wire, WireModel, MoMSolution, dipole, loop, arc, halo, folded_dipole_wire,
    helix_over_ground, lpda_model, rhombic_model,
    top_hat_monopole, bicone_cage, tl_admittance, NetworkSolution,
    solve_network,
    impedance_matrix, solve, input_impedance, far_field, pattern_power,
    radiated_power, directivity, directivity_towards,
    vswr, resonant_scale, antenna_q, vswr_bandwidth,
)

__all__ = [
    "bor", "corner", "dra", "hallen", "horn_pattern", "lens", "loop_modal", "paraboloid", "patch_cavity", "patch_q", "waveguide_slot", "zone_plate",
    "Wire", "WireModel", "MoMSolution", "dipole", "loop", "arc", "halo",
    "folded_dipole_wire", "helix_over_ground", "lpda_model",
    "rhombic_model", "top_hat_monopole", "bicone_cage",
    "tl_admittance", "NetworkSolution", "solve_network",
    "impedance_matrix", "solve", "input_impedance", "far_field",
    "pattern_power", "radiated_power", "directivity", "directivity_towards",
    "vswr", "resonant_scale", "antenna_q", "vswr_bandwidth",
]
