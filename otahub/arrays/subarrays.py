"""Subarray architectures: steering applied per subarray instead of per element.

A large array is often steered by one phase shifter per subarray (or one time
delay per subarray over element-level phase shifters). With phase applied only
at the subarray level, every element of a subarray shares its centre's phase,
and for identical subarrays the pattern factorises exactly:

    AF(u) = S(u) * sum_s exp(j k c_s.(u - u0))

S the (unsteered) subarray pattern and c_s the subarray centres. The second
factor is an array at the subarray period, so it grows grating lobes -
quantisation lobes - at u0 - m lambda/D_sub; the first, still pointing at
broadside, weighs them: the lobe sits S(u_q)/S(u0) below the beam in voltage,
and the beam itself is down by S(u0) - the scan loss of steering this way.

`subarray_steering` builds the excitation; `directivity_toward` evaluates any
complex excitation toward a direction exactly (the planar kernel, with an
element pattern if given). Positions in wavelengths, theta from the normal.
"""
from __future__ import annotations

import math

import numpy as np

__all__ = ["rect_subarray_groups", "subarray_centres", "subarray_steering", "directivity_toward"]


def rect_subarray_groups(nx: int, ny: int, sx: int, sy: int) -> np.ndarray:
    """Group index of each element of an nx x ny rectangular lattice (x-fastest,
    as `planar.rectangular_lattice` orders it) cut into sx x sy blocks."""
    if nx % sx or ny % sy:
        raise ValueError("subarray size must divide the array size")
    i, j = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    return ((i // sx) * (ny // sy) + (j // sy)).ravel()


def subarray_centres(positions, groups) -> np.ndarray:
    pos = np.asarray(positions, dtype=float)
    g = np.asarray(groups)
    return np.array([pos[g == k].mean(axis=0) for k in range(int(g.max()) + 1)])


def subarray_steering(positions, groups, weights=None, scan_theta_deg: float = 0.0, scan_phi_deg: float = 0.0,
                      level: str = "subarray") -> np.ndarray:
    """Complex element excitations that point the beam at (scan_theta, scan_phi):
    level='element' puts the exact phase on every element, 'subarray' gives every
    element its subarray centre's phase (one phase shifter per subarray)."""
    pos = np.asarray(positions, dtype=float)
    w = np.ones(len(pos)) if weights is None else np.asarray(weights, dtype=float)
    st = math.sin(math.radians(scan_theta_deg))
    u0 = np.array([st * math.cos(math.radians(scan_phi_deg)), st * math.sin(math.radians(scan_phi_deg))])
    if level == "element":
        where = pos
    elif level == "subarray":
        where = subarray_centres(pos, groups)[np.asarray(groups)]
    else:
        raise ValueError("level must be 'element' or 'subarray'")
    return w * np.exp(-2j * math.pi * (where @ u0))


def directivity_toward(positions, excitation, theta_deg: float, phi_deg: float = 0.0, element=None) -> float:
    """4 pi U / P toward (theta, phi) for an arbitrary complex excitation."""
    from .elements import isotropic, power_kernel
    el = element or isotropic()
    pos = np.asarray(positions, dtype=float)
    a = np.asarray(excitation, dtype=complex)
    diff = np.round((pos[:, None, :] - pos[None, :, :]).reshape(-1, 2), 9)
    uniq, inv = np.unique(diff, axis=0, return_inverse=True)
    Kd = power_kernel(el, uniq)[inv.ravel()].reshape(len(pos), len(pos))
    den = float(np.real(a @ Kd @ np.conj(a)))
    t, f = math.radians(theta_deg), math.radians(phi_deg)
    u = np.array([math.sin(t) * math.cos(f), math.sin(t) * math.sin(f)])
    af = np.sum(a * np.exp(2j * math.pi * (pos @ u)))
    return 4 * math.pi * float(el(t, f)) * abs(af) ** 2 / den
