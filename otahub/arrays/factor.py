"""Array factors, directivity and grating-lobe limits for linear arrays.

Patterns are returned as :class:`otahub.core.pattern.Pattern` objects so the
directivity, beamwidth and sidelobe machinery already verified against
closed-form cases applies unchanged to arrays.
"""
from __future__ import annotations

import math

import numpy as np

from ..core.pattern import Pattern, make_grid
from .tapers import taper_efficiency


def array_factor(weights, d_over_lambda: float, theta: np.ndarray,
                 scan_deg: float = 90.0) -> np.ndarray:
    """Complex array factor of a linear array along z.

    `theta` is measured from the +z axis. `scan_deg` is the desired main-beam
    direction in the same convention, so 90 degrees is broadside and 0 is
    endfire; the inter-element phase is set to steer there.
    """
    w = np.asarray(weights, dtype=float)
    n = len(w)
    positions = np.arange(n) - (n - 1) / 2.0
    psi_scan = 2.0 * math.pi * d_over_lambda * math.cos(math.radians(scan_deg))
    phase = 2.0 * math.pi * d_over_lambda * np.cos(theta)[:, None] - psi_scan
    return (w[None, :] * np.exp(1j * positions[None, :] * phase)).sum(axis=1)


def array_pattern(weights, d_over_lambda: float, scan_deg: float = 90.0,
                  n_theta: int = 3601) -> Pattern:
    """Power pattern of the array factor alone (isotropic elements)."""
    theta, phi = make_grid(n_theta, 181)
    af = array_factor(weights, d_over_lambda, theta, scan_deg)
    u = np.abs(af) ** 2
    return Pattern(theta, phi, np.repeat(u[:, None], phi.size, axis=1))


def broadside_directivity(weights, d_over_lambda: float) -> float:
    """Directivity of a broadside linear array of isotropic elements.

    For d = lambda/2 and uniform weights this is exactly N. A taper reduces it
    by the taper efficiency, which is the whole cost of low sidelobes.
    """
    w = np.asarray(weights, dtype=float)
    n = len(w)
    positions = (np.arange(n) - (n - 1) / 2.0) * d_over_lambda
    num = w.sum() ** 2
    # sum_m sum_n a_m a_n sinc(2*(z_m - z_n)) with numpy's normalised sinc
    delta = positions[:, None] - positions[None, :]
    den = float((w[:, None] * w[None, :] * np.sinc(2.0 * delta)).sum())
    if den <= 0:
        raise ValueError("degenerate array: no radiated power")
    return float(num / den)


def grating_lobe_free_spacing(scan_deg: float = 90.0) -> float:
    """Largest d/lambda with no grating lobe in real space.

    d/lambda < 1 / (1 + |cos(scan)|) in the theta-from-z convention. Broadside
    allows just under 1; scanning to endfire tightens it to 1/2, which is the
    reason phased arrays are built on a half-wave lattice.
    """
    return 1.0 / (1.0 + abs(math.cos(math.radians(scan_deg))))


def has_grating_lobe(d_over_lambda: float, scan_deg: float = 90.0) -> bool:
    return d_over_lambda >= grating_lobe_free_spacing(scan_deg)


def beam_broadening_factor(weights) -> float:
    """Rough beamwidth penalty of a taper relative to uniform excitation."""
    return 1.0 / math.sqrt(taper_efficiency(weights))


def summarise(weights, d_over_lambda: float, scan_deg: float = 90.0) -> dict:
    """Everything an array designer wants in one call."""
    from ..core.pattern import first_sidelobe_db, hpbw_deg

    pattern = array_pattern(weights, d_over_lambda, scan_deg)
    out = {
        "elements": len(np.asarray(weights)),
        "d_over_lambda": d_over_lambda,
        "scan_deg": scan_deg,
        "taper_efficiency": taper_efficiency(weights),
        "sidelobe_db": first_sidelobe_db(pattern),
        "grating_lobe": has_grating_lobe(d_over_lambda, scan_deg),
        "max_spacing_no_grating": grating_lobe_free_spacing(scan_deg),
    }
    try:
        out["hpbw_deg"] = hpbw_deg(pattern)
    except ValueError:
        out["hpbw_deg"] = float("nan")
    try:
        out["directivity_dbi"] = 10 * math.log10(
            broadside_directivity(weights, d_over_lambda))
    except ValueError:
        out["directivity_dbi"] = float("nan")
    return out
