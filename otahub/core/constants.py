"""Physical constants. SI throughout — the whole engine is SI-internal."""
from __future__ import annotations

import math

#: Speed of light in vacuum [m/s] (exact, SI definition)
C0 = 2.99792458e8
#: Vacuum permeability [H/m] (CODATA 2018; no longer exactly 4e-7*pi)
MU0 = 1.25663706212e-6
#: Vacuum permittivity [F/m]
EPS0 = 1.0 / (MU0 * C0 * C0)
#: Impedance of free space [ohm]
ETA0 = math.sqrt(MU0 / EPS0)  # 376.730313412...

#: Conductivity of common metals [S/m] at 20 C
CONDUCTIVITY = {
    "copper": 5.8e7,
    "annealed_copper": 5.80e7,
    "aluminium": 3.77e7,
    "aluminum": 3.77e7,
    "silver": 6.30e7,
    "gold": 4.10e7,
    "brass": 1.57e7,
    "steel": 1.45e6,
    "pec": float("inf"),
}


def wavelength(f_hz: float, eps_r: float = 1.0, mu_r: float = 1.0) -> float:
    """Wavelength [m] in a medium at frequency `f_hz`."""
    if f_hz <= 0:
        raise ValueError(f"frequency must be positive, got {f_hz}")
    return C0 / (f_hz * math.sqrt(eps_r * mu_r))


def wavenumber(f_hz: float, eps_r: float = 1.0, mu_r: float = 1.0) -> float:
    """Free-space (or in-medium) wavenumber k = 2*pi/lambda [rad/m]."""
    return 2.0 * math.pi / wavelength(f_hz, eps_r, mu_r)


def skin_depth(f_hz: float, sigma: float, mu_r: float = 1.0) -> float:
    """Skin depth [m]. delta = 1/sqrt(pi*f*mu*sigma)."""
    if sigma == float("inf"):
        return 0.0
    return 1.0 / math.sqrt(math.pi * f_hz * MU0 * mu_r * sigma)


def surface_resistance(f_hz: float, sigma: float, mu_r: float = 1.0) -> float:
    """Surface resistance Rs = sqrt(pi*f*mu/sigma) [ohm/square]."""
    if sigma == float("inf"):
        return 0.0
    return math.sqrt(math.pi * f_hz * MU0 * mu_r / sigma)
