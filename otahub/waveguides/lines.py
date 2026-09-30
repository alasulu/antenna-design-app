"""Planar and coaxial transmission lines.

Closed-form synthesis and analysis for the line types an antenna feed network
actually uses: coax, microstrip, stripline and coplanar waveguide.
"""
from __future__ import annotations

import math

from scipy.special import ellipkm1

from ..core.constants import C0, ETA0, surface_resistance

_LOG = math.log


# ------------------------------------------------------------------- coaxial

def coax_impedance(a_inner: float, b_outer: float, eps_r: float = 1.0) -> float:
    """Coaxial characteristic impedance [ohm]: Z0 = eta0/(2*pi*sqrt(eps_r))*ln(b/a)."""
    if not 0 < a_inner < b_outer:
        raise ValueError(f"need 0 < a < b, got a={a_inner}, b={b_outer}")
    return ETA0 / (2.0 * math.pi * math.sqrt(eps_r)) * _LOG(b_outer / a_inner)


def coax_te11_cutoff(a_inner: float, b_outer: float, eps_r: float = 1.0) -> float:
    """Approximate TE11 cutoff [Hz], above which coax stops being single-mode.

    f_c ~ c / (pi * (a + b) * sqrt(eps_r)), from the mean circumference.
    """
    return C0 / (math.pi * (a_inner + b_outer) * math.sqrt(eps_r))


def coax_optimum_ratios() -> dict[str, float]:
    """The classic coax optima, all functions of b/a alone.

    Minimum attenuation at b/a = 3.5911 (Z0 = 76.7 ohm in air), maximum power
    handling at b/a = 1.6487 = sqrt(e) (Z0 = 30 ohm). 50 ohm is the historical
    compromise between the two, not an optimum of anything.
    """
    return {
        "min_attenuation_ratio": 3.5911,
        "min_attenuation_z0_air": ETA0 / (2 * math.pi) * _LOG(3.5911),
        "max_power_ratio": math.sqrt(math.e),
        "max_power_z0_air": ETA0 / (2 * math.pi) * _LOG(math.sqrt(math.e)),
    }


# ---------------------------------------------------------------- microstrip
#
# Hammerstad and Jensen, "Accurate models for microstrip computer-aided
# design", IEEE MTT-S 1980: continuous in w/h (the older two-branch Hammerstad
# forms jump 0.4% at w/h = 1, which leaves some impedances with no width at
# all), Z0 to 0.01% of the quasi-static solution for 0.01 < w/h < 100 and
# eps_eff to 0.2% for eps_r < 128, zero-thickness strip, no dispersion.

def _hj_z01(u: float) -> float:
    """Air-filled microstrip impedance at w/h = u (Hammerstad-Jensen)."""
    f = 6.0 + (2.0 * math.pi - 6.0) * math.exp(-(30.666 / u) ** 0.7528)
    return ETA0 / (2.0 * math.pi) * _LOG(f / u + math.sqrt(1.0 + (2.0 / u) ** 2))


def microstrip_eps_eff(w: float, h: float, eps_r: float) -> float:
    """Effective permittivity (Hammerstad-Jensen, static)."""
    u = w / h
    a = (1.0 + _LOG((u ** 4 + (u / 52.0) ** 2) / (u ** 4 + 0.432)) / 49.0
         + _LOG(1.0 + (u / 18.1) ** 3) / 18.7)
    b = 0.564 * ((eps_r - 0.9) / (eps_r + 3.0)) ** 0.053
    return (eps_r + 1) / 2 + (eps_r - 1) / 2 * (1 + 10 / u) ** (-a * b)


def microstrip_impedance(w: float, h: float, eps_r: float) -> float:
    """Microstrip characteristic impedance [ohm] (Hammerstad-Jensen).

    Z0 = Z01(w/h) / sqrt(eps_eff), continuous and monotonic in w/h; see the
    section note for its accuracy.
    """
    return _hj_z01(w / h) / math.sqrt(microstrip_eps_eff(w, h, eps_r))


def microstrip_width_for(z0: float, h: float, eps_r: float,
                         tol: float = 1e-10, max_iter: int = 200) -> float:
    """Synthesise the strip width [m] giving `z0`, by bisection on w/h.

    Bisection rather than the Wheeler closed form: it inverts exactly the same
    model used for analysis, so width -> impedance -> width round-trips to
    machine precision instead of carrying a second approximation's error.
    """
    lo, hi = 1e-4, 1e3                       # w/h bracket
    if microstrip_impedance(lo * h, h, eps_r) < z0:
        raise ValueError(f"Z0={z0} ohm is above what this substrate can reach")
    if microstrip_impedance(hi * h, h, eps_r) > z0:
        raise ValueError(f"Z0={z0} ohm is below what this substrate can reach")
    for _ in range(max_iter):
        mid = math.sqrt(lo * hi)             # geometric bisection: w/h spans decades
        z = microstrip_impedance(mid * h, h, eps_r)
        if abs(z - z0) < tol * z0:
            return mid * h
        if z > z0:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi) * h


def microstrip_guide_wavelength(f_hz: float, w: float, h: float, eps_r: float) -> float:
    return C0 / (f_hz * math.sqrt(microstrip_eps_eff(w, h, eps_r)))


# ----------------------------------------------------------------- stripline

def _K_ratio(k: float, kp: float) -> float:
    """K(k) / K(k'), the complete elliptic integrals by modulus, with the
    complementary modulus k' = sqrt(1 - k^2) passed in rather than formed.

    scipy's ellipkm1(p) is K at parameter m = 1 - p, so K(k) = ellipkm1(k'^2)
    and K(k') = ellipkm1(k^2): neither side ever computes 1 - k^2, which is
    where a wide stripline (k' -> 1 to the last bit) or a wide CPW centre
    conductor (k -> 1) used to lose every digit and return 0 or infinity.
    """
    return float(ellipkm1(kp * kp)) / float(ellipkm1(k * k))


def stripline_impedance(w: float, b: float, eps_r: float) -> float:
    """Stripline characteristic impedance [ohm], zero-thickness conductor.

    Z0 = (eta0 / (4*sqrt(eps_r))) * K(k)/K(k'), with k = sech(pi*w/(2b)).

    The ratio orientation is fixed by the limits: a wide strip must give a LOW
    impedance, and K(k)/K(k') is the form that shrinks as w grows. Checked
    against the familiar 30*pi/sqrt(eps_r) * b/(w + 0.441b) approximation, the
    two agree to 0.3% at w/b = 0.5.

    Stripline is a true TEM line, so eps_eff is exactly eps_r - no dispersion,
    unlike microstrip.
    """
    if w <= 0 or b <= 0:
        raise ValueError("w and b must be positive")
    arg = math.pi * w / (2.0 * b)
    k, kp = 1.0 / math.cosh(arg), math.tanh(arg)      # k' = tanh, exactly
    return (ETA0 / (4.0 * math.sqrt(eps_r))) * _K_ratio(k, kp)


# -------------------------------------------------------- coplanar waveguide

def cpw_impedance(w: float, s: float, eps_r: float) -> float:
    """CPW characteristic impedance [ohm] on an electrically thick substrate.

    `w` is the centre conductor width, `s` the gap to each ground plane.
    Z0 = (eta0 / (4*sqrt(eps_eff))) * K(k')/K(k), with k = w/(w + 2s).

    Note the ratio is the OPPOSITE way up from stripline: here a wide centre
    conductor drives k toward 1, which sends K(k) to infinity and the
    impedance to zero, as it must.

    eps_eff = (eps_r + 1)/2 because half the field is in air - which is why
    thick-substrate CPW is only weakly dispersive.
    """
    if w <= 0 or s <= 0:
        raise ValueError("w and s must be positive")
    k = w / (w + 2.0 * s)
    kp = 2.0 * math.sqrt(s * (w + s)) / (w + 2.0 * s)   # the same, cancellation-free
    eps_eff = (eps_r + 1.0) / 2.0
    return (ETA0 / (4.0 * math.sqrt(eps_eff))) * _K_ratio(kp, k)


def cpw_eps_eff(eps_r: float) -> float:
    """Quasi-static CPW effective permittivity on a thick substrate."""
    return (eps_r + 1.0) / 2.0


def quarter_wave_transformer(z_source: float, z_load: float) -> float:
    """Impedance [ohm] of the quarter-wave section matching source to load."""
    if z_source <= 0 or z_load <= 0:
        raise ValueError("impedances must be positive")
    return math.sqrt(z_source * z_load)
