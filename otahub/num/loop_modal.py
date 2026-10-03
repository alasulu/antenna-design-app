"""Fourier-mode solution of the circular loop - an arbiter for the MoM.

Derived rather than copied, so it can be read and checked. Write the current on
a loop of radius b (wire radius a) as I(phi) = sum_n I_n exp(j n phi). The
phi-directed field of one mode separates, because every integral around the
ring turns into a Fourier coefficient of the kernel:

    A_phi,n  = (mu b / 8 pi) I_n [ Ghat_{n-1} + Ghat_{n+1} ]
    Phi_n    = -(n I_n / 4 pi eps omega) Ghat_n
    Ghat_m   = integral over psi of exp(-j m psi) exp(-jkR)/R dpsi,
    R(psi)   = sqrt(4 b^2 sin^2(psi/2) + a^2)

A delta gap of 1 V has E_inc,n = 1/(2 pi b) for every n, and imposing
E_phi^scat + E_inc = 0 mode by mode gives

    I_n = 4 k / ( j eta [ k^2 b^2 (Ghat_{n-1} + Ghat_{n+1}) - 2 n^2 Ghat_n ] )

so Y_in = sum_n I_n. The cos(phi) term dominating at one wavelength is what
puts the beam broadside; a uniform current would peak in the loop's own plane.

Why it is worth having: it discretises no geometry at all, so it shares nothing
with `mom.py` but Maxwell. Where the two agree, the answer is the physics
rather than either code. Geometry is in wavelengths, k = 2*pi, as in `mom`.
"""
from __future__ import annotations

import functools
import math

import numpy as np

from ..core.constants import ETA0

__all__ = ["kernel_coefficients", "mode_currents", "input_impedance",
           "current_at", "resonant_circumference", "loop_current_resistance",
           "tuned_loop"]

K = 2.0 * math.pi
#: Default kernel quadrature. R(psi) peaks over a width of about a/b, so a thin
#: wire needs a lot of points; 20000 is converged to five digits at a/lambda =
#: 1e-4 and 2000 would be 25% wrong in the reactance there.
NQ = 20000


@functools.lru_cache(maxsize=512)
def kernel_coefficients(b: float, a: float, nmax: int, nq: int = NQ) -> tuple:
    """Ghat_m for m = 0 .. nmax+1. The integrand is periodic and smooth, so the
    plain trapezoid rule is spectrally accurate once the peak is resolved."""
    psi = np.arange(nq) * 2.0 * math.pi / nq
    r = np.sqrt(4.0 * b * b * np.sin(0.5 * psi) ** 2 + a * a)
    g = np.exp(-1j * K * r) / r
    return tuple((g * np.exp(-1j * m * psi)).sum() * (2.0 * math.pi / nq)
                 for m in range(nmax + 2))


def mode_currents(circumference: float, a: float, nmax: int = 60,
                  loss: float = 0.0) -> dict[int, complex]:
    """I_n for a 1 V delta gap at phi = 0. `loss` is the conductor's series
    resistance per wavelength of wire (ohm), uniform round the loop - a skin-effect
    wire's R_s/(2 pi a) in these units. It enters every mode alike, E_inc,n =
    (Z_n + loss) I_n, because the tangential field at a resistive wire is that
    resistance times its current wherever the current is."""
    b = circumference / (2.0 * math.pi)
    g = kernel_coefficients(b, a, nmax)
    out: dict[int, complex] = {}
    for n in range(-nmax, nmax + 1):
        den = ((K * b) ** 2 * (g[abs(n - 1)] + g[abs(n + 1)])
               - 2.0 * n * n * g[abs(n)])
        z_n = 1j * ETA0 * den / (4.0 * K) / (2.0 * math.pi * b)    # ohm per wavelength
        out[n] = 1.0 / (2.0 * math.pi * b) / (z_n + loss)
    return out


def input_impedance(circumference: float, a: float, nmax: int = 60, loss: float = 0.0) -> complex:
    return 1.0 / sum(mode_currents(circumference, a, nmax, loss).values())


def loop_current_resistance(circumference: float, a: float, nmax: int = 60) -> float:
    """Radiation resistance referred to the loop's mean-square current, 2 P / <|I|^2>.

    The resistance that sets efficiency: a conductor loses r <|I|^2> per unit
    length wherever the current is, so efficiency is R / (R + r 2 pi b) with THIS
    R, not the terminal resistance - the two part as the current stops being
    uniform (at 0.3 wavelengths round the far side carries 1.8 times the feed's
    current, and the terminal resistance is 2.4 times this one). It is also the
    robust one: the delta gap's own capacitance, which the terminal value takes
    up with the mode count, carries almost no current round the loop."""
    modes = mode_currents(circumference, a, nmax)
    y = sum(modes.values())
    return float(y.real / sum(abs(v) ** 2 for v in modes.values()))


def tuned_loop(circumference: float, a: float, rs_ohm: float = 0.0, nmax: int = 60,
               h: float = 1e-4) -> dict:
    """The loop tuned by a series capacitor at its terminals: Q from the impedance's
    slope (Yaghjian and Best's Q_Z, the capacitor's reactance slope included),
    radiation efficiency from the modal powers. `rs_ohm` is the conductor's surface
    resistance at this frequency; it is scaled as sqrt(f) across the derivative, as
    a skin-effect wire's is. Geometry in wavelengths at the design frequency."""
    def z(scale: float) -> complex:
        loss = rs_ohm * math.sqrt(scale) / (2.0 * math.pi * a * scale)
        return input_impedance(circumference * scale, a * scale, nmax, loss)
    z0 = z(1.0)
    dz = (z(1.0 + h) - z(1.0 - h)) / (2.0 * h)                  # dZ/dw times w
    q = abs(dz + 1j * z0.imag) / (2.0 * z0.real)
    modes = mode_currents(circumference, a, nmax, rs_ohm / (2.0 * math.pi * a))
    p_in = 0.5 * (1.0 / z0).real                                # 1 V across the gap
    p_loss = 0.5 * rs_ohm / (2.0 * math.pi * a) * circumference * sum(abs(v) ** 2 for v in modes.values())
    return {"z": z0, "q": float(q), "efficiency": float(1.0 - p_loss / p_in)}


def current_at(circumference: float, a: float, phi, nmax: int = 60):
    """I(phi) for a 1 V delta gap, on an array of angles."""
    modes = mode_currents(circumference, a, nmax)
    phi = np.asarray(phi, dtype=float)
    return sum(v * np.exp(1j * n * phi) for n, v in modes.items())


def resonant_circumference(a: float, nmax: int = 60, steps: int = 24,
                           lo: float = 0.95, hi: float = 1.30):
    """(circumference in wavelengths, resistance in ohm) where X crosses zero.

    A thicker conductor resonates at a LONGER circumference, which is the
    opposite of a dipole and the opposite of what `one_wavelength_circular_loop`
    used to claim. A bracket with no crossing in it is an error, not its end.
    """
    z_lo, z_hi = input_impedance(lo, a, nmax), input_impedance(hi, a, nmax)
    if z_lo.imag == 0.0:
        return lo, z_lo.real
    if z_hi.imag == 0.0:
        return hi, z_hi.real
    below = z_lo.imag < 0
    if (z_hi.imag < 0) == below:
        raise ValueError(f"no resonance between C = {lo} and {hi} wavelengths "
                         f"for a = {a} wavelengths; widen the bracket")
    for _ in range(steps):
        mid = 0.5 * (lo + hi)
        if (input_impedance(mid, a, nmax).imag < 0) == below:
            lo = mid
        else:
            hi = mid
    return lo, input_impedance(lo, a, nmax).real
