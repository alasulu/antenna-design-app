"""Network parameters and mismatch figures.

One-port and two-port conversions between S, Z, Y and ABCD, plus the
mismatch quantities an antenna engineer quotes daily. Everything is SI and
complex-valued; reference impedance defaults to 50 ohm but is never assumed
silently - it is an explicit argument everywhere it matters.
"""
from __future__ import annotations

import cmath
import math
from dataclasses import dataclass

import numpy as np


# ------------------------------------------------------------ one-port math

def reflection_coefficient(z_load: complex, z0: float = 50.0) -> complex:
    """Gamma = (ZL - Z0) / (ZL + Z0). An open circuit (infinite ZL) is +1."""
    z_load = complex(z_load)
    if cmath.isinf(z_load):
        return 1.0 + 0j
    if z_load == -z0:
        raise ValueError("ZL = -Z0 is a singular load")
    return (z_load - z0) / (z_load + z0)


def impedance_from_gamma(gamma: complex, z0: float = 50.0) -> complex:
    """Invert the reflection coefficient back to an impedance."""
    gamma = complex(gamma)
    if abs(1.0 - gamma) < 1e-15:
        return complex(float("inf"), 0.0)
    return z0 * (1.0 + gamma) / (1.0 - gamma)


def vswr(gamma: complex | float) -> float:
    """Voltage standing wave ratio from a reflection coefficient."""
    mag = abs(complex(gamma))
    if mag >= 1.0 - 1e-15:
        return float("inf")
    return (1.0 + mag) / (1.0 - mag)


def gamma_from_vswr(s: float) -> float:
    """Reflection magnitude corresponding to a VSWR."""
    if s < 1.0:
        raise ValueError(f"VSWR cannot be below 1, got {s}")
    if not math.isfinite(s):
        return 1.0
    return (s - 1.0) / (s + 1.0)


def return_loss_db(gamma: complex | float) -> float:
    """Return loss [dB], reported POSITIVE as is conventional.

    A perfect match is infinite return loss. Note the sign convention trap:
    S11 in dB is negative, return loss is its magnitude.
    """
    mag = abs(complex(gamma))
    if mag <= 1e-15:
        return float("inf")
    return -20.0 * math.log10(mag)


def mismatch_loss_db(gamma: complex | float) -> float:
    """Power lost to reflection [dB, positive]: -10*log10(1 - |Gamma|^2)."""
    mag2 = abs(complex(gamma)) ** 2
    if mag2 >= 1.0 - 1e-15:
        return float("inf")
    return -10.0 * math.log10(1.0 - mag2)


def s11_db(z_load: complex, z0: float = 50.0) -> float:
    """S11 in dB (negative for a load that absorbs anything)."""
    return -return_loss_db(reflection_coefficient(z_load, z0))


@dataclass(frozen=True, slots=True)
class MatchQuality:
    """The family of mismatch numbers, computed once and consistent."""
    z_load: complex
    z0: float
    gamma: complex
    vswr: float
    return_loss_db: float
    mismatch_loss_db: float

    def __str__(self) -> str:
        return (f"Z = {self.z_load.real:.2f}{self.z_load.imag:+.2f}j ohm, "
                f"VSWR {self.vswr:.3f}, RL {self.return_loss_db:.2f} dB, "
                f"mismatch loss {self.mismatch_loss_db:.3f} dB")


def match_quality(z_load: complex, z0: float = 50.0) -> MatchQuality:
    gamma = reflection_coefficient(z_load, z0)
    return MatchQuality(complex(z_load), z0, gamma, vswr(gamma),
                        return_loss_db(gamma), mismatch_loss_db(gamma))


# ----------------------------------------------------------- two-port math

def z_to_s(z: np.ndarray, z0: float = 50.0) -> np.ndarray:
    """Two-port Z -> S. S = (Z - Z0*I)(Z + Z0*I)^-1 for real Z0."""
    z = np.asarray(z, dtype=complex)
    _check_2x2(z)
    eye = np.eye(2, dtype=complex) * z0
    return (z - eye) @ np.linalg.inv(z + eye)


def s_to_z(s: np.ndarray, z0: float = 50.0) -> np.ndarray:
    """Two-port S -> Z."""
    s = np.asarray(s, dtype=complex)
    _check_2x2(s)
    eye = np.eye(2, dtype=complex)
    return z0 * (eye + s) @ np.linalg.inv(eye - s)


def z_to_y(z: np.ndarray) -> np.ndarray:
    return np.linalg.inv(np.asarray(z, dtype=complex))


def y_to_z(y: np.ndarray) -> np.ndarray:
    return np.linalg.inv(np.asarray(y, dtype=complex))


def z_to_abcd(z: np.ndarray) -> np.ndarray:
    """Z -> ABCD. Requires Z21 != 0 (a reciprocal, non-degenerate network)."""
    z = np.asarray(z, dtype=complex)
    _check_2x2(z)
    z11, z12, z21, z22 = z[0, 0], z[0, 1], z[1, 0], z[1, 1]
    if abs(z21) < 1e-300:
        raise ValueError("Z21 = 0: the network has no forward transmission")
    return np.array([[z11 / z21, (z11 * z22 - z12 * z21) / z21],
                     [1.0 / z21, z22 / z21]], dtype=complex)


def abcd_to_z(abcd: np.ndarray) -> np.ndarray:
    abcd = np.asarray(abcd, dtype=complex)
    _check_2x2(abcd)
    a, b, c, d = abcd[0, 0], abcd[0, 1], abcd[1, 0], abcd[1, 1]
    if abs(c) < 1e-300:
        raise ValueError("C = 0: the network cannot be expressed in Z parameters")
    return np.array([[a / c, (a * d - b * c) / c],
                     [1.0 / c, d / c]], dtype=complex)


def abcd_to_s(abcd: np.ndarray, z0: float = 50.0) -> np.ndarray:
    abcd = np.asarray(abcd, dtype=complex)
    _check_2x2(abcd)
    a, b, c, d = abcd[0, 0], abcd[0, 1], abcd[1, 0], abcd[1, 1]
    denom = a + b / z0 + c * z0 + d
    return np.array([
        [(a + b / z0 - c * z0 - d) / denom, 2.0 * (a * d - b * c) / denom],
        [2.0 / denom, (-a + b / z0 - c * z0 + d) / denom]], dtype=complex)


def cascade(*networks: np.ndarray) -> np.ndarray:
    """Cascade ABCD matrices left to right. This is why ABCD exists."""
    if not networks:
        raise ValueError("nothing to cascade")
    out = np.eye(2, dtype=complex)
    for net in networks:
        out = out @ np.asarray(net, dtype=complex)
    return out


def series_impedance(z: complex) -> np.ndarray:
    """ABCD of a series element."""
    return np.array([[1.0, complex(z)], [0.0, 1.0]], dtype=complex)


def shunt_admittance(y: complex) -> np.ndarray:
    """ABCD of a shunt element."""
    return np.array([[1.0, 0.0], [complex(y), 1.0]], dtype=complex)


def transmission_line(z0_line: float, electrical_length_rad: float) -> np.ndarray:
    """ABCD of a lossless transmission line of given electrical length."""
    beta_l = float(electrical_length_rad)
    return np.array([[math.cos(beta_l), 1j * z0_line * math.sin(beta_l)],
                     [1j * math.sin(beta_l) / z0_line, math.cos(beta_l)]], dtype=complex)


def input_impedance(z_load: complex, z0_line: float,
                    electrical_length_rad: float) -> complex:
    """Impedance looking into a lossless line terminated in `z_load`.

    Zin = Z0 * (ZL + j*Z0*tan(bl)) / (Z0 + j*ZL*tan(bl)), with the
    quarter-wave and half-wave special cases handled exactly rather than
    through a tan() that blows up. An open-circuit load (infinite ZL) gives
    -j Z0 cot(bl); a load that puts the input at a pole (ZL = j Z0 / tan(bl))
    gives an open circuit, returned as infinity.
    """
    bl = float(electrical_length_rad)
    z_load = complex(z_load)
    inf = complex(float("inf"), 0.0)
    open_load = cmath.isinf(z_load)
    quarter = abs((bl % math.pi) - math.pi / 2) < 1e-12
    if quarter:
        if open_load:
            return 0j
        if z_load == 0:
            return inf
        return z0_line ** 2 / z_load
    if abs(bl % math.pi) < 1e-12:
        return z_load
    t = math.tan(bl)
    if open_load:
        return z0_line / (1j * t)
    den = z0_line + 1j * z_load * t
    if abs(den) <= 1e-12 * (z0_line + abs(z_load * t)):
        return inf
    return z0_line * (z_load + 1j * z0_line * t) / den


def _check_2x2(m: np.ndarray) -> None:
    if m.shape != (2, 2):
        raise ValueError(f"expected a 2x2 network matrix, got shape {m.shape}")
