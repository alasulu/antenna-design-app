"""Impedance matching networks.

L-sections, quarter-wave transformers and single-stub tuners. Every routine
returns *all* the solutions the topology admits, because the choice between
them is an engineering one - bandwidth, component values, board area - that
the tool should not make silently.
"""
from __future__ import annotations

import cmath
import math
from dataclasses import dataclass, field

from ..core.constants import C0
from .network import input_impedance, match_quality, reflection_coefficient


@dataclass(frozen=True, slots=True)
class LumpedElement:
    """A reactance realised as an inductor or capacitor at one frequency."""
    kind: str            # "L" or "C"
    value: float         # henries or farads
    reactance: float     # ohm at the design frequency

    def __str__(self) -> str:
        if self.kind == "L":
            return f"L = {self.value * 1e9:.4g} nH  (jX = {self.reactance:+.2f} ohm)"
        return f"C = {self.value * 1e12:.4g} pF  (jX = {self.reactance:+.2f} ohm)"


def element_from_reactance(x: float, f_hz: float) -> LumpedElement:
    """Positive reactance is an inductor, negative a capacitor."""
    omega = 2.0 * math.pi * f_hz
    if x >= 0:
        return LumpedElement("L", x / omega, x)
    return LumpedElement("C", -1.0 / (omega * x), x)


def element_from_susceptance(b: float, f_hz: float) -> LumpedElement:
    """Positive susceptance is a capacitor, negative an inductor."""
    omega = 2.0 * math.pi * f_hz
    if b >= 0:
        return LumpedElement("C", b / omega, -1.0 / b if b else math.inf)
    return LumpedElement("L", -1.0 / (omega * b), -1.0 / b)


@dataclass(slots=True)
class LSection:
    """One L-network solution."""
    topology: str                 # description of element order
    series: LumpedElement
    shunt: LumpedElement
    f_hz: float
    z_load: complex
    z0: float
    achieved: complex = field(default=0j)

    def __str__(self) -> str:
        return (f"{self.topology}\n    series: {self.series}\n"
                f"    shunt : {self.shunt}")


def l_section(z_load: complex, z0: float = 50.0, f_hz: float = 1e9) -> list[LSection]:
    """Match `z_load` to `z0` with a two-element L network.

    Returns both solutions. Which topology applies depends on whether the load
    lies inside the 1 + jx circle on the Smith chart:

    - RL > Z0: shunt element across the LOAD, then series to the source.
    - RL < Z0: series element at the load, then shunt to the source.

    Verified against Pozar Example 5.1 (200 - j100 to 100 ohm at 500 MHz),
    which yields C = 0.92 pF with L = 38.8 nH, and C = 2.61 pF with L = 46.1 nH.
    """
    z_load = complex(z_load)
    rl, xl = z_load.real, z_load.imag
    if rl <= 0:
        raise ValueError(f"load resistance must be positive, got {rl}")
    out: list[LSection] = []

    if rl >= z0:
        # Shunt across the load, series toward the source.
        disc = rl * rl + xl * xl - z0 * rl
        if disc < 0:
            return []
        # Pozar 5.3a: sqrt(RL/Z0), not its inverse. Getting this upside down
        # still produces plausible component values that do not match.
        root = math.sqrt(rl / z0) * math.sqrt(disc)
        for sign in (+1.0, -1.0):
            b = (xl + sign * root) / (rl * rl + xl * xl)
            if b == 0:
                continue
            x = 1.0 / b + xl * z0 / rl - z0 / (b * rl)
            section = LSection(
                "shunt at load, then series (RL > Z0)",
                element_from_reactance(x, f_hz),
                element_from_susceptance(b, f_hz),
                f_hz, z_load, z0)
            section.achieved = _verify_shunt_first(z_load, b, x)
            out.append(section)
    else:
        # Series at the load, shunt toward the source.
        disc = rl * (z0 - rl)
        if disc < 0:
            return []
        for sign in (+1.0, -1.0):
            # X and B take the SAME sign here. Derived rather than recalled:
            # with X = +sqrt(RL(Z0-RL)) - XL, the residual susceptance that the
            # shunt must cancel works out to +sqrt((Z0-RL)/RL)/Z0.
            x = sign * math.sqrt(disc) - xl
            b = sign * math.sqrt((z0 - rl) / rl) / z0
            if b == 0:
                continue
            section = LSection(
                "series at load, then shunt (RL < Z0)",
                element_from_reactance(x, f_hz),
                element_from_susceptance(b, f_hz),
                f_hz, z_load, z0)
            section.achieved = _verify_series_first(z_load, x, b)
            out.append(section)
    return out


def _verify_shunt_first(z_load: complex, b: float, x: float) -> complex:
    """Forward-evaluate the network so the caller can check it really matches."""
    y = 1.0 / z_load + 1j * b
    return 1.0 / y + 1j * x


def _verify_series_first(z_load: complex, x: float, b: float) -> complex:
    z = z_load + 1j * x
    return 1.0 / (1.0 / z + 1j * b)


def quarter_wave(z_load: float, z0: float = 50.0) -> float:
    """Transformer impedance for a REAL load: Z1 = sqrt(Z0 * ZL).

    A complex load must first be rotated to a real impedance along a line, or
    resonated out; applying this to a reactive load is a common and quiet
    mistake, so it is refused here.
    """
    if isinstance(z_load, complex) and abs(z_load.imag) > 1e-12:
        raise ValueError(
            f"quarter-wave transformer needs a real load; got {z_load}. "
            f"Rotate the load to a real impedance along a line first, or "
            f"resonate the reactance out.")
    zl = float(getattr(z_load, "real", z_load))
    if zl <= 0 or z0 <= 0:
        raise ValueError("impedances must be positive")
    return math.sqrt(z0 * zl)


def quarter_wave_bandwidth(z_load: float, z0: float = 50.0,
                           vswr_max: float = 2.0) -> float:
    """Fractional bandwidth of a single-section quarter-wave transformer.

    The classic result: bandwidth shrinks as the transformation ratio grows,
    which is why wideband transformers are multi-section.
    """
    zl = float(getattr(z_load, "real", z_load))
    gamma_m = (vswr_max - 1.0) / (vswr_max + 1.0)
    if zl == z0:
        return 2.0
    arg = (gamma_m / math.sqrt(1.0 - gamma_m ** 2)) * 2.0 * math.sqrt(z0 * zl) / abs(zl - z0)
    if arg >= 1.0:
        return 2.0
    return 2.0 - (4.0 / math.pi) * math.acos(arg)


@dataclass(slots=True)
class StubSolution:
    """A single-stub tuner solution, in wavelengths."""
    line_length_lambda: float
    stub_length_lambda: float
    stub_kind: str            # "short" or "open"
    stub_placement: str       # "shunt"
    achieved: complex

    def __str__(self) -> str:
        return (f"line d = {self.line_length_lambda:.4f} lambda, "
                f"{self.stub_kind}-circuit shunt stub "
                f"l = {self.stub_length_lambda:.4f} lambda")


def single_stub(z_load: complex, z0: float = 50.0, stub: str = "short",
                samples: int = 20001) -> list[StubSolution]:
    """Shunt single-stub tuner: find every (d, l) pair that matches.

    Solved by scanning the line length for the points where the input
    conductance equals Y0, then choosing the stub length that cancels the
    residual susceptance. A scan rather than the closed form because it is
    transparent, handles every load without case analysis, and the residual
    match is verified numerically before a solution is returned.
    """
    if stub not in ("short", "open"):
        raise ValueError(f"stub must be 'short' or 'open', got {stub!r}")
    z_load = complex(z_load)
    y0 = 1.0 / z0
    solutions: list[StubSolution] = []

    previous = None
    for i in range(samples):
        d = 0.5 * i / (samples - 1)                     # 0 .. 0.5 lambda
        y_in = 1.0 / input_impedance(z_load, z0, 2.0 * math.pi * d)
        residual = y_in.real - y0
        if previous is not None and previous[1] * residual <= 0 and i > 0:
            # Linear interpolation onto the conductance crossing.
            d0, r0 = previous
            d_star = d0 - r0 * (d - d0) / (residual - r0) if residual != r0 else d
            y_star = 1.0 / input_impedance(z_load, z0, 2.0 * math.pi * d_star)
            b_needed = -y_star.imag
            length = _stub_length(b_needed, z0, stub)
            if length is None:
                previous = (d, residual)
                continue
            y_stub = _stub_admittance(length, z0, stub)
            achieved = 1.0 / (y_star + y_stub)
            if abs(reflection_coefficient(achieved, z0)) < 1e-3:
                solutions.append(StubSolution(d_star, length, stub, "shunt", achieved))
        previous = (d, residual)
    return solutions


def _stub_length(b_needed: float, z0: float, kind: str) -> float | None:
    """Stub length in wavelengths giving the required susceptance."""
    y0 = 1.0 / z0
    if kind == "short":
        # Y = -j*Y0*cot(beta*l)  ->  cot(beta l) = -b/Y0
        target = -b_needed / y0
        beta_l = math.atan2(1.0, target)          # cot(x) = target
    else:
        # Y = j*Y0*tan(beta*l)
        target = b_needed / y0
        beta_l = math.atan(target)
    if beta_l < 0:
        beta_l += math.pi
    return beta_l / (2.0 * math.pi)


def _stub_admittance(length_lambda: float, z0: float, kind: str) -> complex:
    beta_l = 2.0 * math.pi * length_lambda
    y0 = 1.0 / z0
    t = math.tan(beta_l)
    if kind == "short":
        if abs(t) < 1e-12:
            return complex(0.0, -math.inf)
        return complex(0.0, -y0 / t)
    return complex(0.0, y0 * t)


def summarise_match(z_load: complex, z0: float = 50.0, f_hz: float = 1e9) -> dict:
    """Every matching option for a load, with the resulting quality."""
    before = match_quality(z_load, z0)
    out = {"before": before, "l_sections": l_section(z_load, z0, f_hz),
           "stubs_short": single_stub(z_load, z0, "short"),
           "stubs_open": single_stub(z_load, z0, "open")}
    try:
        out["quarter_wave"] = quarter_wave(z_load, z0)
    except ValueError as exc:
        out["quarter_wave"] = str(exc)
    return out
