"""Unit handling for the UI edge. The engine core is always SI."""
from __future__ import annotations

import math
from typing import Iterable

# multiplicative factor: value_in_SI = value_in_unit * FACTOR[unit]
LENGTH = {"m": 1.0, "cm": 1e-2, "mm": 1e-3, "um": 1e-6, "in": 0.0254,
          "inch": 0.0254, "mil": 2.54e-5, "ft": 0.3048}
FREQUENCY = {"Hz": 1.0, "kHz": 1e3, "MHz": 1e6, "GHz": 1e9, "THz": 1e12}
ANGLE = {"rad": 1.0, "deg": math.pi / 180.0}
POWER = {"W": 1.0, "mW": 1e-3, "kW": 1e3}
RESISTANCE = {"ohm": 1.0, "kohm": 1e3, "mohm": 1e-3}

_TABLES = (LENGTH, FREQUENCY, ANGLE, POWER, RESISTANCE)


def to_si(value: float, unit: str) -> float:
    """Convert `value` expressed in `unit` into SI."""
    for table in _TABLES:
        if unit in table:
            return value * table[unit]
    if unit in ("", None, "1", "-", "dimensionless"):
        return value
    raise KeyError(f"unknown unit {unit!r}")


def from_si(value: float, unit: str) -> float:
    """Convert an SI `value` into `unit`."""
    for table in _TABLES:
        if unit in table:
            return value / table[unit]
    if unit in ("", None, "1", "-", "dimensionless"):
        return value
    raise KeyError(f"unknown unit {unit!r}")


def convert(value: float, src: str, dst: str) -> float:
    """Convert between two units of the same kind."""
    for table in _TABLES:
        if src in table and dst in table:
            return value * table[src] / table[dst]
    raise KeyError(f"cannot convert {src!r} -> {dst!r} (different kinds or unknown)")


# ---------------------------------------------------------------- dB helpers

def db10(x: float) -> float:
    """Power ratio -> dB."""
    if x <= 0:
        return float("-inf")
    return 10.0 * math.log10(x)


def db20(x: float) -> float:
    """Field/voltage ratio -> dB."""
    if x <= 0:
        return float("-inf")
    return 20.0 * math.log10(x)


def undb10(x_db: float) -> float:
    return 10.0 ** (x_db / 10.0)


def undb20(x_db: float) -> float:
    return 10.0 ** (x_db / 20.0)


def dbi_to_linear(g_dbi: float) -> float:
    return undb10(g_dbi)


def linear_to_dbi(g: float) -> float:
    return db10(g)


def engineering(value: float, unit: str = "", sig: int = 4) -> str:
    """Format a number with an SI prefix, e.g. 2.4e9 -> '2.400 G'."""
    if value == 0 or not math.isfinite(value):
        return f"{value:g} {unit}".strip()
    prefixes = {-12: "p", -9: "n", -6: "u", -3: "m", 0: "", 3: "k", 6: "M", 9: "G", 12: "T"}
    exp = int(math.floor(math.log10(abs(value)) / 3.0) * 3)
    exp = max(-12, min(12, exp))
    scaled = value / (10.0 ** exp)
    return f"{scaled:.{sig}g} {prefixes[exp]}{unit}".strip()
