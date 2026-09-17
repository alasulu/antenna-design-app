"""Touchstone (.sNp) reading and writing.

The exporters send a model out to CST or HFSS. This is the way back: read the
S-parameters a solver or a VNA produced and put them beside what the archetype
predicted.

Two details in the format catch people out, and both are handled here:

* **Two-port files store their matrix in the other order.** A 2-port line runs
  `freq S11 S21 S12 S22` - column-major - while three ports and up run
  row-major, one matrix row per line. That single exception is the oldest trap
  in the format, and getting it wrong silently transposes every two-port file.
* **A frequency point may span any number of lines.** The specification lets
  data wrap, and real files from real instruments do. Parsing line-by-line
  works until it doesn't, so this reads a flat stream of numbers and chunks it
  by the 1 + 2N^2 values each point must contain.

Version 1.1 keyword blocks (`[Version]`, `[Number of Ports]`, and the rest) are
recognised and their port count honoured; the older option-line form is the
default and needs no keywords at all.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .network import s_to_z, z_to_y

__all__ = ["Network", "read_touchstone", "write_touchstone", "compare_to_prediction"]

_FREQ_UNITS = {"hz": 1.0, "khz": 1e3, "mhz": 1e6, "ghz": 1e9, "thz": 1e12}
_FORMATS = ("ma", "db", "ri")
_PARAMS = ("s", "y", "z", "g", "h")


@dataclass(slots=True)
class Network:
    """Sampled network parameters, always stored as S with frequency in Hz."""

    frequency_hz: np.ndarray
    s: np.ndarray                      # (F, N, N) complex
    z0: float = 50.0
    comments: list[str] = field(default_factory=list)
    source: str = ""

    def __post_init__(self) -> None:
        self.frequency_hz = np.asarray(self.frequency_hz, dtype=float)
        self.s = np.asarray(self.s, dtype=complex)
        if self.s.ndim != 3 or self.s.shape[1] != self.s.shape[2]:
            raise ValueError(f"S must be (F, N, N), got {self.s.shape}")
        if self.s.shape[0] != self.frequency_hz.size:
            raise ValueError(
                f"{self.frequency_hz.size} frequencies but {self.s.shape[0]} "
                f"S-matrices")
        if np.any(np.diff(self.frequency_hz) <= 0):
            raise ValueError("frequencies must increase strictly")

    @property
    def n_ports(self) -> int:
        return int(self.s.shape[1])

    @property
    def z(self) -> np.ndarray:
        """Impedance matrices. One-port files give the load impedance."""
        return np.stack([s_to_z(m, self.z0) if m.shape[0] == 2
                         else _s_to_z_general(m, self.z0) for m in self.s])

    @property
    def y(self) -> np.ndarray:
        return np.stack([z_to_y(m) if m.shape[0] == 2 else np.linalg.inv(m)
                         for m in self.z])

    def s_db(self, i: int = 0, j: int | None = None) -> np.ndarray:
        """|S_ij| in dB. With one index, the reflection at that port."""
        j = i if j is None else j
        return 20.0 * np.log10(np.maximum(np.abs(self.s[:, i, j]), 1e-300))

    def vswr(self, port: int = 0) -> np.ndarray:
        g = np.abs(self.s[:, port, port])
        return (1 + g) / np.maximum(1 - g, 1e-300)

    def impedance_at_port(self, port: int = 0) -> np.ndarray:
        """Port impedance seen looking in, from its own reflection alone.

        For a multi-port this ignores what the other ports are doing, which is
        the right reading only when they are terminated in z0 - the usual case
        for a measurement, and worth stating rather than assuming.
        """
        g = self.s[:, port, port]
        return self.z0 * (1 + g) / (1 - g)

    def at(self, frequency_hz: float) -> np.ndarray:
        """S at one frequency, linearly interpolated between samples."""
        f = float(frequency_hz)
        lo, hi = self.frequency_hz[0], self.frequency_hz[-1]
        if not (lo <= f <= hi):
            raise ValueError(
                f"{f/1e9:.6g} GHz is outside the measured span "
                f"{lo/1e9:.6g}-{hi/1e9:.6g} GHz; this does not extrapolate")
        real = np.stack([np.interp(f, self.frequency_hz, self.s[:, i, j].real)
                         for i in range(self.n_ports) for j in range(self.n_ports)])
        imag = np.stack([np.interp(f, self.frequency_hz, self.s[:, i, j].imag)
                         for i in range(self.n_ports) for j in range(self.n_ports)])
        n = self.n_ports
        return (real + 1j * imag).reshape(n, n)

    def resonances(self, port: int = 0, threshold_db: float = -10.0) -> list[float]:
        """Frequencies where the reflection dips to a local minimum below
        `threshold_db`. What you actually look for in a measured antenna."""
        db = self.s_db(port)
        out = []
        for i in range(1, len(db) - 1):
            if db[i] < threshold_db and db[i] <= db[i - 1] and db[i] <= db[i + 1]:
                out.append(float(self.frequency_hz[i]))
        return out


def _s_to_z_general(s: np.ndarray, z0: float) -> np.ndarray:
    """Z = z0 (I + S)(I - S)^-1, for any port count."""
    n = s.shape[0]
    eye = np.eye(n)
    return z0 * (eye + s) @ np.linalg.inv(eye - s)


# --------------------------------------------------------------------- reading

def _strip(line: str) -> tuple[str, str]:
    """Split a line into (data, comment) on the first '!'."""
    idx = line.find("!")
    return (line, "") if idx < 0 else (line[:idx], line[idx + 1:].strip())


def read_touchstone(source, n_ports: int | None = None) -> Network:
    """Read a Touchstone file, or its text, into a :class:`Network`.

    `source` is a path or the file contents. `n_ports` overrides what the
    filename suffix or the data layout implies.
    """
    text, name = _load(source)
    if n_ports is None:
        n_ports = _ports_from_name(name)

    comments: list[str] = []
    option = {"freq": "ghz", "param": "s", "format": "ma", "z0": 50.0}
    seen_option = False
    numbers: list[float] = []

    for raw in text.splitlines():
        data, comment = _strip(raw)
        if comment:
            comments.append(comment)
        data = data.strip()
        if not data:
            continue
        if data.startswith("["):
            key, _, value = data.partition("]")
            key = key[1:].strip().lower()
            value = value.strip()
            if key == "number of ports":
                n_ports = int(float(value))
            elif key == "two-port data order" and value:
                option["order"] = value.strip().lower()
            elif key in ("network data", "end"):
                continue
            continue
        if data.startswith("#"):
            _parse_option(data[1:], option)
            seen_option = True
            continue
        numbers.extend(float(tok) for tok in data.replace(",", " ").split())

    if not numbers:
        raise ValueError(f"{name or 'input'} contains no data points")
    if n_ports is None:
        n_ports = _infer_ports(text, len(numbers), name)

    stride = 1 + 2 * n_ports * n_ports
    if len(numbers) % stride:
        raise ValueError(
            f"{name or 'input'}: {len(numbers)} numbers is not a whole number of "
            f"{n_ports}-port points ({stride} each). Wrong port count, or a "
            f"truncated file.")
    block = np.asarray(numbers, dtype=float).reshape(-1, stride)

    freq = block[:, 0] * _FREQ_UNITS[option["freq"]]
    pairs = block[:, 1:].reshape(len(block), n_ports * n_ports, 2)
    values = _to_complex(pairs, option["format"])

    if n_ports == 2 and option.get("order", "12_21") != "21_12":
        # The historical exception: 2-port files are column-major.
        s = values.reshape(-1, 2, 2).transpose(0, 2, 1)
    else:
        s = values.reshape(-1, n_ports, n_ports)

    if option["param"] != "s":
        s = _convert_to_s(s, option["param"], option["z0"])

    if not seen_option:
        comments.append(
            "no option line found; assumed the specification's defaults "
            "(GHz, S, MA, R 50)")
    return Network(freq, s, option["z0"], comments, str(name))


def _load(source) -> tuple[str, str]:
    if isinstance(source, Path):
        return source.read_text(), source.name
    if isinstance(source, str) and "\n" not in source and len(source) < 4096:
        path = Path(source)
        if path.exists():
            return path.read_text(), path.name
    if isinstance(source, str):
        return source, ""
    raise TypeError(f"cannot read Touchstone data from {type(source).__name__}")


def _ports_from_name(name: str) -> int | None:
    m = re.search(r"\.s(\d+)p$", name, re.IGNORECASE)
    return int(m.group(1)) if m else None


def _first_data_width(text: str) -> int:
    for raw in text.splitlines():
        data, _ = _strip(raw)
        data = data.strip()
        if data and not data.startswith(("#", "[")):
            return len(data.replace(",", " ").split())
    return 0


def _infer_ports(text: str, total: int, name: str) -> int:
    """Port count, from the first line's width and the total value count.

    Neither signal is sufficient alone. One and two ports put the whole matrix
    on one line (3 and 9 values); three ports and up put one matrix ROW on each
    (2N + 1 values). So a nine-value first line is either a 2-port or the first
    row of a 4-port - a real ambiguity in the format, and the reason the .sNp
    suffix exists. The total count breaks the tie, because only one candidate
    divides it evenly.
    """
    width = _first_data_width(text)
    candidates: list[int] = []
    if width == 3:
        candidates = [1]
    elif width == 9:
        candidates = [2, 4]                       # the ambiguous case
    elif width >= 7 and width % 2 == 1:
        candidates = [(width - 1) // 2]
    if not candidates:
        candidates = list(range(1, 9))
    viable = [n for n in candidates if total % (1 + 2 * n * n) == 0]
    if len(viable) == 1:
        return viable[0]
    if not viable:
        raise ValueError(
            f"{name or 'input'}: cannot infer the port count from {total} "
            f"values with a {width}-value first line. Pass n_ports, or name "
            f"the file .sNp.")
    raise ValueError(
        f"{name or 'input'}: {total} values is consistent with "
        f"{' and '.join(str(n) for n in viable)} ports. Pass n_ports, or name "
        f"the file .sNp - this is why the suffix exists.")


def _parse_option(text: str, option: dict) -> None:
    tokens = text.replace(",", " ").split()
    i = 0
    while i < len(tokens):
        tok = tokens[i].lower()
        if tok in _FREQ_UNITS:
            option["freq"] = tok
        elif tok in _PARAMS:
            option["param"] = tok
        elif tok in _FORMATS:
            option["format"] = tok
        elif tok == "r":
            i += 1
            if i < len(tokens):
                option["z0"] = float(tokens[i])
        elif tok.startswith("r") and len(tok) > 1:
            option["z0"] = float(tok[1:])
        else:
            raise ValueError(f"unrecognised option-line token {tokens[i]!r}")
        i += 1


def _to_complex(pairs: np.ndarray, fmt: str) -> np.ndarray:
    a, b = pairs[..., 0], pairs[..., 1]
    if fmt == "ri":
        return a + 1j * b
    if fmt == "ma":
        return a * np.exp(1j * np.radians(b))
    if fmt == "db":
        return 10.0 ** (a / 20.0) * np.exp(1j * np.radians(b))
    raise ValueError(f"unknown data format {fmt!r}")


def _convert_to_s(m: np.ndarray, param: str, z0: float) -> np.ndarray:
    if param == "z":
        eye = np.eye(m.shape[1])
        return np.stack([(x / z0 - eye) @ np.linalg.inv(x / z0 + eye) for x in m])
    if param == "y":
        eye = np.eye(m.shape[1])
        return np.stack([(eye - x * z0) @ np.linalg.inv(eye + x * z0) for x in m])
    raise ValueError(
        f"{param.upper()}-parameter files are recognised but not converted; "
        f"G and H parameters need a two-port-specific transformation this "
        f"module does not implement")


# --------------------------------------------------------------------- writing

def write_touchstone(network: Network, path=None, fmt: str = "ri",
                     freq_unit: str = "ghz") -> str:
    """Render a network as Touchstone text, optionally writing it to `path`."""
    fmt, freq_unit = fmt.lower(), freq_unit.lower()
    if fmt not in _FORMATS:
        raise ValueError(f"format must be one of {_FORMATS}, not {fmt!r}")
    if freq_unit not in _FREQ_UNITS:
        raise ValueError(f"unknown frequency unit {freq_unit!r}")
    n = network.n_ports
    lines = [f"! written by OTA Hub Antenna Toolkit",
             f"! {n}-port, {len(network.frequency_hz)} points"]
    lines += [f"! {c}" for c in network.comments]
    lines.append(f"# {freq_unit.upper()} S {fmt.upper()} R {network.z0:g}")

    scale = _FREQ_UNITS[freq_unit]
    for f, mat in zip(network.frequency_hz, network.s):
        ordered = mat.T if n == 2 else mat          # the 2-port exception, again
        flat = ordered.reshape(-1)
        cells = []
        for value in flat:
            if fmt == "ri":
                cells.append(f"{value.real:< .9g} {value.imag:< .9g}")
            elif fmt == "ma":
                cells.append(f"{abs(value):< .9g} {math.degrees(np.angle(value)):< .9g}")
            else:
                mag = 20 * math.log10(max(abs(value), 1e-300))
                cells.append(f"{mag:< .9g} {math.degrees(np.angle(value)):< .9g}")
        if n <= 2:
            lines.append(f"{f/scale:< .9g} " + " ".join(cells))
        else:
            lines.append(f"{f/scale:< .9g} " + " ".join(cells[:n]))
            for row in range(1, n):
                lines.append("  " + " ".join(cells[row * n:(row + 1) * n]))
    text = "\n".join(lines) + "\n"
    if path is not None:
        Path(path).write_text(text)
    return text


# ----------------------------------------------------------------- comparison

def compare_to_prediction(network: Network, predicted_impedance: complex,
                          frequency_hz: float, port: int = 0) -> dict:
    """Measured against predicted, at one frequency.

    This is the point of reading Touchstone at all: an archetype predicts an
    input impedance, a solver or a VNA produces a file, and the two need to be
    put side by side in the same terms.
    """
    s = network.at(frequency_hz)[port, port]
    z_meas = network.z0 * (1 + s) / (1 - s)
    z_pred = complex(predicted_impedance)
    g_pred = (z_pred - network.z0) / (z_pred + network.z0)

    def _vswr(g):
        m = abs(g)
        return (1 + m) / (1 - m) if m < 1 else float("inf")

    return {
        "frequency_hz": float(frequency_hz),
        "measured_impedance_ohm": complex(z_meas),
        "predicted_impedance_ohm": z_pred,
        "measured_s11_db": 20 * math.log10(max(abs(s), 1e-300)),
        "predicted_s11_db": 20 * math.log10(max(abs(g_pred), 1e-300)),
        "measured_vswr": _vswr(s),
        "predicted_vswr": _vswr(g_pred),
        "resistance_error_pct": (
            (z_meas.real - z_pred.real) / z_pred.real * 100
            if z_pred.real else float("nan")),
        "reactance_error_ohm": float(z_meas.imag - z_pred.imag),
        "impedance_error_ohm": float(abs(z_meas - z_pred)),
    }
