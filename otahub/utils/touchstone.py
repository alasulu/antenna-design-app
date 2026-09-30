"""Touchstone (.sNp) reading and writing.

The exporters send a model out to CST or HFSS. This is the way back: read the
S-parameters a solver or a VNA produced and put them beside what the archetype
predicted.

The details in the format that catch people out, all handled here (section
numbers and wording from the Touchstone 2.1 specification, IBIS 2024):

* **Two-port files store their matrix in the other order.** A Version 1 2-port
  line runs `freq N11 N21 N12 N22` - column-major - while three ports and up
  run row-major. Version 2 files say which with `[Two-Port Data Order]`:
  `21_12` is that same Version 1 order, `12_21` the row-major one.
* **Version 1 Z and Y data are normalised** to the option line's R; Version 2
  Z and Y data are in ohms and siemens. The same numbers mean different
  networks in the two versions.
* **A frequency point may span any number of lines.** Real files wrap, so this
  reads a flat stream of numbers and chunks it by the values each point must
  contain.
* **Two-port files may carry noise parameters after the network data**: five
  values a line (frequency, NFmin in dB, |Gamma_opt|, its angle, Rn). In
  Version 1 files the noise block begins where the frequency stops increasing;
  Version 2 files mark it `[Noise Data]`.

Version 2 keywords are honoured: `[Version]`, `[Number of Ports]`,
`[Two-Port Data Order]`, `[Reference]` (which may span lines),
`[Matrix Format]` Full/Lower/Upper, `[Network Data]`, `[Noise Data]`, `[End]`,
and `[Begin Information]` blocks are skipped. Unequal per-port reference
resistances (Version 2 `[Reference]`, or a Version 1.1 option line with one R
per port) are renormalised to a single reference, since a :class:`Network`
holds one z0; the comments say so. Mixed-mode data is refused with a message.
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
    """Sampled network parameters, always stored as S with frequency in Hz.

    `noise`, for a two-port file that carries it, is a complex (M, 4) array
    with columns frequency (Hz), NFmin (dB), Gamma_opt (about the option
    line's R, as the specification defines it) and Rn (ohm, de-normalised);
    the real columns have zero imaginary part.
    """

    frequency_hz: np.ndarray
    s: np.ndarray                      # (F, N, N) complex
    z0: float = 50.0
    comments: list[str] = field(default_factory=list)
    source: str = ""
    noise: np.ndarray | None = None

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
    filename suffix, the `[Number of Ports]` keyword or the data layout imply.
    """
    text, name = _load(source)
    label = name or "input"
    if n_ports is None:
        n_ports = _ports_from_name(name)

    comments: list[str] = []
    option = {"freq": "ghz", "param": "s", "format": "ma", "z0": [50.0]}
    version = 1
    order: str | None = None
    matrix = "full"
    reference: list[float] | None = None
    seen_option = False
    numbers: list[float] = []
    noise_numbers: list[float] = []
    target = numbers
    pending_reference = False
    in_information = False

    for raw in text.splitlines():
        data, comment = _strip(raw)
        if comment:
            comments.append(comment)
        data = data.strip()
        if not data:
            continue
        if in_information:
            if data.lower().startswith("[end information]"):
                in_information = False
            continue
        if data.startswith("["):
            pending_reference = False
            key, _, value = data.partition("]")
            key = key[1:].strip().lower()
            value = value.strip()
            if key == "version":
                version = 2 if float(value.split()[0]) >= 2.0 else 1
            elif key == "number of ports":
                keyword_ports = int(float(value))
                if n_ports is None:
                    n_ports = keyword_ports
            elif key == "two-port data order":
                order = value.lower()
                if order not in ("12_21", "21_12"):
                    raise ValueError(f"{label}: [Two-Port Data Order] must be "
                                     f"12_21 or 21_12, not {value!r}")
            elif key == "matrix format":
                matrix = value.lower() or "full"
                if matrix not in ("full", "lower", "upper"):
                    raise ValueError(f"{label}: [Matrix Format] must be Full, "
                                     f"Lower or Upper, not {value!r}")
            elif key == "reference":
                reference = [float(t) for t in value.split()]
                pending_reference = True
            elif key == "mixed-mode order":
                raise ValueError(f"{label}: mixed-mode Touchstone data is not "
                                 f"supported by this reader")
            elif key == "begin information":
                in_information = True
            elif key == "noise data":
                target = noise_numbers
            elif key == "end":
                break
            continue
        if data.startswith("#"):
            _parse_option(data[1:], option)
            seen_option = True
            continue
        values = [float(tok) for tok in data.replace(",", " ").split()]
        if pending_reference and reference is not None and (
                n_ports is None or len(reference) < n_ports):
            reference.extend(values)             # [Reference] may span lines
            continue
        pending_reference = False
        target.extend(values)

    if not numbers:
        raise ValueError(f"{label} contains no data points")

    lower_or_upper = version == 2 and matrix != "full"
    if n_ports is None:
        n_ports = _infer_ports(text, numbers, name, noise_numbers or None)
    per_point = n_ports * (n_ports + 1) // 2 if lower_or_upper else n_ports * n_ports
    stride = 1 + 2 * per_point

    if version == 1 and n_ports == 2 and not noise_numbers:
        numbers, noise_numbers = _split_noise(numbers, stride)
    if len(numbers) % stride:
        raise ValueError(
            f"{label}: {len(numbers)} numbers is not a whole number of "
            f"{n_ports}-port points ({stride} each). Wrong port count, or a "
            f"truncated file.")
    block = np.asarray(numbers, dtype=float).reshape(-1, stride)

    freq = block[:, 0] * _FREQ_UNITS[option["freq"]]
    pairs = block[:, 1:].reshape(len(block), per_point, 2)
    values = _to_complex(pairs, option["format"])

    if lower_or_upper:
        m = _from_triangle(values, n_ports, matrix)
    elif n_ports == 2 and (order or "21_12") == "21_12":
        # Version 1's order, and 21_12 in Version 2: N11 N21 N12 N22.
        m = values.reshape(-1, 2, 2).transpose(0, 2, 1)
    else:
        m = values.reshape(-1, n_ports, n_ports)

    refs = reference if (version == 2 and reference) else option["z0"]
    if len(refs) not in (1, n_ports):
        raise ValueError(f"{label}: {len(refs)} reference resistances for "
                         f"{n_ports} ports")
    z0 = float(refs[0])            # all equal, or the one they renormalise to

    param = option["param"]
    if param == "s":
        s = m
        if len(refs) > 1 and len(set(refs)) > 1:
            s = _renormalise(s, np.asarray(refs, dtype=float), z0)
            comments.append(
                f"per-port references {refs} renormalised to one reference, "
                f"{z0:g} ohm")
    else:
        if version == 1:
            # normalised data: z = Z / R (or y = Y R), each port by its own R
            r = np.asarray(refs if len(refs) == n_ports else refs * n_ports,
                           dtype=float)
            sq = np.sqrt(r)
            if param == "z":
                m = m * np.outer(sq, sq)
            elif param == "y":
                m = m / np.outer(sq, sq)
        s = _convert_to_s(m, param, z0)

    noise = None
    if noise_numbers:
        if n_ports != 2:
            raise ValueError(f"{label}: noise data is only allowed in 2-port files")
        if len(noise_numbers) % 5:
            raise ValueError(f"{label}: noise data must be five values a line; "
                             f"got {len(noise_numbers)} values")
        nb = np.asarray(noise_numbers, dtype=float).reshape(-1, 5)
        rn = nb[:, 4] * (z0 if version == 1 else 1.0)
        noise = np.column_stack([nb[:, 0] * _FREQ_UNITS[option["freq"]], nb[:, 1],
                                 nb[:, 2] * np.exp(1j * np.radians(nb[:, 3])), rn])

    if not seen_option:
        comments.append(
            "no option line found; assumed the specification's defaults "
            "(GHz, S, MA, R 50)")
    return Network(freq, s, z0, comments, str(name), noise)


def _split_noise(numbers: list[float], stride: int) -> tuple[list[float], list[float]]:
    """Version 1 two-port noise data follows the network data with no marker:
    its first frequency is at or below the last network frequency."""
    last = None
    for start in range(0, len(numbers), stride):
        f = numbers[start]
        if last is not None and f <= last:
            return numbers[:start], numbers[start:]
        last = f
    return numbers, []


def _from_triangle(values: np.ndarray, n: int, matrix: str) -> np.ndarray:
    """Expand Lower/Upper [Matrix Format] data (row by row) to full symmetric."""
    out = np.zeros((values.shape[0], n, n), dtype=complex)
    k = 0
    for i in range(n):
        cols = range(0, i + 1) if matrix == "lower" else range(i, n)
        for j in cols:
            out[:, i, j] = values[:, k]
            out[:, j, i] = values[:, k]
            k += 1
    return out


def _renormalise(s: np.ndarray, refs: np.ndarray, z0: float) -> np.ndarray:
    """S with real per-port references `refs` to S with one reference `z0`,
    through Z = sqrt(R) (I + S)(I - S)^-1 sqrt(R)."""
    n = refs.size
    eye = np.eye(n)
    sq = np.diag(np.sqrt(refs))
    out = []
    for m in s:
        z = sq @ (eye + m) @ np.linalg.inv(eye - m) @ sq
        out.append((z - z0 * eye) @ np.linalg.inv(z + z0 * eye))
    return np.stack(out)


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


def _data_widths(text: str, count: int = 2) -> list[int]:
    """Value counts of the first `count` data lines."""
    out: list[int] = []
    for raw in text.splitlines():
        data, _ = _strip(raw)
        data = data.strip()
        if data and not data.startswith(("#", "[")):
            out.append(len(data.replace(",", " ").split()))
            if len(out) == count:
                break
    return out + [0] * (count - len(out))


def _infer_ports(text: str, numbers: list[float], name: str,
                 noise: list[float] | None = None) -> int:
    """Port count, from the first line's width and the total value count.

    Neither signal is sufficient alone. One and two ports put the whole matrix
    on one line (3 and 9 values); three ports and up put one matrix ROW on each
    (2N + 1 values, or at most four pairs a line in Version 1). So a nine-value
    first line is either a 2-port or the first row of a 4-port - a real
    ambiguity in the format, and the reason the .sNp suffix exists. The total
    count breaks the tie, because usually only one candidate divides it evenly
    (a 2-port's trailing noise block counted aside).
    """
    width, second = _data_widths(text)
    total = len(numbers)
    candidates: list[int] = []
    if width == 3:
        candidates = [1]
    elif width == 9:
        candidates = [2, 4]                       # the ambiguous case
        # five ports and up wrap a row at four pairs; the continuation line
        # carries no frequency, so its width tells the port count
        candidates += [n for n in range(5, 9) if second == 2 * min(4, n - 4)]
    elif width >= 7 and width % 2 == 1:
        candidates = [(width - 1) // 2]
    if not candidates:
        candidates = list(range(1, 9))

    def fits(n: int) -> bool:
        stride = 1 + 2 * n * n
        if n == 2 and noise is None:
            net, extra = _split_noise(numbers, stride)
            return len(net) % stride == 0 and len(extra) % 5 == 0
        return total % stride == 0

    viable = [n for n in candidates if fits(n)]
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
            refs = []
            while i + 1 < len(tokens) and _is_number(tokens[i + 1]):
                i += 1
                refs.append(float(tokens[i]))
            if refs:                              # Version 1.1: one R per port
                option["z0"] = refs
        elif tok.startswith("r") and len(tok) > 1 and _is_number(tok[1:]):
            option["z0"] = [float(tok[1:])]
        else:
            raise ValueError(f"unrecognised option-line token {tokens[i]!r}")
        i += 1


def _is_number(tok: str) -> bool:
    try:
        float(tok)
    except ValueError:
        return False
    return True


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
    """Z in ohms or Y in siemens (already de-normalised) to S about z0."""
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
        freq = repr(float(f / scale))                # round-trips exactly
        if n <= 2:
            lines.append(f"{freq} " + " ".join(cells))
        else:
            # Version 1: each matrix row starts a line, at most four pairs a line
            first = True
            for row in range(n):
                cols = cells[row * n:(row + 1) * n]
                for k in range(0, n, 4):
                    lead = f"{freq} " if first else "  "
                    lines.append(lead + " ".join(cols[k:k + 4]))
                    first = False
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
