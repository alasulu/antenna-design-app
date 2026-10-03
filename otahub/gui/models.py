"""Qt-side data plumbing for the catalogue.

Kept free of widget code so it can be exercised headlessly. A
QStandardItemModel behind a QSortFilterProxyModel is deliberate: a bespoke
QAbstractItemModel would be more "correct" and considerably more bug-prone for
a two-level tree that never changes shape.
"""
from __future__ import annotations

from PySide6.QtCore import QSortFilterProxyModel, Qt
from PySide6.QtGui import QStandardItem, QStandardItemModel

from ..core.registry import Registry

#: Custom role carrying the archetype key on leaf rows.
KEY_ROLE = int(Qt.ItemDataRole.UserRole) + 1


def build_catalogue_model(registry: Registry) -> QStandardItemModel:
    """Two-level tree: family -> archetype, with the key on the leaf."""
    model = QStandardItemModel()
    model.setHorizontalHeaderLabels(["Archetype"])
    for family in registry.families:
        members = registry.by_family(family)
        parent = QStandardItem(f"{family}  ({len(members)})")
        parent.setEditable(False)
        parent.setSelectable(False)
        for archetype in members:
            label = archetype.name
            if archetype.spec.confidence == "low":
                label += "   ⚠ low confidence"
            leaf = QStandardItem(label)
            leaf.setEditable(False)
            leaf.setData(archetype.key, KEY_ROLE)
            leaf.setToolTip(archetype.spec.summary or archetype.key)
            parent.appendRow(leaf)
        model.appendRow(parent)
    return model


class CatalogueFilter(QSortFilterProxyModel):
    """Free-text filter that keeps a family visible when any child matches."""

    def __init__(self, registry: Registry, parent=None) -> None:
        super().__init__(parent)
        self._registry = registry
        self._needle = ""
        self.setRecursiveFilteringEnabled(True)

    def set_needle(self, text: str) -> None:
        self._needle = text.strip().lower()
        # invalidateFilter/invalidateRowsFilter are both deprecated in Qt 6.9+;
        # invalidate() re-evaluates the whole proxy, which is trivial here.
        self.invalidate()

    def filterAcceptsRow(self, row: int, parent) -> bool:  # noqa: N802 - Qt API
        if not self._needle:
            return True
        index = self.sourceModel().index(row, 0, parent)
        key = index.data(KEY_ROLE)
        if key is None:
            # A family row never matches on its own. Returning False here is
            # what makes the filter work: with recursive filtering enabled Qt
            # keeps a parent whose descendants match, so families reappear
            # exactly when they contain a hit. Delegating to the base class
            # instead accepts every family, since the base regex is empty.
            return False
        archetype = self._registry.get(key)
        if archetype is None:
            return False
        return matches(archetype, self._needle)


def matches(archetype, needle: str) -> bool:
    """Free-text match on an archetype's key, name, family and summary."""
    needle = needle.strip().lower()
    if not needle:
        return True
    haystack = (f"{archetype.key} {archetype.name} {archetype.family} "
                f"{archetype.spec.summary}").lower()
    return needle in haystack


def free_geometry(archetype) -> set[str]:
    """Geometry the synthesis reads but never produces - a loop's turn count and wire
    radius. The design cannot finish without them, so they are inputs, as the CLI's
    "to complete this design, supply" already treats them."""
    from ..core.expr import ExprError, referenced_symbols
    spec = archetype.spec
    produced = {r.output for r in spec.synthesis}
    read: set[str] = set()
    for expr in [r.expr for r in spec.synthesis] + [r.expr for r in spec.analysis]:
        try:
            read |= referenced_symbols(expr, spec.symbols)
        except ExprError:
            continue
    return {p.symbol for p in spec.parameters
            if p.role == "geometry" and p.symbol not in produced and p.symbol in read}


def requirement_fields(archetype) -> list:
    """Parameters the user must or may supply, ordered for a form.

    Requirements first because they have no defaults and the design cannot
    proceed without them; then free geometry (see :func:`free_geometry`), which it
    cannot finish without either; assumptions and materials after, since those
    carry defensible defaults and most users will leave them alone.
    """
    free = free_geometry(archetype)
    order = {"requirement": 0, "geometry": 1, "assumption": 2, "material": 3}
    fields = [p for p in archetype.spec.parameters
              if p.role in order and (p.role != "geometry" or p.symbol in free)]
    return sorted(fields, key=lambda p: (order[p.role], p.unit != "Hz", p.symbol))


def suggested_value(param, values: dict) -> float | None:
    """The middle of a parameter's typical range, evaluated on the design's own
    values - "0.001*s .. 0.05*s" for a square loop's wire radius, the geometric
    mean of the two ends (the ranges span decades). None when the typical text
    gives no number here. A form suggestion only: the engine itself never
    invents an input."""
    import math
    import re
    from ..core.expr import ExprError, evaluate
    text = (param.typical or "").strip()
    if not text:
        return None
    ends = []
    for part in text.split(".."):
        expr = re.sub(r"\blambda\b", "lambda0", part.strip())
        try:
            v = float(evaluate(expr, values))
        except (ExprError, TypeError, ValueError, ZeroDivisionError):
            return None
        if not math.isfinite(v) or v <= 0:
            return None
        ends.append(v)
    return math.sqrt(ends[0] * ends[-1])


def is_primary(param) -> bool:
    """On the form's first page: the requirements and the free geometry."""
    return param.role in ("requirement", "geometry")


def default_frequency(archetype) -> float:
    """A sensible in-band design frequency for pre-filling the form.

    The geometric mean of the validity band, which lands mid-decade on a log
    axis rather than being dragged to the top end by an arithmetic mean over a
    band spanning several decades. Open-ended bands are clamped so the default
    stays somewhere an engineer would actually work.
    """
    import math

    lo, hi = archetype.spec.freq_range_hz
    lo = max(lo, 1e5)
    hi = min(hi if math.isfinite(hi) else 1e11, 1e11)
    if hi <= lo:
        return lo
    return math.sqrt(lo * hi)


def default_for(param, archetype=None) -> str:
    """Pre-fill text for a parameter's form field.

    Frequency gets special handling: specs rarely declare a `typical` for a
    design frequency (there is no universal answer), but an empty frequency
    field means the form opens unable to synthesise anything, which reads as
    broken.

    Matched on the declared UNIT, not on the name. Keying it to "f0" alone
    missed every archetype whose requirement is `f_low` instead - the whole
    wideband family, which opened blank and produced no geometry.
    """
    if (param.unit == "Hz" and param.role == "requirement"
            and archetype is not None):
        try:
            float(param.typical)
        except (TypeError, ValueError):
            return f"{default_frequency(archetype):.6g}"
    try:
        float(param.typical)
    except (TypeError, ValueError):
        return ""
    return str(param.typical)


# ---------------------------------------------------------------- values at the edge of the form
import math as _math
import re as _re

_UNIT_TABLES = (
    {"m": 1.0, "cm": 1e-2, "mm": 1e-3, "um": 1e-6, "µm": 1e-6, "in": 0.0254, "inch": 0.0254, "mil": 2.54e-5, "ft": 0.3048},
    {"Hz": 1.0, "kHz": 1e3, "MHz": 1e6, "GHz": 1e9, "THz": 1e12},
    {"rad": 1.0, "deg": _math.pi / 180.0, "°": _math.pi / 180.0},
    {"W": 1.0, "mW": 1e-3, "kW": 1e3},
    {"ohm": 1.0, "kohm": 1e3, "mohm": 1e-3, "Ω": 1.0, "kΩ": 1e3},
)
_PREFIX = {"p": 1e-12, "n": 1e-9, "u": 1e-6, "µ": 1e-6, "m": 1e-3, "k": 1e3, "M": 1e6, "G": 1e9, "T": 1e12}
_NUMBER = _re.compile(r"^([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)\s*(\S*)$")
ENGINEERING_UNITS = {"m", "Hz", "H", "F", "S"}
UNIT_SHOWN = {"ohm": "Ω", "deg": "°", "m^2": "m²", "mm^3": "mm³", "ohm^2": "Ω²", "-": ""}


def parse_quantity(text: str, unit: str = "") -> float:
    """A form entry in SI, or with a unit: '2.4 GHz', '1.6 mm', '50 ohm', '2.4G'.

    Raises ValueError for anything else, rather than guessing."""
    m = _NUMBER.match(text.strip())
    if not m:
        raise ValueError(f"{text!r} is not a number")
    value, suffix = float(m.group(1)), m.group(2)
    if not suffix or suffix == unit or (suffix == "°" and unit == "deg"):
        return value
    for table in _UNIT_TABLES:
        if suffix in table and unit in table:
            return value * table[suffix] / table[unit]
    if suffix[0] in _PREFIX and suffix[1:] in (unit, ""):
        return value * _PREFIX[suffix[0]]
    raise ValueError(f"{text!r} is not a number in {unit or 'these units'}")


def _g(value: float, sig: int) -> str:
    return f"{value:.{sig}g}"


def _eng(value: float, unit: str, sig: int) -> str:
    """engineering() with a proper micro sign."""
    from ..core.units import engineering
    text = engineering(value, unit, sig)
    return _re.sub(r" u(?=\S)", " µ", text)


def format_input(value: float, unit: str = "") -> str:
    """What the form shows for a value: engineering notation for lengths and
    frequencies, which parse_quantity reads back without loss at 6 figures."""
    if unit in ENGINEERING_UNITS and value and _math.isfinite(value):
        return _eng(value, unit, 6)
    return _g(value, 6)


def shown_unit(unit: str) -> str:
    return UNIT_SHOWN.get(unit, unit or "")


def display_value(key: str, value, unit: str = "", sig: int = 5) -> str | None:
    """A computed value for a table or tile; None when it is NaN (outside the
    range a fit was solved over), so the caller can say so instead of 'nan'."""
    if isinstance(value, bool):
        value = float(value)
    if isinstance(value, complex):
        if _math.isnan(value.real) or _math.isnan(value.imag):
            return None
        sign = "+" if value.imag >= 0 else "−"
        u = shown_unit(unit)
        return f"{value.real:.{min(sig, 4)}g} {sign} j{abs(value.imag):.{min(sig, 4)}g}{(' ' + u) if u else ''}"
    if not isinstance(value, (int, float)):
        return str(value)
    value = float(value)
    if _math.isnan(value):
        return None
    if _math.isinf(value):                  # a perfect null in dB, say: shown as such
        u = shown_unit(unit)
        return ("−∞" if value < 0 else "∞") + ((" " + u) if u and u != "°" else u)
    if ("fractional_bandwidth" in key or "bandwidth_vswr2" in key) and unit in ("-", "") and abs(value) < 5:
        return f"{value * 100:.3g} %"
    if key.endswith("bandwidth_ratio") and unit in ("-", "") and _math.isfinite(value):
        return f"{value:.3g} : 1"
    if unit in ENGINEERING_UNITS and value and _math.isfinite(value):
        return _eng(value, unit, min(sig, 4))
    u = shown_unit(unit)
    text = _g(value, sig)
    if not u:
        return text
    return f"{text}{u}" if u == "°" else f"{text} {u}"


_WORDS = {"dbi": "", "db": "", "dbd": "", "hz": "", "ohm": "", "deg": "", "m": "", "m2": "", "pct": "", "s": "",
          "vswr2": "VSWR 2", "vswr": "VSWR", "hpbw": "HPBW", "fnbw": "FNBW", "q": "Q", "te10": "TE10",
          "tm11": "TM11", "tm10": "TM10", "te11": "TE11", "sll": "sidelobe", "ar": "axial ratio", "f0": "f0",
          "fb": "front-to-back", "xpol": "cross-pol"}


def humanize(key: str) -> str:
    """'fractional_bandwidth_vswr2' -> 'Fractional bandwidth VSWR 2'; unit tokens drop off the end."""
    parts = key.split("_")
    while len(parts) > 1 and _WORDS.get(parts[-1]) == "":
        parts.pop()
    words = [(_WORDS[w] if _WORDS.get(w) else w) for w in parts]
    text = " ".join(words)
    return text[:1].upper() + text[1:]


#: What the headline tiles show, in order: (label, patterns tried in turn).
KEY_FIGURES = (
    ("Gain", (r"^gain_dbi$", r"^realised_gain_dbi$", r"gain_dbi$", r"^directivity_dbi$", r"directivity_dbi$")),
    ("Bandwidth", (r"^fractional_bandwidth_vswr2$", r"^fractional_bandwidth_vswr2_with_surface_waves$",
                   r"^fractional_bandwidth_vswr2_one_mode$", r"^best_stack_bandwidth_vswr2$",
                   r"^fractional_bandwidth_estimate$", r"^bandwidth_ratio$", r"^fractional_bandwidth_3db$",
                   r"^axial_ratio_bandwidth_3db$", r"^design_bandwidth$")),
    ("Impedance", (r"^input_impedance_driving_point_ohm$", r"^input_impedance_ohm$",
                   r"^input_resistance_driving_point_ohm$", r"^input_resistance_ohm$", r"^inset_resistance_ohm$",
                   r"^feed_resistance_ohm$", r"^input_resistance_f_low_ohm$", r"^resonant_resistance_ohm$",
                   r"^edge_resistance_ohm$", r"^feed_impedance_each_element_ohm$",
                   r"^radiation_resistance_ohm$")),        # never a reference (a bicone's, an infinite sheet's)
    ("Beamwidth", (r"^hpbw_e_deg$", r"^hpbw_deg$", r"hpbw.*deg$", r"beamwidth.*deg$")),
)


def key_figures(metrics: dict, units: dict) -> list[tuple[str, str, str]]:
    """Up to four headline figures: (label, value text, metric key). A figure the
    design does not produce, or produces as NaN, is skipped. Patterns are tried in
    order and match whole names, so a loss or mutual resistance never stands in
    for the input impedance."""
    out = []
    for label, patterns in KEY_FIGURES:
        for pat in patterns:
            hit = next((k for k in metrics if _re.search(pat, k)), None)
            if hit is None:
                continue
            text = display_value(hit, metrics[hit], units.get(hit, ""), sig=4)
            if text is None:
                continue
            if label == "Gain" and "directivity" in hit:
                label = "Directivity"
            out.append((label, text, hit))
            break
    return out
