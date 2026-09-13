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
        haystack = (f"{archetype.key} {archetype.name} {archetype.family} "
                    f"{archetype.spec.summary}").lower()
        return self._needle in haystack


def requirement_fields(archetype) -> list:
    """Parameters the user must or may supply, ordered for a form.

    Requirements first because they have no defaults and the design cannot
    proceed without them; assumptions and materials after, since those carry
    defensible defaults and most users will leave them alone.
    """
    order = {"requirement": 0, "assumption": 1, "material": 2}
    fields = [p for p in archetype.spec.parameters if p.role in order]
    return sorted(fields, key=lambda p: (order[p.role], p.symbol != "f0", p.symbol))


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

    Frequency gets special handling: specs rarely declare a `typical` for f0
    (there is no universal answer), but an empty frequency field means the
    form opens unable to synthesise anything, which reads as broken.
    """
    if param.symbol == "f0" and archetype is not None:
        try:
            float(param.typical)
        except (TypeError, ValueError):
            return f"{default_frequency(archetype):.6g}"
    try:
        float(param.typical)
    except (TypeError, ValueError):
        return ""
    return str(param.typical)
