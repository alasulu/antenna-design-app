"""OTA Hub main window.

A navigation rail on the left and four pages: the antenna catalogue, the linear
array synthesiser, the planar array designer and the waveguide calculator. The
catalogue opens on a gallery - every archetype drawn from its own default design -
and a card opens that antenna's design page: requirements on the left, its
drawing, headline figures and the full results on the right, recomputed as the
requirements change. The form is generated from each spec's declared
parameters, so a new archetype in ``specs/`` gets a working page, and a drawing
from its family, with no code change here.
"""
from __future__ import annotations

import math
import re
import sys

import numpy as np
from PySide6.QtCore import QLocale, QPoint, QRect, QRectF, QSize, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QFont, QPainter, QPainterPath, QPen
from PySide6.QtWidgets import (QAbstractItemView, QApplication, QCheckBox, QComboBox, QDoubleSpinBox,
                               QFormLayout, QFrame, QGroupBox, QHBoxLayout, QHeaderView,
                               QLabel, QLayout, QLineEdit, QListWidget, QListWidgetItem, QMainWindow,
                               QPushButton, QScrollArea, QSizePolicy, QSpinBox, QSplitter,
                               QStackedWidget, QTableWidget, QTableWidgetItem, QTabWidget, QTextBrowser,
                               QVBoxLayout, QWidget)

from ..arrays import (TAPERS, grating_lobe_free_spacing_planar,
                      lattice_element_saving, planar_beam_cut,
                      planar_summarise, rectangular_lattice,
                      separable_weights, summarise, triangular_lattice, visible_grating_lobes)
from ..core import pattern as pat
from ..core.registry import Registry, default_registry
from ..core.units import engineering
from ..waveguides.rectangular import WR_SERIES, recommended_band, standard
from . import drawings
from .models import (default_for, display_value, format_input, humanize, is_primary, key_figures,
                     matches, parse_quantity, requirement_fields, shown_unit, suggested_value)
from .model3d import SolidPanel
from .plots import (Canvas, plot_element_layout, plot_hemisphere_cuts,
                    plot_polar, plot_sweep)
from .style import (ACCENT, FAINT, INK, LINE, MUTED, STYLE, SURFACE, WARN, WARN_SOFT,
                    dot_icon, rail_icon)

#: Archetypes whose far-field pattern we can compute from first principles.
#: Anything absent gets an honest message rather than a fabricated plot.
PATTERN_SOURCES = {
    "short_dipole": lambda: (pat.short_dipole(), "short dipole, sin^2"),
    "half_wave_dipole": lambda: (pat.finite_dipole(0.5), "half-wave dipole"),
    "resonant_dipole": lambda: (pat.finite_dipole(0.48), "resonant dipole, 0.48 lambda"),
    "quarter_wave_monopole": lambda: (_upper_half(pat.finite_dipole(0.5)),
                                      "monopole over ground, the upper half space"),
    "folded_dipole": lambda: (pat.finite_dipole(0.5), "folded dipole"),
    "small_circular_loop": lambda: (pat.short_dipole(), "small loop (dual of short dipole)"),
}

def _upper_half(p: "pat.Pattern") -> "pat.Pattern":
    """A pattern over a ground plane: the image makes the field, but nothing radiates
    below the plane, so the power (and the directivity) is the upper half space's."""
    return pat.Pattern(p.theta, p.phi, np.where((p.theta <= math.pi / 2)[:, None], p.U, 0.0))


FAMILY_NAMES = {"uwb": "UWB", "travelling_wave": "Travelling wave"}
#: The gallery's order: the everyday families first, the specialised ones after.
FAMILY_ORDER = ("wire", "loop", "patch", "horn", "reflector", "lens", "slot", "travelling_wave", "uwb", "dielectric")


def ordered_families(registry) -> list[str]:
    rank = {f: i for i, f in enumerate(FAMILY_ORDER)}
    return sorted(registry.families, key=lambda f: (rank.get(f, len(rank)), f))


def family_name(family: str) -> str:
    return FAMILY_NAMES.get(family, family.replace("_", " ").capitalize())


# ---------------------------------------------------------------- small building blocks

def _table(rows: int = 0, cols: int = 2) -> QTableWidget:
    t = QTableWidget(rows, cols)
    t.setHorizontalHeaderLabels(["quantity", "value"])
    t.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
    t.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
    t.horizontalHeader().setDefaultAlignment(Qt.AlignmentFlag.AlignLeft)
    t.verticalHeader().setVisible(False)
    t.setAlternatingRowColors(True)
    t.setShowGrid(False)
    t.setWordWrap(False)
    t.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
    t.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
    t.setFocusPolicy(Qt.FocusPolicy.NoFocus)
    return t


def _fill(table: QTableWidget, pairs: list[tuple]) -> None:
    """Rows of (name, value[, tooltip]). A value of None reads as outside the solved range."""
    table.setRowCount(len(pairs))
    mono = QFont()
    mono.setStyleHint(QFont.StyleHint.Monospace)
    mono.setFamilies(["SF Mono", "Menlo", "Consolas", "DejaVu Sans Mono", "monospace"])
    for row, pair in enumerate(pairs):
        name, value = pair[0], pair[1]
        tip = pair[2] if len(pair) > 2 else ""
        a = QTableWidgetItem(name)
        b = QTableWidgetItem(value if value is not None else "outside the solved range")
        if value is None:
            b.setForeground(QColor(FAINT))
        else:
            b.setFont(mono)
        b.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        if tip:
            a.setToolTip(tip)
            b.setToolTip(tip)
        table.setItem(row, 0, a)
        table.setItem(row, 1, b)


def _format(value, unit: str = "") -> str:
    text = display_value("", value, unit)
    return text if text is not None else "outside the solved range"


def _label(text: str = "", name: str = "", wrap: bool = False) -> QLabel:
    lab = QLabel(text)
    if name:
        lab.setObjectName(name)
    lab.setWordWrap(wrap)
    return lab


def _card(layout: QLayout | None = None, name: str = "card") -> QFrame:
    frame = QFrame()
    frame.setObjectName(name)
    if layout is not None:
        frame.setLayout(layout)
    return frame



def _shown(value, fmt: str, unit: str) -> str:
    """An array figure for a table: a beamwidth the pattern does not have (no half-power
    point, NaN) says so, and a sidelobe of -inf dB is no sidelobe at all."""
    if math.isnan(value):
        return "not defined"
    if math.isinf(value):
        return "none" if value < 0 and "dB" in unit else ("−∞" if value < 0 else "∞") + unit
    return f"{value:{fmt}}{unit}"

class FlowLayout(QLayout):
    """Lays its items out left to right and wraps them, like text."""

    def __init__(self, parent=None, spacing: int = 14, centred: bool = False) -> None:
        super().__init__(parent)
        self._items = []
        self._spacing = spacing
        self.centred = centred

    def addItem(self, item) -> None:  # noqa: N802 - Qt API
        self._items.append(item)

    def count(self) -> int:
        return len(self._items)

    def itemAt(self, i: int):  # noqa: N802 - Qt API
        return self._items[i] if 0 <= i < len(self._items) else None

    def takeAt(self, i: int):  # noqa: N802 - Qt API
        return self._items.pop(i) if 0 <= i < len(self._items) else None

    def expandingDirections(self):  # noqa: N802 - Qt API
        return Qt.Orientation(0)

    def hasHeightForWidth(self) -> bool:  # noqa: N802 - Qt API
        return True

    def heightForWidth(self, width: int) -> int:  # noqa: N802 - Qt API
        return self._arrange(QRect(0, 0, width, 0), True)

    def setGeometry(self, rect) -> None:  # noqa: N802 - Qt API
        super().setGeometry(rect)
        self._arrange(rect, False)

    def sizeHint(self) -> QSize:  # noqa: N802 - Qt API
        return self.minimumSize()

    def minimumSize(self) -> QSize:  # noqa: N802 - Qt API
        size = QSize()
        for item in self._items:
            size = size.expandedTo(item.minimumSize())
        m = self.contentsMargins()
        return size + QSize(m.left() + m.right(), m.top() + m.bottom())

    def _arrange(self, rect, test: bool) -> int:
        m = self.contentsMargins()
        left = rect.x() + m.left()
        right = rect.right() - m.right()
        if self.centred and self._items:
            w = self._items[0].sizeHint().width()
            n = max(1, (right - left + self._spacing) // (w + self._spacing))
            left += max(0, (right - left - (n * w + (n - 1) * self._spacing)) // 2)
        x, y, line = left, rect.y() + m.top(), 0
        for item in self._items:
            if item.widget() is not None and not item.widget().isVisibleTo(item.widget().parentWidget()):
                continue
            hint = item.sizeHint()
            if x + hint.width() > right and line > 0:
                x, y, line = left, y + line + self._spacing, 0
            if not test:
                item.setGeometry(QRect(QPoint(x, y), hint))
            x += hint.width() + self._spacing
            line = max(line, hint.height())
        return y + line - rect.y() + m.bottom()


# ---------------------------------------------------------------- the gallery

def default_values(archetype) -> dict:
    """The form's default requirements, as numbers."""
    out = {}
    for param in requirement_fields(archetype):
        try:
            out[param.symbol] = float(default_for(param, archetype))
        except ValueError:
            pass
    return out


def design_values(design) -> dict:
    values = dict(design.requirements)
    values.update(design.parameters)
    values.update(design.metrics)
    return values


class AntennaCard(QWidget):
    """One archetype in the gallery: its drawing, its name, its family."""

    clicked = Signal(str)
    W, H, PIC = 216, 206, 132

    def __init__(self, archetype, pixmap) -> None:
        super().__init__()
        self.archetype = archetype
        self.pixmap = pixmap
        self.hover = False
        self.setFixedSize(self.W, self.H)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setToolTip(archetype.spec.summary)
        self.setAttribute(Qt.WidgetAttribute.WA_Hover, True)

    def enterEvent(self, event) -> None:  # noqa: N802 - Qt API
        self.hover = True
        self.update()

    def leaveEvent(self, event) -> None:  # noqa: N802 - Qt API
        self.hover = False
        self.update()

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802 - Qt API
        if event.button() == Qt.MouseButton.LeftButton:
            self.clicked.emit(self.archetype.key)

    def paintEvent(self, event) -> None:  # noqa: N802 - Qt API
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        r = QRectF(self.rect()).adjusted(1, 1, -1, -1)
        colour = drawings.family_colour(self.archetype.family)
        path = QPainterPath()
        path.addRoundedRect(r, 12, 12)
        p.fillPath(path, QColor(SURFACE))
        p.save()
        p.setClipPath(path)
        p.drawPixmap(QRectF(r.left(), r.top(), r.width(), self.PIC), self.pixmap, QRectF(self.pixmap.rect()))
        p.fillRect(QRectF(r.left(), r.top() + self.PIC - 3, r.width(), 3), colour)
        p.restore()
        p.setPen(QPen(QColor(ACCENT) if self.hover else QColor(LINE), 1.6 if self.hover else 1.0))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawPath(path)
        name_font = QFont(self.font())
        name_font.setPointSizeF(self.font().pointSizeF() * 1.0)
        name_font.setWeight(QFont.Weight.DemiBold)
        p.setFont(name_font)
        p.setPen(QColor(INK))
        text_rect = QRectF(r.left() + 12, r.top() + self.PIC + 8, r.width() - 24, 38)
        p.drawText(text_rect, Qt.TextFlag.TextWordWrap | Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop,
                   self.archetype.name)
        small = QFont(self.font())
        small.setPointSizeF(self.font().pointSizeF() * 0.86)
        p.setFont(small)
        p.setPen(colour)
        p.drawText(QRectF(r.left() + 12, r.bottom() - 24, r.width() - 24, 18),
                   Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter, family_name(self.archetype.family))
        conf = self.archetype.spec.confidence
        if conf in ("low", "high"):
            label = "low confidence" if conf == "low" else "high confidence"
            fm = p.fontMetrics()
            w = fm.horizontalAdvance(label) + 14
            badge = QRectF(r.right() - 12 - w, r.bottom() - 25, w, 19)
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(QColor(WARN_SOFT) if conf == "low" else QColor("#e0f0e7"))
            p.drawRoundedRect(badge, 9.5, 9.5)
            p.setPen(QColor(WARN) if conf == "low" else QColor("#2d7a55"))
            p.drawText(badge, Qt.AlignmentFlag.AlignCenter, label)
        p.end()


class DrawingView(QWidget):
    """The design's technical drawing, repainted from its current geometry."""

    def __init__(self) -> None:
        super().__init__()
        self.key, self.family, self.values = "", "", {}
        self.setMinimumHeight(300)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

    def show_design(self, key: str, family: str, values: dict) -> None:
        self.key, self.family, self.values = key, family, values
        self.update()

    def paintEvent(self, event) -> None:  # noqa: N802 - Qt API
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        r = QRectF(self.rect())
        path = QPainterPath()
        path.addRoundedRect(r.adjusted(0.5, 0.5, -0.5, -0.5), 12, 12)
        p.fillPath(path, QColor("#fbfcfd"))
        p.setPen(QPen(QColor(LINE), 1))
        p.drawPath(path)
        if self.key:
            drawings.paint(p, r.adjusted(8, 8, -8, -30), self.key, self.family, self.values)
            self._legend(p, r)
        p.end()

    def _legend(self, p: QPainter, r: QRectF) -> None:
        items = (("copper", drawings.COPPER), ("substrate", drawings.SUBSTRATE), ("dielectric", drawings.CERAMIC),
                 ("metal", drawings.METAL), ("feed", drawings.FEED))
        small = QFont(self.font())
        small.setPointSizeF(self.font().pointSizeF() * 0.82)
        p.setFont(small)
        fm = p.fontMetrics()
        x, y = r.left() + 14, r.bottom() - 16
        for name, colour in items:
            p.setPen(QPen(QColor(colour).darker(125), 1))
            p.setBrush(colour)
            p.drawEllipse(QRectF(x, y - 5, 9, 9))
            p.setPen(QColor(MUTED))
            p.drawText(QRectF(x + 13, y - 8, 90, 16), Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter, name)
            x += 13 + fm.horizontalAdvance(name) + 16
        note = "drawn from this design's dimensions"
        p.drawText(QRectF(r.right() - 14 - fm.horizontalAdvance(note), y - 8, fm.horizontalAdvance(note), 16),
                   Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter, note)


class FigureTile(QFrame):
    def __init__(self) -> None:
        super().__init__()
        self.setObjectName("tile")
        lay = QVBoxLayout(self)
        lay.setContentsMargins(14, 10, 14, 10)
        lay.setSpacing(2)
        self.label = _label("", "tilelabel")
        self.value = _label("", "tilevalue")
        self.key = _label("", "tilekey")
        for w in (self.label, self.value, self.key):
            lay.addWidget(w)

    def set(self, label: str, value: str, key: str) -> None:
        self.label.setText(label.upper())
        self.value.setText(value)
        self.key.setText(humanize(key))
        self.setToolTip(key)


# ---------------------------------------------------------------- catalogue

class CatalogueTab(QWidget):
    """The gallery of every archetype, and the design page of the one chosen."""

    def __init__(self, registry: Registry) -> None:
        super().__init__()
        self.registry = registry
        self.current = None
        self._fields: dict[str, QLineEdit] = {}
        self._errors: dict[str, QLabel] = {}
        self._design = None
        self._family = "all"
        # a symbol one spec leaves undeclared is often named in another (t: "thickness parameter")
        self._names: dict[str, str] = {}
        for archetype in registry:
            for param in archetype.spec.parameters:
                if param.name:
                    self._names.setdefault(param.symbol, param.name)
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(260)
        self._timer.timeout.connect(self._synthesise)

        self.pages = QStackedWidget()
        self.pages.addWidget(self._build_gallery())
        self.pages.addWidget(self._build_design_page())
        QVBoxLayout(self).addWidget(self.pages)
        self.layout().setContentsMargins(0, 0, 0, 0)

    # -------------------------------------------------------------- gallery
    def _build_gallery(self) -> QWidget:
        page = QWidget()
        page.setObjectName("page")
        outer = QVBoxLayout(page)
        outer.setContentsMargins(28, 22, 28, 12)
        outer.setSpacing(14)

        head = QHBoxLayout()
        titles = QVBoxLayout()
        titles.setSpacing(2)
        titles.addWidget(_label("Antennas", "h1"))
        titles.addWidget(_label(f"{len(self.registry)} designs in {len(self.registry.families)} families. "
                                "Pick one to size it for your frequency.", "muted"))
        head.addLayout(titles)
        head.addStretch(1)
        self.search = QLineEdit(placeholderText="Search: patch, horn, Yagi, spiral…")
        self.search.setObjectName("search")
        self.search.setClearButtonEnabled(True)
        self.search.textChanged.connect(self._on_search)
        head.addWidget(self.search, 0, Qt.AlignmentFlag.AlignVCenter)
        outer.addLayout(head)

        chip_host = QWidget()
        chips = FlowLayout(chip_host, spacing=6)
        chips.setContentsMargins(0, 0, 0, 0)
        self.chips: dict[str, QPushButton] = {}
        for fam in ["all"] + ordered_families(self.registry):
            n = len(self.registry) if fam == "all" else len(self.registry.by_family(fam))
            b = QPushButton(f"All  {n}" if fam == "all" else f"{family_name(fam)}  {n}")
            b.setObjectName("chip")
            b.setCheckable(True)
            b.setChecked(fam == "all")
            if fam != "all":
                b.setIcon(dot_icon(drawings.FAMILY_COLOURS.get(fam, MUTED)))
            b.clicked.connect(lambda _=False, f=fam: self._on_family(f))
            chips.addWidget(b)
            self.chips[fam] = b
        outer.addWidget(chip_host)

        self.cards: dict[str, AntennaCard] = {}
        grid_host = QWidget()
        self.flow = FlowLayout(grid_host, spacing=16, centred=True)
        self.flow.setContentsMargins(0, 4, 0, 16)
        for family in ordered_families(self.registry):
            for archetype in self.registry.by_family(family):
                card = AntennaCard(archetype, self._thumbnail(archetype))
                card.clicked.connect(self.select_key)
                self.flow.addWidget(card)
                self.cards[archetype.key] = card
        self.empty = _label("No antenna matches that search.", "muted")
        self.empty.hide()
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(grid_host)
        outer.addWidget(self.empty)
        outer.addWidget(scroll, 1)
        return page

    def _thumbnail(self, archetype):
        try:
            values = design_values(archetype.synthesize(**default_values(archetype)))
        except Exception:  # noqa: BLE001 - a card still shows its family's picture
            values = {}
        return drawings.thumbnail(archetype.key, archetype.family, values, AntennaCard.W - 2, AntennaCard.PIC)

    def _apply_filter(self) -> None:
        needle = self.search.text()
        shown = 0
        for key, card in self.cards.items():
            a = self.registry[key]
            ok = (self._family == "all" or a.family == self._family) and matches(a, needle)
            card.setVisible(ok)
            shown += ok
        self.empty.setVisible(shown == 0)
        self.flow.invalidate()

    def _on_search(self, text: str) -> None:
        self._apply_filter()

    def _on_family(self, family: str) -> None:
        self._family = family
        for fam, b in self.chips.items():
            b.setChecked(fam == family)
        self._apply_filter()

    # -------------------------------------------------------------- design page
    def _build_design_page(self) -> QWidget:
        page = QWidget()
        page.setObjectName("page")
        outer = QVBoxLayout(page)
        outer.setContentsMargins(28, 14, 28, 12)
        outer.setSpacing(10)

        back = QPushButton("←  All antennas")
        back.setObjectName("back")
        back.clicked.connect(self.show_gallery)
        outer.addWidget(back, 0, Qt.AlignmentFlag.AlignLeft)

        head = QHBoxLayout()
        head.setSpacing(12)
        self.title = _label("select an archetype", "h1")
        self.family_chip = _label("")
        head.addWidget(self.title)
        head.addWidget(self.family_chip, 0, Qt.AlignmentFlag.AlignVCenter)
        head.addStretch(1)
        outer.addLayout(head)
        self.summary = _label("", "muted", wrap=True)
        self.band = _label("", "faint")
        outer.addWidget(self.summary)
        outer.addWidget(self.band)

        split = QSplitter(Qt.Orientation.Horizontal)
        split.setChildrenCollapsible(False)
        split.addWidget(self._build_form_column())
        split.addWidget(self._build_results_column())
        split.setStretchFactor(0, 0)
        split.setStretchFactor(1, 1)
        split.setSizes([360, 900])
        outer.addWidget(split, 1)
        return page

    def _build_form_column(self) -> QWidget:
        col = QFrame()
        col.setObjectName("card")
        col.setMinimumWidth(300)
        col.setMaximumWidth(460)
        lay = QVBoxLayout(col)
        lay.setContentsMargins(16, 14, 16, 14)
        lay.setSpacing(10)
        lay.addWidget(_label("Requirements", "h2"))
        lay.addWidget(_label("Type values with units if you like: 2.4 GHz, 1.6 mm, 50 ohm. "
                             "Results update as you type.", "faint", wrap=True))
        self.examples = QComboBox()
        self.examples.activated.connect(self._load_example)
        lay.addWidget(self.examples)

        self.form_host = QWidget()
        self.form = QVBoxLayout(self.form_host)
        self.form.setContentsMargins(0, 4, 0, 0)
        self.form.setSpacing(12)
        self.more_button = QPushButton("Show more settings")
        self.more_button.setObjectName("link")
        self.more_button.clicked.connect(self._toggle_more)
        self.more_host = QWidget()
        self.more = QVBoxLayout(self.more_host)
        self.more.setContentsMargins(0, 0, 0, 0)
        self.more.setSpacing(12)
        self.more_host.hide()

        inner = QWidget()
        il = QVBoxLayout(inner)
        il.setContentsMargins(0, 0, 4, 0)
        il.addWidget(self.form_host)
        il.addWidget(self.more_button, 0, Qt.AlignmentFlag.AlignLeft)
        il.addWidget(self.more_host)
        il.addStretch(1)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(inner)
        lay.addWidget(scroll, 1)

        self.synth_button = QPushButton("Recalculate")
        self.synth_button.setObjectName("primary")
        self.synth_button.clicked.connect(self._synthesise)
        self.synth_button.setEnabled(False)
        lay.addWidget(self.synth_button)
        return col

    def _build_results_column(self) -> QWidget:
        col = QWidget()
        lay = QVBoxLayout(col)
        lay.setContentsMargins(8, 0, 0, 0)
        lay.setSpacing(12)
        self.drawing = DrawingView()
        lay.addWidget(self.drawing, 5)

        self.banner = _card(QVBoxLayout(), "banner")
        self.banner_text = _label("", wrap=True)
        self.banner_text.setTextFormat(Qt.TextFormat.RichText)
        self.banner.layout().setContentsMargins(12, 8, 12, 8)
        self.banner.layout().addWidget(self.banner_text)
        self.banner.hide()
        lay.addWidget(self.banner)

        tiles = QHBoxLayout()
        tiles.setSpacing(12)
        self.tiles = [FigureTile() for _ in range(4)]
        for t in self.tiles:
            tiles.addWidget(t)
        lay.addLayout(tiles)

        self.geometry_table = _table()
        self.metrics_table = _table()
        self.notes = QTextBrowser()
        self.notes.setOpenExternalLinks(False)

        self.metric_picker = QComboBox()
        self.metric_picker.currentIndexChanged.connect(
            lambda _i: self._sweep(self.metric_picker.currentData() or ""))
        self.sweep_canvas = Canvas()
        sweep = QWidget()
        sl = QVBoxLayout(sweep)
        sl.setContentsMargins(4, 8, 4, 4)
        row = QHBoxLayout()
        row.addWidget(_label("Quantity", "muted"))
        row.addWidget(self.metric_picker, 1)
        sl.addLayout(row)
        sl.addWidget(self.sweep_canvas, 1)

        self.pattern_canvas = Canvas(polar=True)
        self.solid_panel = SolidPanel()

        self.results = QTabWidget()
        self.results.addTab(self._wrap(self.geometry_table), "Dimensions")
        self.results.addTab(self._wrap(self.metrics_table), "Performance")
        self.results.addTab(self._wrap(sweep), "Sweep")
        self.results.addTab(self._wrap(self.pattern_canvas), "Pattern")
        self.results.addTab(self._wrap(self.solid_panel), "3D model")
        self.results.addTab(self._wrap(self.notes), "Notes and validity")
        lay.addWidget(self.results, 6)
        return col

    @staticmethod
    def _wrap(widget: QWidget) -> QFrame:
        lay = QVBoxLayout()
        lay.setContentsMargins(6, 6, 6, 6)
        lay.addWidget(widget)
        return _card(lay)

    # -------------------------------------------------------------- actions
    def show_gallery(self) -> None:
        self._timer.stop()                 # a pending recalculation belongs to the page being left
        self.pages.setCurrentIndex(0)

    def select_key(self, key: str) -> None:
        """Open an archetype's design page. Separated out so tests can drive it."""
        archetype = self.registry[key]
        self.current = archetype
        spec = archetype.spec
        self.title.setText(spec.name)
        colour = drawings.FAMILY_COLOURS.get(archetype.family, MUTED)
        self.family_chip.setText(family_name(archetype.family))
        self.family_chip.setStyleSheet(f"color: {colour}; border: 1px solid {colour}; border-radius: 10px; "
                                       "padding: 2px 10px; font-weight: 600;")
        self.summary.setText(spec.summary)
        lo, hi = spec.freq_range_hz
        band = f"Valid {engineering(lo, 'Hz')} – {engineering(hi, 'Hz')}"
        if spec.confidence == "low":
            band += "   ·   low confidence: verify in a full-wave solver"
        self.band.setText(band)

        for layout in (self.form, self.more):
            while layout.count():
                w = layout.takeAt(0).widget()
                if w is not None:                 # detach now: deleteLater alone leaves it painted a moment
                    w.hide()
                    w.setParent(None)
                    w.deleteLater()
        self._fields.clear()
        self._errors.clear()
        n_more = 0
        syms = {p.symbol for p in requirement_fields(archetype)}
        primary = {p.symbol for p in requirement_fields(archetype) if is_primary(p)}
        # the known design that fills the most of the page, the first among equals
        fitting = [c for c in spec.known_cases if c.given and set(c.given) <= syms]
        first_case = max(fitting, key=lambda c: len(set(c.given) & primary), default=None)
        for param in requirement_fields(archetype):
            box = QWidget()
            bl = QVBoxLayout(box)
            bl.setContentsMargins(0, 0, 0, 0)
            bl.setSpacing(3)
            top = QHBoxLayout()
            name = param.name or param.symbol
            name = humanize(name) if "_" in name else name[:1].upper() + name[1:]
            top.addWidget(_label(name + (" *" if is_primary(param) else ""), "fieldname"))
            top.addStretch(1)
            meta = param.symbol + (f"  ·  {shown_unit(param.unit)}" if param.unit and param.unit != "-" else "")
            top.addWidget(_label(meta, "fieldmeta"))
            bl.addLayout(top)
            text = default_for(param, archetype)
            if first_case is not None and is_primary(param) and param.symbol in first_case.given:
                text = str(first_case.given[param.symbol])      # a real design, not a band's midpoint
            try:
                text = format_input(float(text), param.unit)
            except ValueError:
                pass
            edit = QLineEdit(text)
            hint = param.description or param.name
            edit.setPlaceholderText("required" if is_primary(param) else "default")
            edit.setToolTip(f"{param.role}: {hint}")
            edit.textEdited.connect(lambda _=None: self._timer.start())
            edit.returnPressed.connect(self._synthesise)
            bl.addWidget(edit)
            err = _label("", "error")
            err.hide()
            bl.addWidget(err)
            if param.description and param.role != "requirement":
                bl.addWidget(_label(param.description, "fieldmeta", wrap=True))
            (self.form if is_primary(param) else self.more).addWidget(box)
            n_more += not is_primary(param)
            self._fields[param.symbol] = edit
            self._errors[param.symbol] = err
        self.more_button.setVisible(n_more > 0)
        self.more_button.setText(f"Show {n_more} more settings" if not self.more_host.isVisible() else "Hide settings")

        self.examples.blockSignals(True)
        self.examples.clear()
        self.examples.addItem("Load a known design…")
        syms = set(self._fields)
        for i, case in enumerate(spec.known_cases):
            if case.given and set(case.given) <= syms:
                self.examples.addItem(self._case_label(case), i)
        self.examples.setVisible(self.examples.count() > 1)
        self.examples.blockSignals(False)

        self._suggest_blanks()
        self.synth_button.setEnabled(True)
        self.pages.setCurrentIndex(1)
        self._clear_results()
        self._update_pattern()
        self._synthesise()

    def _suggest_blanks(self) -> None:
        """Fill a field no known design and no plain default reaches from the middle
        of its spec's typical range, worked out on the design so far - so a loop's
        page opens complete instead of asking for a wire radius. Twice, since one
        suggestion can unlock another's range."""
        params = {p.symbol: p for p in self.current.spec.parameters}
        for _ in range(2):
            blanks = [s for s, e in self._fields.items() if not e.text().strip()]
            if not blanks:
                return
            values, _ = self._read_form()
            try:
                design = self.current.synthesize(**values)
            except Exception:  # noqa: BLE001 - suggestions are a courtesy
                return
            known = design_values(design)
            for sym in blanks:
                v = suggested_value(params[sym], known)
                if v is not None:
                    edit = self._fields[sym]
                    edit.setText(format_input(v, params[sym].unit))
                    edit.setToolTip(edit.toolTip() + f"\nsuggested: the middle of {params[sym].typical}")

    def _case_label(self, case) -> str:
        units = {p.symbol: p.unit for p in self.current.spec.parameters}
        parts = [f"{k} {display_value(k, v, units.get(k, '')) or v}" for k, v in list(case.given.items())[:3]]
        return ", ".join(parts)

    def _load_example(self, index: int) -> None:
        case_index = self.examples.itemData(index)
        if case_index is None:
            return
        case = self.current.spec.known_cases[case_index]
        params = {p.symbol: p for p in self.current.spec.parameters}
        for sym, edit in self._fields.items():
            # every field the case leaves out goes back to its default: a setting
            # changed by hand is not part of the published design
            text = str(case.given[sym]) if sym in case.given else default_for(params[sym], self.current)
            try:
                text = format_input(float(text), params[sym].unit)
            except ValueError:
                pass
            edit.setText(text)
        # a free size the case does not give (a loop's wire radius) gets the same
        # suggestion it has when the page opens, not a blank that stops the design
        self._suggest_blanks()
        self._synthesise()

    def _toggle_more(self) -> None:
        self.more_host.setVisible(not self.more_host.isVisible())
        n = self.more.count()
        self.more_button.setText("Hide settings" if self.more_host.isVisible() else f"Show {n} more settings")

    def _clear_results(self) -> None:
        _fill(self.geometry_table, [])
        _fill(self.metrics_table, [])
        self.notes.setHtml("")
        self.metric_picker.clear()
        self.sweep_canvas.message("synthesise a design to sweep it")
        for t in self.tiles:
            t.hide()
        self.banner.hide()

    def collect(self) -> dict[str, float]:
        """Read the form. Blank fields are omitted, not zeroed; units are understood."""
        out: dict[str, float] = {}
        units = {p.symbol: p.unit for p in self.current.spec.parameters} if self.current else {}
        for symbol, edit in self._fields.items():
            text = edit.text().strip()
            if not text:
                continue
            try:
                out[symbol] = parse_quantity(text, units.get(symbol, ""))
            except ValueError:
                raise ValueError(f"{symbol}: {text!r} is not a number") from None
        return out

    def _read_form(self) -> tuple[dict[str, float], list[str]]:
        """Like collect, but marks bad fields instead of stopping at the first."""
        out, bad = {}, []
        units = {p.symbol: p.unit for p in self.current.spec.parameters}
        for symbol, edit in self._fields.items():
            text = edit.text().strip()
            err = self._errors[symbol]
            ok = True
            if text:
                try:
                    out[symbol] = parse_quantity(text, units.get(symbol, ""))
                except ValueError:
                    ok = False
                    bad.append(symbol)
                    unit = shown_unit(units.get(symbol, ""))
                    err.setText(f"Enter a number{f', in {unit} or with a unit' if unit else ''}.")
            err.setVisible(not ok)
            edit.setProperty("invalid", "false" if ok else "true")
            edit.style().unpolish(edit)
            edit.style().polish(edit)
        return out, bad

    def _synthesise(self) -> None:
        if self.current is None:
            return
        values, bad = self._read_form()
        try:
            design = self.current.synthesize(**values)
        except Exception as exc:  # noqa: BLE001 - surfaced on the page
            self._show_banner([f"This design could not be computed: {exc}"])
            return
        self._design = design
        units = design.units
        spec = self.current.spec
        notes = {r.output: r.notes for r in spec.synthesis}
        notes.update({r.metric: r.notes for r in spec.analysis})
        names = dict(self._names)
        names.update({p.symbol: p.name for p in spec.parameters if p.name})

        _fill(self.geometry_table, [
            (names.get(name) or humanize(name), display_value(name, value, units.get(name, "")),
             f"{name}\n\n{notes.get(name, '')}".strip())
            for name, value in design.parameters.items()
            if name not in ("lambda0", "k0") and name not in self._fields])
        _fill(self.metrics_table, [
            (humanize(name), display_value(name, value, units.get(name, "")),
             f"{name}\n\n{notes.get(name, '')}".strip())
            for name, value in design.metrics.items()])

        figures = key_figures(design.metrics, units)
        for tile, fig in zip(self.tiles, figures):
            tile.set(*fig)
            tile.show()
        for tile in self.tiles[len(figures):]:
            tile.hide()

        self.drawing.show_design(self.current.key, self.current.family, design_values(design))
        self.solid_panel.set_design(design)

        messages = []
        if bad:
            messages.append(f"Check {', '.join(bad)}: not a number this page can read.")
        missing = design.missing_requirements()
        if missing:
            messages.append(f"Add {', '.join(missing)} to complete the design.")
        for w in design.warnings:
            if "low-confidence" not in w:
                messages.append(w)
        if spec.confidence == "low":
            messages.append("Low confidence: treat these numbers as indicative and verify in a full-wave solver.")
        self._show_banner(messages)
        self._write_notes(design, missing)

        numeric = [k for k, v in design.metrics.items()
                   if isinstance(v, (int, float)) and not isinstance(v, bool)
                   and math.isfinite(float(v))]
        keep = self.metric_picker.currentData()
        self.metric_picker.blockSignals(True)
        self.metric_picker.clear()
        for key in sorted(numeric, key=humanize):
            self.metric_picker.addItem(humanize(key), key)
            self.metric_picker.setItemData(self.metric_picker.count() - 1, key, Qt.ItemDataRole.ToolTipRole)
        pick = keep if keep in numeric else next(
            (k for k in numeric if re.search(r"(gain|directivity)_dbi$", k)), None)
        if pick:
            self.metric_picker.setCurrentIndex(self.metric_picker.findData(pick))
        self.metric_picker.blockSignals(False)
        if numeric:
            self._sweep(self.metric_picker.currentData())
        self._update_pattern()

    def _show_banner(self, messages: list[str]) -> None:
        if not messages:
            self.banner.hide()
            return
        self.banner_text.setText("<br>".join(f"•&nbsp;{m}" for m in messages))
        self.banner.show()

    def _write_notes(self, design, missing) -> None:
        parts = [f"<h3 style='margin:4px 0 6px 0'>Valid only under</h3>"]
        if design.notes:
            parts.append("<ul style='margin-top:0'>" + "".join(f"<li>{n}</li>" for n in design.notes) + "</ul>")
        else:
            parts.append(f"<p style='color:{MUTED}'>No stated limits beyond the frequency band.</p>")
        if design.warnings:
            parts.append("<h3 style='margin:12px 0 6px 0'>Warnings</h3><ul style='margin-top:0'>"
                         + "".join(f"<li>{w}</li>" for w in design.warnings) + "</ul>")
        refs = self.current.spec.references
        if refs:
            parts.append("<h3 style='margin:12px 0 6px 0'>References</h3><ol style='margin-top:0'>"
                         + "".join(f"<li style='color:{MUTED}'>{r}</li>" for r in refs) + "</ol>")
        self.notes.setHtml("".join(parts))

    def _sweep(self, metric: str) -> None:
        """Re-synthesise across the validity band and plot one metric."""
        if not metric or self.current is None or self._design is None:
            return
        base = dict(self._design.requirements)
        freq_key = next((p.symbol for p in requirement_fields(self.current)
                         if p.unit == "Hz" and p.role == "requirement" and p.symbol in base), None)
        if freq_key is None:
            self.sweep_canvas.message(
                "this archetype is not parameterised by a frequency, so there is nothing to sweep")
            return
        lo, hi = self.current.spec.freq_range_hz
        centre = base[freq_key]
        lo = max(lo, centre / 4.0)
        hi = min(hi if math.isfinite(hi) else centre * 4.0, centre * 4.0)
        if not hi > lo:
            self.sweep_canvas.message("validity band is too narrow to sweep")
            return
        freqs, values = [], []
        for f in np.geomspace(lo, hi, 160):
            trial = dict(base, **{freq_key: float(f)})
            try:
                result = self.current.synthesize(**trial)
            except Exception:  # noqa: BLE001 - a point that will not solve is skipped
                continue
            value = result.metrics.get(metric)
            if isinstance(value, (int, float)) and math.isfinite(float(value)):
                freqs.append(f)
                values.append(float(value))
        if len(freqs) < 2:
            self.sweep_canvas.message(f"{metric} could not be evaluated across the band")
            return
        plot_sweep(self.sweep_canvas, freqs, values, humanize(metric),
                   shown_unit(self._design.units.get(metric, "")), self.current.spec.freq_range_hz)
        ax = self.sweep_canvas.axes
        ax.axvline(centre / 1e9, color=INK, lw=0.8, ls=":")
        self.sweep_canvas.draw_idle()

    def _update_pattern(self) -> None:
        if self.current is None:
            return
        factory = PATTERN_SOURCES.get(self.current.key)
        if factory is None:
            self.pattern_canvas.message(
                f"No closed-form pattern is implemented for {self.current.name}.\n\n"
                "Only archetypes whose far field this toolkit can compute from first "
                "principles are plotted here; the rest would need a full-wave solve, "
                "and drawing a plausible-looking guess would be worse than drawing "
                "nothing.")
            return
        pattern, label = factory()
        cut = pattern.cut(0.0)
        plot_polar(self.pattern_canvas, pattern.theta, cut,
                   f"{label}\nD = {pattern.directivity_dbi():.2f} dBi")


# ------------------------------------------------------------------- arrays

def _page(title: str, subtitle: str, body: QWidget) -> tuple[QWidget, QVBoxLayout]:
    page = QWidget()
    page.setObjectName("page")
    lay = QVBoxLayout(page)
    lay.setContentsMargins(28, 22, 28, 12)
    lay.setSpacing(12)
    lay.addWidget(_label(title, "h1"))
    lay.addWidget(_label(subtitle, "muted", wrap=True))
    lay.addWidget(body, 1)
    return page, lay


def _results_tabs(*pairs) -> QTabWidget:
    tabs = QTabWidget()
    for widget, name in pairs:
        tabs.addTab(CatalogueTab._wrap(widget), name)
    return tabs


class ArrayTab(QWidget):
    """Linear array taper designer with a live pattern."""

    def __init__(self) -> None:
        super().__init__()
        self.n = QSpinBox(minimum=2, maximum=256, value=16)
        self.taper = QComboBox()
        self.taper.addItems(sorted(TAPERS))
        self.taper.setCurrentText("chebyshev")
        self.sll = QDoubleSpinBox(minimum=-80.0, maximum=-5.0, value=-30.0, singleStep=5.0)
        self.sll.setSuffix(" dB")
        self.spacing = QDoubleSpinBox(minimum=0.05, maximum=5.0, value=0.5, singleStep=0.05)
        self.spacing.setSuffix(" λ")
        self.scan = QDoubleSpinBox(minimum=0.0, maximum=180.0, value=90.0, singleStep=5.0)
        self.scan.setSuffix("°")
        self.log_scale = QCheckBox("log frequency axis")
        self.log_scale.setVisible(False)

        controls = QFormLayout()
        controls.setVerticalSpacing(10)
        controls.addRow("Elements", self.n)
        controls.addRow("Taper", self.taper)
        controls.addRow("Sidelobe level", self.sll)
        controls.addRow("Spacing", self.spacing)
        controls.addRow("Scan angle", self.scan)
        box = QGroupBox("Array")
        box.setLayout(controls)

        self.summary_table = _table()
        self.canvas = Canvas(polar=True)
        self.weights_canvas = Canvas()

        left = QVBoxLayout()
        left.setContentsMargins(0, 0, 0, 0)
        left.addWidget(box)
        results = QGroupBox("Result")
        rl = QVBoxLayout(results)
        rl.addWidget(self.summary_table)
        left.addWidget(results, 1)
        left_widget = QWidget()
        left_widget.setLayout(left)
        left_widget.setMaximumWidth(420)

        right = _results_tabs((self.canvas, "Pattern"), (self.weights_canvas, "Excitation"))

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(left_widget)
        splitter.addWidget(right)
        splitter.setSizes([340, 760])
        page, _ = _page("Linear arrays", "Choose a taper, sidelobe level, spacing and scan: the beam, its "
                        "sidelobes and any grating lobe update as you go.", splitter)
        QVBoxLayout(self).addWidget(page)
        self.layout().setContentsMargins(0, 0, 0, 0)

        for widget in (self.n, self.sll, self.spacing, self.scan):
            widget.valueChanged.connect(self.refresh)
        self.taper.currentTextChanged.connect(self.refresh)
        self.refresh()

    def weights(self):
        name = self.taper.currentText()
        if name in ("chebyshev", "taylor"):
            return TAPERS[name](self.n.value(), self.sll.value())
        return TAPERS[name](self.n.value())

    def refresh(self) -> None:
        uses_sll = self.taper.currentText() in ("chebyshev", "taylor")
        self.sll.setEnabled(uses_sll)
        try:
            w = self.weights()
            s = summarise(w, self.spacing.value(), self.scan.value())
        except Exception as exc:  # noqa: BLE001
            self.canvas.message(f"could not synthesise: {exc}")
            return

        rows = [
            ("elements", str(s["elements"])),
            ("directivity", _shown(s["directivity_dbi"], ".2f", " dBi")),
            ("half-power beamwidth", _shown(s["hpbw_deg"], ".3f", "°")),
            ("first sidelobe", _shown(s["sidelobe_db"], ".2f", " dB")),
            ("taper efficiency", f"{s['taper_efficiency']:.4f}"),
            ("directivity given up", f"{-10 * math.log10(s['taper_efficiency']):.2f} dB"),
            ("max spacing, no grating lobe", f"{s['max_spacing_no_grating']:.4f} λ"),
            ("grating lobe present", "YES" if s["grating_lobe"] else "no"),
        ]
        _fill(self.summary_table, rows)

        from ..arrays import array_pattern
        p = array_pattern(w, self.spacing.value(), self.scan.value())
        plot_polar(self.canvas, p.theta, p.cut(0.0),
                   f"{self.taper.currentText()} taper, {s['elements']} elements\n"
                   f"SLL {s['sidelobe_db']:.1f} dB, HPBW {s['hpbw_deg']:.2f}°")

        self.weights_canvas.clear()
        ax = self.weights_canvas.axes
        ax.set_axis_on()
        markers, stems, base = ax.stem(range(len(w)), w)
        markers.set_color(ACCENT)
        stems.set_color(ACCENT)
        base.set_color(LINE)
        ax.set_xlabel("element")
        ax.set_ylabel("normalised amplitude")
        self.weights_canvas.draw_idle()


class PlanarArrayTab(QWidget):
    """Planar lattice designer: layout, steering, and the cuts through the beam."""

    def __init__(self) -> None:
        super().__init__()
        self.nx = QSpinBox(minimum=1, maximum=64, value=12)
        self.ny = QSpinBox(minimum=1, maximum=64, value=12)
        self.lattice = QComboBox()
        self.lattice.addItems(["rectangular", "triangular"])
        self.taper = QComboBox()
        self.taper.addItems(sorted(TAPERS))
        self.taper.setCurrentText("chebyshev")
        self.sll = QDoubleSpinBox(minimum=-80.0, maximum=-5.0, value=-30.0, singleStep=5.0)
        self.sll.setSuffix(" dB")
        self.spacing = QDoubleSpinBox(minimum=0.05, maximum=3.0, value=0.5, singleStep=0.05)
        self.spacing.setSuffix(" λ")
        self.scan = QDoubleSpinBox(minimum=0.0, maximum=85.0, value=0.0, singleStep=5.0)
        self.scan.setSuffix("° from normal")
        self.scan_phi = QDoubleSpinBox(minimum=0.0, maximum=360.0, value=0.0, singleStep=15.0)
        self.scan_phi.setSuffix("° azimuth")
        self.ground = QCheckBox("ground-plane backed (one-sided)")
        self.ground.setChecked(True)

        controls = QFormLayout()
        controls.setVerticalSpacing(10)
        controls.addRow("Elements along x", self.nx)
        controls.addRow("Rows along y", self.ny)
        controls.addRow("Lattice", self.lattice)
        controls.addRow("Spacing", self.spacing)
        controls.addRow("Taper", self.taper)
        controls.addRow("Sidelobe level", self.sll)
        controls.addRow("Scan angle", self.scan)
        controls.addRow("Scan azimuth", self.scan_phi)
        controls.addRow("", self.ground)
        box = QGroupBox("Planar array")
        box.setLayout(controls)

        self.summary_table = _table()
        self.warning = QLabel("")
        self.warning.setWordWrap(True)
        self.warning.setStyleSheet(f"color: {WARN};")
        self.canvas = Canvas(polar=True)
        self.layout_canvas = Canvas()

        left = QVBoxLayout()
        left.setContentsMargins(0, 0, 0, 0)
        left.addWidget(box)
        results = QGroupBox("Result")
        rl = QVBoxLayout(results)
        rl.addWidget(self.summary_table)
        rl.addWidget(self.warning)
        left.addWidget(results, 1)
        left_widget = QWidget()
        left_widget.setLayout(left)
        left_widget.setMaximumWidth(440)

        right = _results_tabs((self.canvas, "Pattern"), (self.layout_canvas, "Layout"))

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(left_widget)
        splitter.addWidget(right)
        splitter.setSizes([360, 740])
        page, _ = _page("Planar arrays", "A rectangular or triangular lattice, tapered and steered: the two "
                        "principal cuts through the beam and the element layout.", splitter)
        QVBoxLayout(self).addWidget(page)
        self.layout().setContentsMargins(0, 0, 0, 0)

        for widget in (self.nx, self.ny, self.sll, self.spacing, self.scan,
                       self.scan_phi):
            widget.valueChanged.connect(self.refresh)
        for widget in (self.taper, self.lattice):
            widget.currentTextChanged.connect(self.refresh)
        self.ground.toggled.connect(self.refresh)
        self.refresh()

    # -- geometry and excitation --------------------------------------------

    def _line_taper(self, n: int):
        name = self.taper.currentText()
        if name in ("chebyshev", "taylor"):
            return TAPERS[name](n, self.sll.value())
        return TAPERS[name](n)

    def build(self):
        """Positions and weights for the current settings.

        A triangular lattice is not separable, so a taper cannot be applied to
        it the way it can to a rectangular grid. Saying so is better than
        quietly excising the control or, worse, applying a taper that does not
        mean what the label says.
        """
        nx, ny, d = self.nx.value(), self.ny.value(), self.spacing.value()
        if self.lattice.currentText() == "triangular":
            positions = triangular_lattice(nx, ny, d)
            return positions, np.ones(len(positions)), True
        positions = rectangular_lattice(nx, ny, d)
        return positions, separable_weights(self._line_taper(nx),
                                            self._line_taper(ny)), False

    # -- refresh -------------------------------------------------------------

    def refresh(self) -> None:
        lattice = self.lattice.currentText()
        tapered = lattice == "rectangular"
        self.taper.setEnabled(tapered)
        self.sll.setEnabled(tapered and self.taper.currentText() in ("chebyshev", "taylor"))
        try:
            positions, weights, uniform_forced = self.build()
            s = planar_summarise(positions, weights, self.scan.value(),
                                 self.scan_phi.value(),
                                 half_space=self.ground.isChecked())
        except Exception as exc:  # noqa: BLE001
            self.canvas.message(f"could not synthesise: {exc}")
            return

        limit = grating_lobe_free_spacing_planar(self.scan.value(), lattice)
        rows = [
            ("elements", str(s["elements"])),
            ("aperture",
             f"{s['aperture_x_lambda']:.2f} × {s['aperture_y_lambda']:.2f} λ"),
            ("directivity", _shown(s["directivity_dbi"], ".2f", " dBi")),
            ("beamwidth, scan plane", _shown(s["hpbw_scan_plane_deg"], ".3f", "°")),
            ("beamwidth, cross plane", _shown(s["hpbw_cross_plane_deg"], ".3f", "°")),
            ("sidelobe, scan plane", _shown(s["sidelobe_scan_plane_db"], ".2f", " dB")),
            ("sidelobe, cross plane", _shown(s["sidelobe_cross_plane_db"], ".2f", " dB")),
            ("taper efficiency", f"{s['taper_efficiency']:.4f}"),
            ("max spacing, no grating lobe", f"{limit:.4f} λ"),
        ]
        _fill(self.summary_table, rows)

        notes = []
        if visible_grating_lobes(lattice, self.spacing.value(), None, self.scan.value(),
                                 self.scan_phi.value()):
            notes.append(
                f"A grating lobe is in real space: spacing {self.spacing.value():.3f} λ "
                f"steered to {self.scan.value():.0f}° at φ {self.scan_phi.value():.0f}°.")
        elif self.spacing.value() >= limit:
            notes.append(
                f"No grating lobe at this φ, but spacing {self.spacing.value():.3f} λ exceeds the "
                f"{limit:.4f} λ limit for scanning to {self.scan.value():.0f}° in every plane.")
        if uniform_forced and self.taper.currentText() != "uniform":
            notes.append(
                "A triangular lattice is not separable, so the taper does not "
                "apply and the array is excited uniformly.")
        if lattice == "rectangular":
            notes.append(
                f"A triangular lattice would allow "
                f"{grating_lobe_free_spacing_planar(self.scan.value(), 'triangular'):.4f} λ "
                f"and cover the same aperture with "
                f"{lattice_element_saving() * 100:.1f}% fewer elements.")
        if not self.ground.isChecked():
            notes.append(
                "Isotropic elements radiate both ways, so half the power is in "
                "the mirror beam. Tick the ground-plane box for the one-sided "
                "figure.")
        self.warning.setText("  ".join(notes))

        cuts = []
        for label, plane in (("scan plane", "scan"), ("cross plane", "cross")):
            cut = planar_beam_cut(positions, weights, plane, self.scan.value(),
                                  self.scan_phi.value())
            cuts.append((label, cut.theta - math.pi / 2, cut.cut(0.0)))
        plot_hemisphere_cuts(
            self.canvas, cuts,
            f"{s['elements']} elements, {lattice}, "
            f"{'uniform' if uniform_forced else self.taper.currentText()} taper\n"
            f"D {s['directivity_dbi']:.1f} dBi, beam at {self.scan.value():.0f}° "
            f"from the normal")
        plot_element_layout(
            self.layout_canvas, positions, weights,
            f"{lattice} lattice, {self.spacing.value():.3f} λ")


# --------------------------------------------------------------- waveguides

class WaveguideTab(QWidget):
    """WR-series calculator with a dispersion plot."""

    def __init__(self) -> None:
        super().__init__()
        self.guide = QComboBox()
        self.guide.addItems(list(WR_SERIES))
        self.guide.setCurrentText("WR-90")
        self.freq = QDoubleSpinBox(minimum=0.01, maximum=1000.0, value=10.0, decimals=4)
        self.freq.setSuffix(" GHz")

        controls = QFormLayout()
        controls.setVerticalSpacing(10)
        controls.addRow("Guide", self.guide)
        controls.addRow("Frequency", self.freq)
        box = QGroupBox("Waveguide")
        box.setLayout(controls)

        self.props = _table()
        self.modes = _table()
        self.modes.setHorizontalHeaderLabels(["mode", "cutoff"])
        self.canvas = Canvas()

        left = QVBoxLayout()
        left.setContentsMargins(0, 0, 0, 0)
        left.addWidget(box)
        props_box = QGroupBox("Properties")
        QVBoxLayout(props_box).addWidget(self.props)
        modes_box = QGroupBox("Modes")
        QVBoxLayout(modes_box).addWidget(self.modes)
        left.addWidget(props_box, 3)
        left.addWidget(modes_box, 2)
        left_widget = QWidget()
        left_widget.setLayout(left)
        left_widget.setMaximumWidth(460)

        right = _results_tabs((self.canvas, "Dispersion"))

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(left_widget)
        splitter.addWidget(right)
        splitter.setSizes([400, 700])
        page, _ = _page("Waveguides", "Standard WR rectangular guides: cutoffs, the single-mode band, guide "
                        "wavelength, impedance and loss at your frequency.", splitter)
        QVBoxLayout(self).addWidget(page)
        self.layout().setContentsMargins(0, 0, 0, 0)

        self.guide.currentTextChanged.connect(self.refresh)
        self.freq.valueChanged.connect(self.refresh)
        self.refresh()

    def refresh(self) -> None:
        name = self.guide.currentText()
        g = standard(name)
        f0 = self.freq.value() * 1e9
        lo, hi = g.single_mode_band_hz
        rlo, rhi = recommended_band(name)

        rows = [
            ("dimensions", f"{g.a * 1e3:.3f} × {g.b * 1e3:.3f} mm"),
            ("aspect ratio", f"{g.aspect_ratio:.4f}"),
            ("TE10 cutoff", engineering(g.dominant_cutoff_hz, "Hz")),
            ("single-mode band", f"{lo / 1e9:.4f} – {hi / 1e9:.4f} GHz"),
            ("recommended band", f"{rlo / 1e9:.3f} – {rhi / 1e9:.3f} GHz"),
        ]
        if f0 <= g.dominant_cutoff_hz:
            rows.append(("at this frequency", "BELOW CUTOFF — evanescent"))
        else:
            if not g.is_single_mode(f0):
                rows.append(("at this frequency", "⚠ above single-mode band"))
            rows += [
                ("guide wavelength", engineering(g.guide_wavelength(f0), "m")),
                ("wave impedance", f"{g.wave_impedance(f0):.2f} Ω"),
                ("attenuation", f"{g.attenuation_db_per_m(f0):.4f} dB/m"),
                ("power capacity", f"{g.power_capacity_w(f0) / 1e6:.3f} MW"),
            ]
        _fill(self.props, rows)
        _fill(self.modes, [(m.label, f"{m.f_cutoff_hz / 1e9:.4f} GHz")
                           for m in g.modes(4 * g.dominant_cutoff_hz)])

        fc = g.dominant_cutoff_hz
        freqs = np.linspace(fc * 1.001, fc * 3.0, 400)
        lam_g = np.array([g.guide_wavelength(f) for f in freqs])
        lam_0 = 2.99792458e8 / freqs
        self.canvas.clear()
        ax = self.canvas.axes
        ax.set_axis_on()
        ax.plot(freqs / 1e9, lam_g * 1e3, color=ACCENT, lw=2.0, label="guide wavelength")
        ax.plot(freqs / 1e9, lam_0 * 1e3, "--", color=MUTED, lw=1.4, label="free space")
        ax.axvline(fc / 1e9, color="#a93a2e", ls=":", lw=1, label="TE10 cutoff")
        ax.axvspan(rlo / 1e9, rhi / 1e9, alpha=0.10, color="#2d7a55", lw=0, label="recommended")
        ax.axvline(f0 / 1e9, color=INK, ls=":", lw=0.8)
        ax.set_xlabel("frequency (GHz)")
        ax.set_ylabel("wavelength (mm)")
        ax.set_ylim(0, float(lam_0[0]) * 1e3 * 1.5)
        ax.legend(fontsize=8, frameon=False)
        self.canvas.draw_idle()


# ---------------------------------------------------------------------- app

class MainWindow(QMainWindow):
    SECTIONS = (("Catalogue", "Antennas", "antennas"), ("Linear arrays", "Linear arrays", "linear"),
                ("Planar arrays", "Planar arrays", "planar"), ("Waveguides", "Waveguides", "waveguide"))

    def __init__(self, registry: Registry | None = None) -> None:
        super().__init__()
        registry = registry or default_registry()
        self.setWindowTitle(
            f"OTA Hub Antenna Toolkit — {len(registry)} archetypes, "
            f"{len(registry.families)} families")
        self.resize(1360, 860)
        self.setMinimumSize(1060, 680)
        self.setStyleSheet(STYLE)

        # The pages live in a tab widget whose own tab bar is hidden: the rail drives it.
        self.tabs = QTabWidget()
        self.tabs.tabBar().hide()
        self.tabs.setDocumentMode(True)
        self.catalogue = CatalogueTab(registry)
        self.arrays = ArrayTab()
        self.planar = PlanarArrayTab()
        self.waveguides = WaveguideTab()
        for widget, (tab_name, _, _) in zip((self.catalogue, self.arrays, self.planar, self.waveguides), self.SECTIONS):
            self.tabs.addTab(widget, tab_name)

        self.rail = QListWidget()
        self.rail.setObjectName("rail")
        self.rail.setIconSize(QSize(20, 20))
        self.rail.setFixedWidth(196)
        for _, label, icon in self.SECTIONS:
            item = QListWidgetItem(rail_icon(icon), f"  {label}")
            item.setSizeHint(QSize(180, 44))
            self.rail.addItem(item)
        self.rail.currentRowChanged.connect(self._go)
        self.rail.setCurrentRow(0)

        side = QWidget()
        side.setStyleSheet("background: #1c2730;")
        sl = QVBoxLayout(side)
        sl.setContentsMargins(0, 0, 0, 0)
        sl.setSpacing(0)
        sl.addWidget(_label("OTA Hub", "brand"))
        sl.addWidget(_label("Antenna toolkit", "brandsub"))
        sl.addWidget(self.rail, 1)
        side.setFixedWidth(196)

        central = QWidget()
        cl = QHBoxLayout(central)
        cl.setContentsMargins(0, 0, 0, 0)
        cl.setSpacing(0)
        cl.addWidget(side)
        cl.addWidget(self.tabs, 1)
        self.setCentralWidget(central)

        problems = registry.problems()
        if problems:
            self.statusBar().showMessage(
                f"⚠ {sum(len(v) for v in problems.values())} spec problem(s) — run 'doctor'")
        else:
            self.statusBar().showMessage(
                f"{len(registry)} archetypes loaded, all structurally sound")

    def _go(self, row: int) -> None:
        self.catalogue._timer.stop()       # nothing recomputes behind another page
        self.tabs.setCurrentIndex(row)
        if row == 0:
            self.catalogue.show_gallery()


def main(argv: list[str] | None = None) -> int:
    # numbers read the same everywhere in the app: Python formats the tables with a
    # point, so the spin boxes must not take the region's decimal comma
    QLocale.setDefault(QLocale.c())
    app = QApplication(argv if argv is not None else sys.argv)
    app.setApplicationName("OTA Hub Antenna Toolkit")
    window = MainWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
