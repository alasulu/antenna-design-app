"""OTA Hub main window.

Four tabs mirroring the toolkit: the archetype catalogue, the linear array
synthesiser, the planar array designer and the waveguide calculator. The catalogue form is generated from
each spec's declared parameters, so a new archetype in ``specs/`` gets a
working UI with no code change here.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (QApplication, QCheckBox, QComboBox, QDoubleSpinBox,
                               QFormLayout, QGroupBox, QHBoxLayout, QHeaderView,
                               QLabel, QLineEdit, QMainWindow, QMessageBox,
                               QPlainTextEdit, QPushButton, QSpinBox, QSplitter,
                               QTableWidget, QTableWidgetItem, QTabWidget,
                               QTreeView, QVBoxLayout, QWidget)

from ..arrays import (TAPERS, grating_lobe_free_spacing_planar,
                      lattice_element_saving, planar_beam_cut,
                      planar_summarise, rectangular_lattice,
                      separable_weights, summarise, triangular_lattice)
from ..core import pattern as pat
from ..core.registry import Registry, default_registry
from ..core.units import engineering
from ..waveguides.rectangular import WR_SERIES, recommended_band, standard
from .models import KEY_ROLE, CatalogueFilter, build_catalogue_model, default_for, requirement_fields
from .plots import (Canvas, plot_element_layout, plot_hemisphere_cuts,
                    plot_polar, plot_sweep)

#: Archetypes whose far-field pattern we can compute from first principles.
#: Anything absent gets an honest message rather than a fabricated plot.
PATTERN_SOURCES = {
    "short_dipole": lambda: (pat.short_dipole(), "short dipole, sin^2"),
    "half_wave_dipole": lambda: (pat.finite_dipole(0.5), "half-wave dipole"),
    "resonant_dipole": lambda: (pat.finite_dipole(0.48), "resonant dipole, 0.48 lambda"),
    "quarter_wave_monopole": lambda: (pat.finite_dipole(0.5),
                                      "monopole (equivalent dipole with image)"),
    "folded_dipole": lambda: (pat.finite_dipole(0.5), "folded dipole"),
    "small_circular_loop": lambda: (pat.short_dipole(), "small loop (dual of short dipole)"),
}


def _table(rows: int = 0, cols: int = 2) -> QTableWidget:
    t = QTableWidget(rows, cols)
    t.setHorizontalHeaderLabels(["quantity", "value"])
    t.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
    t.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
    t.verticalHeader().setVisible(False)
    t.setAlternatingRowColors(True)
    t.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
    return t


def _fill(table: QTableWidget, pairs: list[tuple[str, str]]) -> None:
    table.setRowCount(len(pairs))
    for row, (name, value) in enumerate(pairs):
        table.setItem(row, 0, QTableWidgetItem(name))
        table.setItem(row, 1, QTableWidgetItem(value))


def _format(value, unit: str = "") -> str:
    if isinstance(value, complex):
        sign = "+" if value.imag >= 0 else "-"
        return f"{value.real:.4g} {sign} j{abs(value.imag):.4g} {unit}".strip()
    if isinstance(value, float):
        if unit in ("m", "Hz") and value != 0 and math.isfinite(value):
            return engineering(value, unit)
        return f"{value:.6g} {unit}".strip()
    return f"{value} {unit}".strip()


# ---------------------------------------------------------------- catalogue

class CatalogueTab(QWidget):
    """Browse archetypes, synthesise a design, sweep it, plot its pattern."""

    def __init__(self, registry: Registry) -> None:
        super().__init__()
        self.registry = registry
        self.current = None
        self._fields: dict[str, QLineEdit] = {}
        self._design = None

        self.search = QLineEdit(placeholderText="search archetypes…")
        self.tree = QTreeView()
        self.source_model = build_catalogue_model(registry)
        self.proxy = CatalogueFilter(registry, self)
        self.proxy.setSourceModel(self.source_model)
        self.tree.setModel(self.proxy)
        self.tree.setHeaderHidden(True)
        self.tree.expandAll()
        self.search.textChanged.connect(self._on_search)
        self.tree.selectionModel().currentChanged.connect(self._on_select)

        left = QVBoxLayout()
        left.addWidget(self.search)
        left.addWidget(self.tree)
        left_widget = QWidget()
        left_widget.setLayout(left)

        self.title = QLabel("select an archetype")
        font = QFont()
        font.setBold(True)
        font.setPointSize(13)
        self.title.setFont(font)
        self.summary = QLabel("")
        self.summary.setWordWrap(True)
        self.band = QLabel("")

        self.form_box = QGroupBox("requirements")
        self.form = QFormLayout(self.form_box)
        self.synth_button = QPushButton("Synthesise")
        self.synth_button.clicked.connect(self._synthesise)
        self.synth_button.setEnabled(False)

        middle = QVBoxLayout()
        for w in (self.title, self.summary, self.band, self.form_box, self.synth_button):
            middle.addWidget(w)
        middle.addStretch(1)
        middle_widget = QWidget()
        middle_widget.setLayout(middle)

        self.geometry_table = _table()
        self.metrics_table = _table()
        self.notes = QPlainTextEdit(readOnly=True)
        design = QWidget()
        dl = QVBoxLayout(design)
        dl.addWidget(QLabel("geometry"))
        dl.addWidget(self.geometry_table)
        dl.addWidget(QLabel("predicted performance"))
        dl.addWidget(self.metrics_table)
        dl.addWidget(QLabel("warnings and validity"))
        dl.addWidget(self.notes)

        self.metric_picker = QComboBox()
        self.metric_picker.currentTextChanged.connect(self._sweep)
        self.sweep_canvas = Canvas()
        sweep = QWidget()
        sl = QVBoxLayout(sweep)
        sl.addWidget(self.metric_picker)
        sl.addWidget(self.sweep_canvas)

        self.pattern_canvas = Canvas(polar=True)

        self.results = QTabWidget()
        self.results.addTab(design, "Design")
        self.results.addTab(sweep, "Sweep")
        self.results.addTab(self.pattern_canvas, "Pattern")

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(left_widget)
        splitter.addWidget(middle_widget)
        splitter.addWidget(self.results)
        splitter.setSizes([260, 320, 620])
        layout = QHBoxLayout(self)
        layout.addWidget(splitter)

    # -------------------------------------------------------------- actions

    def _on_search(self, text: str) -> None:
        self.proxy.set_needle(text)
        self.tree.expandAll()

    def _on_select(self, current, _previous) -> None:
        key = current.data(KEY_ROLE)
        if key is None:
            return
        self.select_key(key)

    def select_key(self, key: str) -> None:
        """Load an archetype into the form. Separated out so tests can drive it."""
        archetype = self.registry[key]
        self.current = archetype
        spec = archetype.spec
        self.title.setText(spec.name)
        self.summary.setText(spec.summary)
        lo, hi = spec.freq_range_hz
        band = f"valid {engineering(lo, 'Hz')} – {engineering(hi, 'Hz')}"
        if spec.confidence == "low":
            band += "     ⚠ low confidence: verify in a full-wave solver"
        self.band.setText(band)

        while self.form.rowCount():
            self.form.removeRow(0)
        self._fields.clear()
        for param in requirement_fields(archetype):
            edit = QLineEdit(default_for(param, archetype))
            hint = param.description or param.name
            edit.setPlaceholderText(f"{hint} [{param.unit}]" if param.unit else hint)
            edit.setToolTip(f"{param.role}: {hint}")
            edit.returnPressed.connect(self._synthesise)
            label = param.symbol + (" *" if param.role == "requirement" else "")
            self.form.addRow(label, edit)
            self._fields[param.symbol] = edit
        self.synth_button.setEnabled(True)
        self._clear_results()
        self._update_pattern()

    def _clear_results(self) -> None:
        _fill(self.geometry_table, [])
        _fill(self.metrics_table, [])
        self.notes.setPlainText("")
        self.metric_picker.clear()
        self.sweep_canvas.message("synthesise a design to sweep it")

    def collect(self) -> dict[str, float]:
        """Read the form. Blank fields are omitted, not zeroed."""
        out: dict[str, float] = {}
        for symbol, edit in self._fields.items():
            text = edit.text().strip()
            if not text:
                continue
            try:
                out[symbol] = float(text)
            except ValueError:
                raise ValueError(f"{symbol}: {text!r} is not a number") from None
        return out

    def _synthesise(self) -> None:
        if self.current is None:
            return
        try:
            values = self.collect()
        except ValueError as exc:
            QMessageBox.warning(self, "Invalid input", str(exc))
            return
        try:
            design = self.current.synthesize(**values)
        except Exception as exc:  # noqa: BLE001 - surfaced to the user
            QMessageBox.critical(self, "Synthesis failed", str(exc))
            return
        self._design = design
        units = design.units

        _fill(self.geometry_table, [
            (name, _format(value, units.get(name, "")))
            for name, value in design.parameters.items() if name not in ("lambda0", "k0")])
        _fill(self.metrics_table, [
            (name, _format(value, units.get(name, "")))
            for name, value in design.metrics.items()])

        lines: list[str] = []
        if design.warnings:
            lines += ["WARNINGS:"] + [f"  ! {w}" for w in design.warnings] + [""]
        missing = design.missing_requirements()
        if missing:
            lines += [f"Supply these to complete the design: {', '.join(missing)}", ""]
        if design.notes:
            lines += ["VALID ONLY UNDER:"] + [f"  - {n}" for n in design.notes]
        self.notes.setPlainText("\n".join(lines) or "no warnings")

        numeric = [k for k, v in design.metrics.items()
                   if isinstance(v, (int, float)) and not isinstance(v, bool)
                   and math.isfinite(float(v))]
        self.metric_picker.blockSignals(True)
        self.metric_picker.clear()
        self.metric_picker.addItems(sorted(numeric))
        self.metric_picker.blockSignals(False)
        if numeric:
            self._sweep(self.metric_picker.currentText())
        self._update_pattern()

    def _sweep(self, metric: str) -> None:
        """Re-synthesise across the validity band and plot one metric."""
        if not metric or self.current is None or self._design is None:
            return
        base = dict(self._design.requirements)
        if "f0" not in base:
            self.sweep_canvas.message(
                "this archetype is not parameterised by f0, so there is nothing to sweep")
            return
        lo, hi = self.current.spec.freq_range_hz
        centre = base["f0"]
        lo = max(lo, centre / 4.0)
        hi = min(hi if math.isfinite(hi) else centre * 4.0, centre * 4.0)
        if not hi > lo:
            self.sweep_canvas.message("validity band is too narrow to sweep")
            return
        freqs, values = [], []
        for f in np.geomspace(lo, hi, 160):
            trial = dict(base, f0=float(f))
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
        plot_sweep(self.sweep_canvas, freqs, values, metric,
                   self._design.units.get(metric, ""), self.current.spec.freq_range_hz)

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
        controls.addRow("elements", self.n)
        controls.addRow("taper", self.taper)
        controls.addRow("sidelobe level", self.sll)
        controls.addRow("spacing", self.spacing)
        controls.addRow("scan angle", self.scan)
        box = QGroupBox("array")
        box.setLayout(controls)

        self.summary_table = _table()
        self.canvas = Canvas(polar=True)
        self.weights_canvas = Canvas()

        left = QVBoxLayout()
        left.addWidget(box)
        left.addWidget(self.summary_table)
        left.addStretch(1)
        left_widget = QWidget()
        left_widget.setLayout(left)

        right = QTabWidget()
        right.addTab(self.canvas, "Pattern")
        right.addTab(self.weights_canvas, "Excitation")

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(left_widget)
        splitter.addWidget(right)
        splitter.setSizes([320, 720])
        QHBoxLayout(self).addWidget(splitter)

        for widget in (self.n, self.sll, self.spacing, self.scan):
            widget.valueChanged.connect(self.refresh)
        self.taper.currentTextChanged.connect(self.refresh)
        self.refresh()

    def weights(self):
        name = self.taper.currentText()
        if name in ("chebyshev", "taylor"):
            return TAPERS[name](self.n.value(), self.sll.value())
        if name == "cosine":
            return TAPERS[name](self.n.value())
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
            ("directivity", f"{s['directivity_dbi']:.2f} dBi"),
            ("half-power beamwidth", f"{s['hpbw_deg']:.3f}°"),
            ("first sidelobe", f"{s['sidelobe_db']:.2f} dB"),
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
        ax.stem(range(len(w)), w)
        ax.set_xlabel("element")
        ax.set_ylabel("normalised amplitude")
        ax.grid(True, alpha=0.3)
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
        controls.addRow("elements along x", self.nx)
        controls.addRow("rows along y", self.ny)
        controls.addRow("lattice", self.lattice)
        controls.addRow("spacing", self.spacing)
        controls.addRow("taper", self.taper)
        controls.addRow("sidelobe level", self.sll)
        controls.addRow("scan angle", self.scan)
        controls.addRow("scan azimuth", self.scan_phi)
        controls.addRow("", self.ground)
        box = QGroupBox("planar array")
        box.setLayout(controls)

        self.summary_table = _table()
        self.warning = QLabel("")
        self.warning.setWordWrap(True)
        self.warning.setStyleSheet("color: #b00020;")
        self.canvas = Canvas(polar=True)
        self.layout_canvas = Canvas()

        left = QVBoxLayout()
        left.addWidget(box)
        left.addWidget(self.summary_table)
        left.addWidget(self.warning)
        left.addStretch(1)
        left_widget = QWidget()
        left_widget.setLayout(left)

        right = QTabWidget()
        right.addTab(self.canvas, "Pattern")
        right.addTab(self.layout_canvas, "Layout")

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(left_widget)
        splitter.addWidget(right)
        splitter.setSizes([340, 700])
        QHBoxLayout(self).addWidget(splitter)

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
            ("directivity", f"{s['directivity_dbi']:.2f} dBi"),
            ("beamwidth, scan plane", f"{s['hpbw_scan_plane_deg']:.3f}°"),
            ("beamwidth, cross plane", f"{s['hpbw_cross_plane_deg']:.3f}°"),
            ("sidelobe, scan plane", f"{s['sidelobe_scan_plane_db']:.2f} dB"),
            ("sidelobe, cross plane", f"{s['sidelobe_cross_plane_db']:.2f} dB"),
            ("taper efficiency", f"{s['taper_efficiency']:.4f}"),
            ("max spacing, no grating lobe", f"{limit:.4f} λ"),
        ]
        _fill(self.summary_table, rows)

        notes = []
        if self.spacing.value() >= limit:
            notes.append(
                f"Spacing {self.spacing.value():.3f} λ exceeds the "
                f"{limit:.4f} λ grating-lobe limit for scanning to "
                f"{self.scan.value():.0f}°, so a grating lobe is in real space.")
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
        controls.addRow("guide", self.guide)
        controls.addRow("frequency", self.freq)
        box = QGroupBox("waveguide")
        box.setLayout(controls)

        self.props = _table()
        self.modes = _table()
        self.modes.setHorizontalHeaderLabels(["mode", "cutoff"])
        self.canvas = Canvas()

        left = QVBoxLayout()
        left.addWidget(box)
        left.addWidget(QLabel("properties"))
        left.addWidget(self.props)
        left.addWidget(QLabel("modes"))
        left.addWidget(self.modes)
        left_widget = QWidget()
        left_widget.setLayout(left)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(left_widget)
        splitter.addWidget(self.canvas)
        splitter.setSizes([400, 640])
        QHBoxLayout(self).addWidget(splitter)

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
        ax.plot(freqs / 1e9, lam_g * 1e3, label="guide wavelength")
        ax.plot(freqs / 1e9, lam_0 * 1e3, "--", label="free space")
        ax.axvline(fc / 1e9, color="r", ls=":", lw=1, label="TE10 cutoff")
        ax.axvspan(rlo / 1e9, rhi / 1e9, alpha=0.08, color="green", label="recommended")
        ax.set_xlabel("frequency (GHz)")
        ax.set_ylabel("wavelength (mm)")
        ax.set_ylim(0, float(lam_0[0]) * 1e3 * 1.5)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        self.canvas.draw_idle()


# ---------------------------------------------------------------------- app

class MainWindow(QMainWindow):
    def __init__(self, registry: Registry | None = None) -> None:
        super().__init__()
        registry = registry or default_registry()
        self.setWindowTitle(
            f"OTA Hub Antenna Toolkit — {len(registry)} archetypes, "
            f"{len(registry.families)} families")
        self.resize(1280, 800)
        self.tabs = QTabWidget()
        self.catalogue = CatalogueTab(registry)
        self.arrays = ArrayTab()
        self.planar = PlanarArrayTab()
        self.waveguides = WaveguideTab()
        self.tabs.addTab(self.catalogue, "Catalogue")
        self.tabs.addTab(self.arrays, "Linear arrays")
        self.tabs.addTab(self.planar, "Planar arrays")
        self.tabs.addTab(self.waveguides, "Waveguides")
        self.setCentralWidget(self.tabs)
        problems = registry.problems()
        if problems:
            self.statusBar().showMessage(
                f"⚠ {sum(len(v) for v in problems.values())} spec problem(s) — run 'doctor'")
        else:
            self.statusBar().showMessage(
                f"{len(registry)} archetypes loaded, all structurally sound")


def main(argv: list[str] | None = None) -> int:
    app = QApplication(argv if argv is not None else sys.argv)
    app.setApplicationName("OTA Hub Antenna Toolkit")
    window = MainWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
