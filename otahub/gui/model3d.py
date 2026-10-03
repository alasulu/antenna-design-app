"""The design page's 3-D model: the solid the STL export writes, drawn before it
is written, with the construction options beside it.

The options are the choices a build needs that the electrical design does not
fix - a wall's thickness, a wire's radius, a board's margin. Each says whether
the predicted performance shown on the other tabs accounts for it, so changing
one shows its effect on the model and, honestly, on the figures (usually none:
"geometry only").
"""
from __future__ import annotations

import math

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QFileDialog, QFormLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QScrollArea,
                               QTextBrowser, QVBoxLayout, QWidget)

from .models import format_input, parse_quantity, shown_unit

COLOURS = {"PEC": "#c47f3f", "ferrite": "#555b62"}
DIELECTRIC = "#7fae8c"


def _colour(material: str) -> tuple[str, float]:
    if material == "PEC":
        return COLOURS["PEC"], 1.0
    if material.startswith("ferrite"):
        return COLOURS["ferrite"], 1.0
    return DIELECTRIC, 0.8


def draw_model(ax, model, empty_text: str = "No 3-D model") -> None:
    """Shade a model's materials into 3-D axes: dielectrics first and translucent,
    then ferrite, then the conductors, so the copper reads through the board."""
    from ..export import mesh
    ax.clear()
    ax.set_axis_off()
    if not model.built:
        ax.text2D(0.5, 0.5, empty_text, ha="center", va="center", transform=ax.transAxes)
        return
    light = np.array([0.4, -0.5, 0.75])
    light /= np.linalg.norm(light)
    # body by body from the bottom up, painted in that order: a ground under a board
    # under a patch reads correctly from above, which per-material depth sorting did not
    ax.computed_zorder = False
    bodies = [b for b in model.bodies if b.material not in ("VOID", "VACUUM")]
    bodies.sort(key=lambda b: (sum(b.solid.bounding_box()[i] for i in (2, 5)) / 2, b.material == "PEC"))
    for body in bodies:
        material, solid = body.material, body.solid
        tris = mesh._triangles(solid)
        n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
        n /= np.maximum(np.linalg.norm(n, axis=1), 1e-30)[:, None]
        colour, alpha = _colour(material)
        base = np.array([int(colour[i:i + 2], 16) / 255 for i in (1, 3, 5)])
        shade = 0.55 + 0.45 * np.abs(n @ light)
        faces = np.clip(base[None, :] * shade[:, None], 0, 1)
        ax.add_collection3d(Poly3DCollection(tris, facecolors=np.column_stack([faces, np.full(len(faces), alpha)]),
                                             edgecolor="none"))
    lo, hi = model.bounds()
    span = hi - lo
    big = float(span.max())
    shown = np.maximum(span, 0.12 * big)          # a flat board still gets some depth on screen
    centre = (lo + hi) / 2
    for setter, c, w in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), centre, shown):
        setter(c - w / 2, c + w / 2)
    ax.set_box_aspect(tuple(shown / big), zoom=1.3)
    ax.view_init(elev=38 if span[2] < 0.15 * big else 22, azim=-58)


class SolidPanel(QWidget):
    """Construction options, the solid model drawn from them, and Save STL."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._design = None
        self._key = None
        self._stale = True
        self._fields: dict[str, QLineEdit] = {}
        self._options = ()

        outer = QHBoxLayout(self)
        outer.setContentsMargins(4, 4, 4, 4)
        outer.setSpacing(12)

        left = QWidget()
        left.setFixedWidth(300)
        ll = QVBoxLayout(left)
        ll.setContentsMargins(0, 0, 0, 0)
        title = QLabel("Construction options")
        title.setObjectName("cardtitle")
        ll.addWidget(title)
        hint = QLabel("What a build needs that the electrical design does not fix. Hover an option to see "
                      "whether the predicted performance on the other tabs accounts for it (most are geometry "
                      "only). Values take units: 1 mm, 35 um, 32 mil.")
        hint.setWordWrap(True)
        hint.setObjectName("muted")
        ll.addWidget(hint)
        self.form_host = QWidget()
        self.form = QFormLayout(self.form_host)
        self.form.setContentsMargins(0, 4, 0, 4)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(self.form_host)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        ll.addWidget(scroll, 1)
        buttons = QHBoxLayout()
        self.reset_button = QPushButton("Defaults")
        self.reset_button.clicked.connect(self._reset)
        self.save_button = QPushButton("Save STL…")
        self.save_button.setObjectName("primary")
        self.save_button.clicked.connect(lambda: self.save_stl())
        self.enlarge_button = QPushButton("Enlarge")
        self.enlarge_button.clicked.connect(self.enlarge)
        buttons.addWidget(self.reset_button)
        buttons.addWidget(self.enlarge_button)
        buttons.addWidget(self.save_button, 1)
        ll.addLayout(buttons)
        self.status = QLabel("")
        self.status.setWordWrap(True)
        self.status.setObjectName("muted")
        ll.addWidget(self.status)
        outer.addWidget(left)

        right = QVBoxLayout()
        self.figure = Figure(figsize=(5, 4))
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.axes = self.figure.add_subplot(111, projection="3d")
        right.addWidget(self.canvas, 5)
        self.summary = QLabel("")
        self.summary.setWordWrap(True)
        right.addWidget(self.summary)
        self.notes = QTextBrowser()
        self.notes.setMaximumHeight(62)
        right.addWidget(self.notes, 1)
        self._model = None
        outer.addLayout(right, 1)

    # ---------------------------------------------------------------- data

    def set_design(self, design) -> None:
        """A new design: rebuild the form for a new archetype, redraw when visible."""
        from ..export import mesh
        self._design = design
        if not mesh.available():
            self.summary.setText("The 3-D model needs the manifold3d package: pip install manifold3d")
            self.save_button.setEnabled(False)
            return
        if design.archetype != self._key:
            self._key = design.archetype
            self._build_form(mesh.construction_options(design.archetype))
        self._stale = True
        if self.isVisible():
            self.refresh()

    def _build_form(self, options) -> None:
        while self.form.rowCount():
            self.form.removeRow(0)
        self._fields.clear()
        self._options = options
        for o in options:
            edit = QLineEdit()
            edit.setPlaceholderText("default")
            edit.setToolTip(f"{o.label}\n\n{o.effect}")
            edit.editingFinished.connect(self.refresh)
            label = QLabel(f"{o.label}" + (f"  [{shown_unit(o.unit)}]" if o.unit not in ("-", "") else ""))
            label.setWordWrap(True)
            label.setToolTip(o.effect)
            self.form.addRow(label, edit)
            self._fields[o.name] = edit

    def _reset(self) -> None:
        for e in self._fields.values():
            e.clear()
        self.refresh()

    def options(self):
        """The options the user typed; a field left blank takes its default. Bad
        entries are marked and left out."""
        from ..export import mesh
        values = {}
        for o in self._options:
            text = self._fields[o.name].text().strip()
            edit = self._fields[o.name]
            if not text:
                edit.setProperty("invalid", "false")
                continue
            try:
                value = parse_quantity(text, "" if o.unit == "-" else o.unit)
                values[o.name] = value
                bad = mesh.option_problems(self._key, {o.name: value})
                edit.setProperty("invalid", "true" if bad else "false")
                edit.setToolTip(f"{o.label}\n\n" + ("; ".join(bad) if bad else o.effect))
            except ValueError:
                edit.setProperty("invalid", "true")
            edit.style().unpolish(edit)
            edit.style().polish(edit)
        return mesh.Options(values)

    def _show_defaults(self, values: dict) -> None:
        for o in self._options:
            v = values.get(o.name)
            shown = ("design default" if v is None or not math.isfinite(v)
                     else format_input(v, "" if o.unit == "-" else o.unit))
            self._fields[o.name].setPlaceholderText(shown)

    # ---------------------------------------------------------------- drawing

    def showEvent(self, event) -> None:
        super().showEvent(event)
        if self._stale and self._design is not None:
            self.refresh()

    def refresh(self) -> None:
        from ..export import mesh
        if self._design is None or not mesh.available():
            return
        opts = self.options()
        self._stale = False
        try:
            self._show_defaults(mesh.option_values(self._design, opts))
            with mesh.preview():
                model = mesh.solid(self._design, opts)
        except Exception as exc:  # noqa: BLE001 - a model that cannot be built must not stop the page
            model = mesh.Solid3D(self._design.archetype, self._design.archetype,
                                 notes=[f"NO 3-D GEOMETRY: {type(exc).__name__}: {exc}"])
        self._draw(model)
        notes = "".join(f"<li>{n}</li>" for n in model.notes)
        self.notes.setHtml(f"<ul style='margin:0'>{notes}</ul>")
        self.save_button.setEnabled(model.built)

    def _draw(self, model) -> None:
        from ..export import mesh
        self._model = model
        draw_model(self.axes, model, "No 3-D model for this design:\nsee the notes below.")
        self.figure.subplots_adjust(0, 0, 1, 1)
        self.canvas.draw_idle()
        if not model.built:
            self.summary.setText("")
            return
        lo, hi = model.bounds()
        span = hi - lo
        mats = ", ".join(mesh.material_label(m) for m in model.materials())
        self.summary.setText(f"<b>{model.title}</b> · {span[0] * 1e3:.4g} × {span[1] * 1e3:.4g} × "
                             f"{span[2] * 1e3:.4g} mm · {mats} · drag to turn")

    def enlarge(self) -> None:
        """The model in a window of its own, at full detail."""
        from PySide6.QtWidgets import QDialog
        from ..export import mesh
        if self._design is None:
            return
        model = mesh.solid(self._design, self.options())
        dlg = QDialog(self)
        dlg.setWindowTitle(f"{model.title} - 3-D model (drag to turn, scroll to zoom)")
        dlg.resize(1000, 760)
        lay = QVBoxLayout(dlg)
        fig = Figure(figsize=(8, 6))
        canvas = FigureCanvasQTAgg(fig)
        ax = fig.add_subplot(111, projection="3d")
        draw_model(ax, model)
        fig.subplots_adjust(0, 0, 1, 1)
        lay.addWidget(canvas)
        self._dialog = dlg                      # kept, so it is not collected while open
        dlg.show()

    # ---------------------------------------------------------------- saving

    def save_stl(self, path: str | None = None) -> list:
        """Write the full-detail model (not the preview) to `path`, asking for one if not given."""
        from ..export import mesh
        if self._design is None:
            return []
        if path is None:
            path, _ = QFileDialog.getSaveFileName(self, "Save 3-D model", f"{self._design.archetype}.stl",
                                                  "STL files (*.stl)")
            if not path:
                return []
        model = mesh.solid(self._design, self.options())
        try:
            files = mesh.write_stl(model, path)
        except OSError as exc:
            self.status.setText(f"Nothing saved: {exc.strerror or exc} ({path}).")
            return []
        if files:
            self.status.setText("Saved, in millimetres: " + ", ".join(f.name for f in files))
        else:
            self.status.setText("Nothing saved: this design has no 3-D model (see the notes).")
        return files
