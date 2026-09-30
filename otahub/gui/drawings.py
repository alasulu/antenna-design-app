"""Technical drawings of every archetype, painted from the design's own geometry.

Each drawing is a small datasheet figure: conductors in copper, substrates and
dielectrics in their own colours, metal sheets in steel, the feed marked in red,
and the dimensions that matter labelled with the values the synthesis produced.
Proportions follow the computed geometry wherever the spec gives it; a quantity
the spec does not produce is drawn at a nominal proportion and left unlabelled,
so no label ever shows a number the design did not compute.

`paint(painter, rect, key, family, values, labels=True)` draws one figure;
`thumbnail(...)` renders the unlabelled version for the gallery.
"""
from __future__ import annotations

import math
from typing import Callable, Mapping

from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import (QBrush, QColor, QFont, QImage, QLinearGradient, QPainter,
                           QPainterPath, QPen, QPixmap, QPolygonF, QRadialGradient)

from ..core.units import engineering

# ---------------------------------------------------------------- materials
COPPER = QColor("#c47f3f")
COPPER_EDGE = QColor("#8d5322")
COPPER_LIGHT = QColor("#e1ad78")
SUBSTRATE = QColor("#d3e2cc")
SUBSTRATE_EDGE = QColor("#7f9c76")
CERAMIC = QColor("#d9e6ef")
CERAMIC_EDGE = QColor("#6f8ba0")
METAL = QColor("#b4bec6")
METAL_DARK = QColor("#7d8a95")
METAL_LIGHT = QColor("#dde3e8")
FERRITE = QColor("#44494f")
GLASS = QColor("#e4eef2")
FEED = QColor("#d1453a")
INK = QColor("#2b3640")
DIM = QColor("#5b6b77")
GRID = QColor("#e9eef1")
RAY = QColor("#e39a55")

FAMILY_COLOURS = {
    "wire": "#b8672e", "loop": "#b98414", "patch": "#3e8a5c", "horn": "#3b6c9c",
    "travelling_wave": "#2b8585", "uwb": "#7457a6", "reflector": "#5b6a79",
    "slot": "#4a5bb0", "lens": "#2e87ad", "dielectric": "#b0526a",
}


def family_colour(family: str) -> QColor:
    return QColor(FAMILY_COLOURS.get(family, "#5b6a79"))


class _Nominal(float):
    """A proportion the design did not compute. Arithmetic keeps it nominal, and no
    label prints it (a dipole without a length was labelled "L = 1 m")."""


def _nominal_op(name):
    def op(self, *args):
        out = getattr(float, name)(self, *args)
        return _Nominal(out) if isinstance(out, float) else out
    return op


for _name in ("__add__", "__radd__", "__sub__", "__rsub__", "__mul__", "__rmul__", "__truediv__",
              "__rtruediv__", "__pow__", "__rpow__", "__neg__", "__pos__", "__abs__"):
    setattr(_Nominal, _name, _nominal_op(_name))

#: in a label, marks a number the design did not compute: Fig drops such labels
_UNSET = "\x00"


def _num(values: Mapping, *names, default=None):
    for n in names:
        v = values.get(n)
        if isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) and v > 0:
            return float(v)
    return _Nominal(default) if isinstance(default, (int, float)) else default


def _len(v) -> str:
    if isinstance(v, _Nominal):
        return _UNSET
    return engineering(v, "m", 3) if v else ""


def _count(v) -> str:
    """A whole count for a label, or the unset marker."""
    return _UNSET if isinstance(v, _Nominal) or v is None else f"{int(round(v))}"


# ---------------------------------------------------------------- the figure: a world box mapped onto a rect
class Fig:
    def __init__(self, p: QPainter, rect: QRectF, labels: bool) -> None:
        self.p, self.r, self.labels = p, rect, labels
        self.drawn: list[str] = []                 # every label painted, for the tests
        self.u = max(0.7, min(2.2, rect.height() / 300.0))
        self.s, self.ox, self.oy = 1.0, 0.0, 0.0
        size = max(8.5, min(13.0, rect.height() / 26.0))
        self.font = QFont()
        self.font.setPointSizeF(size * 0.78)
        self.font.setStyleHint(QFont.StyleHint.Monospace)
        self.font.setFamilies(["SF Mono", "Menlo", "Consolas", "DejaVu Sans Mono", "monospace"])

    def fit(self, x0: float, y0: float, x1: float, y1: float, margin: float = 0.12) -> None:
        w, h = max(x1 - x0, 1e-12), max(y1 - y0, 1e-12)
        m = margin if self.labels else 0.08
        aw, ah = self.r.width() * (1 - 2 * m), self.r.height() * (1 - 2 * m)
        self.s = min(aw / w, ah / h)
        c = self.r.center()
        self.ox = c.x() - self.s * (x0 + x1) / 2
        self.oy = c.y() + self.s * (y0 + y1) / 2

    def pt(self, x: float, y: float) -> QPointF:
        return QPointF(self.ox + self.s * x, self.oy - self.s * y)

    def minw(self, px: float) -> float:
        return px / self.s

    # -------------------------------------------------------------- primitives
    def pen(self, colour: QColor, width: float = 1.0, style=Qt.PenStyle.SolidLine) -> QPen:
        pen = QPen(colour, width * self.u, style)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
        return pen

    def poly(self, pts, fill=None, edge: QColor | None = INK, width: float = 1.2, closed: bool = True) -> None:
        path = QPainterPath()
        q = [self.pt(x, y) for x, y in pts]
        path.moveTo(q[0])
        for a in q[1:]:
            path.lineTo(a)
        if closed:
            path.closeSubpath()
        self.p.setBrush(QBrush(fill) if isinstance(fill, QColor) else (fill or Qt.BrushStyle.NoBrush))
        self.p.setPen(self.pen(edge, width) if edge else Qt.PenStyle.NoPen)
        self.p.drawPath(path)

    def box(self, x, y, w, h, fill=None, edge: QColor | None = INK, width: float = 1.2) -> None:
        self.poly([(x, y), (x + w, y), (x + w, y + h), (x, y + h)], fill, edge, width)

    def circle(self, cx, cy, r, fill=None, edge: QColor | None = INK, width: float = 1.2) -> None:
        self.p.setBrush(QBrush(fill) if isinstance(fill, QColor) else (fill or Qt.BrushStyle.NoBrush))
        self.p.setPen(self.pen(edge, width) if edge else Qt.PenStyle.NoPen)
        c = self.pt(cx, cy)
        self.p.drawEllipse(c, r * self.s, r * self.s)

    def ellipse(self, cx, cy, rx, ry, fill=None, edge: QColor | None = INK, width: float = 1.2) -> None:
        self.p.setBrush(QBrush(fill) if isinstance(fill, QColor) else (fill or Qt.BrushStyle.NoBrush))
        self.p.setPen(self.pen(edge, width) if edge else Qt.PenStyle.NoPen)
        self.p.drawEllipse(self.pt(cx, cy), rx * self.s, ry * self.s)

    def line(self, x1, y1, x2, y2, colour: QColor = INK, width: float = 1.2, style=Qt.PenStyle.SolidLine) -> None:
        self.p.setPen(self.pen(colour, width, style))
        self.p.drawLine(self.pt(x1, y1), self.pt(x2, y2))

    def curve(self, pts, colour: QColor = COPPER, width: float = 2.6) -> None:
        path = QPainterPath()
        q = [self.pt(x, y) for x, y in pts]
        path.moveTo(q[0])
        for a in q[1:]:
            path.lineTo(a)
        self.p.setBrush(Qt.BrushStyle.NoBrush)
        self.p.setPen(self.pen(colour, width))
        self.p.drawPath(path)

    def wire(self, x1, y1, x2, y2, width: float = 3.0) -> None:
        self.line(x1, y1, x2, y2, COPPER_EDGE, width + 1.2)
        self.line(x1, y1, x2, y2, COPPER, width)

    def gradient(self, x0, y0, x1, y1, a: QColor, b: QColor) -> QBrush:
        g = QLinearGradient(self.pt(x0, y0), self.pt(x1, y1))
        g.setColorAt(0.0, a)
        g.setColorAt(1.0, b)
        return QBrush(g)

    def feed(self, x, y, r_px: float = 4.2) -> None:
        c = self.pt(x, y)
        self.p.setPen(QPen(QColor("#ffffff"), 1.6 * self.u))
        self.p.setBrush(QBrush(FEED))
        self.p.drawEllipse(c, r_px * self.u, r_px * self.u)

    def ground(self, x0, x1, y, depth_px: float = 7.0) -> None:
        """A ground plane seen edge-on: a steel bar with hatching below."""
        a, b = self.pt(x0, y), self.pt(x1, y)
        d = depth_px * self.u
        self.p.setPen(Qt.PenStyle.NoPen)
        self.p.setBrush(QBrush(METAL))
        self.p.drawRect(QRectF(a.x(), a.y(), b.x() - a.x(), d * 0.6))
        self.p.setPen(QPen(METAL_DARK, 1.0 * self.u))
        self.p.drawLine(a, b)
        step = 6 * self.u
        x = a.x()
        while x < b.x():
            self.p.drawLine(QPointF(x, a.y() + d * 0.6), QPointF(x - d * 0.6, a.y() + d * 1.3))
            x += step

    def text(self, x, y, s: str, align: str = "c", colour: QColor = DIM, dx: float = 0, dy: float = 0) -> None:
        if not self.labels or not s or _UNSET in s:
            return
        self.drawn.append(s)
        self.p.setFont(self.font)
        self.p.setPen(QPen(colour))
        c = self.pt(x, y)
        fm = self.p.fontMetrics()
        w = fm.horizontalAdvance(s)
        h = fm.ascent()
        if align == "c":
            px, py = c.x() - w / 2, c.y() + h / 2
        elif align == "l":
            px, py = c.x(), c.y() + h / 2
        else:
            px, py = c.x() - w, c.y() + h / 2
        self.p.drawText(QPointF(px + dx * self.u, py + dy * self.u), s)

    def _arrow(self, tip: QPointF, towards: QPointF) -> None:
        v = QPointF(towards.x() - tip.x(), towards.y() - tip.y())
        n = math.hypot(v.x(), v.y()) or 1.0
        ux, uy = v.x() / n, v.y() / n
        s = 6 * self.u
        a = QPointF(tip.x() + ux * s - uy * s * 0.4, tip.y() + uy * s + ux * s * 0.4)
        b = QPointF(tip.x() + ux * s + uy * s * 0.4, tip.y() + uy * s - ux * s * 0.4)
        self.p.setBrush(QBrush(DIM))
        self.p.setPen(Qt.PenStyle.NoPen)
        self.p.drawPolygon(QPolygonF([tip, a, b]))

    def dim(self, x1, y1, x2, y2, label: str, offset_px: float = 0.0, side: int = 1) -> None:
        """A dimension line from (x1, y1) to (x2, y2), pushed off the object by offset_px."""
        if not self.labels or not label or _UNSET in label:
            return
        self.drawn.append(label)
        a, b = self.pt(x1, y1), self.pt(x2, y2)
        vx, vy = b.x() - a.x(), b.y() - a.y()
        n = math.hypot(vx, vy) or 1.0
        nx, ny = vy / n * side, -vx / n * side
        o = offset_px * self.u
        a2, b2 = QPointF(a.x() + nx * o, a.y() + ny * o), QPointF(b.x() + nx * o, b.y() + ny * o)
        if o:
            ext = 4 * self.u
            self.p.setPen(QPen(QColor(DIM.red(), DIM.green(), DIM.blue(), 150), 0.8 * self.u))
            self.p.drawLine(a, QPointF(a2.x() + nx * ext, a2.y() + ny * ext))
            self.p.drawLine(b, QPointF(b2.x() + nx * ext, b2.y() + ny * ext))
        self.p.setPen(QPen(DIM, 0.9 * self.u))
        self.p.drawLine(a2, b2)
        self._arrow(a2, b2)
        self._arrow(b2, a2)
        self.p.setFont(self.font)
        fm = self.p.fontMetrics()
        w, h = fm.horizontalAdvance(label), fm.ascent()
        mx, my = (a2.x() + b2.x()) / 2, (a2.y() + b2.y()) / 2
        pad = 3 * self.u
        if abs(vx) >= abs(vy):                        # horizontal: text on the side the line was pushed to
            tx, ty = mx - w / 2, (my - pad) if ny <= 0 else (my + h + pad)
        else:                                          # vertical: text beside it
            tx, ty = (mx + pad + 2 * self.u) if nx >= 0 else (mx - w - pad - 2 * self.u), my + h / 2
        bg = QRectF(tx - 2 * self.u, ty - h - 1 * self.u, w + 4 * self.u, h + 4 * self.u)
        self.p.setPen(Qt.PenStyle.NoPen)
        self.p.setBrush(QBrush(QColor(255, 255, 255, 215)))
        self.p.drawRoundedRect(bg, 2 * self.u, 2 * self.u)
        self.p.setPen(QPen(INK))
        self.p.drawText(QPointF(tx, ty), label)


# ---------------------------------------------------------------- wire family
def _dipole(f: Fig, v, name="L"):
    L = _num(v, "L", default=1.0)
    g = L * 0.035
    f.fit(-L * 0.45, -L * 0.56, L * 0.45, L * 0.56)
    f.wire(0, g, 0, L / 2, 3.2)
    f.wire(0, -g, 0, -L / 2, 3.2)
    f.line(0, g, 0, -g, FEED, 1.2)
    f.feed(0, 0)
    f.dim(0, -L / 2, 0, L / 2, f"{name} = {_len(L)}", 30, side=1)
    f.text(0, 0, "feed", "l", dx=10)


def _folded_dipole(f: Fig, v):
    L = _num(v, "L", default=1.0)
    d = max(_num(v, "d_sep", default=L * 0.04), L * 0.08)
    f.fit(-L * 0.5, -L * 0.2, L * 0.72, d + L * 0.26)
    r = d / 2
    path = [(-L / 2 + r, d)]
    path += [(L / 2 - r, d)]
    path += [(L / 2 - r + r * math.sin(t), r + r * math.cos(t)) for t in [i * math.pi / 16 for i in range(17)]]
    path += [(0.03 * L, 0)]
    f.curve(path, COPPER_EDGE, 4.2)
    f.curve(path, COPPER, 3.0)
    path2 = [(-0.03 * L, 0), (-L / 2 + r, 0)]
    path2 += [(-L / 2 + r - r * math.sin(t), r - r * math.cos(t)) for t in [i * math.pi / 16 for i in range(17)]]
    f.curve(path2, COPPER_EDGE, 4.2)
    f.curve(path2, COPPER, 3.0)
    f.feed(0, 0)
    f.dim(-L / 2, d, L / 2, d, f"L = {_len(L)}", 26, side=1)
    sep = _num(v, "d_sep")
    if sep:
        f.dim(L / 2, 0, L / 2, d, f"d = {_len(sep)}", 18, side=-1)


def _dipole_over_ground(f: Fig, v):
    L = _num(v, "L", default=1.0)
    h = _num(v, "h", default=L / 2)
    f.fit(-L * 0.7, -h * 0.25, L * 0.7, h * 1.35)
    f.ground(-L * 0.65, L * 0.65, 0)
    g = L * 0.03
    f.wire(-L / 2, h, -g, h)
    f.wire(g, h, L / 2, h)
    f.feed(0, h)
    f.dim(-L / 2, h, L / 2, h, f"L = {_len(L)}", 18, side=1)
    f.dim(L * 0.58, 0, L * 0.58, h, f"h = {_len(h)}", 0, side=-1)


def _monopole(f: Fig, v, kind="plain"):
    h = _num(v, "h", default=1.0)
    f.fit(-h * 0.8, -h * 0.18, h * 0.8, h * 1.12)
    f.ground(-h * 0.75, h * 0.75, 0)
    base = h * 0.03
    if kind == "loaded":
        f.wire(0, base, 0, h * 0.42)
        # the loading coil
        n, r, y0, y1 = 7, h * 0.05, h * 0.42, h * 0.6
        pts = [(r * math.sin(2 * math.pi * n * t / 60), y0 + (y1 - y0) * t / 60) for t in range(61)]
        f.curve(pts, COPPER_EDGE, 3.4)
        f.curve(pts, COPPER, 2.3)
        f.wire(0, h * 0.6, 0, h)
        f.text(-r, (y0 + y1) / 2, "loading coil", "r", dx=-10)
    else:
        f.wire(0, base, 0, h)
    if kind == "hat":
        a = _num(v, "a_hat", default=h * 0.2)
        f.box(-a, h, 2 * a, f.minw(4.5), COPPER, COPPER_EDGE, 1.0)
        f.dim(-a, h, a, h, f"hat ⌀ {_len(2 * a)}" if _num(v, "a_hat") else "", 14, side=1)
    f.feed(0, base)
    f.dim(h * 0.55, 0, h * 0.55, h, f"h = {_len(h)}", 0, side=-1)


def _biconical(f: Fig, v):
    Lc = _num(v, "Lc", default=1.0)
    th = min(_num(v, "theta_h", default=math.radians(30)), math.radians(85))   # radians, from the axis
    rh, hh = Lc * math.sin(th), Lc * math.cos(th)
    gap = _num(v, "feed_gap", default=hh * 0.08) / 2
    f.fit(-rh * 1.4, -hh * 1.15, rh * 1.7, hh * 1.15)
    for sgn in (1, -1):
        f.poly([(0, sgn * gap), (-rh, sgn * (gap + hh)), (rh, sgn * (gap + hh))],
               f.gradient(-rh, 0, rh, 0, COPPER_LIGHT, COPPER), COPPER_EDGE, 1.2)
        f.ellipse(0, sgn * (gap + hh), rh, rh * 0.18, COPPER_LIGHT, COPPER_EDGE, 1.0)
    f.feed(0, 0)
    f.dim(0, gap, rh, gap + hh, f"slant {_len(Lc)}", 22, side=-1)


def _turnstile(f: Fig, v):
    L = _num(v, "L", default=1.0)
    g = L * 0.04
    f.fit(-L * 0.62, -L * 0.62, L * 0.62, L * 0.62)
    for a in (0, 90):
        c, s = math.cos(math.radians(a)), math.sin(math.radians(a))
        f.wire(g * c, g * s, L / 2 * c, L / 2 * s)
        f.wire(-g * c, -g * s, -L / 2 * c, -L / 2 * s)
    path = [(g * 2.5 * math.cos(t), g * 2.5 * math.sin(t)) for t in [math.pi * 0.1 * i / 4 + math.pi * 0.05 for i in range(21)]]
    f.curve(path, FEED, 1.4)
    f.feed(0, 0)
    f.dim(-L / 2, 0, L / 2, 0, f"L = {_len(L)}", 18, side=1)
    f.text(g * 3, -g * 3, "90° phasing line", "l", dx=10, dy=14)


# ---------------------------------------------------------------- loops
def _circle_loop(f: Fig, v, kind="plain"):
    C = _num(v, "C", "circumference_m", default=None)
    D = _num(v, "Dm", default=C / math.pi if C else 1.0)
    r = D / 2
    f.fit(-r * 1.35, -r * 1.3, r * 1.5, r * 1.3)
    gap = 0.08 if kind != "halo" else 0.16
    t0, t1 = -math.pi / 2 + gap, 3 * math.pi / 2 - gap
    pts = [(r * math.cos(t0 + (t1 - t0) * i / 90), r * math.sin(t0 + (t1 - t0) * i / 90)) for i in range(91)]
    f.curve(pts, COPPER_EDGE, 4.2)
    f.curve(pts, COPPER, 3.0)
    if kind == "halo":
        x = r * math.cos(t0)
        f.line(-x, -r - f.minw(6), -x, -r + f.minw(6), INK, 2.0)
        f.line(x, -r - f.minw(6), x, -r + f.minw(6), INK, 2.0)
        f.text(0, -r, "gap capacitor", "c", dy=18)
        f.feed(-r * 0.7, -r * 0.7)
    else:
        f.feed(0, -r)
    f.dim(-r, 0, r, 0, f"⌀ {_len(D)}", 0)


def _square_loop(f: Fig, v, kind="plain"):
    s = _num(v, "s", default=1.0)
    f.fit(-s * 0.72, -s * 0.72, s * 0.85, s * 0.72)
    g = s * 0.05
    pts = [(-g, -s / 2), (-s / 2, -s / 2), (-s / 2, s / 2), (s / 2, s / 2), (s / 2, -s / 2), (g, -s / 2)]
    f.curve(pts, COPPER_EDGE, 4.2)
    f.curve(pts, COPPER, 3.0)
    if kind == "alford":
        for cx, cy in ((-s / 2, s / 2), (s / 2, s / 2), (-s / 2, -s / 2), (s / 2, -s / 2)):
            f.box(cx - f.minw(5), cy - f.minw(5), f.minw(10), f.minw(10), QColor("#ffffff"), INK, 1.2)
        f.text(-s / 2, -s / 2, "corner loads", "l", dx=12, dy=16)
    f.feed(0, -s / 2)
    f.dim(-s / 2, s / 2, s / 2, s / 2, f"side {_len(s)}", 16, side=1)


def _multiturn(f: Fig, v):
    a = _num(v, "a", default=0.1)
    N = int(_num(v, "N", default=6))
    lc = _num(v, "l_coil", default=a * 1.2)
    n_draw = max(3, min(N, 14))
    f.fit(-lc * 0.8, -a * 1.5, lc * 1.6, a * 1.5)
    for i in range(n_draw):
        x = -lc / 2 + lc * (i + 0.5) / n_draw
        f.ellipse(x, 0, lc / n_draw * 0.28, a, None, COPPER_EDGE, 3.4)
        f.ellipse(x, 0, lc / n_draw * 0.28, a, None, COPPER, 2.2)
    f.feed(-lc / 2, -a)
    f.dim(-lc / 2, a, lc / 2, a, f"{N} turns" if _num(v, "N") else "", 16, side=1)
    f.dim(lc / 2, -a, lc / 2, a, f"⌀ {_len(2 * a)}" if _num(v, "a") else "", 16, side=-1)


def _ferrite(f: Fig, v):
    l = _num(v, "l_rod", default=0.1)
    d = _num(v, "d_rod", default=l * 0.08)
    d = max(d, l * 0.07)
    f.fit(-l * 0.6, -l * 0.25, l * 0.6, l * 0.25)
    f.box(-l / 2, -d / 2, l, d, f.gradient(0, d / 2, 0, -d / 2, QColor("#6a7077"), FERRITE), QColor("#2b2f33"), 1.2)
    lc = l * 0.35
    n = 14
    for i in range(n):
        x = -lc / 2 + lc * (i + 0.5) / n
        f.line(x - lc / n * 0.3, -d * 0.62, x + lc / n * 0.3, d * 0.62, COPPER, 2.2)
    f.feed(-lc / 2, -d * 0.62)
    f.dim(-l / 2, -d / 2, l / 2, -d / 2, f"rod {_len(_num(v, 'l_rod'))}" if _num(v, "l_rod") else "", 18, side=-1)
    f.text(0, d / 2, "ferrite rod · coil", "c", dy=-16)


# ---------------------------------------------------------------- patches
def _patch_frame(f: Fig, W, L, extra=1.55):
    """Top view of a patch on its board, with an edge-on strip underneath."""
    sw, sl = W * extra, L * extra
    f.fit(-sw * 0.66, -sl * 0.85, sw * 0.9, sl * 0.62)
    f.box(-sw / 2, -sl / 2, sw, sl, SUBSTRATE, SUBSTRATE_EDGE, 1.0)
    return sw, sl


def _side_strip(f: Fig, v, W, sw, sl, patches=((None, None),)):
    """The edge-on view under the top view: ground, substrate of thickness h, patch."""
    h = _num(v, "h")
    y0 = -sl / 2 - sl * 0.2
    t = max(sl * 0.05, f.minw(6))
    f.box(-sw / 2, y0 - t, sw, t, SUBSTRATE, SUBSTRATE_EDGE, 1.0)
    f.box(-sw / 2, y0 - t - f.minw(3), sw, f.minw(3), METAL, METAL_DARK, 0.8)
    f.box(-W / 2, y0, W, f.minw(3), COPPER, COPPER_EDGE, 0.8)
    if h:
        f.dim(sw / 2, y0 - t, sw / 2, y0, f"h = {_len(h)}", 14, side=-1)
    f.text(-sw / 2, y0 - t / 2, "side", "r", dx=-6)


def _rect_patch(f: Fig, v, kind="plain"):
    W = _num(v, "W", default=1.0)
    L = _num(v, "L", default=W * 0.75)
    if kind == "cp":
        W = L
    sw, sl = _patch_frame(f, W, L)
    c = _num(v, "c_trunc", default=L * 0.12) if kind == "cp" else 0
    if kind == "cp":
        pts = [(-W / 2 + c, L / 2), (W / 2, L / 2), (W / 2, -L / 2 + c), (W / 2 - c, -L / 2), (-W / 2, -L / 2), (-W / 2, L / 2 - c)]
        f.poly(pts, f.gradient(-W / 2, L / 2, W / 2, -L / 2, COPPER_LIGHT, COPPER), COPPER_EDGE, 1.2)
        f.feed(0, -L * 0.22)
        # c_trunc is the cut's LEG along each edge, so dimension a leg, not the diagonal
        f.dim(W / 2 - c, -L / 2, W / 2, -L / 2, f"cut {_len(c)}", 14, side=-1)
    elif kind == "inset":
        y0 = _num(v, "y0", default=L * 0.3)
        wn = W * 0.08
        wf = W * 0.1
        pts = [(-W / 2, -L / 2), (-wf / 2 - wn, -L / 2), (-wf / 2 - wn, -L / 2 + y0), (wf / 2 + wn, -L / 2 + y0),
               (wf / 2 + wn, -L / 2), (W / 2, -L / 2), (W / 2, L / 2), (-W / 2, L / 2)]
        f.poly(pts, f.gradient(-W / 2, L / 2, W / 2, -L / 2, COPPER_LIGHT, COPPER), COPPER_EDGE, 1.2)
        f.box(-wf / 2, -sl / 2, wf, sl / 2 - L / 2 + y0, COPPER, COPPER_EDGE, 1.0)
        f.feed(0, -sl / 2)
        if _num(v, "y0"):
            f.dim(W / 2, -L / 2, W / 2, -L / 2 + y0, f"inset {_len(y0)}", 14, side=-1)
    else:
        f.box(-W / 2, -L / 2, W, L, f.gradient(-W / 2, L / 2, W / 2, -L / 2, COPPER_LIGHT, COPPER), COPPER_EDGE, 1.2)
        if kind == "shorted":
            f.box(-W / 2, L / 2 - f.minw(6), W, f.minw(6), METAL_DARK, INK, 0.8)
            for i in range(7):
                f.circle(-W / 2 + W * (i + 0.5) / 7, L / 2 - f.minw(3), f.minw(1.6), QColor("#ffffff"), None)
            f.text(W / 2, L / 2, "shorting wall", "l", dx=14, dy=6)
            f.feed(0, L * 0.25)
        else:
            f.feed(0, -L * 0.2)
    if kind != "cp":
        f.dim(-W / 2, L / 2, W / 2, L / 2, f"W = {_len(W)}", 16, side=1)
    f.dim(-W / 2, -L / 2, -W / 2, L / 2, f"L = {_len(L)}", 16, side=1)
    _side_strip(f, v, W, sw, sl)


def _circ_patch(f: Fig, v):
    a = _num(v, "a", default=1.0)
    sw, sl = _patch_frame(f, 2 * a, 2 * a)
    f.circle(0, 0, a, f.gradient(-a, a, a, -a, COPPER_LIGHT, COPPER), COPPER_EDGE, 1.2)
    f.feed(a * 0.3, -a * 0.35)
    f.dim(-a, 0, a, 0, f"⌀ {_len(2 * a)}", 0)
    _side_strip(f, v, 2 * a, sw, sl)


def _tri_patch(f: Fig, v):
    s = _num(v, "a_side", default=1.0)
    ht = s * math.sqrt(3) / 2
    sw, sl = _patch_frame(f, s, ht, 1.45)
    pts = [(-s / 2, -ht / 3), (s / 2, -ht / 3), (0, 2 * ht / 3)]
    pts = [(x, y - ht / 6) for x, y in pts]
    f.poly(pts, f.gradient(-s / 2, ht / 2, s / 2, -ht / 2, COPPER_LIGHT, COPPER), COPPER_EDGE, 1.2)
    f.feed(0, ht * 0.05)
    f.dim(-s / 2, -ht / 2, s / 2, -ht / 2, f"side {_len(s)}", 14, side=-1)
    _side_strip(f, v, s, sw, sl)


def _ring_patch(f: Fig, v):
    b = _num(v, "b_out", default=1.0)
    a = _num(v, "a_in", default=b * 0.5)
    sw, sl = _patch_frame(f, 2 * b, 2 * b)
    path = QPainterPath()
    c = f.pt(0, 0)
    path.addEllipse(c, b * f.s, b * f.s)
    path.addEllipse(c, a * f.s, a * f.s)
    f.p.setBrush(f.gradient(-b, b, b, -b, COPPER_LIGHT, COPPER))
    f.p.setPen(f.pen(COPPER_EDGE, 1.2))
    f.p.drawPath(path)
    f.feed((a + b) / 2 * 0.7071, -(a + b) / 2 * 0.7071)
    f.dim(-b, b, b, b, f"outer ⌀ {_len(2 * b)}", 10, side=1)
    f.text(0, 0, f"inner ⌀ {_len(2 * a)}", "c", colour=INK)
    _side_strip(f, v, 2 * b, sw, sl)


def _pifa(f: Fig, v):
    L = _num(v, "L", default=1.0)
    h = _num(v, "h", default=L * 0.25)
    hh = max(h, L * 0.12)
    f.fit(-L * 0.35, -hh * 0.6, L * 1.35, hh * 1.9)
    f.ground(-L * 0.25, L * 1.25, 0)
    f.box(0, hh, L, f.minw(4), COPPER, COPPER_EDGE, 1.0)
    f.box(0, 0, f.minw(4), hh, METAL_DARK, INK, 0.8)
    xf = L * 0.18
    f.line(xf, 0, xf, hh, COPPER_EDGE, 2.2)
    f.feed(xf, 0)
    f.text(0, hh / 2, "short", "r", dx=-6)
    f.dim(0, hh, L, hh, f"L = {_len(L)}", 16, side=1)
    f.dim(L, 0, L, hh, f"h = {_len(h)}" if _num(v, "h") else "", 18, side=-1)


def _stacked(f: Fig, v):
    L = _num(v, "L", default=1.0)
    L2 = _num(v, "L2", default=L * 1.1)
    h = _num(v, "h", default=L * 0.06)
    h2 = _num(v, "h2", default=L * 0.3)
    hd = max(h, L * 0.05)
    span = max(L, L2) * 1.5
    f.fit(-span * 0.85, -hd * 1.4, span * 0.9, hd + h2 + L * 0.3)
    f.ground(-span / 2, span / 2, 0)
    f.box(-span / 2, 0, span, hd, SUBSTRATE, SUBSTRATE_EDGE, 1.0)
    f.box(-L / 2, hd, L, f.minw(3.5), COPPER, COPPER_EDGE, 0.8)
    f.box(-L2 / 2, hd + h2, L2, f.minw(3.5), COPPER, COPPER_EDGE, 0.8)
    for x in (-L2 * 0.4, L2 * 0.4):
        f.line(x, hd + f.minw(2), x, hd + h2, QColor("#b9c4cc"), 1.0, Qt.PenStyle.DashLine)
    xf = _num(v, "feed_offset", default=L * 0.3)
    f.line(xf, 0, xf, hd, COPPER_EDGE, 2.2)
    f.feed(xf, 0)
    f.text(-span / 2, hd + h2 / 2, "air gap", "r", dx=-6)
    f.dim(-L2 / 2, hd + h2, L2 / 2, hd + h2, f"parasitic {_len(L2)}", 14, side=1)
    f.text(0, hd + f.minw(3.5), f"driven {_len(L)}", "c", colour=INK, dy=-12)
    f.dim(span / 2, hd, span / 2, hd + h2, f"gap {_len(h2)}", 12, side=-1)


# ---------------------------------------------------------------- horns and guides
def _guide_and_flare(f: Fig, a_in, a_out, flare_len, guide_len, label_ap, label_len, shade=True, steps=None):
    x0 = -guide_len
    f.fit(x0 - guide_len * 0.1, -a_out * 0.75, flare_len + (flare_len - x0) * 0.32, a_out * 0.75)
    top = [(x0, a_in / 2), (0, a_in / 2), (flare_len, a_out / 2)]
    bot = [(flare_len, -a_out / 2), (0, -a_in / 2), (x0, -a_in / 2)]
    inside = f.gradient(0, a_out / 2, 0, -a_out / 2, METAL_LIGHT, METAL)
    f.poly(top + bot, inside, METAL_DARK, 1.2)
    f.poly([(flare_len, a_out / 2), (flare_len - f.minw(5), a_out / 2), (flare_len - f.minw(5), -a_out / 2),
            (flare_len, -a_out / 2)], METAL_DARK, INK, 0.8)
    if steps:
        for x, r in steps:
            f.line(x, r, x, a_out, INK, 1.0)
    f.line(x0, a_in / 2, x0, -a_in / 2, INK, 2.0)
    f.feed(x0 + guide_len * 0.35, 0)
    f.dim(flare_len, -a_out / 2, flare_len, a_out / 2, label_ap, 16, side=-1)
    if label_len:
        f.dim(0, -a_out / 2, flare_len, -a_out / 2, label_len, 16, side=-1)


def _pyramidal(f: Fig, v):
    a1 = _num(v, "a1", default=1.0)
    b1 = _num(v, "b1", default=a1 * 0.8)
    awg = _num(v, "a_wg", default=a1 * 0.2)
    ln = _num(v, "p_len", "rho_1", default=a1 * 1.2)
    _guide_and_flare(f, awg, a1, ln, ln * 0.3, f"a1 = {_len(a1)}", f"length {_len(ln)}")
    f.text(ln / 2, a1 / 2, f"b1 = {_len(b1)}" if _num(v, "b1") else "", "c", dy=-18)


def _sectoral(f: Fig, v, plane):
    ap = _num(v, "b1" if plane == "e" else "a1", default=1.0)
    guide = _num(v, "b_wg" if plane == "e" else "a_wg", default=ap * 0.15)
    rho = _num(v, "rho", default=ap * 1.5)
    ln = max(rho * 0.9, ap * 0.5)
    name = "b1" if plane == "e" else "a1"
    _guide_and_flare(f, guide, ap, ln, ln * 0.3, f"{name} = {_len(ap)}", "")
    f.text(ln * 0.45, 0, f"{'E' if plane == 'e' else 'H'}-plane flare", "c", colour=INK)


def _conical(f: Fig, v, kind="plain"):
    dm = _num(v, "dm", "D", default=1.0)
    L = _num(v, "L", default=dm * 1.5)            # the SLANT length, apex to rim
    psi = math.asin(min(0.95, dm / (2 * L)))       # the cone's half angle
    ax = L * math.cos(psi)                         # apex to aperture along the axis
    din = dm * 0.25
    ln = ax * (1 - din / dm)                       # the flare, throat to aperture
    _guide_and_flare(f, din, dm, ln, ln * 0.4, f"⌀ {_len(dm)}", "")
    apex = ln - ax
    for sgn in (1, -1):
        f.line(0, sgn * din / 2, apex, 0, DIM, 0.9, Qt.PenStyle.DashLine)
    f.dim(apex, 0, ln, dm / 2, f"slant {_len(L)}", 14, side=1)
    if kind == "corrugated":
        n = 16
        for i in range(1, n):
            x = ln * i / n
            r = din / 2 + (dm - din) / 2 * i / n
            d = max(_num(v, "slot_depth", default=dm * 0.03), dm * 0.025)
            f.line(x, r, x, r - d, INK, 1.0)
            f.line(x, -r, x, -r + d, INK, 1.0)
        f.text(ln / 2, 0, "corrugated wall", "c", colour=INK)
    if kind == "dual":
        ds = _num(v, "d_step", default=din * 1.3)
        f.box(-ln * 0.12, -ds / 2, ln * 0.12, ds, METAL_LIGHT, METAL_DARK, 1.0)
        f.text(-ln * 0.06, ds / 2, "step", "c", dy=-12)


def _diagonal(f: Fig, v):
    a = _num(v, "a_ap", default=1.0)
    f.fit(-a * 1.3, -a * 0.95, a * 1.6, a * 0.95)
    # side view on the left, the diamond aperture seen from the front on the right
    ln = a * 0.9
    f.poly([(-a * 1.1, a * 0.08), (-a * 1.1 + ln, a / 2), (-a * 1.1 + ln, -a / 2), (-a * 1.1, -a * 0.08)],
           f.gradient(0, a / 2, 0, -a / 2, METAL_LIGHT, METAL), METAL_DARK, 1.2)
    c = a * 0.95
    f.poly([(c, a * 0.62), (c + a * 0.62, 0), (c, -a * 0.62), (c - a * 0.62, 0)], METAL_LIGHT, METAL_DARK, 2.0)
    f.line(c - a * 0.3, -a * 0.3, c + a * 0.3, a * 0.3, FEED, 1.6)
    f.text(c, -a * 0.62, "aperture, diagonal E", "c", dy=16)
    f.dim(-a * 1.1 + ln, -a / 2, -a * 1.1 + ln, a / 2, f"side {_len(a)}", 14, side=-1)


def _open_guide(f: Fig, v):
    a = _num(v, "a_wg", default=1.0)
    b = _num(v, "b_wg", default=a / 2)
    f.fit(-a * 1.25, -a * 0.85, a * 1.35, a * 0.8)
    d = a * 0.45
    f.poly([(-a / 2, b / 2), (-a / 2 + d, b / 2 + d * 0.5), (a / 2 + d, b / 2 + d * 0.5), (a / 2, b / 2)], METAL_LIGHT, METAL_DARK, 1.0)
    f.poly([(a / 2, b / 2), (a / 2 + d, b / 2 + d * 0.5), (a / 2 + d, -b / 2 + d * 0.5), (a / 2, -b / 2)], METAL, METAL_DARK, 1.0)
    f.box(-a / 2 - a * 0.25, -b / 2 - a * 0.25, a * 1.5, b + a * 0.5, None, METAL_DARK, 1.0)
    f.box(-a / 2, -b / 2, a, b, QColor("#3b4650"), METAL_DARK, 2.0)
    f.line(0, -b / 2 * 0.8, 0, b / 2 * 0.8, FEED, 1.6)
    f.text(0, 0, "E", "l", colour=QColor("#ffffff"), dx=6)
    f.dim(-a / 2, -b / 2, a / 2, -b / 2, f"a = {_len(a)}", 12 + a * 0.25 * f.s / f.u, side=-1)
    f.dim(-a / 2, -b / 2, -a / 2, b / 2, f"b = {_len(b)}", 12 + a * 0.25 * f.s / f.u, side=1)


# ---------------------------------------------------------------- travelling wave
def _helix(f: Fig, v, kind="axial"):
    D = _num(v, "D_helix", default=1.0)
    S = _num(v, "S", default=D * 0.7)
    N = _num(v, "N", default=6 if kind == "axial" else 4)
    turns = max(1, int(round(N)))
    Lax = turns * S                                # the winding, first turn to last
    n = min(turns, 60)                             # beyond 60 the pitch is drawn wider, the length true
    S = Lax / n
    Dg = _num(v, "D_gnd", default=D * 2.2) if kind == "axial" else D * 1.6
    f.fit(-Dg * 0.75, -Lax * 0.08, Dg * 0.75 + D * 0.3, Lax * 1.08)
    f.ground(-Dg / 2, Dg / 2, 0)
    pts_back, pts_front = [], []
    for i in range(n * 60 + 1):
        t = i / 60.0
        x = D / 2 * math.cos(2 * math.pi * t)
        y = S * 0.4 + t * S
        (pts_front if math.sin(2 * math.pi * t) >= 0 else pts_back).append((i, x, y))

    def runs(pts):
        out, cur, last = [], [], None
        for i, x, y in pts:
            if last is not None and i != last + 1:
                out.append(cur)
                cur = []
            cur.append((x, y))
            last = i
        if cur:
            out.append(cur)
        return out
    for run in runs(pts_back):
        f.curve(run, QColor("#d9b18a"), 2.2)
    for run in runs(pts_front):
        f.curve(run, COPPER_EDGE, 3.4)
        f.curve(run, COPPER, 2.3)
    f.feed(D / 2, S * 0.4)
    f.dim(-D / 2, Lax + S * 0.6, D / 2, Lax + S * 0.6, f"⌀ {_len(D)}", 6, side=1)
    # the winding's true length (nominal when the turn count or pitch is)
    length = _UNSET if isinstance(N, _Nominal) else _len(Lax)
    f.dim(Dg / 2, S * 0.4, Dg / 2, S * 0.4 + Lax, f"{length}, {_count(N)} turns", 10, side=-1)


def _yagi(f: Fig, v):
    boom = _num(v, "boom_length", default=1.0)
    lr = _num(v, "reflector_length", default=boom * 0.4)
    ld = _num(v, "driven_length", default=lr * 0.95)
    lam = lr / 0.5
    n_dir = max(1, min(12, int(round((boom - 0.2 * lam) / (0.3 * lam)))))
    f.fit(-boom * 0.12, -lr * 0.7, boom * 1.12, lr * 0.75)
    f.box(0, -f.minw(2.5), boom, f.minw(5), METAL, METAL_DARK, 0.8)
    xs = [0, 0.2 * lam] + [0.2 * lam + (boom - 0.2 * lam) * (i + 1) / n_dir for i in range(n_dir)]
    lengths = [lr, ld] + [ld * (0.9 - 0.01 * i) for i in range(n_dir)]
    for i, (x, l) in enumerate(zip(xs, lengths)):
        if i == 1:
            f.wire(x, f.minw(4), x, l / 2, 2.6)
            f.wire(x, -f.minw(4), x, -l / 2, 2.6)
            f.feed(x, 0)
        else:
            f.wire(x, -l / 2, x, l / 2, 2.4)
    f.text(0, -lr / 2, "reflector", "c", dy=12)
    f.text(xs[1], ld / 2, "driven", "c", dy=-12)
    f.text(xs[-1], lengths[-1] / 2, "directors →", "r", dy=-12)
    f.dim(0, -lr / 2, boom, -lr / 2, f"boom {_len(boom)}", 30, side=-1)


def _lpda(f: Fig, v):
    Lmax = _num(v, "L_max", default=1.0)
    count = _num(v, "N_elements", default=8)
    N = max(2, min(int(round(count)), 80))
    tau = min(_num(v, "tau", default=0.88), 0.99)
    lens = [Lmax * tau ** i for i in range(N)]
    xs, x = [], 0.0
    sigma = _num(v, "sigma", default=0.16)
    for i, l in enumerate(lens):
        xs.append(x)
        x += 2 * sigma * l
    span = xs[-1]
    f.fit(-span * 0.1 - Lmax * 0.9, -Lmax * 0.62, span * 1.15 + Lmax * 0.1, Lmax * 0.62)
    f.line(xs[0], f.minw(3), xs[-1], f.minw(3), METAL_DARK, 2.2)
    f.line(xs[0], -f.minw(3), xs[-1], -f.minw(3), METAL_DARK, 2.2)
    for i, (x, l) in enumerate(zip(xs, lens)):
        sgn = 1 if i % 2 == 0 else -1
        f.wire(x, sgn * f.minw(3), x, sgn * l / 2, 2.2)
        f.wire(x, -sgn * f.minw(3), x, -sgn * l / 2, 2.2)
    f.feed(xs[-1], 0)
    f.dim(xs[0], -Lmax / 2, xs[0], Lmax / 2, f"longest {_len(Lmax)}", 26, side=1)
    f.text(span / 2, -Lmax / 2, f"{_count(count)} elements", "c", dy=14)


def _rhombic(f: Fig, v, kind="rhombic"):
    L = _num(v, "L", default=1.0)
    a = math.radians(_num(v, "half_angle_deg", default=25))
    if kind == "rhombic":
        dx, dy = L * math.cos(a), L * math.sin(a)
        f.fit(-dx * 0.08, -dy * 1.4, 2 * dx * 1.06, dy * 1.4)
        f.curve([(0, 0), (dx, dy), (2 * dx, 0), (dx, -dy), (0, 0)], COPPER_EDGE, 3.0)
        f.curve([(0, 0), (dx, dy), (2 * dx, 0), (dx, -dy), (0, 0)], COPPER, 2.0)
        f.box(2 * dx - f.minw(4), -f.minw(9), f.minw(8), f.minw(18), QColor("#ffffff"), INK, 1.2)
        f.text(2 * dx, 0, "R load", "l", dx=12)
        f.feed(0, 0)
        f.dim(0, 0, dx, dy, f"leg {_len(L)}", 14, side=1)
    else:
        dx, dy = L * math.cos(a), L * math.sin(a)
        f.fit(-dx * 0.08, -dy * 1.35, dx * 1.12, dy * 1.35)
        f.wire(0, 0, dx, dy, 2.2)
        f.wire(0, 0, dx, -dy, 2.2)
        f.feed(0, 0)
        f.dim(0, 0, dx, dy, f"leg {_len(L)}", 14, side=1)


def _long_wire(f: Fig, v):
    L = _num(v, "L", default=1.0)
    hgt = L * 0.14
    f.fit(-L * 0.08, -hgt * 0.8, L * 1.2, hgt * 2.4)
    f.ground(-L * 0.04, L * 1.04, 0)
    f.wire(0, hgt, L, hgt, 2.0)
    f.line(0, 0, 0, hgt, COPPER_EDGE, 1.6)
    f.line(L, 0, L, hgt, INK, 1.2)
    f.box(L - f.minw(4), hgt * 0.3, f.minw(8), hgt * 0.4, QColor("#ffffff"), INK, 1.2)
    f.feed(0, 0)
    f.text(L, hgt * 0.5, "R load", "l", dx=10)
    f.dim(0, hgt, L, hgt, f"L = {_len(L)}", 14, side=1)


def _leaky(f: Fig, v):
    L = _num(v, "L", default=1.0)
    w = L * 0.12
    f.fit(-L * 0.08, -w * 2.6, L * 1.08, w * 2.2)
    f.box(0, -w / 2, L, w, f.gradient(0, w / 2, 0, -w / 2, METAL_LIGHT, METAL), METAL_DARK, 1.2)
    f.box(0, -w * 0.08, L, w * 0.16, QColor("#2f3a44"), None)
    for i in range(8):
        x = L * (i + 0.5) / 8
        f.line(x, w / 2, x + L * 0.05, w / 2 + w * 1.1, RAY, 1.2)
    f.feed(0, 0)
    f.dim(0, -w / 2, L, -w / 2, f"L = {_len(L)}", 14, side=-1)
    f.text(L / 2, -w / 2, "slot radiating along its length", "c", dy=34)


# ---------------------------------------------------------------- UWB
def _bowtie(f: Fig, v):
    Ls = _num(v, "Lside", default=1.0)
    We = _num(v, "Wend", default=Ls)
    f.fit(-Ls * 1.25, -We * 0.72, Ls * 1.6, We * 0.72)
    g = Ls * 0.04
    f.poly([(-g, 0), (-Ls, We / 2), (-Ls, -We / 2)], f.gradient(-Ls, We / 2, 0, -We / 2, COPPER_LIGHT, COPPER), COPPER_EDGE, 1.2)
    f.poly([(g, 0), (Ls, We / 2), (Ls, -We / 2)], f.gradient(0, We / 2, Ls, -We / 2, COPPER_LIGHT, COPPER), COPPER_EDGE, 1.2)
    f.feed(0, 0)
    f.dim(g, -We / 2, Ls, -We / 2, f"arm {_len(Ls)}", 14, side=-1)
    f.dim(Ls, -We / 2, Ls, We / 2, f"end {_len(We)}", 12, side=-1)


def _spiral(f: Fig, v, kind="arch"):
    ro = _num(v, "r_out", default=1.0)
    ri = min(_num(v, "r_in", "r0", default=ro * 0.08), ro * 0.95)
    f.fit(-ro * 1.2, -ro * 1.2, ro * 1.3, ro * 1.2)
    if kind == "arch":
        turns = _num(v, "n_turns", default=5.0)
    else:                                          # r = r0 exp(a phi): ln(ro/r0) / (2 pi a) turns
        turns = _num(v, "turns_required", default=math.log(ro / ri) / (2 * math.pi * _num(v, "a_growth", default=0.22)))
    for arm in (0, math.pi):
        pts = []
        steps = max(600, int(160 * turns))
        for i in range(steps):
            t = i / (steps - 1)
            if kind == "arch":
                r = ri + (ro - ri) * t
            else:
                r = ri * (ro / ri) ** t
            ang = arm + 2 * math.pi * turns * t
            pts.append((r * math.cos(ang), r * math.sin(ang)))
        f.curve(pts, COPPER_EDGE, 4.4)
        f.curve(pts, COPPER, 3.0)
    f.feed(0, 0)
    f.dim(-ro, -ro, ro, -ro, f"⌀ {_len(2 * ro)}" if _num(v, "r_out") else "", 10, side=-1)


def _discone(f: Fig, v):
    Dd = _num(v, "disc_diameter", default=0.7)
    Db = _num(v, "cone_base_diameter", default=1.0)
    H = _num(v, "cone_height", default=0.85)
    top = _num(v, "cone_top_diameter", default=Db * 0.03)
    f.fit(-Db * 0.72, -H * 1.12, Db * 0.85, H * 0.3)
    f.poly([(-top / 2, -f.minw(4)), (top / 2, -f.minw(4)), (Db / 2, -H), (-Db / 2, -H)],
           f.gradient(-Db / 2, 0, Db / 2, 0, COPPER_LIGHT, COPPER), COPPER_EDGE, 1.2)
    f.box(-Dd / 2, 0, Dd, f.minw(4), COPPER, COPPER_EDGE, 1.0)
    f.feed(0, -f.minw(2))
    f.dim(-Dd / 2, f.minw(4), Dd / 2, f.minw(4), f"disc ⌀ {_len(Dd)}", 12, side=1)
    f.dim(-Db / 2, -H, Db / 2, -H, f"cone ⌀ {_len(Db)}", 14, side=-1)
    f.dim(Db / 2, -H, Db / 2, 0, f"{_len(H)}", 12, side=-1)


def _cone_monopole(f: Fig, v):
    Db = _num(v, "base_diameter", default=1.0)
    H = _num(v, "height", default=0.5)
    f.fit(-Db * 0.7, -H * 0.25, Db * 0.8, H * 1.2)
    f.ground(-Db * 0.66, Db * 0.66, 0)
    f.poly([(0, f.minw(3)), (Db / 2, H), (-Db / 2, H)], f.gradient(-Db / 2, 0, Db / 2, 0, COPPER_LIGHT, COPPER), COPPER_EDGE, 1.2)
    f.ellipse(0, H, Db / 2, Db * 0.08, COPPER_LIGHT, COPPER_EDGE, 1.0)
    f.feed(0, f.minw(2))
    f.dim(-Db / 2, H + Db * 0.08, Db / 2, H + Db * 0.08, f"⌀ {_len(Db)}", 10, side=1)
    f.dim(Db * 0.62, 0, Db * 0.62, H, f"h = {_len(H)}", 0, side=-1)


def _planar_monopole(f: Fig, v, kind="disc"):
    if kind == "disc":
        r = _num(v, "r_disc", default=1.0)
        w, hgt = 2 * r, 2 * r
    else:
        hgt = _num(v, "Lp", default=1.0)
        w = _num(v, "Wp", default=hgt)
    gap = hgt * 0.05
    f.fit(-w * 1.05, -hgt * 0.3, w * 1.15, hgt * 1.2)
    f.ground(-w, w, 0)
    if kind == "disc":
        f.circle(0, gap + r, r, f.gradient(-r, 2 * r, r, 0, COPPER_LIGHT, COPPER), COPPER_EDGE, 1.2)
        f.dim(-r, gap + r, r, gap + r, f"⌀ {_len(2 * r)}", 0)
    else:
        f.box(-w / 2, gap, w, hgt, f.gradient(-w / 2, hgt, w / 2, 0, COPPER_LIGHT, COPPER), COPPER_EDGE, 1.2)
        f.dim(-w / 2, gap + hgt, w / 2, gap + hgt, f"W = {_len(w)}", 12, side=1)
        f.dim(w / 2, gap, w / 2, gap + hgt, f"L = {_len(hgt)}", 12, side=-1)
    f.line(0, 0, 0, gap, COPPER_EDGE, 2.0)
    f.feed(0, 0)


def _vivaldi(f: Fig, v):
    Wt = _num(v, "W_ap", default=1.0)
    L = _num(v, "Lax", default=Wt * 1.4)
    W = max(Wt, L / 4)                             # the board, widened on paper for a slender design
    f.fit(-L * 0.08, -W * 0.85, L * 1.45, W * 0.75)
    f.box(0, -W * 0.62, L, W * 1.24, SUBSTRATE, SUBSTRATE_EDGE, 1.0)
    s0 = Wt * 0.01                                 # the slot's half width at the feed
    k = math.log((Wt / 2) / s0) / (L * 0.8)        # opening to the aperture W_ap, not the board
    for sgn in (1, -1):
        edge = [(L * 0.2 + x, sgn * s0 * math.exp(k * x)) for x in [L * 0.8 * i / 60 for i in range(61)]]
        pts = [(0, sgn * s0), (L * 0.2, sgn * s0)] + edge + [(L, sgn * W * 0.6), (0, sgn * W * 0.6)]
        f.poly(pts, f.gradient(0, sgn * W * 0.6, L, 0, COPPER_LIGHT, COPPER), COPPER_EDGE, 1.0)
    f.feed(L * 0.2, 0)
    f.dim(L, -Wt / 2, L, Wt / 2, f"aperture {_len(Wt)}", 12, side=-1)
    f.dim(0, -W * 0.62, L, -W * 0.62, f"length {_len(_num(v, 'Lax'))}" if _num(v, "Lax") else "", 12, side=-1)


# ---------------------------------------------------------------- reflectors
def _dish(f: Fig, v, kind="prime"):
    D = _num(v, "W", default=1.0) if kind == "cyl" else _num(v, "D", default=1.0)
    F = _num(v, "F", "focal_length_m")
    if F is None:                                  # a dual reflector gives f/D, not F: F = (f/D) D
        F = D * _num(v, "f_over_D", default=0.4)
    depth = D * D / (16 * F)
    if kind == "offset":
        h0 = _num(v, "h0", default=D * 0.6)
        ys = [h0 - D / 2 + D * i / 60 for i in range(61)]
        f.fit(-F * 0.25, -D * 0.2, F * 1.25, h0 + D * 0.7)
        pts = [(y * y / (4 * F), y) for y in ys]
        f.curve(pts, METAL_DARK, 6.0)
        f.curve(pts, METAL, 4.0)
        f.feed(F, 0)
        f.poly([(F, 0), (F + D * 0.12, D * 0.05), (F + D * 0.12, -D * 0.05)], METAL_LIGHT, METAL_DARK, 1.0)
        for y in (ys[0], ys[-1]):
            f.line(F, 0, y * y / (4 * F), y, RAY, 1.0, Qt.PenStyle.DashLine)
        xl = min(p[0] for p in pts) - D * 0.08      # the projected aperture: vertical, rim to rim
        f.dim(xl, ys[0], xl, ys[-1], f"D = {_len(D)}", 14, side=1)
        return
    f.fit(-depth - D * 0.35, -D * 0.62, F + D * 0.3, D * 0.62)
    ys = [-D / 2 + D * i / 80 for i in range(81)]
    pts = [(y * y / (4 * F) - depth, y) for y in ys]
    f.curve(pts, METAL_DARK, 6.5)
    f.curve(pts, METAL, 4.5)
    if kind in ("prime", "cyl"):
        fx = F - depth
        f.poly([(fx, 0), (fx + D * 0.1, D * 0.05), (fx + D * 0.1, -D * 0.05)], METAL_LIGHT, METAL_DARK, 1.0)
        f.feed(fx, 0)
        for y in (-D * 0.5, D * 0.5):
            f.line(-depth + y * y / (4 * F), y, fx + D * 0.1, 0, METAL_DARK, 0.8)
        f.dim(-depth, -D / 2, F - depth, -D / 2, f"F = {_len(F)}", 16, side=-1)
        if kind == "cyl":
            f.text(F * 0.4, D / 2, "parabolic cylinder, line feed", "c", dy=-14)
    else:
        fx = F - depth
        sr = _num(v, "Ds", default=D * 0.2) / 2
        if kind == "cass":
            sx = fx * 0.78
            sub = [(sx + (y / sr) ** 2 * sr * 0.3, y) for y in [-sr + 2 * sr * i / 30 for i in range(31)]]
            f.curve(sub, METAL_DARK, 4.4)
            f.curve(sub, METAL, 3.0)
        else:
            sx = fx * 1.2
            sub = [(sx - (y / sr) ** 2 * sr * 0.35, y) for y in [-sr + 2 * sr * i / 30 for i in range(31)]]
            f.curve(sub, METAL_DARK, 4.4)
            f.curve(sub, METAL, 3.0)
        f.poly([(-depth + D * 0.02, 0), (-depth + D * 0.14, D * 0.045), (-depth + D * 0.14, -D * 0.045)], METAL_LIGHT, METAL_DARK, 1.0)
        f.feed(-depth + D * 0.02, 0)
        f.line(-depth + D * 0.14, D * 0.04, sx, sr * 0.8, RAY, 1.0, Qt.PenStyle.DashLine)
        f.line(-depth + D * 0.14, -D * 0.04, sx, -sr * 0.8, RAY, 1.0, Qt.PenStyle.DashLine)
        f.text(sx, sr, "subreflector", "c", dy=-14)
        f.dim(sx + D * 0.06, -sr, sx + D * 0.06, sr, f"⌀ {_len(2 * sr)}", 8, side=-1)
        f.dim(-depth, -D / 2, fx, -D / 2, f"F = {_len(F)}", 16, side=-1)
    f.dim(-depth - D * 0.04, -D / 2, -depth - D * 0.04, D / 2, f"{'W' if kind == 'cyl' else 'D'} = {_len(D)}", 14, side=1)


def _corner(f: Fig, v, angle):
    S = _num(v, "S", default=0.5)
    Ls = _num(v, "Lside", default=1.0)
    half = math.radians(angle / 2)
    f.fit(-Ls * 0.1, -Ls * math.sin(half) * 1.2, Ls * math.cos(half) * 1.15, Ls * math.sin(half) * 1.2)
    for sgn in (1, -1):
        x, y = Ls * math.cos(half), sgn * Ls * math.sin(half)
        f.line(0, 0, x, y, METAL_DARK, 7.0)
        f.line(0, 0, x, y, METAL, 5.0)
    f.circle(S, 0, f.minw(5), COPPER, COPPER_EDGE, 1.0)
    f.feed(S, 0, 2.5)
    f.text(S, 0, "dipole (end-on)", "l", dx=12)
    f.dim(0, 0, S, 0, f"S = {_len(_num(v, 'S'))}" if _num(v, "S") else "", 0)
    f.dim(0, 0, Ls * math.cos(half), -Ls * math.sin(half), f"side {_len(_num(v, 'Lside'))}" if _num(v, "Lside") else "", 14, side=-1)
    f.text(0, 0, f"{angle}°", "r", dx=-8)


# ---------------------------------------------------------------- slots
def _sheet(f: Fig, x0, y0, w, h):
    f.box(x0, y0, w, h, f.gradient(x0, y0 + h, x0 + w, y0, METAL_LIGHT, METAL), METAL_DARK, 1.2)


def _slot(f: Fig, v, kind="plain"):
    L = _num(v, "L", "slot_length", default=1.0)
    w = max(_num(v, "w", "slot_width", default=L * 0.05), L * 0.05)
    f.fit(-L * 0.78, -L * 0.45, L * 0.78, L * 0.45)
    _sheet(f, -L * 0.7, -L * 0.36, L * 1.4, L * 0.72)
    if kind == "folded":
        f.box(-L / 2, w * 0.9, L, w, QColor("#26313a"), None)
        f.box(-L / 2, -w * 1.9, L, w, QColor("#26313a"), None)
        f.box(-L / 2, -w * 1.9, w, w * 3.8, QColor("#26313a"), None)
        f.box(L / 2 - w, -w * 1.9, w, w * 3.8, QColor("#26313a"), None)
        f.feed(0, -w * 1.4)
    else:
        f.box(-L / 2, -w / 2, L, w, QColor("#26313a"), None)
        f.feed(0, 0)
    f.dim(-L / 2, w * 2, L / 2, w * 2, f"L = {_len(L)}", 14, side=1)
    f.text(0, -L * 0.36, "metal sheet", "c", dy=14)


def _cavity_slot(f: Fig, v):
    L = _num(v, "slot_length", default=1.0)
    d = _num(v, "cavity_depth", default=L * 0.5)
    f.fit(-L * 0.9, -d * 1.3, L * 1.2, L * 0.35)
    f.box(-L * 0.8, 0, L * 1.6, f.minw(5), METAL, METAL_DARK, 1.0)
    f.box(-L / 2, 0, L, f.minw(5), QColor("#26313a"), None)
    f.poly([(-L * 0.6, 0), (-L * 0.6, -d), (L * 0.6, -d), (L * 0.6, 0)], None, METAL_DARK, 3.0, closed=False)
    f.box(-L * 0.6, -d, L * 1.2, d, QColor(215, 222, 228, 120), None)
    f.feed(0, f.minw(2))
    f.dim(-L / 2, f.minw(5), L / 2, f.minw(5), f"slot {_len(L)}", 12, side=1)
    f.dim(L * 0.6, -d, L * 0.6, 0, f"depth {_len(d)}", 12, side=-1)
    f.text(0, -d * 0.65, "cavity", "c", colour=INK)


def _guide_slots(f: Fig, v, kind):
    a = _num(v, "a_wg", default=1.0)
    sl = _num(v, "slot_length", default=a * 0.6)
    sp = _num(v, "spacing", default=sl * 1.4)
    off = _num(v, "offset", "x1", default=a * 0.12)
    count = 1 if kind == "single" else _num(v, "N", default=6)
    n = max(1, min(int(round(count)), 80))
    length = sp * (n + 0.6) if n > 1 else sl * 2.4
    f.fit(-length * 0.08 - max(a * 1.1, length * 0.16), -a * 1.0, length * 1.08, a * 0.95)
    _sheet(f, 0, -a / 2, length, a)
    f.line(0, 0, length, 0, QColor(90, 106, 119, 110), 0.8, Qt.PenStyle.DashLine)
    for i in range(n):
        x = (length - (n - 1) * sp) / 2 + i * sp if n > 1 else length / 2
        o = off * (1 if (kind != "resonant" or i % 2 == 0) else -1)
        f.box(x - sl / 2, o - a * 0.035, sl, a * 0.07, QColor("#26313a"), None)
    f.feed(0, 0)
    f.dim(0, -a / 2, 0, a / 2, f"a = {_len(a)}", 12, side=1)
    if n > 1:
        x0 = (length - (n - 1) * sp) / 2
        f.dim(x0, -a / 2, x0 + sp, -a / 2, f"pitch {_len(sp)}", 14, side=-1)
        f.text(length / 2, -a / 2, f"{_count(count)} slots", "c", dy=34)
    else:
        f.dim(length / 2 - sl / 2, -a / 2, length / 2 + sl / 2, -a / 2, f"slot {_len(sl)}", 14, side=-1)
    f.text(length / 2, a / 2, "broad wall of the guide", "c", dy=-12)


# ---------------------------------------------------------------- lenses and dielectrics
def _rays(f: Fig, x0, x1, ys):
    for y in ys:
        f.line(x0, y, x1, y, RAY, 1.2)
        f._arrow(f.pt(x1, y), f.pt(x0, y))


def _luneburg(f: Fig, v):
    R = _num(v, "R", default=1.0)
    f.fit(-R * 1.6, -R * 1.2, R * 2.3, R * 1.2)
    for i in range(6, 0, -1):
        r = R * i / 6
        c = QColor(214 - 18 * (6 - i), 228 - 12 * (6 - i), 238 - 8 * (6 - i))
        f.circle(0, 0, r, c, QColor(120, 145, 165), 0.6)
    f.poly([(-R, 0), (-R - R * 0.35, R * 0.18), (-R - R * 0.35, -R * 0.18)], METAL_LIGHT, METAL_DARK, 1.0)
    f.feed(-R, 0)
    _rays(f, R * 1.1, R * 2.1, [-R * 0.6, -R * 0.2, R * 0.2, R * 0.6])
    f.dim(-R, -R, R, -R, f"⌀ {_len(2 * R)}", 12, side=-1)
    f.text(0, R, "graded index", "c", dy=-14)


def _dielectric_lens(f: Fig, v, kind="hyper"):
    D = _num(v, "D_ap", default=1.0)
    F = _num(v, "F", default=D)
    f.fit(-F * 1.15, -D * 0.72, D * 1.5, D * 0.72)
    ys = [-D / 2 + D * i / 40 for i in range(41)]
    if kind == "hyper":
        t = D * 0.35
        prof = [(-t * (1 - (y / (D / 2)) ** 2), y) for y in ys]
        f.poly(prof + [(D * 0.06, D / 2), (D * 0.06, -D / 2)], QColor(206, 226, 236), CERAMIC_EDGE, 1.2)
    elif kind == "plates":
        for i in range(9):
            y = -D / 2 + D * i / 8
            w = D * 0.1 + D * 0.3 * (abs(y) / (D / 2)) ** 2
            f.box(0, y - f.minw(1.5), w, f.minw(3), METAL, METAL_DARK, 0.8)
        f.text(D * 0.2, D / 2, "metal plates", "c", dy=-14)
    f.poly([(-F, 0), (-F - D * 0.14, D * 0.07), (-F - D * 0.14, -D * 0.07)], METAL_LIGHT, METAL_DARK, 1.0)
    f.feed(-F, 0)
    for y in (-D * 0.45, 0, D * 0.45):
        f.line(-F, 0, -D * 0.35 * (1 - (y / (D / 2)) ** 2) if kind == "hyper" else 0, y, RAY, 1.0, Qt.PenStyle.DashLine)
    _rays(f, D * 0.5, D * 1.35, [-D * 0.3, 0, D * 0.3])
    f.dim(D * 0.06, -D / 2, D * 0.06, D / 2, f"⌀ {_len(D)}", 14, side=-1)
    f.dim(-F, -D / 2, 0, -D / 2, f"F = {_len(_num(v, 'F'))}" if _num(v, "F") else "", 10, side=-1)


def _zone_plate(f: Fig, v):
    """Zone m ends at r_m = sqrt(m lambda F + (m lambda / 2)^2), m = 1..M; an L-level
    phase plate steps every 2/L of a zone, an amplitude plate (L = 1) blocks the
    even zones, the centre open."""
    ro = _num(v, "r_outer", default=1.0)
    M = max(1, min(int(round(_num(v, "M", default=4))), 40))
    levels = max(1, int(round(_num(v, "phase_levels", default=1))))
    lam, F = _num(v, "lambda0"), _num(v, "F")
    if isinstance(lam, float) and isinstance(F, float) and not isinstance(F, _Nominal):
        r = lambda m: math.sqrt(m * lam * F + (m * lam / 2) ** 2)
    else:
        r = lambda m: ro * math.sqrt(m / M)
    f.fit(-ro * 1.2, -ro * 1.2, ro * 1.25, ro * 1.2)
    sub = 1 if levels == 1 else max(1, levels // 2)
    edges = [r(j / sub) for j in range(1, M * sub + 1)]
    shades = [QColor("#f7fafb"), QColor("#c9dbe4"), QColor("#a9c4d2"), QColor("#8aaec0")]
    for k in range(len(edges) - 1, -1, -1):
        if levels == 1:
            fill = METAL if k % 2 == 1 else QColor("#f7fafb")
        else:
            fill = shades[k % min(levels, 4)]
        f.circle(0, 0, edges[k], fill, METAL_DARK, 0.6)
    f.feed(0, 0, 3.0)
    f.dim(-ro, -ro, ro, -ro, f"⌀ {_len(2 * ro)}", 12, side=-1)
    f.text(0, ro, "alternate zones blocked" if levels == 1 else f"{levels}-level phase steps", "c", dy=-14)


def _dra(f: Fig, v, kind="cyl"):
    if kind == "cyl":
        a = _num(v, "a", default=1.0)
        h = _num(v, "h", default=a)
        w = 2 * a
    elif kind == "hemi":
        a = _num(v, "a", default=1.0)
        w, h = 2 * a, a
    else:
        w = _num(v, "w", default=1.0)
        h = _num(v, "d", default=w * 0.5)
    f.fit(-w * 1.05, -h * 0.5, w * 1.25, h * 1.45)
    f.ground(-w, w, 0)
    if kind == "cyl":
        e = w * 0.12
        f.box(-a, 0, w, h, f.gradient(-a, 0, a, 0, QColor("#eef4f8"), CERAMIC), CERAMIC_EDGE, 1.2)
        f.ellipse(0, h, a, e, QColor("#f2f7fa"), CERAMIC_EDGE, 1.2)
        f.dim(-a, h + e, a, h + e, f"⌀ {_len(2 * a)}", 10, side=1)
        f.dim(a, 0, a, h, f"h = {_len(h)}", 14, side=-1)
    elif kind == "hemi":
        path = QPainterPath()
        c = f.pt(0, 0)
        r = a * f.s
        path.moveTo(QPointF(c.x() - r, c.y()))
        path.arcTo(QRectF(c.x() - r, c.y() - r, 2 * r, 2 * r), 180, -180)
        path.closeSubpath()
        g = QRadialGradient(f.pt(-a * 0.3, a * 0.6), r * 1.2)
        g.setColorAt(0, QColor("#f4f8fb"))
        g.setColorAt(1, CERAMIC)
        f.p.setBrush(QBrush(g))
        f.p.setPen(f.pen(CERAMIC_EDGE, 1.2))
        f.p.drawPath(path)
        f.dim(-a, 0, a, 0, f"⌀ {_len(2 * a)}", 14, side=-1)
    else:
        dd = w * 0.35
        f.poly([(-w / 2, h), (-w / 2 + dd, h + dd * 0.45), (w / 2 + dd, h + dd * 0.45), (w / 2, h)], QColor("#f2f7fa"), CERAMIC_EDGE, 1.0)
        f.poly([(w / 2, 0), (w / 2 + dd, dd * 0.45), (w / 2 + dd, h + dd * 0.45), (w / 2, h)], CERAMIC, CERAMIC_EDGE, 1.0)
        f.box(-w / 2, 0, w, h, f.gradient(-w / 2, 0, w / 2, 0, QColor("#eef4f8"), CERAMIC), CERAMIC_EDGE, 1.2)
        f.dim(-w / 2, 0, w / 2, 0, f"w = {_len(w)}", 20, side=-1)
        f.dim(-w / 2, 0, -w / 2, h, f"d = {_len(h)}", 14, side=1)
    xf = w * 0.28 if kind != "hemi" else w * 0.22
    f.line(xf, 0, xf, h * 0.55, COPPER_EDGE, 2.0)
    f.feed(xf, 0)
    er = v.get("eps_r")
    f.text(-w * 0.12, h * 0.4 if kind != "hemi" else h * 0.3, f"εr {er:g}" if isinstance(er, (int, float)) else "", "c", colour=INK)


# ---------------------------------------------------------------- dispatch
DRAWINGS: dict[str, Callable[[Fig, Mapping], None]] = {
    "half_wave_dipole": _dipole, "resonant_dipole": _dipole, "short_dipole": _dipole,
    "dipole_arbitrary_length": _dipole, "folded_dipole": _folded_dipole, "dipole_over_ground": _dipole_over_ground,
    "quarter_wave_monopole": _monopole, "inductively_loaded_monopole": lambda f, v: _monopole(f, v, "loaded"),
    "top_loaded_monopole": lambda f, v: _monopole(f, v, "hat"), "biconical": _biconical, "turnstile_dipole": _turnstile,
    "small_circular_loop": _circle_loop, "one_wavelength_circular_loop": _circle_loop,
    "halo_loop": lambda f, v: _circle_loop(f, v, "halo"), "small_square_loop": _square_loop,
    "quad_loop_square": _square_loop, "alford_loop": lambda f, v: _square_loop(f, v, "alford"),
    "multiturn_small_loop": _multiturn, "ferrite_rod_loop": _ferrite,
    "rectangular_patch": _rect_patch, "rectangular_patch_inset": lambda f, v: _rect_patch(f, v, "inset"),
    "truncated_corner_cp_patch": lambda f, v: _rect_patch(f, v, "cp"),
    "quarter_wave_shorted_patch": lambda f, v: _rect_patch(f, v, "shorted"),
    "circular_patch": _circ_patch, "triangular_patch": _tri_patch, "annular_ring_patch": _ring_patch,
    "pifa": _pifa, "stacked_patch": _stacked,
    "pyramidal_horn": _pyramidal, "e_plane_sectoral_horn": lambda f, v: _sectoral(f, v, "e"),
    "h_plane_sectoral_horn": lambda f, v: _sectoral(f, v, "h"), "conical_horn": _conical,
    "corrugated_conical_horn": lambda f, v: _conical(f, v, "corrugated"),
    "conical_horn_dual_mode": lambda f, v: _conical(f, v, "dual"), "diagonal_horn": _diagonal,
    "open_ended_waveguide": _open_guide,
    "axial_mode_helix": _helix, "normal_mode_helix": lambda f, v: _helix(f, v, "normal"),
    "yagi_uda": _yagi, "lpda": _lpda, "rhombic": _rhombic,
    "v_antenna_travelling": lambda f, v: _rhombic(f, v, "v"), "long_wire_travelling": _long_wire,
    "leaky_wave_line_source": _leaky,
    "bowtie": _bowtie, "archimedean_spiral": _spiral, "equiangular_spiral": lambda f, v: _spiral(f, v, "equi"),
    "discone": _discone, "conical_monopole": _cone_monopole, "planar_monopole_circular": _planar_monopole,
    "planar_monopole_rectangular": lambda f, v: _planar_monopole(f, v, "rect"), "vivaldi_tsa": _vivaldi,
    "prime_focus_parabolic": _dish, "offset_parabolic": lambda f, v: _dish(f, v, "offset"),
    "cassegrain": lambda f, v: _dish(f, v, "cass"), "gregorian_dual_reflector": lambda f, v: _dish(f, v, "greg"),
    "cylindrical_parabolic": lambda f, v: _dish(f, v, "cyl"),
    "corner_reflector_90": lambda f, v: _corner(f, v, 90), "corner_reflector_60": lambda f, v: _corner(f, v, 60),
    "half_wave_slot": _slot, "folded_slot": lambda f, v: _slot(f, v, "folded"), "cavity_backed_slot": _cavity_slot,
    "waveguide_longitudinal_slot": lambda f, v: _guide_slots(f, v, "single"),
    "waveguide_slot_array_resonant": lambda f, v: _guide_slots(f, v, "resonant"),
    "waveguide_slot_array_travelling_wave": lambda f, v: _guide_slots(f, v, "travelling"),
    "luneburg_lens": _luneburg, "hyperbolic_dielectric_lens": _dielectric_lens,
    "metal_plate_lens": lambda f, v: _dielectric_lens(f, v, "plates"), "fresnel_zone_plate": _zone_plate,
    "cylindrical_dra": _dra, "hemispherical_dra": lambda f, v: _dra(f, v, "hemi"),
    "rectangular_dra": lambda f, v: _dra(f, v, "rect"),
}
FAMILY_FALLBACK = {
    "wire": _dipole, "loop": _circle_loop, "patch": _rect_patch, "horn": _pyramidal,
    "travelling_wave": _yagi, "uwb": _bowtie, "reflector": _dish, "slot": _slot, "lens": _luneburg,
    "dielectric": _dra,
}


def paint(painter: QPainter, rect: QRectF, key: str, family: str, values: Mapping, labels: bool = True,
          grid: bool = True) -> list[str]:
    """Draw one archetype's figure into rect; returns the labels it painted."""
    painter.save()
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
    painter.setClipRect(rect)
    if grid:
        painter.setPen(QPen(GRID, 1))
        step = max(12.0, rect.height() / 22.0)
        x = rect.left() + step
        while x < rect.right():
            painter.drawLine(QPointF(x, rect.top()), QPointF(x, rect.bottom()))
            x += step
        y = rect.top() + step
        while y < rect.bottom():
            painter.drawLine(QPointF(rect.left(), y), QPointF(rect.right(), y))
            y += step
    fig = Fig(painter, rect, labels)
    draw = DRAWINGS.get(key) or FAMILY_FALLBACK.get(family, _dipole)
    try:
        draw(fig, values)
    except Exception:  # noqa: BLE001 - a figure must never take the window down
        fig = Fig(painter, rect, False)
        FAMILY_FALLBACK.get(family, _dipole)(fig, {})
    painter.restore()
    return fig.drawn


def thumbnail(key: str, family: str, values: Mapping, width: int, height: int, dpr: float = 2.0) -> QPixmap:
    """The unlabelled figure on a tint of its family's colour, for the gallery."""
    img = QImage(int(width * dpr), int(height * dpr), QImage.Format.Format_ARGB32_Premultiplied)
    img.setDevicePixelRatio(dpr)
    tint = family_colour(family)
    bg = QColor(tint)
    bg.setAlpha(26)
    img.fill(QColor("#ffffff"))
    p = QPainter(img)
    p.fillRect(QRectF(0, 0, width, height), bg)
    paint(p, QRectF(0, 0, width, height), key, family, values, labels=False, grid=False)
    p.end()
    return QPixmap.fromImage(img)
