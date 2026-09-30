"""The look of the desktop app: one palette, one stylesheet, a few painted icons.

Copper is the accent because copper is what an antenna is made of; each family
carries its own colour (drawings.FAMILY_COLOURS) for chips and card edges. The
navigation rail is dark so the working area, where the drawings and numbers are,
stays the brightest thing on screen.
"""
from __future__ import annotations

from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QIcon, QPainter, QPainterPath, QPen, QPixmap

BG = "#eef1f3"
SURFACE = "#ffffff"
SUBTLE = "#f6f8f9"
LINE = "#d9dfe4"
INK = "#1f2a33"
MUTED = "#5c6b77"
FAINT = "#8a98a3"
ACCENT = "#b8672e"
ACCENT_DARK = "#9a5423"
ACCENT_SOFT = "#f4e6da"
RAIL = "#1c2730"
RAIL_TEXT = "#aebbc5"
OK = "#2d7a55"
OK_SOFT = "#e0f0e7"
WARN = "#96600a"
WARN_SOFT = "#f7ecd6"
BAD = "#a93a2e"

STYLE = f"""
QMainWindow, QWidget#page {{ background: {BG}; }}
QWidget {{ color: {INK}; }}
QToolTip {{ background: {INK}; color: #ffffff; border: none; padding: 6px 8px; }}

/* navigation rail */
QListWidget#rail {{ background: {RAIL}; border: none; outline: 0; padding-top: 10px; }}
QListWidget#rail::item {{ color: {RAIL_TEXT}; padding: 12px 10px; margin: 2px 8px; border-radius: 8px; }}
QListWidget#rail::item:hover {{ background: #26333d; color: #ffffff; }}
QListWidget#rail::item:selected {{ background: #2f3e4a; color: #ffffff; border-left: 3px solid {ACCENT}; }}
QLabel#brand {{ color: #ffffff; font-size: 17px; font-weight: 700; padding: 18px 16px 2px 16px; background: {RAIL}; }}
QLabel#brandsub {{ color: {RAIL_TEXT}; font-size: 11px; padding: 0 16px 10px 16px; background: {RAIL}; }}

/* type */
QLabel#h1 {{ font-size: 24px; font-weight: 700; }}
QLabel#h2 {{ font-size: 15px; font-weight: 650; }}
QLabel#muted {{ color: {MUTED}; }}
QLabel#faint {{ color: {FAINT}; font-size: 11px; }}
QLabel#fieldname {{ font-weight: 550; }}
QLabel#fieldmeta {{ color: {FAINT}; font-size: 11px; }}
QLabel#error {{ color: {BAD}; font-size: 11px; }}

/* cards */
QFrame#card {{ background: {SURFACE}; border: 1px solid {LINE}; border-radius: 12px; }}
QFrame#tile {{ background: {SURFACE}; border: 1px solid {LINE}; border-radius: 12px; }}
QLabel#tilelabel {{ color: {MUTED}; font-size: 11px; font-weight: 600; letter-spacing: 0.5px; }}
QLabel#tilevalue {{ font-size: 22px; font-weight: 650; }}
QLabel#tilekey {{ color: {FAINT}; font-size: 10px; }}
QFrame#banner {{ background: {WARN_SOFT}; border: 1px solid #ead7b0; border-radius: 10px; }}
QFrame#banner QLabel {{ color: #5b3d06; }}

/* controls */
QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {{
    background: {SURFACE}; border: 1px solid #cbd3d9; border-radius: 8px; padding: 6px 9px; min-height: 18px;
    selection-background-color: {ACCENT_SOFT}; selection-color: {INK}; }}
QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus {{ border: 1px solid {ACCENT}; }}
QLineEdit[invalid="true"] {{ border: 1px solid {BAD}; background: #fbeceb; }}
QLineEdit#search {{ border-radius: 18px; padding: 7px 14px; min-width: 240px; }}
QComboBox::drop-down {{ border: none; width: 22px; }}
QPushButton {{ background: {SURFACE}; border: 1px solid #cbd3d9; border-radius: 8px; padding: 7px 14px; }}
QPushButton:hover {{ border-color: {ACCENT}; color: {ACCENT_DARK}; }}
QPushButton#primary {{ background: {ACCENT}; color: #ffffff; border: none; font-weight: 600; }}
QPushButton#primary:hover {{ background: {ACCENT_DARK}; color: #ffffff; }}
QPushButton#back {{ background: transparent; border: none; color: {ACCENT_DARK}; font-weight: 600; padding: 4px 2px; }}
QPushButton#link {{ background: transparent; border: none; color: {ACCENT_DARK}; padding: 2px 0; text-align: left; }}
QPushButton#chip {{ border-radius: 14px; padding: 5px 12px; background: {SURFACE}; color: {MUTED}; }}
QPushButton#chip:checked {{ background: {INK}; color: #ffffff; border-color: {INK}; }}
QCheckBox {{ spacing: 8px; }}

/* tables and tabs */
QTableWidget {{ background: {SURFACE}; border: none; gridline-color: transparent;
    alternate-background-color: {SUBTLE}; selection-background-color: {ACCENT_SOFT}; selection-color: {INK}; }}
QTableWidget::item {{ padding: 6px 8px; border: none; }}
QHeaderView::section {{ background: {SURFACE}; color: {FAINT}; border: none; border-bottom: 1px solid {LINE};
    padding: 6px 8px; font-size: 11px; font-weight: 600; }}
QTabWidget::pane {{ border: none; background: transparent; }}
QTabBar::tab {{ background: transparent; color: {MUTED}; padding: 8px 14px; margin-right: 4px;
    border-bottom: 2px solid transparent; font-weight: 600; }}
QTabBar::tab:selected {{ color: {INK}; border-bottom: 2px solid {ACCENT}; }}
QTabBar::tab:hover {{ color: {INK}; }}
QGroupBox {{ background: {SURFACE}; border: 1px solid {LINE}; border-radius: 12px; margin-top: 14px; padding: 16px 12px 10px 12px;
    font-weight: 650; }}
QGroupBox::title {{ subcontrol-origin: margin; left: 12px; padding: 0 4px; color: {INK}; }}
QPlainTextEdit, QTextBrowser {{ background: {SURFACE}; border: none; }}
QScrollArea {{ background: transparent; border: none; }}
QScrollArea > QWidget > QWidget {{ background: transparent; }}
QSplitter::handle {{ background: transparent; width: 10px; }}
QStatusBar {{ background: {SURFACE}; color: {MUTED}; border-top: 1px solid {LINE}; }}
QScrollBar:vertical {{ background: transparent; width: 10px; margin: 2px; }}
QScrollBar::handle:vertical {{ background: #c3ccd3; border-radius: 4px; min-height: 30px; }}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{ height: 0; }}
"""


def _pen(colour: str, width: float = 1.7) -> QPen:
    pen = QPen(QColor(colour), width)
    pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    return pen


def rail_icon(kind: str, colour: str = "#e6ecf0", size: int = 22) -> QIcon:
    """Line icons for the navigation rail, painted rather than shipped."""
    pm = QPixmap(size * 2, size * 2)
    pm.setDevicePixelRatio(2.0)
    pm.fill(Qt.GlobalColor.transparent)
    p = QPainter(pm)
    p.setRenderHint(QPainter.RenderHint.Antialiasing)
    p.setPen(_pen(colour))
    s = size
    if kind == "antennas":
        p.drawLine(QPointF(s / 2, s * 0.35), QPointF(s / 2, s * 0.92))
        for r in (0.18, 0.32):
            rect = QRectF(s / 2 - s * r, s * 0.35 - s * r, 2 * s * r, 2 * s * r)
            p.drawArc(rect, 30 * 16, 120 * 16)
        p.setBrush(QColor(colour))
        p.drawEllipse(QPointF(s / 2, s * 0.35), 1.8, 1.8)
    elif kind == "linear":
        p.setBrush(QColor(colour))
        for i in range(5):
            p.drawEllipse(QPointF(s * (0.14 + 0.18 * i), s * 0.8), 1.7, 1.7)
        p.setBrush(Qt.BrushStyle.NoBrush)
        path = QPainterPath(QPointF(s * 0.2, s * 0.7))
        path.cubicTo(QPointF(s * 0.3, s * 0.05), QPointF(s * 0.7, s * 0.05), QPointF(s * 0.8, s * 0.7))
        p.drawPath(path)
    elif kind == "planar":
        p.setBrush(QColor(colour))
        for i in range(4):
            for j in range(4):
                p.drawEllipse(QPointF(s * (0.2 + 0.2 * i), s * (0.2 + 0.2 * j)), 1.6, 1.6)
    elif kind == "waveguide":
        p.drawRect(QRectF(s * 0.12, s * 0.38, s * 0.52, s * 0.34))
        p.drawLine(QPointF(s * 0.12, s * 0.38), QPointF(s * 0.34, s * 0.2))
        p.drawLine(QPointF(s * 0.64, s * 0.38), QPointF(s * 0.86, s * 0.2))
        p.drawLine(QPointF(s * 0.64, s * 0.72), QPointF(s * 0.86, s * 0.54))
        p.drawLine(QPointF(s * 0.34, s * 0.2), QPointF(s * 0.86, s * 0.2))
        p.drawLine(QPointF(s * 0.86, s * 0.2), QPointF(s * 0.86, s * 0.54))
    p.end()
    return QIcon(pm)


def dot_icon(colour: str, size: int = 10) -> QIcon:
    pm = QPixmap(size * 2, size * 2)
    pm.setDevicePixelRatio(2.0)
    pm.fill(Qt.GlobalColor.transparent)
    p = QPainter(pm)
    p.setRenderHint(QPainter.RenderHint.Antialiasing)
    p.setPen(Qt.PenStyle.NoPen)
    p.setBrush(QColor(colour))
    p.drawEllipse(QRectF(1, 1, size - 2, size - 2))
    p.end()
    return QIcon(pm)


def mpl_style(fig, ax, polar: bool = False) -> None:
    """Matplotlib in the app's palette: quiet axes, the accent for data."""
    fig.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    ax.tick_params(colors=MUTED, labelsize=8.5)
    for side in ("top", "right"):
        if side in ax.spines and not polar:
            ax.spines[side].set_visible(False)
    for spine in ax.spines.values():
        spine.set_color(LINE)
    ax.xaxis.label.set_color(MUTED)
    ax.yaxis.label.set_color(MUTED)
    ax.title.set_color(INK)
    ax.grid(True, color="#e3e8ec", linewidth=0.8)
