"""Matplotlib canvases embedded in Qt."""
from __future__ import annotations

import math

import matplotlib
matplotlib.use("QtAgg")

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure


class Canvas(FigureCanvasQTAgg):
    """A figure that repaints itself on demand."""

    def __init__(self, polar: bool = False, width: float = 5.0,
                 height: float = 4.0) -> None:
        self.figure = Figure(figsize=(width, height), layout="constrained")
        self.axes = self.figure.add_subplot(111, projection="polar" if polar else None)
        super().__init__(self.figure)

    def clear(self) -> None:
        self.axes.clear()

    def message(self, text: str) -> None:
        """Show an explanatory message instead of an empty, mysterious plot."""
        self.clear()
        self.axes.text(0.5, 0.5, text, ha="center", va="center", wrap=True,
                       transform=self.axes.transAxes, fontsize=9)
        self.axes.set_axis_off()
        self.draw_idle()


def plot_sweep(canvas: Canvas, freqs_hz, values, metric: str,
               unit: str = "", band: tuple[float, float] | None = None) -> None:
    """Metric against frequency, with the validity band shaded."""
    canvas.clear()
    ax = canvas.axes
    ax.set_axis_on()
    ax.plot(np.asarray(freqs_hz) / 1e9, values, lw=1.8)
    ax.set_xlabel("frequency (GHz)")
    ax.set_ylabel(f"{metric} [{unit}]" if unit else metric)
    ax.grid(True, alpha=0.3)
    if band is not None and all(map(math.isfinite, band)):
        lo, hi = band[0] / 1e9, band[1] / 1e9
        xlo, xhi = ax.get_xlim()
        if hi > xlo and lo < xhi:
            ax.axvspan(max(lo, xlo), min(hi, xhi), alpha=0.08, color="green")
    canvas.draw_idle()


def plot_polar(canvas: Canvas, theta_rad, u_linear, title: str = "",
               floor_db: float = -40.0, mirror: bool = True) -> None:
    """Normalised power pattern in dB on a polar axis.

    `theta` spans 0..pi, but every pattern here is azimuthally symmetric about
    the antenna axis, so the cut is mirrored through 2*pi - theta to fill the
    circle. Leaving the lower half blank would read as a one-sided pattern,
    which would be wrong.
    """
    canvas.clear()
    ax = canvas.axes
    ax.set_axis_on()
    u = np.asarray(u_linear, dtype=float)
    peak = u.max()
    if peak <= 0:
        canvas.message("pattern radiates no power")
        return
    db = 10.0 * np.log10(np.maximum(u / peak, 10 ** (floor_db / 10.0)))
    theta = np.asarray(theta_rad, dtype=float)
    if mirror:
        theta = np.concatenate([theta, 2.0 * np.pi - theta[::-1]])
        db = np.concatenate([db, db[::-1]])
    ax.plot(theta, db, lw=1.6)
    ax.set_ylim(floor_db, 0)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.grid(True, alpha=0.35)
    if title:
        ax.set_title(title, fontsize=10)
    canvas.draw_idle()
