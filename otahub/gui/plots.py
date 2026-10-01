"""Matplotlib canvases embedded in Qt."""
from __future__ import annotations

import math

import matplotlib
matplotlib.use("QtAgg")

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure

from .style import ACCENT, INK, MUTED, mpl_style

SECOND = "#3b6c9c"                     # a second trace, in the horn family's steel blue
AMPLITUDE = LinearSegmentedColormap.from_list("copper_ramp", ["#f1ddca", "#d69a66", "#8d4a1c"])


class Canvas(FigureCanvasQTAgg):
    """A figure that repaints itself on demand."""

    def __init__(self, polar: bool = False, width: float = 5.0,
                 height: float = 4.0) -> None:
        self.figure = Figure(figsize=(width, height), layout="constrained")
        self.polar = polar
        self.axes = self.figure.add_subplot(111, projection="polar" if polar else None)
        self.colorbars: list = []
        super().__init__(self.figure)
        mpl_style(self.figure, self.axes, polar)

    def clear(self) -> None:
        # A colour bar from the last draw would otherwise pile up. Colorbar.remove() also
        # unregisters it from its parent axes; removing its axes alone left constrained
        # layout holding an axes with no figure, which raised on the next draw.
        while self.colorbars:
            self.colorbars.pop().remove()
        for extra in self.figure.axes[1:]:
            extra.remove()
        self.axes.clear()
        mpl_style(self.figure, self.axes, self.polar)

    def message(self, text: str) -> None:
        """Show an explanatory message instead of an empty, mysterious plot."""
        self.clear()
        self.axes.text(0.5, 0.5, text, ha="center", va="center", wrap=True,
                       transform=self.axes.transAxes, fontsize=9.5, color=MUTED)
        self.axes.set_axis_off()
        self.draw_idle()


def plot_sweep(canvas: Canvas, freqs_hz, values, metric: str,
               unit: str = "", band: tuple[float, float] | None = None) -> None:
    """Metric against frequency, with the validity band shaded."""
    canvas.clear()
    ax = canvas.axes
    ax.set_axis_on()
    ax.plot(np.asarray(freqs_hz) / 1e9, values, lw=2.0, color=ACCENT)
    ax.set_xlabel("frequency (GHz)")
    ax.set_ylabel(f"{metric} [{unit}]" if unit else metric)
    if band is not None and all(map(math.isfinite, band)):
        lo, hi = band[0] / 1e9, band[1] / 1e9
        xlo, xhi = ax.get_xlim()
        if hi > xlo and lo < xhi:
            ax.axvspan(max(lo, xlo), min(hi, xhi), alpha=0.08, color=ACCENT, lw=0)
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
    ax.plot(theta, db, lw=2.0, color=ACCENT)
    ax.fill(theta, db, color=ACCENT, alpha=0.08)
    ax.set_ylim(floor_db, 0)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    if title:
        ax.set_title(title, fontsize=10, color=INK)
    canvas.draw_idle()


def plot_hemisphere_cuts(canvas: Canvas, cuts, title: str = "",
                         floor_db: float = -40.0) -> None:
    """Pattern cuts of a planar array, plotted over the forward hemisphere only.

    `cuts` is a sequence of (label, psi_rad, u_linear) where psi is the SIGNED
    angle from the array normal, so the trace spans -90 to +90 degrees and the
    normal points up.

    The lower half of the polar plot stays empty on purpose. A planar array
    with a ground plane radiates into one hemisphere, and there is nothing to
    draw below it — unlike a dipole cut, where a blank half really did mean the
    mirroring was missing.
    """
    canvas.clear()
    ax = canvas.axes
    ax.set_axis_on()
    peak = max(float(np.asarray(u).max()) for _, _, u in cuts)
    if peak <= 0:
        canvas.message("pattern radiates no power")
        return
    for (label, psi, u), colour in zip(cuts, (ACCENT, SECOND, MUTED)):
        db = 10.0 * np.log10(np.maximum(np.asarray(u, dtype=float) / peak,
                                        10 ** (floor_db / 10.0)))
        ax.plot(np.asarray(psi, dtype=float), db, lw=1.8, label=label, color=colour)
    ax.set_ylim(floor_db, 0)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_thetamin(-90)
    ax.set_thetamax(90)
    ax.legend(loc="lower center", fontsize=8, ncol=len(cuts), frameon=False,
              bbox_to_anchor=(0.5, -0.18))
    if title:
        ax.set_title(title, fontsize=10, color=INK)
    canvas.draw_idle()


def plot_element_layout(canvas: Canvas, positions, weights=None,
                        title: str = "") -> None:
    """Element positions in wavelengths, sized by excitation amplitude."""
    canvas.clear()
    ax = canvas.axes
    ax.set_axis_on()
    pos = np.asarray(positions, dtype=float)
    if weights is None:
        sizes = np.full(len(pos), 28.0)
        colours = None
    else:
        w = np.abs(np.asarray(weights, dtype=float))
        w = w / (w.max() or 1.0)
        sizes = 8.0 + 42.0 * w
        colours = w
    sc = ax.scatter(pos[:, 0], pos[:, 1], s=sizes, c=colours if colours is not None else ACCENT,
                    cmap=AMPLITUDE if colours is not None else None, edgecolors="none")
    if colours is not None:
        canvas.colorbars.append(canvas.figure.colorbar(sc, ax=ax, label="normalised amplitude"))
    ax.set_xlabel("x (wavelengths)")
    ax.set_ylabel("y (wavelengths)")
    ax.set_aspect("equal", adjustable="datalim")
    if title:
        ax.set_title(title, fontsize=10, color=INK)
    canvas.draw_idle()
