"""Radiation-pattern utilities.

Closed-form element patterns plus numerical integration of directivity and
beamwidth. These are deliberately independent of the spec files: they are
first-principles code that the spec-driven analysis formulas can be checked
against, which is what makes a disagreement meaningful.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

_EPS = 1e-12


@dataclass(slots=True)
class Pattern:
    """Radiation intensity sampled on a theta/phi grid.

    `theta` is polar angle from +z in radians (0..pi); `phi` is azimuth
    (0..2pi). `U` has shape (len(theta), len(phi)) and carries radiation
    intensity in arbitrary but consistent units.
    """
    theta: np.ndarray
    phi: np.ndarray
    U: np.ndarray

    def __post_init__(self) -> None:
        if self.U.shape != (self.theta.size, self.phi.size):
            raise ValueError(
                f"U shape {self.U.shape} does not match grid "
                f"({self.theta.size}, {self.phi.size})"
            )

    # -------------------------------------------------------------- measures

    def radiated_power(self) -> float:
        """P_rad = int int U sin(theta) dtheta dphi."""
        integrand = self.U * np.sin(self.theta)[:, None]
        over_phi = np.trapezoid(integrand, self.phi, axis=1)
        return float(np.trapezoid(over_phi, self.theta))

    def directivity(self) -> float:
        """Peak directivity D = 4*pi*U_max / P_rad."""
        p_rad = self.radiated_power()
        if p_rad <= _EPS:
            raise ValueError("pattern radiates no power; directivity undefined")
        return float(4.0 * math.pi * self.U.max() / p_rad)

    def directivity_dbi(self) -> float:
        return 10.0 * math.log10(self.directivity())

    def cut(self, phi_rad: float) -> np.ndarray:
        """The theta-cut nearest `phi_rad`, normalised to its own peak."""
        idx = int(np.argmin(np.abs(self.phi - (phi_rad % (2 * math.pi)))))
        cut = self.U[:, idx]
        peak = cut.max()
        return cut / peak if peak > _EPS else cut


def make_grid(n_theta: int = 361, n_phi: int = 181) -> tuple[np.ndarray, np.ndarray]:
    return (np.linspace(0.0, math.pi, n_theta), np.linspace(0.0, 2.0 * math.pi, n_phi))


# ------------------------------------------------------------ element patterns

def _broadcast(theta: np.ndarray, phi: np.ndarray, f_theta: np.ndarray) -> Pattern:
    """Build an azimuthally symmetric pattern from a theta-only profile."""
    return Pattern(theta, phi, np.repeat(f_theta[:, None], phi.size, axis=1))


def short_dipole(n_theta: int = 361, n_phi: int = 181) -> Pattern:
    """Hertzian dipole along z: U ~ sin^2(theta). Exact directivity 1.5."""
    theta, phi = make_grid(n_theta, n_phi)
    return _broadcast(theta, phi, np.sin(theta) ** 2)


def finite_dipole(length_over_lambda: float = 0.5,
                  n_theta: int = 721, n_phi: int = 181) -> Pattern:
    """Centre-fed thin dipole of electrical length `length_over_lambda`.

    U ~ |[cos(k L cos(theta)/2) - cos(k L / 2)] / sin(theta)|^2, with the
    removable singularity at theta = 0, pi handled explicitly.
    """
    theta, phi = make_grid(n_theta, n_phi)
    kl_2 = math.pi * length_over_lambda          # k*L/2 with k = 2*pi/lambda
    sin_t = np.sin(theta)
    num = np.cos(kl_2 * np.cos(theta)) - math.cos(kl_2)
    f = np.zeros_like(theta)
    good = np.abs(sin_t) > 1e-9
    f[good] = num[good] / sin_t[good]
    return _broadcast(theta, phi, f ** 2)


def cosine_q(q: float = 1.0, n_theta: int = 361, n_phi: int = 181) -> Pattern:
    """cos^q(theta) *field* taper over the forward hemisphere.

    Radiation intensity is therefore cos^(2q)(theta), giving the standard
    feed-pattern directivity D = 2*(2q + 1).
    """
    theta, phi = make_grid(n_theta, n_phi)
    f = np.where(theta <= math.pi / 2, np.cos(theta), 0.0)
    return _broadcast(theta, phi, f ** (2.0 * q))


def uniform_line_source(length_over_lambda: float, n_theta: int = 2001,
                        n_phi: int = 181) -> Pattern:
    """Uniformly illuminated broadside line source along z.

    U ~ sinc^2((L/lambda) * cos(theta)), so the beam peaks broadside at
    theta = 90 deg. Classic results: first sidelobe -13.26 dB, and
    D -> 2*L/lambda for L >> lambda.
    """
    theta, phi = make_grid(n_theta, n_phi)
    f = np.sinc(length_over_lambda * np.cos(theta))   # numpy sinc = sin(pi x)/(pi x)
    return _broadcast(theta, phi, f ** 2)


# --------------------------------------------------------------- beam measures

def hpbw(pattern: Pattern, phi_rad: float = 0.0) -> float:
    """Half-power beamwidth [rad] of the main lobe in the given phi-cut."""
    cut = pattern.cut(phi_rad)
    theta = pattern.theta
    peak_idx = int(np.argmax(cut))
    half = 0.5

    def _crossing(indices) -> float | None:
        prev_i = peak_idx
        for i in indices:
            if cut[i] < half:
                t0, t1 = theta[prev_i], theta[i]
                c0, c1 = cut[prev_i], cut[i]
                if abs(c0 - c1) < _EPS:
                    return float(t1)
                return float(t0 + (c0 - half) * (t1 - t0) / (c0 - c1))
            prev_i = i
        return None

    lower = _crossing(range(peak_idx, -1, -1))
    upper = _crossing(range(peak_idx, theta.size))
    if lower is None or upper is None:
        raise ValueError("main lobe does not fall to half power inside the grid")
    return upper - lower


def hpbw_deg(pattern: Pattern, phi_rad: float = 0.0) -> float:
    return math.degrees(hpbw(pattern, phi_rad))


def _main_lobe_bounds(cut: np.ndarray, peak_idx: int) -> tuple[int, int]:
    """Index span of the main lobe, bounded by the first null on each side.

    A null is where the pattern stops falling and starts rising again; if it
    never turns, the lobe runs to the edge of the grid.
    """
    def walk(step: int) -> int:
        i = peak_idx
        while 0 <= i + step < cut.size:
            if cut[i + step] > cut[i] + _EPS:   # turned upward: null at i
                return i
            i += step
        return i
    return walk(-1), walk(+1)


def first_sidelobe_db(pattern: Pattern, phi_rad: float = 0.0) -> float:
    """Highest sidelobe relative to the main beam [dB, negative].

    The main lobe is excised between its bounding nulls before searching, so a
    pattern with a second main beam elsewhere in the cut does not masquerade
    as a 0 dB sidelobe.
    """
    cut = pattern.cut(phi_rad)
    peak_idx = int(np.argmax(cut))
    peak = cut[peak_idx]
    if peak <= _EPS:
        return float("-inf")
    lo, hi = _main_lobe_bounds(cut, peak_idx)
    mask = np.ones(cut.size, dtype=bool)
    mask[lo:hi + 1] = False
    if not mask.any():
        return float("-inf")
    sidelobe = cut[mask].max()
    if sidelobe <= _EPS:
        return float("-inf")
    return float(10.0 * math.log10(sidelobe / peak))
