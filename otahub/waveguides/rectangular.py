"""Rectangular metallic waveguide.

Mode cutoffs, propagation constants, wave impedance and attenuation. All SI.
The convention throughout is a >= b, with a the broad wall, so the dominant
mode is TE10 and its cutoff is c/(2a).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterator

from ..core.constants import C0, ETA0, surface_resistance


@dataclass(frozen=True, slots=True)
class Mode:
    """A waveguide mode and its cutoff frequency."""
    family: str          # "TE" or "TM"
    m: int
    n: int
    f_cutoff_hz: float

    @property
    def label(self) -> str:
        return f"{self.family}{self.m}{self.n}"

    def __str__(self) -> str:
        return f"{self.label} @ {self.f_cutoff_hz / 1e9:.4f} GHz"


class RectangularWaveguide:
    """Air- or dielectric-filled rectangular waveguide.

    >>> wr90 = RectangularWaveguide(0.02286, 0.01016)
    >>> round(wr90.cutoff(1, 0) / 1e9, 4)
    6.5566
    """

    def __init__(self, a: float, b: float, eps_r: float = 1.0, mu_r: float = 1.0,
                 sigma: float = 5.8e7, tan_delta: float = 0.0,
                 wall_mu_r: float = 1.0) -> None:
        if a <= 0 or b <= 0:
            raise ValueError(f"guide dimensions must be positive, got a={a}, b={b}")
        if b > a:
            raise ValueError(
                f"convention is a >= b (a is the broad wall); got a={a}, b={b}. "
                f"Swap them so the dominant mode is TE10."
            )
        self.a, self.b = a, b
        self.eps_r, self.mu_r = eps_r, mu_r
        self.sigma, self.tan_delta = sigma, tan_delta
        self.wall_mu_r = wall_mu_r      # the metal's, not the filling's (mu_r)

    # ------------------------------------------------------------- geometry

    @property
    def aspect_ratio(self) -> float:
        return self.a / self.b

    @property
    def _v(self) -> float:
        """Phase velocity of the filling medium."""
        return C0 / math.sqrt(self.eps_r * self.mu_r)

    @property
    def eta(self) -> float:
        """Intrinsic impedance of the filling medium."""
        return ETA0 * math.sqrt(self.mu_r / self.eps_r)

    # --------------------------------------------------------------- modes

    def cutoff(self, m: int, n: int) -> float:
        """Cutoff frequency [Hz] of the (m, n) mode.

        Identical for TE_mn and TM_mn - cutoff depends only on the transverse
        geometry. TM modes additionally require m, n >= 1 to exist at all.
        """
        if m < 0 or n < 0 or (m == 0 and n == 0):
            raise ValueError(f"invalid mode indices ({m}, {n})")
        return (self._v / 2.0) * math.hypot(m / self.a, n / self.b)

    def modes(self, f_max_hz: float, max_index: int = 6) -> list[Mode]:
        """Every mode with a cutoff below `f_max_hz`, ordered by cutoff."""
        found: list[Mode] = []
        for m in range(max_index + 1):
            for n in range(max_index + 1):
                if m == 0 and n == 0:
                    continue
                fc = self.cutoff(m, n)
                if fc > f_max_hz:
                    continue
                found.append(Mode("TE", m, n, fc))
                if m >= 1 and n >= 1:          # TM_mn needs both indices nonzero
                    found.append(Mode("TM", m, n, fc))
        return sorted(found, key=lambda mo: (mo.f_cutoff_hz, mo.family, mo.m, mo.n))

    @property
    def dominant_cutoff_hz(self) -> float:
        return self.cutoff(1, 0)

    @property
    def single_mode_band_hz(self) -> tuple[float, float]:
        """Frequency range over which only TE10 propagates.

        The upper edge is the next mode up, which is TE20 (cutoff 2 x TE10)
        when a >= 2b, and TE01 (cutoff v/2b) when the guide is squatter than
        that. Standard guides use a = 2b exactly, making the two coincide.
        """
        lo = self.cutoff(1, 0)
        hi = min(self.cutoff(2, 0), self.cutoff(0, 1))
        return lo, hi

    def is_single_mode(self, f_hz: float) -> bool:
        lo, hi = self.single_mode_band_hz
        return lo < f_hz < hi

    # --------------------------------------------------------- propagation

    def k(self, f_hz: float) -> float:
        """Wavenumber of an unbounded wave in the filling medium [rad/m]."""
        return 2.0 * math.pi * f_hz / self._v

    def beta(self, f_hz: float, m: int = 1, n: int = 0) -> float:
        """Propagation constant [rad/m]. Zero at and below cutoff."""
        fc = self.cutoff(m, n)
        if f_hz <= fc:
            return 0.0
        return self.k(f_hz) * math.sqrt(1.0 - (fc / f_hz) ** 2)

    def guide_wavelength(self, f_hz: float, m: int = 1, n: int = 0) -> float:
        """Guide wavelength [m]; diverges at cutoff."""
        beta = self.beta(f_hz, m, n)
        if beta == 0.0:
            return math.inf
        return 2.0 * math.pi / beta

    def wave_impedance(self, f_hz: float, m: int = 1, n: int = 0,
                       family: str = "TE") -> float:
        """Wave impedance [ohm]. TE rises above eta near cutoff; TM falls to zero."""
        fc = self.cutoff(m, n)
        if f_hz <= fc:
            return math.inf if family.upper() == "TE" else 0.0
        ratio = math.sqrt(1.0 - (fc / f_hz) ** 2)
        return self.eta / ratio if family.upper() == "TE" else self.eta * ratio

    # --------------------------------------------------------- attenuation

    def conductor_attenuation(self, f_hz: float) -> float:
        """TE10 conductor loss [Np/m].

        alpha_c = Rs / (b*eta*sqrt(1-(fc/f)^2)) * (1 + (2b/a)*(fc/f)^2)
        """
        fc = self.cutoff(1, 0)
        if f_hz <= fc:
            return math.inf
        if self.sigma == math.inf:
            return 0.0
        rs = surface_resistance(f_hz, self.sigma, self.wall_mu_r)
        ratio = (fc / f_hz) ** 2
        return (rs / (self.b * self.eta * math.sqrt(1.0 - ratio))) * (
            1.0 + (2.0 * self.b / self.a) * ratio)

    def dielectric_attenuation(self, f_hz: float) -> float:
        """Dielectric loss [Np/m]: alpha_d = k^2 * tan_delta / (2*beta)."""
        if self.tan_delta == 0.0:
            return 0.0
        beta = self.beta(f_hz)
        if beta == 0.0:
            return math.inf
        return self.k(f_hz) ** 2 * self.tan_delta / (2.0 * beta)

    def attenuation_db_per_m(self, f_hz: float) -> float:
        """Total TE10 attenuation [dB/m]. 1 Np = 8.685889638 dB."""
        total = self.conductor_attenuation(f_hz) + self.dielectric_attenuation(f_hz)
        return total * 8.685889638065035

    def power_capacity_w(self, f_hz: float, e_max_v_per_m: float = 3.0e6) -> float:
        """TE10 breakdown power [W] for a peak field limit.

        P = a*b*E0^2 / (4*eta) * sqrt(1 - (fc/f)^2). The 3 MV/m default is dry
        air at sea level; derate hard for humidity, altitude and surface finish.
        """
        fc = self.cutoff(1, 0)
        if f_hz <= fc:
            return 0.0
        return (self.a * self.b * e_max_v_per_m ** 2 / (4.0 * self.eta)) * math.sqrt(
            1.0 - (fc / f_hz) ** 2)

    def __repr__(self) -> str:
        return (f"RectangularWaveguide(a={self.a * 1e3:.3f}mm, b={self.b * 1e3:.3f}mm, "
                f"TE10={self.dominant_cutoff_hz / 1e9:.4f}GHz)")


#: Standard WR-series guides: designation -> (a, b) in metres.
#: Stored as EXACT conversions of the EIA nominal inch dimensions (1 in =
#: 25.4 mm exactly), because rounding the millimetre values shifts the aspect
#: ratio enough to matter - WR-28 lands on 1.997 instead of 2.000.
#: The WR number is the broad wall in hundredths of an inch (WR-90 = 0.90 in).
_WR_INCHES: dict[str, tuple[float, float]] = {
    "WR-2300": (23.000, 11.500), "WR-1500": (15.000, 7.500),
    "WR-975": (9.750, 4.875), "WR-650": (6.500, 3.250),
    "WR-430": (4.300, 2.150), "WR-284": (2.840, 1.340),
    "WR-187": (1.872, 0.872), "WR-137": (1.372, 0.622),
    "WR-112": (1.122, 0.497), "WR-90": (0.900, 0.400),
    "WR-62": (0.622, 0.311), "WR-42": (0.420, 0.170),
    "WR-28": (0.280, 0.140), "WR-22": (0.224, 0.112),
    "WR-15": (0.148, 0.074), "WR-10": (0.100, 0.050),
}

_INCH = 0.0254

WR_SERIES: dict[str, tuple[float, float]] = {
    name: (a * _INCH, b * _INCH) for name, (a, b) in _WR_INCHES.items()
}


def standard(name: str, **kwargs) -> RectangularWaveguide:
    """Build a standard WR-series guide by name, e.g. ``standard("WR-90")``."""
    key = name.upper().replace("_", "-")
    if not key.startswith("WR-"):
        key = f"WR-{key}"
    try:
        a, b = WR_SERIES[key]
    except KeyError:
        raise KeyError(
            f"unknown guide {name!r}; known: {', '.join(sorted(WR_SERIES))}"
        ) from None
    return RectangularWaveguide(a, b, **kwargs)


def recommended_band(name: str) -> tuple[float, float]:
    """Conventional operating band [Hz]: roughly 1.25x to 1.9x the TE10 cutoff.

    Guides are never run right down to cutoff (dispersion and loss blow up) nor
    right up to TE20 (moding), so the usable band is inset from the
    single-mode band at both ends.
    """
    guide = standard(name)
    fc = guide.dominant_cutoff_hz
    return 1.25 * fc, 1.90 * fc
