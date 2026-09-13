"""Circular metallic waveguide.

Cutoffs follow from the zeros of Bessel functions and their derivatives:
TE_mn from the n-th zero of J'_m, TM_mn from the n-th zero of J_m.
"""
from __future__ import annotations

import math
from functools import lru_cache

from scipy.special import jn_zeros, jnp_zeros

from ..core.constants import C0, ETA0, surface_resistance
from .rectangular import Mode


@lru_cache(maxsize=256)
def te_root(m: int, n: int) -> float:
    """n-th non-trivial zero of J'_m. TE11 gives 1.8412."""
    return float(jnp_zeros(m, n)[-1])


@lru_cache(maxsize=256)
def tm_root(m: int, n: int) -> float:
    """n-th zero of J_m. TM01 gives 2.4049."""
    return float(jn_zeros(m, n)[-1])


class CircularWaveguide:
    """Circular waveguide of radius `a`.

    >>> g = CircularWaveguide(0.01)
    >>> round(g.cutoff("TE", 1, 1) / 1e9, 4)
    8.7876
    """

    def __init__(self, a: float, eps_r: float = 1.0, mu_r: float = 1.0,
                 sigma: float = 5.8e7) -> None:
        if a <= 0:
            raise ValueError(f"radius must be positive, got {a}")
        self.a, self.eps_r, self.mu_r, self.sigma = a, eps_r, mu_r, sigma

    @property
    def _v(self) -> float:
        return C0 / math.sqrt(self.eps_r * self.mu_r)

    @property
    def eta(self) -> float:
        return ETA0 * math.sqrt(self.mu_r / self.eps_r)

    def cutoff(self, family: str, m: int, n: int) -> float:
        """Cutoff frequency [Hz] of TE_mn or TM_mn."""
        fam = family.upper()
        if fam == "TE":
            root = te_root(m, n)
        elif fam == "TM":
            root = tm_root(m, n)
        else:
            raise ValueError(f"family must be 'TE' or 'TM', got {family!r}")
        return root * self._v / (2.0 * math.pi * self.a)

    def modes(self, f_max_hz: float, max_m: int = 5, max_n: int = 4) -> list[Mode]:
        found: list[Mode] = []
        for fam in ("TE", "TM"):
            for m in range(max_m + 1):
                for n in range(1, max_n + 1):
                    fc = self.cutoff(fam, m, n)
                    if fc <= f_max_hz:
                        found.append(Mode(fam, m, n, fc))
        return sorted(found, key=lambda mo: (mo.f_cutoff_hz, mo.family, mo.m, mo.n))

    @property
    def dominant_cutoff_hz(self) -> float:
        """TE11 is always dominant, since 1.8412 is the smallest root."""
        return self.cutoff("TE", 1, 1)

    @property
    def single_mode_band_hz(self) -> tuple[float, float]:
        """TE11 alone propagates up to the TM01 cutoff (2.4049/1.8412 = 1.306x)."""
        return self.cutoff("TE", 1, 1), self.cutoff("TM", 0, 1)

    def beta(self, f_hz: float, family: str = "TE", m: int = 1, n: int = 1) -> float:
        fc = self.cutoff(family, m, n)
        if f_hz <= fc:
            return 0.0
        return (2.0 * math.pi * f_hz / self._v) * math.sqrt(1.0 - (fc / f_hz) ** 2)

    def guide_wavelength(self, f_hz: float, family: str = "TE",
                         m: int = 1, n: int = 1) -> float:
        beta = self.beta(f_hz, family, m, n)
        return math.inf if beta == 0.0 else 2.0 * math.pi / beta

    def wave_impedance(self, f_hz: float, family: str = "TE",
                       m: int = 1, n: int = 1) -> float:
        fc = self.cutoff(family, m, n)
        if f_hz <= fc:
            return math.inf if family.upper() == "TE" else 0.0
        ratio = math.sqrt(1.0 - (fc / f_hz) ** 2)
        return self.eta / ratio if family.upper() == "TE" else self.eta * ratio

    def te01_attenuation(self, f_hz: float) -> float:
        """TE01 conductor loss [Np/m].

        TE01 is the famous low-loss mode: its attenuation FALLS without limit
        as frequency rises, because the wall currents are purely circumferential.
        It is not the dominant mode, so it must be kept pure against mode
        conversion - the reason TE01 circular guide never displaced coax or
        fibre despite the loss advantage.
        """
        fc = self.cutoff("TE", 0, 1)
        if f_hz <= fc:
            return math.inf
        rs = surface_resistance(f_hz, self.sigma, self.mu_r)
        ratio = (fc / f_hz) ** 2
        return (rs / (self.a * self.eta)) * ratio / math.sqrt(1.0 - ratio)

    def __repr__(self) -> str:
        return (f"CircularWaveguide(a={self.a * 1e3:.3f}mm, "
                f"TE11={self.dominant_cutoff_hz / 1e9:.4f}GHz)")
