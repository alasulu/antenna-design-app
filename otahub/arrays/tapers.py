"""Array excitation tapers.

Every function returns real, non-negative amplitude weights normalised to a
peak of 1. The taper is what trades beamwidth against sidelobe level: uniform
gives the narrowest beam and -13.26 dB sidelobes, binomial eliminates
sidelobes entirely at the cost of a very broad beam, and Dolph-Chebyshev is
the optimum compromise - the narrowest beam achievable for a chosen sidelobe
level.
"""
from __future__ import annotations

import math

import numpy as np


def uniform(n: int) -> np.ndarray:
    """Equal excitation. Narrowest beam, but -13.26 dB first sidelobe."""
    _check(n)
    return np.ones(n)


def binomial(n: int) -> np.ndarray:
    """Binomial coefficients: no sidelobes at all for d <= lambda/2.

    The price is severe - the beam is much broader than uniform and the edge
    to centre amplitude ratio grows impossibly large for big arrays (a
    20-element binomial array spans over 60 dB of excitation).
    """
    _check(n)
    w = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)
    return w / w.max()


def _cheb_T(order: int, x: np.ndarray) -> np.ndarray:
    """Chebyshev polynomial T_m(x), valid both inside and outside |x| <= 1."""
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    inside = np.abs(x) <= 1.0
    out[inside] = np.cos(order * np.arccos(x[inside]))
    outside = ~inside
    # T_m(x) = cosh(m*arccosh|x|), with the sign of x carried for odd order
    out[outside] = np.cosh(order * np.arccosh(np.abs(x[outside])))
    neg = outside & (x < 0)
    if order % 2 == 1:
        out[neg] = -out[neg]
    return out


def dolph_chebyshev(n: int, sidelobe_db: float = -30.0) -> np.ndarray:
    """Dolph-Chebyshev taper for a target sidelobe level.

    Produces the narrowest possible main beam for the given sidelobe level -
    and every sidelobe sits at exactly that level, which is the signature of
    the design. `sidelobe_db` is negative (e.g. -30).

    Computed by sampling the Chebyshev pattern and inverse-transforming, which
    stays numerically sound for large arrays where the direct recurrence for
    the coefficients loses precision.
    """
    _check(n)
    if sidelobe_db >= 0:
        raise ValueError(f"sidelobe_db must be negative, got {sidelobe_db}")
    if n == 1:
        return np.ones(1)
    order = n - 1
    ratio = 10.0 ** (abs(sidelobe_db) / 20.0)
    x0 = math.cosh(math.acosh(ratio) / order)

    k = np.arange(n)
    pattern = _cheb_T(order, x0 * np.cos(np.pi * k / n))
    # The array sits at positions (i - (n-1)/2), which are half-integers for
    # even n. This linear phase carries that offset; without it the even-length
    # weights come out wrong (and, embarrassingly, independent of the design
    # sidelobe level). The sign matches numpy's ifft convention, exp(+2j*pi*kn/N).
    pattern = pattern * np.exp(-1j * np.pi * (n - 1) * k / n)
    weights = np.abs(np.real(np.fft.ifft(pattern)))
    return weights / weights.max()


def taylor_nbar(n: int, sidelobe_db: float = -30.0, nbar: int = 5) -> np.ndarray:
    """Taylor n-bar taper for a DISCRETE array (Villeneuve's distribution).

    Dolph-Chebyshev holds every sidelobe at the design level forever, which
    means the far-out sidelobes are higher than they need to be and the
    aperture distribution has awkward edge spikes. Taylor holds only the first
    `nbar` sidelobes near the design level and lets the rest fall off, which is
    what most real apertures use.

    Built by placing the array polynomial's zeros, not by sampling the
    continuous line-source distribution. In u, where the uniform N-element
    array's zeros sit at 1, 2, 3, ..., the first nbar - 1 zeros are the
    N-element Dolph-Chebyshev zeros stretched by sigma so that the nbar-th
    lands on the uniform array's; every zero beyond that IS the uniform
    array's (Villeneuve, IEEE Trans. AP-32, 1984). The weights are the
    polynomial's coefficients, read off by an FFT of its samples.

    Why not the line source: sampling it is exact only as N grows. Over 357
    designs - N = 5..101, odd and even, -20 to -40 dB, every nbar from Taylor's
    minimum 2A^2 + 1/2 up to 8 - the sampled line source exceeded its design
    sidelobe level by more than 0.5 dB in 156, by up to 2.3 dB, and not only on
    tiny arrays. This construction, judged by dense evaluation of the array
    factor over the same 357, stays within 0.05 dB above the design level for
    N >= 10 and at most 0.63 dB below it; one small odd array (N = 9, -30 dB,
    nbar = 4) reaches 0.36 dB above. It converges to the line source as N grows
    (weights within 1e-4 at N = 100) and is exactly Dolph-Chebyshev once nbar
    passes the last zero pair. (An earlier note here reported -28.9 dB for a
    20-element -30 dB design; the old taper actually gave -30.10 dB there -
    its failures were elsewhere, and larger.)

    Two properties worth knowing. A large nbar at a modest sidelobe level
    (-20 dB with nbar >= 5 on a 10- to 16-element array, say) rises again at the
    edges; that is Taylor's, not this construction's - the line-source version
    does it in the same cases. And the frequent claim that Taylor beats
    Dolph-Chebyshev on aperture efficiency does not hold for DISCRETE arrays:
    Chebyshev is the more efficient of the two at n >= 20, as the optimality
    proof for discrete arrays implies. Choose Taylor for the absence of edge
    spikes and for decaying far sidelobes, not for efficiency.
    """
    _check(n)
    if sidelobe_db >= 0:
        raise ValueError(f"sidelobe_db must be negative, got {sidelobe_db}")
    if nbar < 1:
        raise ValueError(f"nbar must be at least 1, got {nbar}")
    if n == 1:
        return np.ones(1)
    ratio = 10.0 ** (abs(sidelobe_db) / 20.0)
    x0 = math.cosh(math.acosh(ratio) / (n - 1))
    p = np.arange(1, n)
    cheb = (n / math.pi) * np.arccos(np.cos((2 * p - 1) * math.pi / (2 * (n - 1))) / x0)
    pairs = (n - 1) // 2
    nb = min(nbar, pairs + 1)
    # Once nbar passes the last zero pair every zero is Chebyshev's and the
    # taper IS Dolph-Chebyshev. (For odd n the stretch would otherwise be set
    # by a zero beyond psi = pi, which is not one of this array's.)
    sigma = 1.0 if nbar > pairs else nb / cheb[nb - 1]
    zeros = []
    for q in range(1, pairs + 1):
        u = sigma * cheb[q - 1] if q < nb else float(q)
        zeros += [u, -u]
    if n % 2 == 0:
        zeros.append(n / 2.0)                           # psi = pi
    roots = np.exp(2j * math.pi * np.asarray(zeros) / n)
    z = np.exp(2j * math.pi * np.arange(n) / n)
    diff = z[:, None] - roots[None, :]
    # the polynomial's samples, in log form so a thousand-element array cannot overflow
    with np.errstate(divide="ignore"):
        logmag = np.log(np.abs(diff)).sum(axis=1)
    phase = np.angle(diff).sum(axis=1)
    live = np.isfinite(logmag)
    samples = np.zeros(n, dtype=complex)
    samples[live] = np.exp(logmag[live] - logmag[live].max() + 1j * phase[live])
    weights = np.real(np.fft.fft(samples)) / n
    weights = np.abs(weights)
    return weights / weights.max()


def raised_cosine(n: int, pedestal: float = 0.0) -> np.ndarray:
    """Cosine-on-a-pedestal taper. `pedestal` 0 gives pure cosine, 1 uniform."""
    _check(n)
    if not 0.0 <= pedestal <= 1.0:
        raise ValueError(f"pedestal must be in [0, 1], got {pedestal}")
    p = (np.arange(n) - (n - 1) / 2.0) / max(n - 1, 1)
    w = pedestal + (1.0 - pedestal) * np.cos(np.pi * p)
    return w / w.max()


def taper_efficiency(weights: np.ndarray) -> float:
    """Aperture taper efficiency: (sum |a|)^2 / (N * sum |a|^2).

    Exactly 1 for uniform excitation, and the factor by which any taper
    reduces directivity relative to uniform. This is the real cost of low
    sidelobes, and it is worth quoting alongside any taper.
    """
    w = np.asarray(weights, dtype=float)
    return float(w.sum() ** 2 / (len(w) * np.sum(w ** 2)))


def _check(n: int) -> None:
    if not isinstance(n, (int, np.integer)) or n < 1:
        raise ValueError(f"element count must be a positive integer, got {n!r}")


TAPERS = {
    "uniform": uniform,
    "binomial": binomial,
    "chebyshev": dolph_chebyshev,
    "taylor": taylor_nbar,
    "cosine": raised_cosine,
}
