"""Dielectric resonators on a ground plane - an arbiter for the DRA specs.

Two independent routes to a resonator's complex natural frequency, from which
the resonance is the real part and the radiation Q is Re/(2 Im):

``sphere_pole``  the exact magnetic-dipole (TE1) pole of a dielectric sphere,
                 which by image theory is the hemispherical DRA's TE111 mode.
                 Elementary functions and a complex Newton iteration.

``brick`` and ``hemisphere``
                 a small 3-D FDTD solver: Yee grid, one quarter of the upper
                 half space, CPML on the outer faces, a pulse to ring the
                 resonator and a matrix pencil fitted to the ringdown. The
                 symmetry planes are chosen so that only modes whose magnetic
                 dipole lies along x, parallel to the ground, can exist:
                 x = 0 is a magnetic wall, y = 0 and z = 0 (the ground) are
                 electric walls. The dielectric-waveguide model and Mie theory
                 enter nowhere.

The FDTD is checked against ``sphere_pole`` on a staircased hemisphere before
it is trusted with a brick. Grid units throughout: cell size 1, c = 1.
"""
from __future__ import annotations

import math

import numpy as np

__all__ = ["sphere_pole", "sphere_resonance", "matrix_pencil", "ringdown", "brick", "hemisphere"]


# --------------------------------------------------------------- exact sphere

def _psi(z):
    return np.sin(z) / z - np.cos(z)


def _dpsi(z):
    return np.cos(z) / z - np.sin(z) / z ** 2 + np.sin(z)


def _xi(z):      # z*h1^(2)(z): outgoing for exp(+j w t)
    return _psi(z) - 1j * (-np.cos(z) / z - np.sin(z))


def _dxi(z):
    return _dpsi(z) - 1j * (np.sin(z) / z + np.cos(z) / z ** 2 - np.cos(z))


def _te1(x: complex, m: float) -> complex:
    """Vanishes at a TE1 (magnetic dipole) natural frequency x = k0*a."""
    return m * _dpsi(m * x) * _xi(x) - _dxi(x) * _psi(m * x)


def sphere_pole(eps_r: float, guess: complex | None = None) -> complex:
    """Complex k0*a of the lowest TE1 mode of a dielectric sphere of radius a.

    Im > 0 is decay for exp(+j w t). As eps_r grows, sqrt(eps_r)*Re -> pi.
    """
    m = math.sqrt(eps_r)
    if guess is None:
        x = math.pi / m * 0.95
        guess = complex(x, x / (2.0 * (0.8 * eps_r + 0.012 * eps_r ** 2)))
    x = complex(guess)
    for _ in range(100):
        h = 1e-7 * abs(x)
        f = _te1(x, m)
        df = (_te1(x + h, m) - _te1(x - h, m)) / (2 * h)
        step = f / df
        x -= step
        if abs(step) < 1e-15 * abs(x):
            break
    else:  # pragma: no cover
        raise RuntimeError("sphere pole did not converge")
    if not (0.8 < m * x.real / math.pi < 1.05 and x.imag > 0):
        raise RuntimeError(f"converged on the wrong root: {x}")
    return x


def sphere_resonance(eps_r: float) -> tuple[float, float]:
    """(k0*a, radiation Q) of the hemispherical DRA's TE111 mode."""
    x = sphere_pole(eps_r)
    return x.real, x.real / (2.0 * x.imag)


# ----------------------------------------------------------------------- FDTD

def _cpml(n: int, npml: int, half: bool, dt: float, m: int = 3, kappa_max: float = 5.0):
    """CPML coefficients at nodes 0..n-1 (shifted by 1/2 if half). The layer
    fills the last npml cells before the outer PEC wall."""
    last = n if half else n - 1
    pos = np.arange(n) + (0.5 if half else 0.0)
    rho = np.clip((pos - (last - npml)) / npml, 0.0, 1.0)
    sig = 0.8 * (m + 1) * rho ** m
    kap = 1.0 + (kappa_max - 1.0) * rho ** m
    b = np.exp(-sig / kap * dt)
    c = np.where(sig > 0, (b - 1.0) / kap, 0.0)
    start = int(np.argmax(rho > 0)) if np.any(rho > 0) else n
    return kap, b, c, start


class _Yee:
    """Quarter-space Yee grid: PMC at x = 0, PEC at y = 0 and z = 0, CPML outside."""

    def __init__(self, shape, eps_cell, npml: int, dtype=np.float32):
        nx, ny, nz = shape
        self.dt = dt = 0.99 / math.sqrt(3.0)
        z = lambda *s: np.zeros(s, dtype)
        self.Ex, self.Ey, self.Ez = z(nx, ny + 1, nz + 1), z(nx + 1, ny, nz + 1), z(nx + 1, ny + 1, nz)
        self.Hx, self.Hy, self.Hz = z(nx + 1, ny, nz), z(nx, ny + 1, nz), z(nx, ny, nz + 1)
        # cell permittivity with one mirrored layer below each symmetry plane
        c = [np.arange(-1, n) + 0.5 for n in shape]
        C = np.broadcast_to(eps_cell(c[0][:, None, None], c[1][None, :, None], c[2][None, None, :]),
                            (nx + 1, ny + 1, nz + 1)).astype(float)
        # an edge sees the mean of the four cells around it: the arithmetic mean
        # a tangential field wants, and exact for a normal one on a grid-aligned face
        ex = 0.25 * (C[1:, :-1, :-1] + C[1:, 1:, :-1] + C[1:, :-1, 1:] + C[1:, 1:, 1:])
        ey = 0.25 * (C[:-1, 1:, :-1] + C[1:, 1:, :-1] + C[:-1, 1:, 1:] + C[1:, 1:, 1:])
        ez = 0.25 * (C[:-1, :-1, 1:] + C[1:, :-1, 1:] + C[:-1, 1:, 1:] + C[1:, 1:, 1:])
        self.cEx = (dt / ex[:, 1:, 1:]).astype(dtype)
        self.cEy = (dt / ey[:, :, 1:]).astype(dtype)
        self.cEz = (dt / ez[:, 1:, :]).astype(dtype)
        self.p_int = [_cpml(n + 1, npml, False, dt) for n in shape]
        self.p_half = [_cpml(n, npml, True, dt) for n in shape]
        self.psi: dict[str, np.ndarray] = {}

    def _d(self, key, deriv, axis, prof, offset):
        """Stretched derivative: deriv/kappa plus the CPML memory term."""
        kap, b, c, start = prof
        n = deriv.shape[axis]
        idx = np.arange(n) + offset
        shp = [1, 1, 1]
        shp[axis] = n
        out = deriv / kap[idx].reshape(shp).astype(deriv.dtype)
        s = max(start - offset, 0)
        if s < n:
            sl = [slice(None)] * 3
            sl[axis] = slice(s, None)
            sl = tuple(sl)
            shp[axis] = n - s
            psi = self.psi.get(key)
            if psi is None:
                psi = self.psi[key] = np.zeros(deriv[sl].shape, deriv.dtype)
            psi *= b[idx[s:]].reshape(shp).astype(deriv.dtype)
            psi += c[idx[s:]].reshape(shp).astype(deriv.dtype) * deriv[sl]
            out[sl] += psi
        return out

    def step(self):
        dt, (px, py, pz), (hx, hy, hz) = self.dt, self.p_half, self.p_int
        Ex, Ey, Ez = self.Ex, self.Ey, self.Ez
        self.Hx -= dt * (self._d("hxy", Ez[:, 1:] - Ez[:, :-1], 1, py, 0)
                         - self._d("hxz", Ey[:, :, 1:] - Ey[:, :, :-1], 2, pz, 0))
        self.Hy -= dt * (self._d("hyz", Ex[:, :, 1:] - Ex[:, :, :-1], 2, pz, 0)
                         - self._d("hyx", Ez[1:] - Ez[:-1], 0, px, 0))
        self.Hz -= dt * (self._d("hzx", Ey[1:] - Ey[:-1], 0, px, 0)
                         - self._d("hzy", Ex[:, 1:] - Ex[:, :-1], 1, py, 0))
        Hx, Hy, Hz = self.Hx, self.Hy, self.Hz
        self.Ex[:, 1:-1, 1:-1] += self.cEx * (
            self._d("exy", Hz[:, 1:, 1:-1] - Hz[:, :-1, 1:-1], 1, hy, 1)
            - self._d("exz", Hy[:, 1:-1, 1:] - Hy[:, 1:-1, :-1], 2, hz, 1))
        # tangential H is odd through the magnetic wall at x = 0
        hzm = np.concatenate([-Hz[:1, :, 1:-1], Hz[:, :, 1:-1]])
        self.Ey[:-1, :, 1:-1] += self.cEy * (
            self._d("eyz", Hx[:-1, :, 1:] - Hx[:-1, :, :-1], 2, hz, 1)
            - self._d("eyx", hzm[1:] - hzm[:-1], 0, hx, 0))
        hym = np.concatenate([-Hy[:1, 1:-1, :], Hy[:, 1:-1, :]])
        self.Ez[:-1, 1:-1, :] += self.cEz * (
            self._d("ezx", hym[1:] - hym[:-1], 0, hx, 0)
            - self._d("ezy", Hx[:-1, 1:] - Hx[:-1, :-1], 1, hy, 1))


def matrix_pencil(y, dt: float, modes: int = 16):
    """Exponents s_k and amplitudes a_k with y(n dt) ~ sum a_k exp(s_k n dt)."""
    y = np.asarray(y, dtype=complex)
    n = len(y)
    L = n // 3
    Y = np.array([y[i:i + L + 1] for i in range(n - L)])
    _, S, Vh = np.linalg.svd(Y, full_matrices=False)
    M = min(modes, int(np.sum(S > S[0] * 1e-10)))
    V = Vh[:M].conj().T
    z = np.linalg.eigvals(np.linalg.pinv(V[:-1]) @ V[1:])
    a = np.linalg.lstsq(np.vander(z, n, increasing=True).T, y, rcond=None)[0]
    return np.log(z) / dt, a


def ringdown(eps_cell, body: tuple[int, int, int], f_guess: float, source_z: int,
             air: int = 14, npml: int = 12, periods: float = 45.0) -> tuple[float, float]:
    """(frequency, Q) of the dominant x-dipole mode of the body.

    body is the resonator's extent in cells along x, y, z from the symmetry
    planes; f_guess (cycles per unit time) centres the exciting pulse and only
    needs to be within some tens of percent.
    """
    shape = tuple(n + air + npml for n in body)
    g = _Yee(shape, eps_cell, npml)
    tau = 1.0 / f_guess
    t0 = 3.0 * tau
    steps = int((t0 + 3.0 * tau + periods / f_guess) / g.dt)
    start = int((t0 + 3.0 * tau) / g.dt)
    probe = np.empty(steps)
    for n in range(steps):
        g.step()
        t = (n + 1) * g.dt
        # a y-directed current where the mode's E field circulates about x
        g.Ey[1, 0, source_z] += math.exp(-((t - t0) / tau) ** 2) * math.sin(2 * math.pi * f_guess * (t - t0))
        probe[n] = g.Hx[0, 0, 0]
    dec = max(1, int(1.0 / (20.0 * f_guess * g.dt)))
    s, a = matrix_pencil(probe[start::dec], dec * g.dt)
    f = s.imag / (2 * math.pi)
    q = s.imag / (-2.0 * s.real)
    # The magnetic dipole is the lowest mode this symmetry admits, but not
    # always the loudest at the probe: a long bar rings a high-Q higher mode
    # just as hard. So take the lowest frequency that carries real amplitude.
    ok = (f > 0.3 * f_guess) & (q > 0.5) & (np.abs(a) > 0.1 * np.abs(a[f > 0]).max())
    k = int(np.flatnonzero(ok)[np.argmin(f[ok])])
    return f[k], q[k]


def brick(eps_r: float, half_x: int, half_y: int, height: int, k0h_guess: float,
          **kw) -> tuple[float, float]:
    """(k0*height, Q) of the brick |x| < half_x, |y| < half_y, 0 < z < height
    on a ground plane, for the mode whose magnetic dipole lies along x."""
    def eps_cell(x, y, z):
        return 1.0 + (eps_r - 1.0) * ((np.abs(x) < half_x) & (np.abs(y) < half_y) & (np.abs(z) < height))
    f, q = ringdown(eps_cell, (half_x, half_y, height), k0h_guess / (2 * math.pi * height),
                    max(1, int(0.7 * height)), **kw)
    return 2 * math.pi * f * height, q


def hemisphere(eps_r: float, radius: float, **kw) -> tuple[float, float]:
    """(k0*a, Q) of a staircased hemisphere of the given radius in cells; each
    cell's permittivity is its volume fraction of dielectric, from 4^3 samples."""
    sub = (np.arange(4) + 0.5) / 4 - 0.5

    def eps_cell(x, y, z):
        inside = sum(((x + p) ** 2 + (y + q) ** 2 + (z + r) ** 2 < radius ** 2).astype(float)
                     for p in sub for q in sub for r in sub)
        return 1.0 + (eps_r - 1.0) * inside / 64.0
    n = int(math.ceil(radius))
    f, q = ringdown(eps_cell, (n, n, n), 0.94 * math.sqrt(10.0 / eps_r) / (2 * math.pi * radius),
                    max(1, int(0.6 * radius)), **kw)
    return 2 * math.pi * f * radius, q
