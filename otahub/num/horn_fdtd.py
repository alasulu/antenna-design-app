"""A sectoral horn by FDTD, and a closed near-to-far-field transform to give its gain.

The horn's gain across its unflared, waveguide-sized dimension is where the
aperture models disagree (section 4 of the handover): the conventional
aperture-power form, a Huygens aperture and the same aperture in a ground plane
sit up to a decibel apart. Nothing but a full-wave solve of the horn itself
decides between them. This one models the horn as it is built: a TE10 guide from
a back short, the flare, thin walls, open space all round.

Grid: `dra._Yee`, quarter space - a magnetic wall at x = 0 and electric walls at
y = 0 and z = 0 - which is exactly the TE10 horn's symmetry, E along y even
about both mid-planes, the z = 0 wall being the guide's back short. The horn then
sits on an infinite plane at its throat, and that plane is not harmless: it
mirrors the radiation off the thin walls' outer edges back into the beam, and put
the H-plane horn 4 dB below every aperture model. `free=True` uses `_YeeFree`
instead - CPML at both ends of z, the back short a plate over the guide's own
cross-section, the box closed behind the horn - which is the horn as built, and
what the specs use. Units: cell 1, c 1, H scaled by eta0. Walls are staircased
onto cell faces; tangential E on them is zeroed after every update.

`BoxTransform` records the tangential fields on the faces x = X, y = Y and z = Z
of a box round the horn, DFTs them at f, and radiates the equivalent currents
J = n x H and M = -n x E together with their images through the three walls
(PMC: J tangential kept, normal reversed, M the other way; PEC: J tangential
reversed, M tangential kept). Directivity is taken over the physical half space
z > 0, or the whole sphere on the free grid. It is checked against a short dipole
over the back plane, whose directivity follows from the array factor exactly, and
a lone short dipole on the free grid.
"""
from __future__ import annotations

import math

import numpy as np

from .dra import _Yee

__all__ = ["BoxTransform", "dipole_over_plane", "dipole_over_plane_exact", "sectoral_horn"]

# image signs: for a reflection through the plane normal to axis a, the (x, y, z) sign of J and M
_PMC = {"J": lambda a: [(-1 if c == a else 1) for c in range(3)], "M": lambda a: [(1 if c == a else -1) for c in range(3)]}
_PEC = {"J": lambda a: [(1 if c == a else -1) for c in range(3)], "M": lambda a: [(-1 if c == a else 1) for c in range(3)]}
_WALLS = ((0, _PMC), (1, _PEC), (2, _PEC))


class _YeeFree(_Yee):
    """Quarter-space in x and y (PMC at x = 0, PEC at y = 0), CPML at both ends of z: no wall behind."""

    def __init__(self, shape, eps_cell, npml: int, dtype=np.float32):
        super().__init__(shape, eps_cell, npml, dtype)
        from .patch_fdtd import _cpml_both
        nz = shape[2]
        self.p_int[2] = _cpml_both(nz + 1, npml, False, self.dt)
        self.p_half[2] = _cpml_both(nz, npml, True, self.dt)


class BoxTransform:
    """DFT at f of the tangential fields on the faces x = i0, y = j0, z = k0 of the box
    [0, i0] x [0, j0] x [0, k0] - the other faces lie on the walls. With j_lo given the
    grid has no wall at y = 0 (the half-space grid): the box is [0, i0] x [j_lo, j0] x
    [0, k0], with a face at y = j_lo too and no image through y. With k_lo given, likewise
    in z (the free grid): a face at z = k_lo, no image through z, the whole sphere.

    With k_plane given, z = k_plane is an infinite PEC plane with an opening - a flange
    round an aperture, the source behind it. Then the only equivalent source is the
    plane itself: M = -z x E over [0, i0] x [0, j0] of it, outward normal +z into the
    half space in front, imaged through the plane (which doubles M and cancels J), and
    E is zero on the metal, so only the opening counts and the box's size does not. A
    box round the space in front instead encloses no source, and one cutting through
    the flange read 0.36 dB differently as its size changed."""

    def __init__(self, i0: int, j0: int, k0: int, f: float, dt: float, j_lo: int | None = None,
                 k_lo: int | None = None, k_plane: int | None = None):
        if k_lo is not None and k_plane is not None:
            raise ValueError("k_lo (no plane) and k_plane (a plane) exclude each other")
        self.i0, self.j0, self.k0, self.f, self.dt, self.j_lo, self.k_lo = i0, j0, k0, f, dt, j_lo, k_lo
        self.k_plane = k_plane
        self.acc: dict[str, np.ndarray] = {}

    @property
    def _k_first(self) -> int:
        """The box's first z cell: the free grid's k_lo, a flange plane, or the wall at 0."""
        return self.k_lo if self.k_lo is not None else (self.k_plane or 0)

    def _add(self, key, value, w):
        if key in self.acc:
            self.acc[key] += value * w
        else:
            self.acc[key] = value.astype(complex) * w

    def record(self, g, n: int) -> None:
        """Call after the n-th g.step(): E is at (n + 1) dt, H at (n + 1/2) dt."""
        i0, j0, k0 = self.i0, self.j0, self.k0
        a = self.j_lo or 0                               # first y cell of the box
        c = self._k_first                                # first z cell of the box
        we = np.exp(-2j * math.pi * self.f * (n + 1) * self.dt) * self.dt
        wh = np.exp(-2j * math.pi * self.f * (n + 0.5) * self.dt) * self.dt
        Ex, Ey, Ez, Hx, Hy, Hz = g.Ex, g.Ey, g.Ez, g.Hx, g.Hy, g.Hz
        # x = i0 face, points (i0, j + 1/2, k + 1/2)
        self._add("x_Ey", 0.5 * (Ey[i0, a:j0, c:k0] + Ey[i0, a:j0, c + 1:k0 + 1]), we)
        self._add("x_Ez", 0.5 * (Ez[i0, a:j0, c:k0] + Ez[i0, a + 1:j0 + 1, c:k0]), we)
        self._add("x_Hy", 0.25 * (Hy[i0 - 1, a:j0, c:k0] + Hy[i0, a:j0, c:k0] + Hy[i0 - 1, a + 1:j0 + 1, c:k0] + Hy[i0, a + 1:j0 + 1, c:k0]), wh)
        self._add("x_Hz", 0.25 * (Hz[i0 - 1, a:j0, c:k0] + Hz[i0, a:j0, c:k0] + Hz[i0 - 1, a:j0, c + 1:k0 + 1] + Hz[i0, a:j0, c + 1:k0 + 1]), wh)
        # y faces, points (i + 1/2, y, k + 1/2)
        for key, jy in (("y", j0), ("ylo", self.j_lo)):
            if jy is None:
                continue
            self._add(key + "_Ex", 0.5 * (Ex[:i0, jy, c:k0] + Ex[:i0, jy, c + 1:k0 + 1]), we)
            self._add(key + "_Ez", 0.5 * (Ez[:i0, jy, c:k0] + Ez[1:i0 + 1, jy, c:k0]), we)
            self._add(key + "_Hx", 0.25 * (Hx[:i0, jy - 1, c:k0] + Hx[:i0, jy, c:k0] + Hx[1:i0 + 1, jy - 1, c:k0] + Hx[1:i0 + 1, jy, c:k0]), wh)
            self._add(key + "_Hz", 0.25 * (Hz[:i0, jy - 1, c:k0] + Hz[:i0, jy, c:k0] + Hz[:i0, jy - 1, c + 1:k0 + 1] + Hz[:i0, jy, c + 1:k0 + 1]), wh)
        # z faces, points (i + 1/2, j + 1/2, z)
        for key, kz in (("z", k0), ("zlo", self.k_lo if self.k_lo is not None else self.k_plane)):
            if kz is None:
                continue
            self._add(key + "_Ex", 0.5 * (Ex[:i0, a:j0, kz] + Ex[:i0, a + 1:j0 + 1, kz]), we)
            self._add(key + "_Ey", 0.5 * (Ey[:i0, a:j0, kz] + Ey[1:i0 + 1, a:j0, kz]), we)
            self._add(key + "_Hx", 0.25 * (Hx[:i0, a:j0, kz - 1] + Hx[:i0, a:j0, kz] + Hx[1:i0 + 1, a:j0, kz - 1] + Hx[1:i0 + 1, a:j0, kz]), wh)
            self._add(key + "_Hy", 0.25 * (Hy[:i0, a:j0, kz - 1] + Hy[:i0, a:j0, kz] + Hy[:i0, a + 1:j0 + 1, kz - 1] + Hy[:i0, a + 1:j0 + 1, kz]), wh)

    def _currents(self):
        """[(points (P, 3), J (P, 3), M (P, 3))] on the three faces, unit cell areas."""
        a, i0, j0, k0 = self.acc, self.i0, self.j0, self.k0
        lo = self.j_lo or 0
        klo = self._k_first
        if self.k_plane is not None:                     # the aperture plane alone, normal +z
            ii, jj = np.meshgrid(np.arange(i0) + 0.5, np.arange(lo, j0) + 0.5, indexing="ij")
            P = np.stack([ii.ravel(), jj.ravel(), np.full(ii.size, float(klo))], 1)
            z = np.zeros(ii.size, complex)
            return [(P, np.stack([-a["zlo_Hy"].ravel(), a["zlo_Hx"].ravel(), z], 1),
                     np.stack([a["zlo_Ey"].ravel(), -a["zlo_Ex"].ravel(), z], 1))]
        out = []
        jj, kk = np.meshgrid(np.arange(lo, j0) + 0.5, np.arange(klo, k0) + 0.5, indexing="ij")
        P = np.stack([np.full(jj.size, float(i0)), jj.ravel(), kk.ravel()], 1)
        z = np.zeros(jj.size, complex)
        out.append((P, np.stack([z, -a["x_Hz"].ravel(), a["x_Hy"].ravel()], 1), np.stack([z, a["x_Ez"].ravel(), -a["x_Ey"].ravel()], 1)))
        ii, kk = np.meshgrid(np.arange(i0) + 0.5, np.arange(klo, k0) + 0.5, indexing="ij")
        P = np.stack([ii.ravel(), np.full(ii.size, float(j0)), kk.ravel()], 1)
        z = np.zeros(ii.size, complex)
        out.append((P, np.stack([a["y_Hz"].ravel(), z, -a["y_Hx"].ravel()], 1), np.stack([-a["y_Ez"].ravel(), z, a["y_Ex"].ravel()], 1)))
        if self.j_lo is not None:                        # the y = j_lo face, outward normal -y
            P = np.stack([ii.ravel(), np.full(ii.size, float(self.j_lo)), kk.ravel()], 1)
            out.append((P, np.stack([-a["ylo_Hz"].ravel(), z, a["ylo_Hx"].ravel()], 1),
                        np.stack([a["ylo_Ez"].ravel(), z, -a["ylo_Ex"].ravel()], 1)))
        ii, jj = np.meshgrid(np.arange(i0) + 0.5, np.arange(lo, j0) + 0.5, indexing="ij")
        P = np.stack([ii.ravel(), jj.ravel(), np.full(ii.size, float(k0))], 1)
        z = np.zeros(ii.size, complex)
        out.append((P, np.stack([-a["z_Hy"].ravel(), a["z_Hx"].ravel(), z], 1), np.stack([a["z_Ey"].ravel(), -a["z_Ex"].ravel(), z], 1)))
        if self.k_lo is not None or self.k_plane is not None:   # the z = k_lo face, outward normal -z
            P = np.stack([ii.ravel(), jj.ravel(), np.full(ii.size, float(klo))], 1)
            out.append((P, np.stack([a["zlo_Hy"].ravel(), -a["zlo_Hx"].ravel(), z], 1),
                        np.stack([-a["zlo_Ey"].ravel(), a["zlo_Ex"].ravel(), z], 1)))
        return out

    def intensity(self, theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
        """|E_theta|^2 + |E_phi|^2 up to a constant, from the faces and their seven images."""
        th, ph = np.broadcast_arrays(np.asarray(theta, float), np.asarray(phi, float))
        st, ct, sp, cp = np.sin(th).ravel(), np.cos(th).ravel(), np.sin(ph).ravel(), np.cos(ph).ravel()
        rhat = np.stack([st * cp, st * sp, ct], 1)
        k = 2 * math.pi * self.f
        N = np.zeros((len(rhat), 3), complex)
        L = np.zeros((len(rhat), 3), complex)
        faces = self._currents()
        if self.k_plane is not None:                     # image through the flange's plane, not z = 0
            shift = np.array([0.0, 0.0, float(self.k_plane)])
            faces = [(P - shift, J, M) for P, J, M in faces]
        npts = sum(len(P) for P, _, _ in faces)
        step = max(1, int(4e6 / max(npts, 1)))                   # directions per chunk: ~64 MB of phases
        for P, J, M in faces:
            for mask in range(8):
                if self.j_lo is not None and mask & 2:          # no wall at y = 0 on the half-space grid
                    continue
                if self.k_lo is not None and mask & 4:          # no wall at z = 0 on the free grid
                    continue
                sP = np.ones(3)
                sJ = np.ones(3)
                sM = np.ones(3)
                for bit, (axis, rule) in enumerate(_WALLS):
                    if mask >> bit & 1:
                        sP[axis] = -1
                        sJ *= rule["J"](axis)
                        sM *= rule["M"](axis)
                Pm, Jm, Mm = P * sP, J * sJ, M * sM
                for d0 in range(0, len(rhat), step):
                    phase = np.exp(1j * k * (rhat[d0:d0 + step] @ Pm.T))     # (dirs, points)
                    N[d0:d0 + step] += phase @ Jm
                    L[d0:d0 + step] += phase @ Mm
        thv = np.stack([ct * cp, ct * sp, -st], 1)
        phv = np.stack([-sp, cp, np.zeros_like(cp)], 1)
        Nt, Np = np.sum(N * thv, 1), np.sum(N * phv, 1)
        Lt, Lp = np.sum(L * thv, 1), np.sum(L * phv, 1)
        return (np.abs(Lp + Nt) ** 2 + np.abs(Lt - Np) ** 2).reshape(th.shape)

    def peak(self, n_theta: int = 60, n_phi: int = 72) -> tuple[float, float, float]:
        """(peak directivity, its theta, its phi) over the physical half space z > 0 (the
        whole sphere on the free grid)."""
        x, w = np.polynomial.legendre.leggauss(n_theta)
        span = math.pi if self.k_lo is not None else 0.5 * math.pi          # the whole sphere, or the half above the wall
        th = 0.5 * span * (x + 1)
        wt = 0.5 * span * w * np.sin(th)
        ph = 2 * math.pi * (np.arange(n_phi) + 0.5) / n_phi
        TH, PH = np.meshgrid(th, ph, indexing="ij")
        U = self.intensity(TH, PH)
        P = float(np.sum(U * wt[:, None]) * 2 * math.pi / n_phi)
        k = np.unravel_index(np.argmax(U), U.shape)
        return 4 * math.pi * float(U[k]) / P, float(TH[k]), float(PH[k])

    def directivity(self, theta: float = 0.0, phi: float = 0.0, n_theta: int = 60, n_phi: int = 72) -> float:
        """Toward (theta, phi), over the physical half space z > 0 (the whole sphere on the free grid)."""
        x, w = np.polynomial.legendre.leggauss(n_theta)
        span = math.pi if self.k_lo is not None else 0.5 * math.pi          # the whole sphere, or the half above the wall
        th = 0.5 * span * (x + 1)
        wt = 0.5 * span * w * np.sin(th)
        ph = 2 * math.pi * (np.arange(n_phi) + 0.5) / n_phi
        TH, PH = np.meshgrid(th, ph, indexing="ij")
        U = self.intensity(TH, PH)
        P = float(np.sum(U * wt[:, None]) * 2 * math.pi / n_phi)
        return 4 * math.pi * float(self.intensity(np.array([theta]), np.array([phi]))[0]) / P


def _pulse(t, f, tau, t0):
    return math.exp(-((t - t0) / tau) ** 2) * math.sin(2 * math.pi * f * (t - t0))


def dipole_over_plane(height: float, cells_per_lambda: int = 20, margin: float = 0.4, npml: int = 12,
                      free: bool = False) -> float:
    """Broadside directivity of a short y-directed source `height` wavelengths above the
    z = 0 plane, at the corner of the two symmetry walls, by the box transform - to be
    compared with the exact array-factor value. With `free` there is no plane (the free
    grid, the box closed below the source too): a lone short dipole, 1.5 exactly."""
    N = cells_per_lambda
    f = 1.0 / N
    m = int(round(margin * N))
    kh = (npml + 4 + m) if free else int(round(height * N))
    i0, j0, k0 = m, m, kh + m
    shape = (i0 + 4 + npml, j0 + 4 + npml, k0 + 4 + npml)
    g = (_YeeFree if free else _Yee)(shape, lambda x, y, z: 1.0 + 0 * (x + y + z), npml)
    box = BoxTransform(i0, j0, k0, f, g.dt, k_lo=(kh - m) if free else None)
    tau, t0 = 1.0 / f, 3.0 / f
    steps = int((t0 + 3 * tau + 12 / f) / g.dt)
    for n in range(steps):
        g.step()
        g.Ey[0, 0, kh] += _pulse((n + 1) * g.dt, f, tau, t0)
        box.record(g, n)
    return box.directivity()


def dipole_over_plane_exact(height: float, n: int = 400) -> float:
    """y-directed element pattern (1 - sin^2 th sin^2 ph), array factor 2 sin(kh cos th)."""
    x, w = np.polynomial.legendre.leggauss(n)
    th = 0.25 * math.pi * (x + 1)
    wt = 0.25 * math.pi * w * np.sin(th)
    ph = 2 * math.pi * (np.arange(2 * n) + 0.5) / (2 * n)
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    kh = 2 * math.pi * height
    U = (1 - np.sin(TH) ** 2 * np.sin(PH) ** 2) * np.sin(kh * np.cos(TH)) ** 2
    P = float(np.sum(U * wt[:, None]) * 2 * math.pi / (2 * n))
    return 4 * math.pi * math.sin(kh) ** 2 / P


def _wall_masks(shape, XA, YB, k_ap, kz0: int = 0):
    """Zero-tangential-E masks (Ex, Ey, Ez) for a rectangular guide and flare of half-width
    XA[k] (x node of the side wall) and half-height YB[k] (y node of the top wall) in the
    layer z in [kz0 + k, kz0 + k + 1], for k < k_ap; walls are zero-thickness sheets, open
    beyond. With kz0 > 0 the guide is closed at z = kz0 by a plate over its own cross-section."""
    nx, ny, nz = shape
    mx = np.zeros((nx, ny + 1, nz + 1), bool)
    my = np.zeros((nx + 1, ny, nz + 1), bool)
    mz = np.zeros((nx + 1, ny + 1, nz), bool)
    if kz0:
        mx[:XA[0], :YB[0] + 1, kz0] = True
        my[:XA[0] + 1, :YB[0], kz0] = True
    for kk in range(k_ap):
        k = kz0 + kk
        xa, yb = XA[kk], YB[kk]
        # side wall x = xa: Ey on its two bounding z-lines, Ez across it, for y up to yb
        my[xa, :yb, k] = True
        my[xa, :yb, k + 1] = True
        mz[xa, :yb + 1, k] = True
        # top wall y = yb: Ex on its bounding z-lines, Ez across it, for x up to xa
        mx[:xa, yb, k] = True
        mx[:xa, yb, k + 1] = True
        mz[:xa + 1, yb, k] = True
        if kk > 0 and (XA[kk - 1] != xa or YB[kk - 1] != yb):
            # the step face at z = k between the two layers' outlines
            x_lo, x_hi = sorted((XA[kk - 1], xa))
            y_lo, y_hi = sorted((YB[kk - 1], yb))
            if x_hi > x_lo:
                mx[x_lo:x_hi, :min(YB[kk - 1], yb) + 1, k] = True
                my[x_lo:x_hi + 1, :min(YB[kk - 1], yb), k] = True
            if y_hi > y_lo:
                mx[:min(XA[kk - 1], xa), y_lo:y_hi + 1, k] = True
                my[:min(XA[kk - 1], xa) + 1, y_lo:y_hi, k] = True
    return mx, my, mz


def sectoral_horn(a_wg: float, b_wg: float, A: float, B: float, flare_len: float, cells_per_lambda: int = 30,
                  guide_len: float = 1.0, margin: float = 0.35, npml: int = 12, periods: float = 40.0,
                  flange: bool = False, free: bool = False) -> dict:
    """Boresight directivity of a horn (lengths in wavelengths at f): a guide a_wg x b_wg
    from a back short at z = 0, `guide_len` long, flaring linearly over `flare_len` to an
    aperture A x B (A = a_wg for an E-plane horn, B = b_wg for an H-plane one). TE10 is
    launched a quarter guide wavelength in front of the short."""
    N = cells_per_lambda
    f = 1.0 / N
    kg = int(round(guide_len * N))
    kf = int(round(flare_len * N))
    k_ap = kg + kf
    xa0, yb0 = int(round(a_wg / 2 * N)), int(round(b_wg / 2 * N))
    xa1, yb1 = int(round(A / 2 * N)), int(round(B / 2 * N))
    XA = [xa0] * kg + [int(round(xa0 + (xa1 - xa0) * (k - kg + 0.5) / kf)) for k in range(kg, k_ap)]
    YB = [yb0] * kg + [int(round(yb0 + (yb1 - yb0) * (k - kg + 0.5) / kf)) for k in range(kg, k_ap)]
    m = int(round(margin * N))
    kz0 = (npml + 4 + m) if free else 0          # free: room behind the back short for the box and the CPML
    i0, j0, k0 = xa1 + m, yb1 + m, kz0 + k_ap + m
    shape = (i0 + 4 + npml, j0 + 4 + npml, k0 + 4 + npml)
    g = (_YeeFree if free else _Yee)(shape, lambda x, y, z: 1.0 + 0 * (x + y + z), npml)
    mx, my, mz = _wall_masks(shape, XA, YB, k_ap, kz0)
    if flange:                                   # a conducting plane round the aperture, out through the CPML
        kf_ = kz0 + k_ap
        mx[:, :, kf_] |= (np.arange(shape[0])[:, None] + 0.5 > xa1) | (np.arange(shape[1] + 1)[None, :] >= yb1)
        my[:, :, kf_] |= (np.arange(shape[0] + 1)[:, None] >= xa1) | (np.arange(shape[1])[None, :] + 0.5 > yb1)
        # the box closes on the flange and images through it: the plane is infinite
        box = BoxTransform(i0, j0, k0, f, g.dt, k_plane=kf_)
    else:
        box = BoxTransform(i0, j0, k0, f, g.dt, k_lo=(kz0 - m) if free else None)
    lam_g = 1.0 / math.sqrt(1.0 - (1.0 / (2.0 * a_wg)) ** 2)            # guide wavelength, in wavelengths
    ks = max(1, int(round(0.25 * lam_g * N)))
    xs = np.arange(xa0 + 1)
    prof = np.cos(math.pi * xs / (2 * xa0))[:, None] * np.ones((1, yb0))   # TE10 across x, uniform in y
    tau, t0 = 1.5 / f, 4.5 / f
    steps = int((t0 + 3 * tau + periods / f) / g.dt)
    for n in range(steps):
        g.step()
        g.Ex[mx] = 0.0
        g.Ey[my] = 0.0
        g.Ez[mz] = 0.0
        g.Ey[:xa0 + 1, :yb0, kz0 + ks] += prof * _pulse((n + 1) * g.dt, f, tau, t0)
        box.record(g, n)
    return dict(directivity=box.directivity(), cells_per_lambda=N, box=(i0, j0, k0), steps=steps)
