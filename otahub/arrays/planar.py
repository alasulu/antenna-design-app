"""Planar arrays: rectangular and triangular lattices, steering, directivity.

Convention differs deliberately from the linear module. A linear array lies
along z and theta is measured from that axis, so broadside is 90 degrees. A
planar array lies in the x-y plane and theta is measured from the array NORMAL,
so broadside is 0 degrees and the array plane is 90. Both are the usual choice
for their geometry; mixing them up is the obvious trap, so every function here
says which it means.

Positions and spacings are in wavelengths throughout.
"""
from __future__ import annotations

import math

import numpy as np

from ..core.pattern import Pattern
from .elements import ElementPattern, power_kernel

__all__ = [
    "rectangular_lattice", "triangular_lattice", "separable_weights",
    "steering_phase", "planar_array_factor", "planar_directivity",
    "planar_pattern_cut", "planar_beam_cut", "grating_lobe_free_spacing_planar",
    "lattice_element_saving", "planar_summarise",
]


# ------------------------------------------------------------------ lattices

def rectangular_lattice(nx: int, ny: int, dx: float,
                        dy: float | None = None) -> np.ndarray:
    """`nx` by `ny` elements on a rectangular grid, centred on the origin.

    Returns an (nx*ny, 2) array of (x, y) in wavelengths, ordered x-fastest so
    it matches :func:`separable_weights`.
    """
    if nx < 1 or ny < 1:
        raise ValueError("nx and ny must be at least 1")
    if dx <= 0 or (dy is not None and dy <= 0):
        raise ValueError("spacings must be positive")
    dy = dx if dy is None else dy
    xs = (np.arange(nx) - (nx - 1) / 2.0) * dx
    ys = (np.arange(ny) - (ny - 1) / 2.0) * dy
    xx, yy = np.meshgrid(xs, ys, indexing="ij")
    return np.column_stack([xx.ravel(), yy.ravel()])


def triangular_lattice(nx: int, ny: int, s: float) -> np.ndarray:
    """Equilateral triangular lattice of nearest-neighbour spacing `s`.

    Alternate rows are offset by s/2 and the row pitch is s*sqrt(3)/2, so every
    element has six neighbours at exactly `s`. For the same grating-lobe-free
    scan volume this lattice needs 13.4% fewer elements than a rectangular one
    — see :func:`lattice_element_saving`.
    """
    if nx < 1 or ny < 1:
        raise ValueError("nx and ny must be at least 1")
    if s <= 0:
        raise ValueError("spacing must be positive")
    pitch = s * math.sqrt(3) / 2.0
    pts = []
    for j in range(ny):
        offset = s / 2.0 if j % 2 else 0.0
        for i in range(nx):
            pts.append(((i - (nx - 1) / 2.0) * s + offset,
                        (j - (ny - 1) / 2.0) * pitch))
    return np.asarray(pts, dtype=float)


def separable_weights(wx, wy) -> np.ndarray:
    """Outer product of two linear tapers, flattened to match the lattice order.

    Separable excitation is what makes a rectangular array tractable: the
    pattern is the product of the two linear array factors, so a Chebyshev
    taper in each axis gives the design sidelobe level along both principal
    planes (though not along the diagonals, where the two products multiply).
    """
    wx = np.asarray(wx, dtype=float)
    wy = np.asarray(wy, dtype=float)
    return np.outer(wx, wy).ravel()


# ------------------------------------------------------------------ patterns

def steering_phase(positions, scan_theta_deg: float = 0.0,
                   scan_phi_deg: float = 0.0) -> np.ndarray:
    """Element phases that point the beam at (scan_theta, scan_phi).

    theta is measured from the array normal, so 0 is broadside.
    """
    pos = np.asarray(positions, dtype=float)
    st = math.sin(math.radians(scan_theta_deg))
    u0 = st * math.cos(math.radians(scan_phi_deg))
    v0 = st * math.sin(math.radians(scan_phi_deg))
    return np.exp(-2j * math.pi * (pos[:, 0] * u0 + pos[:, 1] * v0))


def planar_array_factor(positions, weights, theta, phi,
                        scan_theta_deg: float = 0.0,
                        scan_phi_deg: float = 0.0) -> np.ndarray:
    """Complex array factor at the given angles, theta from the array normal."""
    pos = np.asarray(positions, dtype=float)
    w = np.asarray(weights, dtype=complex) * steering_phase(
        pos, scan_theta_deg, scan_phi_deg)
    theta = np.asarray(theta, dtype=float)
    phi = np.asarray(phi, dtype=float)
    u = np.sin(theta) * np.cos(phi)
    v = np.sin(theta) * np.sin(phi)
    shape = np.broadcast(u, v).shape
    out = np.zeros(shape, dtype=complex)
    for (x, y), wi in zip(pos, w):
        out = out + wi * np.exp(2j * math.pi * (x * u + y * v))
    return out


def planar_directivity(positions, weights, scan_theta_deg: float = 0.0,
                       scan_phi_deg: float = 0.0,
                       half_space: bool = False,
                       element: ElementPattern | None = None) -> float:
    """Directivity toward the scan direction on an arbitrary planar layout.

    Exact, not integrated numerically. The average of exp(j*k.d) over the
    sphere is sin(kd)/(kd), so the radiated power is a double sum over
    sinc(2*|r_m - r_n|/lambda) — the same kernel the linear module uses,
    with the separation taken in the plane rather than along a line.

    With `half_space` the result is doubled, which is the ideal ground-plane
    case: the array can only radiate upward, so the same beam carries all the
    power. Without it, isotropic elements radiate equally up and down and half
    the power goes into the mirror beam.

    With an `element` pattern (see `elements`) the kernel is that pattern's
    power integral instead of sinc, still exact, and the element decides which
    hemispheres radiate - so `half_space` does not apply. The figure is toward
    the scan direction; a sloping element pattern can pull the peak slightly
    off it (toward broadside), which this deliberately does not chase.
    """
    pos = np.asarray(positions, dtype=float)
    w0 = np.asarray(weights, dtype=complex)
    w = w0 * steering_phase(pos, scan_theta_deg, scan_phi_deg)
    if element is not None:
        if half_space:
            raise ValueError("an element pattern sets its own hemisphere; half_space does not apply")
        diff = np.round((pos[:, None, :] - pos[None, :, :]).reshape(-1, 2), 9)
        uniq, inv = np.unique(diff, axis=0, return_inverse=True)
        Kd = power_kernel(element, uniq)[inv.ravel()].reshape(len(pos), len(pos))
        den = float(np.real((w[:, None] * np.conj(w)[None, :] * Kd).sum()))
        if den <= 0:
            raise ValueError("degenerate array: no radiated power")
        pk = float(element(math.radians(scan_theta_deg), math.radians(scan_phi_deg)))
        return 4.0 * math.pi * pk * abs(w0.sum()) ** 2 / den
    sep = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1)
    # The peak is at the SCAN direction, where the steering phase cancels and
    # the elements add as their bare weights. Summing the steered weights
    # instead would evaluate the pattern at broadside, which after steering is
    # near a null - and would report a scanned array as having lost 25 dB.
    num = abs(w0.sum()) ** 2
    den = float(np.real((w[:, None] * np.conj(w)[None, :] * np.sinc(2.0 * sep)).sum()))
    if den <= 0:
        raise ValueError("degenerate array: no radiated power")
    return float(num / den) * (2.0 if half_space else 1.0)


def planar_pattern_cut(positions, weights, phi_deg: float = 0.0,
                       scan_theta_deg: float = 0.0, scan_phi_deg: float = 0.0,
                       n_theta: int = 2001, element: ElementPattern | None = None) -> Pattern:
    """One azimuth cut, as a :class:`Pattern` so the beamwidth and sidelobe
    machinery applies.

    The returned Pattern's theta axis is the angle in the cut plane measured
    from one horizon to the other, so BROADSIDE SITS AT 90 DEGREES — matching
    the linear module's convention, which is what the beamwidth and sidelobe
    finders expect. Only those readings are meaningful on a single cut; use
    :func:`planar_directivity` for directivity, which needs the whole sphere.
    """
    # psi runs -90 to +90 from the array normal, negative psi meaning the
    # opposite azimuth: sin(-psi) flips the direction through the origin in the
    # (u, v) plane, which is exactly phi + 180. Mapping psi onto the Pattern's
    # theta axis as psi + 90 puts BROADSIDE IN THE MIDDLE of the grid, where
    # the beamwidth and sidelobe finders can see both flanks of the main lobe.
    # Running theta from 0 instead leaves a broadside beam jammed against the
    # edge, and they return NaN and 0 dB respectively.
    psi = np.linspace(-math.pi / 2, math.pi / 2, n_theta)
    phi_r = math.radians(phi_deg)
    pos = np.asarray(positions, dtype=float)
    w = np.asarray(weights, dtype=complex) * steering_phase(
        pos, scan_theta_deg, scan_phi_deg)
    u = np.sin(psi) * math.cos(phi_r)
    v = np.sin(psi) * math.sin(phi_r)
    af = np.zeros(psi.shape, dtype=complex)
    for (x, y), wi in zip(pos, w):
        af = af + wi * np.exp(2j * math.pi * (x * u + y * v))
    power = np.abs(af) ** 2
    if element is not None:
        power = power * element(np.abs(psi), np.where(psi < 0, phi_r + math.pi, phi_r))
    theta_axis = psi + math.pi / 2
    return Pattern(theta_axis, np.array([phi_r]), power[:, None])


def _beam_basis(scan_theta_deg: float, scan_phi_deg: float):
    """Orthonormal frame around the beam: (beam, in-scan-plane, cross-plane)."""
    t = math.radians(scan_theta_deg)
    f = math.radians(scan_phi_deg)
    beam = np.array([math.sin(t) * math.cos(f), math.sin(t) * math.sin(f), math.cos(t)])
    z = np.array([0.0, 0.0, 1.0])
    cross = np.cross(z, beam)
    norm = np.linalg.norm(cross)
    if norm < 1e-12:
        # Broadside: the beam IS the normal, so the scan plane is undefined by
        # the geometry alone. Fall back on the requested azimuth, which makes
        # this case agree with the fixed-azimuth principal-plane cut.
        cross = np.array([-math.sin(f), math.cos(f), 0.0])
        in_plane = np.array([math.cos(f), math.sin(f), 0.0])
    else:
        cross = cross / norm
        in_plane = np.cross(beam, cross)
        in_plane = in_plane / np.linalg.norm(in_plane)
    return beam, in_plane, cross


def planar_beam_cut(positions, weights, plane: str = "scan",
                    scan_theta_deg: float = 0.0, scan_phi_deg: float = 0.0,
                    n_theta: int = 2001, element: ElementPattern | None = None) -> Pattern:
    """Great-circle cut THROUGH the beam, in the scan plane or across it.

    A fixed-azimuth cut only contains the beam when the array is unscanned or
    scanned along that very azimuth. Steer to (60 deg, phi = 0) and the phi =
    90 cut misses the beam entirely, so its "first sidelobe" is measured
    against whatever ridge happens to be highest there - a reading that looks
    like a catastrophic 0 dB sidelobe and means nothing.

    These cuts are great circles containing the beam: `plane="scan"` is the one
    that also contains the array normal, `plane="cross"` the one perpendicular
    to it. The beam sits at theta = 90 degrees on the returned Pattern, so the
    beamwidth and sidelobe finders see both of its flanks.
    """
    kind = plane.strip().lower()
    if kind not in ("scan", "cross"):
        raise ValueError(f"plane must be 'scan' or 'cross', not {plane!r}")

    if kind == "scan":
        # Sweep the angle from the NORMAL, signed, so the cut stays in the
        # forward hemisphere. Sweeping a great circle outward from the beam
        # instead runs past the array plane and into the mirror beam, which for
        # isotropic elements is at full amplitude - and the sidelobe finder
        # then dutifully reports a 0 dB sidelobe at every spacing. A real
        # planar array has a ground plane or directive elements and radiates
        # into this hemisphere only.
        f = math.radians(scan_phi_deg)
        psi = np.linspace(-math.pi / 2, math.pi / 2, n_theta)
        dirs = np.column_stack([
            np.sin(psi) * math.cos(f),
            np.sin(psi) * math.sin(f),
            np.cos(psi),
        ])
        # theta axis for the Pattern: psi + 90, so the beam lands at
        # 90 + scan_theta and both of its flanks are on the grid.
        theta_axis = psi + math.pi / 2
    else:
        # The cross-plane great circle is perpendicular to the scan plane, so
        # its z component is cos(psi)*cos(scan_theta): never negative for
        # |psi| <= 90, and the mirror beam is never reached.
        beam, _in_plane, cross = _beam_basis(scan_theta_deg, scan_phi_deg)
        psi = np.linspace(-math.pi / 2, math.pi / 2, n_theta)
        dirs = (np.cos(psi)[:, None] * beam[None, :]
                + np.sin(psi)[:, None] * cross[None, :])
        theta_axis = psi + math.pi / 2

    pos = np.asarray(positions, dtype=float)
    w = np.asarray(weights, dtype=complex) * steering_phase(
        pos, scan_theta_deg, scan_phi_deg)
    af = np.zeros(psi.shape, dtype=complex)
    for (x, y), wi in zip(pos, w):
        af = af + wi * np.exp(2j * math.pi * (x * dirs[:, 0] + y * dirs[:, 1]))
    power = np.abs(af) ** 2
    if element is not None:
        power = power * element(np.arccos(np.clip(dirs[:, 2], -1.0, 1.0)), np.arctan2(dirs[:, 1], dirs[:, 0]))
    return Pattern(theta_axis, np.array([0.0]), power[:, None])

# ------------------------------------------------------------------ lattices, limits

def grating_lobe_free_spacing_planar(scan_max_deg: float = 0.0,
                                     lattice: str = "rectangular") -> float:
    """Largest spacing with no grating lobe for scanning anywhere within a cone.

    The grating lobes of a lattice sit at its reciprocal-lattice points, so the
    condition is that the shortest reciprocal vector exceeds (1 + sin(scan))*k0.
    For a rectangular lattice that shortest vector is 2*pi/d, giving the
    familiar d <= lambda/(1 + sin(scan)). For an equilateral triangular lattice
    it is 4*pi/(sqrt(3)*s), so s <= 2*lambda/(sqrt(3)*(1 + sin(scan))) — a
    15.47% longer spacing, and the reason triangular lattices are standard in
    wide-scan arrays.

    Returns dx for a rectangular lattice, or the nearest-neighbour spacing s
    for a triangular one.
    """
    limit = 1.0 + abs(math.sin(math.radians(scan_max_deg)))
    kind = lattice.strip().lower()
    # Exact names, not prefixes: startswith("hex") happily accepted
    # "hexagonalish" and returned a confident number for it.
    if kind in ("rect", "rectangular", "square"):
        return 1.0 / limit
    if kind in ("tri", "triangular", "hex", "hexagonal"):
        return 2.0 / (math.sqrt(3.0) * limit)
    raise ValueError(
        f"unknown lattice {lattice!r}; use 'rectangular' or 'triangular'")


def lattice_element_saving() -> float:
    """Fraction of elements a triangular lattice saves over a rectangular one.

    Exactly 1 - sqrt(3)/2 = 0.1340. The triangular unit cell is 2/sqrt(3) larger
    in area for the same grating-lobe-free scan volume, so the same aperture
    needs 13.40% fewer elements — which at millimetre wave is 13.40% fewer
    transmit-receive modules, and the dominant cost of the array.
    """
    return 1.0 - math.sqrt(3.0) / 2.0


def planar_summarise(positions, weights, scan_theta_deg: float = 0.0,
                     scan_phi_deg: float = 0.0,
                     half_space: bool = False,
                     element: ElementPattern | None = None) -> dict:
    """Everything a planar-array designer wants in one call."""
    from ..core.pattern import first_sidelobe_db, hpbw_deg
    from .tapers import taper_efficiency

    pos = np.asarray(positions, dtype=float)
    w = np.asarray(weights, dtype=float)
    out: dict = {
        "elements": len(pos),
        "scan_theta_deg": scan_theta_deg,
        "scan_phi_deg": scan_phi_deg,
        "aperture_x_lambda": float(np.ptp(pos[:, 0])),
        "aperture_y_lambda": float(np.ptp(pos[:, 1])),
        "taper_efficiency": taper_efficiency(w),
    }
    d = planar_directivity(pos, w, scan_theta_deg, scan_phi_deg, half_space, element)
    out["directivity"] = d
    out["directivity_dbi"] = 10.0 * math.log10(d)
    out["half_space"] = half_space
    out["element"] = "isotropic" if element is None else element.name
    if element is not None:
        from .elements import element_directivity
        out["element_directivity_dbi"] = 10.0 * math.log10(element_directivity(
            element, math.radians(scan_theta_deg), math.radians(scan_phi_deg)))
    # Cuts through the BEAM, not at fixed azimuth: a scanned array's beam is
    # not in the phi = 90 plane, and measuring a sidelobe against a cut that
    # misses the beam reports nonsense.
    for label, plane in (("scan_plane", "scan"), ("cross_plane", "cross")):
        cut = planar_beam_cut(pos, w, plane, scan_theta_deg, scan_phi_deg, element=element)
        try:
            out[f"hpbw_{label}_deg"] = hpbw_deg(cut, 0.0)
        except ValueError:
            out[f"hpbw_{label}_deg"] = float("nan")
        try:
            out[f"sidelobe_{label}_db"] = first_sidelobe_db(cut, 0.0)
        except ValueError:
            out[f"sidelobe_{label}_db"] = float("nan")
    return out
