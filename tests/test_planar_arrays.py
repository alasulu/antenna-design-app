"""Planar array checks.

The decisive tests here are cross-checks rather than quoted numbers: a
separable rectangular array must reproduce the linear module exactly in its
principal planes, the exact directivity kernel must agree with brute-force
integration over the sphere, and the triangular lattice's advantage must fall
out of the reciprocal lattice rather than being asserted.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from otahub import arrays as A


def _brute_force_directivity(pos, w, n=721):
    """4*pi*U_max / P_rad by direct integration, to check the closed kernel."""
    th = np.linspace(1e-9, np.pi - 1e-9, n)
    ph = np.linspace(0, 2 * np.pi, 2 * n)
    T, P = np.meshgrid(th, ph, indexing="ij")
    u = np.sin(T) * np.cos(P)
    v = np.sin(T) * np.sin(P)
    af = np.zeros(T.shape, dtype=complex)
    for (x, y), wi in zip(np.asarray(pos), np.asarray(w, dtype=complex)):
        af += wi * np.exp(2j * np.pi * (x * u + y * v))
    power = np.abs(af) ** 2
    p_rad = np.trapezoid(np.trapezoid(power * np.sin(T), ph, axis=1), th, axis=0)
    return 4 * np.pi * power.max() / p_rad


# ----------------------------------------------------------------- lattices

def test_rectangular_lattice_is_centred_and_ordered():
    pos = A.rectangular_lattice(3, 2, 0.5, 0.25)
    assert pos.shape == (6, 2)
    assert np.allclose(pos.mean(axis=0), 0.0)
    assert np.allclose(np.unique(pos[:, 0]), [-0.5, 0.0, 0.5])
    assert np.allclose(np.unique(pos[:, 1]), [-0.125, 0.125])


def test_triangular_lattice_neighbours_are_all_the_same_distance():
    """Every interior element must have six neighbours at exactly s."""
    s = 0.6
    pos = A.triangular_lattice(7, 7, s)
    centre = np.argmin(np.linalg.norm(pos, axis=1))
    d = np.linalg.norm(pos - pos[centre], axis=1)
    near = np.sort(d)[1:7]
    assert np.allclose(near, s, rtol=1e-12), f"neighbour distances {near}"


@pytest.mark.parametrize("bad", [(0, 4, 0.5), (4, 0, 0.5), (4, 4, 0.0), (4, 4, -1.0)])
def test_lattices_reject_degenerate_input(bad):
    nx, ny, d = bad
    with pytest.raises(ValueError):
        A.rectangular_lattice(nx, ny, d)
    with pytest.raises(ValueError):
        A.triangular_lattice(nx, ny, d)


# -------------------------------------------------------------- directivity

@pytest.mark.parametrize("nx,ny,d", [(2, 2, 0.5), (4, 4, 0.5), (4, 4, 0.7), (6, 3, 0.55)])
def test_directivity_kernel_matches_brute_force_integration(nx, ny, d):
    pos = A.rectangular_lattice(nx, ny, d)
    w = np.ones(nx * ny)
    assert A.planar_directivity(pos, w) == pytest.approx(
        _brute_force_directivity(pos, w), rel=2e-3)


def test_triangular_directivity_also_matches_brute_force():
    pos = A.triangular_lattice(4, 4, 0.5)
    w = np.ones(len(pos))
    assert A.planar_directivity(pos, w) == pytest.approx(
        _brute_force_directivity(pos, w), rel=2e-3)


def test_single_row_reproduces_the_linear_module_exactly():
    """A one-row planar array IS a linear array; the two kernels must agree."""
    for n, d in ((8, 0.5), (16, 0.5), (10, 0.75)):
        pos = A.rectangular_lattice(1, n, d)
        planar = A.planar_directivity(pos, np.ones(n))
        linear = A.broadside_directivity(A.uniform(n), d)
        assert planar == pytest.approx(linear, rel=1e-12)


def test_uniform_half_wave_array_approaches_two_pi_times_the_aperture():
    """Isotropic elements radiate both ways, so a big planar array reaches
    2*pi*A/lambda^2 — half the 4*pi*A/lambda^2 a one-sided aperture gets."""
    pos = A.rectangular_lattice(20, 20, 0.5)
    d = A.planar_directivity(pos, np.ones(400))
    area = (19 * 0.5) ** 2                     # aperture between outer elements
    assert d == pytest.approx(2 * math.pi * area, rel=0.10)
    # and with a ground plane it is exactly twice that
    assert A.planar_directivity(pos, np.ones(400), half_space=True) == pytest.approx(
        2 * d, rel=1e-12)


def test_scan_loss_follows_the_projected_aperture():
    pos = A.rectangular_lattice(16, 16, 0.5)
    w = np.ones(256)
    broadside = A.planar_directivity(pos, w)
    for angle in (15.0, 30.0, 45.0):
        scanned = A.planar_directivity(pos, w, scan_theta_deg=angle)
        assert scanned / broadside == pytest.approx(
            math.cos(math.radians(angle)), rel=0.05), f"at {angle} deg"


def test_steering_does_not_change_the_peak_value():
    """The peak is at the scan angle, where the steering phase cancels. A
    directivity that collapses when the beam is steered means the numerator is
    being evaluated at broadside instead."""
    pos = A.rectangular_lattice(8, 8, 0.5)
    w = np.ones(64)
    assert A.planar_directivity(pos, w, scan_theta_deg=45.0) > 0.5 * A.planar_directivity(pos, w)


# ------------------------------------------------------------- separability

def test_principal_planes_reproduce_the_linear_array_exactly():
    """A separable rectangular array's principal-plane cuts ARE the linear
    array's pattern. This is the property that makes planar design tractable."""
    from otahub.core.pattern import first_sidelobe_db, hpbw_deg

    for n in (8, 16):
        pos = A.rectangular_lattice(n, n, 0.5)
        w = np.ones(n * n)
        cut = A.planar_pattern_cut(pos, w, phi_deg=0.0)
        linear = A.summarise(A.uniform(n), 0.5)
        # Not bit-identical: the two cuts use different theta grids, so the
        # half-power crossing is interpolated at a different pitch. 1e-4
        # relative is still four orders tighter than any real discrepancy.
        assert hpbw_deg(cut, 0.0) == pytest.approx(linear["hpbw_deg"], rel=1e-4)
        assert first_sidelobe_db(cut, 0.0) == pytest.approx(
            linear["sidelobe_db"], abs=0.02)


@pytest.mark.parametrize("sll", [-20.0, -30.0, -40.0])
def test_separable_chebyshev_holds_its_design_level_in_both_planes(sll):
    from otahub.core.pattern import first_sidelobe_db

    n = 12
    wx = A.dolph_chebyshev(n, sll)
    pos = A.rectangular_lattice(n, n, 0.5)
    w = A.separable_weights(wx, wx)
    for phi in (0.0, 90.0):
        cut = A.planar_pattern_cut(pos, w, phi_deg=phi)
        assert first_sidelobe_db(cut, math.radians(phi)) == pytest.approx(sll, abs=0.05)


def test_diagonal_plane_sidelobes_are_lower_than_the_principal_planes():
    """The two array factors multiply off the principal planes, so the
    diagonal sidelobes sit at roughly twice the design level in dB."""
    from otahub.core.pattern import first_sidelobe_db

    n = 12
    wx = A.dolph_chebyshev(n, -25.0)
    pos = A.rectangular_lattice(n, n, 0.5)
    w = A.separable_weights(wx, wx)
    principal = first_sidelobe_db(A.planar_pattern_cut(pos, w, 0.0), 0.0)
    diagonal = first_sidelobe_db(A.planar_pattern_cut(pos, w, 45.0), math.radians(45.0))
    assert diagonal < principal - 10.0, f"principal {principal}, diagonal {diagonal}"


# ------------------------------------------------------------ grating lobes

@pytest.mark.parametrize("scan", [0.0, 30.0, 45.0, 60.0])
def test_grating_lobe_spacing_matches_the_reciprocal_lattice(scan):
    """The limit is that the shortest reciprocal-lattice vector exceeds
    (1 + sin(scan))*k0. Recomputing it from the primitive vectors is an
    independent route to the same number."""
    limit = 1.0 + math.sin(math.radians(scan))

    def shortest_reciprocal(a1, a2):
        B = 2 * np.pi * np.linalg.inv(np.array([a1, a2]).T)   # rows are b1, b2
        return min(np.linalg.norm(m * B[0] + n * B[1])
                   for m in range(-3, 4) for n in range(-3, 4) if (m, n) != (0, 0))

    d = A.grating_lobe_free_spacing_planar(scan, "rectangular")
    assert shortest_reciprocal([d, 0], [0, d]) / (2 * np.pi) == pytest.approx(limit, rel=1e-9)

    s = A.grating_lobe_free_spacing_planar(scan, "triangular")
    assert shortest_reciprocal(
        [s, 0], [s / 2, s * math.sqrt(3) / 2]) / (2 * np.pi) == pytest.approx(limit, rel=1e-9)


def test_triangular_lattice_saves_thirteen_percent_of_the_elements():
    saving = A.lattice_element_saving()
    assert saving == pytest.approx(1 - math.sqrt(3) / 2, rel=1e-12)
    assert saving == pytest.approx(0.1340, abs=5e-4)
    # and it falls out of the cell areas, at every scan angle
    for scan in (0.0, 30.0, 60.0):
        d = A.grating_lobe_free_spacing_planar(scan, "rectangular")
        s = A.grating_lobe_free_spacing_planar(scan, "triangular")
        cell_rect = d * d
        cell_tri = math.sqrt(3) / 2 * s * s
        assert cell_tri / cell_rect == pytest.approx(2 / math.sqrt(3), rel=1e-12)
        assert 1 - cell_rect / cell_tri == pytest.approx(saving, rel=1e-12)


def test_triangular_lattice_at_its_limit_really_has_no_grating_lobe():
    """Brute force: scan to the design limit and check nothing far from the
    beam comes back up to full amplitude."""
    scan = 45.0
    s = A.grating_lobe_free_spacing_planar(scan, "triangular")
    pos = A.triangular_lattice(16, 16, s)
    w = np.ones(len(pos))
    psi = np.linspace(0, np.pi / 2, 300)
    phi = np.linspace(0, 2 * np.pi, 480)
    T, P = np.meshgrid(psi, phi, indexing="ij")
    af = A.planar_array_factor(pos, w, T, P, scan_theta_deg=scan)
    mag = np.abs(af) / len(pos)
    u = np.sin(T) * np.cos(P)
    v = np.sin(T) * np.sin(P)
    u0 = math.sin(math.radians(scan))
    far = np.hypot(u - u0, v) > 0.3
    assert 20 * math.log10(mag[far].max()) < -8.0


def test_unknown_lattice_is_rejected():
    with pytest.raises(ValueError):
        A.grating_lobe_free_spacing_planar(0.0, "hexagonalish")


# ---------------------------------------------------------------- summarise

def test_summarise_reports_a_usable_set():
    pos = A.rectangular_lattice(8, 8, 0.5)
    s = A.planar_summarise(pos, np.ones(64))
    assert s["elements"] == 64
    assert s["aperture_x_lambda"] == pytest.approx(3.5)
    assert s["directivity_dbi"] == pytest.approx(19.7368, abs=0.01)
    assert s["hpbw_scan_plane_deg"] == pytest.approx(s["hpbw_cross_plane_deg"], rel=1e-9)
    assert math.isfinite(s["sidelobe_scan_plane_db"])


# ------------------------------------------------------- cuts through the beam

def _hemisphere_peak_away_from_beam(pos, w, scan_deg, exclude=0.3):
    """Largest normalised |AF| in the forward hemisphere, well away from the beam."""
    psi = np.linspace(0, np.pi / 2, 420)
    phi = np.linspace(0, 2 * np.pi, 720)
    T, P = np.meshgrid(psi, phi, indexing="ij")
    af = A.planar_array_factor(pos, w, T, P, scan_theta_deg=scan_deg)
    mag = np.abs(af) / np.sum(np.abs(w))
    u = np.sin(T) * np.cos(P)
    v = np.sin(T) * np.sin(P)
    u0 = math.sin(math.radians(scan_deg))
    far = np.hypot(u - u0, v) > exclude
    return mag[far].max()


def test_grating_lobe_appears_just_past_the_limit_and_not_before():
    """The converse of the limit test: step across it and a full-amplitude lobe
    must appear somewhere in the hemisphere."""
    scan = 45.0
    limit = A.grating_lobe_free_spacing_planar(scan, "triangular")
    inside = A.triangular_lattice(14, 14, limit * 0.90)
    outside = A.triangular_lattice(14, 14, limit * 1.25)
    lo = _hemisphere_peak_away_from_beam(inside, np.ones(len(inside)), scan)
    hi = _hemisphere_peak_away_from_beam(outside, np.ones(len(outside)), scan)
    assert 20 * math.log10(lo) < -8.0, "grating lobe inside the limit"
    assert 20 * math.log10(hi) > -3.0, "no grating lobe beyond the limit"


def test_scan_plane_cut_stays_in_the_forward_hemisphere():
    """Sweeping a great circle outward from a scanned beam runs past the array
    plane into the mirror beam, which is at full amplitude for isotropic
    elements. The sidelobe reading then pins to 0 dB at every spacing."""
    from otahub.core.pattern import first_sidelobe_db

    pos = A.triangular_lattice(12, 12, 0.55)
    w = np.ones(len(pos))
    for scan in (0.0, 30.0, 60.0):
        cut = A.planar_beam_cut(pos, w, "scan", scan_theta_deg=scan)
        sll = first_sidelobe_db(cut, 0.0)
        assert sll < -10.0, f"scan {scan} deg reported {sll:.2f} dB"


def test_scan_plane_beam_broadens_as_one_over_cos_and_cross_plane_does_not():
    pos = A.rectangular_lattice(16, 16, 0.5)
    w = np.ones(256)
    base = A.planar_summarise(pos, w)
    for scan in (30.0, 45.0):
        s = A.planar_summarise(pos, w, scan_theta_deg=scan)
        assert s["hpbw_scan_plane_deg"] == pytest.approx(
            base["hpbw_scan_plane_deg"] / math.cos(math.radians(scan)), rel=0.05)
        assert s["hpbw_cross_plane_deg"] == pytest.approx(
            base["hpbw_cross_plane_deg"], rel=1e-3)


def test_broadside_beam_cut_equals_the_principal_azimuth_cut():
    """At broadside the scan plane is undefined by geometry, so it falls back
    on the requested azimuth. The two routes must then agree."""
    from otahub.core.pattern import first_sidelobe_db, hpbw_deg

    pos = A.rectangular_lattice(12, 12, 0.5)
    w = np.ones(144)
    a = A.planar_beam_cut(pos, w, "scan", 0.0, 0.0)
    b = A.planar_pattern_cut(pos, w, phi_deg=0.0)
    assert hpbw_deg(a, 0.0) == pytest.approx(hpbw_deg(b, 0.0), rel=1e-9)
    assert first_sidelobe_db(a, 0.0) == pytest.approx(first_sidelobe_db(b, 0.0), abs=1e-9)


def test_beam_cut_rejects_an_unknown_plane():
    pos = A.rectangular_lattice(4, 4, 0.5)
    with pytest.raises(ValueError):
        A.planar_beam_cut(pos, np.ones(16), plane="diagonal")
