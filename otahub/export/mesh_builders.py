"""Solid models for the 30 archetypes the CST/HFSS exporters build no geometry
for: horns, reflectors, lenses, the travelling-wave and frequency-independent
antennas, and the patches and loop that need shapes those exporters lack.

Each builder works from the design's own dimensions - requirements, parameters
and metrics - and takes what the spec does not give (a wall thickness, an
element radius, a feed horn) from a declared construction option, with a
default the user can change (`mesh.Opt`). Every option says whether the
predicted performance accounts for it.
"""
from __future__ import annotations

import math

import numpy as np

from ..core.archetype import DesignResult
from .base import _dielectric
from .mesh import (Opt, Options, Solid3D, box, frustum, hull, m3, mesh_builder, option_values, place, plate,
                   revolve, rod, sphere, tube)

C0 = 2.99792458e8
GEO = "geometry only: the predicted performance does not depend on it"
NOT_IN = "geometry only: the predicted performance assumes an idealised version of it"


class _B:
    """A design's values and its construction options, with somewhere to put the solid."""

    def __init__(self, design: DesignResult, opts: Options, title: str):
        self.v = {**design.requirements, **design.parameters, **design.metrics}
        f = self.v.get("f0") or self.v.get("f_low") or 1e9
        self.v.setdefault("lambda0", C0 / float(f))
        self.lam = float(self.v["lambda0"])
        self.o = option_values(design, opts)
        self.copper = self.o["copper"]
        self.s = Solid3D(design.archetype, title)
        self.design = design

    def num(self, *names, default=None) -> float:
        for n in names:
            x = self.v.get(n)
            if isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(float(x)):
                return float(x)
        if default is None:
            raise KeyError(f"the design gives none of {', '.join(names)}")
        return float(default)

    def margin(self, default: float) -> float:
        m = self.o.get("margin")
        return default if m is None or not math.isfinite(m) else m

    def note(self, text: str) -> None:
        self.s.notes.append(text)

    def opt_note(self, *names: str) -> None:
        """Say which construction options this model used, and at what values."""
        parts = []
        for n in names:
            val = self.o[n]
            parts.append(f"{n} = {val * 1e3:.4g} mm" if abs(val) < 50 and n not in ("count", "turns", "dk", "tand")
                         else f"{n} = {val:.4g}")
        self.note("Construction options: " + ", ".join(parts) + " (set with --opt NAME=VALUE).")


def _lam(v) -> float:
    return float(v.get("lambda0") or C0 / float(v.get("f0") or v.get("f_low") or 1e9))


# ---------------------------------------------------------------- horns

def _rect(a, b, z):
    return [(-a / 2, -b / 2, z), (a / 2, -b / 2, z), (a / 2, b / 2, z), (-a / 2, b / 2, z)]


def _horn(B: _B, a_g, b_g, a1, b1, flare, name="horn"):
    """A hollow rectangular horn: guide (a_g x b_g) from z = -guide_length to 0,
    flare to the aperture (a1 x b1) at z = flare, walls `wall` thick."""
    t, lg = B.o["wall"], B.o["guide_length"]
    outer = hull(_rect(a_g + 2 * t, b_g + 2 * t, 0.0) + _rect(a1 + 2 * t, b1 + 2 * t, flare))
    # the cutter runs a wall's thickness past both ends, its sections extrapolated along
    # the flare, so the hollow is a_g x b_g at the throat and a1 x b1 at the aperture exactly
    ext = lambda u0, u1, z: u0 + (u1 - u0) * z / flare
    inner = hull(_rect(ext(a_g, a1, -t), ext(b_g, b1, -t), -t)
                 + _rect(ext(a_g, a1, flare + t), ext(b_g, b1, flare + t), flare + t))
    guide = box((-a_g / 2 - t, a_g / 2 + t), (-b_g / 2 - t, b_g / 2 + t), (-lg, 0.0)) - \
        box((-a_g / 2, a_g / 2), (-b_g / 2, b_g / 2), (-lg - t, t))
    return (outer - inner) + guide


_HORN_OPTS = (
    Opt("wall", "wall thickness", "m", lambda v: _lam(v) / 30, GEO),
    Opt("guide_length", "feed waveguide length", "m", lambda v: _lam(v),
        "geometry only: the predictions take the aperture field as the horn's; a short guide lets "
        "evanescent modes reach the port in a solver"),
)


@mesh_builder("pyramidal_horn", options=_HORN_OPTS)
def _pyramidal(design, opts):
    B = _B(design, opts, "Pyramidal horn")
    a_g, b_g, a1, b1, ln = (B.num(n) for n in ("a_wg", "b_wg", "a1", "b1", "p_len"))
    B.s.add("horn", "PEC", _horn(B, a_g, b_g, a1, b1, ln))
    B.note(f"Aperture {a1 * 1e3:.4g} x {b1 * 1e3:.4g} mm at z = {ln * 1e3:.4g} mm, broad wall along x; "
           "the guide ends open at z = -guide_length, where a waveguide port goes.")
    B.opt_note("wall", "guide_length")
    return B.s


@mesh_builder("e_plane_sectoral_horn", options=_HORN_OPTS + (
        Opt("b_wg", "feed guide narrow wall", "m", lambda v: float(v["a_wg"]) * 0.4444,
            "geometry only: the E-plane horn's design takes the narrow wall at the throat as given by b1's flare"),))
def _e_sectoral(design, opts):
    B = _B(design, opts, "E-plane sectoral horn")
    a, b1, rho = B.num("a_wg"), B.num("b1"), B.num("rho")
    b_g = B.o["b_wg"]
    flare = rho * (1 - b_g / b1)                        # throat to aperture, from the apex distance
    B.s.add("horn", "PEC", _horn(B, a, b_g, a, b1, flare))
    B.note(f"Flared in the E-plane (y) only, from b = {b_g * 1e3:.4g} mm to b1 = {b1 * 1e3:.4g} mm over "
           f"{flare * 1e3:.4g} mm: rho = {rho * 1e3:.4g} mm is apex to aperture.")
    B.opt_note("wall", "guide_length", "b_wg")
    return B.s


@mesh_builder("h_plane_sectoral_horn", options=_HORN_OPTS + (
        Opt("a_wg", "feed guide broad wall", "m", lambda v: float(v["b_wg"]) * 2.25,
            "geometry only: the H-plane horn's design takes the broad wall at the throat from a1's flare"),))
def _h_sectoral(design, opts):
    B = _B(design, opts, "H-plane sectoral horn")
    b, a1, rho = B.num("b_wg"), B.num("a1"), B.num("rho")
    a_g = B.o["a_wg"]
    flare = rho * (1 - a_g / a1)
    B.s.add("horn", "PEC", _horn(B, a_g, b, a1, b, flare))
    B.note(f"Flared in the H-plane (x) only, from a = {a_g * 1e3:.4g} mm to a1 = {a1 * 1e3:.4g} mm over "
           f"{flare * 1e3:.4g} mm.")
    B.opt_note("wall", "guide_length", "a_wg")
    return B.s


@mesh_builder("diagonal_horn", options=_HORN_OPTS + (
        Opt("throat", "square feed guide side", "m", lambda v: 0.75 * _lam(v),
            "geometry only: the predictions are the diagonal aperture's own field"),))
def _diagonal(design, opts):
    B = _B(design, opts, "Diagonal horn")
    a, R = B.num("a_ap"), B.num("R_axial")
    s = B.o["throat"]
    flare = R * (1 - s / a)
    B.s.add("horn", "PEC", _horn(B, s, s, a, a, flare).rotate((0.0, 0.0, 45.0)))
    B.note(f"Square aperture {a * 1e3:.4g} mm, turned 45 degrees so the E-field (along y in the guide) "
           "runs corner to corner of the aperture diagonal.")
    B.opt_note("wall", "guide_length", "throat")
    return B.s


def _cone_profile(B: _B, r_in, R, flare, extra=()):
    """Inner wall (r_in at z=0 to R at z=flare), the guide below, walls `wall` thick."""
    t, lg = B.o["wall"], B.o["guide_length"]
    return [(r_in, -lg), (r_in, 0.0), *extra, (R, flare), (R + t, flare), (r_in + t, -0.0), (r_in + t, -lg)]


_CONE_OPTS = _HORN_OPTS + (
    Opt("throat", "circular feed guide diameter", "m", lambda v: 0.7 * _lam(v),
        "geometry only: the predictions take the TE11 (or HE11) aperture field as the horn's"),)


def _cone_geometry(B: _B):
    dm, L = B.num("dm"), B.num("L")
    psi = math.asin(min(0.95, dm / (2 * L)))              # half angle; L is the slant, apex to rim
    return dm, L, psi


@mesh_builder("conical_horn", options=_CONE_OPTS)
def _conical(design, opts):
    B = _B(design, opts, "Conical horn")
    dm, L, psi = _cone_geometry(B)
    d_in = min(B.o["throat"], dm * 0.9)
    flare = (dm - d_in) / 2 / math.tan(psi)
    B.s.add("horn", "PEC", revolve(_cone_profile(B, d_in / 2, dm / 2, flare)))
    B.note(f"Aperture {dm * 1e3:.4g} mm across, half angle {math.degrees(psi):.3g} deg (slant {L * 1e3:.4g} mm "
           f"apex to rim), fed by a {d_in * 1e3:.4g} mm guide.")
    B.opt_note("wall", "guide_length", "throat")
    return B.s


@mesh_builder("corrugated_conical_horn", options=_CONE_OPTS + (
        Opt("pitch", "corrugation pitch", "m", lambda v: _lam(v) / 8,
            "geometry only: the predictions assume the balanced hybrid (HE11) field ideal corrugations make"),
        Opt("tooth", "tooth width as a fraction of the pitch", "-", lambda v: 0.25, GEO)))
def _corrugated(design, opts):
    B = _B(design, opts, "Corrugated conical horn")
    dm, L, psi = _cone_geometry(B)
    depth = B.num("slot_depth", default=B.lam / 4)
    d_in = min(B.o["throat"], dm * 0.9)
    r_in, R = d_in / 2, dm / 2
    flare = (R - r_in) / math.tan(psi)
    t, lg, p = B.o["wall"], B.o["guide_length"], B.o["pitch"]
    if not p > 0 or flare / p > 2000:
        raise ValueError(f"a {p * 1e3:.4g} mm pitch puts {flare / p if p > 0 else float('inf'):.4g} "
                         f"corrugations on the {flare * 1e3:.4g} mm flare; at most 2000 are built")
    wt = p * B.o["tooth"]
    r = lambda z: r_in + (R - r_in) * z / flare
    prof = [(r_in, -lg), (r_in, 0.0)]
    z = 0.0
    while z + p <= flare + 1e-12:
        prof += [(r(z + wt), z + wt), (r(z + wt) + depth, z + wt), (r(z + p) + depth, z + p), (r(z + p), z + p)]
        z += p
    prof += [(R, flare), (R + depth + t, flare), (r_in + depth + t, 0.0), (r_in + t, 0.0), (r_in + t, -lg)]
    B.s.add("horn", "PEC", revolve(prof))
    n = int(flare // p)
    B.note(f"{n} corrugations {depth * 1e3:.4g} mm deep (the design's slot depth) at a {p * 1e3:.4g} mm pitch "
           f"along a {flare * 1e3:.4g} mm flare.")
    B.opt_note("wall", "guide_length", "throat", "pitch")
    return B.s


@mesh_builder("conical_horn_dual_mode", options=_HORN_OPTS)
def _potter(design, opts):
    B = _B(design, opts, "Dual-mode (Potter) horn")
    dm, L, psi = _cone_geometry(B)
    d_in = B.num("d_in_over_lambda") * B.lam
    d_s = B.num("d_step")
    lp = B.num("phasing_length_m")
    r_in, r_s, R = d_in / 2, d_s / 2, dm / 2
    flare = (R - r_s) / math.tan(psi)
    t, lg = B.o["wall"], B.o["guide_length"]
    prof = [(r_in, -lg), (r_in, 0.0), (r_s, 0.0), (r_s, lp), (R, lp + flare), (R + t, lp + flare),
            (r_s + t, lp), (r_s + t, -t), (r_in + t, -t), (r_in + t, -lg)]
    B.s.add("horn", "PEC", revolve(prof))
    B.note(f"Input guide {d_in * 1e3:.4g} mm, step to {d_s * 1e3:.4g} mm at z = 0, phasing section "
           f"{lp * 1e3:.4g} mm, then the cone to {dm * 1e3:.4g} mm - the spec's first-order step and phasing; "
           "`otahub potter` solves them jointly.")
    B.opt_note("wall", "guide_length")
    return B.s


# ---------------------------------------------------------------- reflectors

def _feed_horn(B: _B, at_z: float, diameter: float, up=True):
    """A small conical feed horn with its aperture at z = at_z, looking +z (or -z)."""
    t = B.o["thickness"]
    lh = 1.5 * diameter
    r0 = 0.35 * diameter
    s = -1 if up else 1
    prof = [(r0, s * lh), (diameter / 2, 0.0), (diameter / 2 + t, 0.0), (r0 + t, s * lh)]
    if not up:
        prof = prof[::-1]
    return revolve(prof).translate((0.0, 0.0, at_z))


def _dish_shell(F: float, r_max: float, t: float, n: int = 120):
    rs = np.linspace(0.0, r_max, n)
    front = [(r, r * r / (4 * F)) for r in rs]
    back = [(r, r * r / (4 * F) - t) for r in rs[::-1]]
    return revolve(front + back, seg=96)


_DISH_OPTS = (
    Opt("thickness", "reflector thickness", "m", lambda v: max(_lam(v) / 10, float(v.get("D") or v.get("W") or 1) / 400),
        GEO),
    Opt("feed_diameter", "feed horn aperture", "m", lambda v: 2 * _lam(v),
        "geometry only: the predictions take a cos^n feed pattern set by the edge taper, not this horn"),
)


_PRIME_OPTS = (_DISH_OPTS[0], Opt("feed_diameter", "feed horn aperture", "m",
                                   lambda v: float(v.get("d_blockage") or 0.0) or 2 * _lam(v),
                                   "the predicted blockage is d_blockage, a central disc: this sets the drawn feed, "
                                   "which defaults to it"))


@mesh_builder("prime_focus_parabolic", options=_PRIME_OPTS + (
        Opt("struts", "number of feed struts", "-", lambda v: 4.0,
            "geometry only: strut blockage is not in the predictions (d_blockage is a central disc)"),
        Opt("strut_radius", "feed strut radius", "m", lambda v: _lam(v) / 4, GEO)))
def _prime(design, opts):
    B = _B(design, opts, "Prime-focus paraboloid")
    D = B.num("D")
    F = B.num("focal_length_m", default=B.num("f_over_D", default=0.4) * D)
    t = B.o["thickness"]
    B.s.add("dish", "PEC", _dish_shell(F, D / 2, t))
    fd = B.o["feed_diameter"]
    B.s.add("feed", "PEC", _feed_horn(B, F, fd, up=False))
    n = int(round(B.o["struts"]))
    zr = D * D / (16 * F)
    for k in range(n):
        a = 2 * math.pi * k / max(n, 1)
        rim = (0.45 * D * math.cos(a), 0.45 * D * math.sin(a), (0.45 * D) ** 2 / (4 * F))
        top = (fd / 2 * math.cos(a), fd / 2 * math.sin(a), F + 1.5 * fd)
        B.s.add(f"strut{k}", "PEC", rod(rim, top, B.o["strut_radius"], 12))
    B.note(f"Paraboloid D = {D * 1e3:.4g} mm, F = {F * 1e3:.4g} mm (depth {zr * 1e3:.4g} mm), vertex at the "
           "origin opening toward +z; the feed horn's aperture is at the focus, looking down.")
    B.opt_note("thickness", "feed_diameter", "strut_radius")
    return B.s


@mesh_builder("offset_parabolic", options=_DISH_OPTS)
def _offset(design, opts):
    B = _B(design, opts, "Offset paraboloid")
    D, F, h0 = B.num("D"), B.num("F"), B.num("h0")
    t = B.o["thickness"]
    parent = _dish_shell(F, h0 + D / 2 * 1.02, t, 200)
    zmax = (h0 + D / 2) ** 2 / (4 * F) + 1.0
    cut = rod((0.0, h0, -t - 1.0), (0.0, h0, zmax), D / 2, 96)
    B.s.add("dish", "PEC", parent ^ cut)
    centre = np.array([0.0, h0, h0 * h0 / (4 * F)]) - np.array([0.0, 0.0, F])
    horn = _feed_horn(B, 0.0, B.o["feed_diameter"], up=True)
    B.s.add("feed", "PEC", place(horn, (0.0, 0.0, F), centre / np.linalg.norm(centre)))
    B.note(f"The part of the parent paraboloid (F = {F * 1e3:.4g} mm) whose projected aperture is a circle "
           f"of D = {D * 1e3:.4g} mm centred {h0 * 1e3:.4g} mm off the axis; the feed at the focus looks at "
           "its centre.")
    B.opt_note("thickness", "feed_diameter")
    return B.s


def _dual(B: _B, kind: str):
    D = B.num("D")
    F = B.num("F", default=B.num("f_over_D") * D)
    M = max(B.num("magnification", default=4.0), 1.01)
    rs = B.num("Ds", default=0.2 * D) / 2
    zr = D * D / (16 * F)
    psi0 = math.atan2(D / 2, F - zr)                     # F1 to the rim, from the -z axis
    if kind == "cass":                                   # between F1 and the dish
        P = (rs, F - rs / math.tan(psi0))
    else:                                                # beyond F1, on the far rim's ray
        P = (rs, F + rs / math.tan(psi0))
    th_e = 2 * math.atan(D / (4 * M * F))
    f2 = P[1] - rs / math.tan(th_e)
    c = (F - f2) / 2
    d1 = math.hypot(P[0], P[1] - F)
    d2 = math.hypot(P[0], P[1] - f2)
    a = (d2 - d1) / 2 if kind == "cass" else (d2 + d1) / 2
    pts = []
    for i in range(61):
        phi = th_e * i / 60
        r = ((c * c - a * a) / (c * math.cos(phi) - a) if kind == "cass"
             else (a * a - c * c) / (a - c * math.cos(phi)))
        pts.append((r * math.sin(phi), f2 + r * math.cos(phi)))
    t = B.o["thickness"]
    sub = revolve(pts + [(x, z + t) for x, z in pts[::-1]], seg=96)
    B.s.add("main", "PEC", _dish_shell(F, D / 2, t))
    B.s.add("subreflector", "PEC", sub)
    B.s.add("feed", "PEC", _feed_horn(B, f2, B.o["feed_diameter"], up=True))
    B.note(f"Main paraboloid D = {D * 1e3:.4g} mm, F = {F * 1e3:.4g} mm; {'hyperboloidal' if kind == 'cass' else 'ellipsoidal'} "
           f"subreflector {2 * rs * 1e3:.4g} mm across, eccentricity {c / a:.4g} (magnification {M:g}); the feed "
           f"horn's aperture at the second focus, z = {f2 * 1e3:.4g} mm.")
    B.opt_note("thickness", "feed_diameter")
    return B.s


@mesh_builder("cassegrain", options=_DISH_OPTS)
def _cass(design, opts):
    return _dual(_B(design, opts, "Cassegrain dual reflector"), "cass")


@mesh_builder("gregorian_dual_reflector", options=_DISH_OPTS)
def _greg(design, opts):
    return _dual(_B(design, opts, "Gregorian dual reflector"), "greg")


@mesh_builder("cylindrical_parabolic", options=_DISH_OPTS[:1] + (
        Opt("feed_radius", "line feed radius", "m", lambda v: _lam(v) / 20,
            "geometry only: the predictions take a cos^n line-source feed"),))
def _cylinder(design, opts):
    B = _B(design, opts, "Parabolic cylinder")
    W, Lc, F = B.num("W"), B.num("Lc"), B.num("F")
    t = B.o["thickness"]
    xs = np.linspace(-W / 2, W / 2, 121)
    prof = [(x, x * x / (4 * F)) for x in xs] + [(x, x * x / (4 * F) - t) for x in xs[::-1]]
    sheet = m3.CrossSection([np.asarray(prof)], m3.FillRule.EvenOdd).extrude(Lc).translate((0.0, 0.0, -Lc / 2))
    B.s.add("reflector", "PEC", sheet.rotate((90.0, 0.0, 0.0)))
    B.s.add("line_feed", "PEC", rod((0.0, -Lc / 2, F), (0.0, Lc / 2, F), B.o["feed_radius"], 24))
    B.note(f"Parabolic in x-z (W = {W * 1e3:.4g} mm, F = {F * 1e3:.4g} mm), {Lc * 1e3:.4g} mm long along y; "
           "the line feed lies on the focal line.")
    B.opt_note("thickness", "feed_radius")
    return B.s


def _corner(design, opts, angle):
    B = _B(design, opts, f"{angle}-degree corner reflector")
    S, Ld, Ls = B.num("S"), B.num("Ld"), B.num("Lside")
    a = B.num("aw_over_lambda", default=0.001) * B.lam
    H, t, g = B.o["plate_height"], B.o["thickness"], B.o["feed_gap"]
    for sgn in (1, -1):
        plate_ = box((0.0, Ls), (-t / 2, t / 2), (-H / 2, H / 2)).rotate((0.0, 0.0, sgn * angle / 2))
        B.s.add(f"plate{'+' if sgn > 0 else '-'}", "PEC", plate_)
    B.s.add("dipole_top", "PEC", rod((S, 0.0, g / 2), (S, 0.0, Ld / 2), a, 16))
    B.s.add("dipole_bottom", "PEC", rod((S, 0.0, -Ld / 2), (S, 0.0, -g / 2), a, 16))
    B.note(f"Two plates {Ls * 1e3:.4g} mm wide from the apex line (the z axis) at +-{angle / 2:g} deg; a "
           f"{Ld * 1e3:.4g} mm dipole parallel to the apex, {S * 1e3:.4g} mm out on the bisector, fed across a "
           f"{g * 1e3:.3g} mm gap.")
    B.opt_note("plate_height", "thickness", "feed_gap")
    return B.s


_CORNER_OPTS = (
    Opt("plate_height", "plate length along the dipole", "m", lambda v: max(1.2 * float(v["Ld"]), 0.6 * _lam(v)),
        "geometry only: the predictions assume plates long enough to act infinite along the dipole"),
    Opt("thickness", "plate thickness", "m", lambda v: _lam(v) / 100, GEO),
    Opt("feed_gap", "dipole feed gap", "m", lambda v: _lam(v) / 100, GEO),
)


@mesh_builder("corner_reflector_90", options=_CORNER_OPTS)
def _corner90(design, opts):
    return _corner(design, opts, 90)


@mesh_builder("corner_reflector_60", options=_CORNER_OPTS)
def _corner60(design, opts):
    return _corner(design, opts, 60)


# ---------------------------------------------------------------- lenses

def _lens_curve(F, n, theta_rim, k=80):
    """Surface facing a feed at the origin: r = F(n-1)/(n cos - 1) (n > 1, hyperbola)
    or F(1-n)/(1 - n cos) (n < 1, ellipse), as (rho, z)."""
    pts = []
    for i in range(k + 1):
        th = theta_rim * i / k
        r = F * (n - 1) / (n * math.cos(th) - 1) if n > 1 else F * (1 - n) / (1 - n * math.cos(th))
        pts.append((r * math.sin(th), r * math.cos(th)))
    return pts


_LENS_FEED = Opt("feed_diameter", "feed horn aperture", "m", lambda v: 1.5 * _lam(v),
                 "geometry only: the predictions take a cos^n feed set by the edge taper")


@mesh_builder("hyperbolic_dielectric_lens", options=(
        Opt("edge_thickness", "lens thickness at the rim", "m", lambda v: _lam(v),
            "geometry only: the predictions' path lengths do not depend on it"),
        _LENS_FEED, Opt("thickness", "feed horn wall", "m", lambda v: _lam(v) / 20, GEO)))
def _hyper_lens(design, opts):
    B = _B(design, opts, "Plano-hyperbolic dielectric lens")
    F, n = B.num("F"), B.num("n_index")
    th = math.radians(B.num("theta_rim_deg"))
    pts = _lens_curve(F, n, th)
    z_back = pts[-1][1] + B.o["edge_thickness"]
    body = revolve([(0.0, F)] + pts[1:] + [(pts[-1][0], z_back), (0.0, z_back)], seg=96)
    eps = B.num("eps_r", default=n * n)
    B.s.add("lens", _dielectric(eps, design), body)
    B.s.add("feed", "PEC", _feed_horn(B, 0.0, B.o["feed_diameter"], up=True))
    B.note(f"Hyperbolic face toward the feed at the origin (vertex {F * 1e3:.4g} mm away, eccentricity n = "
           f"{n:.4g}), flat back at z = {z_back * 1e3:.4g} mm, {2 * pts[-1][0] * 1e3:.4g} mm across.")
    B.opt_note("edge_thickness", "feed_diameter")
    return B.s


@mesh_builder("metal_plate_lens", options=(
        Opt("centre_thickness", "lens thickness on the axis", "m", lambda v: _lam(v) / 4,
            "geometry only: the predictions' path lengths do not depend on it"),
        Opt("plate_thickness", "plate thickness", "m", lambda v: _lam(v) / 50,
            "geometry only: the predicted index assumes infinitely thin plates"),
        _LENS_FEED, Opt("thickness", "feed horn wall", "m", lambda v: _lam(v) / 20, GEO)))
def _plate_lens(design, opts):
    B = _B(design, opts, "Metal-plate lens")
    F, n, s = B.num("F"), B.num("n_index"), B.num("plate_spacing")
    th = math.radians(B.num("theta_rim_deg"))
    pts = _lens_curve(F, n, th)
    rho = pts[-1][0]
    z_back = F + B.o["centre_thickness"]
    body = revolve([(0.0, F)] + pts[1:] + [(rho, z_back), (0.0, z_back)], seg=96)
    tp = B.o["plate_thickness"]
    k = int(math.floor(rho / s))
    plates = []
    for i in range(-k, k + 1):
        y = i * s
        plates.append(body ^ box((-rho - 1, rho + 1), (y - tp / 2, y + tp / 2), (pts[-1][1] - 1, z_back + 1)))
    from .mesh import _union
    B.s.add("plates", "PEC", _union(plates))
    B.s.add("feed", "PEC", _feed_horn(B, 0.0, B.o["feed_diameter"], up=True))
    B.note(f"{2 * k + 1} parallel plates {s * 1e3:.4g} mm apart (normal to y: the E-field lies in them), cut to "
           f"the n = {n:.4g} lens: elliptical face toward the feed, flat back; thicker at the rim, as an index "
           "below one needs. Unzoned.")
    B.opt_note("centre_thickness", "plate_thickness", "feed_diameter")
    return B.s


# ---------------------------------------------------------------- travelling wave

@mesh_builder("yagi_uda", options=(
        Opt("element_radius", "element radius", "m", lambda v: 0.00425 * _lam(v),
            "the NBS designs behind the predictions use d/lambda = 0.0085, this default"),
        Opt("director_length", "director length", "m", lambda v: 0.425 * _lam(v),
            "geometry only: the predicted gain is the NBS design's, whose directors are 0.42-0.43 lambda"),
        Opt("boom_radius", "boom radius", "m", lambda v: 0.01 * _lam(v),
            "geometry only: the predictions assume a non-conducting or well-isolated boom"),
        Opt("feed_gap", "driven element feed gap", "m", lambda v: 0.01 * _lam(v), GEO)))
def _yagi(design, opts):
    B = _B(design, opts, "Yagi-Uda array")
    boom, lr, ld = B.num("boom_length"), B.num("reflector_length"), B.num("driven_length")
    lam = B.lam
    n_dir = max(1, int(round(B.num("director_count_nbs", default=(boom - 0.2 * lam) / (0.3 * lam)))))
    a, g = B.o["element_radius"], B.o["feed_gap"]
    xs = [0.0, 0.2 * lam] + [0.2 * lam + (boom - 0.2 * lam) * (i + 1) / n_dir for i in range(n_dir)]
    zb = -B.o["boom_radius"] - 2 * a                    # just below the elements
    B.s.add("boom", "PEC", rod((0.0, 0.0, zb), (boom, 0.0, zb), B.o["boom_radius"], 24))
    B.s.add("reflector", "PEC", rod((0.0, -lr / 2, 0.0), (0.0, lr / 2, 0.0), a, 16))
    B.s.add("driven_a", "PEC", rod((xs[1], g / 2, 0.0), (xs[1], ld / 2, 0.0), a, 16))
    B.s.add("driven_b", "PEC", rod((xs[1], -ld / 2, 0.0), (xs[1], -g / 2, 0.0), a, 16))
    for i, x in enumerate(xs[2:]):
        l = B.o["director_length"]
        B.s.add(f"director{i + 1}", "PEC", rod((x, -l / 2, 0.0), (x, l / 2, 0.0), a, 16))
    B.note(f"Boom {boom * 1e3:.4g} mm along x with the reflector, the split driven element 0.2 lambda ahead of "
           f"it and the NBS design's {n_dir} directors evenly over the rest; elements along y, the boom just "
           "below them.")
    B.opt_note("element_radius", "director_length", "boom_radius", "feed_gap")
    return B.s


@mesh_builder("lpda", options=(
        Opt("boom_spacing", "gap between the two booms", "m", lambda v: float(v["L_max"]) / 40,
            "geometry only: the predictions assume the feeder's impedance Carrel's design calls for"),
        Opt("boom_radius", "boom radius", "m", lambda v: float(v["L_max"]) / 150, GEO)))
def _lpda(design, opts):
    B = _B(design, opts, "Log-periodic dipole array")
    Lmax, tau, sigma = B.num("L_max"), B.num("tau"), B.num("sigma")
    N = max(2, int(round(B.num("N_elements"))))
    ld = B.num("length_over_diameter", default=125.0)
    gap, rb = B.o["boom_spacing"], B.o["boom_radius"]
    lens = [Lmax * tau ** i for i in range(N)]
    xs, x = [], 0.0
    for l in lens:
        xs.append(x)
        x += 2 * sigma * l
    for zb, name in ((gap / 2 + rb, "boom_top"), (-gap / 2 - rb, "boom_bottom")):
        B.s.add(name, "PEC", rod((xs[0], 0.0, zb), (xs[-1], 0.0, zb), rb, 24))
    for i, (x, l) in enumerate(zip(xs, lens)):
        r = l / ld / 2
        sgn = 1 if i % 2 == 0 else -1                  # alternate halves cross over between the booms
        B.s.add(f"el{i}_a", "PEC", rod((x, 0.0, gap / 2 + rb), (x, sgn * l / 2, gap / 2 + rb), r, 12))
        B.s.add(f"el{i}_b", "PEC", rod((x, 0.0, -gap / 2 - rb), (x, -sgn * l / 2, -gap / 2 - rb), r, 12))
    B.note(f"{N} dipoles from {lens[0] * 1e3:.4g} mm down to {lens[-1] * 1e3:.4g} mm (tau {tau:g}, sigma "
           f"{sigma:.4g}), each half on its own boom, alternating sides element to element; fed at the short end.")
    B.opt_note("boom_spacing", "boom_radius")
    return B.s


def _helix_path(D, S, turns, z0=0.0, per_turn=48):
    t = np.linspace(0.0, turns, int(turns * per_turn) + 1)
    return np.column_stack([D / 2 * np.cos(2 * np.pi * t), D / 2 * np.sin(2 * np.pi * t), z0 + S * t])


@mesh_builder("axial_mode_helix", options=(
        Opt("wire_radius", "helix wire radius", "m", lambda v: _lam(v) / 200,
            "geometry only: the predicted gain and beam are for a thin wire"),
        Opt("ground_thickness", "ground plane thickness", "m", lambda v: _lam(v) / 100, GEO)))
def _axial_helix(design, opts):
    B = _B(design, opts, "Axial-mode helix")
    D, S, N, Dg = B.num("D_helix"), B.num("S"), B.num("N"), B.num("D_gnd")
    a, tg = B.o["wire_radius"], B.o["ground_thickness"]
    z0 = S / 2
    B.s.add("ground", "PEC", rod((0.0, 0.0, -tg), (0.0, 0.0, 0.0), Dg / 2, 96))
    B.s.add("helix", "PEC", tube(_helix_path(D, S, N, z0), a))
    B.s.add("feed_lead", "PEC", rod((D / 2, 0.0, S / 8), (D / 2, 0.0, z0 + a), a, 16))
    B.note(f"{N:g} turns, {D * 1e3:.4g} mm diameter, {S * 1e3:.4g} mm pitch, starting {z0 * 1e3:.4g} mm above a "
           f"{Dg * 1e3:.4g} mm ground disc; the feed gap is the {S / 8 * 1e3:.3g} mm between the ground and the "
           "lead.")
    B.opt_note("wire_radius", "ground_thickness")
    return B.s


@mesh_builder("normal_mode_helix", options=(
        Opt("wire_radius", "helix wire radius", "m", lambda v: float(v["D_helix"]) / 40,
            "geometry only: the predicted axial ratio is the thin-wire model's"),))
def _normal_helix(design, opts):
    B = _B(design, opts, "Normal-mode helix")
    D, S, N = B.num("D_helix"), B.num("S"), B.num("N")
    a = B.o["wire_radius"]
    B.s.add("helix", "PEC", tube(_helix_path(D, S, N, 0.0), a))
    B.note(f"{N:g} turns, {D * 1e3:.4g} mm diameter, {S * 1e3:.4g} mm pitch, {N * S * 1e3:.4g} mm long along z, "
           "in free space as the spec models it; fed at either end.")
    B.opt_note("wire_radius")
    return B.s


def _wire_poly(B: _B, pts, a, name):
    for i, (p, q) in enumerate(zip(pts, pts[1:])):
        B.s.add(f"{name}{i}", "PEC", rod(p, q, a, 16))
    for i, p in enumerate(pts[1:-1]):
        B.s.add(f"{name}_joint{i}", "PEC", sphere(p, a, 16))


@mesh_builder("v_antenna_travelling", options=(
        Opt("wire_radius", "wire radius", "m", lambda v: _lam(v) / 1000, GEO),
        Opt("feed_gap", "feed gap at the apex", "m", lambda v: _lam(v) / 100, GEO)))
def _vee(design, opts):
    B = _B(design, opts, "V antenna")
    L, half = B.num("L"), math.radians(B.num("half_angle_deg"))
    a, g = B.o["wire_radius"], B.o["feed_gap"]
    for sgn in (1, -1):
        p0 = (0.0, sgn * g / 2, 0.0)
        B.s.add(f"leg{'+' if sgn > 0 else '-'}", "PEC",
                rod(p0, (L * math.cos(half), sgn * (g / 2 + L * math.sin(half)), 0.0), a, 16))
    B.note(f"Two {L * 1e3:.4g} mm legs at +-{math.degrees(half):.4g} deg from the x axis, fed across the apex gap.")
    B.opt_note("wire_radius", "feed_gap")
    return B.s


@mesh_builder("rhombic", options=(
        Opt("feed_gap", "feed and load gaps", "m", lambda v: _lam(v) / 100, GEO),))
def _rhombic(design, opts):
    B = _B(design, opts, "Rhombic antenna")
    L, half = B.num("L"), math.radians(B.num("half_angle_deg"))
    a, g = B.num("aw", default=B.lam / 1000), B.o["feed_gap"]
    dx, dy = L * math.cos(half), L * math.sin(half)
    for sgn in (1, -1):
        _wire_poly(B, [(0.0, sgn * g / 2, 0.0), (dx, sgn * dy, 0.0), (2 * dx, sgn * g / 2, 0.0)], a,
                   f"side{'+' if sgn > 0 else '-'}")
    B.note(f"Four {L * 1e3:.4g} mm legs, half angle {math.degrees(half):.4g} deg, {2 * dx * 1e3:.4g} mm long; "
           "fed across the gap at x = 0, the terminating resistor goes across the gap at the far end.")
    B.opt_note("feed_gap")
    return B.s


@mesh_builder("long_wire_travelling", options=(
        Opt("height", "height above ground", "m", lambda v: float(v.get("height") or _lam(v) / 4),
            "defaults to the design's own height, which sets its predicted termination resistance; change it "
            "there (--set height=...) for the predictions to follow - this option only places the wire"),
        Opt("ground_width", "ground plane width", "m", lambda v: _lam(v), "geometry only: an infinite ground is the "
            "usual assumption"),
        Opt("feed_gap", "feed and load gaps", "m", lambda v: _lam(v) / 100, GEO)))
def _long_wire(design, opts):
    B = _B(design, opts, "Terminated long wire")
    L = B.num("L")
    a = B.num("aw", default=B.lam / 1000)
    h, wg, g = B.o["height"], B.o["ground_width"], B.o["feed_gap"]
    B.s.add("ground", "PEC", box((-wg / 2, L + wg / 2), (-wg / 2, wg / 2), (-B.lam / 200, 0.0)))
    _wire_poly(B, [(0.0, 0.0, g), (0.0, 0.0, h), (L, 0.0, h), (L, 0.0, g)], a, "wire")
    B.note(f"A {L * 1e3:.4g} mm wire {h * 1e3:.4g} mm above a finite ground, fed across the gap at x = 0 and "
           "loaded across the gap at x = L. The spec's height sets its termination resistance; its pattern "
           "and directivity are the wire's alone, without the ground.")
    B.opt_note("height", "ground_width", "feed_gap")
    return B.s


# ---------------------------------------------------------------- printed / frequency independent

_SUB_OPTS = (
    Opt("substrate_h", "substrate thickness (0 = free-standing copper)", "m", lambda v: 0.0,
        "geometry only: the predictions are for the antenna in free space"),
    Opt("dk", "substrate permittivity (Dk)", "-", lambda v: 4.4,
        "geometry only: a substrate lowers the frequencies the predictions give"),
    Opt("tand", "substrate loss tangent", "-", lambda v: 0.0, GEO),
)


def _board(B: _B, outline_lo, outline_hi, default_margin):
    """A substrate under the copper (z from -h to 0) when one is asked for."""
    h = B.o.get("substrate_h", 0.0)
    if not h or h <= 0:
        return
    m = B.margin(default_margin)
    tand = B.o.get("tand", 0.0)
    mat = f"eps_r={B.o['dk']:g}" + (f";tand={tand:g}" if tand > 0 else "")
    B.s.add("substrate", mat, box((outline_lo[0] - m, outline_hi[0] + m), (outline_lo[1] - m, outline_hi[1] + m),
                                  (-h, 0.0)))
    B.note(f"On a {h * 1e3:.4g} mm substrate, Dk {B.o['dk']:g}" + (f", tan d {tand:g}" if tand > 0 else "")
           + f", {m * 1e3:.4g} mm margin.")


@mesh_builder("bowtie", options=(
        Opt("feed_gap", "feed gap between the arms", "m", lambda v: 0.08 * float(v["Lside"]), GEO),) + _SUB_OPTS)
def _bowtie(design, opts):
    B = _B(design, opts, "Bowtie")
    Ls, We = B.num("Lside"), B.num("Wend")
    g = B.o["feed_gap"] / 2
    for sgn in (1, -1):
        tri = [(sgn * g, 0.0), (sgn * (g + Ls), We / 2), (sgn * (g + Ls), -We / 2)]
        if sgn < 0:
            tri = tri[::-1]
        B.s.add(f"arm{'+' if sgn > 0 else '-'}", "PEC", plate(tri, 0.0, B.copper))
    _board(B, (-g - Ls, -We / 2), (g + Ls, We / 2), 0.1 * Ls)
    B.note(f"Two triangular arms {Ls * 1e3:.4g} mm long, {We * 1e3:.4g} mm across the ends, fed across the "
           "central gap; copper on z = 0.")
    B.opt_note("feed_gap", "copper", "substrate_h")
    return B.s


def _arm_polygon(centre, width):
    """A strip of the given width along a 2-D centreline (offset along its normals)."""
    c = np.asarray(centre)
    d = np.gradient(c, axis=0)
    nrm = np.column_stack([-d[:, 1], d[:, 0]])
    nrm /= np.linalg.norm(nrm, axis=1)[:, None]
    w = np.broadcast_to(np.asarray(width, float), (len(c),))[:, None]
    return np.vstack([c + nrm * w / 2, (c - nrm * w / 2)[::-1]])


@mesh_builder("archimedean_spiral", options=_SUB_OPTS)
def _arch_spiral(design, opts):
    B = _B(design, opts, "Archimedean spiral")
    ri, ro = B.num("r_in"), B.num("r_out")
    turns = B.num("n_turns", default=5.0)
    pitch = (ro - ri) / turns                         # one arm's growth per turn
    w = pitch / 4                                     # two arms, two gaps: self-complementary
    for arm in (0.0, math.pi):
        phi = np.linspace(0.0, 2 * np.pi * turns, int(turns * 120) + 1)
        r = ri + pitch * phi / (2 * np.pi)
        centre = np.column_stack([r * np.cos(phi + arm), r * np.sin(phi + arm)])
        B.s.add(f"arm{int(arm > 0) + 1}", "PEC", plate(_arm_polygon(centre, w), 0.0, B.copper))
    _board(B, (-ro, -ro), (ro, ro), 0.1 * ro)
    B.note(f"Two arms, {turns:g} turns from r = {ri * 1e3:.4g} to {ro * 1e3:.4g} mm, each {w * 1e3:.3g} mm wide "
           "with equal gaps (self-complementary); fed between the arms' inner ends.")
    B.opt_note("copper", "substrate_h")
    return B.s


@mesh_builder("equiangular_spiral", options=_SUB_OPTS)
def _equi_spiral(design, opts):
    B = _B(design, opts, "Equiangular spiral")
    r0, ro = B.num("r0"), B.num("r_out")
    a = B.num("a_growth", default=0.221)
    delta = math.pi / 2                               # arm and gap equal: self-complementary
    pmax = math.log(ro / r0) / a
    for arm in (0.0, math.pi):
        phi = np.linspace(0.0, pmax, 600)
        outer = np.column_stack([r0 * np.exp(a * phi) * np.cos(phi + arm), r0 * np.exp(a * phi) * np.sin(phi + arm)])
        inner = np.column_stack([r0 * np.exp(a * (phi - delta)) * np.cos(phi + arm),
                                 r0 * np.exp(a * (phi - delta)) * np.sin(phi + arm)])
        B.s.add(f"arm{int(arm > 0) + 1}", "PEC", plate(np.vstack([outer, inner[::-1]]), 0.0, B.copper))
    _board(B, (-ro, -ro), (ro, ro), 0.1 * ro)
    B.note(f"Two arms r = r0 exp(a phi) with a = {a:.4g}, from r0 = {r0 * 1e3:.4g} to {ro * 1e3:.4g} mm "
           f"({pmax / (2 * math.pi):.3g} turns), each arm the region between two such curves 90 deg apart.")
    B.opt_note("copper", "substrate_h")
    return B.s


@mesh_builder("vivaldi_tsa", options=(
        Opt("substrate_h", "substrate thickness", "m", lambda v: 0.813e-3,
            "geometry only: the predictions do not model the board (its Dk lowers the band edge)"),
        Opt("dk", "substrate permittivity (Dk)", "-", lambda v: 3.55, NOT_IN),
        Opt("tand", "substrate loss tangent", "-", lambda v: 0.0027, GEO),
        Opt("back_length", "board behind the taper (slotline and cavity)", "m", lambda v: 0.2 * float(v["Lax"]),
            "geometry only: the predictions assume an ideal broadband slotline termination"),
        Opt("cavity_radius", "slotline cavity radius", "m", lambda v: float(v["W_ap"]) / 12,
            "geometry only: the predictions assume an ideal broadband slotline termination")))
def _vivaldi(design, opts):
    B = _B(design, opts, "Vivaldi tapered-slot antenna")
    L, Wap, ws = B.num("Lax"), B.num("W_ap"), B.num("w_slot", default=3e-4)
    R = B.num("R_open", default=10.0)
    m = B.margin(0.1 * Wap)
    W = Wap + 2 * m
    back = B.o["back_length"]
    x0 = -back / 2                                    # the cavity's centre, behind the taper
    s0 = ws / 2
    # the spec's taper, y = C1 exp(R z) + C2 from the throat (z = 0, y = s0) to the
    # aperture (z = Lax, y = W_ap/2): R_open sets how fast it opens, the ends fix C1, C2
    c1 = (Wap / 2 - s0) / math.expm1(R * L)
    c2 = s0 - c1
    xs = np.linspace(0.0, L, 160)
    top = [(x, c1 * math.exp(R * x) + c2) for x in xs]
    slot = [(x0, s0), *top, (L * 1.001, top[-1][1]), (L * 1.001, -top[-1][1]),
            *[(x, -y) for x, y in top[::-1]], (x0, -s0)]
    metal = plate([(-back, -W / 2), (L, -W / 2), (L, W / 2), (-back, W / 2)], 0.0, B.copper) \
        - plate(slot, -B.copper, 3 * B.copper) - rod((x0, 0.0, -B.copper), (x0, 0.0, 2 * B.copper), B.o["cavity_radius"], 48)
    B.s.add("metal", "PEC", metal)
    tand = B.o["tand"]
    B.s.add("substrate", f"eps_r={B.o['dk']:g}" + (f";tand={tand:g}" if tand > 0 else ""),
            box((-back, L), (-W / 2, W / 2), (-B.o["substrate_h"], 0.0)))
    B.note(f"Coplanar Vivaldi: the taper y = C1 exp(R z) + C2 with the design's R_open = {R:g} /m opens from a "
           f"{ws * 1e3:.3g} mm slot at z = 0 to {Wap * 1e3:.4g} mm at z = Lax = {L * 1e3:.4g} mm (growth "
           f"{math.exp(R * L):.4g}, the spec's taper ratio); behind it a slotline to a circular cavity. Board "
           f"{(L + back) * 1e3:.4g} x {W * 1e3:.4g} mm; the feed crosses the slot at z = 0.")
    B.opt_note("substrate_h", "dk", "back_length", "cavity_radius", "copper")
    return B.s


# ---------------------------------------------------------------- patches and a loop

def _patch_on_board(B: _B, outline, probe_xy, size):
    eps, h = B.num("eps_r"), B.num("h")
    m = B.margin(0.6 * size)
    xs, ys = [p[0] for p in outline], [p[1] for p in outline]
    lo, hi = (min(xs) - m, min(ys) - m), (max(xs) + m, max(ys) + m)
    rp = B.o["probe_radius"]
    hole = 2.3 * rp                                    # a 50-ohm coax outer for a PTFE-filled line
    px, py = probe_xy
    B.s.add("substrate", _dielectric(eps, B.design), box((lo[0], hi[0]), (lo[1], hi[1]), (0.0, h)))
    B.s.add("ground", "PEC", box((lo[0], hi[0]), (lo[1], hi[1]), (-B.copper, 0.0))
            - rod((px, py, -2 * B.copper), (px, py, B.copper), hole, 32))
    B.s.add("patch", "PEC", plate(outline, h, B.copper))
    B.s.add("probe", "PEC", rod((px, py, -B.copper), (px, py, h), rp, 24))
    B.note(f"On {h * 1e3:.4g} mm of eps_r {eps:g} (the design's own substrate), board margin {m * 1e3:.4g} mm; "
           f"copper on both faces. The probe passes through a {2 * hole * 1e3:.3g} mm clearance hole in the "
           "ground (a coaxial feed): the port goes between the probe and the hole's rim.")


_PATCH_OPTS = (
    Opt("probe_radius", "feed probe radius", "m", lambda v: 0.65e-3,
        "geometry only: the predictions do not include the probe's reactance"),
    Opt("probe_offset", "feed probe offset from the centre", "m", lambda v: float("nan"),
        "geometry only: the spec gives no feed position - tune it for 50 ohm in a solver"),
)


@mesh_builder("triangular_patch", options=_PATCH_OPTS)
def _tri_patch(design, opts):
    B = _B(design, opts, "Equilateral triangular patch")
    s = B.num("a_side")
    ht = s * math.sqrt(3) / 2
    outline = [(-s / 2, -ht / 3), (s / 2, -ht / 3), (0.0, 2 * ht / 3)]     # centroid at the origin
    off = B.o["probe_offset"]
    off = 0.2 * ht if not math.isfinite(off) else off
    _patch_on_board(B, outline, (0.0, -off), s)
    B.note(f"Side {s * 1e3:.4g} mm, centroid at the origin, a base edge along x; probe {off * 1e3:.4g} mm "
           "below the centroid on the altitude.")
    B.opt_note("probe_radius", "copper")
    return B.s


@mesh_builder("truncated_corner_cp_patch", options=_PATCH_OPTS)
def _cp_patch(design, opts):
    B = _B(design, opts, "Truncated-corner CP patch")
    L, c = B.num("L"), B.num("c_trunc")
    h2 = L / 2
    outline = [(-h2 + c, h2), (h2, h2), (h2, -h2 + c), (h2 - c, -h2), (-h2, -h2), (-h2, h2 - c)]
    off = B.o["probe_offset"]
    off = 0.22 * L if not math.isfinite(off) else off
    _patch_on_board(B, outline, (0.0, -off), L)
    B.note(f"Square {L * 1e3:.4g} mm with two opposite corners cut back {c * 1e3:.3g} mm along each edge; probe "
           f"{off * 1e3:.4g} mm off centre on the y axis.")
    B.opt_note("probe_radius", "copper")
    return B.s


@mesh_builder("ferrite_rod_loop", options=())
def _ferrite(design, opts):
    B = _B(design, opts, "Ferrite-rod loop")
    l, d = B.num("l_rod"), B.num("d_rod")
    N = B.num("N", default=60)
    lc = B.num("l_coil", default=0.35 * l)
    b = B.num("b", default=2e-4)
    mu = B.num("mu_i", default=125)
    B.s.add("rod", f"ferrite_mu_i={mu:g}", rod((-l / 2, 0.0, 0.0), (l / 2, 0.0, 0.0), d / 2, 48))
    per_layer = max(1, int(lc // (2.05 * b)))          # turns that fit side by side along l_coil
    layers = int(math.ceil(N / per_layer))
    left = N
    for k in range(layers):
        n_k = min(per_layer, left)
        left -= n_k
        rc = d / 2 + b * 1.05 + k * 2.1 * b
        length = lc * n_k / per_layer if layers > 1 else lc
        pitch = length / n_k
        t = np.linspace(0.0, n_k, int(n_k * 24) + 1)
        path = np.column_stack([-length / 2 + pitch * t, rc * np.cos(2 * np.pi * t), rc * np.sin(2 * np.pi * t)])
        B.s.add(f"coil{k + 1}", "PEC", tube(path, b, 10))
    B.note(f"Ferrite rod {l * 1e3:.4g} x {d * 1e3:.4g} mm (initial permeability {mu:g}) along x, wound with "
           f"{N:g} turns of {2 * b * 1e3:.3g} mm wire over the design's l_coil = {lc * 1e3:.4g} mm at its middle"
           + (f", in {layers} layers since they do not fit in one" if layers > 1 else "")
           + ". The ferrite file carries the rod; give it the material's permeability and loss in the solver.")
    return B.s
