"""Solid 3-D models of every archetype as watertight triangle meshes, written as STL.

The CST and HFSS exporters describe geometry to a solver's own modeller. This
module builds the solids itself, on manifold3d's mesh kernel - closed,
consistently oriented meshes and real booleans - so a design can go to any CAD
tool, 3-D viewer, printer or solver that imports STL. Every one of the 72
archetypes has a model: the 42 with a CST/HFSS builder are converted from that
model, the other 30 are built here (`mesh_builders`).

Conventions:
- Metres internally; the STL is written in MILLIMETRES, the usual CAD unit -
  choose millimetres on import.
- Zero-thickness conductors (patches, ground planes, strips) become sheets
  `SHEET` thick (35 um, one-ounce copper), since a mesh cannot hold a surface.
- One STL per material: conductor, each dielectric, ferrite. Where a conductor
  and a dielectric overlap, the conductor wins and the dielectric is cut, so
  the files can be imported together without intersecting.
- Air regions and cutting tools (materials VOID and VACUUM) are not written.
- What a builder had to choose, because the spec does not give it (a wall
  thickness, a feed guide's length, a wire radius), is listed in `notes`.

STEP is not offered: a B-rep writer needs a CAD kernel (OpenCascade, several
hundred megabytes); STL carries the same geometry as facets and is read by
CST, HFSS, FreeCAD, Fusion 360 and SolidWorks alike.
"""
from __future__ import annotations

import math
import struct
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Sequence

import numpy as np

try:                                   # an optional dependency: pip install manifold3d
    import manifold3d as m3
except ImportError:                    # pragma: no cover - reported when used
    m3 = None

from ..core.archetype import DesignResult

__all__ = ["SHEET", "Opt", "COMMON", "Options", "preview", "option_values", "Body", "Solid3D", "solid", "write_stl", "MESH_BUILDERS",
           "MESH_OPTIONS", "construction_options", "available"]

#: thickness given to a zero-thickness conductor by default [m]: one-ounce copper
SHEET = 35e-6
#: facets around a full circle
SEG = 48
#: scales every facet count: 1 for files, lower for a quick on-screen preview
_DETAIL = [1.0]


class preview:
    """Within `with preview():` circles get fewer facets - for drawing on screen,
    not for files."""

    def __init__(self, detail: float = 0.35):
        self.detail = detail

    def __enter__(self):
        self.saved = _DETAIL[0]
        _DETAIL[0] = self.detail
        return self

    def __exit__(self, *exc):
        _DETAIL[0] = self.saved


def _seg(n: int) -> int:
    return max(8, int(round(n * _DETAIL[0])))
#: materials that are air or cutting tools, never written
_NOT_SOLID = ("VOID", "VACUUM")


def available() -> bool:
    return m3 is not None


def _need() -> None:
    if m3 is None:
        raise ImportError("3-D export needs the manifold3d package: pip install manifold3d")


# ---------------------------------------------------------------- options

@dataclass(frozen=True)
class Opt:
    """A construction choice the electrical design does not fix - a wall's
    thickness, a wire's radius, a feed guide's length - with its default.

    `default` takes the design's values (requirements, parameters, metrics and
    `lambda0`) and returns the value in SI units. `effect` says whether the
    predicted performance accounts for it, so a user changing it knows what the
    figures beside the model do and do not reflect."""
    name: str
    label: str
    unit: str                                  # "m", "deg" or "-"
    default: Callable[[dict], float]
    effect: str = "geometry only: the predicted performance assumes the default"


#: options every antenna has
COMMON = (
    Opt("copper", "conductor sheet thickness", "m", lambda v: SHEET,
        "geometry only: the predictions treat printed conductors as infinitely thin"),
    Opt("margin", "board / ground margin beyond the antenna", "m", lambda v: float("nan"),
        "geometry only: the predictions assume an infinite board and ground where the spec does; "
        "blank keeps each design's own default"),
)


@dataclass
class Options:
    """Values the user set for construction options, by name, in SI units.
    Anything not set takes the option's default."""
    values: dict[str, float] = field(default_factory=dict)

    @property
    def copper(self) -> float:
        return float(self.values.get("copper", SHEET))

    @property
    def margin(self) -> float | None:
        m = self.values.get("margin")
        return None if m is None or not math.isfinite(m) else float(m)


# ---------------------------------------------------------------- the model

@dataclass
class Body:
    name: str
    material: str
    solid: object                      # a manifold3d.Manifold


@dataclass
class Solid3D:
    """A design as solid bodies, each with a material."""
    archetype: str
    title: str
    bodies: list[Body] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def add(self, name: str, material: str, solid) -> None:
        if solid is not None and not solid.is_empty():
            self.bodies.append(Body(name, material, solid))

    @property
    def built(self) -> bool:
        return bool(self.bodies)

    def materials(self) -> dict[str, object]:
        """One solid per material, the conductor taking precedence over anything it overlaps."""
        groups: dict[str, list] = {}
        for b in self.bodies:
            if b.material in _NOT_SOLID:
                continue
            groups.setdefault(b.material, []).append(b.solid)
        out = {k: _union(v) for k, v in groups.items()}
        metal = out.get("PEC")
        if metal is not None:
            for k in list(out):
                # cut only where they truly overlap: subtracting copper that merely touches
                # a board face to face leaves degenerate slivers along the shared plane
                if k != "PEC" and (out[k] ^ metal).volume() > 1e-9 * out[k].volume():
                    out[k] = out[k] - metal
        return {k: v for k, v in out.items() if not v.is_empty()}

    def bounds(self) -> tuple[np.ndarray, np.ndarray]:
        lo, hi = np.full(3, np.inf), np.full(3, -np.inf)
        for b in self.bodies:
            if b.material in _NOT_SOLID:
                continue
            bb = b.solid.bounding_box()
            lo, hi = np.minimum(lo, bb[:3]), np.maximum(hi, bb[3:])
        return lo, hi


def _union(solids: Sequence):
    solids = [s for s in solids if not s.is_empty()]
    if not solids:
        return m3.Manifold()
    if len(solids) == 1:
        return solids[0]
    return m3.Manifold.batch_boolean(solids, m3.OpType.Add)


# ---------------------------------------------------------------- primitives

def _frame(d) -> np.ndarray:
    """A rotation taking +z onto the unit vector d."""
    d = np.asarray(d, float)
    d = d / np.linalg.norm(d)
    z = np.array([0.0, 0.0, 1.0])
    c = float(np.dot(z, d))
    if c > 1 - 1e-12:
        return np.eye(3)
    if c < -1 + 1e-12:
        return np.diag([1.0, -1.0, -1.0])
    v = np.cross(z, d)
    s = np.linalg.norm(v)
    k = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + k + k @ k * ((1 - c) / (s * s))


def place(solid, origin=(0.0, 0.0, 0.0), direction=(0.0, 0.0, 1.0)):
    """Turn a solid built along +z onto `direction`, then move it to `origin`."""
    r = _frame(direction)
    m = np.hstack([r, np.asarray(origin, float).reshape(3, 1)])
    return solid.transform(m)


def box(x, y, z, min_thickness: float = SHEET):
    """An axis-aligned block; a ZERO extent (a sheet) becomes `min_thickness`, centred
    on its plane. A finite extent, however small, is kept: a 70 um slot stays 70 um."""
    _need()
    lo, size = [], []
    for a, b in (x, y, z):
        a, b = min(a, b), max(a, b)
        if b - a <= 1e-12 * max(1.0, abs(a), abs(b)):
            mid = (a + b) / 2
            a, b = mid - min_thickness / 2, mid + min_thickness / 2
        lo.append(a)
        size.append(b - a)
    return m3.Manifold.cube(tuple(size)).translate(tuple(lo))


def rod(p0, p1, r: float, seg: int = SEG):
    """A cylinder of radius r from p0 to p1."""
    return frustum(p0, p1, r, r, seg)


def frustum(p0, p1, r0: float, r1: float, seg: int = SEG):
    _need()
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    length = float(np.linalg.norm(p1 - p0))
    c = m3.Manifold.cylinder(length, max(r0, 0.0), max(r1, 0.0), _seg(seg))
    return place(c, p0, p1 - p0)


def sphere(centre, r: float, seg: int = SEG):
    _need()
    return m3.Manifold.sphere(r, _seg(seg)).translate(tuple(float(v) for v in centre))


def torus(centre, axis, R: float, r: float, seg: int = SEG):
    _need()
    ring = m3.CrossSection.circle(r, max(8, _seg(seg) // 2)).translate((R, 0.0)).revolve(_seg(seg))
    return place(ring, centre, axis)


def plate(points, z0: float = 0.0, thickness: float = SHEET):
    """A flat polygon (x, y points, holes allowed as further loops) extruded in z."""
    _need()
    loops = points if isinstance(points[0][0], (list, tuple, np.ndarray)) else [points]
    cs = m3.CrossSection([np.asarray(p, float) for p in loops], m3.FillRule.EvenOdd)
    return cs.extrude(thickness).translate((0.0, 0.0, z0))


def revolve(profile, seg: int = SEG):
    """A solid of revolution about z from a closed (r, z) outline, r >= 0."""
    _need()
    return m3.CrossSection([np.asarray(profile, float)], m3.FillRule.EvenOdd).revolve(_seg(seg))


def hull(points):
    _need()
    return m3.Manifold.hull_points(np.asarray(points, float))


def tube(path, r: float, seg: int = 16):
    """A wire of radius r swept along a polyline, capped at both ends (a helix, a spiral)."""
    _need()
    pts = np.asarray(path, float)
    seg = max(6, _seg(seg))
    n = len(pts)
    t = np.gradient(pts, axis=0)
    t /= np.linalg.norm(t, axis=1)[:, None]
    # parallel-transported normals, so the rings do not twist
    a = np.cross(t[0], [0.0, 0.0, 1.0])
    if np.linalg.norm(a) < 1e-6:
        a = np.cross(t[0], [1.0, 0.0, 0.0])
    nrm = [a / np.linalg.norm(a)]
    for i in range(1, n):
        v = nrm[-1] - np.dot(nrm[-1], t[i]) * t[i]
        nrm.append(v / np.linalg.norm(v))
    nrm = np.array(nrm)
    bin_ = np.cross(t, nrm)
    ang = 2 * np.pi * np.arange(seg) / seg
    ring = (np.cos(ang)[None, :, None] * nrm[:, None, :] + np.sin(ang)[None, :, None] * bin_[:, None, :])
    verts = (pts[:, None, :] + r * ring).reshape(-1, 3)
    verts = np.vstack([verts, pts[0], pts[-1]])
    tri = []
    for i in range(n - 1):
        for j in range(seg):
            a0, a1 = i * seg + j, i * seg + (j + 1) % seg
            b0, b1 = a0 + seg, a1 + seg
            tri += [(a0, a1, b1), (a0, b1, b0)]
    c0, c1 = n * seg, n * seg + 1
    for j in range(seg):
        tri.append((c0, (j + 1) % seg, j))
        tri.append((c1, (n - 1) * seg + j, (n - 1) * seg + (j + 1) % seg))
    mesh = m3.Mesh(np.ascontiguousarray(verts, dtype=np.float32), np.ascontiguousarray(tri, dtype=np.uint32))
    out = m3.Manifold(mesh)
    if out.is_empty() or out.volume() < 0:
        raise ValueError("a swept wire folded onto itself: its radius is too large for its bends")
    return out


# ---------------------------------------------------------------- from a CST/HFSS model

def _from_model(model, opts: Options | None = None) -> Solid3D:
    """The neutral export model's solids and booleans, as meshes."""
    from .base import Brick, Cone, Cylinder, Sphere, Subtract, Torus, Unite
    opts = opts or Options()
    t_cu = opts.copper
    if opts.margin is not None:
        model = _with_margin(model, opts.margin, t_cu)
    axes = {"x": (1.0, 0.0, 0.0), "y": (0.0, 1.0, 0.0), "z": (0.0, 0.0, 1.0)}

    def point(axis, along, centre):
        a, b = centre
        return {"x": (along, a, b), "y": (a, along, b), "z": (a, b, along)}[axis]

    def turn(solid, deg):
        return solid.rotate((0.0, 0.0, deg)) if deg else solid

    tools = {t for op in model.operations for t in op.tools}
    made: dict[str, object] = {}
    for s in model.solids:
        # a cutting sheet must pass right through the sheet it cuts
        thick = 4 * t_cu if (s.name in tools and s.material == "VOID") else t_cu
        if isinstance(s, Brick):
            made[s.name] = box(s.x, s.y, s.z, thick)
        elif isinstance(s, Cylinder):
            a, b = s.span
            if abs(b - a) <= 1e-12 * max(1.0, abs(a), abs(b)):   # a flat disc: a patch, a ground
                mid = (a + b) / 2
                a, b = mid - thick / 2, mid + thick / 2
            made[s.name] = turn(rod(point(s.axis, a, s.centre), point(s.axis, b, s.centre), s.radius),
                                getattr(s, "rotate_z", 0.0))
        elif isinstance(s, Cone):
            made[s.name] = frustum(point(s.axis, s.span[0], s.centre), point(s.axis, s.span[1], s.centre),
                                   s.radius_start, s.radius_end)
        elif isinstance(s, Sphere):
            made[s.name] = sphere(s.centre, s.radius)
        elif isinstance(s, Torus):
            made[s.name] = torus(s.centre, axes[s.axis], s.major_radius, s.minor_radius)
    gone: set[str] = set()
    for op in model.operations:
        if isinstance(op, Subtract):
            made[op.target] = made[op.target] - _union([made[t] for t in op.tools])
            if not op.keep_tools:
                gone.update(op.tools)
        elif isinstance(op, Unite):
            made[op.target] = _union([made[op.target]] + [made[t] for t in op.tools])
            gone.update(op.tools)
    out = Solid3D(model.archetype, model.title, notes=list(model.notes))
    for s in model.solids:
        if s.name not in gone and s.material not in _NOT_SOLID:
            out.add(s.name, s.material, made[s.name])
    return out


_BOARD_NAMES = ("substrate", "ground", "ground_plane", "gnd")


def _with_margin(model, margin: float, t_cu: float = SHEET):
    """The model with its board and ground plane resized to the antenna's own
    outline plus `margin` on every side. The outline is the plan view of every other
    solid as actually built - rotated wires, horizontal hats and all - not of the
    primitives' nominal extents."""
    from dataclasses import replace
    from .base import Brick, Cylinder
    board = [s for s in model.solids if s.name in _BOARD_NAMES]
    rest = [s for s in model.solids if s.name not in _BOARD_NAMES and s.material not in _NOT_SOLID]
    if not board or not rest:
        return model
    bare = replace(model, solids=[s for s in model.solids if s.name not in _BOARD_NAMES],
                   operations=[op for op in model.operations
                               if op.target not in _BOARD_NAMES and not set(op.tools) & set(_BOARD_NAMES)])
    lo, hi = np.full(2, np.inf), np.full(2, -np.inf)
    for b in _from_model(bare, Options({"copper": t_cu})).bodies:
        bb = np.asarray(b.solid.bounding_box(), float)
        lo, hi = np.minimum(lo, bb[:2]), np.maximum(hi, bb[3:5])
    if not np.all(np.isfinite(lo)):
        return model
    new = []
    for s in model.solids:
        if s in board and isinstance(s, Brick):
            s = replace(s, x=(lo[0] - margin, hi[0] + margin), y=(lo[1] - margin, hi[1] + margin))
        elif s in board and isinstance(s, Cylinder) and s.axis == "z":
            corners = [math.hypot(x - s.centre[0], y - s.centre[1]) for x in (lo[0], hi[0]) for y in (lo[1], hi[1])]
            s = replace(s, radius=max(corners) + margin)
        new.append(s)
    out = replace(model, solids=new)
    kept = [n for n in model.notes if not (n.startswith("Substrate extended") or "ground plane extends" in n)]
    out.notes = kept + [f"Board and ground resized to the antenna's outline plus {margin * 1e3:.4g} mm "
                        "(the margin option); the predictions assume an infinite board and ground."]
    return out


#: archetype key -> builder(design, options) -> Solid3D, for the archetypes the CST/HFSS
#: exporters build no geometry for (filled by mesh_builders)
MESH_BUILDERS: dict[str, Callable[[DesignResult], Solid3D]] = {}


#: archetype key -> its own construction options (beyond COMMON)
MESH_OPTIONS: dict[str, tuple[Opt, ...]] = {}


def mesh_builder(*keys: str, options: Sequence[Opt] = ()):
    def wrap(fn):
        for k in keys:
            MESH_BUILDERS[k] = fn
            MESH_OPTIONS[k] = tuple(options)
        return fn
    return wrap


def construction_options(key: str) -> tuple[Opt, ...]:
    """Every option an archetype's solid model takes: its own, then the common ones."""
    from . import mesh_builders  # noqa: F401
    return tuple(MESH_OPTIONS.get(key, ())) + COMMON


def option_values(design: DesignResult, options: "Options | None" = None) -> dict[str, float]:
    """Each construction option's value for this design: the user's, else the default."""
    v = {**design.requirements, **design.parameters, **design.metrics}
    v.setdefault("lambda0", 2.99792458e8 / float(v.get("f0") or v.get("f_low") or 1e9))
    given = options.values if options else {}
    out = {}
    for o in construction_options(design.archetype):
        out[o.name] = float(given[o.name]) if o.name in given else float(o.default(v))
    return out


def solid(design: DesignResult, options: Options | None = None) -> Solid3D:
    """The design as solid bodies. Empty (with a note saying why) when a dimension
    it needs is unavailable."""
    _need()
    from . import mesh_builders  # noqa: F401  (registers the builders)
    from .base import BUILDERS, build
    opts = options or Options()
    key = design.archetype
    if key in MESH_BUILDERS:
        try:
            out = MESH_BUILDERS[key](design, opts)
        except (KeyError, ValueError) as exc:
            return Solid3D(key, key, notes=[f"NO 3-D GEOMETRY: {exc}."])
        return out if _finite(out) else Solid3D(key, out.title, notes=[
            "NO 3-D GEOMETRY: a dimension this design needs is unavailable (NaN)."])
    if key in BUILDERS:
        model = build(design)
        if not model.built_geometry:
            return Solid3D(key, model.title, notes=list(model.notes))
        out = _from_model(model, opts)
        if opts.copper != SHEET:
            out.notes.append(f"Sheet conductors {opts.copper * 1e6:.4g} um thick.")
        else:
            out.notes.append("Sheet conductors (patches, grounds, strips) 35 um thick, one-ounce copper; set --copper to change.")
        return out
    return Solid3D(key, key, notes=[f"NO 3-D GEOMETRY: no builder for {key!r}."])


def _finite(s: Solid3D) -> bool:
    for b in s.bodies:
        bb = np.asarray(b.solid.bounding_box(), float)
        if not np.all(np.isfinite(bb)):
            return False
    return True


# ---------------------------------------------------------------- STL

def material_label(material: str) -> str:
    if material == "PEC":
        return "conductor"
    if material.startswith("eps_r="):
        return "dielectric_" + material.replace("eps_r=", "er").replace(";tand=", "_tand").replace(";", "_")
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in material)


def _triangles(solid) -> np.ndarray:
    mesh = solid.to_mesh()
    v = np.asarray(mesh.vert_properties, float)[:, :3]
    t = np.asarray(mesh.tri_verts, np.int64)
    return v[t]                                        # (n, 3, 3)


def _stl_bytes(tris: np.ndarray, header: str) -> bytes:
    tris = tris * 1000.0                               # metres -> millimetres
    n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    norm = np.linalg.norm(n, axis=1)
    n = np.where(norm[:, None] > 0, n / np.where(norm > 0, norm, 1)[:, None], 0.0)
    rec = np.zeros(len(tris), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
    rec["n"], rec["v"] = n, tris
    head = header.encode("ascii", "replace")[:80].ljust(80, b" ")
    return head + struct.pack("<I", len(tris)) + rec.tobytes()


def write_stl(model: Solid3D, path) -> list[Path]:
    """Write `path` (every body, for viewing) and, when there is more than one
    material, `<stem>_<material>.stl` per material for import into a solver.
    Returns the files written; nothing when the model has no geometry."""
    path = Path(path)
    if not model.built:
        return []
    groups = model.materials()
    written = []
    allt = np.concatenate([_triangles(s) for s in groups.values()])
    path.write_bytes(_stl_bytes(allt, f"{model.archetype} - all bodies - millimetres"))
    written.append(path)
    if len(groups) > 1:
        for mat, s in groups.items():
            p = path.with_name(f"{path.stem}_{material_label(mat)}.stl")
            p.write_bytes(_stl_bytes(_triangles(s), f"{model.archetype} - {mat} - millimetres"))
            written.append(p)
    return written
