"""Geometry checks on exported models: where a port's ends lie, and which solids
intersect. Axis-aligned primitives only, plus cylinders turned about z - which is
every shape the builders make."""
from __future__ import annotations

import math

from otahub.export.base import Brick, Cone, Cylinder, Sphere, Subtract, Torus, Unite


def _distance(pt, s) -> float:
    """From a point to a solid, 0 inside."""
    x, y, z = pt
    if isinstance(s, Brick):
        d = [max(lo - v, 0.0, v - hi) for v, (lo, hi) in zip(pt, (s.x, s.y, s.z))]
        return math.sqrt(sum(c * c for c in d))
    if isinstance(s, (Cylinder, Cone)):
        if getattr(s, "rotate_z", 0.0):
            a = math.radians(-s.rotate_z)
            x, y = x * math.cos(a) - y * math.sin(a), x * math.sin(a) + y * math.cos(a)
        along = {"x": x, "y": y, "z": z}[s.axis]
        u, v = {"x": (y, z), "y": (x, z), "z": (x, y)}[s.axis]
        lo, hi = min(s.span), max(s.span)
        r = s.radius if isinstance(s, Cylinder) else max(s.radius_start, s.radius_end)
        rho = math.hypot(u - s.centre[0], v - s.centre[1])
        return math.hypot(max(rho - r, 0.0), along - min(max(along, lo), hi))
    if isinstance(s, Torus) and s.axis == "z":
        cx, cy, cz = s.centre
        return max(math.hypot(math.hypot(x - cx, y - cy) - s.major_radius, z - cz) - s.minor_radius, 0.0)
    if isinstance(s, Sphere):
        return max(math.dist(pt, s.centre) - s.radius, 0.0)
    return math.inf


def port_problems(model) -> list:
    """Ports with an end touching no conductor."""
    removed = {t for op in model.operations if isinstance(op, Subtract) and not op.keep_tools for t in op.tools}
    pec = [s for s in model.solids if s.material == "PEC" and s.name not in removed]
    out = []
    for port in model.ports:
        tol = max(1e-6, 1e-3 * math.dist(port.start, port.end))
        for end in (port.start, port.end):
            if not any(_distance(end, s) <= tol for s in pec):
                out.append((port.name, end))
    return out


def _box(s):
    if isinstance(s, Brick):
        return (*s.x, *s.y, *s.z)
    if isinstance(s, (Cylinder, Cone)):
        r = s.radius if isinstance(s, Cylinder) else max(s.radius_start, s.radius_end)
        lo, hi = min(s.span), max(s.span)
        c0, c1 = s.centre
        if getattr(s, "rotate_z", 0.0) and s.axis == "x":
            a = math.radians(s.rotate_z)
            xs, ys = [lo * math.cos(a), hi * math.cos(a)], [lo * math.sin(a), hi * math.sin(a)]
            return (min(xs) - r, max(xs) + r, min(ys) - r, max(ys) + r, c1 - r, c1 + r)
        return {"x": (lo, hi, c0 - r, c0 + r, c1 - r, c1 + r), "y": (c0 - r, c0 + r, lo, hi, c1 - r, c1 + r),
                "z": (c0 - r, c0 + r, c1 - r, c1 + r, lo, hi)}[s.axis]
    if isinstance(s, Torus):
        R, (cx, cy, cz) = s.major_radius + s.minor_radius, s.centre
        return (cx - R, cx + R, cy - R, cy + R, cz - s.minor_radius, cz + s.minor_radius)
    if isinstance(s, Sphere):
        (cx, cy, cz), R = s.centre, s.radius
        return (cx - R, cx + R, cy - R, cy + R, cz - R, cz + R)
    return None


def overlaps(model) -> list:
    """Pairs of solids whose volumes intersect and that no boolean resolves - by
    bounding box, so a hollow shell can look intersecting when it is not."""
    removed = {t for op in model.operations if isinstance(op, Subtract) and not op.keep_tools for t in op.tools}
    resolved = [frozenset((op.target, *op.tools)) for op in model.operations
                if isinstance(op, Unite) or (isinstance(op, Subtract) and op.keep_tools)]
    solids = [s for s in model.solids if s.name not in removed]
    out = []
    for i, a in enumerate(solids):
        for b in solids[i + 1:]:
            ba, bb = _box(a), _box(b)
            if ba is None or bb is None:
                continue
            scale = max(max(ba[1::2]) - min(ba[0::2]), 1e-12)
            gaps = [min(ba[2 * k + 1], bb[2 * k + 1]) - max(ba[2 * k], bb[2 * k]) for k in range(3)]
            if all(g > 1e-6 * scale for g in gaps) and not any({a.name, b.name} <= r for r in resolved):
                out.append((a.name, b.name))
    return out
