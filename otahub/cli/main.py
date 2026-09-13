"""Command line front end. ``python -m otahub.cli.main --help``"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from typing import Any

from ..core.registry import Registry, default_registry
from ..core.units import engineering, to_si

_QTY = re.compile(r"^\s*([+-]?[0-9.]+(?:[eE][+-]?[0-9]+)?)\s*([A-Za-z]*)\s*$")
_SUFFIX_SI = {
    "GHz": 1e9, "MHz": 1e6, "kHz": 1e3, "Hz": 1.0,
    "mm": 1e-3, "cm": 1e-2, "m": 1.0, "um": 1e-6, "mil": 2.54e-5, "in": 0.0254,
    "ohm": 1.0, "deg": math.pi / 180.0, "rad": 1.0,
}


def parse_quantity(text: str) -> float:
    """Parse ``2.4GHz`` / ``1.6mm`` / ``3.5`` into an SI float."""
    m = _QTY.match(text)
    if not m:
        raise argparse.ArgumentTypeError(f"cannot parse quantity {text!r}")
    value, suffix = float(m.group(1)), m.group(2)
    if not suffix:
        return value
    for unit, factor in _SUFFIX_SI.items():          # exact, case-sensitive first
        if suffix == unit:
            return value * factor
    for unit, factor in _SUFFIX_SI.items():          # then case-insensitive
        if suffix.lower() == unit.lower():
            return value * factor
    raise argparse.ArgumentTypeError(f"unknown unit {suffix!r} in {text!r}")


def _kv(pairs: list[str]) -> dict[str, float]:
    out: dict[str, float] = {}
    for item in pairs or []:
        if "=" not in item:
            raise argparse.ArgumentTypeError(f"expected name=value, got {item!r}")
        name, _, raw = item.partition("=")
        out[name.strip()] = parse_quantity(raw)
    return out


def _looks_numeric(text: str) -> bool:
    """True when a spec's `typical` is a single number we can paste into a flag."""
    try:
        float(text)
    except (TypeError, ValueError):
        return False
    return True


def _fmt(value: Any, unit: str = "") -> str:
    if isinstance(value, complex):
        sign = "+" if value.imag >= 0 else "-"
        return f"{value.real:.4g} {sign} j{abs(value.imag):.4g} {unit}".strip()
    if isinstance(value, float):
        if unit in ("m", "Hz") and value != 0 and math.isfinite(value):
            return engineering(value, unit)
        return f"{value:.6g} {unit}".strip()
    return f"{value} {unit}".strip()


# --------------------------------------------------------------- subcommands

def cmd_list(args: argparse.Namespace, reg: Registry) -> int:
    items = list(reg)
    if args.family:
        items = reg.by_family(args.family)
    if args.search:
        items = [a for a in reg.search(args.search) if a in items]
    if args.covering is not None:
        cover = set(a.key for a in reg.covering(args.covering))
        items = [a for a in items if a.key in cover]
    if not items:
        print("no archetypes match.", file=sys.stderr)
        return 1
    width = max(len(a.key) for a in items)
    current = None
    for a in sorted(items, key=lambda x: (x.family, x.key)):
        if a.family != current:
            current = a.family
            print(f"\n{current.upper()}")
        flag = "  [low confidence]" if a.spec.confidence == "low" else ""
        print(f"  {a.key:<{width}}  {a.name}{flag}")
    print(f"\n{len(items)} archetype(s), {len(reg.families)} family(ies).")
    return 0


def cmd_show(args: argparse.Namespace, reg: Registry) -> int:
    a = reg[args.key]
    s = a.spec
    lo, hi = s.freq_range_hz
    print(f"{s.name}  [{s.key}]")
    print(f"family      : {s.family}")
    if s.summary:
        print(f"summary     : {s.summary}")
    print(f"valid band  : {engineering(lo,'Hz')} .. {engineering(hi,'Hz')}")
    if s.confidence == "low":
        print("confidence  : LOW — verify in a full-wave solver")
    if s.parameters:
        print("\nparameters:")
        for p in s.parameters:
            typ = f"  (typical {p.typical})" if p.typical else ""
            print(f"  {p.symbol:<12} {p.unit:<6} {p.role:<12} {p.description}{typ}")
    if s.synthesis:
        print("\nsynthesis:")
        for r in s.synthesis:
            print(f"  {r.output:<12} = {r.expr}")
            if r.notes:
                print(f"               ^ {r.notes}")
    if s.analysis:
        print("\nanalysis:")
        for r in s.analysis:
            print(f"  {r.metric:<24} = {r.expr}   [{r.units_out}]")
    if s.validity:
        print("\nvalidity:")
        for v in s.validity:
            print(f"  - {v}")
    if s.references:
        print("\nreferences:")
        for r in s.references:
            print(f"  - {r}")
    probs = s.problems()
    if probs:
        print("\nPROBLEMS:")
        for p in probs:
            print(f"  ! {p}")
    return 0


def cmd_synth(args: argparse.Namespace, reg: Registry) -> int:
    a = reg[args.key]
    reqs = _kv(args.set)
    if args.f0 is not None:
        reqs["f0"] = args.f0
    if not reqs:
        print("give at least --f0 (e.g. --f0 2.4GHz)", file=sys.stderr)
        return 2
    try:
        design = a.synthesize(**reqs)
    except Exception as exc:  # noqa: BLE001 - CLI boundary
        print(f"synthesis failed: {exc}", file=sys.stderr)
        return 1
    if args.json:
        def _enc(o: Any) -> Any:
            if isinstance(o, complex):
                return {"real": o.real, "imag": o.imag}
            if hasattr(o, "tolist"):
                return o.tolist()
            return str(o)
        print(json.dumps(design.as_dict(), indent=2, default=_enc))
        return 0
    print(f"{a.name}  [{a.key}]")
    print(f"requirements: " + ", ".join(f"{k}={_fmt(v)}" for k, v in reqs.items()))
    print("\ngeometry:")
    for k, v in design.parameters.items():
        if k in ("lambda0", "k0"):
            continue
        print(f"  {k:<22} {_fmt(v, design.units.get(k,''))}")
    if design.metrics:
        print("\npredicted performance:")
        for k, v in design.metrics.items():
            print(f"  {k:<22} {_fmt(v, design.units.get(k,''))}")
    if design.warnings:
        print("\nwarnings:")
        for w in design.warnings:
            print(f"  ! {w}")

    missing = design.missing_requirements()
    if missing:
        typicals = {p.symbol: p.typical for p in a.spec.parameters}
        print("\nto complete this design, supply:")
        for sym in missing:
            hint = typicals.get(sym, "")
            print(f"  {sym:<20} {('typical ' + hint) if hint else ''}")
        suggestion = " ".join(
            f"--set {s}={typicals[s]}" for s in missing
            if _looks_numeric(typicals.get(s, ""))
        )
        if suggestion:
            freq = f"--f0 {reqs['f0']:.6g}" if "f0" in reqs else ""
            print(f"\n  python OTA_Hub_AntennaToolkit.py synth {a.key} {freq} {suggestion}".rstrip())
    if design.notes:
        print("\nvalid only under:")
        for n in design.notes:
            print(f"  - {n}")
    return 0


def cmd_check(args: argparse.Namespace, reg: Registry) -> int:
    rows: list[dict[str, Any]] = []
    for a in reg:
        if args.key and a.key != args.key:
            continue
        rows.extend(a.check_known_cases())
    if not rows:
        print("no known cases to check.", file=sys.stderr)
        return 1
    failed = [r for r in rows if not r["passed"]]
    for r in rows:
        if r["passed"] and not args.verbose:
            continue
        mark = "PASS" if r["passed"] else "FAIL"
        err = f"{r['error_pct']:.2f}%" if r["error_pct"] is not None else "n/a"
        print(f"[{mark}] {r['archetype']}.{r['quantity']}: "
              f"expected {r['expected']}, got {r['actual']} (err {err}, "
              f"tol {r['tol_pct']}%) {r['detail']}")
    print(f"\n{len(rows) - len(failed)}/{len(rows)} known cases pass.")
    return 1 if failed else 0


def cmd_doctor(args: argparse.Namespace, reg: Registry) -> int:
    probs = reg.problems()
    if not probs:
        print(f"all {len(reg)} archetype(s) structurally sound.")
        return 0
    for key, items in sorted(probs.items()):
        print(f"\n{key}:")
        for p in items:
            print(f"  ! {p}")
    total = sum(len(v) for v in probs.values())
    print(f"\n{total} problem(s) across {len(probs)} archetype(s)/file(s).")
    return 1


# --------------------------------------------------------------------- entry

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="otahub", description="OTA Hub Antenna Toolkit — synthesise and analyse antennas.")
    p.add_argument("--spec-dir", default=None, help="override the spec directory")
    sub = p.add_subparsers(dest="command", required=True)

    pl = sub.add_parser("list", help="list archetypes")
    pl.add_argument("--family"); pl.add_argument("--search")
    pl.add_argument("--covering", type=parse_quantity, metavar="FREQ",
                    help="only those valid at this frequency, e.g. 2.4GHz")
    pl.set_defaults(func=cmd_list)

    ps = sub.add_parser("show", help="show one archetype in full")
    ps.add_argument("key"); ps.set_defaults(func=cmd_show)

    py = sub.add_parser("synth", help="synthesise a design")
    py.add_argument("key")
    py.add_argument("--f0", type=parse_quantity, help="design frequency, e.g. 2.4GHz")
    py.add_argument("--set", action="append", metavar="NAME=VALUE",
                    help="extra requirement, repeatable (e.g. --set eps_r=4.4)")
    py.add_argument("--json", action="store_true")
    py.set_defaults(func=cmd_synth)

    pc = sub.add_parser("check", help="run the citable known cases")
    pc.add_argument("--key"); pc.add_argument("-v", "--verbose", action="store_true")
    pc.set_defaults(func=cmd_check)

    pd = sub.add_parser("doctor", help="report structural faults in the specs")
    pd.set_defaults(func=cmd_doctor)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    reg = Registry.load(args.spec_dir) if args.spec_dir else default_registry()
    if not len(reg) and not reg.load_errors:
        print("no specs found — is the spec directory populated?", file=sys.stderr)
        return 1
    return args.func(args, reg)


if __name__ == "__main__":
    raise SystemExit(main())
