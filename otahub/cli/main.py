"""Command line front end. ``python -m otahub.cli.main --help``"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from typing import Any

from pathlib import Path

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


def cmd_export(args: argparse.Namespace, reg: Registry) -> int:
    from ..export import build
    from ..export import cst as cst_backend
    from ..export import hfss as hfss_backend

    a = reg[args.key]
    reqs = _kv(args.set)
    if args.f0 is not None:
        reqs["f0"] = args.f0
    if not reqs:
        print("give at least --f0 (e.g. --f0 2.4GHz)", file=sys.stderr)
        return 2
    try:
        design = a.synthesize(**reqs)
    except Exception as exc:  # noqa: BLE001
        print(f"synthesis failed: {exc}", file=sys.stderr)
        return 1

    model = build(design)
    backend = {"cst": cst_backend, "hfss": hfss_backend}[args.format]
    text = backend.render(model)

    if args.output:
        path = Path(args.output)
        path.write_text(text)
        print(f"wrote {path}  ({len(text.splitlines())} lines)")
    else:
        print(text)

    if not model.built_geometry:
        print(f"\nNOTE: no solid geometry was generated for {a.key!r}; the file "
              f"defines parameters only.", file=sys.stderr)
    return 0


def cmd_match(args: argparse.Namespace, reg: Registry) -> int:
    from ..utils import matching as M
    from ..utils.network import match_quality

    z_load = complex(args.r, args.x)
    before = match_quality(z_load, args.z0)
    print(f"load {z_load.real:.4g}{z_load.imag:+.4g}j ohm into {args.z0:.4g} ohm")
    print(f"  unmatched: {before}")

    sections = M.l_section(z_load, args.z0, args.f0 or 1e9)
    if sections:
        print(f"\nL-section solutions at {engineering(args.f0 or 1e9, 'Hz')}:")
        for i, s in enumerate(sections, 1):
            print(f"  [{i}] {s.topology}")
            print(f"      series: {s.series}")
            print(f"      shunt : {s.shunt}")
            print(f"      achieved Zin = {s.achieved.real:.6g}{s.achieved.imag:+.6g}j ohm")
    else:
        print("\nno L-section solution exists for this load")

    for kind in ("short", "open"):
        stubs = M.single_stub(z_load, args.z0, kind)
        if stubs:
            print(f"\nshunt {kind}-circuit stub tuners:")
            for s in stubs:
                print(f"  {s}")

    try:
        z1 = M.quarter_wave(z_load, args.z0)
        bw = M.quarter_wave_bandwidth(z_load.real, args.z0)
        print(f"\nquarter-wave transformer: Z1 = {z1:.4g} ohm, "
              f"fractional bandwidth (VSWR<2) = {bw:.3f}")
    except ValueError as exc:
        print(f"\nquarter-wave transformer: not applicable - {exc}")
    return 0


def cmd_gui(args: argparse.Namespace, reg: Registry) -> int:
    try:
        from ..gui.app import main as gui_main
    except ImportError as exc:
        print(f"the GUI needs PySide6: {exc}\n\n"
              f"  python -m pip install PySide6", file=sys.stderr)
        return 1
    return gui_main([sys.argv[0]])


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



# ------------------------------------------------------- waveguides & arrays

def cmd_guide(args: argparse.Namespace, reg: Registry) -> int:
    from ..waveguides.rectangular import WR_SERIES, recommended_band, standard

    if args.name.lower() in ("list", "all"):
        print(f"{'guide':<10}{'a (mm)':>10}{'b (mm)':>10}{'TE10 (GHz)':>12}"
              f"{'band (GHz)':>18}")
        for name in WR_SERIES:
            g = standard(name)
            lo, hi = recommended_band(name)
            print(f"{name:<10}{g.a*1e3:>10.3f}{g.b*1e3:>10.3f}"
                  f"{g.dominant_cutoff_hz/1e9:>12.4f}"
                  f"{lo/1e9:>9.2f}-{hi/1e9:<8.2f}")
        return 0

    try:
        g = standard(args.name)
    except KeyError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    lo, hi = g.single_mode_band_hz
    rlo, rhi = recommended_band(args.name)
    print(f"{args.name.upper()}  {g.a*1e3:.3f} x {g.b*1e3:.3f} mm  "
          f"(aspect {g.aspect_ratio:.3f})")
    print(f"TE10 cutoff     : {engineering(g.dominant_cutoff_hz, 'Hz')}")
    print(f"single-mode band: {lo/1e9:.4f} - {hi/1e9:.4f} GHz")
    print(f"recommended band: {rlo/1e9:.3f} - {rhi/1e9:.3f} GHz")

    if args.f0 is not None:
        f0 = args.f0
        print(f"\nat {engineering(f0, 'Hz')}:")
        if f0 <= g.dominant_cutoff_hz:
            print("  BELOW CUTOFF - the guide is evanescent here, nothing propagates.")
            return 0
        if not g.is_single_mode(f0):
            print("  ! WARNING: above the single-mode band, higher modes propagate")
        print(f"  guide wavelength  {engineering(g.guide_wavelength(f0), 'm')}")
        print(f"  wave impedance    {g.wave_impedance(f0):.2f} ohm")
        print(f"  attenuation       {g.attenuation_db_per_m(f0):.4f} dB/m")
        print(f"  power capacity    {g.power_capacity_w(f0)/1e6:.3f} MW (dry air, derate hard)")
    print("\nmodes:")
    for mode in g.modes(args.f_max if args.f_max else 3 * g.dominant_cutoff_hz):
        print(f"  {mode}")
    return 0


def cmd_array(args: argparse.Namespace, reg: Registry) -> int:
    from ..arrays import TAPERS, summarise

    taper = args.taper.lower()
    if taper not in TAPERS:
        print(f"unknown taper {taper!r}; choose from {', '.join(TAPERS)}", file=sys.stderr)
        return 1
    # "30 dB sidelobes" and "-30 dB sidelobes" mean the same thing to everyone
    # who says it out loud, and a sidelobe ABOVE the main beam is not a thing,
    # so accept either sign rather than raising on the commoner spelling.
    sll = -abs(args.sll)
    if taper in ("chebyshev", "taylor"):
        try:
            weights = TAPERS[taper](args.n, sll)
        except ValueError as exc:
            print(f"cannot synthesise that taper: {exc}", file=sys.stderr)
            return 1
    elif taper == "cosine":
        weights = TAPERS[taper](args.n, args.pedestal)
    else:
        weights = TAPERS[taper](args.n)

    s = summarise(weights, args.d, args.scan)
    design = f", {sll:.1f} dB design" if taper in ("chebyshev", "taylor") else ""
    print(f"{args.n}-element linear array, {taper} taper{design}, d = {args.d:.3f} lambda, "
          f"scan {args.scan:.1f} deg")
    print(f"\n  directivity        {s['directivity_dbi']:.2f} dBi")
    print(f"  half-power beam    {s['hpbw_deg']:.3f} deg")
    print(f"  first sidelobe     {s['sidelobe_db']:.2f} dB")
    print(f"  taper efficiency   {s['taper_efficiency']:.4f}"
          f"  ({-10*math.log10(s['taper_efficiency']):.2f} dB of directivity given up)")
    print(f"  max spacing        {s['max_spacing_no_grating']:.4f} lambda "
          f"before a grating lobe at this scan angle")
    if s["grating_lobe"]:
        print("  ! GRATING LOBE present at this spacing and scan angle")
    print("\nexcitation (normalised amplitude):")
    for i, w in enumerate(weights):
        bar = "#" * int(round(w * 40))
        print(f"  {i:>3} {w:7.4f} {bar}")
    return 0


def cmd_planar(args: argparse.Namespace, reg: Registry) -> int:
    from ..arrays import (TAPERS, dolph_chebyshev, grating_lobe_free_spacing_planar,
                          lattice_element_saving, planar_summarise,
                          rectangular_lattice, separable_weights, taylor_nbar,
                          triangular_lattice, uniform)

    lattice = args.lattice.strip().lower()
    taper = args.taper.lower()
    if taper not in TAPERS:
        print(f"unknown taper {taper!r}; choose from {', '.join(TAPERS)}", file=sys.stderr)
        return 1
    sll = -abs(args.sll)

    def line_taper(n):
        if taper in ("chebyshev", "taylor"):
            return TAPERS[taper](n, sll)
        if taper == "cosine":
            return TAPERS[taper](n, args.pedestal)
        return TAPERS[taper](n)

    if lattice in ("rect", "rectangular", "square"):
        positions = rectangular_lattice(args.nx, args.ny, args.d, args.dy or args.d)
        try:
            weights = separable_weights(line_taper(args.nx), line_taper(args.ny))
        except ValueError as exc:
            print(f"cannot synthesise that taper: {exc}", file=sys.stderr)
            return 1
    elif lattice in ("tri", "triangular", "hex", "hexagonal"):
        positions = triangular_lattice(args.nx, args.ny, args.d)
        if taper != "uniform":
            print("note: a triangular lattice is not separable, so the taper is "
                  "ignored and the array is excited uniformly", file=sys.stderr)
        weights = [1.0] * len(positions)
    else:
        print(f"unknown lattice {args.lattice!r}; use rectangular or triangular",
              file=sys.stderr)
        return 1

    from ..arrays import elements as elem
    element = None
    spec = (args.element or "isotropic").strip().lower()
    kind, _, val = spec.partition(":")
    try:
        if kind == "cos":
            element = elem.cosine(float(val) if val else 1.0)
        elif kind in ("dipole-x", "dipole-y"):
            element = elem.short_dipole(kind[-1], float(val) if val else None)
        elif kind != "isotropic":
            raise ValueError(f"unknown element {args.element!r}")
    except ValueError as exc:
        print(f"{exc}; use isotropic, cos[:q], dipole-x[:height] or dipole-y[:height]", file=sys.stderr)
        return 1
    if element is not None and args.ground_plane:
        print("note: the element pattern sets the hemisphere itself; --ground-plane is ignored", file=sys.stderr)
    s = planar_summarise(positions, weights, args.scan, args.scan_phi,
                         half_space=args.ground_plane and element is None, element=element)
    label = "triangular" if lattice.startswith(("tri", "hex")) else "rectangular"
    print(f"{s['elements']}-element {label} planar array, {taper} taper, "
          f"d = {args.d:.3f} lambda")
    print(f"  aperture           {s['aperture_x_lambda']:.3f} x "
          f"{s['aperture_y_lambda']:.3f} lambda")
    print(f"  beam               theta {s['scan_theta_deg']:.1f} deg from the "
          f"normal, phi {s['scan_phi_deg']:.1f} deg")
    if element is not None:
        print(f"\n  element            {s['element']}, {s['element_directivity_dbi']:.2f} dBi on its own "
              f"toward the beam")
        print(f"  directivity        {s['directivity_dbi']:.2f} dBi toward the beam, element pattern included")
    else:
        print(f"\n  directivity        {s['directivity_dbi']:.2f} dBi"
              + ("  (ground-plane backed)" if args.ground_plane else
                 "  (isotropic elements, so half the power goes into the mirror beam;"
                 " pass --ground-plane for the one-sided figure)"))
    print(f"  half-power beam    {s['hpbw_scan_plane_deg']:.3f} deg in the scan plane, "
          f"{s['hpbw_cross_plane_deg']:.3f} deg across it")
    print(f"  first sidelobe     {s['sidelobe_scan_plane_db']:.2f} dB in the scan plane, "
          f"{s['sidelobe_cross_plane_db']:.2f} dB across it")
    print(f"  taper efficiency   {s['taper_efficiency']:.4f}")

    limit = grating_lobe_free_spacing_planar(abs(args.scan), label)
    print(f"\n  grating-lobe limit {limit:.4f} lambda for scanning to "
          f"{abs(args.scan):.1f} deg")
    if args.d >= limit:
        print(f"  ! SPACING {args.d:.4f} EXCEEDS THAT LIMIT - a grating lobe is in real space")
    other = "triangular" if label == "rectangular" else "rectangular"
    other_limit = grating_lobe_free_spacing_planar(abs(args.scan), other)
    print(f"  a {other} lattice would allow {other_limit:.4f} lambda")
    if label == "rectangular":
        print(f"  and cover the same aperture with "
              f"{lattice_element_saving()*100:.2f}% fewer elements")
    return 0


def cmd_line(args: argparse.Namespace, reg: Registry) -> int:
    from ..waveguides import lines as L

    kind = args.kind.lower()
    if kind == "microstrip":
        if args.z0 is not None:
            w = L.microstrip_width_for(args.z0, args.h, args.eps_r)
            print(f"microstrip on eps_r={args.eps_r}, h={engineering(args.h,'m')}")
            print(f"  target Z0     {args.z0:.2f} ohm")
            print(f"  strip width   {engineering(w, 'm')}  (w/h = {w/args.h:.4f})")
        elif args.w is not None:
            w = args.w
            print(f"  Z0            {L.microstrip_impedance(w, args.h, args.eps_r):.3f} ohm")
        else:
            print("give either --z0 (synthesis) or --w (analysis)", file=sys.stderr)
            return 2
        ee = L.microstrip_eps_eff(w, args.h, args.eps_r)
        print(f"  eps_eff       {ee:.4f}")
        if args.f0:
            lam = L.microstrip_guide_wavelength(args.f0, w, args.h, args.eps_r)
            print(f"  guide wavelen {engineering(lam, 'm')} at {engineering(args.f0,'Hz')}")
            print(f"  quarter wave  {engineering(lam/4, 'm')}")
        return 0
    if kind == "coax":
        if args.z0 is None:
            print("coax needs --z0", file=sys.stderr)
            return 2
        ratio = math.exp(args.z0 * 2 * math.pi * math.sqrt(args.eps_r) / 376.730313412)
        print(f"coax for {args.z0:.2f} ohm in eps_r={args.eps_r}: b/a = {ratio:.4f}")
        opt = L.coax_optimum_ratios()
        print(f"  minimum attenuation would be b/a = {opt['min_attenuation_ratio']:.4f} "
              f"({opt['min_attenuation_z0_air']:.1f} ohm in air)")
        print(f"  maximum power     would be b/a = {opt['max_power_ratio']:.4f} "
              f"({opt['max_power_z0_air']:.1f} ohm in air)")
        print("  50 ohm is the historical compromise between those, not an optimum")
        return 0
    print(f"unknown line kind {kind!r}; choose microstrip or coax", file=sys.stderr)
    return 1


# --------------------------------------------------------------------- entry

def cmd_touchstone(args: argparse.Namespace, reg: Registry) -> int:
    from ..utils.touchstone import compare_to_prediction, read_touchstone

    try:
        net = read_touchstone(Path(args.path), n_ports=args.ports)
    except (OSError, ValueError) as exc:
        print(f"could not read {args.path}: {exc}", file=sys.stderr)
        return 1

    lo, hi = net.frequency_hz[0], net.frequency_hz[-1]
    print(f"{args.path}: {net.n_ports}-port, {len(net.frequency_hz)} points, "
          f"{lo/1e9:.4g} to {hi/1e9:.4g} GHz, reference {net.z0:g} ohm")
    for note in net.comments[:4]:
        print(f"  ! {note}")

    port = args.port
    db = net.s_db(port)
    best = int(db.argmin())
    print(f"\n  best match at port {port + 1}: {db[best]:.2f} dB at "
          f"{net.frequency_hz[best]/1e9:.6g} GHz")
    print(f"  impedance there      {net.impedance_at_port(port)[best]:.4g} ohm")
    print(f"  VSWR there           {net.vswr(port)[best]:.4f}")
    dips = net.resonances(port, args.threshold)
    if dips:
        print(f"  resonances below {args.threshold:.0f} dB: "
              + ", ".join(f"{d/1e9:.6g} GHz" for d in dips))
    else:
        print(f"  no resonance dips below {args.threshold:.0f} dB")

    if args.compare:
        key = args.compare
        if key not in reg:
            print(f"unknown archetype {key!r}", file=sys.stderr)
            return 1
        overrides = dict(pair.split("=", 1) for pair in (args.set or []))
        given = {k: parse_quantity(v) for k, v in overrides.items()}
        f_cmp = args.at if args.at else net.frequency_hz[best]
        given.setdefault("f0", f_cmp)
        try:
            design = reg[key].synthesize(**given)
        except Exception as exc:  # noqa: BLE001
            print(f"could not synthesise {key}: {exc}", file=sys.stderr)
            return 1
        r = design.get("input_resistance_ohm") or design.get("radiation_resistance_ohm")
        x = design.get("input_reactance_ohm") or 0.0
        if r is None:
            print(f"{key} predicts no input resistance to compare against",
                  file=sys.stderr)
            return 1
        out = compare_to_prediction(net, complex(float(r), float(x)), f_cmp, port)
        print(f"\n  against {key} at {f_cmp/1e9:.6g} GHz:")
        print(f"    predicted  {out['predicted_impedance_ohm']:.4g} ohm, "
              f"S11 {out['predicted_s11_db']:.2f} dB, VSWR {out['predicted_vswr']:.3f}")
        print(f"    measured   {out['measured_impedance_ohm']:.4g} ohm, "
              f"S11 {out['measured_s11_db']:.2f} dB, VSWR {out['measured_vswr']:.3f}")
        print(f"    resistance off by {out['resistance_error_pct']:+.2f}%, "
              f"reactance by {out['reactance_error_ohm']:+.3g} ohm")
    return 0


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

    pe = sub.add_parser("export", help="export a design to CST or HFSS")
    pe.add_argument("key")
    pe.add_argument("--format", choices=["cst", "hfss"], default="cst")
    pe.add_argument("--f0", type=parse_quantity, help="design frequency")
    pe.add_argument("--set", action="append", metavar="NAME=VALUE")
    pe.add_argument("-o", "--output", help="write to a file instead of stdout")
    pe.set_defaults(func=cmd_export)

    pm = sub.add_parser("match", help="impedance matching networks for a load")
    pm.add_argument("--r", type=float, required=True, help="load resistance [ohm]")
    pm.add_argument("--x", type=float, default=0.0, help="load reactance [ohm]")
    pm.add_argument("--z0", type=float, default=50.0, help="reference impedance")
    pm.add_argument("--f0", type=parse_quantity, help="frequency for component values")
    pm.set_defaults(func=cmd_match)

    pui = sub.add_parser("gui", help="launch the graphical interface")
    pui.set_defaults(func=cmd_gui)

    pg = sub.add_parser("guide", help="analyse a rectangular waveguide")
    pg.add_argument("name", help="WR designation, or 'list' for the whole series")
    pg.add_argument("--f0", type=parse_quantity, help="operating frequency, e.g. 10GHz")
    pg.add_argument("--f-max", type=parse_quantity, help="highest cutoff to list")
    pg.set_defaults(func=cmd_guide)

    pa = sub.add_parser("array", help="synthesise a linear array taper")
    pa.add_argument("-n", type=int, default=16, help="element count")
    pa.add_argument("--taper", default="chebyshev",
                    help="uniform | binomial | chebyshev | taylor | cosine")
    pa.add_argument("--sll", type=float, default=-30.0,
                    help="design sidelobe level in dB below the main beam; "
                         "either sign accepted, so 30 and -30 both mean 30 dB down")
    pa.add_argument("--d", type=float, default=0.5, help="spacing in wavelengths")
    pa.add_argument("--scan", type=float, default=90.0,
                    help="beam direction from the array axis, 90 = broadside")
    pa.add_argument("--pedestal", type=float, default=0.0, help="cosine taper pedestal")
    pa.set_defaults(func=cmd_array)

    pp = sub.add_parser("planar", help="synthesise a planar array")
    pp.add_argument("--nx", type=int, default=8, help="elements along x")
    pp.add_argument("--ny", type=int, default=8, help="elements along y (rows)")
    pp.add_argument("--d", type=float, default=0.5,
                    help="element spacing in wavelengths; for a triangular "
                         "lattice this is the nearest-neighbour spacing")
    pp.add_argument("--dy", type=float, default=None,
                    help="y spacing if it differs from --d (rectangular only)")
    pp.add_argument("--lattice", default="rectangular",
                    help="rectangular | triangular")
    pp.add_argument("--taper", default="uniform",
                    help="uniform | binomial | chebyshev | taylor | cosine, "
                         "applied separably along each axis")
    pp.add_argument("--sll", type=float, default=-30.0,
                    help="design sidelobe level in dB below the main beam; "
                         "either sign accepted")
    pp.add_argument("--pedestal", type=float, default=0.0, help="cosine taper pedestal")
    pp.add_argument("--scan", type=float, default=0.0,
                    help="beam angle from the array NORMAL in degrees; 0 is broadside")
    pp.add_argument("--scan-phi", dest="scan_phi", type=float, default=0.0,
                    help="azimuth of the scan direction in degrees")
    pp.add_argument("--ground-plane", dest="ground_plane", action="store_true",
                    help="report the one-sided directivity of a backed array")
    pp.add_argument("--element", default="isotropic",
                    help="element pattern: isotropic | cos[:q] (cos^q power, ground-backed; "
                         "q = 1 is the ideal matched-array element) | dipole-x[:height] | "
                         "dipole-y[:height] (a short dipole, optionally that many wavelengths "
                         "over ground)")
    pp.set_defaults(func=cmd_planar)

    pt = sub.add_parser("touchstone",
                        help="read measured or simulated S-parameters")
    pt.add_argument("path", help="a .sNp file")
    pt.add_argument("--ports", type=int, default=None,
                    help="port count, if the filename does not say")
    pt.add_argument("--port", type=int, default=0,
                    help="which port to report, 0-based")
    pt.add_argument("--threshold", type=float, default=-10.0,
                    help="return-loss threshold for calling something a resonance")
    pt.add_argument("--compare", default=None,
                    help="archetype key to compare the measurement against")
    pt.add_argument("--at", type=parse_quantity, default=None,
                    help="frequency for the comparison; defaults to the best match")
    pt.add_argument("--set", action="append", metavar="NAME=VALUE",
                    help="requirement for the compared archetype")
    pt.set_defaults(func=cmd_touchstone)

    pn = sub.add_parser("line", help="synthesise a transmission line")
    pn.add_argument("kind", help="microstrip | coax")
    pn.add_argument("--z0", type=float, help="target impedance [ohm]")
    pn.add_argument("--w", type=parse_quantity, help="strip width, for analysis")
    pn.add_argument("--h", type=parse_quantity, default=1.6e-3, help="substrate thickness")
    pn.add_argument("--eps-r", dest="eps_r", type=float, default=4.4)
    pn.add_argument("--f0", type=parse_quantity, help="frequency for guide wavelength")
    pn.set_defaults(func=cmd_line)
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
