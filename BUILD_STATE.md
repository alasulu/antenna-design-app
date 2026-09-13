# OTA Hub Antenna Toolkit — Build State

**Single source of truth for scheduled resume runs.** Read this first, update it last.

- Project root: `/Users/macbookair/Desktop/Antenna Design App`
- Python: `/Users/macbookair/miniconda3/bin/python` (numpy 2.2.6, scipy 1.16.3, matplotlib 3.10.8)
- Build order chosen by user: **breadth of archetypes first** (CLI only in S1; GUI in S3)
- Disk: ~30 GB free as of 2026-09-13 20:30. The earlier 3.4 GB squeeze has cleared,
  so PySide6 is no longer a risk.

## Session plan — SCHEDULE HAS SLIPPED, READ THIS

S1 consumed both the interactive session and scheduled run 1 and 2, because every
spec-writing subagent died: 16 agent launches across three token windows produced
one usable file. The specs were ultimately written directly in the main loop, which
worked. **Do not delegate spec authoring to subagents.**

| # | Scope | Status |
|---|-------|--------|
| S1 | Core engine, registry, 40 archetypes, CLI, tests | **DONE** (2026-09-13 22:30) |
| S2 | Waveguides module + arrays (layouts, tapers, array factor) | PENDING |
| S3 | PySide6 GUI shell: catalog browser, param panel, plots | PENDING |
| S4 | Exporters (CST VBA, HFSS script), matching/network utils | PENDING |

Only ONE scheduled run remains (2026-09-14 01:12). It will pick up S2. **S3 and S4
will still be outstanding when it halts to await the user's `Continue`.** Say so
plainly rather than implying the build is finished.

## Invariants
- `otahub.core` must stay GUI-free and import-clean under plain numpy/scipy.
- Every archetype: declarative JSON in `specs/`, never Python.
- Every formula carries a `reference`, and every archetype at least one
  `known_cases` entry that becomes a live pytest case.
- Expectations of zero MUST use `tol_abs`; relative error is undefined there.
- Anything indicative rather than derived must say so in `notes`.

## Completed in S1
- `otahub/core/`: spec model, whitelisted AST evaluator (with Si/Ci/Bessel/elliptic),
  partial synthesis solver, registry, first-principles pattern maths.
- `otahub/cli/` + `OTA_Hub_AntennaToolkit.py` launcher: list / show / synth / check / doctor.
- 8 families, 40 archetypes: wire, loop, patch, horn, reflector, travelling_wave, uwb, slot.
- 241/241 known cases pass. Full suite green.

## Engine bugs found and fixed in S1
- All-or-nothing synthesis discarded whole designs over one unresolvable optional
  rule (79 of 100 known cases failed from this alone).
- Relative comparison against an expected zero reported a 3 ohm residual as 348% error.
- `db10()` used where a plain `log10()` fit was intended, inflating Yagi gain tenfold.
- `cosine_q` malformed expression; sidelobe finder mistook a second main beam for a sidelobe.

## Spec errors caught by the known-case harness
- quad_loop_square loss resistance out by 30.9x; halo_loop by 3.3x.
- Open-ended WR-90 gain 4.20 dBi, not 6.17. Ruze at 100 GHz 0.2057, not 0.2088.
- DRA radiation Q scales as eps_r^1.3 (2.46x per doubling), not 2x.
- Biconical Zc at 10 deg is 375.50 ohm, not 417.46.

## Next action
S2: waveguides (rectangular, circular, ridged, SIW; cutoff, modes, impedance) and
arrays (linear/planar layouts, uniform/binomial/Chebyshev/Taylor tapers, array factor,
grating-lobe limits). `otahub/core/pattern.py` already has the pattern machinery to
build array factors on.
