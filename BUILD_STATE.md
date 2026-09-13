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
| S2 | Waveguides module + arrays (layouts, tapers, array factor) | **DONE** (2026-09-14) |
| S3 | PySide6 GUI shell: catalog browser, param panel, plots | PENDING |
| S4 | Exporters (CST VBA, HFSS script), matching/network utils | PENDING |

The user typed `Continue` before the 01:12 run fired, and S2 was done in that
session. The 01:12 run, if it still fires, should pick up **S3 (GUI)**. S4
(exporters) remains outstanding after that.

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

## Completed in S2
- `otahub/waveguides/`: rectangular (modes, cutoff, guide wavelength, wave impedance,
  conductor and dielectric attenuation, power capacity, exact WR-series table),
  circular (Bessel-root cutoffs, TE01 low-loss mode), lines (coax, microstrip with
  exact-inverse synthesis, stripline, CPW, quarter-wave transformer).
- `otahub/arrays/`: uniform, binomial, Dolph-Chebyshev, Taylor n-bar and raised-cosine
  tapers; array factor, broadside directivity, beam steering, grating-lobe limits.
- CLI gained `guide`, `array` and `line` subcommands.
- Verified: WR-90 cutoff 6.557 GHz, attenuation 0.108 dB/m and power 1.05 MW all match
  datasheets; the 1.25x-1.90x band heuristic reproduces the published 8.2-12.4 GHz.
  Dolph-Chebyshev measured sidelobes equal the design level to 0.000 dB for every
  array size and level tested. Uniform d=lambda/2 directivity is exactly N.

## Bugs found in S2
- Dolph-Chebyshev even-length weights were independent of the design sidelobe level
  (missing half-sample phase term), and then rotated by one (ifft sign convention).
- `stripline_impedance` and `cpw_impedance` shipped dead `if False else` branches that
  double-applied the elliptic ratio. The limits settled the orientation: stripline
  needs K(k)/K(k'), CPW the inverse.
- WR table held rounded millimetres, perturbing WR-28's aspect to 1.997; now exact
  inch conversions.

## Corrections to my own assumptions in S2
- Not all WR guides are 2:1 - WR-90 is 2.25, WR-42 is 2.47.
- The WR number is hundredths of an inch, not mils, and is a ROUNDED nominal
  (WR-22 is really 0.224 in).
- Taylor is NOT more aperture-efficient than Dolph-Chebyshev for discrete arrays;
  measured, Chebyshev wins for n >= 20, consistent with its optimality proof.
  Taylor's real advantages are no edge spikes and decaying far sidelobes.

## Next action
S3: PySide6 GUI shell - catalogue browser over the 40 archetypes, parameter panel
driven by each spec's declared parameters, and pattern/impedance plots. Disk now has
~30 GB free so PySide6 installs without trouble. Keep `otahub.core` GUI-free.
