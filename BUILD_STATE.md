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
| S3 | PySide6 GUI shell: catalogue browser, param panel, plots | **DONE** (2026-09-14) |
| S4 | Exporters (CST VBA, HFSS script), matching/network utils | **DONE** (2026-09-14) |

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

## Completed in S3
- PySide6 6.11.2 installed (~2 GB; 28 GB free remains).
- `otahub/gui/`: three tabs - Catalogue (tree + search + generated form + Design/Sweep/
  Pattern results), Arrays (live taper designer), Waveguides (WR calculator with a
  dispersion plot). `otahub.core` remains GUI-free.
- The catalogue form is generated from each spec's declared parameters, so a new
  archetype in `specs/` gets a working UI with no GUI code change.
- `gui` subcommand added to the CLI and launcher.
- 40/40 archetypes verified to produce a non-empty design from form defaults.
- GUI tests run headless under QT_QPA_PLATFORM=offscreen.

## Bugs found in S3
- Catalogue filter delegated family rows to the base class, whose empty regex accepts
  everything, so every family survived any search. Family rows must return False and
  be pulled back in by recursive filtering.
- Form opened with an empty f0 (no spec declares a typical for it), so synthesis
  produced zero geometry and looked broken. Frequency is now prefilled from the
  geometric mean of each archetype's validity band.
- Polar plots left the lower half blank; patterns here are azimuthally symmetric, so
  a blank half reads as a one-sided pattern. Now mirrored.

## Completed in S4
- `otahub/utils/`: network parameters (S/Z/Y/ABCD, cascading, line transforms) and
  matching (L-section, quarter-wave, single-stub). Pozar 5.1 and 5.2 reproduced exactly.
- `otahub/export/`: neutral geometry IR with two backends (CST VBA, HFSS IronPython).
  Builders for 7 archetypes; everything else exports parameters only and says so.
- CLI gained `export` and `match`.

## Bugs found in S4
- L-section used sqrt(Z0/RL) where Pozar 5.3a has sqrt(RL/Z0) -> matched to 181 ohm.
- I then flipped a sign in the other branch that was already correct; deriving it
  directly showed X and B take the SAME sign there.
- Parameter export guessed lengths from magnitude, writing 319 ohm as 319105 mm.
  Now driven by each spec's declared unit.
- DesignResult.get() could not see requirements; the engine did not record units for
  the lambda0 and k0 it derives itself. Both fixed.

## Completed in S5 (catalogue breadth)

Ran on the user's `Continue` after the three scheduled runs. Scope followed the
original build-order choice: breadth of archetypes first.

- 40 -> 72 archetypes, 8 -> 10 families. 925 tests, 423/423 known cases.
- New families: `lens` (plano-hyperbolic dielectric, Luneburg, Fresnel zone plate,
  metal-plate) and `dielectric` (hemispherical and rectangular DRAs, plus
  `cylindrical_dra` moved out of `reflector` where it never belonged).
- Extended: wire +5, travelling_wave +4, reflector +4, patch +5, slot +3, uwb +3,
  horn +2.
- Evaluator gained Bessel functions of the second kind, needed for the annular
  ring's exact resonance condition.

## Errors found in SHIPPED specs during S5

- `corner_reflector_90` had AF = 4*sin^2(kS), putting a null at S = lambda/2 where
  the optimum actually is. Image theory gives 2*[1-cos(kS)]; the null is at
  S = lambda. It also used a field factor as a power factor and ignored the mutual
  coupling that sets the feed resistance. Rebuilt; now returns the textbook ~12 dBi.
- `short_dipole` labelled 20*pi^2*(L/lambda)^2 as the uniform-current value. It is
  the triangular one. The second metric then quartered an already-triangular value,
  under-reporting a real short dipole by 4x.
- `cassegrain` carried `magnification = 1.0` as a hard-coded placeholder.

## Test-guard bugs found in S5

- `_bounds_for` classified by substring before suffix, so `efficiency_db` was
  bounded to [0,1] and any negative decibel value failed.
- The dBi floor rejected genuine pattern nulls; a dipole lambda/2 over ground has
  an exact zenith null and reported -310 dBi.
- The degree bound assumed angles are unsigned. A beam angle measured from
  broadside is negative when the beam scans the other way.

## CLI bug found in S5

- `array --sll 30`, which is how everyone says "30 dB sidelobes", hit an unhandled
  ValueError. Either sign is now accepted, and an unsynthesisable taper reports
  cleanly instead of raising.

## Lessons that generalise (S5)

- **Hand arithmetic was the single biggest source of false failures.** The formulas
  were nearly always right; my expected values were not. Roughly a dozen known
  cases failed on my arithmetic alone. Compute expectations with a short
  independent script instead - once I did that for the lens family, it landed with
  zero failures on the first run.
- **String-interpolating expressions needs parentheses.** Splicing the guide
  wavelength in unparenthesised made `(c/f0)/LAMG` parse as `1/sqrt(...)`,
  inverting lambda0/lambda_g and inflating a conductance by 76%.
- **Check which branch a root-finder landed on.** The hemispherical DRA took three
  attempts: the first two converged on higher-order modes and produced a Q that
  fell with permittivity, which is physically backwards. Scanning for the FIRST
  peak of the Mie coefficient fixed it, and the result then independently
  reproduced the published eps_r^1.3 scaling.
- **Physical sanity checks catch what algebra does not.** A V antenna modelled with
  both legs carrying outward current has an exact null on its axis. It is a flared
  transmission line: the return conductor's current runs against propagation.

## Dimensional audit added in S5

`tests/test_scale_invariance.py`. Scale every frequency by S, every length by
1/S and conductivity by S: the antenna is electrically identical, so every
quantity must follow the power of S its declared unit implies. 793 quantities
across all 72 archetypes, with no reference data at all.

Found two defects in specs that had passed their own cited cases since S1:
`circular_patch` carried an intermediate radius in centimetres while declaring
it dimensionless, and `half_wave_slot` had an `m2` unit typo. Verified to fire
by injecting a fixed 1 mm into a dipole length.

Worth knowing: the first version of the harness had its own bug - lower-casing
units folded siemens ("S") together with seconds ("s"), giving conductance the
time scaling rule. Unit strings need case-sensitive matching.

## Planar arrays added in S5

`otahub/arrays/planar.py` plus a `planar` CLI subcommand. Rectangular and
equilateral-triangular lattices, separable tapers, steering, scan loss, and
directivity for an arbitrary element layout.

Directivity is exact, not integrated: the sphere average of exp(j k.d) is
sin(kd)/(kd), so radiated power is a double sum over sinc(2|r_m - r_n|/lambda).
It reproduces the linear module to machine precision for a single row.

The triangular lattice's 13.40% element saving is derived from the reciprocal
lattice (shortest vector 4*pi/(sqrt3*s) against 2*pi/d), not asserted, and a
test recomputes it from the primitive vectors.

Three bugs in my own new code, all caught by writing the checks first:

- The steered directivity summed the STEERED weights for the peak. The peak is
  at the scan angle, where the steering phase cancels, so it must be the bare
  weights - the wrong version reported a 45-degree scan as -4.7 dBi instead of
  18.2.
- The pattern cut put broadside at theta = 0, the grid edge, so the beamwidth
  finder returned NaN and the sidelobe finder 0 dB. Broadside belongs in the
  middle of the cut where both flanks are visible.
- The scan-plane cut swept a great circle outward from the beam, which runs
  past the array plane into the MIRROR beam - full amplitude for isotropic
  elements - and pinned the sidelobe reading to 0 dB at every spacing. It now
  sweeps the signed angle from the normal and stays in the forward hemisphere.

The last one is the same failure mode as the S1 sidelobe bug: a second main
beam mistaken for a sidelobe.

## Export builders extended in S5

7 archetypes had geometry builders against a catalogue of 40; the catalogue
then grew to 72 and the gap widened. Now 18 of 72.

Three primitives added to the neutral IR and rendered in both backends: Cone
(truncated, so one radius may be zero), Sphere, and Subtract as a boolean
operation applied after the solids exist.

New builders: folded_dipole, dipole_over_ground, turnstile_dipole (two ports),
quarter_wave_shorted_patch, rectangular_dra, cylindrical_dra,
hemispherical_dra (sphere minus a half-space, since neither tool has a
hemisphere primitive), conical_monopole, biconical, discone.

Two tests worth keeping: booleans must only name solids that exist, and they
must be emitted after every solid they operate on. Both are silent failures in
the simulator otherwise - a subtract on a missing object is a runtime error
inside the script, not something the exporter would notice.

Still parameters-only: horns (need a loft or truncated pyramid), loops (need a
torus), and Yagi-Uda - that last one for a different reason. The spec gives
boom length, reflector and driven lengths and a director count, but not the
individual director lengths, so the geometry is genuinely underdetermined and
building it would mean inventing dimensions.

## GUI planar tab added in S5

A fourth tab, and the linear one renamed so the two are distinguishable.
Controls for lattice, spacing, taper, steering in both angles, and a
ground-plane switch; a layout scatter sized by excitation amplitude; and the
two cuts through the beam on one polar plot.

The tab says what it cannot do rather than hiding it: a triangular lattice is
not separable, so the taper control disables itself and a note explains why.
Ignoring the taper silently would have been the easy option and the wrong one.

The polar plot spans the forward hemisphere only, and that empty lower half is
correct here - unlike the S3 dipole plots, where a blank half meant the
mirroring was missing.

## GUI bug found in the S5 verification pass

Eight archetypes opened in the GUI with a blank frequency field and produced no
geometry - the whole wideband family, whose requirement is `f_low` rather than
`f0`. The S3 fix that prefills a design frequency keyed on the NAME "f0", so it
never reached them. Now keyed on the declared unit.

It survived because the test summing geometry and metric rows let them pass on
their constant metrics alone - a fixed bandwidth ratio and directivity - while
every dimension came out blank. The test now checks geometry specifically.

This is the same failure the S3 notes describe ("form opened with an empty f0,
so synthesis produced zero geometry and looked broken"), in a place the S3 fix
did not reach.

## BUILD STATE

S1-S5 done. 72 archetypes, 10 families, 1235 tests, 423/423 known cases.
See `docs/HANDOVER.md` for what is verified, what is not, and next steps.

The largest untested surface is unchanged: neither exporter has been run against a
real CST or HFSS installation. Structurally validated only.
