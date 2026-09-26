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

## Inset-patch discrepancy resolved (post-S5)

Open since S1 and recorded in the spec's own validity block: the exact integral
of Balanis 14-18a was said to give 212.5 ohm against the published 228.35, a
6.9% disagreement.

There was no disagreement. Direct quadrature of 14-12 and 14-18a reproduces
Example 14.2 to 0.05% - G1 = 1.5747e-3 S against the published 1.57e-3, and
G12 = 6.1651e-4 against 6.1683e-4. The 212.5 ohm figure was simply a bad
integration.

The spec's 1.7% agreement was two errors cancelling: G1 used the (1/90)(W/lam)^2
small-width branch, 10% high, while G12 = G1*J0(k0*L) was 28% low. Both fixed.
G1 is now the exact closed form; G12/G1 is a polynomial fitted to the exact
integral over k0L in [0.8, 3.3] and k0W in [1.15, 3.25], the range eps_r from 1
to 12 produces. Example 14.2 now reproduces to 0.066%, and two further known
cases pin the fit at the other end of the design space.

Worth remembering: an agreement that looks acceptable can be two errors
cancelling, and the only way to tell is to check each term on its own.

## Cross-consistency audit added (post-S5)

tests/test_cross_consistency.py: 23 relationships between archetypes that must
agree because the same physics reaches two specs by different routes. Monopole
against dipole, folded against plain, slot against its complementary dipole,
conical monopole against biconical, slot array against single slot, Potter horn
cutoffs against the circular-guide Bessel zeros, and the three DRA shapes -
exact Mie, published fits, magnetic-wall algebra - against each other.

All clean on the first run. Verified to fire by injecting a 3% drift into
half_wave_dipole's radiation resistance: six tests across three families
failed, which is the blast radius that makes it worth having. The known-case
harness flags the one archetype; this shows how far the error reaches.

Together with the dimensional audit the catalogue now has two systematic checks
that need no external reference. Neither covers an archetype with no sibling
and no dimensional quirk; those still need reading against their source.

## Rectangular DRA resonance corrected (post-S5)

The all-magnetic-wall model is gone, replaced by the dielectric waveguide model:
magnetic walls on the four sides, the open top treated properly, so kz solves
kz*tan(kz*h) = sqrt((eps_r-1)*k0^2 - kz^2).

Non-dimensionalising first made it tractable. With u = kz*h and
P = pi^2*(1/aw^2 + 1/aL^2), the condition collapses to
u*tan(u) = sqrt(((eps_r-1)*P - u^2)/eps_r) - a function of P and eps_r alone,
with the aspect ratios entering nowhere else. Fitting u over P in [2,20] and
eps_r in [6,50] gives 0.032% error in the resulting size.

I had the limit backwards at first and it is worth recording why. I assumed the
magnetic wall was the eps_r -> infinity limit. It is not: k0 shrinks with eps_r
too, so the right-hand side tends to sqrt(P), which is finite, and u settles
below pi/2 for any permittivity. The magnetic wall is a cruder model, not a
limiting case - so the oversize does not vanish at high eps_r, and the spec now
says so explicitly.

The validation is the useful part. Correcting the resonance brings the brick
from 1.41x the volume of the hemispherical and cylindrical archetypes down to
1.06x, holding within 1.14x across eps_r from 8 to 40. Those two come from
exact Mie theory and from published curve fits respectively - three unrelated
models converging. The cross-consistency guard on DRA volumes was loosened to
2.5x to accommodate the old brick, which made it nearly useless; it is now 1.25x.

Q is still borrowed from the hemisphere. That is the one substantive
approximation left in this family.

## Touchstone import added (post-S5)

otahub/utils/touchstone.py plus a `touchstone` CLI subcommand. The exporters
send a model out to CST or HFSS; this is the way back.

Reads v1.0 and v1.1, MA/DB/RI formats, every frequency unit, S and Z
parameters, any port count. G and H files are recognised and refused rather
than mis-converted. Writes too, and the round trip is tested at 1, 2 and 3
ports in all three formats.

The two traps, both silent failures, both tested:

- Two-port files are COLUMN-major (freq S11 S21 S12 S22) while three ports and
  up are row-major. Reading a 2-port row-major transposes it, swapping forward
  gain with reverse isolation. The test data is deliberately asymmetric so a
  transpose cannot hide.
- A frequency point may wrap across any number of lines, so the reader works on
  a flat stream of values chunked by 1 + 2N^2, not line by line.

Port-count inference turned out to need both signals. A nine-value first line
is either a 2-port (whole matrix on one line) or the first row of a 4-port -
a real ambiguity in the format, and why the .sNp suffix exists. The reader
breaks the tie on the total value count where only one candidate divides
evenly, and says so plainly where both do.

## Horns rebuilt on exact theory (post-S5)

Found by surveying the catalogue for weak verification: the two sectoral horns
sat at the bottom, one known case each asserting one quantity out of four
produced.

They pinned the aperture efficiency at 0.65 and 0.64 and synthesised the
aperture so the horn was always optimum by construction. Self-consistent, and
useless for analysing a horn that already exists. Worse, with efficiency fixed
the gain rose without limit as the flare grew - the old spec would have
recommended over-flaring indefinitely.

Balanis gives the directivity exactly in Fresnel integrals. Rearranged:

  eta_E = (8/pi^2)*[C(q)^2 + S(q)^2]/q^2,   q = b1/sqrt(2*lam*rho)
  eta_H = (lam*rho/a1^2)*{[C(u)-C(v)]^2 + [S(u)-S(v)]^2}
  eta_P = (pi^2/8)*eta_E*eta_H

The 8/pi^2 is the TE10 cosine taper across the unflared plane. In the pyramidal
form it appears in both sectoral efficiencies and must be removed once - that
factor is exactly what my first independent check was missing, which showed up
as a CONSTANT 0.912 dB offset at every parameter. A constant offset means a
constant factor, and pi^2/8 = 1.2337 is 0.912 dB.

Optimum values 0.64870, 0.64276, 0.51440 against the pinned 0.65, 0.64, 0.51.
Fresnel integrals added to the evaluator whitelist (scipy returns (S, C) in that
order, which is the easy thing to reverse).

Verified three ways, agreeing to 1e-14: the closed form, Balanis 13-19 and
13-41 as published, and direct aperture integration. The spec keeps the
published directivity expression as a live metric so the two routes are checked
against each other on every run.

## circular_patch directivity integrated (post-S5)

Next on the weakness survey after the horns. It carried db10(6.3) with no
derivation at all. The TM110 far fields are closed in J0 -/+ J2 (Balanis
14-79/14-80), so the directivity is one quadrature over the upper hemisphere.

The hard-coded value was 14% high on eps_r 2.2, 54% high on FR-4 and 84% high
on eps_r 10.2 - worst exactly where circular patches get used, since high
permittivity is why you pick one. True values: 5.50, 4.08, 3.43.

Everything depends on k0*a_e alone, and the resonance condition k*a_e = 1.8412
ties that to eps_r, so it is a one-variable fit. Worst error 0.19%. Also added
the feed-position law R(rho0)/R_edge = [J1(k*rho0)/J1(k*a_e)]^2, which is
permittivity-independent for the same reason.

The radiation conductance is derived from the same fields but is NOT checked
against any published example, and says so in its notes. The directivity and
the feed law are checked; those are the numbers to lean on.

Bug in my own work worth recording: I built the fit expression by string
replacement and `.replace("**0", "")` turned `k0ae**0` into `k0ae`, making the
constant term linear. The known cases caught it immediately - three failures at
18%, 12% and 48%. Build expressions term by term, not by patching strings.

## Conical and corrugated horns integrated (post-S5)

The last two pinned efficiencies in the horn family: 0.51 smooth, 0.69
corrugated. Both now integrated from their aperture fields - TE11 for the
smooth cone, the balanced hybrid HE11 (a J0(2.405*rho/a) taper) for the
corrugated one - with quadratic phase across the flare.

Two validations fell straight out. The uniform-phase limits are 0.836829 and
0.691660, against 0.836 and 0.69 in the literature. The corrugated horn's 0.69
is therefore DERIVED here, not quoted.

The 0.51 puzzle resolved rather than being papered over. The integral gives
0.53847 at the sqrt(3*lambda*L) flare, not 0.51. The two are both right: 0.51
is the efficiency at the TRUE maximum-gain flare, s = 0.3908, while the
sqrt(3*lambda*L) rule puts s at exactly 0.375. They differ by 4% in efficiency
and 0.0078 dB in gain - which is what "optimum" means, the peak being flat. The
old spec paired one convention's efficiency with the other's aperture, costing
0.24 dB.

Found while doing it: gain_advantage_over_smooth_db shipped as
db10(0.69/0.51) = +1.31 dB, comparing the corrugated horn at ZERO phase error
against the smooth horn at its optimum. Like for like at the same flare the
corrugated horn is slightly BEHIND, -0.17 dB at the optimum, because the J0
taper is heavier. It pulls ahead only when over-flared, where its taper
tolerates phase error better. Corrugation buys pattern symmetry and
cross-polarisation, not gain, and the spec now says so.

Also: my uniform-phase limit case used flare = 0.0001, which describes a horn
a fraction of a millimetre across and a 109,000 degree beamwidth. The
plausibility guard caught it. Limit checks still have to describe real objects.

## Rectangular patch directivity integrated (post-S5)

Found by sweeping the catalogue for directivities that are bare constants. Of
33 such expressions most are legitimately exact - a short dipole really is 1.5,
a ground-plane DRA really is 3.0 - but `6.6` was shared by THREE patch specs
and was a flat number for every substrate.

Derived from the two-slot model instead. Each radiating edge is a uniform
magnetic line current of length W, giving a single slot D1 = (k0*W)^2/I1 over
the half space, and the pair D = 2*D1/(1 + G12/G1). Both pieces already existed
in the inset patch's conductance work.

  eps_r  2.2   true 5.44 (7.36 dBi)   the 6.6 was 21% high
  eps_r  4.4   true 4.06 (6.08 dBi)   63% high
  eps_r 10.2   true 3.47 (5.40 dBi)   90% high

and LOW for a wide patch on thin low-permittivity board, where it passes 7.8.

Two checks: direct 2-D pattern integration agrees to 0.000% at five aperture
sizes, and the narrow-slot limit gives D1 = 3.00001 - a magnetic dipole's 1.5
doubled by the ground plane.

A convention split worth knowing, and now stated in both specs: the slot
SEPARATION for the pattern is the EFFECTIVE length L + 2*dL, because the
equivalent currents sit at the fringing edges, while the mutual conductance
that sets the input resistance uses the PHYSICAL L - which is what reproduces
Balanis Example 14.2. Balanis uses both in one chapter. Using L_eff for the
impedance instead would move it from 228 to 261 ohm.

A test had to be retired. It asserted that a circular patch is less directive
than a rectangular one, which was not physics but an artifact of comparing two
invented constants, 6.3 against 6.6. With both integrated they agree to within
2% and cross over around eps_r 6. Asserting an ordering there would be
asserting numerical noise; the test now checks that they track each other and
that both fall with permittivity.

## Shorted patch: what shorting actually costs (post-S5)

quarter_wave_shorted_patch has ONE radiating edge, so its directivity is exactly
the single-slot term the two-slot model is built from, (k0*W)^2/I1. Dropped
straight in from the previous iteration's work.

The more interesting finding is the claim the spec carried alongside it: that
one slot instead of two costs "about 3 dB". Computed, the cost is 2.16 dB on
eps_r 2.2, 1.06 dB on FR-4 and 0.50 dB on eps_r 10.2. Always LESS than 3, and
varying by 1.65 dB across ordinary substrates.

The reason: the two slots sit well under a half wavelength apart - 0.36
wavelengths at eps_r 2.2, 0.17 at eps_r 10.2 - so they never arrayed perfectly.
On high-permittivity board they are nearly coincident and the second slot was
adding almost nothing, so removing it costs almost nothing. Shorting a patch is
much cheaper on high-permittivity board than the rule of thumb suggests.

## Loop geometry builders (post-S5)

18 -> 24 builders. The loop family was the largest with none at all, and it
turned out to need only one new primitive: a torus, which both CST and HFSS
expose directly. That is the difference from the horns, where a loft would mean
guessing at the API's STRUCTURE (CST picks faces by id, HFSS wants a polyline
then a loft) rather than just its parameter names.

Built: small_circular_loop, one_wavelength_circular_loop, small_square_loop,
quad_loop_square, alford_loop, halo_loop. Square loops needed no new primitive
at all - four cylinders, one of them split for the feed.

The two backends specify a torus differently and the difference is silent: CST
takes inner and outer radii measured from the axis, HFSS takes major and minor.
Feed either one the other's numbers and you get a ring of the wrong size with
nothing in the file looking wrong. A test pins both forms.

Feed gaps are cut with a real boolean rather than drawn, because a port across
an unbroken ring shorts itself out. The halo is the exception worth noting: its
gap is a DESIGN parameter, since the tip capacitance sets the resonance, where
every other loop's gap is invented by the exporter and says so.

Still nothing for reflector, slot, travelling-wave or lens.

## Slot geometry builders (post-S5)

24 -> 30 builders, and the slot family is now complete at 6/6. No new
primitives at all: every one is bricks with boolean cuts, both already in use.

A slot has to be SUBTRACTED from its ground plane. Drawn as a separate solid it
is not a slot, it is a plate with a bar lying on it - and nothing in the
exported file would look wrong. Four tests pin that, including that the
waveguide cuts span the full wall thickness rather than leaving a blind hole.

The resonant array's alternating offsets get their own test, because the
alternation IS the mechanism: it undoes the 180 degrees of propagation phase
between slots half a guide wavelength apart. Build them all on one side and the
array splits into two beams off broadside instead of one on it - a mistake that
produces a perfectly plausible-looking model.

Each builder states what it chose that the spec could not supply: slot widths
(the specs give only lengths), cavity lateral dimensions, and for the
travelling-wave array the fact that its offsets are UNIFORM where a real one
tapers them along the guide.

Left without builders: reflector, travelling-wave and lens. Parabolic dishes
need a swept profile and horns need a loft, both of which mean guessing at API
structure rather than parameter names.

## Seven more geometry builders (post-S5)

30 -> 37. Found by asking which of the remaining archetypes are buildable with
the primitives ALREADY in the IR, rather than assuming the empty families all
needed new ones. Most did not.

long_wire_travelling, leaky_wave_line_source, planar_monopole_rectangular,
planar_monopole_circular, annular_ring_patch, pifa, stacked_patch. Bricks,
cylinders and boolean cuts throughout.

Two worth noting. The long wire exports TWO ports, the second standing in for
the TERMINATION at a non-50-ohm impedance - a travelling-wave wire with an open
far end is a standing-wave wire and its pattern splits, so leaving that implicit
would produce a model that runs and answers the wrong question. And the annular
ring is CUT from a disc rather than drawn as two, for the same reason the slots
are.

New test: an archetype whose spec is flagged low confidence must carry that
warning into the exported file. pifa, stacked_patch, planar_monopole_circular
and halo_loop all do. Without it a shaky model arrives in the solver looking as
solid as any other.

Caught while adding the cases: I passed eps_r to the PIFA, which declares no
such parameter. The unit test on exported parameters caught it - an undeclared
requirement reaches the exporter with no unit and cannot be classified.

Left: reflector 0/7 and lens 0/4, plus the horns. All need swept profiles or
lofts, which means guessing at API structure rather than parameter names.

## Lens geometry builders (post-S5)

37 -> 39. The lens family looked like it would need swept profiles; two of the
four turned out to be concentric shells, which Sphere, Cylinder and Subtract
already cover.

fresnel_zone_plate builds the metal rings of a Soret plate, with zone edges
from the EXACT r(m) = sqrt(m*lambda*F + (m*lambda/2)^2) rather than the
sqrt(m*lambda*F) approximation - the difference matters when F is only a few
wavelengths. A test pins the builder's radii against the spec's own r_first and
r_outer, so the two independent copies of that formula cannot drift apart.

luneburg_lens builds the eight concentric shells a real one is made from, with
boundaries at equal steps in r^2 rather than in r. Permittivity is exactly
linear in r^2, so that spacing gives every shell the same permittivity span;
equal steps in r would not. A test asserts the r^2 spacing directly.

Both lenses export NO port and say why: they are illuminated by a separate
feed, and a Luneburg lens can carry several at once, which is the reason to
build one. A model with no port and no explanation reads as unfinished.

Left at zero: reflector 0/7, plus the hyperbolic lens and the horns - all
needing swept profiles or lofts.

Worth recording about the session itself: repeated classifier timeouts dropped
several large heredocs before they ran. Checking with grep rather than assuming
the edit landed, and falling back to the Write tool for the new test file, kept
the work from being silently lost.

## Wire family complete, loop nearly (post-S5)

39 -> 43. wire is 11/11 and loop 7/8, both with primitives already in the IR.

dipole_arbitrary_length needed no new builder at all - it is the same two arms
and a gap, so it joined the existing dipole decorator stack rather than getting
a copy that would have to be kept in step. A test asserts the two produce
identical geometry at 0.5 lambda.

Three builders export an extra port, and in each case leaving it out would give
a model that runs cleanly and answers a different question:

- inductively_loaded_monopole's port 2 IS the loading coil's gap. Left open the
  whip is just a short whip and the spec's resonance never happens. The note
  names loading_inductance_H and coil_q so the reader knows what to put there.
- long_wire_travelling's port 2 is the termination, from the earlier round.
- multiturn_small_loop draws N separate rings, which solved as-drawn is one
  driven ring and N-1 parasitic ones. The spec's N-squared radiation resistance
  depends entirely on the turns being in series, so the note says CONNECT THE
  TURNS in as many words.

Two honesty notes went into the specs' exported files. top_loaded_monopole's
hat radius and its beta_top are BOTH inputs to the spec and are not derived
from one another, so the model does not imply the drawn hat produces the
assumed current distribution. multiturn_small_loop's pitch comes from the wire
radius rather than the spec's l_coil, so a loosely wound coil will not match.

ferrite_rod_loop stays unbuilt: its geometry does not resolve without
ld_target and L_target, and it needs a ferrite material model besides.

## The annular ring's directivity is two-dimensional

The ring carried a flat `db10(5.0)`. It is now integrated, and the integration
turned up something that changed how it had to be shipped.

The ring radiates from BOTH edge walls. Each is a phi-directed magnetic ring
current with a cos(phi) dependence - the same source that gives the circular
patch its J0 -/+ J2 pattern - so the far field is the superposition of two such
rings, weighted by the cavity eigenfunction

    R(rho) = J1(k*rho)*Y1'(k*a) - Y1(k*rho)*J1'(k*a)

at each radius, and OPPOSITE in sign, because the two outward normals oppose.
That sign is the whole character of the answer: the two rings partly cancel, so
the ring sits just below a solid disc of the same outer radius rather than
above it. Quadrature and an independent 2-D angular grid agree to four decimals.

The flat 5.0 was wrong in both directions - 5% LOW at eps_r 2.2 and 48% HIGH at
eps_r 10.2. A constant cannot be merely conservative when the true value crosses
it.

The fit is in the two RADII, not in the ratio b/a, and that was not the first
choice. A 1-D fit in k0*b failed no matter the degree, and the diagnostic said
why: D is not a single-valued function of k0*b. Two designs reaching the same
k0*b from different (eps_r, b/a) pairs differ by as much as 13.5%, so a
one-variable fit has to average away a real 13% spread. In (k0*a_in, k0*b_out)
a degree-3 fit already beats degree-6 in (k0*b, ratio) - 0.042% against 0.028%
for twice the terms. Shipped: degree 4, 15 terms, worst error 0.0046% over
eps_r 2-13 and b/a 1.2-3. The validity block says the fit diverges rather than
degrades outside that box, because a bivariate polynomial does.

`tests/test_ring_directivity.py` pins both facts: that the shipped fit tracks
its own quadrature, and that the 2-D-ness is real - two designs sharing a k0*b
to within 10% whose directivities still differ by more than 5%. If someone
later "simplifies" this to a function of the ratio, that test fails.

One cross-check earned its place: the ring must land below the circular patch
on the same board, but not far below - two independently derived patterns
agreeing to within 10% is worth more than either alone.

Remaining bare constants after this: `triangular_patch` (5.0) and `biconical`
(1.6409).

## The triangle's mode, found rather than quoted

`triangular_patch` carried `db10(5.0)` with the note "roughly 1 dB below a
rectangular one because its aperture is smaller. Placeholder, not computed."
Both halves of that sentence turned out to be wrong.

The cavity model needs the dominant Neumann eigenfunction of an equilateral
triangle. Rather than copy a closed form out of a handbook and hope the
placement of the axes matched, it was found here: the six plane waves of
magnitude k = 4*pi/(3a) sitting on the hexagonal star span a two-dimensional
null space of the magnetic-wall condition - two-dimensional because the mode IS
degenerate, which is the reason a triangle can be made to radiate circularly -
and the member symmetric about a median is the one a probe on that median
excites. The Neumann residual on the walls comes back at 2e-14.

All THREE walls radiate, as magnetic line currents M = 2*Ez*(z_hat x n_hat).
Ez restricted to a straight side is a sum of twelve exponentials, so each
wall's far-field integral is elementary and exact - no edge discretisation at
all. The first version sampled the walls numerically and took 84 seconds per
sweep; the closed form takes 40 for the whole design space and agrees with it
to five digits, which is also a check on both.

Two things fell out that are worth more than the number itself:

- The fringing correction cancels exactly out of k0*a_eff. The synthesis sets
  a_eff = 2c/(3*f0*sqrt(eps_r)), so k0*a_eff = 4*pi/(3*sqrt(eps_r)) whatever
  the substrate height. Directivity is therefore a function of ONE variable,
  and eps_r alone fixes it. That is shipped as its own metric so the claim is
  checkable rather than assumed.
- The small-patch limit is exactly 3. A patch much smaller than a wavelength is
  a horizontal magnetic dipole over a ground plane, and the raw integral
  returns 3.0000039 at k0*a_eff = 0.0042. Nothing in the chain was set up to
  make that happen, so it verifies the mode, the wall currents and the
  hemisphere integral at once - with no reference data. The shipped form is
  D = 3 + (k0*a_eff)^2*g(k0*a_eff) so that the limit survives outside the fit
  range, where a bare polynomial would wander off. Degree 5 in g, worst error
  0.028% over eps_r 1 to 21.

The flat 5.0 was 43% LOW on an air substrate and 47% HIGH on eps_r = 10.2.

And the old note's reasoning was backwards. The triangle does not run 1 dB
below a rectangular patch - it lands within 0.1 dB of one at every permittivity
tested, +0.09 dB on air and -0.04 dB on eps_r 10.2. So the 36% board-area
saving costs essentially nothing in directivity, which is a better argument for
the shape than the spec was making for it. The summary line was corrected to
say so. Three patch shapes now agree within 5% at eps_r 2.2 and all converge on
3.0 as the substrate gets denser, by three independently derived patterns.

Remaining bare constant: `biconical` (1.6409).

## A full-wave solver, and a reference that was not what it looked like

A sweep of the whole catalogue for metrics whose expression is a bare number
found 66 of them, 37 not labelled indicative. Most are exactly right and should
be constants - a short dipole's 1.5, a half-wave dipole's 78.08 degree
beamwidth, an offset parabolic's blockage efficiency of 1. A handful are not:
`one_wavelength_circular_loop` asserts 3.4 dBi and 100 ohm, `quad_loop_square`
3.3 dBi and 120 ohm, and neither has a closed form to check against or a
citable number precise enough to use. Deriving them needs a full-wave solver,
which is also what HANDOVER section 6 item 4 - validate a low-confidence
archetype end to end - had been blocked on since session 1.

So `otahub/num/mom.py`: thin-wire method of moments, EFIE in mixed-potential
form, rooftop basis, Galerkin testing. The kernel's singular part is removed
analytically rather than quadratured, which matters more than it sounds: with
a = 0.001 lambda and segments of 0.006 lambda the kernel peak is six times
narrower than a segment, and plain quadrature loses the self term silently.

The afternoon went on one thing. The solver returned 86 ohm for a half-wave
dipole where every textbook says 73.08 + j42.52, and the discrepancy GREW under
mesh refinement - the signature of a convergence bug. It was not one.

- Driving the matrix with a prescribed sinusoid returned 73.083 + j42.494. So
  the matrix was right, and the solve was the suspect.
- Raising the quadrature from 24 points to 160 changed the answer by nothing
  whatsoever. So it was not quadrature either.
- An independent Hallen solve - a different integral equation, no divergence
  term, collocation instead of Galerkin - returned 86.6 ohm.

Two formulations sharing no algebra agreed with each other and disagreed with
the constant. 73.08 + j42.52 is the INDUCED-EMF value for an ASSUMED sinusoidal
current on a VANISHINGLY THIN wire. It is not the driving-point impedance of a
delta-gap-fed wire of finite radius, and the real current is visibly fatter than
a sinusoid near the ends - 0.385 against 0.309 at z = 0.2 lambda. The lesson is
that a famous number can be a different quantity wearing the same name, and the
way to find out is a second formulation rather than a finer mesh.

Six checks now licence the solver, chosen so each touches a different part:
the matrix alone against the induced-EMF integral; the far-field code alone
against the same number by a different route; the pattern integral against the
circuit power, to 1 part in 1e5; D = 1.5 for a short dipole and a small loop;
Rr -> 20*pi^2*(C/lambda)^4 as a loop shrinks; and the folded dipole's 4.014:1,
which is also the only geometry with corners and two parallel conductors.

Two real bugs found on the way, both in my own code:

- `dipole()` defaulted to an ODD segment count, and an odd count has no node at
  the centre - so the feed sat half a segment off and broke the symmetry of the
  problem. It now rounds up, and a test asserts the feed node is at z = 0 for
  odd and even requests alike.
- The far-field routine looped in Python over every angle, which made the
  validation suite take 18 seconds. Vectorising the segment phasor and
  accumulating per segment rather than per basis function took it to 1.4. The
  matrix build got the same treatment - one pass over segment pairs computing
  the four linear moments, instead of recomputing each pair up to four times -
  which took a 60-segment solve from about 10 seconds to 0.02.

What it does not do, and these are limits rather than approximations that wash
out: no ground plane, no dielectric, no loss, no junction of more than two
wires, and the reduced rather than the exact kernel.

## The full-wave loop, and three numbers that were all wrong

With the solver in hand, `one_wavelength_circular_loop` was the first place to
point it: the spec asserted a fixed 1.09*lambda resonance, 100 ohm and 3.4 dBi,
none of them derived from anything.

Because a single solver proves nothing, the arbiter came first. The circular
loop has a classical Fourier-mode solution - expand I(phi) in exp(j n phi), and
the phi-directed EFIE separates mode by mode into

    I_n = 4 k V / ( j eta [ k^2 b^2 (G_{n-1} + G_{n+1}) - 2 n^2 G_n ] )

with G_m the Fourier coefficients of exp(-jkR)/R around the ring. It was
derived here rather than copied, and it discretises no geometry at all, so it
shares nothing with the method of moments but Maxwell. The two agree to 0.25%
on impedance and to four decimals on the current distribution, which is a
stronger statement than either alone.

All three asserted numbers are wrong, and not in a way a tolerance would hide:

- The resonant circumference is NOT fixed. It runs 1.0384 lambda on 1e-4 lambda
  wire to 1.0937 on 3e-3 - the flat 1.09 is 2.6 to 5.0% long.
- The driving-point resistance is 138 to 151 ohm. The asserted 100 ohm is 28 to
  30% LOW. Sanity: the square quad loop, whose ~120 ohm is well attested, comes
  out at 126-136 ohm from the same solver, and a circle encloses more area than
  a square of equal perimeter, so the ordering is right.
- Broadside directivity is 3.61 to 3.77 dBi. The asserted 3.4 is 4.7 to 6.2%
  low, and the gain over a dipole is 1.46-1.62 dBd rather than 1.25.

The sharpest finding is the one that was written down twice. The spec carried a
thickness correction `1.09 - 0.0045*(20 - Omega)` AND a validity note saying
"thicker conductors resonate at a shorter circumference". Both solvers say the
opposite: thicker resonates LONGER. The old correction was +6.0% at thin wire
and -5.8% at thick, crossing the truth in the middle, which is how it could
look plausible on a spot check. `test_a_thicker_conductor_resonates_at_a_longer
_circumference` pins the direction now.

The thickness parameter itself was wrong too - defined as 2*ln(2*pi*C/b) where
the loop's Omega is 2*ln(2*pi*a/b) = 2*ln(C/b), an extra 2*pi inside the
logarithm putting it 3.676 high.

Everything now keys off that one corrected number, with degree-3 fits at 0.061%
(circumference), 0.209% (resistance) and 0.051% (directivity).

Two checks worth keeping for their own sake. The far-field routine reproduces
the EXACT uniform-current loop directivity, 2*k*a*J1(k*a)^2 / int J2, to 0.01%
on this same ring geometry - a validation of the radiation integral on a curved
path rather than a straight one. And broadside being the main beam is not
geometry, it is the current: with a UNIFORM current a one-wavelength loop peaks
in its own plane instead, 90 degrees away. The test asserts the beam direction
for that reason.

The dimensional audit earned its keep again, on this session's own work: it
rejected `dBd` as a unit with no scaling rule the moment it appeared. It is a
legitimate unit - a ratio of two directivities, invariant for the same reason
dBi is - so the rule went into the table rather than the unit coming out of the
spec.

Still constants in the loop family: `quad_loop_square` (3.3 dBi, 120 ohm, and a
1.0218*lambda perimeter which the solver says is not resonant at all - it reads
107 - 108j there), `halo_loop` and `alford_loop`. The quad is next; the halo
needs a lumped capacitive gap, which this solver cannot yet model.

## The quad loop, verified by continuity rather than by a second formula

`quad_loop_square` asserted 1.0218*lambda, 120 ohm and 3.3 dBi. A square has no
modal solution to check it against, so the question was how to verify it at all
without simply trusting one solver.

The answer was continuity. A regular N-gon of a given perimeter tends to the
circle of the same perimeter, and the circle HAS an independent Fourier-mode
solution. The N-gon uses exactly the code path the square does - corners,
non-collinear segments, a closed wrap - so if the sequence lands on the modal
circle, that whole path is anchored to something outside this solver, and the
square is simply N = 4 of it. It lands: 4 sides gives 114.3 - 69.0j, 8 gives
129.4 - 38.1j, 48 gives 136.4 - 21.3j, 96 gives 136.5 - 21.3j, against the
modal circle's 137.1 - 19.4j. The test asserts the approach is monotone as well
as that the limit is right.

Of the three asserted numbers, one was wrong, one was close and one was right:

- The perimeter is wrong for the stated context. 1.0218*lambda is the ham-radio
  1005/f(MHz) feet rule - an empirical figure for insulated HF wire near ground
  inside a multi-element array, all of which shorten the resonance. An isolated
  bare-wire loop is not resonant there at all: it reads 107 - 108j at
  b = 8e-4 lambda. Free-space resonance is 1.0516 to 1.1362 lambda, 3.5 to 6.8%
  longer. The number was not invented, it was just carrying a context the spec
  never stated - so the spec now says which is which, keeps the empirical value
  as the nominal, and derives the free-space one beside it.
- 120 ohm is 4 to 9% low; the truth is 124 to 140 ohm. Much closer than the
  circular loop's 100 ohm was.
- 3.3 dBi is essentially EXACT - 3.299 on 1e-4 lambda wire - and only 2.1% low
  on the fattest conductor tested. It now has a test of its own so a later
  refactor cannot quietly move it. After a run of constants that were 30% and
  48% out, it is worth recording that one of them was simply right.

The cross-check that costs nothing and says the most: a square encloses less
area than a circle of the same perimeter, so it must resonate LONGER and show a
LOWER resistance and a LOWER directivity. All three orderings come out that way
across every wire gauge, from two geometries fitted independently.

The modal solver moved out of the scratchpad into `otahub/num/loop_modal.py`,
beside the MoM, because it is a reference model and not a test fixture - the
quad tests anchor to it as well as the circular ones.

Remaining constants in the loop family: `halo_loop` (1.0 dBi, 15 ohm, 1.5 dB
azimuth ripple) and `alford_loop` (0.5 dB ripple). The halo is a bent dipole
with a lumped capacitive gap, which this solver cannot model yet - adding a
lumped load to the EFIE is a small, well-defined piece of work and the obvious
next step.

## The halo, and an engine bug that answered the wrong question

`halo_loop` needed two things the solver did not have: bent open wires and
lumped loading. Both went in, and both were checked before being used. At zero
bend `arc` reproduces `dipole` to every digit, so following the bend round from
straight to almost-closed never leaves ground the dipole already validated
twice. And a huge series load at one point of a ring agrees with a real gap cut
there to four figures on resistance - two unrelated ways of breaking a ring,
and they meet.

The geometry was wrong before any of the physics was. The spec synthesised a
0.47*lambda conductor beside a 0.5*lambda ring and a 0.015*lambda gap, and
those do not add up: 0.47 + 0.015 = 0.485. The gap had been subtracted twice.
Conductor + gap = circumference is now an identity the spec computes and a test
asserts exactly.

The spec's REASONING was consistently sound, which is worth saying because it
has not always been. It said the resistance is about a fifth of a dipole's
because the bent halves' currents partly cancel - and following the bend
continuously shows exactly that, 69 ohm straight falling to 11 ohm almost
closed. It said the polarisation is horizontal, and that turns out to be exact
rather than approximate: the horizon field is 100.000% E_phi at every azimuth.
It said the peak is in the plane of the ring, and the horizon maximum IS the
global maximum to 0.006%.

Its numbers were another matter:

- Resistance is 9.6 to 13.9 ohm, not 15 - 8 to 56% high, and 32% high at a
  typical gap and tubing.
- Peak directivity is 1.21 to 1.27 dBi against an assumed 1.0. Close, and
  inside the 0 to 2 dBi its own note claimed.
- Azimuth ripple is 2.7 to 3.3 dB, roughly DOUBLE the asserted 1.5, so the
  worst-case azimuth gain is -1.4 to -2.1 dBi rather than -0.5. That is 1 to
  1.6 dB optimistic in precisely the number the spec's note says an
  omnidirectional link budget runs on. A halo is less omnidirectional than it
  is usually sold as.

One thing worth knowing that is the opposite of the intuition: a gap capacitor
makes a halo SMALLER, not larger. It end-loads the dipole. A 0.5*lambda ring
cannot be series-resonated by a capacitor opposite the feed at all - the only
zero crossing on that branch is a parallel antiresonance around 25 kilohm. That
took a wrong turn to establish: the first tuning sweep assumed the input
reactance fell as gap capacitance was added, and it rises, from the closed
ring's -3475j towards the open gap's +88j. Scanning instead of assuming
monotonicity found the pole in between.

## An engine bug: the spec's default silently beat the caller

Sweeping the halo's gap did nothing at all, which is how this surfaced. The
resolver ran every synthesis rule unconditionally, so a value the CALLER
supplied was overwritten by the spec's own nominal for it. Ask for a
0.03*lambda gap, get a design sheet for a 0.015*lambda one - and nothing warns,
because every number on that sheet is perfectly self-consistent. It simply
answers a question nobody asked, which is the quietest kind of wrong this
project keeps running into.

A supplied value is a requirement and now outranks the default. The change
touched every archetype in the catalogue and broke nothing: 487/487 known cases
still pass. Two tests in `test_engine_behaviour.py` pin it, including that
overriding one output does not strand the rules downstream of it.

Two of my OWN test claims failed on first run and were the things that were
wrong, not the code: a wide gap on thin wire resonates very slightly ABOVE
0.5*lambda, so "the resonant ring is always smaller than nominal" was an
overreach; and the loaded branch dips below 50 ohm, so a threshold I picked by
eye was simply wrong. Both assertions were narrowed to what is true, and the
spec sentence that overstated the same thing was corrected with them.

The loop family is now fully derived apart from `alford_loop`'s 0.5 dB ripple.

## Bandwidth that is computed, and a 4:1 that needed qualifying

Fractional bandwidth is the metric this catalogue has most often carried as
"indicative", because the closed forms reach it only through an ASSUMED Q. The
solver can do better. Geometry here is in wavelengths, so scaling every
dimension is the same thing as moving the frequency, and a sweep costs nothing
beyond rebuilding the model.

It is computed twice, by routes sharing no algebra: directly, by walking
outwards from resonance until VSWR crosses 2 and then bisecting, and through
the antenna Q of Yaghjian and Best, Q = (w0/2R0)|dZ/dw|. They agree to 2.0% at
worst, and they part company in the right direction - the gap widens as the
bandwidth grows, which is exactly how a narrowband approximation should fail.
Neither would be worth putting in a spec on its own.

`folded_dipole` was the first customer, and the headline is a correction to a
correction. The classic 4:1 step-up holds, and holds tightly - 3.949 to 4.021
across every spacing and gauge - but only between structures AT THEIR OWN
RESONANCE. At a fixed frequency it is not 4 and not even constant: it runs 2.2
to 4.9 across a 20% span of length, because the transmission-line mode does not
scale with the radiating one. My own first reading of the sweep took the
same-frequency ratio for the law and briefly had the spec wrong; the two
comparisons are different questions and the spec now says which is which.

What else was found:

- Resonant length is 0.449 to 0.474 lambda, not the flat 0.48 synthesised.
  Wider spacing shortens it markedly.
- Resistance is 284 to 290 ohm, where the spec computed 4 x 73.079 = 292.3.
  The 1 to 3% overshoot is the induced-EMF 73.079 standing in for a real
  dipole's 71.9 ohm resonant resistance - the same confusion that cost an
  afternoon when the MoM first disagreed with the textbook.
- Bandwidth is 10.7% to 18.7%, against a flat 0.10. And it is NOT "roughly
  twice a plain thin dipole's" as the spec said: against a plain dipole of the
  same gauge it is 1.36 to 1.72 times, never 2.
- The pattern is not quite "indistinguishable from a plain dipole" either. The
  two conductors are a short end-fire pair, so the directivity is 1.1% up on
  the closest spacing and 7.0% up on the widest - 0.05 to 0.29 dB, small but
  systematic rather than noise.

Two cross-consistency tests failed on the rewrite, and both were encoding the
error rather than catching it: they compared the folded dipole's resistance
against the induced-EMF 73.079 and demanded exactly 4.00 and 9.00. Against that
constant the real ratios are 3.91 and 8.81. One now asserts 4 within 3% AND
that it lands below 4, with the reason; the other compares each archetype's own
N-scaling instead of dividing a derived number by an asserted one, which is
what it had been doing.

## A low-confidence archetype promoted, for the first time

`halo_loop` has come off the low-confidence list. That is section 6 item 4 of
the handover discharged - "validate a low-confidence archetype end to end and
either promote it or record why it cannot be" - which had been marked
impossible since session 1 for want of a solver.

The remaining eight cannot follow it, and the reason is the same for nearly all
of them: `cassegrain`, `conical_horn_dual_mode`, `pifa`,
`planar_monopole_circular`, `stacked_patch`, `vivaldi_tsa` and
`waveguide_longitudinal_slot` are patches, horns, reflectors and slots, and a
thin-wire solver cannot reach any of them. `ferrite_rod_loop` needs a ferrite
material model besides. That is a limit worth stating plainly rather than
leaving as an open task that looks actionable and is not.

## Four asserted Q values, three of them derived and one explained

A sweep for bandwidth metrics found four reaching their answer through
`1.0/(Q_meas*sqrt(2))` with Q simply asserted, all in the loop family. With
last round's machinery three of them are solvable.

The circular loop got three routes, because it is the one shape here with a
modal solution: the MoM's direct VSWR walk, the MoM's Q-derivative, and the
Fourier-mode solver's own walk. They agree to 0.38% between the two SOLVERS and
0.76% including the Q route. That is worth more than any single number.

- `one_wavelength_circular_loop`: 7.6% to 16.7%, against 5.89% from an assumed
  Q of 12. Real Q is 9.3 down to 4.3, so the old figure was up to 65% low.
- `quad_loop_square`: 7.0% to 15.6%, against 5.05% from an assumed Q of 14.
  Real Q 10.2 down to 4.6, up to 67% low.
- `halo_loop`: 1.70% to 2.98%, against 2.83% from an assumed Q of 25. This one
  erred the OTHER way - the halo is genuinely higher-Q than 25 over most of its
  range, up to 42, so the assumption promised bandwidth the antenna does not
  have. Worth noting that two assumptions in the same family, written by the
  same hand, were wrong in opposite directions; there was no systematic bias to
  correct for, just three guesses.

`Q_meas` is gone as an input from all three. A number the solver can produce
should not be something the user is asked to supply.

## The one that could not be derived, and why that is the interesting one

`alford_loop` keeps its assumed Q, and its notes now say so in as many words.
Its bandwidth is set by the folded, capacitively loaded corner sections that
force the current uniform - and that network is not part of the geometry this
spec describes. Fitting a number to the bare square would have looked exactly
like the other three and meant nothing.

What CAN be checked about it was, and both results were worth having:

- The equal-area circular equivalence its whole model rests on is far better
  than the "good to a few percent" its own validity note claimed. Against a
  uniform-current solve of the actual square it is within 0.73% on radiation
  resistance and 0.08% on directivity over perimeters 0.5 to 1.0 lambda,
  worsening monotonically with size exactly as it should.
- The asserted 0.5 dB azimuth ripple is not a property of the shape at all.
  With the current actually uniform, a square loop's horizon pattern is flat to
  0.001 dB at half a wavelength and 0.010 dB at a full one - fifty times
  smaller. So the ENTIRE ripple budget is how well the corner loading does its
  job. That is a more useful thing for a builder to know than the 0.5 was, and
  it only turned up because the ideal case was computed rather than assumed to
  be roughly the real one.

## The reference dipoles: which numbers are definitions and which are the wire

`half_wave_dipole` and `resonant_dipole` are what everything else is measured
against, and both were built entirely on Balanis 4-70 and 4-79. Those are
exact - for an ASSUMED sinusoidal current on a vanishingly thin wire - and this
build had already found, twice, that such a current is not what a delta-gap-fed
wire of finite radius carries. This round took the lesson back to the source.

The decisive check was a one-liner. Hand the MoM the sinusoidal current at the
resonant dipole's own length and it reproduces the spec's resistance to 0.3%.
So the disagreement with the driving point is the real current's shape and
nothing else - not mesh, not quadrature, not the gap model.

What was found:

- `resonant_dipole` presented the induced-EMF radiation resistance as what the
  antenna presents: "55-68 ohm across practical thicknesses". A real resonant
  wire sits at about 72 ohm and stays there - within about an ohm of 72.5 from
  aw = 1e-5 to 2.7e-3 wavelengths. The closed form is 6% low on the thinnest
  wire and 22% low on the fattest, because it falls away with thickness and the
  wire does not. The folk figure "a resonant dipole is about 72 ohm" was right
  all along; the formula quoted in its place was not.
- Its `vswr_in_50_ohm` inherited the error: about 1.26:1 at aw = 1e-3 where a
  real one sits at 1.46:1. Optimistic in the first number a builder checks.
- Its reactance was "zero by construction", which was true of the model: the
  length fit was tuned to the closed-form reactance. A real wire at that length
  shows -3 to +6 ohm. Recorded, and bounded in a test, but deliberately NOT
  fitted - a quantity that crosses zero has no meaningful relative error.
- It borrowed the half-wave directivity 1.6409. The shortened dipole's own is
  1.636 to 1.638. Its note had the physics right and the number from the wrong
  antenna.
- Bandwidth is now computed: 5.8% on 1e-5 lambda wire to 17.4% on 5e-3, a
  factor of three - the quantitative form of "fat dipoles are broadband".

`half_wave_dipole` was handled differently, on purpose. Its 73.08 + j42.52 is
the definition dBd rests on and half the field quotes, and replacing it would
break that. So it stays exactly as it was, now named for what it is, and a test
pins it. The driving-point impedance sits beside it: 78.1 + j44.3 ohm on
1e-5 lambda wire rising to 94.7 + j45.8 on 2.7e-3, tending to 73.08 only as the
wire vanishes - and slowly, logarithmically.

One correction to my own work before committing. The notes first said the
driving point was "7 to 22% higher" - a range copied from the resonance
comparison. Evaluated at the spec's actual length it is 6 to 28% higher, which
is also to say the closed form is 6 to 22% LOW: the two percentages measure
the same gap from opposite ends, and I had conflated them. Caught by printing
what the spec produces before writing the tests, rather than after.

## A ten-fold error in the loaded whip, from one missing division

The loaded monopole's entire design chain - reactance, then the coil that
cancels it, then that coil's loss, then the efficiency - hangs on one number:
the reactance of a short wire. Image theory lets the ground-plane-free solver
reach it exactly, since a monopole on an infinite PEC plane has half a
dipole's impedance.

At h = 0.05 lambda the spec said -j92 ohm. The method of moments said -j931,
Hallen's equation -j920, and Schelkunoff's transmission-line model - cruder -
-j1091. Three independent routes an order of magnitude from the spec.

The cause was diagnosable exactly, not just detectable. Balanis 4-79 is
referred to the current MAXIMUM of the assumed sinusoid, and on a short wire
that maximum lies far outside the antenna. The feed-point reactance is the
current-maximum one divided by sin^2(kL/2), which at L = 0.1 lambda is 0.0955 -
a factor of 10.5. Dividing the spec's value by it lands within 0.7 to 5.7% of
the MoM at every length and radius tried.

And the omission had a paper trail. `dipole_arbitrary_length` already carried
`input_resistance_ohm = radiation_resistance_current_max_ohm / sin(k0*L/2)**2`
- correctly referred - with `input_reactance_ohm` on the very next line NOT
referred. Its R and X described two different points on the same antenna. The
loaded monopole then copied the reactance, and even its own known case recorded
the feed-referred RESISTANCE (0.9994 ohm) as a cross-check while asserting the
unreferred reactance beside it. The half-wave dipole and quarter-wave monopole
were never affected, because at exactly lambda/2 the referral is the identity -
which is presumably how it survived. A test now asserts that R and X carry the
identical 1/sin^2 factor, so the two cannot drift apart again.

The consequence for a builder, 1.5 m whip at 10 MHz:

- reactance -j1163 ohm, not -j111 (the MoM says -j1140, within 2%);
- loading coil 18.5 uH, not 1.77 - wound to the old figure, the whip would
  still be strongly capacitive, a hundred half-bandwidths off resonance;
- coil loss 5.81 ohm, 5.9 times the radiation resistance, not the 0.555 ohm
  "over half the radiation resistance" the case used to teach;
- radiation efficiency 14% (-8.4 dB), not 61% (-2.1 dB). The spec overstated
  a short loaded whip's efficiency by a factor of 4.3.

`turnstile_dipole` had the same omission for both R and X; at its default
0.4788 lambda it moves the impedance only 0.45%, but it would not have stayed
small at any other length, so it was fixed too.

## Knowing where the solver stops

Checking the whip's RESISTANCE turned up a limit of the solver rather than of
the spec, and it is recorded in `mom.py` because the next person needs it. On a
short wire the delta gap puts a local current excess on the feed node - about
8% at 0.1 lambda, 17% at 0.04 - over the smooth distribution, and it grows as
the mesh resolves it. That is the whole of the resistance drift seen under
refinement: referencing the radiated power to the smooth current instead gives
a mesh-independent answer. But that smooth answer lands 6-8% ABOVE the
triangular closed form rather than on it - which could be the genuine
finite-radius correction or could be the extrapolation. So the solver
brackets short-wire resistance at roughly -10% / +7% and cannot settle it.

The closed form sits inside the bracket and agrees with the sinusoidal-current
input resistance to 1.3%, so it stays, and the spec says why. The reactance,
which the gap barely touches, was the 10x error and is unambiguous by every
method. Claiming the resistance was wrong as well would have been easy and
unsupported.

Hallen's solver moved into `otahub/num/hallen.py`, beside the MoM and the modal
loop solver, now that two test files arbitrate with it.

## The helix: a correction worse than the formula it corrected

The previous round's lesson was that a textbook formula can be badly wrong
outside the conditions it was derived for, while still looking textbook-correct.
The travelling-wave family carries two formulas the literature has long called
optimistic, Kraus's helix gain and Carrel's LPDA curves, so they were next. The
LPDA needs a transmission-line feeder network the solver does not have yet. The
helix did not.

The spec already knew Kraus overestimates, and said so plainly. It then told
designers to "use gain_corrected_dbi for design": 8.3 + 10 log10(C^2 N S) +
20 log10(1 + N/10), which its note described as "a lower constant with a
sub-linear term in N". The added term is super-linear. The formula exceeds
Kraus for every N >= 5, by 2.5 dB at ten turns, 6 at twenty and 8.5 at thirty,
where it rates a 30-turn helix at 28.75 dBi, and all of that is in exactly the
regime where the same note says Kraus is already too high. No solver was needed
for that, only arithmetic on the spec's own words. The metric is gone rather
than relabelled: the spec had been recommending it for design.

What the right answer IS needed the solver. A helix over an infinite PEC ground
plane is one continuous wire by image theory: image helix, short vertical feed
straddling the plane, real helix. Current continuity through the bends produces
the image-current rules with nothing special, and the far field comes out
mirror-symmetric to 1e-7, which is the check that the construction is right.
The directivity was then found a second way that never touches the pattern
integral: illuminate with a circularly polarised plane wave, take the
open-circuit voltage, and use reciprocity. The two agree to 0.02 dB at ten
turns. Power balance holds to 1e-4.

What was found:

- Kraus's error grows with BOTH length and circumference. At C = lambda and 13
  degrees: 1.6 dB high at 3 turns, 3.8 at 10, 4.4 at 20. Across the design core
  it is 0.9 to 5.4 dB high from five turns up, and roughly right only in the
  short, small corner (3 turns at C = 0.9 lambda).
- Pitch runs the OPPOSITE way. Steeper pitch lengthens the helix, so Kraus
  rises. But it moves the phase velocity off the end-fire optimum, so the real
  gain falls, by about 0.5 dB per degree.
- At the band edge long helices fall off a cliff: 11.8 dBi at 20 turns and
  C = 1.2 lambda, where Kraus says 20.5. No polynomial here follows that cliff,
  so the fit covers only the smooth design core (C 0.9 to 1.1) and the whole
  135-point grid ships in the spec's tables instead, so the edge stays visible.
- Kraus's beamwidth is too narrow for the same reason: 34 degrees for the
  canonical helix where the solver finds 46.
- Input resistance and axial ratio are NOT fitted, deliberately. Feed height
  barely moves the gain (0.1 dB) but swings the terminal resistance from 70 to
  150 ohm, and the axial ratio varies irregularly with the standing wave off
  the open end. Both remain indicative, now saying why.

Two corrections to my own work before committing, both of the same kind. The
first draft of the Kraus note read "-0.1 to 5.4 dB HIGH", a range taken
blindly from a min and max over the whole core, and it contradicted itself.
The second draft claimed Kraus was "about right for a three-turn helix", which
holds only at C = 0.9 lambda: at C = lambda a 3-turn helix is already 1.6 dB
high. Both were caught by a test that asserted the claim. The note now quotes
the grid directly rather than paraphrasing it. Three rounds running, a range
written into a note has needed checking against the numbers it summarises,
which is a good enough reason to keep doing it.

Wire size matters more than expected: about 1 dB lower at a = 0.002 lambda and
0.8 dB higher at 0.01, at ten turns, by moving the phase velocity. The fit is at
0.005 lambda and the validity block says so.

## Feed networks, and an LPDA whose number held up

The LPDA was left over from the helix round: its directivity is a reading of
Carrel's 1961 chart, a chart the literature has long called optimistic, and it
could not be checked because the solver had no way to model the feeder line
joining the elements. That is now in `mom.py`: a non-radiating network coupled
to the wires in admittance form, the way NEC's TL cards work. At every port the
antenna draws Y_ant V, with Y_ant taken from the inverse impedance matrix; the
network adds its own admittance; Kirchhoff's current law closes the system.

Three checks licensed it before it touched an LPDA, each against something that
shares no code with the coupling:

- a dipole behind a length of line reproduces the textbook impedance
  transformation to machine precision;
- two dipoles tied in parallel through a stiff network match the plain solver
  driving both gaps at once, to 1e-6;
- a lossless feeder delivers exactly the power the pattern integral says is
  radiated, to 1e-5 - which would fail at once if the coupling were wrong.

Then the LPDA, with every element solved and the feeder as a transposed line.
An LPDA repeats in frequency with period tau, so a single frequency samples one
phase of a ripple; each tau was scanned over a full period, with the active
region held seven elements from the back and twelve from the front.

**The spec's number held.** I went in expecting Carrel's reported optimism and
did not find it: along the optimum-sigma line the spec's reading is within
0.6 dB of the solved log-period mean at every tau, from 0.57 dB high at 0.8 to
0.55 dB low at 0.95. Like the quad loop's 3.3 dBi, it is worth recording
plainly that a number was right. What the single number hides is the ripple -
2.0 dB peak-to-peak at tau = 0.8, 0.7 at 0.9 - so the solved mean, minimum and
maximum now ship as a table.

**What the spec was missing entirely** was the input resistance, and with it
the feeder impedance, which is how an LPDA is matched. Directivity also moves
with the feeder - 0.7 dB between 100 and 200 ohm, 0.9 dB across 50 to 200 - so
it was never really a function of tau alone. The standard relation R = Z0/sqrt(1 + Z0/(4 sigma' Za))
is added, with its inverse for the feeder that gives 50 ohm, and checked
against the solver: -6% to +11% over feeders of 50 to 200 ohm and element
ratios of 60 to 250.

**One claim in the spec was wrong.** Its validity said that without the feeder
transposition "the array fires backwards". Solved over a full period, it does
not do that reliably. It loses its front-to-back ratio altogether - within a
few dB either way - and its input impedance swings tenfold inside one period,
37 to 479 ohm at tau = 0.9, against 12 to 33 dB front-to-back and 72 to 99 ohm
when transposed. So the transposition is what makes it frequency-independent,
which is the whole point of the design. The test I first wrote asserted the
spec's version, and failed; it now asserts what the solver shows, as the spec
does.

Two claims of my own were checked against the data before they went in, as
has become routine: "reactance under about 15 ohm" was 17.5 at worst, and is
quoted as that. The tau = 0.95 point is the least converged - a high-tau array
has a wide active region, and enlarging the solved structure moved it 0.6 dB -
and the validity block says so.

## The rhombic, and the parameter it was missing

The rhombic's and V's directivities were marked "derived here" in session 5,
from an idealised model: an unattenuated travelling wave on every leg and a
perfect termination. Of the travelling-wave family, only the rhombic can be
solved faithfully in free space. A V's and a long wire's terminations need a
return path to ground that free space does not have, but a rhombic closes on
itself: four legs, the feed at one acute vertex, the resistor at the other. It
went in as one closed wire with the resistor as a lumped load, and power into
the resistor balances to 1e-4.

Mesh mattered more here than anywhere before. At 12 segments per wavelength
the axial directivity read 0.35 dB low, which first made the idealised model
look 0.6 to 1.2 dB optimistic. At 28 per wavelength it is within about 0.05 dB
of converged, and the real gap is much smaller:

- The idealised model is 0.24 to 0.45 dB HIGH across leg lengths 2 to 12
  lambda, growing slowly with length, because the real current decays as it
  radiates and the model's does not. Modest, but systematic.

The bigger finding was a parameter the spec did not have. The wire radius sets
three things while barely touching the directivity (0.33 dB across 2.5
decades):

- The optimal termination - the resistance that maximises front-to-back - runs
  from 878 ohm at a = 1e-5 lambda to 273 at 3e-3, linear in ln(lambda/a) at
  about 106 ohm per neper, which is how a two-wire line's impedance goes. It
  hardly depends on leg length (596 to 635 ohm from 2 to 12 lambda). The flat
  600 ohm the spec carried is right to 5% at a = 1e-4, typical HF wire, but 32%
  low on the thinnest wire and 55% high at 1e-3.
- The share of power the termination burns runs from 26% on fat wire to 56% on
  the thinnest. The old "about half" is right for HF wire and wrong for fat
  wire, and the gain sits 1.3 to 3.5 dB below the directivity, not a flat 3.
- `aw` is now a parameter, and the termination, radiation efficiency and gain
  follow it.

And the leg angle: the alignment angle, where each leg's cone lines up with the
axis, is right for an unattenuated wave. With the real current the axial
directivity peaks near 27 degrees at 4 lambda, 0.34 dB above the alignment
angle's 24.9 - a WIDER diamond. The V antenna's spec says its own optimum is "a
few degrees tighter"; that could not be checked, since the V cannot be solved
in free space, and the rhombic points the other way. Worth solving the V over
ground one day before trusting either claim.

The idealised directivity is kept as `directivity_travelling_wave_model`,
because the V still uses the same model and `gain_over_v_antenna_db` should
compare like with like. The old known cases asserted the idealised figure under
`directivity_linear` and now assert it under its own name. Two of those had
been passing only because a relative tolerance on a dB value is loose: 15.93
against 15.63 dBi is under 2%. They were moved to the right metric rather than
left passing by accident.

A first draft of the termination note said 600 ohm was "right to within about
7% for thin wire". That was true only at a = 1e-4; on the thinnest wire it is
32% low. The note now computes those percentages from the solved data rather
than paraphrasing them - the same correction as the last three rounds.

The suite has grown to 105 seconds, from 77, mostly the rhombic's fine-mesh
checks. Worth marking the heaviest MoM tests slow if it keeps climbing.

## Junctions, so wires that touch are connected

Until now two wires that merely touched were not connected at all: basis
functions lived inside single wires, and a node shared by several wire ends
carried nothing. That put a whole class of real antennas out of reach - a top
hat of radial wires, a ground plane of radials, a discone built as a wire cage -
and `top_loaded_monopole` has carried a documented honesty problem since
session 5 for exactly that reason: its hat radius and its top current are
independent inputs, because nothing could derive one from the other.

The change was to the core model, so it was done in the order that keeps it
safe. First the basis layout was generalised: each function now stores its two
halves with their orientation, so a wire may END at a node or START there and
still carry current through it. A half whose current runs against its
segment's direction simply carries negative weights, and its divergence comes
out the same either way. Before touching junctions at all, seven reference
impedances were recorded from the old code - dipole, loop, halo, helix, rhombic,
folded dipole, LPDA - and the refactored code reproduces every one of them with
a relative difference of exactly zero.

Then the junctions: wherever K open-wire ends coincide, K - 1 junction
functions carry current from the first wire into each of the others, so
Kirchhoff's current law holds by construction. A feed at a junction drives all
of them together, and its terminal current is their sum.

Checked, each against something outside the new code:

- a dipole split at its centre is the same dipole to 1e-14, in all four
  combinations of which way each arm runs - the orientation logic in isolation;
- a bent wire split into three pieces, the middle one reversed, matches the
  single polyline to 7e-15;
- a symmetric T splits its current equally between its arms to 4e-15, balances
  power to 1e-5, and keeps the impedance matrix exactly symmetric;
- a ground-plane vertical on four quarter-wave radials, fed at its five-wire
  junction, reproduces the two familiar figures - 22.5 ohm with the radials
  flat, 52.2 with them drooped 45 degrees, against the usual ~22 and ~50.

The last one carried a small lesson worth keeping. The first comparison was at
exactly a quarter wavelength, which gave 66 ohm at 45 degrees and looked wrong.
It was not wrong, it was not resonant: that antenna carried +41j of reactance,
and the rule of thumb is about resonant antennas. Resonating each droop first
gave 52.2. A comparison with a published or remembered figure is only as good
as the match between the conditions it was quoted under and the ones computed.

The capability is in; putting it to work on `top_loaded_monopole`, the discone
and the conical monopole is the next round.

## The top hat: a documented gap closed, and a fit that was declined

`top_loaded_monopole` had said since session 5, in its own validity block,
that "beta_top is an input, not a prediction - relating hat size to beta_top
needs a numerical solve this spec does not attempt". With junctions it could
be attempted: a vertical with a hat of radial wires, joined where they meet,
over infinite ground by image theory.

**The model was right; its input was not.** The spec's radiation resistance,
160 pi^2 (h(1+beta)/2)^2, is within 0.1 to 2.4% of the solved radiated power
for every hat tried - provided beta is the true one. And the spec's defaults
disagreed with each other: a 0.01 lambda hat with eight radials gives beta =
0.52, not the 0.6 carried beside it, so the default radiation resistance was
11% high. beta is now derived from the hat and can still be overridden; every
existing case that set it explicitly passes unchanged, which is the engine fix
from the halo round earning its keep.

Two things about the fit are worth recording.

- **A variable was missing.** The first fit used three dimensionless groups -
  hat radius over height, height over wire radius, radial count - and could not
  get below about 0.1 in beta, whatever form or restriction was tried, including
  a physics-shaped predictor with the hat's coupling to ground. There are four
  groups, not three: the ELECTRICAL height h/lambda had been left out, on the
  quiet assumption that a short antenna is quasi-static. At 0.05 to 0.1 lambda
  it is not. Adding it took the worst error from 0.1 to 0.017, and 0.027 on
  held-out cases.
- **Radial count matters far more than expected.** On a 0.01 lambda hat, beta
  is 0.41 with four radials, 0.61 with sixteen and 0.67 with thirty-two, still
  rising. A solid disc beats any radial count the fit covers, and the spec says
  so; radial count is now a parameter.

The capture of the base current needed care: on a wire this short the delta gap
over-reads the feed node, so the base current is taken from a straight-line fit
to the vertical's current away from the feed. That definition is also the one
that reproduces the radiation resistance to 2.4%, which is the reason to trust
it.

**The reactance was declined, deliberately.** It is the number a builder needs
most - it sizes any remaining loading coil - and the spec has never computed
it. Three forms were tried: log reactance (fails where it crosses zero at
resonance, 480% error), a transmission-line angle X = -Z_v cot(theta) that
passes smoothly through resonance (17% held-out), and reactance relative to the
bare whip as a function of beta (0.2 spread - radials change capacitance and
current taper differently). 17% on a high-Q short antenna would size a coil
that misses resonance by several bandwidths, so no formula went in. The solved
values ship as a table instead, including the hats large enough to resonate the
whip by themselves, which the beta fit also excludes because the linear taper
stops applying there. Declining a fit that looks respectable is the right call
when the error would land in exactly the quantity the user acts on.

The radial meshing needed one fix first: tying radial segments to the
vertical's segment length gave a wide hat on a short whip several thousand
segments. Capping at 12 per radial keeps the length ratio across the junction
at 5 or less and moved beta by 0.0002.

## The cones: an angle bug hidden by its own cross-checks

The plan was to use the new junctions on the cone family - a biconical, and a
conical monopole as half of one by image theory - built as wire cages: two
cones of radial wires joined at their apexes to a short feed wire. Reading the
existing cross-check between the two first turned up something better than the
directivity.

`test_conical_monopole_is_half_a_biconical` compared a biconical at theta_h =
10 degrees with a conical monopole at 5. The discone check paired a 30 degree
discone with a 60 degree biconical. The reason: `biconical` documents theta_h as
the HALF angle ("half angle measured from the axis; 30 deg = 0.5236 rad") and
computes Zc = (eta0/pi) ln cot(theta_h/4) - Kraus's formula for the FULL cone
angle. So it reported the impedance of a cone half as wide as the one described:
243 ohm for a 30 degree cone that presents 158, and 188 for the classic 47
degree cone - the one that gives a 50 ohm monopole and a 100 ohm bicone - that
presents 100. Its siblings used the half angle correctly. The two
cross-consistency tests existed to catch exactly this, and had been written
with a doubled angle on one side so that they passed. They now compare cones of
the same half angle, and the formula is fixed. The solver agrees on the
direction: a 47 degree cage's resistance oscillates between 82 and 177 ohm over
the band sampled, mean 128 - nowhere near 188.

That impedance is only an ordering, though, because a wire cage's impedance does
not converge at affordable wire counts: it was still moving 3 to 5% between 24
and 32 wires, and 32 already costs minutes. Directivity does converge - 0.2%
between 16 and 24 wires and between two segment lengths - so directivity is
what the cages were used for:

- `biconical` carried the half-wave dipole's 1.6409 as a flat constant. At the
  band edge (slant = lambda/4) a narrow cone slightly exceeds it (1.70 at 5
  degrees) and a wide one falls well below (1.45 at 47, 1.33 at 65): up to 23%
  high. Its note had the right idea - "close to a dipole for narrow cones" -
  and now has numbers. Fitted over 5 to 65 degrees to 0.05%.
- `conical_monopole` carried 3.0, which is the SHORT monopole's value. A
  quarter-wave cone is 3.39 when narrow and 2.66 when wide; at the default 47
  degrees it is 2.90, close to 3.0 by coincidence. It is now exactly twice the
  biconical, as image theory requires, and the cross-check asserts that too.
- Across the band the two flares go opposite ways: a 30 degree cone GAINS
  directivity up to 0.6 lambda slant, a 47 degree cone LOSES it (2.47 over
  ground by twice the band-edge frequency) as its beam lifts off the horizon.
  Both specs carry the solved values as a table rather than pretend to a
  single number.

The bandwidth ratios stay indicative, now with the reason stated: a thin-wire
cage cannot represent a solid cone at the top of a decade, where the gap
between wires stops being small against the wavelength.

The suite now takes 132 seconds, from 107 two rounds ago, and the growth is
all solver tests. The next housekeeping item is a slow marker on the heaviest
of them, so the default run stays quick.

## Auditing the cross-checks, and a limit of the solver found on the way

The cone round found two cross-checks that had been written to agree rather
than to check. That warranted reading all of `test_cross_consistency.py` with
one question: does each relationship compare the same physical antenna on both
sides? Most do - same frequency, same geometry, or a reduction within one
archetype. Three were worth a closer look.

- The slot-array check asserts equality with a single slot at a hand-typed
  3 mm offset where the array derives 3.06 mm. Sound after all: both sides
  compute the same Stevenson coefficient, which depends on the guide and the
  frequency and not on the offset.
- The corner reflector's self-resistance check compares like with like.
- **The Babinet check was a tautology.** "Slot times its complementary dipole
  is eta squared over four" - but the slot spec DEFINES its impedance as
  eta0^2 / (4 Z_dipole) from dipole constants it carries itself, so the
  product is eta0^2/4 by construction. It checked nothing about Babinet, and
  Babinet cannot be checked here at all, since the solver models wires, not a
  slot in a screen. What it could check - that the slot's borrowed constants
  agree with what the dipole spec computes - it now checks, under a name that
  says so.

Following that thread turned up a real inconsistency. Babinet is exact, so a
slot resonates where its complementary dipole does - and that dipole's radius
is a quarter of the slot width, 0.006 lambda at the default. `half_wave_slot`
uses 0.4785 lambda, which is a THIN dipole's resonant length; the catalogue's
own `resonant_dipole` fit at that radius says 0.463, and the solver says 0.465.
The slot's borrowed resonant resistance is thin-wire too: about 74 ohm would put
the resonant slot near 480 rather than 530. It is NOT fixed yet, and says so:
the dipole fits stop at 5e-3 lambda, and extending them is the next round's job.
Until then a strict expected failure holds the gap open in the cross-checks and
the spec's validity block names it.

**The solver limit.** Checking that complementary dipole directly, the answer
would not converge - 79 ohm at 60 segments, 113 at 100 - and the reason is the
reduced thin-wire kernel. Mapped: on a 0.006 lambda wire the resonant
resistance holds at 72-74 ohm while segments are 3 to 8 radii long, drifts to
79 at 1.3 radii, and reads 115 - with the resonance moved 6% - at 0.8 radii. A
finer mesh on a fat wire is not more accurate; it leaves the approximation. The
rule, keep segments at least about three radii long, is now in the module
docstring, and WireModel warns when any segment is shorter than its own radius.
The whole suite runs clean with that warning promoted to an error, so nothing
already built was relying on the broken regime.

Housekeeping, overdue: the suite had reached 133 seconds, all of the growth
solver applications. Those ten modules now carry a `slow` marker, so
`pytest -m "not slow"` runs the other 1862 tests in about 25 seconds while the
default run still includes everything. The solver's own licensing tests - the
MoM validation, junctions, reference dipoles, feed referral - stay in the quick
set, because everything else rests on them.

## The slot's resonance, and the kernel it needed

Last round left `half_wave_slot` held open: it resonated at a thin dipole's
0.4785 lambda whatever its width, when Babinet says it resonates where its
complementary dipole does - a wire of radius w/4, 0.0058 lambda at the default
w/L = 0.05. Closing it needed that fat dipole solved honestly first, and the
reduced kernel could not do it: its answer depended on the mesh.

**The exact kernel.** `mom` now carries it as an option, `exact=True`: current
and observer both spread round the circumference, R(phi) = sqrt(v^2 + 4a^2
sin^2(phi/2)) averaged over phi, applied to the self and near collinear pairs
where the reduced kernel goes wrong. The log singularity at phi = 0 is removed
by phi = pi t^2; 24 nodes are converged to 1 part in 10^6 (74.9224 against
74.9225 ohm at 48). On the 0.006 lambda wire the resonance now holds from 3.7
radii per segment down to 0.8 - reactance within 0.6 ohm at a fixed length,
where the reduced kernel falls 16 ohm - and the remaining 3% creep in
resistance is the delta gap's, not the kernel's. On thin wire the two kernels
differ by 0.35%. It is OFF by default, and the default path is pinned bit for
bit, so nothing derived with the reduced kernel moves. What it does not cure:
at a = 0.015 lambda the delta gap itself diverges with either kernel. That is
the next solver limit, and it is recorded as a next step, a finite-gap feed.

**The slot, fixed.** Length and resonant resistance now come from
`resonant_dipole`'s laws at the complementary radius, and `Rd_res` becomes
derived (supplying it still overrides, and a known case pins that the old 67
reproduces the old 529.6). At the default width: 0.4637 lambda and a 468 ohm
slot, against the old 0.4785 and 529.6 - the old resistance was 13% high. The
laws were fitted only to 5e-3 lambda, so they were checked past that against
two solvers that share nothing, each at meshes it has converged on: exact-kernel
MoM (42 segments up; 22 is still 0.1% out) and Hallen (two radii per segment
or more, which its own docstring sets). They give 0.4645-0.4654 lambda and a
471-479 ohm slot; the law is 0.2-0.4% short on length and 0.6-2.5% high on
resistance. At w/L = 0.02 it is inside 0.1% on length. An honest first pass
claimed "within 0.4%" off a coarse mesh; the convergence scan put the worst
case at 0.44% on that mesh and 0.38% on converged ones, and the spec text was
rewritten from the scan, not the first pass.

What the old constants hid was the DIRECTION: a wider slot has a fatter
complement, which resonates shorter at MORE resistance, so the slot - its
inverse - sits LOWER. The spec now shows 482 ohm at w/L = 0.02 and 468 at 0.05,
and a cross-check asserts that ordering. `folded_slot` follows the plain slot,
so its two-slot figure drops from 132 to 117 ohm and three slots reach 52.

The strict expected failure is now an ordinary pass at two widths, on length
and resistance both. Past w/L = 0.05 the spec says it is extrapolating: there
the complementary wire is too fat for a delta gap to be solved honestly, and
Hallen already drifts at 0.05 under a fine mesh.

## The CP patch that would have radiated linear polarisation

HANDOVER section 6 says to pick the next audit by coverage - quantities a spec's
known cases assert, over quantities it produces. The bottom of that ranking was
`truncated_corner_cp_patch`, asserting 3 of 16. Nothing in the repo checked
the rest, and a single-feed CP patch is unforgiving: its whole axial-ratio band
is a fraction of a percent wide.

The arbiter is new: `otahub/num/patch_cavity.py`, a cavity model of the actual
truncated outline. Neumann modes by linear finite elements - the mesh follows
the cut exactly and is split symmetrically, so the uncut square stays exactly
degenerate and any split is the cut's - a probe feed, one loss Q for every
mode, and the broadside polarisation summed over 30 modes. It assumes none of
the rules it was brought in to check. Converged to 0.004% between 240 and 480
cells per side.

What it found, in order of damage:

- **The feed rule was the nearly-square patch's.** The spec said feed on a
  diagonal. Cutting a square's corners splits it into DIAGONAL modes, and on a
  diagonal one of them is identically zero: the cavity model gives about
  3000 dB of axial ratio there. Linear polarisation. The feed goes on a
  centreline, and the other centreline reverses the handedness.
- **The square was sized for the wrong frequency.** Only one mode moves when
  the corners are cut, so the CP centre sits about 0.5/Q above the uncut
  square's resonance. The spec sized the square for f0; in its own 2.4 GHz
  example that put the CP centre 0.79% high, outside a +-0.29% axial-ratio
  band - 7.9 dB of axial ratio at f0. The square is now sized for f_sq, and the
  cavity model gives 0.006 dB at f0.
- **The mode split was halved.** f0/(2Q) is each mode's offset from the
  centre; the split is f0/Q. The cavity model puts perfect CP at Q = 1/split to
  within 0.5%.
- **The classic cut is first-order.** dS/S = 1/(2Q) comes from split = 2 dS/S;
  the solved ratio falls to 1.93 by c/L = 0.14, so the classic cut is 0.2% short
  in c at Q = 1150 and 2.6% short at Q = 12.
- **Both bandwidths were one mode's.** Two staggered modes have
  |Gamma| = x^2/(4 + x^2), a VSWR-2 band of sqrt(2)/Q - twice one mode's - and
  a 3 dB axial-ratio band of exactly 0.347107/Q, where the spec had an
  "indicative" 0.35 x one mode's band, 29% low. The cavity model reproduces
  both to 4 digits.

The cut and the offset are fitted in u = sqrt(1/(2Q)) over Q = 12..1156 and
tested on three meshes the fit never saw (0.03%). Known cases come from a
separate script: the patch formulas re-implemented, the corrections read by
spline off the 480-cell scan rather than from the 240-cell polynomial the spec
carries. Q0 is now a derived parameter that can be supplied, because the
thin-substrate formula only knows radiation loss, and a lossier patch needs a
bigger cut.

In passing: `otahub.num.__all__` promised `bicone_cage` without importing it, so
`from otahub.num import *` raised. Fixed. And two entries under "known open
discrepancies" in the handover were stale - the helix gain and the LPDA
directivity were both solved in earlier rounds - and now say so.

## A horn that could not be fitted to its waveguide

Next on the coverage ranking was `pyramidal_horn`, 4 of 18 quantities asserted.
Its efficiency had been rebuilt on exact Fresnel theory two sessions ago, so the
question was what else it claimed. The answer: a1, b1 and ONE apex distance,
and no feed waveguide anywhere in the spec.

A pyramidal horn is buildable only if its E-plane and H-plane flares reach the
feed guide at the same axial length. With one apex for both planes that needs
a guide shaped like the aperture, a/b = sqrt(3/2); every standard guide is
about 2:1. On WR-90 at 20 dBi the old design's two flares met the guide 15 mm
apart - 172.8 mm and 157.7 mm behind the aperture. Two smaller errors rode
along: the note called the single distance a slant length while every formula
used it as the axial one, and the validity blamed the single apex for a
0.037 dB gain offset that was really the 0.51 sizing efficiency.

The textbook fix (Balanis's design procedure) was checked before adopting it,
and not adopted. It writes the optimum proportions in slant lengths and sizes
with an implied efficiency of 0.5105; solved exactly and integrated over the
aperture, its horns come in 0.04-0.19 dB short at 20-25 dBi and 0.8 dB short at
15. Written instead in the AXIAL apex distances the aperture phase depends on,
with the optimum's exact efficiency, the same condition is a quartic in
sqrt(rho1/lambda), and the horn is on target: s = 1/4, t = 3/8 exactly.
`np.roots` joined the expression whitelist so the spec solves it directly; the
root nearest the point-feed limit is the physical one, checked real.

Checks, none using the spec's algebra: the flares rebuilt as straight lines
from the output dimensions meet the guide together to 1e-9 m; direct 2-D
aperture integration puts the gain on target to 1e-3 dB from 12 to 30 dBi on
WR-90, WR-28 and WR-187; the quartic root matches bisection on the unsquared
condition to 1e-9; a guide shaped like the aperture recovers the old single
apex exactly. The beamwidth constants 54.1 and 78.1, marked indicative, were
integrated too: 53.9-54.4 and 78.0-78.9 from 15 to 25 dBi, drifting 2% by
12 dBi. The feed guide defaults to WR-90's proportions scaled to lambda, and
`axial_length_m` is now the horn's machined length rather than the apex
distance.

## Stevenson's slot conductance was upside down

The handover had carried "Stevenson's g1 should be cross-checked against the
source" since session 5, and `waveguide_longitudinal_slot`'s own note said
its wavelength ratio "is inverted relative to some printings". Reading the
spec against itself settled which way round it had gone: the note quoted
g1 = 2.09 (lambda_g/lambda)(a/b) cos^2(pi lambda/(2 lambda_g)), and the
expression computed 2.09 (lambda/lambda_g)... - a factor (lambda_g/lambda)^2,
1.76 on WR-90 at 10 GHz. `waveguide_slot_array_resonant` used the same
expression and derived its slot offsets from it.

Picking a printing would not have been a check, so the conductance was derived
instead, by quadrature, in `otahub/num/waveguide_slot.py`. The slot's aperture
field is a magnetic current; reciprocity gives the TE10 waves it launches, which
are equal both ways because it couples through H_z - a shunt element. The same
current, doubled by the ground plane, radiates into the half space outside, and
that integral is the slot's one-sided conductance (2 x 73.08/eta^2, to 1e-5).
Power balance for a shunt conductance closes it: g = 2 E_s^2 a b/(Z_TE V^2 G_ext).
Nothing in that uses Stevenson's algebra. It reproduces the lambda_g/lambda form
to +0.03% - the rounding of 2.09 against the derived 2.0893 - at every
frequency and offset tried, 8.2 to 12.4 GHz on WR-90 and 26.5 to 40 GHz on
WR-28, and a real 1.5 mm slot width changes it by 0.35%. The inverted form was
28% low at 12.4 GHz, 43% at 10 and 64% at 8.2, and went to zero toward cutoff
where the right one rises with the guide's wave impedance.

What it did to a design: the twelve-slot WR-90 array at 10 GHz placed its slots
3.06 mm off the centreline. By the first-principles model those slots sum to a
conductance of 1.76, not 1 - VSWR 1.76, 7.5% of the power reflected, at a feed
the spec called matched. The offsets are now 2.28 mm and sum to 1.0000.

Both specs now write g1 as it is derived, 4 eta0/(pi^2 x 73.08) = 2.0893, with
lambda_g/lambda expressed as 1/sqrt(1 - (lambda/2a)^2) so the ratio cannot be
turned over by a typo again. `waveguide_longitudinal_slot` stays low
confidence, and says why: the formula is now right, but Stevenson's model -
a half-wave cosine slot in a thin wall, blind to the guide's evanescent modes -
has still not been checked against a real guide.

## The Taylor taper: a wrong example, and a real problem beside it

The handover's open discrepancies still carried "Taylor taper realises its
design sidelobe level to about 1 dB for small arrays (-28.9 dB measured for a
20-element -30 dB design)", repeated in the function's own docstring. Measured
again - dense evaluation of the array factor, and the package's own
`array_pattern` and `first_sidelobe_db` - that taper gives -30.10 dB there. The
example was wrong. The problem it gestured at was real and worse: the taper was
Taylor's continuous line-source distribution sampled at the elements, and over
357 designs (N = 5..101 odd and even, -20 to -40 dB, every nbar from Taylor's
minimum 2A^2 + 1/2 to 8) it overshot the design sidelobe level by more than
0.5 dB in 156, by up to 2.3 dB, and not only on tiny arrays.

`taylor_nbar` is now Villeneuve's discrete distribution, built by placing the
zeros of the array polynomial: the first nbar - 1 are the N-element
Dolph-Chebyshev zeros stretched so the nbar-th lands on the uniform array's,
and every one after that IS the uniform array's. The weights are the
polynomial's coefficients, from an FFT of its samples taken in log form so a
thousand-element array cannot overflow. Over the same 357 designs it sits
within 0.05 dB above the design level for N >= 10 and at most 0.63 dB below;
one 9-element case reaches 0.36 dB above. Odd N needed care the first draft did
not have: once nbar passes the last zero pair the stretch would have been set
by a zero beyond psi = pi, so there it is now exactly Dolph-Chebyshev.

Checked by routes the construction does not use: realised sidelobes by dense
array-factor evaluation; nulls at the uniform array's positions to 1e-12;
equality with the separate, pattern-sampling Dolph-Chebyshev implementation to
1e-12; and convergence to Taylor's line source re-derived from the textbook
formula, 5e-6 at N = 400. Edge brightening at large nbar and modest sidelobe
levels is Taylor's own - the old taper did it in 78 of 252 cases, this one in
77 - and the docstring now says so. Chebyshev remains the more efficient of the
two for N >= 20.

## Corner reflectors: the directivity was right, the feed impedance was not the feed's

Next on the coverage ranking were the two corner reflectors, 3 of 12 quantities
asserted each. They are the one family on that list the solver can reach
exactly: a dipole in an infinite PEC corner IS its image set, and solving the
set as one method-of-moments problem - every dipole driven with its image sign -
gives the real current rather than an assumed sinusoid. That is now
`otahub/num/corner.py`. Its images conserve power to 1e-5: the image system
radiates exactly as many times the real feed's power as it has wedges.

**The directivity held.** Array factor times dipole directivity times a
resistance ratio, with induced-EMF mutual terms, agrees with the exact solution
to 0.01-0.04 dB from S = 0.1 to 0.75 lambda in the 90-degree corner and 0.35 to
0.8 in the 60-degree one, and the validity note's "nearly flat, 12.6 down to
11.8 dBi" is right with the real current too. Worth recording, because it is
the number most people take from a corner reflector.

**The impedance did not.** The quoted "input resistance" was the induced-EMF
resistance of a sinusoid on an exactly half-wave element: what the directivity
formula needs, and exactly right for that, but not what the feed sees. The
real element presents 6 to 41% more resistance - 151 ohm where the spec said
126, at S = lambda/2 - and a reactance the spec never reported: +49 ohm there,
+107 at S = 0.25 lambda, and at S = 0.1 lambda about 1 ohm of resistance under
40 of reactance. The element that would resonate in the corner is anywhere from
0.42 to 0.52 lambda long depending on spacing and radius.

A closed-form repair was tried first and failed honestly: the isolated MoM
dipole plus induced-EMF image terms, even with a fitted complex scale, is 9-79%
off at its worst point depending on radius and corner - close-range coupling
with the real current does not follow the sinusoidal formulas. So the impedance is solved: 444 image-MoM
solves over S = 0.1..1.0 lambda and element radius 1e-4..5e-3 lambda, fitted
in spacing to degree 12 and log-radius to degree 3. Grid residual under
1.3 ohm; within 2.2 ohm of fresh solves at radii and spacings the fit never saw,
which is inside the solver's own 2-3 ohm mesh uncertainty on the fattest wire.
Degree 10 in spacing had left 5 ohm; degree 14 bought nothing.

Both specs now carry `input_resistance_ohm` and `input_reactance_ohm` from that
fit, with a new `aw_over_lambda` (default 1e-3); the old figure survives as
`induced_emf_resistance_ohm`, named for what it is, and still feeds the
directivity. Known cases come from direct solves, not the fit.

## A dish budget that counted spillover twice

With the coverage ranking run down to archetypes that need a full-wave solver,
the next source of work was the list of metrics the specs themselves label
"indicative", read for any the arbiters built since could now derive. The
prime-focus dish was the most used of them.

Its efficiency budget asserted illumination efficiency 0.82, spillover 0.85,
a beamwidth of 70 lambda/D and a first sidelobe of -24 dB - four numbers that
all follow from one feed pattern and one f/D. And the first two did not fit
together: 0.82 is the textbook optimum of the PRODUCT of taper and spillover
efficiency for a cos^n feed, so multiplying it by a further 0.85 counted
spillover twice. At the default design that cost 0.74 dB of gain.

They are now derived, from the edge taper (a new input, default -11 dB) and
f/D, through Silver's model: a feed with power pattern 2(n+1)cos^n, n set so
the rim sits at the edge taper including spreading loss. Spillover is exactly
1 - cos^(n+1) of the rim angle. The taper efficiency is the feed-angle integral
evaluated in the spec itself, on 4001 points. Beamwidth and first sidelobe are
fitted to the Hankel transform of the aperture field over edge taper -20..-3 dB
and f/D 0.3..1.0. At f/D = 0.4 and -11 dB: taper 0.886, spillover 0.933,
beamwidth 66.5 lambda/D rather than 70, first sidelobe -25.2 dB rather than -24.

The arbiter, `otahub/num/paraboloid.py`, works from the aperture field directly
and shares no algebra with the feed-angle form; the two agree to five figures,
and a cos^2 feed reproduces the classic 0.829 at a 66-degree rim angle and
-10.9 dB. The "-11 dB optimum" holds for f/D 0.4 to 0.6 (best edge taper -10.4
to -10.9 dB) and moves to about -9 dB at f/D = 0.3.

One measurement bug of my own on the way: the first sidelobe fit would not
converge below 1.3 dB, because the search window caught the second sidelobe
once heavy taper had pushed the first below it. Measured properly - the lobe
between the first and second nulls - it fits to 0.06 dB. Past about -17 dB of
edge taper the second sidelobe really is the higher one at larger f/D, and the
spec now says so.

`offset_parabolic` carries the same 0.82 x 0.85 defaults. It needs the offset
geometry's own edge-taper relation rather than a copy, so it is recorded as the
next step rather than patched in passing.

## The offset dish: same double count, and a gain advantage that was a loss

`offset_parabolic` carried the prime-focus defaults, 0.82 x 0.85 with spillover
counted twice, plus a 70 lambda/D beamwidth and a "gain advantage over prime
focus" of +0.09 dB - the value of the 10% blockage it removes, and nothing else.

One geometric fact makes half of it exact. Seen from the focus, the rim of an
offset paraboloid is a circular cone about the bisector of the upper and lower
rim angles. That was checked ray by ray, not assumed: every rim ray lies at
theta* from the axis to 1e-14 rad. So a cos^n feed aimed along it has spillover
exactly 1 - cos^(n+1)(theta*), as the symmetric dish does.

The other half does not carry over. The lower rim is much closer to the feed
than the upper, so spherical spreading tilts the aperture illumination, and
the symmetric formula with theta* overstates the taper efficiency - by 1% at
long focal lengths and up to 20% at parent f/D = 0.3. The taper efficiency and
the beamwidth are therefore fitted to a 2-D integration over the projected
aperture: 864 points over feed taper -20..-3 dB, parent f/D 0.3..1.0 and h0/D
0.55..1.2, held-out errors 0.0013 and 0.03 lambda/D. The two principal planes
differ by under 2%. The offset integrator, set to zero offset, reproduces the
prime-focus arbiter to 1e-9 - two codes written differently agreeing where the
geometries coincide.

What it changed: the input is now the feed's own taper at the rim (default
-10 dB, beside a best of -9.6 dB for the default geometry), since a single
"edge taper including spreading" is ambiguous when the two rims spread
differently. The default design gains 0.48 dB over the old budget. And the gain
comparison, done honestly - against a symmetric dish with the same feed and rim
angle and a 10% blockage - comes out behind wherever the focal length is short
or moderate: up to 0.75 dB at parent f/D = 0.3, 0.09-0.22 dB at the default
0.5, about even at 0.7-0.8, and ahead by at most 0.07 dB at f/D = 1. A first
draft said "negative over the practical range" from three designs; scanning all
864 showed the long-focal-length exception, and the spec says so now. The old
note already said the offset's prize is sidelobes and noise temperature, not
gain; the metric now agrees with its note.

One tolerance lesson on the way: a near-zero dB difference (-0.043 dB) failed a
1% relative check by 0.001 dB. Those expectations now sit in their own case
with an absolute tolerance, as the house rule for zeros already says.

## The parabolic cylinder, one dimension down

`cylindrical_parabolic` asserted taper and spillover of 0.81 and 0.88 and
borrowed the dish's 70 lambda/W for its focusing-plane beamwidth. It is the
dish's analysis one dimension down: a line feed launches a CYLINDRICAL wave,
spreading as 1/sqrt(r) rather than 1/r, so the aperture field across the
focusing plane is sqrt(cos^n(theta)/r). The taper efficiency is then
[int cos^(n/2)/sqrt(1 + cos) dtheta]^2 / (tan(theta0/2) int cos^n dtheta),
which the spec evaluates in its own expression, and the spillover is the
fraction of the feed's transverse power the reflector intercepts. A direct
integration over the aperture coordinate agrees to six figures, and the nearly
uniform limit reproduces the uniform line source's 50.8 lambda/W and -13.3 dB.

Two things came out that the dish did not predict. The old budget understated
the cylinder by a full decibel - the product reaches 0.900 at f/W = 0.4 against
the asserted 0.713. And the dish's -11 dB rule does not carry over: a line
aperture wants a LIGHTER taper, -6.5 to -8.5 dB for f/W 0.3-0.7, so the default
is -8 dB. There the focusing-plane beam is 57.9 lambda/W, and the default fan
beam 3.4:1 rather than 4.1:1.

The sidelobe needed the same care as the dish's, for a different reason. Under
heavy taper the FIRST sidelobe does not merely fall; it collapses, to -51 dB and
-61 dB, as its two nulls merge and the lobe between them vanishes. That is real,
but it is not what a radar designer needs and it will not fit. The spec reports
the PEAK sidelobe instead, which a later lobe takes over smoothly: fitted to
0.3 dB at the corner where the lobes hand over and within 0.08 dB elsewhere on
held-out points.

The dual reflectors are what is left of the efficiency-pair family; they need
a subreflector model before the same derivation applies.

## The dual reflectors, through a traced equivalent paraboloid

The last reflectors with asserted efficiency pairs were the Cassegrain and the
Gregorian: taper 0.85 and spillover 0.90, blockage from the uniform-illumination
(1 - (Ds/D)^2)^2, beamwidth 70 lambda/D. Both specs already described
themselves through an equivalent paraboloid of focal length M*F; nothing had
checked that the equivalence holds. It was checked first, by tracing rays from
the feed off the real hyperboloid and ellipsoid and out through the prime focus:
tan(psi/2)/tan(theta/2) equals M = (e+1)/(e-1) or (1+e)/(1-e) to nine figures,
at M from 1.5 to 9. In geometric optics the dual reflector illuminates its
aperture exactly as the equivalent paraboloid would.

That licensed the prime-focus derivation at the feed's rim angle: n from the
edge taper, spillover past the subreflector exactly 1 - cos^(n+1)(theta_f), taper
efficiency from the feed-angle integral. It also exposed the blockage. The
uniform rule understates it, because the subreflector shadows the brightest
part of a tapered aperture; the geometric-optics blockage is the square of the
unshadowed fraction of the aperture FIELD - 0.965 at Ds/D = 0.1 where the rule
said 0.980, 0.930 at 0.15 against 0.956. Two existing known cases asserted the
uniform 0.9801 and now carry the traced value. The default budget rises 0.21 dB
net: better taper and spillover, worse blockage.

The arbiter works from the traced rays rather than the equivalence: each feed
angle mapped to an aperture radius through the real mirrors, power conserved
along ray tubes. It agrees with the spec's integrals to five figures, the
blockage to 1e-4. Beamwidth and PEAK sidelobe are fitted to the Hankel
transform of the blocked aperture over edge taper, equivalent f/D 0.8-3.0 and
Ds/D 0-0.2 (0.002 lambda/D; the sidelobe to 0.21 dB on held-out points, 0.65 at
the corners where the peak hops between lobes). With no blockage that transform
reproduces the prime-focus arbiter's, and Cassegrain and Gregorian give the same
numbers, as geometric optics says they must.

The scale-invariance audit caught one thing of mine on the way: the
blockage fit variable was centred to be exactly zero at the default Ds/D, and
the audit, rightly, will not accept floating-point noise as a dimensionless
quantity that should not change. It is now Ds/D in tenths, 1 at the default.

That closes the reflector family's asserted efficiencies. The lenses carry
theirs still, and are the next place for the same method.

## Lenses: one flattered by 0.80, one undersold by 0.75

The lenses were the last archetypes with a flat asserted illumination
efficiency: 0.80 for the hyperbolic dielectric lens, 0.75 for the metal-plate
lens, 0.65 for the Luneburg, each silently including spillover, beside a flat
70 lambda/D. A lens decides its own illumination: a ray leaving the feed at
theta refracts at the curved face and crosses the flat one at
rho = |n-1| F sin(theta)/|n cos(theta) - 1|, and power conserved along those ray
tubes turns the feed pattern into the aperture field. The mapping was traced,
not assumed - Snell's law at the actual face, with the normal taken from the
curve: rays leave parallel to 1e-9 and land on the formula to 1e-15.

It cuts both ways. A HYPERBOLA tapers the aperture on top of the feed: -8 dB at
the feed is -22.4 dB at the rim of the default Rexolite lens, and with its rim
angle at 40 of a possible 51 degrees the best taper-spillover product it can
reach is 0.52 - against the asserted 0.80, some 2 dB of gain the spec did not
have. Pulled in to 25 degrees it reaches 0.73. An ELLIPSE (the metal-plate
lens, n < 1) brightens the rim instead: -12 dB at the feed is -8.5 dB at the
rim, and the product reaches 0.89 to 0.95, which the flat 0.75 undersold. So the
two want opposite feeds - broad for the dielectric (-6.7 to -9.3 dB is best),
narrow for the metal plates (-11.8 to -13.6) - and the defaults now say so. The
ellipse also has a limit its spec never stated: the aperture radius peaks at
cos(theta) = n and folds back beyond it, so the rim angle must stay below
acos(n); that is now a reported metric.

Beamwidth is fitted to the Hankel transform of the ray-mapped aperture: 61 to
112 lambda/D for the dielectric lens, where heavy lens taper broadens the beam
a long way past 70, and 54 to 74 for the metal plates. The Luneburg lens keeps
its 0.65, deliberately: in geometric optics a cos(theta) feed gives it a
perfectly uniform aperture, so what keeps real ones near 0.65 is construction
loss, which ray optics cannot see.

## The rectangular DRA was sized for a mode the ground plane forbids

The open item was the rectangular DRA's radiation Q, borrowed from the
hemisphere. Reading the spec to derive it turned up a bigger problem in the
resonance it would have been derived from. The dielectric waveguide model had
been solved with its transcendental along the HEIGHT: magnetic walls on the four
sides, the top face open. That is the even mode of the imaged resonator, 2d
tall - a VERTICAL magnetic dipole, whose tangential E field is largest exactly
where the ground plane sits. A ground plane shorts it out. The physical
broadside mode has its dipole horizontal, and the DWM for it runs the
transcendental along the dipole, with magnetic walls across the width and on
the top face. The two orientations give the same answer only when the imaged
brick is a cube, w = L = 2d - which is the default, and why every known case
passed. A w = 3d, L = 1.5d brick came out 12.6% too big, a cube standing on the
ground 15.9%.

Before rewriting anything I wanted a model that shared nothing with the DWM, and
the toolkit had none for dielectrics: the MoM is thin-wire only. So
`otahub/num/dra.py` is a small 3-D FDTD solver: a quarter of the upper half
space, symmetry planes chosen so that only horizontal-magnetic-dipole modes
exist, CPML outside, a pulse to ring the resonator and a matrix pencil on the
ringdown for the complex frequency. Its check is the exact TE1 pole of a
dielectric sphere, which a staircased hemisphere reproduces to 0.15% in
frequency and 0.9% in Q at eps_r = 10 and 50.

That check found the second problem. The hemispherical spec - "the reference
the others are judged against" - had taken its resonance and Q from the peak
and half-power width of the Mie coefficient |b1| on the real frequency axis.
At Q near 9 the scattering background drags that peak off the pole: at
eps_r = 10 it was 1.2% high in k0*a and 12% LOW in Q (8.04 against 9.15). The
FDTD saw the pole, not the peak; Mongia & Bhartia's published fit sits within
about 1% of the pole too. Both fits are now to the complex pole, 0.23% over eps_r
6-100, in forms that keep the right limits: sqrt(eps_r)*k0*a -> pi, and Q
growing as eps_r^1.5. The old note that its "eps_r^1.32" reproduced the
published eps_r^1.3 law does not survive; the pole is no power law at all.

Then the brick. With the DWM the right way round, 150 ringdowns over eps_r
6-50 and both aspect ratios 1-3 showed it is no better than "a few percent":
10.6% high at eps_r = 6 on long narrow bricks, drifting to 6.2% LOW at
eps_r = 50 on tall ones - where the whole VSWR 2 bandwidth is 1.2%. Its error
does not shrink with permittivity; the magnetic walls it assumes are an
approximation everywhere, not a limit. The DWM's own Q (Mongia & Ittipiboon's
stored energy over magnetic-dipole radiation, which I derived and checked by
direct integration of the fields first) runs from 33% low to 26% high, so it is
not what the spec now uses. The borrowed hemisphere Q had been 6% low at the
default and anywhere from 18% low to 78% high across the design space. Both k0*d
and Q are fitted to the ringdowns instead: cubics in 1/sqrt(eps_r), ln(L/d)
and ln(w/d), within 0.15% and 0.7%, and 0.13% and 0.65% on 33 held-out runs at
other permittivities, other aspect ratios and finer grids. The DWM stays in the
spec as a reported comparison, along with the other broadside mode's frequency
(dipole along w) - the same fit with the aspects exchanged. The default brick
moves 0.3% (3.7441 mm) and its Q from 7.96 to 8.45.

One error of mine on the way, caught because the data would not fit: on long
narrow bars a higher-order mode rang louder at the probe than the dipole mode,
and the first survey took it on seven of 183 runs - resonances 20-30% off with
Q from 14 to 109. The dipole mode is the lowest this symmetry admits, so the
solver now takes the lowest mode with real amplitude, and those seven were
re-run. The exporter's probe note was wrong as well: it still warned of a
"magnetic-wall resonance 10-20% high", two corrections stale. It now says which
mode the probe excites.

The same solver took a first look at `cylindrical_dra`, which rests on
published Mongia & Bhartia fits: five single runs put its Q 12-19% low at
eps_r = 10 and its resonance 6.4% off at a/h = 0.5, while holding at
eps_r = 30. That is the next item - a cylinder is one permittivity map away.

## The cylinder's published fits, and a limit of the ringdown

`cylindrical_dra` was the last DRA resting on someone else's numbers: Mongia &
Bhartia's curve fits for k0*a and Q, "about 2%" over 0.4 < a/h < 6. The FDTD
that checked the brick needed only a cylinder's permittivity map (area
fractions across the axis, the top face grid-aligned); at 16 cells of radius a
finer grid moves k0*a by 0.16% and Q by 0.3% at most.

Sixty-nine ringdowns over eps_r 6-50 and a/h 0.4-6 say the fits hold where they
were presumably made - within about 3% for a/h 0.75-2 at eps_r >= 20 - and
nowhere else. The resonance fit is 6-9% LOW on every puck with a/h <= 0.5, at every
permittivity, so a tall puck built to it resonates that much high. The Q fit
is 6-26% low at eps_r 12 and below, where its eps_r^1.3 scaling has run out,
and up to 16% high on the tallest pucks. At the default, eps_r = 10 and a/h =
1, it read 7.60 against 8.63, overstating the bandwidth by 14%. Both are now
fitted to 48 ringdowns (0.34% and 0.42%; 0.16% and 0.49% on 12 held out), with
the published fits kept as reported comparisons. The three DRA shapes now
agree on how Q scales from eps_r = 10 to 20 - 2.23-2.25 times, exact pole,
brick and puck alike - where the old test allowed 15% around the 2^1.3 the
published fit assumed; it now holds them to 2% of each other.

The survey also found where the ringdown stops. On flat pucks at low
permittivity (a/h = 6, eps_r <= 12) the HE11 mode has Q under 2, and by the
time the window opens - six pulse widths after excitation - it has gone; the
solver reported a higher mode, about 9% off in k0*a with Q of 6-14 where the
published fit says 1.0-2.3. A shorter pulse does not rescue it: it rings every
higher mode and they bury the one it was meant to catch. So the fit stops at
a/h = 4, where the lowest Q in the domain is 1.9 and extraction still holds to
a few percent; the runs beyond are kept in the test data, unfitted, so the
limit is on record. At eps_r >= 15 the flat pucks ring cleanly and the
published fit is within 2.3% in k0*a there, but the domain was kept
rectangular rather than stitched.

The DRA family has no curve-fit physics left: exact pole for the hemisphere,
ringdowns for the brick and the puck. What none of them models is the feed,
a finite ground and dielectric loss.

## The zone plate's gain had no feed in it

`fresnel_zone_plate` computed its gain as the first-order grating efficiency -
1/pi^2, 4/pi^2, 8/pi^2 for opaque, phase-reversing and four-level plates - times
(pi D/lambda)^2: a perfectly lit aperture with no feed in sight. The lens and
dish rounds had each found that the feed's taper and spillover are where most
of an aperture budget goes, so this one was worth doing the same way.

The arbiter, `otahub/num/zone_plate.py`, is scalar Kirchhoff diffraction done
bluntly: a 2(n+1)cos^n feed at the focus, spherical spreading to the plate, each
point's transmission from its own path length, obliquity (1 + cos)/2, and the
radiation integral summed over a Cartesian grid on the plate, gain referred to
the feed's total power. The spec evaluates the same physics as a feed-angle
integral on 8001 points; the two agree to 0.015%. Under near-uniform light on
many zones the integral falls back to 1/pi^2, 4/pi^2 and 8/pi^2, so the old
numbers were the right limit of the right model - just not the case a feed
produces.

For a phase plate the budget factorises exactly, to better than 1%: grating
efficiency times the feed's taper and spillover. That factor is 0.72-0.81 even
at the best feed taper, so the old gain was 0.9-1.4 dB optimistic; at the
default design, phase reversal, 29.13 dBi rather than 30.11. The opaque plate
does NOT factorise, and that was the surprise. Its unfocused zero order passes
the open zones as well, and under a tapered feed the bright central zone, which
is open, carries it into the beam: with few zones and a heavy taper the opaque
plate beats 1/pi^2 (0.115 at best on the default four zones), with many zones
it falls short, so its old gain was anywhere from 0.9 dB pessimistic to 1.3 dB
optimistic. Its best feed taper is heavier too, -6 to -18 dB against -6 to -11
for a phase plate.

The zero order also shapes the opaque plate's beam, and makes it depend on the
parity of M: an odd count ends on an open zone, an even one on a blocked one,
and their beams differ by up to 14%. Beamwidth is therefore fitted per plate
type and, for the opaque plate, per parity - quadratics in feed taper, 1/M and
1 - cos(rim angle) over 1260 patterns, 0.79% worst and 0.18% on held-out 2-D
patterns. A phase plate's beam is the tapered aperture's, 61-92 lambda/D; the
opaque plate's reaches 116. The old flat 70 lambda/D suited a phase plate under
a -10 dB feed; the default opaque plate's beam is 16% wider.

The "INDICATIVE 1/M" bandwidth held up better than anything else in the spec.
Swept in frequency with the plate fixed, M times the 1 dB gain bandwidth runs
0.89-1.37, so 1/M was the right size; it is now computed the same way (0.70%
fit, 0.30% held out), and the 3 dB band is 1.69-1.81 times wider. One detail on
the way: a four-level plate whose steps each compensate at the lower edge of
their band peaks above f0, 6% at four zones.

## A solver for solid cones, and a discone whose feed was part of the answer

The discone carried three conventions: a quarter-wave slant at f_low, a decade
of VSWR under 2 ("indicative"), and a flat 2.2 directivity. The wire-cage model
the cones round built could not settle any of them - a cage's impedance does not
converge in the wire count - so this round built `otahub/num/bor.py`: the method
of moments for a body of revolution driven symmetrically. The current is a ring
current along the generating curve, rooftop basis, Galerkin testing, the wire
kernel replaced by its exact ring averages; for each azimuth a straight
generating segment is a straight line in space, so the static part integrates in
closed form and phi = pi u^2 removes the log singularity left behind. The feed is
a finite gap - a uniform field over a stretch of the curve - so a fat body gets a
physical feed rather than a delta gap whose capacitance diverges.

It was checked four ways before it touched the discone: a thin tube against the
wire code's exact kernel (0.3%), input power against radiated power through a
separate far-field integral (1e-5), a 47-degree bicone against Schelkunoff's
99.9 ohm (within 1% in resistance), and - for the disc's purely radial current -
a cone over discs of 1.5, 3 and 5 wavelengths converging on image theory's half
bicone (7.7%, 3.7%, 2.1%).

Then the discone, and a lesson in how much the feed decides. The first survey
said the decade did not exist - about 3:1 at 30 degrees - and the second, with
the feed made to match the spec's gap rule, said the cutoff moved up by a fifth
or more.
Both were artefacts of the feed. My first geometry put the ideal apex a gap below
the disc and then truncated the cone, which quietly made the real spacing nearly
four times the stated one; the second took the spec's own rule at its word, and that
rule implies a cone top a third of the base RADIUS - a capacitor plate facing
the disc. With a coax-sized top (2% of the base diameter, Nail's spacing of 0.3
times it), converged in the feed region, the picture reversed:

- The decade is real for cones of 30 degrees or more at every disc ratio
  surveyed (VSWR under 2 right through 10 f_low, worst 1.48-2.00 between 2 and 10
  f_low) and not for narrower ones: 3.1:1 to 5.1:1 at 20-25 degrees, the one
  exception being 25 degrees with a 0.6 disc. The
  indicative figure was right for the default and has limits now.
- The quarter wave is short everywhere: the VSWR-2 cutoff sits at 0.257-0.395
  wavelengths of slant, 12% longer than a quarter wave at the default and 58% on
  a 50-degree cone. Wide cones match 50 ohm better once they work, and start
  working later.
- Directivity at f_low is 1.24-1.54, 1.48 (1.71 dBi) at the default, where the
  spec said 2.2; it climbs to 2.35 at 3 f_low as the beam tilts 50 degrees below
  the horizon.

The spec now tabulates all four over half angle 20-50 degrees and disc ratio
0.6-0.9 (the low cutoff extrapolated from two meshes, 0.9% on six held-out
designs), states the feed as geometry - `cone_top_diameter` and the spacing -
and the exporter truncates the cone there. A test rebuilds the antenna from the
spec's own dimensions and solves it live: VSWR 2 at f_low, as the table says.

The survey also cost a restart: the machine hibernated mid-run and the
28-design job had to be redone. It is recorded in tests/data so nothing needs
re-solving to check the tables.

## The bicone's decade holds, the monopole's only sometimes

`biconical` and `conical_monopole` sized their cones with a quarter-wave slant
and asserted an "indicative" decade, with directivities fitted from a wire cage.
`otahub/num/bor.py`, built for the discone, solves both: a conical monopole on an
infinite plane is exactly half a bicone by image theory - half the impedance,
twice the directivity - so one survey of solid bicones, 5 to 65 degrees, swept to
ten times the cutoff and extrapolated from two meshes, serves both. Both carry
the discone's coax-sized feed (cone top 2% of the base diameter), stated as
geometry, and the exporters now truncate the cones there.

What decides the answer is the reference. Matched to its own characteristic
impedance - the best any line can do - a bicone holds VSWR under 2 for more than
a decade at every angle from 10 to 65 degrees, from a slant of 0.22 to 0.29
wavelengths: the quarter wave is right to 15% there, and the indicative decade
was right. A 5-degree cone breaks off at 1.26:1 before settling again, and a
65-degree one needs 0.33 wavelengths.

A conical monopole is fed from 50 ohm, and there the cone angle is everything.
From 30 to 55 degrees it holds a decade (from 0.22-0.23 wavelengths at 35-47
degrees - a quarter wave is 8% long at the default 47). Below that the cone's
impedance climbs away from 50 ohm, and the cutoff with it, steeply: 0.22 at 34
degrees, 0.30 at 33, 0.35 at 32, 0.40 at 30. The early near-match at a fifth of
a wavelength stalls just above VSWR 2 - 2.04 at 32 degrees, 2.08 at 30 - so the
50 ohm band only begins once the cone is long. A table at 5-degree steps would
have interpolated straight across that, so it carries extra nodes at 32, 32.5
and 33. At 15-20 degrees the band is 1.2-1.5:1, below 15 there is none, and at 65
degrees (27 ohm) it is 1.6:1.

Directivity is now quoted at each antenna's own cutoff. The old cage fits were
evaluated at a quarter wave and ran 0.7-2.8% above the solid cone even there -
a cage is not a surface. Because the two specs now define their cutoffs against
different references, "monopole directivity = twice the bicone's" no longer holds
at every angle; it holds exactly where the cutoffs coincide, at 47 degrees, whose
monopole Zc is 50 ohm, and the cross-check now tests that instead of an identity
the definitions no longer imply.

One more survivor of the old angle bug: the biconical exporter still halved
theta_h, drawing cones half as wide as the design. The spec was fixed several
rounds ago; the builder was not, and nothing tested the angle it drew. It is
fixed, with a test on the flare.

## Fat dipoles: where the feed takes over

HANDOVER item 9 asked for a finite-gap feed so the slot family could be checked
past w/L = 0.05, where the complementary wire (a = w/4) is too fat for a delta
gap. `otahub/num/bor.py` already has one: a straight conductor is a body of
revolution, a solid tube with a gap of physical length. So the check was run
there - tube dipoles of radius 0.002 to 0.025 wavelengths, each with gaps of
a/2, a, 2a and 4a - and the useful result is the spread, not a number.

Past a = 0.005 wavelengths the resonance belongs to the feed as much as to the
wire. The resonant length moves +-0.3% across those gaps at a = 0.005, +-0.9% at
0.01 and +-2-4% at 0.015-0.025; the resistance at resonance spans 72-74, 73-78
and 74-96 ohm. At a = 0.025 with a gap of a/2 the reactance tops out at -6 ohm:
the gap capacitance wins and there is no resonance at any length. A delta gap is
the limit of shrinking that gap, which is why it never converged.

`resonant_dipole`'s length law lies inside the spread at every radius tested, out
to 0.025 wavelengths - better than its own validity claimed. Its resistance law,
fitted to delta-gap solutions, sits at the small-gap end of the spread, 2-4%
above the finite gaps at 0.005. The slot family now says what that means for it:
checked to within the feed's +-0.9% up to w/L of about 0.09, and past that its
resonant length and resistance mean something only with the feed modelled - at
w/L = 0.2 the complementary dipole spans +-2-4% and 74-96 ohm, a 370-480 ohm slot.
What stays unchecked is the strip-to-tube equivalence a = w/4 itself for wide
strips; that needs a strip solver the toolkit does not have.

## Patch bandwidths from a solved Q, and an air patch that had none

Three patch shapes carried the same bandwidth: 3.771 (eps_r - 1)/eps_r^2 x
(h/lambda0)(W/L). It is Jackson and Alexopoulos's, simplified twice - their
substrate factor c1 = 1 - 1/eps_r + 2/(5 eps_r^2) cut to 1 - 1/eps_r, their
pattern factor p set to 1 - and the circular and triangular patches borrowed it
with W/L = 1, marked "indicative for this shape". The corner-truncated CP patch
took its Q0 from it, and Q0 sets how big the corners are cut.

`otahub/num/patch_q.py` computes the radiation Q directly: the stored energy of
each shape's own cavity mode against the space wave of the patch's surface
current lying on its grounded substrate, through the transverse-equivalent-
network slab factors - the physics Jackson's formula was derived from, before
its thin-substrate expansion. Two checks that share nothing with it: a
vanishing patch on a thin substrate reproduces Jackson's c1 to 0.02%, and on
air its broadside directivity reproduces the cavity model's edge-current
integrals - an electric surface current against magnetic edge currents - to
0.06% for the disc and 0.003% for the triangle. Against Jackson's full formula
it runs 3-10% high, most of it in p, whose Taylor series overstates the pattern
factor of an electrically large (low-permittivity) patch; the rest is substrate
thickness, which the thin-substrate c1 leaves out.

What the solved Q says about the formulas:
- The simplified formula gives an AIR patch zero bandwidth, through its
  (eps_r - 1) factor - the wideband case. It reads 5-45% low at eps_r 1.25-1.5
  and 0-33% high from 2 upward on the rectangle.
- The disc was lucky: within 13% from eps_r 2 to 12, but 27-51% low at 1.25-1.5
  and zero on air.
- The triangle was 40-64% too broad from eps_r 2 up. It stores more energy for
  each watt it radiates than the square it was borrowed from.
- The CP patch's Q0 was 5-11% low (62.3 against 59.0 on the 2.4 GHz design), so
  every corner cut was 2.5-5.5% too large. The cut follows the solved Q now;
  its cavity-model test needed its meshes re-matched to the smaller cut, and the
  design is still circular at f0 on them.

All four are fitted as ln(Q (h/lambda0) c1/eps_r) in 1/eps_r and
h sqrt(eps_r)/lambda0 - Jackson's thin-substrate law carries the bulk, so twelve
terms reach 0.6% (0.25% held out) - over eps_r 1-12, h/lambda0 0.003-0.05 and
h sqrt(eps_r)/lambda0 up to 0.1, where the cavity model still holds. Radiation
only: dielectric and conductor loss would widen every band, and surface waves
are not counted. The quarter-wave shorted patch keeps its indicative halving:
its shorting wall carries a vertical current this model does not have.

## Horn beamwidths that let the flare vary, and beams that break up

Every horn in the family quoted its beamwidth as K lambda/D with K fixed - 54 and
78 for the E- and H-plane sectoral horns, 60 and 70 for the smooth conical horn,
66 for the corrugated one - "indicative, near the optimum flare only", while the
same specs take the flare as an input and compute the efficiency exactly at any
value of it. `otahub/num/horn_pattern.py` does the beam the same way: the
aperture field (uniform, cosine, the circular guide's TE11, the corrugated
horn's J0 HE11) with its quadratic phase error, Huygens obliquity, a Fourier
integral by quadrature and a root-find for the half-power point. Its E-plane
pattern matches the closed form in Fresnel integrals, which it does not use, to
1e-10.

What it found:
- The sectoral coefficients are right exactly where they were derived: 53.8-
  53.9 and 77.6-77.9 at the optimum flare. Away from it they are not - 50.7 and
  68 with no phase error, 63.7 at s = 0.4, 107 at t = 0.6 - and the H-plane horn
  at the flare of 1.4 one of its own known cases uses has a beam nearly twice 78.
- Over-flared E-plane horns break up: at an edge phase error of 0.47 the
  half-power point jumps from 87 to 126 lambda/b1 as the main lobe splits, for
  every aperture size tried. The smooth conical horn's E plane collapses the
  same way past s = 0.55 - at the 1.3 flare one of its known cases uses it is
  156 lambda/dm against the 60 quoted. Past those points the spec now returns
  NaN rather than a number.
- The conical coefficients were 7-8% narrow at the default flare (64.3 and
  75.9, E narrower than H as the note said), and the corrugated horn's 66 was
  wrong at any flare: a J0 field tapering to zero at the rim gives 75.7 with no
  phase error at all and 80.2 at the default. The one figure for both planes
  survives; its value did not.

Each plane is fitted as the half-power point U = D sin(theta_h) in phase error
and 1/D^2, so HPBW = 2 asin(U/D) keeps the sine mapping exact for small
apertures - 21 terms, 0.5% over 1.5 to 40 wavelengths. The pyramidal horn,
designed at the optimum by construction, keeps its 54.1 and 78.1: aperture
integration with obliquity agrees to 0.5%.


## Two-mode horns: a dual-mode "advantage" that is a loss, and cross-polar lobes 3.5 dB worse

The diagonal horn and Potter's dual-mode conical horn were the last horns with
fixed figures. The dual-mode horn carried an INDICATIVE 0.62 efficiency (so +0.85
dB over a smooth horn), 68 lambda/D "equal in E and H" and -30 dB cross-polar;
the diagonal horn 58 lambda/a and cross-polar lobes "typically around -19 dB, not
modelled". `otahub/num/horn_pattern.py` now transforms BOTH components of a
two-mode aperture, which gives Ludwig-3 co- and cross-polar patterns from one
integral. Its TE11 limit reproduces the scalar TE11 routine to 1e-9 and the
textbook 0.837; the diagonal horn's efficiency matches a closed form in Fresnel
integrals to 1e-6.

What it found:
- Dual-mode efficiency is 0.506 at the default TM11 fraction of 0.15, not 0.62.
  The 0.62 sat between the smooth horn at this phase error (0.51) and the
  corrugated horn at none (0.69) - two numbers that were not comparable. Against
  the smooth horn at the same 3/8 phase error it is a 0.27 dB LOSS: the TM11
  tapers the E-plane rim, which equalises the beams and cancels the cross-polar
  field, and a tapered aperture has less gain.
- Beams are equal at a TM11 power fraction of 0.097, at 78 lambda/dm, not 68; at
  the default 0.15 the E plane is already 4% wider (82.0 against 78.6).
- Cross-polar nulls at a fraction near 0.13 (-34.1 to -34.5 dB depending on
  size); -32.2 dB at the default, so the -30 was about right. The two design
  goals pull a third apart.
- The diagonal horn's cross-polar lobes, in the planes 45 degrees from the
  polarisation, are -15.5 dB for a long horn, -14.7 at an edge phase error of
  1/4 and -12.4 at 1/2 - not -19. Small horns are a little cleaner (-16.1 dB at
  2 wavelengths), so the fit carries 1/D^2.
- Its 58 lambda/a beam was right for a long horn (58.3), and E and H are equal
  by symmetry, exactly; the beam widens to 61.2 at an edge error of 1/4, 81.9 at
  1/2. Its efficiency now follows the flare, exactly - the E-plane sectoral phase
  factor times the H-plane sectoral efficiency - and so does its gain.
- Its notes put the phase error (a/2)^2/(2R) "at the aperture corner"; that is
  the edge centre, and the corner has twice it, so R >= a^2/lambda holds the
  corners to a QUARTER wave, not an eighth.

Dual-mode efficiency is exact as a ratio of four mode integrals; beams are fitted
in sqrt(fraction) and 1/D^2 to 0.002%; cross-polar is tabulated at 3.2, 5.5 and
12 wavelengths and interpolated in 1/D^2 (0.03 dB on eight held-out solutions);
the diagonal horn's beam fit holds 0.07% and its cross-polar fit 0.003 dB. NaN
outside the surveyed range. The Potter horn stays low confidence: the step that
sets the TM11 fraction needs mode matching, so the fraction is an input, and the
bandwidth stays an estimate.


## The shorted patch is not half as wide: its side walls radiate

`quarter_wave_shorted_patch` carried half the rectangular patch's simplified
bandwidth formula, "INDICATIVE - this halving has not been computed". The
reasoning was that Q doubles when one of two slots is removed. It does not: the
stored energy halves along with the patch, so even in the slot picture Q only
rises by 1 + G12/G1. And the open edge is not the only radiator - the side
edges carry sin(pi x/2L) fields that add in phase, and the shorting wall a
vertical current down to the ground.

`otahub/num/patch_q.py` now carries that wall current: K_z equal to the patch
current arriving at the wall, uniform over the height, radiating TM through the
slab with a factor G j sin(theta)/(k0 N^2). The factor was derived from the
fields of a TM plane wave in the slab and is checked against a direct
boundary-value solve to 1e-8; on air it is exactly a vertical current and its
image. The whole shorted patch on air reproduces the cavity model's magnetic
currents on its three open walls, written independently, to 0.004%. Without
the wall current the Q is 3.4 times too high.

What it found, as the ratio of the shorted patch's bandwidth to the full patch's
on the same board:
- 1.83-1.85 on AIR: the shorted patch is nearly twice as broad, not half. The
  old formula gave it zero.
- About 1 near eps_r 2; 0.66-0.72 on FR-4; 0.56-0.57 at eps_r 12. It approaches
  a half as the patch shrinks and the slot picture takes over, and never gets
  there inside the domain.
- The old formula was 2.4-5.5 times too narrow at eps_r 1.25-1.5 and 3-53%
  narrow at eps_r 2-3.5. From eps_r 6 up its two errors - the halving and the
  simplified Jackson formula - roughly cancel, to between 16% narrow and 5% wide.

The Q is fitted in the same form as the other patches (12 terms, 0.39%, 0.23%
held out) over eps_r 1-12, h/lambda0 0.003-0.05, h sqrt(eps_r)/lambda0 <= 0.1.
The spec now also reports the full patch's Q on the same board and the
bandwidth relative to it.

Found on the way, and recorded rather than fixed: patch DIRECTIVITY. The
shorted patch's single-slot figure ignores the side walls that dominate its
bandwidth - the complete free-space cavity model gives 5-26% less - and the
rectangular patch's two-slot figure ignores its side walls too, which raise it
2-10%. The substrate (Jackson's space wave, `patch_q`'s slab model) pushes
both up by as much as 13-15% at eps_r 10.2. Both specs' directivity notes now
say so; HANDOVER section 6 item 12 describes the family-wide revision.


## Patch directivity, settled full-wave: every rectangular figure was 7.5-13% low

The shorted-patch round left patch directivity in question. The rectangular
family carried the two-slot formula; circular, triangular and annular carried
free-space edge-current models; `patch_q`'s cavity current radiated through the
slab disagreed with all of them. Deciding which is right needed a solver that
does not assume the current, so this round built one.

`otahub/num/patch_sdm.py` is a spectral-domain method of moments for a
rectangular patch on a grounded slab: the exact spectral Green's function,
entire-domain edge-conditioned currents in both directions, Galerkin testing, a
contour over the surface-wave poles, reactance extrapolated from two
truncations. It was checked before being believed: a vanishing dipole's space
and surface waves against Jackson and Alexopoulos's thin-slab closed forms
(0.4-2.4% and 0.002-0.9%, residues growing with thickness as they should);
microstrip effective permittivity against Kirschning and Jansen (0.1-0.9%); the
radiated power from the spectral reaction against a far-field integral (1 part
in 10^4). Two things were caught on the way: an azimuthal quadrature with half
the samples a product of two spectra needs (the reactance then wandered with
truncation), and a Maxwell-only transverse profile that put a very wide patch's
resonance 3.5% off - the current across a wide patch is flatter than the edge
factor alone.

What it found, directivity of the spec's own rectangular designs (largest basis):
- air 9.53, eps_r 2.2 5.65, FR-4 4.39, eps_r 10.2 3.99; `patch_q` within 0.3-0.9%
  of every one.
- The two-slot formula: 8.58, 5.22, 4.06, 3.47 - 10%, 7.5%, 7.5% and 13% low. It
  leaves out the non-radiating edges (2-10% of the answer, most on low eps_r)
  and the substrate (up to 22% on electrically thick board).
- The circular and triangular free-space models are exact on thin board but up
  to 18% low on thick; the CP square's two-slot figure 6-10% low; the shorted
  patch's single slot 5-35% high in the thin limit, and low on thick board.
- The triangle is 0.07-0.33 dB below the rectangle, not "within 0.1 dB" - that
  agreement was the rectangle's missing side walls.

Six specs (rectangular, inset, CP, circular, triangular, shorted) now carry
directivity as a thin-substrate value times a substrate factor. The thin value
is exact as the board vanishes and goes to 3 as the patch shrinks: the circular
and triangular keep their edge-current models for it, the rest fit patch_q as 3
plus a series in 1/sqrt(eps_r) (0.074%). The factor is patch_q over thin,
fitted with no constant term so it is exactly 1 on thin board (0.23%, 0.16% held
out), and the directivity is NaN past h sqrt(eps_r)/lambda0 = 0.1. The circular
patch's edge resistance stays on the free-space model and now says so.

Not settled, and recorded as HANDOVER section 6 item 13: the full-wave solver
puts the Hammerstad-length patch 1-2% low in frequency on thin board and about
4% low at eps_r 10.2, and its radiation Q 7-10% below patch_q's on thick FR-4
and 10.2 - but at eps_r 10.2 those move 1.7% and 14% across bases, too much to
redesign on.


## The textbook patch resonates low: 0.5-7%, the 2.4 GHz FR-4 patch at 2.33 GHz

Last round's spectral MoM settled directivity but left resonance and Q open:
they moved 1.7% and 14% with the basis. The cause was the basis. Sines along the
current vanish linearly at the open ends; the real current vanishes as the
square root of the distance. With U_n sqrt(1-t^2) along the current (closed-form
spectra in Bessel functions) the resonance moves 1e-4 from the default basis to
the next. The solver was also made about thirteen times faster: one quadrant by
symmetry, the two truncations sharing their common range, a secant search in
place of a frequency scan.

Then it was checked against something that shares nothing with it but the
geometry. `otahub/num/patch_fdtd.py` rings the patch down in `dra`'s Yee grid - a
PEC sheet on a slab running into the CPML - at three cell sizes of an exactly
representable geometry, extrapolated at the first order the sheet edge gives
(successive differences in ratio 2.03 and 2.04). On a thick eps_r 10.2 patch:
FDTD 0.95854 f0 and total Q 46.90, the MoM 0.95866 and 46.80. On FR-4: 0.99121
and 60.77 against 0.99141 and 60.86. The sine basis had put the first at
0.951-0.986.

What the full-wave solver found on the rectangular spec's own designs (68, and
8 held out, over eps_r 1-12 and h/lambda0 0.003-0.05):
- The transmission-line design with Hammerstad's fringing extension resonates
  LOW: 0.54-0.96% on the thinnest boards, 2.2-2.8% at h = 0.01 lambda0, up to 7%
  on thick high-permittivity board. The classic 2.4 GHz FR-4 patch (L = 29.42 mm)
  resonates at 2.332 GHz; 28.52 mm resonates at 2.4. The spec's validity allowed
  "1-3%".
- The cavity-current radiation Q (`patch_q`), which the bandwidths came from,
  runs 1.6-28% HIGH - least on thin board, most on thick - so every rectangular
  bandwidth was that much too narrow.
- The surface-wave efficiency matches Jackson and Alexopoulos's thin-slab closed
  form to 0.0-0.5% on thin board and falls 7.5% below it at eps_r 10.2,
  h = 0.03 lambda0; 0.66 at eps_r 6, h = 0.04 lambda0.
- Directivity at the full-wave-resonant design is within 0.02-2.5% of the
  cavity-current value the specs carry.

`rectangular_patch` keeps L as the textbook design (Balanis's worked example
still reproduces) and now reports full_wave_resonance_hz, L_full_wave_m (the one
to build), a full-wave radiation_q, surface_wave_efficiency and a bandwidth that
counts the surface wave; the inset patch carries the two geometric ones. Fits:
0.20%, 0.21%, 0.9% and 0.11% (0.13%, 0.12%, 0.37%, 0.13% held out). The circular,
triangular and shorted patches' cavity-current Q, and the CP patch's Q0 and so
its corner cut, are now labelled as not checked full-wave. HANDOVER section 6
item 13 lists what building L_full_wave_m would take.


## The rectangular patch now builds the length that resonates

Last round showed the textbook length resonates 0.5-7% low and left the spec
reporting it rather than building it. `rectangular_patch` and
`rectangular_patch_inset` now synthesise L as the full-wave resonant length (the
fitted ratio to the textbook length, 0.21%, from the survey the FDTD checked),
and keep the transmission-line design as `L_textbook`, so every worked example
still reproduces - under that name. Balanis's Example 14.1 patch, solved directly:
the textbook 0.9053 cm resonates at 9.534 GHz, 4.7% low; 0.8556 cm resonates at
10 GHz. The inset patch keeps its transmission-line edge resistance and feed
FRACTION on the textbook length they were derived for, and cuts the notch that
fraction of the full-wave L: Example 14.2's inset is 0.2952 cm deep on the patch
that resonates, not 0.3124 cm on one that resonates at 9.53 GHz. The exporter
builds from L, so it now builds the patch that resonates.

Found on the way, an engine bug: an expression guarded with np.where came back
as a 0-d array rather than a number. The exporter refused it outright, and the
plausibility test skipped it - so metrics that go NaN outside their fitted
domains (horn beamwidths past break-up, patch figures on too thick a board) had
never been seen by it. The evaluator now unwraps 0-d results; the plausibility
test lets deliberate NaN through, because an expected value that comes out NaN
is caught by the reference match. The rectangular patch's surface-wave
efficiency is capped at 1, which its fit overshot on air by parts in 10^4.

Still on textbook lengths, and said so in HANDOVER section 6 item 13: the CP
square (the MoM handles a square, and its Q0-driven corner cut is also off), the
shorted patch (needs a shorting wall the MoM lacks) and the stacked patch.


## The CP square, full-wave: a smaller square with a bigger cut

`truncated_corner_cp_patch` took its square side from the same textbook chain
that resonates 0.5-7% low, and its corner cut from Q0 - the cavity-current
radiation Q, which runs high and leaves the surface wave out. A single-feed CP
patch lives on the ratio of its mode split to its modes' damping, and the
surface wave damps each mode as surely as radiation, so Q0 is now the full-wave
radiation Q times the surface-wave efficiency; the square side is the one that
resonates at f_sq full-wave; the cut-to-split relation stays the cavity model's,
which `patch_cavity` checked on the truncated outline and the MoM cannot model.

The MoM was checked on a square first: an FDTD ringdown of a thick eps_r 10.2
square, extrapolated in cell size, gives 0.98403 f0 and Q 57.56 against its
0.98396 and 57.46.

What changed, on the spec's own designs (68 squares and 8 held out):
- The old Q0 was 1-99% above the full-wave total - up to 28% from the radiation
  Q itself, the rest from the surface wave it left out, which nearly halves the
  Q of a square on thick high-permittivity board. Every cut was too small and
  every axial-ratio band too narrow.
- The textbook square is 0.54-7.1% too large.
- The classic 2.4 GHz, eps_r 2.2, 1.6 mm design: Q0 53.9 (was 62.2), square
  40.70 mm (was 41.77), cut 3.975 mm (was 3.793). A 1.575 GHz FR-4 GPS patch:
  Q0 118.2 (was 141.9), square 44.66 mm (was 45.54), cut 2.931 mm (was 2.725).

Fits: side ratio 0.21% (0.12% held out), radiation Q 0.87% (0.67%), efficiency
0.11% (0.13%). The cavity-model design test re-matched its meshes to the new cut
(246 cells, a 24-cell cut, c/L to 6e-5) and still finds under 0.1 dB of axial
ratio at f0, and over 7 dB for the classic cut; the first-order offset test now
allows the cut's second-order term, 5% at this Q0. A supplied Q0 still
overrides, for loss the full-wave Q does not include.


## The Potter horn's step, solved - and a flare that moves the phase 20 degrees

`conical_horn_dual_mode` had its aperture derived for a given TM11 share but took
that share as an assumption: the step that launches it was "indicative, 1.3
wavelengths". `otahub/num/waveguide_step.py` solves the step by mode matching
over the m = 1 modes (TE1n and TM1n, counts in proportion to the radii), and
`otahub/num/bor_fdtd.py` checks it: an FDTD at azimuthal order 1 on a (rho, z)
grid, which shares nothing with it. Taken to zero cell size (order 1.4, set by
the step's corner) the FDTD gives a TM11 share of 0.17682 at 41.11 degrees where
mode matching gives 0.17671 at 41.15.

What the step does:
- It launches TM11 20-60 degrees out of phase with TE11, not in phase; the
  phasing section has to make that up.
- The input guide matters. From a 1.0-wavelength guide a share of 0.15 needs a
  step only 2.7% above TM11 cutoff, and at 0.97 f0 TM11 no longer propagates;
  from 1.1 the step is 11% clear, the share holds to 0.002 over +-2%, and 0.1% of
  the power is reflected. The spec now takes the input guide as an input (1.1
  typical) and synthesises the step for the share asked (tabulated, within 0.002
  of it; NaN where the guide cannot launch it).

What the flare does, found by cascading the mode matching over a stepped cone
(generalised scattering matrices) and checked by an FDTD of the whole
10-wavelength horn in its own staircase: the cascade's aperture projections sit
where the FDTD converges (0.367 at 37.4 degrees; FDTD 0.355 and 0.363 at 60 and
120 cells per wavelength). A plain cone barely converts (0.3% of the power), but
coherently: the Potter horn's aperture phase lands at 15 degrees where the
first-order chain - local propagation constants, no conversion - puts it at -7,
and its frequency slope is 11.4 degrees per percent, not 7.7. So the phasing
length the spec now gives is first-order and labelled INDICATIVE, with that
error stated.

Cross-polar bandwidth, derived rather than asserted: the aperture model's
tolerance to phase error is lopsided - TM11 leading spoils the null fast (+5.5
degrees at a share of 0.13), lagging slowly (-19) - so a design can centre its
band by phasing about 6 degrees long. With the first-order slope the -30 dB band
of a 10-wavelength horn at 0.15 is about 2.5%, 1.8% with the cascade's; the -25
dB band about 6.4% and 4.5%. Longer horns are narrower (-17 degrees per percent
at 33 wavelengths). The 6% the spec carried is roughly the -25 dB figure.

The horn stays low confidence, for a narrower reason than before: the step is
solved, the phasing is not. HANDOVER section 6 item 11 describes the cascade-based
design that would finish it.

## Finish line

Set by the user on 2026-09-26, through the VS Code session. It governs every
round from here on; the house rules are unchanged.

1. **Scope freeze.** Pick only from items already open in HANDOVER section 6 as
   it stood at commit 53df1ab. The order after the Potter horn round:
   - item 13: the shorted patch, and the stacked patch's length, full-wave;
   - item 12: the annular ring mode in patch_q, and PIFA;
   - item 9: the wide-strip equivalence;
   - item 5: array element-pattern embedding;
   - item 4: taking one low-confidence archetype end to end.
   Anything newly discovered goes into section 6 as future work, not a new round.
2. **Hard stop at 2026-09-27 20:00 local time.** No new round starts after it. A
   round in progress at the deadline is finished only if it can be green (suite,
   doctor, check) and committed within about an hour; otherwise its work goes on a
   local branch `wip/<topic>` with a commit message saying what is unfinished, and
   master returns to the last green commit.
3. **Stop early** if section 6 runs out first: every remaining item needs
   CST/HFSS, measurement data, or a new solver judged out of scope.
4. **Wrap-up round, the last one:** rewrite HANDOVER section 6 as the final
   future-work list, grouped by what each item needs (CST/HFSS installation,
   measurement data, or a new solver); update the README and HANDOVER counts; run
   the full suite, doctor and check; add a final "Development closed" section to
   this file, last as always; commit and tag that commit v1.0 (local tag, no
   remote, do not push); stop the loop with ScheduleWakeup stop: true and do not
   re-arm it.


## The shorted and stacked patches: a wall that pulls lower, a stack sized backwards

HANDOVER section 6 item 13, the first round under the finish line. The spectral
MoM has neither a shorting wall nor a second metal layer, so both are done with
the FDTD. `otahub/num/patch_fdtd.shorted_ringdown` runs a quarter-wave patch
with a full-width shorting wall in a half-space grid: CPML on both ends of the
length, since a shorted patch has no mirror plane along it. On a symmetric
(unshorted) patch the half-space grid gives the quarter-space answer exactly.

The shorted patch, textbook length, at FDTD cell sizes taken to zero:
- eps_r 10.2, h = 0.0245 lambda0: 8.4% low (the full patch on the same board
  5.9%); three grids, first-order convergence (difference ratio 2.03).
- eps_r 2.2, h = 0.02 lambda0: 5.9% low (full patch 3.6%); two grids.
- FR-4, h = 0.0128 lambda0: about 4% low (full patch 2.8%).
The shorting wall adds its own lowering - the inductance of the current it
carries to ground - which the transmission-line model has no term for, and it
does not scale simply with h across these three boards. The spec now builds the
textbook length times the full patch's full-wave ratio (keeping L_textbook),
which closes most of the gap: the FDTD puts the corrected patch 1-3% low, and
the spec says so and labels L an estimate. Its cavity-current radiation Q,
times the surface-wave efficiency, runs about 10% above the FDTD's total Q on
two boards - the rectangular patch's direction too.

The stacked patch: the spec sized the parasitic at 0.95 of the driven patch,
because "sitting in a lower-permittivity medium, it is usually made slightly
SMALLER for the same resonance". That is backwards: in air a patch needs to be
longer. An FDTD of the default stack (1.6 mm eps_r 2.2 at 2.4 GHz, a 0.03
wavelength air gap) finds its coupled modes at 0.99 and 1.37 f0 - nothing near
f0 but the first. Growing the parasitic lowers both (1.16: 0.95 / 1.20; 1.36:
0.87 / 1.13), but at this gap they stay 25-38% apart for every ratio tried: too
strongly coupled to straddle f0, so the default design is a single patch with
extra height, not a double-tuned one. The driven patch now takes the full-wave
length; the stack stays low confidence, and a two-layer solver with a feed
model is future work.

## BUILD STATE

S1-S5 done. 72 archetypes, 10 families, 2698 tests, 936/936 known cases.
Numerical arbiters now cover wires (MoM), dielectric resonators (FDTD), lenses and
zone plates (ray tracing, Kirchhoff), solid bodies of revolution (`bor`) and
printed patches on a grounded slab (`patch_sdm`, full-wave, checked by `patch_fdtd`) and
waveguide steps and stepped horns (`waveguide_step`, checked by `bor_fdtd`).
The cone family (discone, biconical, conical monopole) is solved end to end, and
fat dipoles to 0.025 wavelengths with a finite gap; patch bandwidths come from a
solved radiation Q, the shorted patch's with its wall current; horn beamwidths,
and the two-mode horns' efficiency and cross-polar levels, from aperture
integration; patch directivity from the cavity current through the slab, and
the rectangular patch's resonance, Q and surface waves full-wave - and it
builds the length that resonates, and the CP square takes its side and Q0
full-wave; the Potter horn's step is solved. Development runs to the finish line
above (hard stop 2026-09-27 20:00): item 13 is done; next, item 12 (the annular
ring mode in patch_q, and PIFA), then items 9, 5 and 4, then the wrap-up round
that tags v1.0.
See `docs/HANDOVER.md` for what is verified, what is not, and next steps.

The largest untested surface is unchanged: neither exporter has been run against a
real CST or HFSS installation. Structurally validated only.
