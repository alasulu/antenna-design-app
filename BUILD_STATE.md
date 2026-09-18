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

## BUILD STATE

S1-S5 done. 72 archetypes, 10 families, 1521 tests, 426/426 known cases.
See `docs/HANDOVER.md` for what is verified, what is not, and next steps.

The largest untested surface is unchanged: neither exporter has been run against a
real CST or HFSS installation. Structurally validated only.
