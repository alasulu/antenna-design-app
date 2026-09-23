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

## BUILD STATE

S1-S5 done. 72 archetypes, 10 families, 1781 tests, 479/479 known cases.
See `docs/HANDOVER.md` for what is verified, what is not, and next steps.

The largest untested surface is unchanged: neither exporter has been run against a
real CST or HFSS installation. Structurally validated only.
