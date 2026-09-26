# OTA Hub Antenna Toolkit — handover

An open reimplementation of the Antenna Magus workflow: state electrical
requirements, get a parameterised geometry, its predicted performance, and an
export to CST Studio or Ansys HFSS.

Built across five sessions on 2026-09-13/14. **72 archetypes, 10 families,
19,172 lines of Python, 23,984 lines of spec data, 2745 tests, 946/946 citable
known cases passing.**

---

## 1. What exists

| Area | Module | Contents |
|---|---|---|
| Engine | `otahub/core/` | Spec model, whitelisted AST evaluator, partial synthesis solver, registry, first-principles pattern maths |
| Catalogue | `specs/*.json` | 72 archetypes across wire (11), patch (9), loop (8), horn (8), travelling-wave (8), UWB (8), reflector (7), slot (6), lens (4), dielectric (3) |
| Reference solvers | `otahub/num/` | Thin-wire method of moments (EFIE, mixed potential, rooftop basis, Galerkin) with lumped loading, bent wires, junctions of any number of wires, computed VSWR bandwidth, feed networks (transmission lines coupled in admittance form), plus a Fourier-mode solution of the circular loop and a Hallén solver as a second opinion on straight wires. Independent full-wave checks on the closed forms |
| Arrays | `otahub/arrays/` | Uniform, binomial, Dolph-Chebyshev, Taylor n-bar, raised-cosine tapers; linear array factor, steering, grating-lobe limits; planar rectangular and triangular lattices with exact directivity, scan loss and beam-following cuts |
| Waveguides | `otahub/waveguides/` | Rectangular and circular guides, exact WR-series table, coax, microstrip, stripline, CPW |
| Utilities | `otahub/utils/` | S/Z/Y/ABCD conversion and cascading; L-section, quarter-wave and single-stub matching; Touchstone read/write and comparison against a prediction |
| Export | `otahub/export/` | Neutral geometry IR rendered to CST VBA and HFSS IronPython |
| Interfaces | `otahub/cli/`, `otahub/gui/` | 12 CLI subcommands; PySide6 GUI with catalogue, linear-array, planar-array and waveguide tabs |

```bash
python -m pytest -m "not slow"      # the quick loop, ~30 s
python OTA_Hub_AntennaToolkit.py gui
python OTA_Hub_AntennaToolkit.py --help
python -m pytest tests/ -q
```

---

## 2. The central design decision

**Archetypes are data, not code.** Each is a JSON document declaring its
parameters, synthesis formulas, analysis formulas, validity ranges and
citations. Formulas are strings evaluated by a whitelisted AST evaluator.

Three consequences, in order of importance:

1. **Every cited number is an executable test.** Each archetype's
   `known_cases` entries become pytest cases automatically. A formula that
   drifts from its reference fails the build rather than quietly producing a
   plausible, wrong antenna.
2. **Adding an antenna needs no Python.** A spec file is the whole
   contribution, and it gets a CLI entry and a GUI form for free.
3. **Nothing in a spec can execute.** Imports, attribute escapes,
   comprehensions and lambdas are rejected — enforced by test.

---

## 3. What is verified, and how

### Reproduced exactly from published worked examples

| Source | Result |
|---|---|
| Balanis Ex 14.1 | Patch W = 1.186 cm, L = 0.906 cm |
| Balanis Ex 14.2 | Inset edge resistance 228.35 Ω (to 1.7%) |
| Balanis Ex 14.4 | Circular patch a = 0.525 cm |
| Balanis Ch. 4 | Dipole R_r = 73.08 Ω, X = +42.52 Ω at λ/2 |
| Pozar Ex 5.1 | L-section C = 0.92 pF / L = 38.98 nH, C = 2.60 pF / L = 46.14 nH |
| Pozar Ex 5.2 | Stub tuner d = 0.110λ, l = 0.095λ; d = 0.259λ, l = 0.405λ |
| Viezbicke NBS TN688 | Yagi gains 7.1–14.2 dBi, fit residual ≤ 0.18 dB |
| WR-90 datasheet | Cutoff 6.557 GHz, 0.108 dB/m, 1.05 MW, band 8.2–12.4 GHz |

### A full-wave reference of its own

`otahub/num/mom.py` is a thin-wire method of moments — EFIE in mixed-potential
form, rooftop basis, Galerkin testing, the singular part of the kernel removed
analytically rather than quadratured. It exists because several archetypes
(one-wavelength loop, quad loop, folded dipole) had no closed form to check
against and no published number precise enough to cite, and because §6 item 4
— validate a low-confidence archetype end to end — was previously marked
impossible for want of a solver.

It is validated by six independent checks, in `tests/test_mom.py`:

| Check | Result |
|---|---|
| Matrix alone, driven with a prescribed sinusoid (no solve) | 73.083 + j42.494 Ω against the induced-EMF 73.0796 + j42.5152 |
| Far-field code alone, same prescribed current (no matrix) | Rr = 73.088 Ω, D = 1.6405 against 1.64093 |
| Pattern integral vs circuit power 0.5·Re(V I\*) | agree to 1 part in 10⁵ |
| Short dipole and small loop directivity | 1.4996 and 1.4971 against the exact 1.5 |
| Small loop Rr vs 20π²(C/λ)⁴ | ratio → 0.998 at C = 0.025λ |
| Folded dipole vs plain dipole at resonance | 4.014 against the classic 4 |

**Limits that bite:** no ground plane (image theory stands in), no
dielectric, no loss, junctions only at wire ends, a delta-gap feed only, and
the exact kernel only on straight runs and only on request (`exact=True`).
The module docstring records the measured edge of each: the delta gap reads
short-wire resistance ~10% low; the reduced kernel needs segments of at least
~3 radii, which the exact kernel lifts; and on a wire fatter than ~0.01 lambda
the delta gap itself has no converged answer with either kernel.

### A second one, for dielectric bodies

The MoM has no dielectric, so the DRA family had nothing full-wave behind it.
`otahub/num/dra.py` is a small 3-D FDTD solver for a resonator on a ground
plane: one quarter of the upper half space (x = 0 a magnetic wall, y = 0 and
the ground electric walls, which admit only the horizontal magnetic dipole
mode), CPML outside, a pulse to ring the resonator and a matrix pencil on the
ringdown for the complex frequency - resonance and Q in one number. Beside it
is the exact TE₁ pole of a dielectric sphere, which is the check it must pass:
a hemisphere 10-12 cells in radius lands on the pole to 0.15% in frequency
and 0.9% in Q at εr = 10 and 50, and the air margin and excitation move
nothing in the fourth figure. A grid-aligned brick has no staircasing at all;
refining it from 12 to 16 or 18 cells per height moves k₀d by 0.1% and Q by
0.2% at most; a staircased puck 16 cells in radius, by 0.16% and 0.3%.
`tests/data/` holds the 183 brick and 69 puck ringdowns the fits came from and
were checked against.

**Limits that bite:** 10-60 s a run in numpy, so fits are made offline and the
suite runs eight live checks; only modes with a horizontal magnetic dipole;
lossless dielectric; infinite ground; no feed. It reports the LOWEST mode with
real amplitude at the probe, because on long narrow bars a higher-order mode
can ring louder - the first survey took that mode on seven bars before the
selection rule was fixed. And it cannot see a mode with Q below about 2: the
ringdown is read from six pulse widths on, by which time such a mode has gone,
and a shorter pulse rings so many higher modes that they bury it instead. On a
flat low-permittivity puck (a/h = 6, εr <= 12) it therefore reports a higher
mode; `tests/data/cylindrical_dra_fdtd.json` keeps those runs, unfitted.

### A third, for solid surfaces of revolution

A wire cage cannot stand in for a solid cone: its impedance was still moving
3-5% between 24 and 32 wires. `otahub/num/bor.py` solves a rotationally
symmetric surface driven symmetrically - a cone, a disc, a tube, any polyline
generating curve off the axis - as the one-dimensional problem it is: EFIE in
mixed-potential form, rooftop basis on the TOTAL ring current along the curve,
Galerkin testing, and the wire kernel replaced by its exact ring averages. For
each azimuth a straight generating segment is a straight line in space, so the
static 1/R is integrated along it in closed form and the log singularity left in
azimuth is removed by phi = pi u^2. The feed is a FINITE gap, a uniform field
over a stretch of the curve, so a fat body has a physical feed instead of a delta
gap whose capacitance diverges.

Checked four ways that share nothing: a thin tube against the wire code's exact
kernel, which models the same tube (0.3%); input power against power radiated
through an independent far-field integral (1 part in 10^5); a long 47-degree
bicone against Schelkunoff's characteristic impedance (within 1%); and a cone
over a disc of 1.5, 3 and 5 wavelengths converging on image theory's half
bicone (7.7%, 3.7%, 2.1%), which is the check on the disc's purely radial
current.

**Limits that bite:** m = 0 only - symmetric excitation, so no tilted beams and
no cross-polar; perfectly conducting, zero-thickness sheets; a surface cannot
reach the axis (stop it at a small radius, where a feed tube goes).

### A fourth, for printed patches

`otahub/num/patch_sdm.py` is a spectral-domain method of moments for a rectangular
patch on an infinite grounded slab: the exact spectral Green's function from the
slab's transverse equivalent network, entire-domain currents in both directions
that meet the edge conditions EXACTLY (U_n sqrt(1-t^2) along the current, so it
vanishes as the square root of the distance to an open edge; T_j/sqrt(1-t^2)
across it), Galerkin testing, spectral integrals on a contour lifted over the
surface-wave poles, one quadrant by symmetry, and the slowly converging reactance
extrapolated from two truncations (every element converges as 1/k_max; the
extrapolation matches one from twice the truncation to 0.02%). The patch's mode is
the characteristic mode that radiates most; its resonance is where the eigenvalue
crosses zero. A first version with sines along the current - linear, not square
root, at the open ends - wandered +-1.5% in resonance with the basis; this one
moves 1e-4 from its default set to the next.

`otahub/num/patch_fdtd.py` checks it with nothing in common but the geometry: a
3-D FDTD ringdown of the patch as a PEC sheet on a slab that runs into the CPML
(`dra`'s Yee grid), at three cell sizes of an exactly representable geometry,
extrapolated at the first order the sheet edge gives (ratio of successive
differences 2.03 and 2.04). Resonance and total Q agree to 0.013% and 0.2% on a
thick eps_r 10.2 patch, 0.02% and 0.15% on an FR-4 one, and 0.007% and 0.17% on
a thick eps_r 10.2 square; doubling the air margin moves the FDTD not at all.

Also checked against answers it shares no code with: a vanishing current
element's space- and surface-wave power against Jackson and Alexopoulos's
thin-slab closed forms; an infinite microstrip's effective permittivity against
Kirschning and Jansen (0.1-0.9%); the radiated power from the spectral reaction
against a far-field integral (1 part in 10^4); and the patch's surface-wave
efficiency against Jackson's closed form, to 0.0-0.5% on thin board.

**Limits that bite:** rectangles only (the TM10 class), no feed model,
lossless, infinite substrate and ground; FDTD is affordable only on moderately
thick boards (the FR-4 check took 80 minutes at its finest grid). What the MoM
cannot hold, the FDTD can: `patch_fdtd.shorted_ringdown` runs a quarter-wave
patch with a full shorting wall in a half-space grid (checked against the
quarter-space grid on a symmetric patch: identical), and two sheets make a
stacked patch. On air the same grid is a PIFA: two plates put L + h at 0.2428 and
0.2423 of their resonant wavelength, 3.1-6.7% below the unshorted plate of twice
the length solved by the MoM.

### A fifth, for waveguide steps and flared horns

`otahub/num/waveguide_step.py` solves a step between coaxial circular guides by
mode matching over the m = 1 modes (TE1n and TM1n, in numbers proportional to
the radii), and cascades such steps with uniform sections through generalised
scattering matrices, so a stepped or flared horn is a staircase of them.
`otahub/num/bor_fdtd.py` checks it with nothing in common: an FDTD at azimuthal
order 1 on a (rho, z) grid, walls of any radius profile on grid lines, CPML ends,
a TE11 pulse in, TE11 and TM11 read out by projection. On a single step the two
agree to 0.06% in the TM11 power share and 0.04 degrees in its phase once the FDTD
is taken to zero cell size (order about 1.4, set by the step's corner). On a whole
10-wavelength Potter horn in the FDTD's own staircase, the cascade's aperture
projections sit where the FDTD is converging (0.367 at 37.4 degrees against
0.355 and 0.363 at 60 and 120 cells per wavelength).

**Limits that bite:** m = 1 only; perfectly conducting walls; the cascade models
a staircase (fine steps approximate a smooth cone); the FDTD's finest horn run
took 11 minutes.

### Independent first-principles checks

`otahub/core/pattern.py` is deliberately spec-independent, so it acts as an
external yardstick rather than restating what the specs claim: short dipole
D = 1.500000, half-wave D = 1.640922 and 78.08° HPBW, cos^q feed D = 2(2q+1),
uniform line source −13.263 dB sidelobe and 50.8 λ/L beamwidth.

Dolph-Chebyshev **measured** sidelobes equal their design level to 0.000 dB
across n = 4…21 and −20/−30/−40 dB. Uniform half-wave array directivity is
exactly N. Binomial arrays have no sidelobes at all.

---

### The dimensional audit — verification that needs no reference

`tests/test_scale_invariance.py` exploits the fact that Maxwell's equations have
no preferred length. Scale every frequency by S, every length by 1/S and
conductivity by S, and the antenna is electrically identical: directivity,
beamwidths, impedances and efficiency unchanged, lengths down by S, areas by S².
Every quantity must follow the power of S its declared unit implies.

**793 quantities across all 72 archetypes**, checked without a single reference
number. This is the counterweight to the `known_cases` harness, which can only
check what a human chose to assert — and which the same human wrote the formulas
for. A companion test asserts that every declared unit has a known scaling rule,
which is how unit typos surface.

It found, on specs that had been passing their own cited cases since session 1:

- `circular_patch` carried its intermediate radius `F` in **centimetres while
  declaring it dimensionless**, with factors of 100 threaded through the
  synthesis to compensate. It worked, and it violated the SI-internal rule the
  spec contract opens with. Rewritten in SI; Balanis Example 14.4 still
  reproduces exactly.
- `half_wave_slot` declared a slot area as `m2` rather than `m^2`, so nothing
  downstream — including the exporter's unit-driven classifier — recognised it
  as an area.

Injecting a deliberate dimensional error (a fixed 1 mm added to a dipole length)
confirms the audit fires, and reports the downstream impedance error at 89.6%
where the known case sees only 1.63% on the length.

---

### Results derived here rather than quoted

Session 5 added archetypes whose key numbers are not in any table I could cite,
so they were computed and then checked against an independent model before
being written into a spec.

| Archetype | Result | Independent check |
|---|---|---|
| `long_wire_travelling` | Rr = (η/2π)[γ + ln(2kl) − Ci(2kl) − 1 + sin(2kl)/(2kl)], integrated from the pattern | Matches quadrature at every length tested; reduces to 80π²(l/λ)² as l → 0 |
| `dipole_over_ground` | Zenith directivity from image theory with exact mutual impedance | Hemisphere integration agrees to five digits; mutual Z reproduces −12.5 −j29.9 Ω at d = λ/2 |
| `turnstile_dipole` | On-axis D equals a single dipole's; element plane exactly 3 dB down | Spherical integration of the summed-power pattern |
| `v_antenna_travelling`, `rhombic` | Axial directivity fitted to a four-leg travelling-wave model | 1.3% max fit error over 1.5–12 λ; axial lobe confirmed to be the global peak at the design angle |
| `corner_reflector_90`, `corner_reflector_60` | Image array factors; and now the driving-point impedance of the half-wave element, fitted over spacing and radius | Summed field leaves ~1e-15 tangential E on the plates. `otahub/num/corner.py` solves dipole plus images as one MoM problem: directivity within 0.04 dB of the specs everywhere tested, power balance to 1e-5, impedance fit within 2.2 ohm of fresh solves off its grid |
| `diagonal_horn` | Aperture efficiency 8/π² = 0.8106, and at any flare exactly, as the E-plane sectoral phase factor times the H-plane sectoral efficiency in Fresnel integrals; beamwidth and cross-polar lobes against phase error and size | Aperture integration on a 2001² grid: 0.8110. `otahub/num/horn_pattern.py` transforms both field components: the Fresnel product to 1e-6; beam fit 0.07% (0.04% held out), cross-polar fit 0.003 dB. The lobes are −15.5 dB, not the −19 the notes claimed, and the edge-centre/corner phase wording was a factor of 2 out |
| `annular_ring_patch` | Cubic correction to the narrow-ring rule | Bisection on the exact Bessel cross-product; 0.20% error against 2.71% uncorrected |
| `annular_ring_patch` | Directivity from the TWO edge walls' magnetic ring currents, fitted in (k₀a, k₀b) | Quadrature and an independent 2-D angular grid agree to four decimals; replaces a hard-coded 5.0 that was 5% low at εr = 2.2 and 48% high at εr = 10.2 |
| `annular_ring_patch` | TM11 radiation Q and substrate directivity factor, from `patch_q.annulus` through the slab | Meets the edge ring currents' Q and directivity to within 0.7% on thin board; fits 0.63% (Q) and 0.12% (factor) on held-out designs |
| `one_wavelength_circular_loop` | Resonant circumference, driving-point resistance and broadside directivity, all as functions of the conductor thickness | A method-of-moments solve and an independent Fourier-mode solution of the same loop agree to 0.25% on impedance and four decimals on the current; replaces a fixed 1.09λ / 100 Ω / 3.4 dBi that were 2.6–5% long, 28–30% low and 5–6% low |
| `quad_loop_square` | Resonant perimeter, resistance and directivity vs conductor thickness | Anchored to the circle's modal solution through an N-gon sequence (96-sided polygon 136.46 Ω against the modal circle's 137.07); replaces a 1.0218λ perimeter that is not resonant in free space at all |
| `halo_loop` | Resonant ring size, resistance, peak directivity and azimuth ripple | Bent-wire MoM anchored by continuity to the straight dipole; the ripple is 2.7–3.3 dB against an asserted 1.5, and the spec's conductor length subtracted the gap twice |
| `folded_dipole` | Resonant length, resistance, directivity and a COMPUTED VSWR bandwidth | Bandwidth found twice, directly and via the Yaghjian–Best antenna Q, agreeing to 2%; the 4:1 step-up confirmed to 3.949–4.021 between resonances, and shown NOT to hold at fixed frequency |
| Loop family bandwidth and Q | VSWR-2 windows and antenna Q for the full-wave circular loop, the quad and the halo | Three routes on the circular loop — MoM walk, MoM Q-derivative, and the independent modal solver's own walk — agreeing to under 1%; replaces four asserted Q values, three of which were 23–67% wrong |
| `alford_loop` | Equal-area circle equivalence, and the true ripple of a uniform-current square | Within 0.73% on resistance and 0.08% on directivity, far better than the “few percent” claimed; ideal ripple 0.01 dB against an asserted 0.5 |
| `resonant_dipole`, `half_wave_dipole` | Driving-point impedance, bandwidth, Q and directivity of the real wire, BESIDE the induced-EMF closed forms | Handing the MoM the assumed sinusoidal current reproduces the closed form to 0.3%, so the 6–28% gap is the current's shape, not the mesh; the canonical 73.08 Ω is kept and pinned as a definition |
| `axial_mode_helix` | Gain and beamwidth over (turns, C/λ, pitch), plus the full solved grid as a table | Image-theory helix, mirror-symmetric to 1e-7; directivity by reciprocity agrees with the pattern integral to 0.02 dB; replaces a “correction” that was worse than the formula it corrected |
| `lpda` | Directivity and its periodic ripple over one log-period; input resistance and the feeder that sets it | Every element solved, feeder as a transposed transmission line; the spec's Carrel reading holds to ±0.6 dB, and the missing feeder design relation is added and checked to within 11% |
| `rhombic` | Directivity, the optimal termination and the power it burns, solved with the resistor in circuit | The idealised travelling-wave model is 0.24–0.45 dB high; the termination is set by the wire radius (878 Ω at a = 10⁻⁵λ to 273 at 3×10⁻³), which the spec had left out |
| `top_loaded_monopole` | The top-current ratio derived from the hat geometry, which the spec used to take as an unrelated input | Junction-connected radial hat over image ground; fitted to 0.017 in beta (0.027 held out); the reactance ships as a solved table because no fit was good enough to size a coil |
| `biconical`, `conical_monopole` | Directivity at the band edge against flare, and across the band; the biconical's angle convention | 16-wire cages, directivity converged to 0.2%; the biconical used the full-angle impedance formula on a half-angle input (188 Ω for a cone that presents 100) |
| `half_wave_slot`, `folded_slot` | Resonant length and resonant resistance from the complementary dipole of radius w/4, through `resonant_dipole`'s laws: 0.4637λ and 468 Ω at the default w/L = 0.05, replacing a thin dipole's 0.4785λ and a flat-67-Ω 529.6 | Exact-kernel MoM and Hallén's equation, each on meshes it has converged on: 0.4645–0.4654λ and a 471–479 Ω slot; at w/L = 0.02 the law is inside 0.1% on length |
| `truncated_corner_cp_patch` | Corner cut, the frequency the square is sized for, mode split, and the axial-ratio and VSWR bands, all from a cavity model of the real outline | `otahub/num/patch_cavity.py`: Neumann modes on linear triangles, probe feed, broadside polarisation over 30 modes, converged to 0.004% between 240 and 480 cells per side; fits hold to 0.03% on meshes they never saw. The old design had 7.9 dB of axial ratio at its own f0 |
| `pyramidal_horn` | Two apex distances solved so both flares meet a real feed guide at one length (a quartic in √(ρ₁/λ)), with the optimum proportions in axial distances and the exact optimum efficiency | Flare geometry rebuilt from the output meets the guide to 1e-9 m; direct 2-D aperture integration puts the gain on target to 1e-3 dB from 12 to 30 dBi on three guides; root matches bisection on the unsquared condition |
| `waveguide_longitudinal_slot`, `waveguide_slot_array_resonant` | Stevenson's shunt conductance with the wavelength ratio the right way up, and the array offsets it sets | `otahub/num/waveguide_slot.py` derives g by reciprocity, a half-space radiation integral and power balance, all by quadrature: +0.03% (the rounding of 2.09) from 8.2 to 40 GHz on two guides; twelve slots at the new offset sum to unity by that model |
| Taylor taper (`otahub.arrays`) | Villeneuve's discrete n̄ distribution by zero placement, replacing a sampled line source that overshot its sidelobe level by up to 2.3 dB | Realised sidelobes by dense array-factor evaluation over 357 designs; nulls land where placed to 1e-12; equals the separate Dolph-Chebyshev implementation past the last zero pair to 1e-12; converges to Taylor's textbook line source (5e-6 at N = 400) |
| `prime_focus_parabolic` | Taper and spillover efficiency, beamwidth and first sidelobe from the edge taper and f/D (Silver's cos^n feed), replacing four separately asserted numbers | `otahub/num/paraboloid.py` integrates the aperture field directly: efficiencies to 1e-5, beamwidth to 0.03 lambda/D and first sidelobe to 0.1 dB off the fit grid; a cos^2 feed reproduces the classic 0.829 optimum at 66 deg |
| `offset_parabolic` | Taper and spillover efficiency and beamwidth from the feed taper and offset geometry; the gain comparison with a prime-focus dish, honestly signed | The rim-cone circularity that makes spillover exact checked ray by ray (1e-14 rad); taper and beamwidth by 2-D integration over the projected aperture, held-out errors 0.0013 and 0.03 lambda/D; the offset integrator reproduces the prime-focus one at zero offset to 1e-9 |
| `cylindrical_parabolic` | Taper and spillover under a line feed (cylindrical spreading), focusing-plane beamwidth and peak sidelobe, from the edge taper and f/W | `otahub/num/paraboloid.py` integrates over the aperture coordinate: efficiencies to 1e-5, beamwidth to 0.05 lambda/W and peak sidelobe to 0.15 dB off the fit grid; the uniform limit reproduces 50.8 lambda/W and -13.3 dB |
| `cassegrain`, `gregorian_dual_reflector` | Taper and spillover through the equivalent paraboloid, field-weighted subreflector blockage, beamwidth and peak sidelobe of the blocked aperture | The equivalent paraboloid traced ray by ray through the real hyperboloid and ellipsoid (magnification to 1e-9); efficiencies from the traced ray-tube mapping to 1e-5; blocked-aperture transform off the fit grid to 0.02 lambda/D and 0.3 dB |
| `hyperbolic_dielectric_lens`, `metal_plate_lens` | Taper and spillover efficiency, beamwidth and aperture edge illumination from the feed taper and the lens's own ray mapping; the metal-plate lens's fold-back limit acos(n) | `otahub/num/lens.py` traces rays with Snell's law at the actual face (parallel to 1e-9, on the closed-form mapping to 1e-12) and integrates the traced ray tubes: efficiencies to 2e-4, beamwidth to 0.3 lambda/D, edge illumination to 0.02 dB |
| Sectoral, conical and corrugated horns (beamwidths) | HPBW of each principal plane at any flare, from the aperture field with its quadratic phase error; NaN past the point where the beam breaks up | `otahub/num/horn_pattern.py` by quadrature, matching the E-plane pattern's Fresnel-integral closed form to 1e-10; fits to 0.5% (0.5% held out) over apertures of 1.5-40 wavelengths |
| `conical_horn_dual_mode` | Efficiency (exact, as a ratio of four mode integrals), both beamwidths and the peak cross-polar level against the TM11 power fraction; the fractions that equalise the beams (0.097) and null the cross-polar field (0.13) | `otahub/num/horn_pattern.py` with TE11 + TM11 in phase at the aperture: its TE11 limit reproduces the scalar routine to 1e-9 and the textbook 0.837; beam fits to 0.002%, cross-polar table to 0.03 dB on held-out sizes. Efficiency 0.506 against an asserted 0.62, so the "+0.85 dB over a smooth horn" is a 0.27 dB loss; beams 78-82 lambda/D, not 68. The step is now solved by mode matching (`waveguide_step`, FDTD-checked by `bor_fdtd` to 0.06%); the phasing length and cross-polar bandwidth are first-order - the flare's mode conversion moves the aperture phase ~20 degrees on a 10-wavelength horn (stepped-cone cascade and whole-horn FDTD agree), and the -30 dB band is about 2% there, not 6% |
| Rectangular, circular, triangular and corner-truncated CP patches | Radiation Q and VSWR-2 bandwidth, each shape from its own cavity mode; the CP patch's Q0 and therefore its cut | `otahub/num/patch_q.py`: stored energy of the mode against the space wave of the patch current on its grounded substrate. Reproduces Jackson's c1 to 0.02% for a vanishing patch, and on air the cavity model's edge-current directivities to 0.06%; fits to 0.6% (0.25% held out) |
| Rectangular, inset, CP, circular, triangular and shorted patches (directivity) | Broadside directivity as the thin-substrate value (all walls radiating) times a substrate factor, NaN past h sqrt(eps_r)/lambda0 = 0.1 | `otahub/num/patch_sdm.py` full-wave: `patch_q`'s directivity within 0.3-0.9% on air, eps_r 2.2, FR-4 and 10.2. The two-slot formula was 7.5-13% low (side walls and substrate); the circular and triangular free-space edge models up to 18% low on thick board; the shorted patch's single slot 5-35% high |
| `rectangular_patch`, `rectangular_patch_inset` (resonance, Q, surface waves) | L is now the length that resonates at f0 full-wave (the textbook design kept as L_textbook, where it resonates reported); full-wave radiation Q, surface-wave efficiency and a bandwidth counting it; the inset depth is the transmission-line fraction of the full-wave L | `patch_sdm` against `patch_fdtd` (0.013-0.02% in frequency, 0.15-0.2% in Q). The Hammerstad design resonates 0.5-7% LOW (the classic 2.4 GHz FR-4 patch at 2.33 GHz); the cavity-current Q ran 1.6-28% high; efficiency matches Jackson's thin-slab form to 0.5% |
| `truncated_corner_cp_patch` (square side, Q0) | Q0 = full-wave radiation Q times surface-wave efficiency; the square side resonant at f_sq full-wave (textbook kept as L_textbook); the cut-to-split relation stays the cavity model's | `patch_sdm` on squares, `patch_fdtd` on one (0.007%, 0.17%). The old cavity-current radiation Q0 was 1-99% above the total - the surface wave it left out dominates on thick high-permittivity board - so the cut was too small; the textbook square 0.5-7% too large |
| `quarter_wave_shorted_patch`, `stacked_patch` (lengths) | Both now build the rectangular patch's full-wave length ratio (textbook kept as L_textbook) | FDTD spot checks: the textbook shorted patch resonates 4-8% low, more than the full patch (the wall's inductance); the corrected one still 1-3% low. The default stack's two modes sit at 0.99 and 1.37 f0 and stay 25-38% apart at any size ratio at its 0.03-wavelength gap; the spec's 'make the parasitic smaller' was backwards |
| `quarter_wave_shorted_patch` | Radiation Q and VSWR-2 bandwidth with the shorting wall's vertical current; the bandwidth relative to the full patch on the same board | `otahub/num/patch_q.py`: the wall's slab factor against a plane-wave boundary-value solve (1e-8); on air the shorted patch against the cavity model's magnetic currents on its three open walls, written independently (0.004%). Not half the full patch's bandwidth: 1.84 times it on air, equal near eps_r 2, 0.57 at eps_r 12. The side walls radiate a third of the power the old one-slot picture left out |
| `biconical` | Cutoff slant (VSWR 2 against its own Zc), continuous bandwidth, worst in-band VSWR and directivity at the cutoff, tabulated over 5-65 deg | `otahub/num/bor.py`, cutoff extrapolated from two meshes; the spec's own dimensions solved live sit at VSWR 2 at f0; the old wire-cage directivity sits 0.7-2.8% above the solid cone at a quarter wave |
| `conical_monopole` | The same in 50 ohm, from the bicone by image theory, over 15-65 deg with extra nodes where the cutoff climbs steeply (32-34 deg) | As above; at 47 deg, where the monopole's Zc is 50 ohm, the two specs' cutoffs coincide and image theory holds to 5e-4 |
| `discone` | Low cutoff (slant at VSWR 2 in 50 ohm), continuous VSWR-2 bandwidth, worst in-band VSWR and directivity at f_low, tabulated over half angle 20-50 deg and disc ratio 0.6-0.9 with a stated coax-sized feed | `otahub/num/bor.py`, the exact axisymmetric surface solution: low cutoff extrapolated from two meshes (0.4%), 0.9% on six held-out designs; the spec's own dimensions, rebuilt and solved live, sit at VSWR 2 at f_low |
| `fresnel_zone_plate` | Aperture efficiency and gain from a cos^n feed through the plate's own zones (Kirchhoff, in-spec feed-angle integral); spillover; HPBW and 1 dB gain bandwidth fitted per plate type and, for the opaque plate, per parity of M | `otahub/num/zone_plate.py` integrates on a Cartesian grid over the plate: efficiency to 0.015%, the fits to 0.18% (beam) and 0.30% (bandwidth) on held-out 2-D runs; many zones under uniform light recover 1/π², 4/π², 8/π² |
| `hemispherical_dra` | k₀a and Q from the exact complex TE₁ pole of the equivalent sphere, fitted to 0.23% over εr 6–100 | Newton on the characteristic equation; `otahub/num/dra.py`'s FDTD ringdown of a staircased hemisphere lands on the pole to 0.15% in frequency and 0.9% in Q; Mongia & Bhartia's published fit is within 1.04%. The peak of the Mie coefficient b₁ on the real axis, used before, was 1.2% high and its half-power Q 12% low at εr = 10 |
| `cylindrical_dra` | k₀a and radiation Q of the HE₁₁ mode over εr 6–50 and a/h 0.4–4 | Fitted to 48 FDTD ringdowns (`otahub/num/dra.py`): 0.34% in k₀a, 0.42% in Q; 0.16% / 0.49% on 12 held out. Mongia & Bhartia's published fits kept as a comparison: 9.4% low to 8.8% high in resonance, 26% low to 16% high in Q |
| `rectangular_dra` | k₀d and radiation Q of the broadside mode (magnetic dipole along L), the other mode's frequency | Fitted to 150 FDTD ringdowns (`otahub/num/dra.py`) over εr 6–50 and both aspects 1–3: 0.15% in k₀d, 0.7% in Q, and within 0.13% / 0.65% on 33 held-out ones. The DWM, run the right way round, is kept as a reported comparison: 10.6% high to 6.2% low |
| `rectangular_patch_inset` | G1 exact in Si(X); G12/G1 fitted to eq. 14-18a | Quadrature reproduces Balanis Example 14.2 to 0.05% |
| Sectoral and pyramidal horns | Aperture efficiency exact in Fresnel integrals | Closed form, Balanis 13-19/13-41 as published, and direct aperture integration all agree to 1e-14 |
| `circular_patch` | Directivity integrated from the TM110 fields in J0 ∓ J2 | Quadrature and an independent 2-D angular grid agree to 0.0000%; replaces a hard-coded 6.3 that was 84% high on εr = 10.2 |
| Conical and corrugated horns | Efficiency integrated from the TE₁₁ and HE₁₁ aperture fields | Uniform-phase limits reproduce the published 0.836 and 0.69; the latter derived from a J₀ taper rather than quoted |
| Rectangular patches (3 specs) | Directivity from the two-slot model, `2·D₁/(1 + G₁₂/G₁)` | Matches direct 2-D pattern integration to 0.000%; the narrow-slot limit gives exactly 3.0, a magnetic dipole doubled by the ground plane |
| `triangular_patch` | Directivity from the triangle's OWN Neumann eigenfunction, three radiating walls | D → 3.0000039 as the patch shrinks — a horizontal magnetic dipole over a ground plane — which nothing in the derivation was arranged to produce; replaces a hard-coded 5.0 that was 43% low on air and 47% high on εr = 10.2 |

**These fits are only as good as the model behind them.** Each is a
closed-form or ray-optics idealisation, not a full-wave result, and the
validity block on each archetype says where it stops.

---

## 4. What is NOT verified — read this before trusting a number

### Eight archetypes are marked low confidence

`cassegrain`, `conical_horn_dual_mode`, `ferrite_rod_loop`, `pifa`,
`planar_monopole_circular`, `stacked_patch`, `vivaldi_tsa`,
`waveguide_longitudinal_slot`.

**`halo_loop` was promoted out of this list**, which is §6 item 4 discharged
for the first time. Every number it asserts is now solved rather than assumed,
by a method of moments whose bent-wire path is anchored by continuity to the
straight dipole. The remaining eight cannot follow it yet for a reason that is
the same in every case: they are patches, horns, reflectors or slots, and a
thin-wire solver cannot reach them. `ferrite_rod_loop` additionally needs a
ferrite material model.

Two of those are new in session 5 and both are honest about why:
`stacked_patch`'s bandwidth multiplier is an expectation drawn from published
designs rather than a computed result. `conical_horn_dual_mode` is most of the
way out: its aperture is derived for a given TM11 share (0.62 was 0.506, the gain
"advantage" a loss), and its step is now SOLVED by mode matching and checked by
FDTD. What keeps it low is the phasing: the flare converts between TE11 and TM11
enough to move the aperture phase about 20 degrees on a 10-wavelength horn, which
a first-order phasing length cannot see, so that length and the bandwidth are
first-order and labelled so. The stepped-cone cascade that does see it exists
(§3); a spec-level design built on it is the way out.

`pifa` stays in the list with its length now FDTD-calibrated (§3): what keeps it
there is that a real PIFA's bandwidth and pattern belong to its ground plane,
which neither this toolkit nor the spec models.

They announce themselves in `list`, `show`, the GUI, and in every design they
produce. Treat their numbers as indicative and verify in a full-wave solver.

### Known open discrepancies

- ~~**`half_wave_slot` resonates at a thin dipole's length.**~~ **Resolved.**
  By Babinet a slot resonates where its complementary dipole does, and that
  dipole's radius is a quarter of the slot width. The slot carried a thin
  dipole's 0.4785 lambda and a flat 67 ohm whatever its width; it now takes
  both from `resonant_dipole` at radius w/4 - 0.4637 lambda and a 468 ohm
  resonant slot at the default w/L = 0.05, against the old 529.6, which was
  13% high. `folded_slot` follows, so its two-slot figure drops from 132 to
  117 ohm. Checking the complementary dipole needed an exact kernel first:
  at that radius the reduced kernel's answer depended on the mesh. Exact-kernel
  MoM and Hallen's equation agree with the laws to 0.4% on length and 2.5% on
  resistance. Past w/L = 0.05 it is extrapolation, and the spec says so.

- ~~**Inset patch mutual conductance.**~~ **Resolved — and the discrepancy was
  never real.** Direct quadrature of eqs. 14-12 and 14-18a reproduces Balanis
  Example 14.2 to 0.05%: G1 = 1.5747e-3 S against the published 1.57e-3, and
  G12 = 6.1651e-4 against 6.1683e-4. The recorded "212.5 Ω by direct
  integration" was a bad integration, not a disagreement with the textbook.

  The spec's old 1.7% agreement was two errors cancelling: `G1` used the
  `(1/90)(W/λ₀)²` small-width branch, which runs 10% high, while
  `G12 = G1·J0(k0·L)` runs 28% low. `G1` is now the exact closed form
  (`I₁ = −2 + cos X + X·Si(X) + sin X/X`, X = k0W, verified against quadrature
  to eight figures), and `G12/G1` is a polynomial fitted to the exact integral
  over the whole practical design space — k0L ∈ [0.8, 3.3], k0W ∈ [1.15, 3.25],
  which is what εr from 1 to 12 produces. Worst fit error 0.0021 in the ratio,
  under 0.1% in the resistance. Example 14.2 now reproduces to **0.066%**.
- ~~**Axial-mode helix gain.**~~ **Resolved** by a full-wave solve (§3): the
  shipped gain is fitted to the solved helix, within 0.35 dB over the design
  core, and the Kraus formula is kept beside it with its measured error.
- ~~**Stevenson's g1**~~ **Resolved, and it was wrong.** The spec's
  expression carried the wavelength ratio inverted, lambda/lambda_g, while its
  own note quoted lambda_g/lambda. A first-principles derivation by quadrature
  (`otahub/num/waveguide_slot.py`: reciprocity for the TE10 the slot excites,
  its one-sided radiation over a half space, power balance for a shunt element)
  reproduces lambda_g/lambda to the rounding of the 2.09 across WR-90 and WR-28.
  The inverted form was 28-64% low; the resonant array's WR-90 offsets were
  3.06 mm where 2.28 mm is right, and the old ones summed to a conductance of
  1.76 - VSWR 1.76 at a feed meant to be matched. Stevenson's model itself
  (half-wave cosine slot, thin wall) is still not checked against a real guide.
- **LPDA directivity** collapses a two-dimensional Carrel chart onto the
  optimum-σ line. Checked since against a full-wave solve of every element: the
  reading holds to ±0.6 dB of the log-period mean, and the periodic ripple
  around that mean is reported separately.
- ~~**Taylor taper**~~ **Resolved - and the recorded example was wrong.** The
  "-28.9 dB for a 20-element -30 dB design" does not reproduce: that taper gave
  -30.10 dB. Its real failures were elsewhere and larger - the sampled
  line-source distribution overshot the design sidelobe level by more than
  0.5 dB in 156 of 357 designs, by up to 2.3 dB. `taylor_nbar` is now
  Villeneuve's discrete distribution, built by placing the array polynomial's
  zeros: within 0.05 dB above design for N >= 10 (0.36 dB in one 9-element
  case), the uniform array's nulls exactly from the nbar-th on, exactly
  Dolph-Chebyshev past the last zero pair, and the line source again as N grows.
- **The DRA family models the isolated resonator.** All three shapes are now
  backed by an exact pole or by FDTD ringdowns, so what remains unverified is
  what none of them models: the feed (probe, slot or microstrip, which pulls
  the resonance a few percent and sets the match), a finite ground plane, and
  dielectric loss. The cylinder's fits stop at a/h = 4: flatter, low-εr pucks
  have Q under 2, which the ringdown cannot isolate.
- **The cone family is solved for one feed.** Discone, biconical and conical
  monopole are now solved as bodies of revolution (§3), each with a stated
  coax-sized feed (cone top 2% of the base diameter). Their band edges move with
  the feed, so a very different feed is outside what was solved, and the conical
  monopole assumes an infinite ground plane.
- **The Fresnel zone plate is scalar Kirchhoff diffraction through an
  infinitely thin plate.** Its gain, beam and bandwidth now come from the feed
  through the plate, but a real phase plate's thickness, its shadowing at
  oblique incidence and polarisation are not modelled.

### Metrics labelled "indicative"

Bandwidth used to be the largest single group of these. Four loop archetypes
reached it only through an asserted Q; three are now solved, and the fourth,
`alford_loop`, keeps its assumption **and says so in the metric's own notes**,
because its bandwidth is set by a corner-loading network the spec does not
describe. That is the pattern to follow: derive it, or state plainly that it
cannot be derived and why — never quietly fit something that looks derived.

Several archetypes carry bandwidth, front-to-back or directivity figures that
are engineering rules of thumb rather than derived results. Every one says so
in its `notes`. They exist so a design sheet is complete, not because they are
predictions.

### Exporters build geometry for 43 of 72 archetypes

Dipoles (3), monopole, folded dipole, dipole over ground, turnstile, four
patch variants (seven now), six loops, six slots, three dielectric resonators,
two planar monopoles, a terminated long wire, a leaky-wave line source,
biconical, conical monopole, discone, and open-ended waveguide. **The other 35
export parameters only and say so.** By family: dielectric 3/3, slot 6/6,
patch 7/9, loop 7/8, wire 11/11, uwb 4/8, travelling-wave 2/8, horn 1/8, lens 2/4, and
nothing yet for reflector.
Where a builder must choose something the spec cannot supply — feed gap, inset
notch width, a finite ground plane standing in for an infinite one — the
choice is stated in the exported file's header.

Neither exporter has been run against a real CST or HFSS installation. The
scripts are structurally validated (balanced blocks, parseable Python, correct
units) but **not execution-tested**. That is the single largest untested
surface in the project.

---

## 5. Engine bugs found during the build

**A supplied value was silently overwritten by the spec's own default.** Many
archetypes synthesise a nominal geometry — a halo's gap, a loop's circumference
— so they are usable from a frequency alone. The resolver ran those rules
unconditionally, so a caller who supplied the value watched it be replaced and
got a design sheet for a different antenna. Nothing warned: every number on the
sheet was self-consistent, it just answered a question nobody had asked. A
supplied value is a requirement and now outranks the default. Found while
sweeping a halo's gap, which refused to move.

Recorded because they show what the test harness is for.

| Bug | Consequence |
|---|---|
| All-or-nothing synthesis | 79 of 100 known cases failed because one unresolvable optional rule discarded the whole design |
| Relative error against expected zero | A 3 Ω residual reactance scored as a 348% failure; `tol_abs` added |
| `db10()` used where `log10()` was meant | A 4.2λ Yagi reported **53 dBi** — more gain than a 30 m dish |
| Inverted `sqrt(RL/Z0)` in the L-section | Matched to 181 Ω instead of 100 Ω |
| Length guessed from magnitude on export | A 319 Ω resistance written as `319105 mm` |
| Chebyshev even-length phase term | Weights came out independent of the design sidelobe level |
| Stripline/CPW elliptic ratio applied twice | Wrong impedance; limits settled the correct orientation each way |
| Catalogue filter delegated parent rows | Every family survived every search |

Spec errors the harness caught include a loop loss resistance out by **30.9×**,
another by 3.3×, and several of my own arithmetic slips (WR-90 open-ended gain,
Ruze at 100 GHz, DRA Q scaling, biconical Z_c).

### A CP patch that would have radiated linear polarisation

`truncated_corner_cp_patch` asserted 3 of the 16 quantities it produced, the
lowest coverage in the catalogue, and a cavity model of its real outline found
three errors, any one of which would have spoiled a built antenna:

- **The feed rule was the wrong antenna's.** It said to feed on a diagonal.
  That is the NEARLY-SQUARE patch's rule; cutting a square's corners splits it
  into DIAGONAL modes, and on a diagonal one of them is identically zero. The
  cavity model gives an axial ratio of about 3000 dB there - linear
  polarisation. The feed belongs on a centreline.
- **The mode split was half what CP needs**, f0/(2Q): that is each mode's
  offset from the centre, not their separation.
- **The square was sized for the wrong frequency.** Cutting the corners lifts
  one mode and leaves the other, so the CP centre sits about 0.5/Q ABOVE the
  uncut square's resonance. Sizing the square for f0 put the CP centre 0.79%
  high in the spec's own example, outside its ±0.29% axial-ratio band: 7.9 dB
  of axial ratio at f0.

Also corrected: the classic cut dS/S = 1/(2Q) is first-order and under-cuts by
up to 2.6% in c; the impedance band of two staggered modes is twice a single
mode's; and the axial-ratio band, 0.347/Q exactly in the two-mode model, was an
"indicative" figure 29% low. Every one was checked in the cavity model, which
assumes none of the rules it checks.

### A dish budget that counted spillover twice

`prime_focus_parabolic` asserted illumination efficiency 0.82, spillover 0.85,
beamwidth 70 lambda/D and first sidelobe -24 dB, separately - though all four
follow from one feed pattern and f/D. Worse, 0.82 is the classic optimum of the
PRODUCT of taper and spillover, so multiplying it by 0.85 counted spillover
twice: 0.74 dB of gain lost at the default design. At f/D = 0.4 and -11 dB edge
taper Silver's cos^n model gives taper 0.886 and spillover 0.933 (product
0.827), beamwidth 66.5 lambda/D and first sidelobe -25.2 dB. All four are now
derived from the edge taper, and a cos^2 feed reproduces the textbook 0.829.

### A horn that could not be fitted to its waveguide

`pyramidal_horn` produced a1, b1 and one apex distance, and no feed guide at
all. One apex for both planes fits only a guide with the aperture's own
sqrt(3/2) aspect ratio; standard guides are about 2:1, so on WR-90 at 20 dBi its
two flares reached the guide 15 mm apart. Its note also called that distance a
slant length while every formula used it as an axial one, and blamed the
single apex for a gain offset that was really the 0.51 sizing efficiency.

The design now takes the guide (default WR-90 proportions scaled to lambda)
and solves for two apex distances so the flares meet it together. Written in
axial distances with the exact optimum efficiency, the condition is a quartic
and the gain lands on target. The textbook procedure writes the same condition
in slant lengths with an implied 0.5105 efficiency; integrated over the
aperture, its horns fall 0.04-0.19 dB short at 20-25 dBi and 0.8 dB at 15 dBi.
`np.roots` joined the expression whitelist to solve the quartic in the spec.

### The reference that was not what it looked like

73.08 + j42.52 Ω is quoted everywhere as "the" half-wave dipole impedance, and
the MoM appeared to be 18% wrong against it — 86 Ω at a/λ = 0.001, drifting
upward under mesh refinement, which looks exactly like a convergence bug. It is
not one. That figure is the **induced-EMF** value for an *assumed sinusoidal
current* on a *vanishingly thin* wire. The driving-point impedance of a
delta-gap-fed wire of finite radius is a different quantity, and the real
current is visibly fatter than a sinusoid near the ends.

Three things settled it, and it is worth recording which, because the first two
were not enough on their own:

1. driving the matrix with the sinusoid reproduced 73.083 + j42.494 — so the
   matrix was right and the solve was the suspect;
2. raising the quadrature from 24 points to 160 changed the answer by nothing
   at all — so it was not a quadrature error either;
3. an independent **Hallén** solve — different integral equation, no divergence
   term, collocation instead of Galerkin — returned 86.6 Ω. Two formulations
   sharing no algebra agreed with each other and disagreed with the constant.

`test_delta_gap_impedance_is_not_the_induced_emf_value` pins this so the next
reader does not spend the same afternoon on it.

That lesson has now been carried back to the two archetypes it came from.
`half_wave_dipole` keeps 73.08 + j42.52 exactly, because it is the definition
dBd rests on, and gains driving-point metrics beside it. `resonant_dipole` had
been presenting the closed form as what the antenna presents — "55-68 ohm
across practical thicknesses" — when a real resonant wire sits at about 72 ohm
whatever its gauge. Its 50 ohm VSWR moved from about 1.26:1 to 1.46:1.

### Two cross-checks written to agree instead of to check

`biconical` documented `theta_h` as the half angle and computed its impedance
with cot(theta_h/4) - Kraus's formula for the FULL cone angle - so it reported
the impedance of a cone half as wide: 188 ohm for the 47 degree cone that
presents 100, 243 for a 30 degree cone that presents 158. Its siblings
`conical_monopole` and `discone` used the half angle correctly. Two
cross-consistency tests existed to catch exactly this kind of disagreement,
and both had been written with a doubled angle on the biconical side - 10
against 5, 60 against 30 - which made them pass. A cross-check is only a check
if it compares like with like; these compared whatever made the numbers match.

### The helix "correction" went the wrong way

`axial_mode_helix` knew Kraus's gain formula overestimates and told designers to
use `gain_corrected_dbi` instead: 8.3 + 10 log10(C^2 N S) + 20 log10(1 + N/10).
Its note called that "sub-linear in N". It is super-linear, and it exceeds
Kraus for every N >= 5 — by 2.5 dB at ten turns, 6 at twenty, 8.5 at thirty,
rating a 30-turn helix at 28.75 dBi — in exactly the regime where Kraus is
already too high. That needed arithmetic, not a solver. The metric is removed.

The replacement is solved: a helix over an infinite ground plane by image
theory (`mom.helix_over_ground`), cross-checked by receive-mode reciprocity to
0.02 dB. At C = lambda and 13 degrees Kraus is 1.6 dB high at 3 turns, 3.8 at
10 and 4.4 at 20. Kraus's gain and beamwidth stay, under `_kraus` names, as the
figures everyone quotes.

### A ten-fold error in the loaded whip, from a missing referral

Both induced-EMF closed forms, Balanis 4-70 and 4-79, are referred to the
current MAXIMUM of the assumed sinusoid. For any length other than lambda/2 that
point is not the feed, and the feed-point value is the current-maximum one
divided by sin^2(kL/2). `dipole_arbitrary_length` applied that to its
resistance and not to the reactance on the next line, so its R and X described
different points on one antenna. `inductively_loaded_monopole` copied the
reactance without it, and for a short whip that is catastrophic: at
h = 0.05 lambda, sin^2(kL/2) = 0.0955.

| Quantity, 1.5 m whip at 10 MHz | Was | Is |
|---|---|---|
| Unloaded reactance | -j111 Ω | -j1163 Ω (MoM -j1140) |
| Loading coil | 1.77 µH | 18.5 µH |
| Coil loss, Q = 200 | 0.555 Ω, "over half" Rr | 5.81 Ω, 5.9 × Rr |
| Radiation efficiency | 61% (−2.1 dB) | 14% (−8.4 dB) |

A coil wound to the old figure would have left the whip still strongly
capacitive, about a hundred half-bandwidths from resonance. The half-wave
dipole and quarter-wave monopole were never affected, because at exactly
lambda/2 the referral is the identity — which is presumably how it went
unnoticed.

### Errors found in already-shipped specs during session 5

These had been passing their own tests since session 1, which is the point: a
known case only checks what it asserts.

| Spec | Error |
|---|---|
| `corner_reflector_90` | Array factor `4·sin²(kS)` put a **null at S = λ/2, where the optimum actually is**. Image theory gives `2[1−cos(kS)]`; the null is at S = λ. It also used a field factor as a power factor and ignored the mutual coupling that sets the feed resistance. Rebuilt as a four-element image array; it now returns the textbook ~12 dBi |
| `short_dipole` | `20π²(L/λ)²` was labelled the uniform-current value. It is the triangular one. The second metric then quartered an already-triangular value, **under-reporting a real short dipole by 4×** |
| `cassegrain` | `magnification` was a hard-coded `1.0` placeholder |

### Test-guard and CLI bugs found in session 5

| Bug | Consequence |
|---|---|
| `_bounds_for` matched substring before suffix | `efficiency_db` was bounded to [0,1], failing any negative decibel value |
| dBi floor rejected genuine nulls | A dipole λ/2 over ground has an exact zenith null and reported −310 dBi |
| Degree bound assumed unsigned angles | A beam angle measured from broadside is negative when it scans the other way |
| `array --sll 30` raised | The commoner spelling of "30 dB sidelobes" hit an unhandled `ValueError`; either sign is now accepted |
| `np.where` results stayed 0-d arrays | A spec output guarded with `np.where` came back as an array, not a number: the exporter refused it ("has none of ('L',)"), and the plausibility test skipped it - so metrics that return NaN outside their fitted domain had never been seen by it. The evaluator now unwraps 0-d results; the plausibility test lets deliberate NaN through, since an expected NaN is caught by the reference match |

---

## 6. Recommended next steps

> **Finish line (set by the user, 2026-09-26).** Scope is frozen to the items open
> here as of commit 53df1ab; development stops at 2026-09-27 20:00 local time,
> followed by a wrap-up round that rewrites this section as the final future-work
> list and tags v1.0. The full terms are in BUILD_STATE.md, section "Finish line".
> Anything newly found goes into this section as future work, not a new round.

1. **Execution-test the exporters** against real CST and HFSS installations.
   This is the biggest gap, and session 5 did not touch it — the catalogue grew
   by 32 archetypes while the exporter still builds geometry for 7. Start with
   `half_wave_dipole` (simplest geometry, strongest analytical reference:
   73.08 + j42.52 Ω at λ/2).
2. ~~**Close the inset-patch discrepancy**~~ — done. The published figure was
   right; the "direct integration" that disagreed with it was wrong. See §4.
3. **Add geometry builders** for horns, which need a loft or truncated-pyramid
   primitive — the one shape where I would be guessing at the API's structure
   rather than just its parameter names, so it has been left alone. The
   pyramidal horn's geometry is at least fully determined now: guide, aperture
   and a length at which both flares meet the guide. The loop
   family is now done (a torus turned out to be a direct primitive in both
   tools). Yagi-Uda remains blocked for a different reason: the
   spec gives boom length, reflector and driven lengths and a director count,
   but not individual director lengths, so its geometry is genuinely
   underdetermined and building it would mean inventing dimensions.
4. **Validate a low-confidence archetype** end to end and either promote it or
   record why it cannot be.
5. ~~**Planar arrays**~~ — done. Rectangular and triangular lattices, separable
   tapers, steering and exact directivity are in `otahub/arrays/planar.py`,
   with a `planar` CLI subcommand. Still missing: circular and thinned
   layouts, subarray architectures, and element-pattern embedding (the module
   assumes isotropic elements, so real gains need the element pattern folded
   in separately).
6. ~~**Touchstone import**~~ — done. `otahub/utils/touchstone.py` plus a
   `touchstone` CLI subcommand: reads v1.0 and v1.1 files in MA/DB/RI, any
   frequency unit, S and Z parameters, with the two-port column-major exception
   and arbitrary line wrapping handled. G and H parameter files are recognised
   and refused rather than mis-converted.
7. **Re-audit the session 1–4 specs** — continuing. A survey by how much each
   archetype actually asserts put the sectoral horns at the bottom (one known
   case asserting one quantity out of four produced); they have since been
   rebuilt on exact theory. The same survey is the way to pick what to look at
   next. Partly done otherwise: Two systematic audits now
   run over the whole catalogue: the dimensional one (793 quantities against
   scale invariance) and the cross-consistency one (23 relationships between
   archetypes that must agree). Both are clean. What neither covers is an
   archetype with no sibling and no dimensional quirk — those still need
   reading against their source.
   Every archetype passes the cases it declares; that is not the same as being
   right. `corner_reflector_90` had a null where its optimum is and said so
   confidently for four sessions. The archetypes carrying a single known case
   are the place to start. Ranked by quantities asserted over quantities
   produced, the bottom was `truncated_corner_cp_patch` (3 of 16), and auditing
   it found three design-breaking errors (§5). `pyramidal_horn` was next and
   turned out not to be buildable (§5). The ranking undercounts archetypes
   already verified another way (`annular_ring_patch`, the inset patch); of the
   ones it rightly flags, the waveguide slot family was next: Stevenson's g1
   turned out inverted in the spec (§4). The corner reflectors followed: their
   directivity held to 0.04 dB against an exact image solution, but the quoted
   "input resistance" was an induced-EMF figure the feed never sees. That leaves
   `stacked_patch`, which needs a full-wave solve this toolkit does not have.
   After it, the "indicative" list: `prime_focus_parabolic` was done from it
   (§5), and `offset_parabolic` after it - which needed its own 2-D
   integration, because the tilted spreading loss makes the symmetric formula
   up to 20% optimistic. `cylindrical_parabolic` followed, with a line feed's
   cylindrical wave, and the dual reflectors after it, through an equivalent
   paraboloid verified by ray tracing. That closes the reflector family's
   asserted efficiency pairs. The hyperbolic and metal-plate lenses followed,
   through a Snell's-law ray trace, and the Fresnel zone plate through Kirchhoff
   diffraction (its gain had left the feed out altogether). `luneburg_lens`
   keeps its 0.65 on purpose: in geometric optics a cos(theta) feed gives it a
   perfectly uniform aperture, so its shortfall is construction - stepped shells, material loss - which ray
   optics cannot see and which needs measured or full-wave data.
8. ~~**Replace `rectangular_dra`'s borrowed Q**~~ — done, and the resonance
   with it: the DWM had been solved along the height, which is the mode of a
   VERTICAL magnetic dipole a ground plane shorts out. Deriving Q from the DWM
   fields was tried and rejected: against a 3-D FDTD ringdown
   (`otahub/num/dra.py`, itself checked on the exact hemisphere pole) the DWM
   Q runs from 33% low to 26% high, and the DWM resonance from 10.6% high to
   6.2% low. Both are now fitted to 150 ringdowns. `cylindrical_dra` followed:
   its published Mongia & Bhartia fits are 6-9% low in resonance on every
   puck with a/h <= 0.5 and 6-26% low in Q at εr <= 12, and are now replaced
   by fits to 48 ringdowns, kept as a reported comparison. The DRA family has
   no curve-fit physics left in it.

9. ~~**A finite-gap feed for the MoM.**~~ - answered with `otahub/num/bor.py`
   instead: a straight conductor solved as a solid tube with a gap of physical
   length. The answer is that past a = 0.005 wavelengths the resonance belongs
   to the feed as much as the wire - gaps of a/2 to 4a spread the resonant length
   +-0.3% at 0.005, +-0.9% at 0.01, +-2-4% at 0.015-0.025 wavelengths, and at
   0.025 a gap of a/2 leaves no resonance at all. `resonant_dipole`'s length law
   lies inside that spread throughout; the slot family is now checked to w/L of
   about 0.09 and says plainly that past it L and R need a feed model. The wire
   MoM itself still has only the delta gap. What remains unchecked is the strip
   to tube equivalence a = w/4 for WIDE strips, which needs a strip solver.
10. ~~**`conical_monopole` and `biconical` bandwidths**~~ - done with
   `otahub/num/bor.py`. Matched to its own Zc a bicone holds a decade from 10 to
   65 degrees; a conical monopole in 50 ohm only from 30 to 55. The biconical
   exporter was found still halving theta_h - the spec's old full-angle bug,
   surviving in the builder - and is fixed.
11. ~~**The diagonal and dual-mode horns' fixed figures**~~ - done with the
   two-component aperture transform in `otahub/num/horn_pattern.py`. The Potter
   horn's step is now solved too (§3, fifth reference). **Left:** design the
   phasing length with the stepped-cone cascade rather than first-order - the
   flare moves the aperture phase ~20 degrees and steepens its frequency slope
   by a third - and tabulate it over horn length, input guide and share; that is
   what would take the horn out of low confidence.
12. ~~**Patch directivity, family-wide.**~~ - done, with a full-wave arbiter
   (`otahub/num/patch_sdm.py`, §3) that settled which model is right: the cavity
   current radiated through the slab (`patch_q`) holds to 0.3-0.9%; the two-slot
   formula (7.5-13% low) and the free-space edge models (thin-substrate only)
   do not. Six specs now carry thin-substrate directivity times a substrate
   factor. The annular ring followed: `patch_q.annulus` is its TM11 mode (it
   meets the two edge walls' ring currents to within 0.7% on thin board), so the
   ring now carries a substrate factor (up to 1.20 on thick board) and a radiation
   Q it never had. PIFA's length is FDTD-calibrated on two air plates (L + h =
   0.2425 wavelengths: the lambda/4 rule resonates 3% low); its bandwidth
   formula is 14-43% too wide against the FDTD's Q, and its directivity stays a
   placeholder. **Left:** a far-field transform in the half-space FDTD, which
   would give the PIFA a directivity - on a finite ground plane, which is what
   sets a real PIFA's pattern.
13. ~~**Patch resonance and Q, full-wave.**~~ - done for the rectangle, with an
   edge-exact basis and an FDTD check (§3); `rectangular_patch` and the inset
   patch build the full-wave length, and the CP square takes its side and Q0
   full-wave (its cut-to-split relation stays the cavity model's - the MoM does
   not model the cut). The shorted and stacked patches take the rectangle's
   length correction, checked by FDTD spot checks (§3). **Left, as future work:**
   a shorted-patch survey with the half-space FDTD (hours per board) to replace
   the borrowed correction; a two-layer spectral solver with a feed model to
   design the stack's gap and size ratio; and a full-wave solver of their own
   shape for the circular and triangular patches.
---

## 7. Conventions worth preserving

- `otahub.core` stays GUI-free and importable on bare numpy/scipy.
- Everything is SI internally; conversion happens only at the UI edge.
- Unknown symbols are errors, never zero — a typo must not silently synthesise
  a shorter antenna.
- Material properties and declared assumptions may default; **requirements
  never do**, because inventing a design target changes the antenna.
- Every archetype needs at least one citable `known_case`, or nothing about it
  is verifiable.
- Expectations of zero must use `tol_abs`, **in their own known case** — the
  tolerance applies to every expectation in the case, and will otherwise drag a
  good non-zero value into a 1e-9 comparison.
- Anything indicative rather than derived must say so in its `notes`.
- Compute `known_cases` expectations with a short independent script, never by
  hand. Across session 5 the formulas were nearly always right and my hand
  arithmetic was wrong about a dozen times; the one family whose expectations
  were computed programmatically landed with zero failures on the first run.
- Prefer expressions that stay correct in their limits over expressions that
  need a guard. Clamping an `acos` argument at 0 rather than −1 made the
  travelling-wave directivity return 1.5 for a short wire — the exact
  uniform-current dipole value — instead of zero.
