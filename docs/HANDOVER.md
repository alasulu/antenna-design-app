# OTA Hub Antenna Toolkit — handover

An open reimplementation of the Antenna Magus workflow: state electrical
requirements, get a parameterised geometry, its predicted performance, and an
export to CST Studio or Ansys HFSS.

Built across five sessions on 2026-09-13/14. **72 archetypes, 10 families,
6,031 lines of Python, 14,835 lines of spec data, 1262 tests, 426/426 citable
known cases passing.**

---

## 1. What exists

| Area | Module | Contents |
|---|---|---|
| Engine | `otahub/core/` | Spec model, whitelisted AST evaluator, partial synthesis solver, registry, first-principles pattern maths |
| Catalogue | `specs/*.json` | 72 archetypes across wire (11), patch (9), loop (8), horn (8), travelling-wave (8), UWB (8), reflector (7), slot (6), lens (4), dielectric (3) |
| Arrays | `otahub/arrays/` | Uniform, binomial, Dolph-Chebyshev, Taylor n-bar, raised-cosine tapers; linear array factor, steering, grating-lobe limits; planar rectangular and triangular lattices with exact directivity, scan loss and beam-following cuts |
| Waveguides | `otahub/waveguides/` | Rectangular and circular guides, exact WR-series table, coax, microstrip, stripline, CPW |
| Utilities | `otahub/utils/` | S/Z/Y/ABCD conversion and cascading; L-section, quarter-wave and single-stub matching |
| Export | `otahub/export/` | Neutral geometry IR rendered to CST VBA and HFSS IronPython |
| Interfaces | `otahub/cli/`, `otahub/gui/` | 12 CLI subcommands; PySide6 GUI with catalogue, linear-array, planar-array and waveguide tabs |

```bash
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
| `corner_reflector_90`, `corner_reflector_60` | Image array factors | Summed field leaves ~1e-15 tangential E on the plates |
| `diagonal_horn` | Aperture efficiency 8/π² = 0.8106 | Aperture integration on a 2001² grid: 0.8110 |
| `annular_ring_patch` | Cubic correction to the narrow-ring rule | Bisection on the exact Bessel cross-product; 0.20% error against 2.71% uncorrected |
| `hemispherical_dra` | k₀a = 2.900 εr^−0.484, Q = 0.380 εr^1.321 | First peak of the Mie magnetic-dipole coefficient; Q exponent independently reproduces the published εr^1.3 |

**These fits are only as good as the model behind them.** Each is a
closed-form or ray-optics idealisation, not a full-wave result, and the
validity block on each archetype says where it stops.

---

## 4. What is NOT verified — read this before trusting a number

### Nine archetypes are marked low confidence

`cassegrain`, `conical_horn_dual_mode`, `ferrite_rod_loop`, `halo_loop`,
`pifa`, `planar_monopole_circular`, `stacked_patch`, `vivaldi_tsa`,
`waveguide_longitudinal_slot`.

Two of those are new in session 5 and both are honest about why:
`stacked_patch`'s bandwidth multiplier is an expectation drawn from published
designs rather than a computed result, and `conical_horn_dual_mode`'s
efficiency, cross-polar level and beamwidth constant are placed by analogy
with its neighbouring horns. In the Potter horn only the two mode cutoff
diameters are exact — they are Bessel zeros.

They announce themselves in `list`, `show`, the GUI, and in every design they
produce. Treat their numbers as indicative and verify in a full-wave solver.

### Known open discrepancies

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
- **Axial-mode helix gain.** The Kraus formula is known to overestimate for
  large N, and the published corrections disagree by 1–2 dB. Both the classic
  and a corrected variant are exposed; neither should be trusted to better
  than a couple of dB.
- **Stevenson's g1** for the waveguide shunt slot should be cross-checked
  against the source before committing an array to fabrication.
- **LPDA directivity** collapses a two-dimensional Carrel chart onto the
  optimum-σ line; expect ~1 dB error.
- **Taylor taper** realises its design sidelobe level to about 1 dB for small
  arrays (−28.9 dB measured for a 20-element −30 dB design).
- **`rectangular_dra`'s radiation Q is borrowed**, not derived: it is the
  hemispherical DRA's exact result reused. Shape matters less than permittivity
  here — the cylindrical archetype's independent fit sits within 6% of the
  hemisphere's at εr = 10 — but aspect ratio moves a rectangular DRA's Q by
  considerably more than that.
- **`rectangular_dra`'s resonance uses an all-magnetic-wall model**, whose
  algebra is exact and whose physical assumption is not: it predicts f₀ high by
  roughly 10–20%. The spec says so in a dedicated note.
- **`discone` and `conical_monopole` rest on engineering conventions** — the
  quarter-wavelength slant and the decade bandwidth figure — not on derivations.
- **Fresnel zone plate efficiencies** (1/π², 4/π², 8/π²) are the standard
  first-order grating results. They ignore a real phase plate's finite thickness
  and its shadowing at angle.

### Metrics labelled "indicative"

Several archetypes carry bandwidth, front-to-back or directivity figures that
are engineering rules of thumb rather than derived results. Every one says so
in its `notes`. They exist so a design sheet is complete, not because they are
predictions.

### Exporters build geometry for 18 of 72 archetypes

Dipoles (3), monopole, folded dipole, dipole over ground, turnstile, four
patch variants, three dielectric resonators, biconical, conical monopole,
discone, and open-ended waveguide. **The other 54 export parameters only and
say so.**
Where a builder must choose something the spec cannot supply — feed gap, inset
notch width, a finite ground plane standing in for an infinite one — the
choice is stated in the exported file's header.

Neither exporter has been run against a real CST or HFSS installation. The
scripts are structurally validated (balanced blocks, parseable Python, correct
units) but **not execution-tested**. That is the single largest untested
surface in the project.

---

## 5. Engine bugs found during the build

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

---

## 6. Recommended next steps

1. **Execution-test the exporters** against real CST and HFSS installations.
   This is the biggest gap, and session 5 did not touch it — the catalogue grew
   by 32 archetypes while the exporter still builds geometry for 7. Start with
   `half_wave_dipole` (simplest geometry, strongest analytical reference:
   73.08 + j42.52 Ω at λ/2).
2. ~~**Close the inset-patch discrepancy**~~ — done. The published figure was
   right; the "direct integration" that disagreed with it was wrong. See §4.
3. **Add geometry builders** for horns, which need a loft or truncated-pyramid
   primitive neither backend abstraction has yet, and for the loop family,
   which needs a torus. Yagi-Uda remains blocked for a different reason: the
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
6. **Touchstone import** so measured or simulated S-parameters can be read
   back and compared against predictions.
7. **Re-audit the session 1–4 specs** — partly done. Two systematic audits now
   run over the whole catalogue: the dimensional one (793 quantities against
   scale invariance) and the cross-consistency one (23 relationships between
   archetypes that must agree). Both are clean. What neither covers is an
   archetype with no sibling and no dimensional quirk — those still need
   reading against their source.
   Every archetype passes the cases it declares; that is not the same as being
   right. `corner_reflector_90` had a null where its optimum is and said so
   confidently for four sessions. The archetypes carrying a single known case
   are the place to start.
8. **Replace `rectangular_dra`'s borrowed Q** with a proper solve, and its
   magnetic-wall resonance with the dielectric-waveguide transcendental. The
   machinery used for the hemispherical DRA transfers directly.

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
