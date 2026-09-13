# OTA Hub Antenna Toolkit — handover

An open reimplementation of the Antenna Magus workflow: state electrical
requirements, get a parameterised geometry, its predicted performance, and an
export to CST Studio or Ansys HFSS.

Built across four sessions on 2026-09-13/14. **40 archetypes, 8 families,
5,977 lines of Python, 7,948 lines of spec data, 644 tests, 241/241 citable
known cases passing.**

---

## 1. What exists

| Area | Module | Contents |
|---|---|---|
| Engine | `otahub/core/` | Spec model, whitelisted AST evaluator, partial synthesis solver, registry, first-principles pattern maths |
| Catalogue | `specs/*.json` | 40 archetypes across wire, loop, patch, horn, reflector, travelling-wave, UWB, slot |
| Arrays | `otahub/arrays/` | Uniform, binomial, Dolph-Chebyshev, Taylor n-bar, raised-cosine tapers; array factor, steering, grating-lobe limits |
| Waveguides | `otahub/waveguides/` | Rectangular and circular guides, exact WR-series table, coax, microstrip, stripline, CPW |
| Utilities | `otahub/utils/` | S/Z/Y/ABCD conversion and cascading; L-section, quarter-wave and single-stub matching |
| Export | `otahub/export/` | Neutral geometry IR rendered to CST VBA and HFSS IronPython |
| Interfaces | `otahub/cli/`, `otahub/gui/` | 11 CLI subcommands; PySide6 GUI with catalogue, array and waveguide tabs |

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

## 4. What is NOT verified — read this before trusting a number

### Seven archetypes are marked low confidence

`ferrite_rod_loop`, `halo_loop`, `pifa`, `cassegrain`,
`waveguide_longitudinal_slot`, `planar_monopole_circular`, `vivaldi_tsa`.

They announce themselves in `list`, `show`, the GUI, and in every design they
produce. Treat their numbers as indicative and verify in a full-wave solver.

### Known open discrepancies

- **Inset patch mutual conductance.** Approximating G12 as `G1·J0(k0·L)`
  lands within 1.7% of Balanis Example 14.2, but numerically integrating
  eq. 14-18a directly gives 212.5 Ω against the published 228.35 Ω — a 6.9%
  disagreement that is recorded in the spec and unresolved. Inset depths carry
  real uncertainty.
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

### Metrics labelled "indicative"

Several archetypes carry bandwidth, front-to-back or directivity figures that
are engineering rules of thumb rather than derived results. Every one says so
in its `notes`. They exist so a design sheet is complete, not because they are
predictions.

### Exporters build geometry for only 7 archetypes

Dipoles (3), monopole, rectangular patch (2 variants), circular patch, and
open-ended waveguide. **Everything else exports parameters only and says so.**
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

---

## 6. Recommended next steps

1. **Execution-test the exporters** against real CST and HFSS installations.
   This is the biggest gap. Start with `half_wave_dipole` (simplest geometry,
   strongest analytical reference: 73.08 + j42.52 Ω at λ/2).
2. **Close the inset-patch discrepancy** — decide whether the 228.35 Ω
   published figure or the 212.5 Ω direct integration is right.
3. **Add geometry builders** for horns and Yagi-Uda; both have unambiguous
   constructions and are common export targets.
4. **Validate a low-confidence archetype** end to end and either promote it or
   record why it cannot be.
5. **Planar arrays** — the array module is linear-only; rectangular and
   triangular lattices are the obvious extension, and the pattern machinery
   already supports them.
6. **Touchstone import** so measured or simulated S-parameters can be read
   back and compared against predictions.

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
- Expectations of zero must use `tol_abs`.
- Anything indicative rather than derived must say so in its `notes`.
