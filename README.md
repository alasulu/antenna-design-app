# OTA Hub Antenna Toolkit

Antenna, array and waveguide synthesis with export to full-wave simulators —
an open reimplementation of the Antenna Magus workflow: state your electrical
requirements, get a parameterised geometry and its predicted performance.

```bash
python OTA_Hub_AntennaToolkit.py gui                        # graphical interface
python OTA_Hub_AntennaToolkit.py list
python OTA_Hub_AntennaToolkit.py synth half_wave_dipole --f0 2.4GHz
python OTA_Hub_AntennaToolkit.py check                      # every citable known case

python OTA_Hub_AntennaToolkit.py guide WR-90 --f0 10GHz     # waveguides
python OTA_Hub_AntennaToolkit.py array -n 16 --taper chebyshev --sll -30
python OTA_Hub_AntennaToolkit.py planar --nx 16 --ny 16 --taper chebyshev --sll 30 --scan 45
python OTA_Hub_AntennaToolkit.py planar --nx 12 --ny 12 --lattice triangular --d 0.62
python OTA_Hub_AntennaToolkit.py line microstrip --z0 50 --h 1.6mm --eps-r 4.4
python OTA_Hub_AntennaToolkit.py match --r 200 --x -100 --z0 100 --f0 500MHz
python OTA_Hub_AntennaToolkit.py touchstone measured.s1p --compare half_wave_dipole

python OTA_Hub_AntennaToolkit.py export rectangular_patch_inset \
    --f0 2.4GHz --set eps_r=4.4 --set h=0.0016 --format cst -o patch.bas
```

## The catalogue

72 archetypes across 10 families, every one carrying citable `known_cases` that
run as tests:

| Family | n | Contents |
|---|---|---|
| `wire` | 11 | dipoles (including arbitrary-length and over-ground, solved exactly in Si/Ci), monopoles, loaded and top-loaded verticals, folded dipole, biconical, turnstile |
| `patch` | 9 | rectangular, inset-fed, circular, triangular, annular ring, quarter-wave shorted, truncated-corner CP, PIFA, stacked |
| `loop` | 8 | small circular and square, multi-turn, ferrite rod, resonant loop, quad, halo, Alford |
| `horn` | 8 | pyramidal, E- and H-plane sectoral, conical, corrugated, diagonal, dual-mode (Potter), open-ended guide |
| `travelling_wave` | 8 | Yagi-Uda, LPDA, axial and normal-mode helix, terminated long wire, V, rhombic, leaky-wave line source |
| `uwb` | 8 | bowtie, planar monopoles, Archimedean and equiangular spirals, Vivaldi, discone, conical monopole |
| `reflector` | 7 | prime-focus and offset parabolic, Cassegrain, Gregorian, parabolic cylinder, 90- and 60-degree corner |
| `slot` | 6 | half-wave, cavity-backed, folded, waveguide longitudinal, resonant and travelling-wave slot arrays |
| `lens` | 4 | plano-hyperbolic dielectric, Luneburg, Fresnel zone plate, metal-plate |
| `dielectric` | 3 | cylindrical, hemispherical and rectangular resonator antennas |

## The central design decision: archetypes are data, not code

An antenna archetype is a JSON document under [`specs/`](specs/) — its
parameters, synthesis formulas, analysis formulas, validity ranges, and
citations. Formulas are *strings*, evaluated by a whitelisted AST evaluator in
[`otahub/core/expr.py`](otahub/core/expr.py).

This buys three things:

1. **Adding an antenna never means writing Python.** A spec file is the whole
   contribution.
2. **Nothing in a spec can execute.** Imports, attribute escapes,
   comprehensions and lambdas are rejected — enforced by test, not convention.
3. **Every cited number is a test.** Each archetype's `known_cases` entries
   become pytest cases automatically, so a formula that drifts away from its
   reference fails the build instead of quietly producing a plausible, wrong
   antenna.

That third point is the one that matters for an engineering tool. Antenna
design formulas are full of coefficients that look right and are not; the
only defence is making the textbook value an executable assertion.

## Honesty guarantees

The engine is built to refuse rather than guess:

- **Unknown symbols are errors, never zero.** A typo in a formula fails loudly;
  it does not silently synthesise a shorter antenna.
- **Synthesis is partial, and says so.** Rules it cannot resolve become
  warnings listing exactly which requirements would unlock them
  (`DesignResult.unresolved`). `strict=True` demands a complete geometry.
- **Material properties default; requirements never do.** Copper and
  `mu_r = 1` are defensible assumptions. A design target is not, so the engine
  will not invent one.
- **Out-of-band use warns.** Every archetype carries an honest validity band;
  synthesising outside it produces a warning rather than a confident number.
- **Low-confidence specs announce themselves** in `list`, `show` and in every
  design they produce.

## Reading measurements back

`export` sends a model out to a solver. `touchstone` is the way back: read the
S-parameters a solver or a VNA produced and put them beside what the archetype
predicted.

```bash
python OTA_Hub_AntennaToolkit.py touchstone patch.s1p --compare rectangular_patch_inset \
    --set eps_r=4.4 --set h=1.6e-3
```

It finds the resonances, reports impedance and VSWR there, and states the
prediction and the measurement in the same terms.

Two details in the Touchstone format are silent failures if you get them wrong,
and both are tested here. **Two-port files store their matrix column-major**
(`freq S11 S21 S12 S22`) while three ports and up are row-major — reading a
2-port row-major transposes it, swapping forward gain with reverse isolation.
And **a frequency point may wrap across any number of lines**, so the reader
works on a flat stream of values rather than line by line. A nine-value first
line is genuinely ambiguous between a 2-port and a 4-port; the reader resolves
it from the total value count where it can and says so plainly where it cannot,
rather than guessing.

## Verification

`otahub/core/pattern.py` is deliberately independent of the spec files —
first-principles pattern integration, so it can serve as an external yardstick
rather than restating what the specs claim. It reproduces:

| Quantity | Computed | Reference |
|---|---|---|
| Short dipole D | 1.500000 | 1.5 exact |
| Half-wave dipole D | 1.640922 | 1.6409 (Balanis) |
| Half-wave HPBW | 78.08° | 78° |
| cos^q feed D | 6.00011 (q=1) | 2(2q+1) = 6 |
| Uniform line source SLL | −13.263 dB | −13.26 dB |
| Uniform line source HPBW | 5.078° | 50.8 λ/L = 5.08° |

```bash
python -m pytest tests/ -q
python OTA_Hub_AntennaToolkit.py doctor   # structural faults in the specs
```

## Layout

| Path | Contents |
|---|---|
| `otahub/core/` | Engine: spec model, evaluator, synthesis solver, pattern maths. **GUI-free by invariant** — it must stay importable on a bare numpy/scipy install. |
| `otahub/arrays/` | Tapers (uniform, binomial, Dolph-Chebyshev, Taylor), array factor, grating-lobe limits; planar rectangular and triangular lattices with steering, scan loss and beam cuts |
| `otahub/waveguides/` | Rectangular and circular guides, coax, microstrip, stripline, CPW |
| `otahub/utils/` | Network parameters (S/Z/Y/ABCD), matching networks, Touchstone import/export |
| `otahub/export/` | CST Studio VBA and Ansys HFSS script generation |
| `otahub/cli/` | Command line front end |
| `otahub/gui/` | PySide6 interface: catalogue, linear arrays, planar arrays, waveguides |
| `specs/` | Archetype definitions — see [`specs/SPEC_FORMAT.md`](specs/SPEC_FORMAT.md) |
| `tests/` | Hand-written engine tests plus cases generated from every spec |

## Requirements

Python 3.10+, numpy, scipy, matplotlib. `pytest` for the suite; `PySide6` for the GUI.

## Export to simulators

`export` emits a parameterised CST VBA macro or HFSS IronPython script. Lengths
become named variables in millimetres, frequencies in GHz, and everything else
keeps its own units — driven by each spec's declared unit rather than guessed
from magnitude.

Geometry is built for **43 of the 72 archetypes** — those whose construction is
unambiguous from their primary dimensions:

| Group | Archetypes |
|---|---|
| Wire | every one: dipoles (half-wave, resonant, short, arbitrary length), monopoles (quarter-wave, top-loaded, inductively loaded), folded dipole, dipole over ground, turnstile, biconical |
| Patch | rectangular, inset-fed, circular, quarter-wave shorted, annular ring, PIFA, stacked |
| Dielectric resonator | rectangular, cylindrical, hemispherical |
| Wideband | biconical, conical monopole, discone, planar monopoles (disc and rectangular) |
| Loop | small circular and square, one-wavelength circular, multi-turn, quad, Alford, halo |
| Slot | half-wave, folded, cavity-backed, waveguide longitudinal, resonant and travelling-wave arrays |
| Travelling wave | terminated long wire, leaky-wave line source |
| Lens | Fresnel zone plate, Luneburg |
| Aperture | open-ended waveguide |

**Everything else exports its parameters and states plainly that no solid
geometry was generated.** A half-built model that looks
finished is worse than none, so the exporter refuses to guess — and where it
does make a choice the spec cannot supply (a feed gap, an inset notch width, a
finite ground plane standing in for an infinite one) it says so in the file
header.

## Verification

```bash
python -m pytest tests/ -q                # 1747 tests
python OTA_Hub_AntennaToolkit.py check    # every archetype against its citations
python OTA_Hub_AntennaToolkit.py doctor   # structural faults in the specs
```

Worked examples reproduced exactly: Balanis 14.1 (patch W/L), 14.2 (inset
resistance), 14.4 (circular patch radius); Pozar 5.1 (L-section) and 5.2
(single-stub tuner); Kraus helix; Viezbicke NBS Yagi gains; WR-90 datasheet
cutoff, attenuation and power. Dolph-Chebyshev sidelobes match their design
level to 0.000 dB.

Several archetypes carry results derived here rather than quoted, each checked
against an independent numerical model before being written into a spec:

| Result | How it was obtained | Check |
|---|---|---|
| Travelling-wave wire radiation resistance | Integrated the pattern analytically | Reduces to 80π²(l/λ)² as l → 0 |
| Dipole over ground | Image theory with exact mutual impedance | Hemisphere integration, to 5 digits |
| Turnstile directivity | Summed-power spherical integration | 1.64092 on axis, exactly a single dipole |
| V and rhombic directivity | Four-leg travelling-wave model, fitted | 1.3% max fit error over 1.5–12 λ |
| Corner reflector image sets | Boundary condition on the plates | 1e-15 residual tangential E |
| Diagonal horn aperture efficiency | Aperture integration | 8/π² = 0.8106 analytic vs 0.8110 numeric |
| Annular ring resonance | Bisection on the Bessel cross-product | Fit error 0.20%, vs 2.71% for the textbook narrow-ring rule |
| Hemispherical DRA resonance and Q | Mie magnetic-dipole resonance | Q ∝ εr^1.32, independently reproducing the published εr^1.3 |

### Planar arrays

`otahub/arrays/planar.py` covers rectangular and equilateral-triangular
lattices, separable tapers, beam steering, and directivity for an arbitrary
element layout.

Directivity is exact rather than integrated: the average of `exp(j k·d)` over
the sphere is `sin(kd)/(kd)`, so the radiated power is a double sum over
`sinc(2|rₘ − rₙ|/λ)`. For a single row it reproduces the linear module to
machine precision, and it agrees with brute-force spherical integration for
every planar layout tested.

The triangular lattice's advantage falls out of the reciprocal lattice rather
than being asserted: grating lobes sit at reciprocal-lattice points, the
shortest such vector is `4π/(√3 s)` against a rectangular lattice's `2π/d`, so
the triangular unit cell is `2/√3` larger in area — **13.40% fewer elements**
for the same grating-lobe-free scan volume. A test recomputes that from the
primitive vectors, and another confirms by brute force that a grating lobe
appears just past the limit and not before.

### Horns solved rather than assumed

The sectoral and pyramidal horns carried a pinned aperture efficiency — 0.65,
0.64, 0.51 — and synthesised the aperture so the horn was always at its optimum
flare by construction. Self-consistent, but it could not analyse a horn someone
already had, and worse, **its gain rose without limit as the flare grew**: with
efficiency fixed, over-flaring always looked better.

Balanis gives the directivity exactly, in Fresnel integrals. Rearranged, the
efficiency is closed-form for *any* flare:

```
η_E = (8/π²)·[C(q)² + S(q)²]/q²          q = b₁/√(2λρ)
η_H = (λρ/a₁²)·{[C(u) − C(v)]² + [S(u) − S(v)]²}
η_P = (π²/8)·η_E·η_H
```

The 8/π² is the TE₁₀ cosine taper across the unflared plane; in the pyramidal
form it appears in both sectoral efficiencies and has to be removed once. At the
optimum these give 0.64870, 0.64276 and 0.51440 against the pinned 0.65, 0.64
and 0.51 — and away from it they give the right answer instead of a badly wrong
one. Gain now peaks at the optimum flare, which is the property that makes it
optimum.

### Cross-consistency between archetypes

`tests/test_cross_consistency.py` checks archetypes against **each other**, at
the places where the same physics reaches two specs by different routes: a
quarter-wave monopole is exactly half a half-wave dipole, a folded dipole four
times one, a turnstile's on-axis directivity equal to a single dipole's, a slot
times its complementary dipole equal to η₀²/4, a conical monopole exactly half a
biconical, a slot array's guide conductance identical to the single-slot spec's,
the Potter horn's mode cutoffs equal to the circular-guide Bessel zeros, and the
three dielectric-resonator shapes — solved by exact Mie theory, by published
curve fits and by magnetic-wall algebra respectively — agreeing on volume within
a factor of 1.41.

It exists because neither a citation nor a unit check catches the failure mode
that matters most here. `corner_reflector_90` shipped for four sessions with its
null where its optimum is, passing its own cited cases throughout;
`rectangular_patch_inset` agreed with Balanis to 1.7% while both of its
conductance terms were wrong in cancelling directions. A second independent
route to the same number catches both. Injecting a 3% drift into one shared
archetype fires six of these tests across three families.

### The dimensional audit

`tests/test_scale_invariance.py` checks every archetype against physics rather
than against a reference. Maxwell's equations have no preferred length: multiply
every frequency by S, divide every length by S, scale conductivity by S, and the
antenna's electrical behaviour must be *identical* — same directivity, same
beamwidths, same impedances, same efficiency — while lengths shrink by S, areas
by S², and so on.

793 quantities across all 72 archetypes are checked this way. It needs no
textbook, which is what makes it worth running over formulas whose known cases
were written by the same hand. It found a patch radius carried in centimetres
while labelled dimensionless, and a `m2` unit typo, on specs that had been
passing their own cited cases since the first session.

Cross-checks between independently written specs are the useful kind, and
several hold to machine precision: the slot array's guide conductance equals
`waveguide_longitudinal_slot`'s; `conical_monopole`'s impedance is exactly half
`biconical`'s at the same flare angle; the Potter horn's mode cutoff diameters
agree with the circular waveguide module's Bessel zeros.
