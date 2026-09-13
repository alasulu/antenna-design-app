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
python OTA_Hub_AntennaToolkit.py line microstrip --z0 50 --h 1.6mm --eps-r 4.4
python OTA_Hub_AntennaToolkit.py match --r 200 --x -100 --z0 100 --f0 500MHz

python OTA_Hub_AntennaToolkit.py export rectangular_patch_inset \
    --f0 2.4GHz --set eps_r=4.4 --set h=0.0016 --format cst -o patch.bas
```

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
| `otahub/arrays/` | Tapers (uniform, binomial, Dolph-Chebyshev, Taylor), array factor, grating-lobe limits |
| `otahub/waveguides/` | Rectangular and circular guides, coax, microstrip, stripline, CPW |
| `otahub/utils/` | Network parameters (S/Z/Y/ABCD) and matching networks |
| `otahub/export/` | CST Studio VBA and Ansys HFSS script generation |
| `otahub/cli/` | Command line front end |
| `otahub/gui/` | PySide6 interface: catalogue, arrays, waveguides |
| `specs/` | Archetype definitions — see [`specs/SPEC_FORMAT.md`](specs/SPEC_FORMAT.md) |
| `tests/` | Hand-written engine tests plus cases generated from every spec |

## Requirements

Python 3.10+, numpy, scipy, matplotlib. `pytest` for the suite; `PySide6` for
the GUI once it lands.

## Export to simulators

`export` emits a parameterised CST VBA macro or HFSS IronPython script. Lengths
become named variables in millimetres, frequencies in GHz, and everything else
keeps its own units — driven by each spec's declared unit rather than guessed
from magnitude.

Geometry is built only for archetypes whose construction is unambiguous from
their primary dimensions: dipoles, monopoles, rectangular and circular patches,
and open-ended waveguide. **Everything else exports its parameters and states
plainly that no solid geometry was generated.** A half-built model that looks
finished is worse than none, so the exporter refuses to guess — and where it
does make a choice the spec cannot supply (a feed gap, an inset notch width, a
finite ground plane standing in for an infinite one) it says so in the file
header.

## Verification

```bash
python -m pytest tests/ -q                # 400+ tests
python OTA_Hub_AntennaToolkit.py check    # every archetype against its citations
python OTA_Hub_AntennaToolkit.py doctor   # structural faults in the specs
```

Worked examples reproduced exactly: Balanis 14.1 (patch W/L), 14.2 (inset
resistance), 14.4 (circular patch radius); Pozar 5.1 (L-section) and 5.2
(single-stub tuner); Kraus helix; Viezbicke NBS Yagi gains; WR-90 datasheet
cutoff, attenuation and power. Dolph-Chebyshev sidelobes match their design
level to 0.000 dB.
