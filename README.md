# OTA Hub Antenna Toolkit

Antenna, array and waveguide synthesis with export to full-wave simulators —
an open reimplementation of the Antenna Magus workflow: state your electrical
requirements, get a parameterised geometry and its predicted performance.

```bash
python OTA_Hub_AntennaToolkit.py list
python OTA_Hub_AntennaToolkit.py synth half_wave_dipole --f0 2.4GHz
python OTA_Hub_AntennaToolkit.py check      # run every citable known case
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
| `otahub/arrays/`, `otahub/waveguides/` | Array and waveguide synthesis |
| `otahub/export/` | CST Studio VBA and Ansys HFSS script generation |
| `otahub/cli/` | Command line front end |
| `specs/` | Archetype definitions — see [`specs/SPEC_FORMAT.md`](specs/SPEC_FORMAT.md) |
| `tests/` | Hand-written engine tests plus cases generated from every spec |

## Requirements

Python 3.10+, numpy, scipy, matplotlib. `pytest` for the suite; `PySide6` for
the GUI once it lands.
