# Archetype spec contract

Every family agent writes `specs/<family>.json` matching this shape EXACTLY.
The implementation is generated against it, so unit discipline is not optional.

```json
{
  "family": "wire",
  "archetypes": [
    {
      "key": "half_wave_dipole",
      "name": "Half-wave dipole",
      "summary": "one sentence",
      "freq_range_hz": [3e6, 3e11],
      "parameters": [
        {"symbol": "L", "name": "total_length", "unit": "m", "role": "geometry",
         "description": "tip-to-tip length", "typical": "0.47*lambda"}
      ],
      "synthesis": [
        {"output": "L", "expr": "0.4788 * c / f0", "units_out": "m",
         "notes": "includes end-effect shortening for thin wire",
         "depends_on": ["f0"]}
      ],
      "analysis": [
        {"metric": "directivity_dbi", "expr": "2.15", "units_out": "dBi",
         "notes": "exact value 1.643 linear"},
        {"metric": "input_impedance_ohm", "expr": "73.1 + 42.5j", "units_out": "ohm",
         "notes": "at exactly L=lambda/2; resonant length is ~0.47 lambda where X->0"}
      ],
      "validity": ["thin wire, a << lambda", "free space, no ground plane"],
      "known_cases": [
        {"given": {"f0": 300e6}, "expect": {"L_m": 0.4788}, "tol_pct": 2,
         "source": "Balanis 4th ed. Ch.4"}
      ],
      "references": ["Balanis, Antenna Theory 4th ed., Ch. 4.6"]
    }
  ]
}
```

## Hard rules
- **SI internally.** All `expr` must evaluate in SI (m, Hz, ohm, W). Convert at the UI edge only.
- `expr` is a Python expression over the parameter symbols plus: `c` (2.99792458e8), `pi`, `eps0`, `mu0`, `eta0` (376.730313412), and numpy as `np`. No statements, no imports.
- `depends_on` lists the symbols the expression reads. It must be accurate — the engine topologically sorts on it.
- Every archetype needs >=1 `known_cases` entry with a real numeric expectation. These become pytest cases. An archetype with no verifiable known case is worthless to us.
- Parameter `role` is one of: `requirement` (no default, ever), `geometry`, `material` or `assumption` (both auto-default from a numeric `typical`), `derived`.
- `freq_range_hz` is the honest validity band, not the band someone could force it into.
- Prefer closed-form engineering formulas with stated accuracy over hand-waving. If a quantity genuinely needs a numerical solve, say so in `notes` and give the defining equation.
- `tol_abs` applies to EVERY expectation in its case, not just the zero one. Put a
  zero expectation in its own `known_cases` entry, or it will drag a perfectly good
  non-zero expectation into a 1e-9 comparison and fail it.
- Prefer expressions that stay correct in their limits. Clamping an `acos` argument at
  0 rather than -1 made the travelling-wave directivity return 1.5 for a short wire —
  the exact uniform-current dipole value — instead of zero. A formula that degrades
  gracefully outside its band is worth more than one that needs a guard.
- Compute the numbers for `known_cases` with a short independent script, not by hand.
  Hand arithmetic was the single largest source of false failures while this catalogue
  was being written; the formulas were nearly always right and the expectations were not.
- If you are unsure of a coefficient, mark the archetype `"confidence": "low"` and say why in `notes`. Do NOT invent precision.
