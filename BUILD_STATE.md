# OTA Hub Antenna Toolkit — Build State

**Single source of truth for scheduled resume runs.** Read this first, update it last.

- Project root: `/Users/macbookair/Desktop/Antenna Design App`
- Python: `/Users/macbookair/miniconda3/bin/python` (numpy 2.2.6, scipy 1.16.3, matplotlib 3.10.8)
- Build order chosen by user: **breadth of archetypes first** (CLI only in S1; GUI in S3)
- Disk constraint: ~3.4 GB free. Do NOT install large deps without checking `df -h` first.

## Session plan

| # | When | Scope | Status |
|---|------|-------|--------|
| S1 | 2026-09-13 (interactive) | Core engine, archetype registry, ~40 antenna archetypes, CLI, tests | IN PROGRESS |
| S2 | 2026-09-13 15:12 CEST | Waveguides module + arrays (layouts, tapers, array factor) | PENDING |
| S3 | 2026-09-13 20:12 CEST | PySide6 GUI shell: catalog browser, param panel, plots | PENDING |
| S4 | 2026-09-14 01:12 CEST | Exporters (CST VBA, HFSS script), matching/network utils, polish | PENDING |
| — | after S4 | HALT. Await user typing `Continue`. | — |

## Invariants
- `otahub.core` must stay GUI-free and import-clean under plain numpy/scipy.
- Every archetype: `synthesize(requirements) -> params`, `analyze(params) -> metrics`, plus validity ranges.
- Every formula carries a `reference` string and a tested limiting case.

## Completed
(nothing yet)

## Next action
Implement core engine + archetype registry from `specs/`.
