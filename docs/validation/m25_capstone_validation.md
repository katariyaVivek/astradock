# M25 Full-System Verification & Capstone Mission — Validation Report

## Scope

Capstone demo (`tools/capstone_demo.cpp`: 9-phase mission + 100-run Monte
Carlo), verification tests (`tests/cpp/test_capstone.cpp`: isolation audit,
interface audit, frozen benchmarks), telemetry (`data/m25_*.csv`), figure
(`artifacts/figures/m25_capstone.png`), lesson (`docs/lessons/025_capstone.md`).

## Mission results

- Nominal (seed 42): converged; final 255 m; Δv 3.0 m/s; attitude errors
  sub-degree by docking phase; GNSS dropout 900–960 s ridden through
  (NIS gap visible in telemetry, monitor coasts, no trigger pileup).
- Monte Carlo: 100/100 converged; p50 255 m; p95 256 m; success 1.000;
  zero aborts/divergences/timeouts.

## Audits

- Isolation: declaration scan over 6 headers passes (1146 assertions incl.
  per-line checks); shared-type rationale documented.
- Interface: handedness, SI, JD, sign, ordering, double cover all pass.
- Benchmarks: 11 frozen rows pass against documented bounds.

## Regression

New demo + tests only. Full suite 311/311 (308 + 3 capstone). Ruff clean.

## Performance (M25G)

Full Release build + 311-test suite + capstone demo (101 mission runs) in one
session: no optimization work needed (header-only, seconds-scale builds,
minute-scale missions). Profiling unnecessary — documented as measured, not
assumed.

## Reproducibility (M25F)

Clean configure + build + ctest from this session's tree; Catch2 pinned
v3.7.1 with FetchContent fallback; Python extras declared in pyproject
(analysis, ml, dev). Data/artifacts git-ignored and regenerable via demos.

## Limitations (project-close statement)

Capstone holds at the rendezvous ball (docking contact stays in the M19 rig
with docking-grade nav); single orbit regime; static cooperative target;
educational atmosphere; advisory-only ML. Every limitation points at its
milestone report — nothing is hidden, nothing is overstated.
