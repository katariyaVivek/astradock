# M21 Mission Simulation & Monte Carlo Framework — Validation Report

## Scope

`cpp/mission/mission_framework.hpp` (config, hash, identity, timeline,
statistics, regression missions), `tests/cpp/test_mission.cpp` (6 cases / 55
assertions), demo (`tools/mission_demo.cpp`: 4 missions × 20 config-driven ECI
runs), oracle (`python/audit/independent_mission_reference.py`), telemetry
(`data/m21_*.csv` + scorecard), figures (`artifacts/figures/m21_*`), lesson
(`docs/lessons/021_mission_framework.md`).

## Analytical checks (all pass)

- Hash deterministic; sensitive to all 6 field groups; invalid configs throw.
- Seed streams distinct + stable; version pinned "0.1.0".
- Timeline boundaries (300/2400/2940 s); empty/singleton behavior.
- Uniform 1..100 grid: mean 50.5, sample std, p5 = 5.95, p50 = 50.5,
  p95 = 95.05, success 0.95, mean Δv 5.05 — all exact.
- Outcome counts 1/1/1/1; pass/fail both directions; NaN rejected.

## Independent verification

Pure-Python oracle (seeds, percentiles, Welford, policy); 4 spot checks pass.

## Scenario validation

| Mission | Result | p50 | p95 | Verdict |
|---|---|---|---|---|
| nominal | 20/20 | 260 m | 261 m | PASS |
| dropout | 20/20 | 260 m | 260 m | PASS |
| hot | 20/20 | 260 m | 261 m | PASS |
| docking-grade | 20/20 | 60 m | 60 m | PASS |

Every CSV carries scenario_id + base seed + config hash + version in a `#`
header comment. Distinct configs hash distinctly (verified in output).

## Regression

New `cpp/mission/` + tests/demo only. Full suite 290/290. Ruff clean.

## Limitations

Structs, not files (no cross-language exchange); percentiles stored not
streamed (fine to N ~ 1e5); 20 demo runs per mission (framework supports any
N via mc_runs); no parallel execution (determinism first).
