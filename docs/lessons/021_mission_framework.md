# Lesson 021 — Mission Simulation & Monte Carlo Framework

## 1. Physical Problem & Motivation

M13–M20 demos hard-code their scenarios: orbit, noise, faults, waypoints,
gains, seeds scattered across `main()` bodies. Nothing is comparable or
auditable across runs — change one gain and no record says which runs used
what. M21 centralizes the scenario as DATA (`MissionConfig`: every knob, one
struct), stamps every run with identity (scenario ID, seed, config hash,
version), drives execution through phase timelines, and scores seeded Monte
Carlo with percentiles and failure classification.

## 2. Design Decisions

- **Structs, not YAML**: a file format buys a parser dependency for zero
  physics benefit. The struct IS the schema, printable to CSV headers.
- **FNV-1a config hash** over ordered fields: same config → same hash on any
  platform; any field change flips it. Catches "same seed, different config",
  the silent Monte Carlo killer.
- **Per-run seeds** `base * stride + index`: independent streams, one base
  seed reproduces everything.
- **Percentiles by sort** (N ≤ thousands — storing doubles is honest), mean /
  variance by Welford (online, exact).
- **Outcomes classified** (converged / diverged / aborted / timeout) by the
  scoring harness; the summary counts each exactly once.

## 3. Components

`MissionConfig` + `validate()`; `config_hash()`; `RunIdentity`;
`derive_seed()`; `active_phase()` + `default_rendezvous_timeline()`;
`ScoredRun` + `summarize_runs()` (mean/std/p5/p50/p95/success rate);
`RegressionMission` + `passes()` + `default_regression_missions()` (4 missions).

## 4. Verification

6 cases / 55 assertions: hash determinism + per-group sensitivity, validation
rejection, seed streams, timeline boundaries, uniform-grid statistics
(mean 50.5, std, p5/p50/p95 exact), outcome counts, pass/fail directions.
Oracle (seeds, percentiles, Welford, policy). Demo: 4 missions × 20 runs, all
PASS — nominal/dropout/hot converge 20/20 (p95 ≈ 260 m), docking-grade 20/20
(p95 ≈ 60 m). Scorecard + dispersion figures.

## What you should now understand

1. Why does "same seed" not imply "same experiment"?
2. What failure does the config hash prevent?
3. Why Welford for moments but sorting for percentiles?
4. Why classify outcomes instead of averaging errors?
5. When would a file-format config justify its parser?
6. What does p95 tell you that the mean hides?
