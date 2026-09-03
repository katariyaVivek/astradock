# M17 Integrated Closed-Loop GNC — Validation Report

## Scope

`cpp/gnc/closed_loop.hpp` (estimate-based tick, timeline, metrics),
`tests/cpp/test_closed_loop.cpp` (7 cases), demo
(`tools/closed_loop_demo.cpp`: hold + 60° maneuver + 100-run Monte Carlo),
oracle (`python/audit/independent_closed_loop_reference.py`), telemetry
(`data/m17_*.csv`), figures (`artifacts/figures/m17_*`), lesson
(`docs/lessons/017_closed_loop_gnc.md`).

## Analytical / wiring checks (all pass)

- Estimate-based tick ≡ M16A law on identical inputs to 1e-12.
- Saturation clamp + flag; timeline selection incl. pre-switch hold.
- Static audit: no simulated-state-typed declaration in the seam headers.
- Perfect-estimate loop detumbles the master-spec tumble, settles < 90 s.
- MEKF-in-loop: attitude error shrinks below init and 0.01 rad; bias < 1 mrad/s.
- Metrics latching verified; mismatched/empty columns rejected.
- Full 6-DOF wheel-actuated hold: error < 0.57°, rate < 0.29°/s, orbit radius
  unchanged to 1e-6 (decoupling preserved through actuation).

## Independent verification

Pure-Python oracle (tick, timeline, settle); 4 spot cross-checks pass. No
production code imported.

## Scenario validation

Seed 42, 20 Hz loop, 0.2 Nm wheel envelope, 1 mrad sensors:
- M17A hold: settles 80.8 s (0.57° bound), final 0.038°, estimator RMS 0.19°.
- M17B 60° yaw at t = 40 s: post-slew settle +31.6 s, final 0.042°.
- M17D Monte Carlo (100 seeds): 100/100 converge; worst final 0.125°;
  worst rate 0.0158°/s; saturated fraction ≈ 0.39 early-transient (flagged,
  never clipped silently). No runs hidden.

## Regression

New `cpp/gnc/` + tests/demo only. Full suite 263/263. Ruff clean.

## Bugs found (fixed)

- Catch2 `&&` inside assertion (MSVC static assert) — parenthesized.
- Static audit opened headers relative to CWD: fails under ctest
  (WORKING_DIRECTORY = build). Fixed with source/build candidate paths.
- Ruff import sort in the oracle — fixed.

## Limitations

Attitude-only; no translation loop (M18+); no sensor dropouts in scenarios
(M20); instantaneous timeline switches; single bus / single orbit regime.
