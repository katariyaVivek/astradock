# M18 Rendezvous & Proximity Operations — Validation Report

## Scope

`cpp/rendezvous/rendezvous_guidance.hpp` (sequence, profile, predicates, tick),
`tests/cpp/test_rendezvous.cpp` (7 cases / 47 assertions), demo
(`tools/rendezvous_demo.cpp`: CW sweeps + ECI truth ± dropout), oracle
(`python/audit/independent_rendezvous_reference.py`), telemetry
(`data/m18_*.csv`), figures (`artifacts/figures/m18_*`), lesson
(`docs/lessons/018_rendezvous.md`).

## Analytical checks (all pass)

- Lateral √(30²+40²) = 50; closing sign convention ±; lateral velocity ignored.
- Keep-out / corridor / speed-warning predicates incl. strict-`<` boundary.
- Capture advance + monotonic legs + mission completion.
- Profile: 1.0 m/s cruise at 4000 m; √(2·10⁻³·4) = 0.089 m/s braking at 4 m.
- Sequence marches −5000→−50 m with decreasing limits; final hold outside
  keep-out. CW-plant tracking reaches each leg (±15 m).

## Independent verification

Pure-Python oracle (lateral, closing, profile, capture); 5 spot checks pass.

## Scenario validation

- CW sweep nominal: complete, Δv proxy 16.70 m/s. 5% under-burn: complete,
  16.71 m/s (feedback absorbs the error — robustness, not luck).
- Corridor stress (60 m lateral start): complete, zero flags (offset beyond
  corridor range is correctly silent; tracker nulls it on approach).
- ECI truth + GNSS relative nav: complete, final 55.5 m; est-vs-scored gap
  ≤ 3.6 m; capture on estimate (4.9 m in a 5 m ball, truth 0.45 m outside).
- 60 s GNSS dropout (t = 200–260): complete, final 55.6 m — thrust-aware
  predict bridges the outage; no divergence.

## Development bugs found (all fixed, documented as findings)

1. `chaser_eci_from_relative` with zero velocity builds a flyby
   (−ω×ρ ≈ 1.3 m/s), not a hold — same-velocity ECI construction used instead.
2. Unmodeled-thrust EKF lag → under-braking → waypoint flyby — fixed by
   thrust-aware mean correction.
3. Off-axis init + tiny ball + fast flyby is geometrically impossible —
   on-axis start (realistic far-field condition).
4. Target/estimate epoch mismatch (~7.6 km bias) — epoch-aligned sequencing.

## Regression

New `cpp/rendezvous/` + tests/demo only. Full suite 270/270. Ruff clean.

## Limitations

Static LVLH waypoints (short arcs); cooperative static target; no abort
trajectory (M19); 5 m GNSS noise bounds capture accuracy (M19 needs better
relative nav); no plume/obstacle modeling.
