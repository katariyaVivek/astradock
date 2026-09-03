# M19 Autonomous Docking — Validation Report

## Scope

`cpp/docking/docking.hpp` (ports, alignment, contact, acceptance, abort),
`tests/cpp/test_docking.cpp` (7 cases / 44 assertions), demo
(`tools/docking_demo.cpp`: nominal, fast, lateral, tilted, dropout), oracle
(`python/audit/independent_docking_reference.py`), telemetry
(`data/m19_*.csv`), figures (`artifacts/figures/m19_*`), lesson
(`docs/lessons/019_docking.md`).

## Analytical checks (all pass)

- Port position: identity offset; 90° yaw maps +X into +Y; invalid rejection.
- Alignment: 5° tilt reads 5°; double cover invariant; envelope/abort lines.
- Contact: 0.05 m @ 0.1 m/s → 60 N; separating fast unloads to 0 (no pull);
  rest preload 20 N.
- Acceptance: perfect mates latch; each single failure (lateral, closing,
  tilt, rate) rejects independently.
- Aborts: lateral / closing / attitude / timeout each fire with reason codes.
- Latch banking: acceptance without banked time does not latch.
- Axial sign: +0.5 apart (no contact), −0.05 penetration with k·pen + c·rate.

## Independent verification

Pure-Python oracle (ports, alignment, contact, acceptance); 4 spot checks pass.

## Scenario validation

- Nominal (50 m, on-axis): latches, no abort; settled lateral 0.06 m.
- Fast (3× cruise, envelope-scaled): latches — correctly, since acceptance
  enforces the envelope every tick; overspeed can never latch, only abort.
- 2 m lateral: tracker recovers, latches (final 0.14 m).
- 20° tilt: aborts `attitude_exceeded` on tick 1 (deterministic, instant).
- Dropout (60 s): latches (thrust-aware predict bridges).

## Development bugs (fixed, kept as findings)

1. Axial sign triple-site (separation vs penetration, closing = −ds/dt).
2. Port-offset lateral vs range-gated abort (envelopes live at their range).
3. Coarse nav + contact-unaware guidance → crush (sensor requirement +
   contact freeze).

## Regression

New `cpp/docking/` + tests/demo only. Full suite 277/277. Ruff clean.

## Limitations

Rigid ports; axial-only penalty contact (lateral contact aborts, no sliding);
relative rate is a frame-mismatched norm bound; static target; retreat is M21.
