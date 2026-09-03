# M15 Rotating-Frame Kinematics & Relative Dynamics — Validation Report

## Scope

`cpp/frames/lvlh_rate.hpp` (LVLH rate, angular acceleration, transport theorem),
`cpp/relative/relative_state.hpp` (relative state, CW prediction/propagation),
tests (`tests/cpp/test_relative_dynamics.cpp`, 11 cases / 69 assertions), demo
(`tools/relative_dynamics_demo.cpp`), oracle
(`python/audit/independent_relative_reference.py`), telemetry
(`data/m15_*.csv`), figures (`artifacts/figures/m15_*`), lesson
(`docs/lessons/015_relative_dynamics.md`).

## Analytical checks (all pass)

- Circular 500 km: `omega = +Z n`, matches `compute_circular_orbit_reference` to 1e-12.
- `omega = h/r^2` with 500 m/s radial-rate independence; inclined direction ∥ h.
- Transport hand case (`w = z, v = x` → `w x v = y`); inverse split recovers zero.
- Coriolis `2z x y = -2x`; centrifugal `z x (z x x) = -x`; Euler `3z x x = +3y`.
- Full composition `[−1.9, −7.0, 0]` term-by-term to 1e-12.
- Circular `alpha = 0` to 1e-15; rate-consistency residual < 1e-9 (also eccentric).
- ECI↔LVLH round trip to 1e-12 (pos) / 1e-9 (vel); zero-separation exact.
- CW: drift hold, radial secular growth, cross-track oscillator, bounded ellipse
  closes over one period (`x → 100 m`, `y → 0`), `Phi(0) = I`.
- Forced −Y push moves chaser −Y vs unforced arc; unforced RK4-vs-closed-form 1e-6.

## Independent verification

Pure-Python oracle (rate, basis, relative state, CW prediction); rate, ellipse
closure, and relative-state spot checks pass. No production code imported.

## Scenario validation (M15E breakdown)

Bounded-ellipse CW vs two-body ECI truth, 500 km circular reference:

| sep | 1-orbit err | 3-orbit err |
|---|---|---|
| 100 m | 0.027 m | 0.082 m |
| 500 m | 0.69 m | 2.06 m |
| 1 km | 2.74 m | 8.22 m |
| 5 km | 68.5 m | 205.5 m |
| 10 km | 274 m | 822 m |
| 25 km | 1713 m | 5139 m |

Error ∝ sep² (100× sep → 10⁴× err) and ∝ horizon (3 orbits ≈ 3× 1 orbit).
Operational reading: CW valid sub-percent to ~1 km; degraded at 5 km (1.4%);
invalid as a predictor at ≥10 km. 5 km overlay: 68.2 m final error.

Two review-driven corrections during development (wrong expected sign in the
along-track-offset velocity test; mislabeled derivative slots in the CW RHS
test) were root-caused to test-side errors — implementation matched the
independent oracle in both cases.

## Regression

No existing-file behavior changed (new `cpp/frames/lvlh_rate.hpp`,
`cpp/relative/` only; M06 APIs untouched). Full suite 245/245.

## Limitations

Circular-reference CW only (no eccentric, J2-differential, or drag-differential
terms); short-to-moderate horizons; LVLH undefined for collinear states.
