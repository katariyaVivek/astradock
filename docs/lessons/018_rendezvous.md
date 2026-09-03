# Lesson 018 — Rendezvous & Proximity Operations

## 1. Physical Problem & Motivation

M17 closed the attitude loop. M18 flies the *relative trajectory*: from a 6 km
offset to a 50 m hold point along the approach axis, tracking waypoints with
the M16 translation law while safety geometry watches every tick. Two findings
from development carry the teaching weight:

1. **Known-actuation awareness**: the spacecraft commanded its thrust, so the
   navigation predict must propagate it. A two-body-only EKF lags under
   sustained thrust (0.05–0.2 m/s² vs 1e-3 process noise) and guidance
   under-brakes into a flyby. The demo folds commanded thrust into the filter
   mean by exact double-integral correction.
2. **Epoch alignment**: forming the LVLH relative state from a freshly stepped
   target and a stale chaser injects a ~7.6 km bias (one tick of orbital
   motion). Scenario sequencing propagates both vehicles and the filter to the
   new epoch *before* differencing — a mission-software bug class, caught by
   estimate-vs-scored telemetry, not by unit tests.

## 2. Guidance Construction and Safety

Per leg A→B: line-of-sight unit, closing profile `min(cruise, √(2·a_brake·d))`
(cruise far, constant-deceleration braking near), first-order velocity tracking
`a_track = (v_des − v)/τ`, then the M16 PD + CW feedforward around it. Safety
on the estimated state: keep-out sphere (25 m, hard), corridor (50 m lateral
within 1000 m), closing-speed warning (3 m/s). Legs advance monotonically on
capture; retreat is abort (M19). Waypoints static in LVLH (short-arc validity).

## 3. State Variables, Inputs, Outputs

Estimated LVLH relative state in; LVLH acceleration command out (allocated to
body thrusters via live attitude downstream). Telemetry: range, lateral error,
closing speed, accel, leg, flags, estimate-vs-scored gap, Δv proxy ∫|a|dt.

## 4. Verification

7 cases / 47 assertions: geometry hand values, predicate boundaries (exactly at
keep-out radius = no violation, strict `<`), capture advance + completion,
cruise/braking branches, decreasing sequence speeds, CW-plant tracking to each
leg, invalid rejection. Independent Python oracle (4 functions). Demo: CW sweep
nominal + 5% under-burn (both complete, Δv ≈ 16.7 m/s) + corridor stress; ECI
truth with GNSS relative nav completes to 55.5 m (estimate-captured at 4.9 m
vs 5 m ball, truth 0.45 m outside — sub-meter, docking needs better nav);
60 s dropout completes to 55.6 m. Est-vs-scored gap ≤ 3.6 m throughout.

## What you should now understand

1. Why must the navigation predict know the commanded thrust?
2. What happens when relative states mix epochs?
3. Why does the closing profile use `√(2·a·d)` near the waypoint?
4. Why is capture declared on the estimate, and what does that cost truth?
5. When does the corridor constraint apply, and when is it silent?
6. Why did 100 m off-axis + 10 m ball + 1 m/s flyby fail by geometry?
