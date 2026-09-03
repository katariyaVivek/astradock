# Lesson 017 — Integrated Closed-Loop GNC

## 1. Physical Problem & Motivation

M16 proved each control law converges fed perfect state. Flight never offers
perfect state — it offers a navigation estimate with noise, bias, latency, and
dropouts. M17 wires the full loop for the first time:

```text
bus -> sensors -> navigation -> guidance -> control -> actuators -> bus
```

with the inviolable rule: **guidance and control see estimates, never the
simulated state**. The controller input struct carries only `estimated_*` and
`reference_*` fields; a simulated-state-typed member would be a compile-time
shape break, and a static source audit test scans the seam headers for such
declarations. Truth appears in exactly two roles: inside sensor models
(measurement generation, M12) and inside telemetry scoring (harness-side
columns explicitly named `scored_*`), never as a controller input.

## 2. Loop Ticks and State Variables

Per 20 Hz tick: gyro sample + MEKF predict; star tracker (10 Hz) + MEKF update;
timeline reference (hold, or 60° yaw step at t = 40 s); estimate-based PD torque
(M16A law on estimated attitude/rate); per-axis saturation vs the 0.2 Nm wheel
envelope (M16E accounting); wheel-assembly step (M14, first-order motor lag
τ = 20 ms); achieved torque into `rk4_step_spacecraft` (M10). Telemetry per
tick: simulated/reference/estimate errors, rate, desired/achieved torque,
saturation, sensor validity, NIS, wheel speed, total wheel momentum, covariance
(via filter). Units/frames inherit M13/M16/M14 (SI, BODY torque, body-to-ECI).

## 3. Governing Equations

No new control math — M16A/M16E/M14/M13 equations composed. New elements are
the seam (`EstimateBasedAttitudeCommand`, `step_estimate_based_attitude_control`,
`AttitudeReferenceTimeline`) and the scoring (`summarize_closed_loop_run` with
latching settle: first index from which all subsequent errors hold the bound).

## 4. Assumptions and Failure Modes

Attitude-only loop (translation is M18+); single bus, circular 500 km orbit;
gyro bias estimated online from zero init; star tracker always available in
these scenarios (dropout robustness is FDIR/M20 scope); timeline switches are
instantaneous (slew transient is the controller's measured job). Failure
accounting: saturation flagged per tick; estimator divergence would appear as
NIS/NEES blowup in telemetry (monitored, not papered over).

## 5. Verification

7 cases: tick≡M16 law on identical inputs; saturation clamp; timeline selection
+ empty rejection; static state-isolation audit (declaration scan, build-dir
aware); perfect-estimate loop detumble + hold; MEKF-in-loop convergence
(attitude + bias); metrics latching + rejection; full 6-DOF wheel-actuated hold
with orbit untouched (decoupling). Independent Python oracle cross-checks tick,
timeline, settle. Demo: hold settles 80.8 s / 0.038°; slew settles +31.6 s /
0.042°; Monte Carlo 100/100 converge, worst 0.125°, worst rate 0.0158°/s.
Three figures.

## What you should now understand

1. Why must the controller struct carry estimates, not just "not truth"?
2. What two roles may simulated state play, and which role is forbidden?
3. Why does the audit scan declarations, not words?
4. Why settle latching instead of first-crossing?
5. Why did closed-loop settle slower (80.8 s) than open-loop PD (60 s)?
6. What would NIS blowup during a run tell you?
