# Milestone M13 Validation Report — IMU-Aided Translational Propagation (M13B)

## 1. Executive Summary

M13B extends the M13A translational EKF with **IMU-aided prediction**: between
GNSS fixes the filter now dead-reckons on measured specific force instead of
re-integrating its own dynamics model. The central reconstruction chain

```text
a_I = C_I_B(q_hat) (f_m - b_a) + g(r_hat)
```

is implemented as one pure, independently audited function; prediction runs at
the 100 Hz IMU rate, GNSS corrects at 1 Hz. Verification follows the mandated
three layers; all previous M13A behavior is preserved bitwise.

**Scope statement:** M13A and M13B are implemented and validated. M13C
(attitude error-state estimation, range measurement update, integrated
multi-rate filter) and M13D are NOT implemented; no claims are made for them.
ROADMAP carries M13 as IN PROGRESS.

## 2. Known-Attitude Assumption (Explicit Staging)

- The attitude used for BODY -> ECI transformation is **externally supplied**
  through an explicit navigation-state type (`NavigationAttitudeEstimate`).
  It is treated as an estimate; no API accepts `SpacecraftState`, sensor
  truth containers, or any truth state.
- Scenarios propagate a reference attitude by pure quaternion kinematics in
  lock-step with truth (zero-error idealization); diagnostic tests inject
  deliberate attitude error to quantify sensitivity.
- The translational error-state covariance **excludes attitude columns**:
  F is still [[0, I], [G, 0]]. This is exact for the known-attitude
  idealization and benign in free fall (dC/d(delta-theta) * f = 0 when
  f ≈ 0), but understates uncertainty under sustained non-gravitational
  force with imperfect attitude. M13C closes this coupling.
- **No attitude estimation. No bias estimation.** Gyro samples flow through
  the pipeline and timestamps but are not consumed inside the filter;
  accelerometer bias is known configuration, deliberately mismatched in
  dedicated experiments to motivate bias states.

## 3. Measurement Model

Reuses the M12 sensor models unchanged:

```text
f_m = f_true + b_a + n_a      (BODY; f_true = non-gravitational acceleration)
w_m = w_true + b_g + n_g      (BODY; carried, not yet used in-filter)
z   = [r; v]_true + b + n     (ECI GNSS fix, H = I6)
```

Nominal configuration: sigma_a = 1e-3 m/s², sigma_r = 10 m,
sigma_v = 0.05 m/s; IMU 100 Hz, GNSS 1 Hz; dropout window [600, 700] s.

## 4. Propagation Equations

Mean (RK4 across each interval, ZOH on measured force and supplied
attitude; gravity varies continuously):

```text
dr/dt = v
dv/dt = C_I_B(q_hat(t_k)) (f_m - b_a) + g(r)
```

Covariance:

```text
Phi = I + F dt,  F = [[0, I], [G(r_hat), 0]],  G = -mu/r^3 (I - 3 rhat rhat^T)
P_prior = Phi P Phi^T + Q
```

Numerical policy unchanged from M13A: Cholesky solves (no explicit
inverses), Joseph-form update, relative symmetry tolerance, hard failure on
non-finite/non-SPD states. `predict_with_imu` degenerates **bitwise** to the
M13A `predict` when f_m − b_a = 0 (verified over 300 steps).

## 5. Process-Noise Model

Q derives from the configured accelerometer noise via the exact continuous
white-noise-acceleration discretization (`imu_prediction_process_noise`):

```text
Q = sigma_f^2 [ dt^3/3 I   dt^2/2 I ]
              [ dt^2/2 I   dt     I ]        (state order [r(0:2), v(3:5)])
```

Assumptions documented in-code: isotropic per-axis noise; unmodeled-
perturbation margin intentionally zero because non-gravitational accelerations
now enter THROUGH the measurement; RSS extension point documented; arbitrary Q
inflation forbidden. A wiring test proves the filter uses the IMU-driven Q and
not the M13A baseline sigma when both differ.

## 6. Jacobians

The acceleration map a(r) = C(q)(f − b) + g(r) has gravity as its only
position dependence, so da/dr = G exactly regardless of attitude or force.

| Audit | Result |
| --- | --- |
| C++ central finite differences across 4 attitudes x 2 LEO positions, step sweep h in {1e6..0.1} m | worst best-relative error < 1e-10 |
| Python oracle FD audit (independent NumPy) | 5.87e-11 |
| d(a)/d(v) | identically 0 (no velocity dependence) |

Step-sweep behavior documents truncation O(h²) vs roundoff eps/h tradeoff;
optimal steps near h ~ 1e3–1e4 m, consistent with M13A findings. Attitude
columns are excluded by design (Section 2).

## 7. Layer 1 — Analytical Tests

| Case | Result |
| --- | --- |
| DCM direction: q(90° about z) applied to body x̂ | maps to ECI ŷ exactly; catches conjugation errors |
| Known-attitude reconstruction vs independent textbook DCM (4 cases incl. identity, 120° about [1,1,1], small tilt) | < 1e-12 relative |
| predict_with_imu(f=0,b=0) vs M13A predict | bitwise equal over 300 steps |
| First-order small-dt limit (dt=1e-4) | dr = dt·v, dv = dt·a₀ within 1e-6 relative |
| Covariance wiring: IMU-sigma Q vs M13A-sigma Q | matches accelerometer-driven Q to < 1e-18 |
| GNSS update after IMU prediction vs direct ekf_update algebra | identical posterior |
| Bias drift laws | position linear-in-bias (ratio 5.000000 for 5x bias), quadratic-in-time (ratio ~4 at half horizon); flat-space 0.5bt² matched to ~1.7% |
| Attitude dependency | delta-a = 2 sin(delta-theta/2) |f| exact to 1e-12; 10 mrad on 0.1 m/s² gives 1e-3 m/s² corruption |

## 8. Free-Fall Orbital Test

Ideal circular 500 km orbit, ideal IMU (zero noise/bias), tumbling known
attitude, 100 Hz propagation for 60 s (6000 steps): the IMU-aided estimate
matches the open-loop two-body reference to **< 1e-6 m / 1e-8 m/s** at all
checkpoints (residuals are timestamp-subtraction roundoff that the filter
correctly honors). In free fall the chain reduces to a_I ≈ g(r̂): gravity
compensation verified end-to-end against the M04/M10 reference trajectory.

## 9. Non-Gravitational Force Test (Drag)

M11 drag environment (500 kg, Cd 2.2, 2 m², rho_ref 4e-12 kg/m³ @500 km):
truth propagates WITH drag (~0.89 µm/s²); a high-grade accelerometer
(sigma_f = 1e-8 m/s² — resolving drag-class forces demands ng-class
instrumentation, an explicit lesson) feeds the estimator.

| Metric | Value |
| --- | --- |
| Measured \|f(t0)\| vs true drag magnitude | within noise (norm preserved by rotation) |
| IMU-aided prediction error after 400 s | < 1 m (sub-metre class) |
| Gravity-only coasting error after 400 s | 76.19 m (quadratic divergence) |

The drag information reaches the estimator ONLY through the simulated
accelerometer; nothing injects the known drag acceleration into the filter.

## 10. GNSS+IMU vs GNSS-Only (Monte Carlo, 100 seeds x 300 s)

Identical truth and identical GNSS realizations per seed; both filters
consume the same fixes. Initial error randomized from P0-scale Gaussians.

| Metric (mean across seeds) | GNSS-only | IMU+GNSS |
| --- | --- | --- |
| EKF position RMSE (2nd half) | 1.247 m | 1.248 m |
| Raw GNSS position RMSE | 3.995 m | 3.995 m |
| EKF velocity RMSE | 0.009 m/s | 0.009 m/s |
| Mean converged NIS | 5.979 | 5.979 |
| Mean NEES (2nd half) | 4.266 | 4.293 |
| Runs beating raw GNSS | 100/100 | 100/100 |
| Median final NEES (IMU+GNSS) | — | 3.773 (chi-square(6) band [1.64, 14.45]) |

Interpretation (honest reporting): with continuous 1 Hz fixes and a
zero-mismatch model, both predictors carry equivalent information and
performance remains similar BY EXPECTATION — the IMU adds no information the
dynamics model does not already provide in this artificial configuration.
The IMU advantage is demonstrated where the model alone fails: dropouts
(Section 11) and drag (Section 9).

## 11. GNSS Dropout Experiment

Nominal telemetry run (1200 s, dropout [600, 700] s):

| Metric | Value |
| --- | --- |
| Initial position error | 61 644.1 m |
| Final position error | 0.88 m |
| Final velocity error | 0.0044 m/s |
| Position RMSE (2nd half) | 1.354 m (raw GNSS: 3.975 m) |
| Velocity RMSE (2nd half) | 0.0096 m/s (raw: 0.281 m/s) |
| Mean NIS | 5.84 (chi-square(6): mean 6) |
| Mean NEES (2nd half) | 4.35 |
| Max position error during outage | 1.56 m (dead reckoning holds) |
| Max velocity error during outage | 0.0022 m/s |
| Position sigma pre / late-dropout / post-recovery | 1.217 m -> 2.259 m -> 1.293 m |
| Covariance trace growth during outage | strictly monotone (asserted every step) |

During the outage only IMU prediction runs; the estimate stays within metres
of truth purely by inertial integration, and the first returning fix snaps
both estimate and covariance back.

## 12. Bias Sensitivity Study

Prediction-only, identity attitude, t = 300 s, uncalibrated true bias along
body x (data/m13b_bias_sensitivity.csv):

| true bias | pos drift | vel drift | 0.5·b·t² |
|---|---|---|---|
| 0 | 0.000 m | 0.000 m/s | 0 |
| 1e-4 m/s² | 4.578 m | 0.031 m/s | 4.50 m |
| 5e-4 m/s² | 22.890 m | 0.155 m/s | 22.50 m |

Position drift grows linearly with bias and quadratically with time; the
~1.7% excess above the flat-space law is radial gravity-gradient
amplification. This quantifies why M13C/D require bias states.

## 13. Covariance Behavior

- Dropout: monotone growth asserted step-by-step (trace strictly increasing).
- Recovery: contraction below late-outage values within one update cycle;
  not every component shrinks monotonically (documented; geometry and Q shape
  P's evolution).
- Growth magnitude consistent with theory: velocity random walk
  sigma_v(T) = sigma_f sqrt(T), position variance ~ sigma_f² T³/3.

## 14. Independent Python Oracle

`python/audit/independent_imu_ekf_reference.py` re-implements quaternion DCMs,
the inertial reconstruction equation, RK4 over the IMU ODE, gravity Jacobian,
Phi, Q, Cholesky solves, Joseph update, and NIS from the specification only
(no AstraDock import). It replays data/m13b_oracle_segment.csv (6002 full-rate
records, 61 GNSS fixes) and audits rotation direction and Jacobians
independently:

| Comparison | Result |
| --- | --- |
| Rotation-direction audit (90° z) | PASS |
| Jacobian FD audit (worst best relative) | 5.87e-11 |
| Max position discrepancy (last 100 steps) | **4.44e-16 m** |
| Max velocity discrepancy (last 100 steps) | 5.68e-14 m/s |
| Python mean NIS | 5.35 (61 updates, short segment) |

Two independent implementations of the same mathematics agree at machine
precision.

## 15. Numerical Sensitivity

- Timestep refinement: the RK4 mean shows fourth-order accuracy (first-order
  analytic limits hold at 1e-6 relative for dt = 1e-4; free-fall nm-level
  agreement at dt = 10 ms); the covariance discretization remains
  first-order in Phi by design (omega*dt ~ 1e-3 makes neglected terms O(1e-6)).
- Timestamp jitter (+-1 ms deterministic wobble): intervals honored exactly;
  final-state deviation bounded by |v|·delta-T as predicted (>0, <10 m for
  the tested profile).
- Per-step Q scalings verified exactly (dt³/dt²/dt blocks).

## 16. Failure Cases

| Case | Behavior |
| --- | --- |
| NaN force / NaN gyro / NaN dt | domain_error thrown, never propagated |
| dt <= 0 (out-of-order timestamps) | rejected loudly |
| Non-unit attitude quaternion (norm 2) | rejected |
| NaN configured bias | rejected at validation |
| IMU unavailable | harness falls back to documented dynamics-only predict(dt); no invented measurements; output finite and reference-consistent |
| Large uncalibrated bias | drifts quadratically but never diverges to NaN/Inf |
| Incorrect assumed attitude | bounded, quantified error (Section 7); no instability |

## 17. Truth Non-Interference & Determinism

- No estimator API accepts any truth container (header-level architecture;
  imu_prediction.hpp/translational_ekf.hpp include no spacecraft/sensor
  truth types). Regression scenario confirms truth telemetry bitwise
  untouched after estimation.
- Full demo rerun produces bitwise-identical CSVs (telemetry, oracle
  segment, bias study, Monte Carlo) and console output.

## 18. Build and Test Scorecard

```text
Previous tests:            178 / 178 passing
New M13B tests:             19 / 19 passing (1 697 646 assertions)
Total CTest:               197 / 197 passing
Clean C++20 build:         zero errors, zero warnings (/W4 permissive-)
Ruff:                      PASS (all python/ sources clean)
Deterministic repeat run:  bitwise identical outputs
Truth non-interference:    PASS
Independent oracle:        PASS (< 1e-15 m disagreement)
Scope guard:               PASS (see Section 20)
```

## 19. Generated Artifacts

```text
data/m13b_imu_ekf_telemetry.csv       12 001 rows @10 Hz incl. all GNSS epochs
data/m13b_oracle_segment.csv          6 002 rows @100 Hz full-rate replay segment
data/m13b_bias_sensitivity.csv        3 bias scenarios + analytic references
data/m13b_imu_ekf_monte_carlo.csv     200 rows (100 seeds x 2 cases)
artifacts/figures/m13b_position.png            (Plot A)
artifacts/figures/m13b_velocity.png            (Plot B)
artifacts/figures/m13b_dropout_covariance.png  (Plot C)
artifacts/figures/m13b_specific_force.png      (Plot D)
artifacts/figures/m13b_error_envelopes.png     (Plot E)
artifacts/figures/m13b_innovations.png         (Plot F)
docs/lessons/014_imu_aided_navigation.md
```

## 20. Scope Guard Compliance

Implemented in M13B: IMU-aided translational prediction with externally
supplied attitude, gravity compensation, multi-rate loop, dropout/bias/
attitude-dependency experiments, extended oracle and plots.

NOT implemented and NOT claimed: attitude error-state estimation, quaternion
EKF, star-tracker update, gyro or accel bias ESTIMATION, range measurement
update, integrated attitude/navigation filter, outlier gating beyond basic
NIS diagnostics, fault detection, guidance, control, rendezvous, docking,
computer vision, machine learning. Those belong to M13C/D/M14+.

## 21. Known Limitations

- Attitude is externally supplied and assumed perfect in nominal scenarios;
  translational covariance ignores attitude uncertainty (exact only for the
  known-attitude idealization; benign in free fall, not under sustained thrust).
- Accelerometer bias must be perfectly calibrated; unestimated bias drifts
  the solution quadratically (quantified in Section 12).
- Drag-class specific forces require ng-class accelerometers; the nominal
  1e-3 m/s² unit cannot resolve them (Lesson 014 discusses instrument classes).
- Monte Carlo randomizes sensor noise and initial error; attitude kept known
  per spec.
- No innovation gating (deliberate; deferred).
