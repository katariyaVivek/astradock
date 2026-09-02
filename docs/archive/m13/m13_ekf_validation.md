# Milestone M13 Validation Report — State Estimation & Extended Kalman Filter (M13A)

## 1. Executive Summary

M13A delivers the first AstraDock state-estimation system: a 6-state
translational Extended Kalman Filter estimating ECI position and velocity from
simulated M12 GNSS measurements. The estimator consumes measurements only;
truth is never visible to the filter. Verification follows the mandated three
layers: analytical Kalman cases, simulation/statistical consistency, and an
independent Python oracle.

**Scope statement:** M13A is implemented and validated. M13B (IMU-specific-
force propagation) and M13C (attitude error-state estimation, range update,
integrated multi-sensor filter) are NOT yet implemented; no claims are made
for them. ROADMAP carries M13 as IN PROGRESS.

## 2. Estimator Design

```text
State ordering (ECI, SI): x = [rx ry rz vx vy vz]  (m, m/s)
Covariance P: 6x6, same ordering (m^2, m^2/s^2, m^2/s cross terms)
Dynamics: two-body point-mass gravity (identical to truth model -> zero
          mismatch baseline)
Mean propagation: RK4 (bitwise-identical to M04 open-loop propagation)
Covariance transition: Phi = I + F dt, F = [[0,I],[G,0]],
          G = -mu/r^3 (I - 3 rhat rhat^T), first-order discretization
Process noise: continuous white-noise acceleration, sigma_a = 1e-3 m/s^2,
          Q = sigma_a^2 [[dt^3/3 I, dt^2/2 I],[dt^2/2 I, dt I]]
Measurement: GNSS z = [r; v] + n, H = I6, isotropic R
Update: innovation y, S = HPH^T+R, K = PH^T S^-1 via Cholesky solve,
          Joseph-form posterior, documented symmetrization
Numerical policy: NaN/Inf rejection, relative symmetry tolerance,
          negative-diagonal rejection, non-SPD S raises (no silent repair)
```

Scenario: circular 500 km orbit, RK4 dt=1 s; GNSS 1 Hz, sigma_r = 10 m,
sigma_v = 0.05 m/s; deliberate initialization error ~61.6 km / 6.16 m/s;
P0 = diag(25 km)^2 x3, (5 m/s)^2 x3.

## 3. Layer 1 — Analytical Verification

| Case | Result |
| --- | --- |
| Scalar KF update (x=0,P=1,z=2,R=1) | K=0.5, x+=1, P+=0.5, NIS=2 exact |
| Two-state scalar-measurement update | K=[0.5,0.25], x+=[1,1.5], P+=[[0.5,0.25],[0.25,0.875]] to 1e-15 |
| Exact covariance prediction (CV model) | matches hand value [[6.4167,5],[5,11]] to 1e-15 |
| Cholesky SPD solve hand case | x=[-0.125,1.75] exact; non-PD/singular throw |
| Joseph vs short form agreement | equal in exact arithmetic (verified) |

## 4. Jacobian Audits

Gravity Jacobian G = da/dr verified against independent central finite
differences at 4 non-axis-aligned LEO positions with step sweep
h in {1e6 ... 0.1} m:

| Audit | Worst best relative error |
| --- | --- |
| C++ test suite (Frobenius relative) | < 1e-10 (asserted), typical ~1e-11 |
| Python oracle (independent NumPy FD) | 5.87e-11 |

Step sweep documents the truncation O(h^2) vs roundoff eps/h tradeoff;
optimal steps fall near h ~ 1e3-1e4 m. Analytical structure checks: radial
eigenvalue +2mu/r^3, transverse -mu/r^3, body-diagonal zero-diagonal case,
exact symmetry, tracelessness to 1e-18 (floating-point unit-vector limit).

GNSS H matrix = I6 verified exactly (identity model).

## 5. Layer 2 — Simulation & Statistical Consistency

### Nominal convergence run (6000 s, dropout [3000,3900] s)

| Metric | Value |
| --- | --- |
| Initial position error | 61 644.1 m |
| Final position error | 1.28 m |
| Final velocity error | 0.0080 m/s |
| EKF position RMSE (2nd half) | 1.38 m |
| Raw GNSS position RMSE | 4.00 m |
| EKF velocity RMSE (2nd half) | 0.0094 m/s |
| Raw GNSS velocity RMSE | 0.282 m/s |
| Mean NIS (valid updates) | 5.94 (chi-square(6): mean 6) |
| Mean NEES (2nd half) | 4.68 (chi-square(6) band [1.64,14.45]) |
| Position sigma growth during dropout | 0.70 m -> 18.59 m (monotone) |
| Max position error during dropout | 15.9 m (no divergence) |

The EKF beats raw GNSS by ~2.9x in position and ~30x in velocity RMSE with a
zero-mismatch model and 1 Hz updates.

### Monte Carlo study (100 seeds x 2000 s)

| Metric | Value |
| --- | --- |
| Median final position error | 1.07 m |
| Median final NEES | 3.12 (inside chi-square(6) 95% band) |
| Mean converged NIS across runs | 5.96 |
| Mean EKF / raw position RMSE | 1.17 m / 3.99 m |
| Mean EKF / raw velocity RMSE | 0.0083 / 0.282 m/s |
| Runs beating raw GNSS (position) | 100 / 100 |

Interpretation: the filter is consistent (NIS ~ 6, NEES within band) and
mildly conservative rather than overconfident — acceptable and safe for GNC.
Seed variation covers sensor-noise realization; initialization offset is held
fixed by design and its dispersion is exercised separately in the unit suite.

### Robustness and contract tests

| Test | Result |
| --- | --- |
| GNSS dropout covariance growth monotone over 900 s | PASS |
| Recovery shrinks trace below 1.10x pre-dropout within 300 s | PASS |
| No divergence during outage (< 500 m bound; actual 15.9 m) | PASS |
| Bitwise determinism (same seed -> identical telemetry) | PASS |
| Truth non-interference (estimator never mutates inputs; const-by-design + regression scenario identical reruns) | PASS |
| Noise limiting cases (R->0 snaps, R->large rejects; Q small/large gain ordering) | PASS |
| Invalid numerical states throw (NaN measurement/state/Phi, negative variance, singular S) | PASS |
| Predict-only mean bitwise equals open-loop RK4 (600 steps) | PASS |

## 6. Layer 3 — Independent Python Oracle

`python/audit/independent_ekf_reference.py` re-implements RK4, the gravity
Jacobian, Phi, Q, Cholesky solves, the Joseph update, and NIS from the
mathematical specification only (no AstraDock import), replays the C++
telemetry measurements, and compares full estimate histories:

| Comparison | Max discrepancy (last 100 steps) |
| --- | --- |
| Position | 8.88e-16 m |
| Velocity | 6.07e-18 m/s |
| Python mean NIS (converged) | 5.93 (vs C++ 5.94) |

Agreement is at floating-point roundoff level: two independent
implementations of the same mathematics produce the same filter.

## 7. Generated Artifacts

```text
data/m13_ekf_telemetry.csv        6001 samples, truth/measurement/estimate/
                                  covariance sigma/innovations/NIS/NEES
data/m13_ekf_monte_carlo.csv      100 per-seed summary rows
artifacts/figures/m13_position.png
artifacts/figures/m13_position_error.png
artifacts/figures/m13_velocity.png
artifacts/figures/m13_covariance.png
artifacts/figures/m13_innovations.png
```

## 8. Build and Test Scorecard

```text
Previous tests:            159 / 159 passing
New M13A tests:             18 / 18 passing (9721 assertions)
Total:                     177 / 177 passing
Clean C++20 build:         zero errors, zero warnings (/W4 or -Wall -Wextra
                           -Wpedantic -Wconversion -Wsign-conversion)
Ruff:                      PASS (python/audit + python/analysis clean)
Deterministic repeat run:  identical outputs
Truth non-interference:    PASS
Independent oracle:        PASS (< 1e-15 m disagreement)
```

## 9. Scope Guard Compliance

Implemented in M13A: translational estimation only. NOT implemented and NOT
claimed: guidance, control, actuators, rendezvous/docking logic, computer
vision, ML, reinforcement learning, optimal control, particle filters,
smoothing, factor graphs, SLAM, IMU-aided propagation (M13B), attitude
estimation (M13C), range update (M13C), outlier gating (optional future work).

## 10. Known Limitations

- Zero process-model mismatch by construction; J2/drag mismatch robustness is
  future work once M13B/C exist.
- First-order Phi discretization (valid: omega*dt ~ 1e-3).
- Single-rate GNSS-only updating; multi-rate fusion arrives with M13B/C.
- Monte Carlo randomizes sensor noise only, not initialization error.
- No innovation gating implemented (deliberate; see spec section 47).
