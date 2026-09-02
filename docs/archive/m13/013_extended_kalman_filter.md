# Lesson 013 — State Estimation & the Extended Kalman Filter

## 1. The Physical Problem: Truth vs Measurement vs Estimate

In every previous milestone we had direct access to the spacecraft's **truth
state**. This is the simulator's private knowledge. A real navigation computer
never has it. It receives only **imperfect sensor measurements**:

```text
GNSS          position/velocity corrupted by noise (and bias)
Accelerometer specific force f = a - g, corrupted by bias and noise
Gyroscope     body rates corrupted by bias and noise
Star tracker  attitude corrupted by small random rotations
```

From these signals the computer must **infer** position, velocity, attitude,
angular velocity, and sensor biases. That inferred belief is the **estimate**.
Three different quantities must never be confused:

```text
Truth       = physical reality (simulator-only, never enters the filter)
Measurement = what sensors reported (noisy, biased, rate-limited)
Estimate    = the filter's probabilistic belief about truth
```

M13's central engineering rule:

> **The estimator must never see truth.**

Truth-dependent statistics (estimation error, NEES) are computed only by
evaluation harnesses outside the estimation module.

## 2. Why a Linear Kalman Filter Is Not Enough

The classical Kalman filter is optimal only for **linear** systems driven by
zero-mean Gaussian noise. Spacecraft navigation is nonlinear:

- two-body gravity `a = -mu * r / |r|^3` is nonlinear in position,
- quaternion attitude kinematics are nonlinear,
- range measurement `rho = ||r_target - r_sc||` is nonlinear in position,
- star-tracker measurements act on SO(3), not on a flat vector space.

The **Extended Kalman Filter** keeps the KF recursion but linearizes the
nonlinear models about the current estimate using Jacobians:

```text
x_{k+1} = f(x_k) + w_k        w ~ N(0, Q)   process model + noise
z_k     = h(x_k) + v_k        v ~ N(0, R)   measurement model + noise
```

## 3. Covariance Semantics

The filter maintains a Gaussian belief `N(x_hat, P)`:

```text
P = E[(x_true - x_hat)(x_true - x_hat)^T]
```

- **Diagonal entries** are variances; `sigma_i = sqrt(P_ii)` is 1-sigma
  uncertainty of state component i.
- **Off-diagonal entries** are correlations between component uncertainties.
- P must remain **symmetric** (`P = P^T`). Roundoff erodes this, so AstraDock
  symmetrizes `P <- (P + P^T)/2` after justified operations and *documents*
  each place it does so. Arbitrary entry clamping is forbidden.
- Physical covariance must be **positive semidefinite**. A negative diagonal
  variance means the filter is broken; the M13 implementation throws rather
  than repairs silently.
- Symmetry tolerance must be **relative** to matrix scale: position variances
  reach `(7.5e4 m)^2 ~ 5.6e9 m^2` while velocity variances are O(10 m^2/s^2),
  so an absolute tolerance would reject physically valid matrices.

## 4. Q Is Model Uncertainty; R Is Sensor Uncertainty

```text
Q = process noise   = "how wrong can the dynamics model be between
                       measurements?"  (unmodeled accelerations)
R = measurement noise = "how noisy is this sensor?"  (GNSS jitter)
```

They are not interchangeable tuning knobs. Large Q makes the filter distrust
its own prediction and follow measurements. Large R makes it coast on its
model. The M13A Q uses a continuous white-noise acceleration model with
`sigma_a` in m/s^2, discretized exactly for the double integrator:

```text
Q = sigma_a^2 * [ dt^3/3 I   dt^2/2 I ]
                [ dt^2/2 I   dt     I ]
```

with `sigma_a = 1e-3 m/s^2`, a reasonable LEO unmodeled-acceleration scale.

## 5. Prediction

```text
x_prior = f(x_posterior)                 mean: full RK4 through real dynamics
P_prior = Phi P Phi^T + Q                covariance: linearized transition
Phi ~= I + F dt                          first-order discretization
F = [[0, I], [d a/d r, 0]]               continuous dynamics Jacobian
```

Validity of `Phi ~= I + F dt`: requires `omega*dt << 1`. In LEO at dt = 1 s,
`omega*dt ~ 1e-3`, so neglected terms are O(1e-6). Documented approximation;
a second-order term could be added later without changing the architecture.

Note the deliberate asymmetry: the **mean** uses the exact same RK4 machinery
as the truth propagator (verified bitwise against M04 in the tests), while the
**covariance** accepts first-order discretization error.

## 6. The Gravity Jacobian

For `a(r) = -mu r / |r|^3`:

```text
d a_i / d r_j = -mu ( delta_ij / r^3 - 3 r_i r_j / r^5 )
G = d a / d r = -mu/r^3 (I - 3 rhat rhat^T)
```

Physical reading of its eigenstructure:

- radial direction: eigenvalue `+2 mu/r^3` — radial separation grows,
- transverse directions: eigenvalue `-mu/r^3` — tangential separation oscillates,
- traceless (`trace(G) = 0`): Laplace's equation away from the origin.

This analytical Jacobian is audited against central finite differences over
multiple non-axis-aligned positions and a step sweep (truncation O(h^2) vs
roundoff ~ eps/h); best relative error observed: 5.9e-11. Finite-differencing
the production function alone would prove nothing.

## 7. Measurement Update

```text
y      = z - h(x_prior)                    innovation: sensor saw minus
                                           filter expected
S      = H P H^T + R                       innovation covariance
K      = P H^T S^{-1}                      Kalman gain (via Cholesky solve,
                                           never explicit inverse)
x_post = x_prior + K y
P_post = (I-KH) P (I-KH)^T + K R K^T       Joseph form
```

Why Joseph form: the short form `(I-KH)P` is algebraically equal but can lose
symmetry/positive semidefiniteness when K carries roundoff. Joseph preserves
both by construction. For GNSS in M13A, `h(x) = Hx = x`, so `H = I_6`.

## 8. NIS — Innovation Consistency (online)

```text
NIS = y^T S^{-1} y ~ chi-square(m)
```

For df = 6 (GNSS position+velocity): mean 6, std sqrt(12).
Measured mean over converged windows: 5.94 (C++), 5.93 (Python oracle).
A single number available onboard without truth; the primary health signal.

## 9. NEES — Estimation Consistency (evaluation only)

```text
NEES = e^T P^{-1} e,   e = x_true - x_hat ~ chi-square(n=6)
```

> A filter can have low error while still being statistically overconfident.

Measured: mean NEES 4.68 (nominal run second half), Monte Carlo median 3.12 —
within the chi-square(6) central band [1.64, 14.45], mildly conservative.
Computed in `estimation/diagnostics.hpp`, explicitly outside the filter.

## 10. Multi-Rate Handling and Dropouts

The filter supports predict-only and update-only steps; any sensor that
produces a valid timestamped sample can update between predictions. During a
900 s GNSS dropout in the demo run:

```text
position sigma: 0.70 m -> 18.59 m   (pure prediction, monotone growth)
max position error during outage: 15.9 m (no divergence)
after recovery: information reabsorbed, RMSE returns < 25 m
```

Covariance growth during information starvation and shrinkage after recovery
are qualitative behaviors every navigation engineer must recognize.

## 11. Why an Accurate Sensor Does Not Guarantee an Accurate Estimate

Filter performance depends on:

```text
model accuracy         (does f match reality?)
Q justification        (is sigma_a honest?)
R calibration          (does quoted noise match the true sensor?)
initial covariance P0  (honest initial uncertainty?)
initial state error    (how far off did we start?)
sensor timing          (timestamps, dropouts)
sensor bias            (modeled? estimated?)
Jacobian correctness   (F and H verified independently?)
```

A perfect sensor with a bad covariance or wrong Jacobian still produces a bad
filter — plausible trajectories are not proof of correctness. That is why M13
requires analytical tests + Jacobian audits + statistical consistency +
Monte Carlo + an independent Python reference before claiming anything.

## 12. Why We Do Not Put the Quaternion Directly Into a Naive EKF

A quaternion stores 4 numbers but represents only 3 rotational DOF under the
nonlinear unit constraint `||q|| = 1`. An additive update
`q_plus = q_minus + K*y` violates the constraint and produces invalid
attitudes. The principled approach is an **error-state** filter: carry the
unit quaternion as the nominal state, keep a 3-component small-angle error
`delta_theta` in the covariance, correct multiplicatively via
`q_plus = delta_q (x) q_nominal` with `delta_q ~ [1, delta_theta/2]`, then
reset the error state and transform the covariance. This belongs to M13C and
is documented here because it drives the whole state design.

## What You Should Now Understand

1. Why must truth, measurement, and estimate stay strictly separate?
2. What breaks in a linear KF applied to orbital dynamics?
3. What do Q and R mean physically, and why can they not be traded freely?
4. Why must covariance stay symmetric and positive semidefinite, and what is
   the policy when roundoff attacks it?
5. Why is the gravity Jacobian traceless and why is its radial eigenvalue
   positive?
6. What does the innovation measure and why is NIS the online consistency
   statistic?
7. Why does the Joseph form beat the short-form covariance update?
8. Why is NEES evaluation-only while NIS is flight-usable?
9. What happens to covariance during a GNSS outage and after recovery?
10. Why will M13C use an error-state quaternion instead of 4 additive states?
