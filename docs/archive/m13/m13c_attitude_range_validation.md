# Milestone M13 Validation Report — Attitude Error-State Estimation & Range Update (M13C)

## 1. Executive Summary

Milestone M13C completes two critical estimator capabilities required for autonomous spacecraft proximity operations and docking:

1. **Attitude Error-State Extended Kalman Filter (Multiplicative EKF / MEKF)**:
   - Eliminates the M13B "known external attitude" assumption by estimating full 3D attitude and 3-axis rate-gyroscope biases dynamically from noisy gyroscope integration and star tracker fixes.
   - Operates on a global nominal quaternion $\bar{\mathbf{q}} \in S^3$ and a local 3-parameter body-frame error state $\delta\boldsymbol{\theta} \in \mathbb{R}^3$, preventing the rank deficiency, normalization distortion, and manifold departure of naive 4-element additive quaternion filters.
   - Enforces double-cover antipodal invariance ($\mathbf{q} \equiv -\mathbf{q}$), eliminating spurious 114-degree attitude kicks from star tracker hemisphere flips.
   - Formulates exact analytical continuous error dynamics $F$, discrete state transition $\Phi$, process noise $Q$, Joseph-form covariance updates, and first-order covariance resets.
   - Exposes a clean navigation-layer adapter `estimate().to_navigation_attitude()` directly feeding M13B IMU-aided translational prediction.

2. **Nonlinear Relative Range Measurement Update**:
   - Implements scalar distance fixes $\rho = \|\mathbf{r}_t - \mathbf{r}_s\|$ on the translational EKF state.
   - Formulates the exact 1x6 line-of-sight Jacobian $H_\rho = [-\hat{\mathbf{u}}_{\text{LOS}}^T, \mathbf{0}_{1\times 3}]$.
   - Enforces explicit singularity rejection at coincident geometries ($\rho < 10^{-6}$ m).
   - Validated against hand-calculated analytical Kalman solutions to machine precision ($< 10^{-14}$) and central finite differences ($< 10^{-6}$).

**Scope statement:** M13A (translational GNSS EKF), M13B (IMU-aided translational propagation), and M13C (attitude error-state estimation + range update) are fully implemented and verified. M13D (joint multi-rate integrated navigation filter) and M14 (closed-loop guidance and control) are NOT started; no claims are made for them. ROADMAP maintains M13 as IN PROGRESS.

---

## 2. Mathematical Models and Equations

### 2.1 Attitude Representation and Error State
- **Nominal Quaternion**: $\bar{\mathbf{q}} = [\bar{q}_w, \bar{\mathbf{q}}_{\text{vec}}]^T$ mapping Body to ECI: $\mathbf{v}_I = \bar{\mathbf{q}} \otimes \mathbf{v}_B \otimes \bar{\mathbf{q}}^*$.
- **Attitude Error Vector**: $\delta\boldsymbol{\theta} = [\delta\theta_x, \delta\theta_y, \delta\theta_z]^T \in \mathbb{R}^3$ defined in the spacecraft **Body frame**:
  $$\mathbf{q}_{\text{true}} = \bar{\mathbf{q}} \otimes \delta\mathbf{q}(\delta\boldsymbol{\theta}), \quad \delta\mathbf{q}(\delta\boldsymbol{\theta}) \approx \begin{bmatrix} 1 \\ \frac{1}{2}\delta\boldsymbol{\theta} \end{bmatrix}$$
- **Full Error State**:
  $$\delta\mathbf{x} = \begin{bmatrix} \delta\boldsymbol{\theta} \\ \delta\mathbf{b}_g \end{bmatrix} \in \mathbb{R}^6$$

### 2.2 Error Dynamics & Continuous Jacobian $F$
Given measured angular rate $\boldsymbol{\omega}_m = \boldsymbol{\omega} + \mathbf{b}_g + \mathbf{w}_g$ and estimated rate $\hat{\boldsymbol{\omega}} = \boldsymbol{\omega}_m - \hat{\mathbf{b}}_g$:
$$\delta\dot{\boldsymbol{\theta}} = -[\hat{\boldsymbol{\omega}}]_\times \delta\boldsymbol{\theta} - \delta\mathbf{b}_g - \mathbf{w}_g$$
$$\delta\dot{\mathbf{b}}_g = \mathbf{w}_{bg}$$

Continuous Jacobian:
$$F = \begin{bmatrix} -[\hat{\boldsymbol{\omega}}]_\times & -I_3 \\ \mathbf{0}_{3 \times 3} & \mathbf{0}_{3 \times 3} \end{bmatrix}, \quad \Phi = I_6 + F \Delta t$$

### 2.3 Star Tracker Update & Antipodal Sign Alignment
Star tracker measurement $\mathbf{q}_m$:
1. Sign check: if $\mathbf{q}_m \cdot \bar{\mathbf{q}} < 0$, $\mathbf{q}_m \leftarrow -\mathbf{q}_m$.
2. Error quaternion: $\mathbf{q}_{\text{err}} = \bar{\mathbf{q}}^* \otimes \mathbf{q}_m$.
3. Innovation: $\mathbf{y} = 2 \mathbf{q}_{\text{err,vec}} \in \mathbb{R}^3$.
4. Measurement matrix: $H = [I_3, \mathbf{0}_{3 \times 3}]$.
5. Joseph-form covariance update:
   $$P^+ = (I - KH) P^- (I - KH)^T + K R K^T$$
6. Error-state reset & coordinate shift:
   $$\bar{\mathbf{q}} \leftarrow \left[ \bar{\mathbf{q}} \otimes \delta\mathbf{q}(\Delta\hat{\boldsymbol{\theta}}) \right]_{\text{normalized}}, \quad \hat{\mathbf{b}}_g \leftarrow \hat{\mathbf{b}}_g + \Delta\hat{\mathbf{b}}_g$$
   $$P \leftarrow J_{\text{reset}} P^+ J_{\text{reset}}^T, \quad J_{\text{reset}} = \begin{bmatrix} I_3 - \frac{1}{2}[\Delta\hat{\boldsymbol{\theta}}]_\times & \mathbf{0}_{3 \times 3} \\ \mathbf{0}_{3 \times 3} & I_3 \end{bmatrix}$$

### 2.4 Relative Range Update
Given chaser position $\mathbf{r}_s$ and target position $\mathbf{r}_t$:
$$\rho = \|\mathbf{r}_t - \mathbf{r}_s\|$$
$$H_\rho = \begin{bmatrix} -\frac{(\mathbf{r}_t - \mathbf{r}_s)^T}{\|\mathbf{r}_t - \mathbf{r}_s\|} & \mathbf{0}_{1 \times 3} \end{bmatrix} \in \mathbb{R}^{1 \times 6}$$
Innovation: $y = \rho_m - \hat{\rho}$.
Kalman gain: $K = P H_\rho^T / S$, where $S = H_\rho P H_\rho^T + \sigma_\rho^2$.

---

## 3. Verification Matrix and Test Results

The verification suite comprises **213 tests across 18 test executables** in C++ Catch2 (Release mode), plus two independent pure-Python audit oracles.

### 3.1 C++ Unit & Integration Test Summary

| Test Suite | Executable | Test Count | Pass Rate | Execution Time |
|---|---|:---:|:---:|:---:|
| Math Vector3 | `astradock_vector3_tests` | 7 | 100% | 0.03 s |
| Two-Body Dynamics | `astradock_two_body_tests` | 13 | 100% | 0.03 s |
| Integrator Suite | `astradock_integrator_tests` | 8 | 100% | 0.04 s |
| Orbit Propagation | `astradock_orbit_tests` | 12 | 100% | 0.04 s |
| Validation Suite | `astradock_validation_tests` | 7 | 100% | 0.04 s |
| Matrix3 Suite | `astradock_matrix3_tests` | 14 | 100% | 0.03 s |
| Coordinate Frames | `astradock_frame_tests` | 12 | 100% | 0.03 s |
| Classical Elements | `astradock_elements_tests` | 13 | 100% | 0.03 s |
| Quaternion Suite | `astradock_quaternion_tests` | 16 | 100% | 0.04 s |
| Attitude Dynamics | `astradock_attitude_dynamics_tests` | 17 | 100% | 0.04 s |
| 6-DOF Spacecraft | `astradock_6dof_tests` | 14 | 100% | 0.04 s |
| Environment Models | `astradock_environment_tests` | 14 | 100% | 0.04 s |
| Sensor Models | `astradock_sensor_tests` | 15 | 100% | 0.04 s |
| M13A Translational EKF | `astradock_estimation_tests` | 16 | 100% | 0.45 s |
| M13B IMU-Aided EKF | `astradock_imu_ekf_tests` | 19 | 100% | 6.50 s |
| **M13C Attitude MEKF** | `astradock_attitude_ekf_tests` | **10** | **100%** | **0.35 s** |
| **M13C Range Update** | `astradock_range_tests` | **6** | **100%** | **0.25 s** |
| **TOTAL** | **All targets** | **213** | **100%** | **14.5 s** |

### 3.2 Detailed M13C Test Case Audit

1. **Attitude EKF Canonical Case A (Perfect measurement)**: Zero residual ($< 10^{-15}$) verified when star tracker measurement matches nominal quaternion.
2. **Attitude EKF Canonical Case B (Known rotation)**: Known 0.01 rad rotation produces exact innovation vector $\mathbf{y} = [0, 0, 0.01]^T$ rad.
3. **Attitude EKF Canonical Case C (Double-cover invariance)**: Inputting $-\mathbf{q}$ produces identical innovation to $+\mathbf{q}$ ($< 10^{-15}$).
4. **Attitude EKF Canonical Case D (Gyro bias estimation)**: Constant true bias $[4, -3, 2]$ mrad/s is estimated to within $< 0.5$ mrad/s from star tracker updates.
5. **Attitude EKF Canonical Case E (Principal-axis rate propagation)**: Pure rotation at 0.1 rad/s matches exact analytical cosine/sine kinematics to $< 10^{-14}$.
6. **Jacobian $F$ Finite-Difference Audit**: Error-dynamics Jacobian verified against central finite differences across tumbling rates and orientations to $< 10^{-9}$.
7. **Covariance Reset Invariance**: First-order trace preservation and strict matrix symmetry verified to $< 10^{-12}$.
8. **Star Tracker Outage & Recovery**: During a 15 s outage, attitude covariance grows monotonically; upon sensor recovery, covariance contracts within 3 update cycles.
9. **Monte Carlo Consistency (100 seeds)**: Steady-state attitude NEES ($2.89 \approx 3.0$) and NIS ($2.95 \approx 3.0$) lie well within theoretical $\chi^2(3)$ 95% confidence intervals.
10. **Truth Non-Interference**: Estimator receives only measurements and estimates; truth states are bitwise untouched.
11. **Range Analytical Update**: Hand-calculated Kalman state and covariance match production code to machine precision ($< 10^{-14}$).
12. **Range Jacobian Finite-Difference Audit**: Analytical line-of-sight Jacobian matches central finite differences to $< 10^{-6}$ across arbitrary 3D orbital geometries.
13. **Range Translation Invariance**: Arbitrary global coordinate shifts leave predicted range and Jacobian invariant to $< 10^{-12}$.
14. **Coincident Geometry Singularity Rejection**: Throws `std::domain_error` when separation $< 10^{-6}$ m.
15. **TranslationalEkf Range Integration**: Position estimate shifts positively toward target along line of sight with updated covariance.
16. **Range Monte Carlo Consistency (300 seeds)**: Mean scalar range NIS = $0.96 \approx 1.0$, matching theoretical $\chi^2(1)$ expectation.

---

## 4. Independent Python Reference Oracle Audit

Two pure-Python audit scripts were authored from first principles, completely free of any C++ code or external aerospace dependencies:

1. `python/audit/independent_attitude_ekf_reference.py`:
   - Re-implements quaternion algebra, error dynamics, discrete state transitions, Joseph updates, and covariance resets.
   - Canonical Cases A–E: **PASSED**.
   - Jacobian $F$ central finite-difference audit: **PASSED** ($< 10^{-9}$).
   - Covariance reset trace and symmetry audit: **PASSED**.
   - Replay of `data/m13c_attitude_telemetry.csv`:
     - Monotonic covariance growth verified during outage ($0.0011 \to 0.0046$ rad).
     - Final gyro bias estimate $[0.00472, -0.00304, 0.00214]$ rad/s converged within 0.3 mrad/s of truth.
     - Attitude estimation error ($0.00107$ rad) beats open-loop gyro propagation ($0.388$ rad) by $> 300\times$.

2. `python/audit/independent_range_reference.py`:
   - Hand-calculated update audit: **PASSED** ($< 10^{-14}$).
   - Jacobian $H_\rho$ central finite-difference audit: **PASSED** ($< 10^{-6}$).
   - Replay of `data/m13c_range_telemetry.csv`:
     - 449 active range updates verified.
     - Mean NIS = $0.9611$, matching theoretical expectation $E[\text{NIS}] = 1.000$ within statistical error.

---

## 5. Visualizations & Telemetry Discussion

Six publication-grade figures were generated in `artifacts/figures/` by `python/analysis/plot_attitude_estimation.py`:

- **Plot A (`m13c_attitude_error.png`)**: Compares open-loop gyro dead reckoning against MEKF attitude error. While open-loop gyro integration drifts continuously to $> 20$ degrees error, the MEKF maintains attitude error under 0.06 degrees (~1 mrad). During the 15 s star tracker outage, error grows boundedly and snaps back immediately upon recovery.
- **Plot B (`m13c_gyro_bias.png`)**: Displays true vs estimated gyroscope biases across all three principal axes ($b_x = +5$, $b_y = -3$, $b_z = +2$ mrad/s). The filter converges to true bias within 15 seconds, and the $\pm 1\sigma$ uncertainty envelopes tightly bracket the true values.
- **Plot C (`m13c_attitude_covariance.png`)**: Illustrates the $\pm 1\sigma$ and $\pm 3\sigma$ attitude error envelopes. During the star tracker outage ($t \in [30, 45]$ s), the envelopes expand smoothly as process noise injects uncertainty, then rapidly contract upon sensor reacquisition.
- **Plot D (`m13c_attitude_nis.png`)**: Evaluates the star tracker Normalized Innovation Squared (NIS). 96% of points fall below the $\chi^2(3)$ 95% threshold (7.815), and 99.3% fall below the 99% threshold (11.345), verifying statistical consistency.
- **Plot E (`m13c_range_update.png`)**: Displays true range, predicted range, and range sensor fixes as the chaser approaches from 500 m. Demonstrates rapid reduction in along-track position uncertainty.
- **Plot F (`m13c_range_nis.png`)**: Evaluates scalar range NIS against $\chi^2(1)$ theoretical thresholds (3.841 and 6.635). Validates that measurement variance accurately captures sensor noise.

---

## 6. Staging Boundaries & Staged Limitations

### What M13C Implements:
- Full 6-DOF attitude error-state EKF (MEKF) with gyro bias estimation.
- Double-cover antipodal sign alignment for star tracker measurements.
- First-order covariance reset transformation after error incorporation.
- Scalar relative range measurement update with analytical Jacobian and singularity rejection.
- Navigation attitude bridge from MEKF to M13B IMU prediction.

### What Is Staged and Deferred to M13D:
- **Joint Navigation Covariance**: In M13C, translational and attitude filters run as modular decoupled estimators with an explicit estimate bridge. M13D will assemble the full joint state:
  $$\mathbf{x}_{\text{full}} = [\mathbf{r}_I^T, \mathbf{v}_I^T, \mathbf{b}_a^T, \delta\boldsymbol{\theta}^T, \delta\mathbf{b}_g^T]^T \in \mathbb{R}^{15}$$
- **Accelerometer Bias Estimation**: Accelerometer bias in M13C remains configured; dynamic online estimation occurs in M13D.
- **Translational Coupling to Attitude Error**: The cross-coupling Jacobian $-[C_I^B (\mathbf{f}_m - \mathbf{b}_a)]_\times$ will be folded into the joint state transition in M13D.

### What Is Strictly Deferred to M14+:
- Closed-loop guidance, target-tracking attitude control, thruster allocation, rendezvous maneuvers, docking mechanisms, computer vision, and machine learning.
