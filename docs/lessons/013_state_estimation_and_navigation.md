# Lesson 013 — Spacecraft State Estimation & Integrated Navigation

## 1. Physical Problem & Motivation

A spacecraft cannot navigate by dead reckoning alone, nor can it trust raw sensor readings directly. Sensors are noisy, biased, rate-limited, asynchronous, and subject to environmental dropouts:
- **GNSS receivers** yield low-frequency ($1\,\text{Hz}$) absolute position and velocity fixes in ECI, but are noisy ($\sigma_r \approx 5\text{--}10\,\text{m}$, $\sigma_v \approx 0.05\,\text{m/s}$) and subject to antenna occlusion or orbit geometry dropouts.
- **Inertial Measurement Units (IMUs)** sample at high frequency ($100\,\text{Hz}$), providing specific force $\mathbf{f}_B = \mathbf{a}_B - \mathbf{g}_B$ and angular velocity $\boldsymbol{\omega}_B$, but carry sensor biases ($\mathbf{b}_a, \mathbf{b}_g$) and white noise that cause unconstrained quadratic and cubic drift over time.
- **Star trackers** deliver high-accuracy attitude fixes ($10\,\text{Hz}$, sub-milliradian), but can be blinded by the Sun, Earth albedo, or plume contamination.
- **Relative range sensors** provide high-accuracy line-of-sight distance fixes to cooperative targets, but measure only one spatial dimension.

The central problem of spacecraft navigation is **state estimation**: inferring the true continuous physical state of the spacecraft by fusing heterogeneous, multi-rate, noisy measurements with a physical dynamics model.

---

## 2. Why Separate Filters Fail: The Necessity of Cross-Covariance

A naive architectural approach divides estimation into two decoupled filters:
1. A **translational filter** estimating position and velocity $[\mathbf{r}_I, \mathbf{v}_I]$.
2. An **attitude filter** estimating orientation $\mathbf{q}$.

This decoupling fundamentally violates orbital physics. In an inertial frame, the specific force measured in the spacecraft body frame ($\mathbf{f}_B = \mathbf{f}_m - \mathbf{b}_a$) is transformed into inertial acceleration via the estimated attitude quaternion $\mathbf{q}$:

$$\mathbf{a}_I = C_I^B(\mathbf{q}) (\mathbf{f}_m - \mathbf{b}_a) + \mathbf{g}(\mathbf{r}_I)$$

When the spacecraft undergoes non-gravitational acceleration (e.g. atmospheric drag, solar radiation pressure, RCS attitude control pulses, or orbit transfer thrust), an attitude error $\delta\boldsymbol{\theta}$ rotates the large specific force vector into an erroneous inertial direction:

$$\delta\mathbf{a}_I \approx - C_I^B(\mathbf{q}) [\mathbf{f}_B]_\times \delta\boldsymbol{\theta}$$

This generates a direct physical coupling:
$$\delta\boldsymbol{\theta} \implies \delta\mathbf{a}_I \implies \delta\mathbf{v}_I \implies \delta\mathbf{r}_I$$

Furthermore, an uncalibrated accelerometer bias $\delta\mathbf{b}_a$ creates a continuous translational acceleration error:
$$\delta\mathbf{a}_I \approx - C_I^B(\mathbf{q}) \delta\mathbf{b}_a$$

If attitude and translation are handled by separate filters with block-diagonal covariances:
1. The translational filter understates its covariance during maneuvers because it assumes the attitude used to project specific force is perfectly known with zero variance.
2. High-accuracy GNSS updates cannot refine attitude estimates or calibrate accelerometer biases.
3. Star tracker dropouts cannot correctly inflate position and velocity uncertainty envelopes.

**The $15 \times 15$ joint error covariance matrix is the mathematical conduit that maintains these physical correlations.**

---

## 3. State Representation: Nominal vs Error State

To respect the geometric constraints of rotation on the 3-sphere $S^3 \subset \mathbb{R}^4$ while maintaining an unconstrained linear Gaussian estimation error, AstraDock employs a **Multiplicative Extended Kalman Filter (MEKF)** formulation.

### The 16-Component Nominal State
The nominal state vector $\mathbf{x}_{\text{nom}}$ stores the best current physical estimates:
$$\mathbf{x}_{\text{nom}} = \begin{bmatrix} \mathbf{r}_I \\ \mathbf{v}_I \\ \mathbf{q} \\ \mathbf{b}_a \\ \mathbf{b}_g \end{bmatrix} \in \mathbb{R}^3 \times \mathbb{R}^3 \times S^3 \times \mathbb{R}^3 \times \mathbb{R}^3$$
- $\mathbf{r}_I \in \mathbb{R}^3$: Spacecraft center-of-mass position in ECI ($\text{m}$).
- $\mathbf{v}_I \in \mathbb{R}^3$: Spacecraft inertial velocity in ECI ($\text{m/s}$).
- $\mathbf{q} \in S^3$: Unit quaternion representing orientation of the spacecraft Body frame relative to ECI ($[w, x, y, z]$, scalar-first, $\|\mathbf{q}\| = 1$).
- $\mathbf{b}_a \in \mathbb{R}^3$: Accelerometer bias in Spacecraft **Body frame** ($\text{m/s}^2$).
- $\mathbf{b}_g \in \mathbb{R}^3$: Rate-gyroscope bias in Spacecraft **Body frame** ($\text{rad/s}$).

### The 15-Component Error State
The estimation error $\delta\mathbf{x} \in \mathbb{R}^{15}$ is strictly unconstrained:
$$\delta\mathbf{x} = \begin{bmatrix} \delta\mathbf{r} \\ \delta\mathbf{v} \\ \delta\boldsymbol{\theta} \\ \delta\mathbf{b}_a \\ \delta\mathbf{b}_g \end{bmatrix} \in \mathbb{R}^{15}$$
where:
$$\mathbf{r}_{\text{true}} = \hat{\mathbf{r}} + \delta\mathbf{r}$$
$$\mathbf{v}_{\text{true}} = \hat{\mathbf{v}} + \delta\mathbf{v}$$
$$\mathbf{q}_{\text{true}} = \hat{\mathbf{q}} \otimes \delta\mathbf{q}(\delta\boldsymbol{\theta}) \approx \hat{\mathbf{q}} \otimes \begin{bmatrix} 1 \\ \frac{1}{2}\delta\boldsymbol{\theta} \end{bmatrix}$$
$$\mathbf{b}_{a,\text{true}} = \hat{\mathbf{b}}_a + \delta\mathbf{b}_a$$
$$\mathbf{b}_{g,\text{true}} = \hat{\mathbf{b}}_g + \delta\mathbf{b}_g$$

### State Index Ordering (0-indexed)
```text
Indices  0..2  : delta_r_x,   delta_r_y,   delta_r_z    (Position error in ECI, m)
Indices  3..5  : delta_v_x,   delta_v_y,   delta_v_z    (Velocity error in ECI, m/s)
Indices  6..8  : delta_th_x,  delta_th_y,  delta_th_z   (Attitude error vector in BODY, rad)
Indices  9..11 : delta_ba_x,  delta_ba_y,  delta_ba_z   (Accel bias error in BODY, m/s^2)
Indices 12..14 : delta_bg_x,  delta_bg_y,  delta_bg_z   (Gyro bias error in BODY, rad/s)
```

The error-state covariance is a dense $15 \times 15$ symmetric positive-definite matrix:
$$P = E\left[\delta\mathbf{x} \delta\mathbf{x}^T\right] \in \mathbb{R}^{15 \times 15}$$

---

## 4. Continuous Error Dynamics & Analytical $15 \times 15$ Jacobian

The continuous-time linearized error dynamics take the form:
$$\delta\dot{\mathbf{x}}(t) = F(t) \delta\mathbf{x}(t) + W(t) \mathbf{w}(t)$$

### Block Structure of $F \in \mathbb{R}^{15 \times 15}$
$$
F = \begin{bmatrix}
\mathbf{0}_{3\times 3} & I_3 & \mathbf{0}_{3\times 3} & \mathbf{0}_{3\times 3} & \mathbf{0}_{3\times 3} \\
G(\hat{\mathbf{r}}) & \mathbf{0}_{3\times 3} & -C_I^B(\hat{\mathbf{q}}) [\hat{\mathbf{f}}_B]_\times & -C_I^B(\hat{\mathbf{q}}) & \mathbf{0}_{3\times 3} \\
\mathbf{0}_{3\times 3} & \mathbf{0}_{3\times 3} & -[\hat{\boldsymbol{\omega}}_B]_\times & \mathbf{0}_{3\times 3} & -I_3 \\
\mathbf{0}_{3\times 3} & \mathbf{0}_{3\times 3} & \mathbf{0}_{3\times 3} & \mathbf{0}_{3\times 3} & \mathbf{0}_{3\times 3} \\
\mathbf{0}_{3\times 3} & \mathbf{0}_{3\times 3} & \mathbf{0}_{3\times 3} & \mathbf{0}_{3\times 3} & \mathbf{0}_{3\times 3}
\end{bmatrix}
$$

where:
1. **Gravity Gradient Tensor**:
   $$G(\hat{\mathbf{r}}) = \left.\frac{\partial \mathbf{g}}{\partial \mathbf{r}}\right|_{\hat{\mathbf{r}}} = -\frac{\mu}{r^3}\left(I_3 - 3\frac{\hat{\mathbf{r}}\hat{\mathbf{r}}^T}{r^2}\right)$$
2. **Specific Force Cross-Coupling**:
   $$\frac{\partial \delta\dot{\mathbf{v}}}{\partial \delta\boldsymbol{\theta}} = - C_I^B(\hat{\mathbf{q}}) [\hat{\mathbf{f}}_B]_\times, \quad \hat{\mathbf{f}}_B = \mathbf{f}_m - \hat{\mathbf{b}}_a$$
3. **Accelerometer Bias Coupling**:
   $$\frac{\partial \delta\dot{\mathbf{v}}}{\partial \delta\mathbf{b}_a} = - C_I^B(\hat{\mathbf{q}})$$
4. **Attitude Rate Error Dynamics**:
   $$\frac{\partial \delta\dot{\boldsymbol{\theta}}}{\partial \delta\boldsymbol{\theta}} = -[\hat{\boldsymbol{\omega}}_B]_\times, \quad \hat{\boldsymbol{\omega}}_B = \boldsymbol{\omega}_m - \hat{\mathbf{b}}_g$$
5. **Gyro Bias Coupling**:
   $$\frac{\partial \delta\dot{\boldsymbol{\theta}}}{\partial \delta\mathbf{b}_g} = -I_3$$

---

## 5. Discrete Propagation & Covariance Prediction

Between measurement arrivals, the high-rate IMU ($100\,\text{Hz}$) drives time updates across interval $\Delta t = t_k - t_{k-1}$ computed from actual timestamps:

1. **Nominal State Kinematics**:
   $$\mathbf{a}_I = C_I^B(\hat{\mathbf{q}})(\mathbf{f}_m - \hat{\mathbf{b}}_a) + \mathbf{g}(\hat{\mathbf{r}})$$
   $$\hat{\mathbf{r}}_{k} = \hat{\mathbf{r}}_{k-1} + \hat{\mathbf{v}}_{k-1}\Delta t + \frac{1}{2}\mathbf{a}_I \Delta t^2$$
   $$\hat{\mathbf{v}}_{k} = \hat{\mathbf{v}}_{k-1} + \mathbf{a}_I \Delta t$$
   $$\hat{\mathbf{q}}_{k} = \left[\hat{\mathbf{q}}_{k-1} \otimes \Delta\mathbf{q}\left((\boldsymbol{\omega}_m - \hat{\mathbf{b}}_g)\Delta t\right)\right]_{\text{normalized}}$$

2. **Discrete Transition**:
   $$\Phi_k = I_{15} + F_k \Delta t$$

3. **Discrete Process Noise Covariance $Q \in \mathbb{R}^{15 \times 15}$**:
   Driven by sensor white noise and bias random walk spectral densities:
   - Position-velocity sub-blocks: $Q_{rr} = \frac{1}{3}\sigma_a^2 \Delta t^3 I_3$, $Q_{rv} = \frac{1}{2}\sigma_a^2 \Delta t^2 I_3$, $Q_{vv} = \sigma_a^2 \Delta t I_3$.
   - Attitude sub-block: $Q_{\theta\theta} = \sigma_g^2 \Delta t I_3$.
   - Bias random walk sub-blocks: $Q_{ba} = \sigma_{ba}^2 \Delta t I_3$, $Q_{bg} = \sigma_{bg}^2 \Delta t I_3$.

4. **Covariance Propagation**:
   $$P_k^- = \Phi_k P_{k-1}^+ \Phi_k^T + Q_k$$
   symmetrized via $P_k^- \leftarrow \frac{1}{2}\left(P_k^- + (P_k^-)^T\right)$.

---

## 6. Asynchronous Measurement Updates & Covariance Resets

The estimator processes measurements sequentially as they arrive on the multi-rate timeline:

### 1. GNSS Measurement Update ($1\,\text{Hz}$)
- **Measurement**: $\mathbf{z}_{\text{GNSS}} = [\mathbf{r}_m^T, \mathbf{v}_m^T]^T \in \mathbb{R}^6$.
- **Jacobian**:
  $$H_{\text{GNSS}} = \begin{bmatrix} I_3 & \mathbf{0}_{3\times 3} & \mathbf{0}_{3\times 3} & \mathbf{0}_{3\times 3} & \mathbf{0}_{3\times 3} \\ \mathbf{0}_{3\times 3} & I_3 & \mathbf{0}_{3\times 3} & \mathbf{0}_{3\times 3} & \mathbf{0}_{3\times 3} \end{bmatrix} \in \mathbb{R}^{6 \times 15}$$

### 2. Star Tracker Measurement Update ($10\,\text{Hz}$)
- **Measurement**: $\mathbf{q}_m \in S^3$.
- **Antipodal Double-Cover Sign Alignment**:
  $$\text{If } \mathbf{q}_m \cdot \hat{\mathbf{q}} < 0 \implies \mathbf{q}_m \leftarrow -\mathbf{q}_m$$
  This prevents catastrophic $114^\circ$ spurious innovations caused by the $SO(3)$ double cover.
- **Physical Error Vector Innovation**:
  $$\mathbf{q}_{\text{err}} = \hat{\mathbf{q}}^* \otimes \mathbf{q}_m, \quad \mathbf{y}_{\text{ST}} = 2 \mathbf{q}_{\text{err,vec}} \in \mathbb{R}^3$$
- **Jacobian**:
  $$H_{\text{ST}} = \begin{bmatrix} \mathbf{0}_{3\times 3} & \mathbf{0}_{3\times 3} & I_3 & \mathbf{0}_{3\times 3} & \mathbf{0}_{3\times 3} \end{bmatrix} \in \mathbb{R}^{3 \times 15}$$

### 3. Relative Range Measurement Update ($10\,\text{Hz}$)
- **Measurement**: $\rho_m \in \mathbb{R}^+$.
- **Line-of-Sight Unit Vector**:
  $$\hat{\rho} = \|\mathbf{r}_{\text{target}} - \hat{\mathbf{r}}\|, \quad \hat{\mathbf{u}}_{\text{LOS}} = \frac{\mathbf{r}_{\text{target}} - \hat{\mathbf{r}}}{\hat{\rho}}$$
  If $\hat{\rho} < 10^{-6}\,\text{m}$, reject geometry as singular (`domain_error`).
- **Jacobian**:
  $$H_{\text{range}} = \begin{bmatrix} -\hat{\mathbf{u}}_{\text{LOS}}^T & \mathbf{0}_{1\times 3} & \mathbf{0}_{1\times 3} & \mathbf{0}_{1\times 3} & \mathbf{0}_{1\times 3} \end{bmatrix} \in \mathbb{R}^{1 \times 15}$$

### Joseph-Form Covariance Update
For any measurement model $(H, R, \mathbf{y})$:
$$S = H P^- H^T + R$$
$$K = P^- H^T S^{-1} \quad (\text{via Cholesky factorization of } S)$$
$$\Delta\hat{\mathbf{x}} = K \mathbf{y} \in \mathbb{R}^{15}$$
$$P^+ = (I_{15} - K H) P^- (I_{15} - K H)^T + K R K^T$$

### State Correction and Covariance Reset
The correction $\Delta\hat{\mathbf{x}}$ is folded into the nominal states:
$$\hat{\mathbf{r}} \leftarrow \hat{\mathbf{r}} + \Delta\hat{\mathbf{x}}_{0..2}, \quad \hat{\mathbf{v}} \leftarrow \hat{\mathbf{v}} + \Delta\hat{\mathbf{x}}_{3..5}$$
$$\hat{\mathbf{q}} \leftarrow \left[\hat{\mathbf{q}} \otimes \delta\mathbf{q}(\Delta\hat{\boldsymbol{\theta}})\right]_{\text{normalized}}, \quad \Delta\hat{\boldsymbol{\theta}} = \Delta\hat{\mathbf{x}}_{6..8}$$
$$\hat{\mathbf{b}}_a \leftarrow \hat{\mathbf{b}}_a + \Delta\hat{\mathbf{x}}_{9..11}, \quad \hat{\mathbf{b}}_g \leftarrow \hat{\mathbf{b}}_g + \Delta\hat{\mathbf{x}}_{12..14}$$

Because the attitude correction $\Delta\hat{\boldsymbol{\theta}}$ resets the error state to zero, the coordinate frame of the attitude error rotates, requiring a first-order covariance reset transformation:
$$J_{\text{reset}} = \text{diag}\left(I_3, I_3, I_3 - \frac{1}{2}[\Delta\hat{\boldsymbol{\theta}}]_\times, I_3, I_3\right) \in \mathbb{R}^{15 \times 15}$$
$$P \leftarrow J_{\text{reset}} P^+ J_{\text{reset}}^T$$

---

## 7. Consistency Diagnostics: NEES and NIS

Estimator correctness is evaluated using two rigorous statistical metrics:

1. **Normalized Innovation Squared (NIS)**:
   $$\text{NIS} = \mathbf{y}^T S^{-1} \mathbf{y} \sim \chi^2(m)$$
   Evaluated online per sensor:
   - GNSS: $m = 6$, $E[\text{NIS}] = 6.0$.
   - Star Tracker: $m = 3$, $E[\text{NIS}] = 3.0$.
   - Range: $m = 1$, $E[\text{NIS}] = 1.0$.

2. **Normalized Estimation Error Squared (NEES)**:
   $$\text{NEES}_{15} = \delta\mathbf{x}^T P^{-1} \delta\mathbf{x} \sim \chi^2(15)$$
   Evaluated in validation harnesses where truth is known. The theoretical expectation is $E[\text{NEES}_{15}] = 15.0$ with two-tailed 95% confidence bounds $[6.26, 27.49]$.

---

## 8. Summary of Milestones M13A–M13D

| Sub-Stage | Focus Area | Key Architectural Contribution |
|---|---|---|
| **M13A** | Translational GNSS EKF | Linearized gravity gradient $G(r)$, Cholesky solves, Joseph update, NIS/NEES diagnostics. |
| **M13B** | IMU-Aided Propagation | Translating specific force to inertial acceleration, gravity compensation, timestamp discipline. |
| **M13C** | Attitude Error-State & Range | MEKF on $S^3$, double-cover sign alignment, covariance reset, scalar nonlinear range update. |
| **M13D** | Integrated Multi-Rate Navigation | Unified 15-state EKF, cross-covariance coupling ($\delta\boldsymbol{\theta} \to \delta\mathbf{v}$, $\delta\mathbf{b}_a \to \delta\mathbf{v}$), online bias estimation. |
