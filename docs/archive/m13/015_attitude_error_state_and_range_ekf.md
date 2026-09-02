# Lesson 015 — Attitude Error-State Estimation and Nonlinear Range Updates

## 1. The Physical Problem: Why We Cannot Dock Blind or Oriented by Guesswork

In Milestone M13B, we learned how an IMU provides inertial propagation for translational states. But M13B made an explicit operational assumption: it received the spacecraft attitude from an external estimate.

During rendezvous, proximity operations, and docking, a spacecraft cannot rely on an external attitude crutch. Consider what happens if attitude knowledge degrades:

1. **IMU Specific Force Misdirection**: Accelerometers measure specific force $\mathbf{f}_B$ in the spacecraft **body frame**. To recover inertial acceleration $\mathbf{a}_I$, this vector must be rotated into ECI:
   $$\mathbf{a}_I = C_I^B(\hat{\mathbf{q}}) (\mathbf{f}_m - \mathbf{b}_a) + \mathbf{g}(\hat{\mathbf{r}})$$
   An attitude error $\delta\boldsymbol{\theta}$ of just 1 degree (17.5 mrad) during a 1 m/s² docking burn rotates 17.5 mm/s² of thrust into the wrong inertial direction. Within 60 seconds, this accumulates into over 30 metres of along-track or cross-track navigation drift!
2. **Sensor Pointing and Line of Sight**: Rendezvous range sensors, optical cameras, and star trackers are rigidly mounted to the vehicle structure. If the vehicle attitude estimate is corrupted, the estimated line-of-sight vector to the target is erroneous, causing docking sensors to misidentify target features or lose lock completely.
3. **Range Alone Is Not Enough**: Scalar range $\rho = \|\mathbf{r}_t - \mathbf{r}_s\|$ measures distance, but distance alone provides zero angular information. To use range measurements effectively, the filter must update position along the line of sight while preserving observability and covariance geometry.

Milestone M13C solves these two fundamental aerospace estimation challenges:
- Estimating full 3D attitude and rate-gyro biases using an **Error-State Extended Kalman Filter (Multiplicative EKF / MEKF)**.
- Updating translational states with **nonlinear scalar relative-range measurements**.

---

## 2. Why a Naive 4-Element Additive EKF Fails for Quaternions

A developer new to aerospace estimation might naively write an EKF where the state vector simply includes the four quaternion components:
$$\mathbf{x}_{\text{naive}} = \begin{bmatrix} q_w \\ q_x \\ q_y \\ q_z \end{bmatrix} \in \mathbb{R}^4$$

This approach is fundamentally flawed for four mathematical and physical reasons:

### A. The Unit Norm Constraint and Covariance Singularity
Quaternions representing physical rotations in $SO(3)$ must have unit magnitude:
$$\mathbf{q}^T \mathbf{q} = q_w^2 + q_x^2 + q_y^2 + q_z^2 = 1$$
This means true rotations lie on the 3-dimensional unit sphere $S^3$ embedded in $\mathbb{R}^4$. The state space has only **3 degrees of freedom**.
If you allocate a $4 \times 4$ covariance matrix $P_q$ for a 4-element state, the covariance in the direction normal to the hypersphere ($\nabla(\mathbf{q}^T\mathbf{q}) = 2\mathbf{q}$) must be identically zero. A standard Kalman filter operating in $\mathbb{R}^4$ does not preserve this geometric constraint. Roundoff errors quickly make $P_q$ lose positive definiteness or develop spurious rank deficiency, leading to numerical divergence.

### B. Additive Updates Depart from $S^3$
In a standard EKF, the measurement update adds a correction vector:
$$\hat{\mathbf{q}}^+ = \hat{\mathbf{q}}^- + K \mathbf{y}$$
The sum of two unit quaternions is **never** a unit quaternion (unless the correction is zero). The updated state $\hat{\mathbf{q}}^+$ is pushed off the manifold $S^3$.

### C. Normalization Destroys Uncertainty Consistency
If you attempt to "fix" the norm violation by normalizing after every update:
$$\hat{\mathbf{q}}^+ \leftarrow \frac{\hat{\mathbf{q}}^+}{\|\hat{\mathbf{q}}^+\|}$$
you are performing an ad-hoc nonlinear projection that was never accounted for in the Kalman gain derivation! The computed covariance $P^+$ reflects the uncertainty of the unnormalized point, not the projected point. Over hundreds of filter cycles, the covariance ceases to represent the true state error dispersion.

### D. The Antipodal ($q \equiv -q$) Double-Cover Ambiguity
Quaternions form a 2-to-1 double cover of the rotation group $SO(3)$. For any physical orientation, the quaternions $\mathbf{q}$ and $-\mathbf{q}$ describe the **identical physical rotation matrix**:
$$C(\mathbf{q}) = C(-\mathbf{q})$$
If a star tracker outputs $-\hat{\mathbf{q}}$, a naive additive innovation produces:
$$\mathbf{y}_{\text{naive}} = (-\hat{\mathbf{q}}) - \hat{\mathbf{q}} = -2\hat{\mathbf{q}}$$
The filter sees a catastrophic innovation of magnitude 2.0 (114 degrees of error!) when the physical measurement error is actually **zero**! The filter will violently kick the spacecraft attitude estimate, causing total estimator divergence.

---

## 3. The Multiplicative Extended Kalman Filter (MEKF)

The aerospace standard for attitude estimation is the **Multiplicative Extended Kalman Filter (MEKF)**.

Instead of representing attitude error additively in 4 dimensions, the MEKF decomposes the total attitude into:
1. A **global nominal quaternion** $\bar{\mathbf{q}} \in S^3 \subset \mathbb{R}^4$ that tracks the large-angle rotational trajectory.
2. A **local error-state vector** $\delta\boldsymbol{\theta} \in \mathbb{R}^3$ expressed in the spacecraft **body frame**.

```text
               q_true = q_nom ⊗ δq(δθ)
```

where $\delta\mathbf{q}(\delta\boldsymbol{\theta})$ is the local attitude error quaternion. For small angular errors $\|\delta\boldsymbol{\theta}\| \ll 1$:
$$\delta\mathbf{q}(\delta\boldsymbol{\theta}) \approx \begin{bmatrix} 1 \\ \frac{1}{2}\delta\boldsymbol{\theta} \end{bmatrix}$$

Because $\delta\boldsymbol{\theta}$ has dimension 3, it is completely unconstrained:
- No norm constraint exists on $\delta\boldsymbol{\theta}$.
- Its covariance $P_{\theta\theta}$ is a non-singular $3 \times 3$ positive-definite matrix representing physical angular uncertainty along the spacecraft body axes (roll, pitch, yaw).
- The state space matches the physical degrees of freedom of $SO(3)$ exactly!

### Error State Definition
Augmenting the 3 attitude error angles with 3 gyroscope bias errors gives the 6-state error vector:
$$\delta\mathbf{x} = \begin{bmatrix} \delta\boldsymbol{\theta} \\ \delta\mathbf{b}_g \end{bmatrix} \in \mathbb{R}^6$$

---

## 4. Error Dynamics and Continuous Jacobian $F$

The true angular rate $\boldsymbol{\omega}$ is related to the measured rate $\boldsymbol{\omega}_m$ by:
$$\boldsymbol{\omega}_m = \boldsymbol{\omega} + \mathbf{b}_g + \mathbf{w}_g$$
where $\mathbf{b}_g$ is the true gyro bias and $\mathbf{w}_g$ is gyro white noise.

The estimated rate used to propagate the nominal quaternion is:
$$\hat{\boldsymbol{\omega}} = \boldsymbol{\omega}_m - \hat{\mathbf{b}}_g$$

Subtracting kinematics yields the continuous-time error dynamics:
$$\delta\dot{\boldsymbol{\theta}} = -[\hat{\boldsymbol{\omega}}]_\times \delta\boldsymbol{\theta} - \delta\mathbf{b}_g - \mathbf{w}_g$$
$$\delta\dot{\mathbf{b}}_g = \mathbf{w}_{bg}$$

where $[\hat{\boldsymbol{\omega}}]_\times$ is the $3 \times 3$ skew-symmetric cross-product matrix:
$$[\hat{\boldsymbol{\omega}}]_\times = \begin{bmatrix} 0 & -\hat{\omega}_z & \hat{\omega}_y \\ \hat{\omega}_z & 0 & -\hat{\omega}_x \\ -\hat{\omega}_y & \hat{\omega}_x & 0 \end{bmatrix}$$

The continuous error-dynamics Jacobian $F \in \mathbb{R}^{6 \times 6}$ is analytical:
$$F = \begin{bmatrix} -[\hat{\boldsymbol{\omega}}]_\times & -I_3 \\ \mathbf{0}_{3 \times 3} & \mathbf{0}_{3 \times 3} \end{bmatrix}$$

Notice the physical elegance:
- The upper-left block $-[\hat{\boldsymbol{\omega}}]_\times$ describes the geometric rotation of the error frame due to body tumbling.
- The upper-right block $-I_3$ directly couples gyro bias errors into attitude angle divergence ($\delta\dot{\boldsymbol{\theta}} \sim -\delta\mathbf{b}_g$).
- The bottom rows reflect the random-walk assumption for gyro bias drift.

For propagation across interval $\Delta t$, the discrete state transition matrix is:
$$\Phi = I_6 + F \Delta t$$

---

## 5. Star Tracker Update with Double-Cover Sign Alignment

When a star tracker delivers a measured orientation quaternion $\mathbf{q}_m$, the MEKF processes it as follows:

### Step 1: Antipodal Sign Alignment
We compute the 4D dot product between the measured quaternion and the nominal quaternion:
$$s = \mathbf{q}_m \cdot \bar{\mathbf{q}} = q_{m,w}\bar{q}_w + q_{m,x}\bar{q}_x + q_{m,y}\bar{q}_y + q_{m,z}\bar{q}_z$$
If $s < 0$, the measurement is in the opposite hemisphere of $S^3$. We simply negate the measurement:
$$\mathbf{q}_m \leftarrow -\mathbf{q}_m$$
This costs a single dot product and guarantees that antipodal pairs produce identical physical residuals.

### Step 2: Physical Attitude Error Quaternion
The error between the nominal quaternion and the star tracker measurement is extracted by quaternion multiplication:
$$\mathbf{q}_{\text{err}} = \bar{\mathbf{q}}^* \otimes \mathbf{q}_m$$
Writing $\mathbf{q}_{\text{err}} = [q_{\text{err},w}, \mathbf{q}_{\text{err,vec}}]^T$:
$$\mathbf{y} = 2 \mathbf{q}_{\text{err,vec}} \in \mathbb{R}^3$$
For small angles, $q_{\text{err},w} \approx 1$ and $\mathbf{q}_{\text{err,vec}} \approx \frac{1}{2}\delta\boldsymbol{\theta}$, so $\mathbf{y}$ is the physical attitude error innovation in radians!

### Step 3: Measurement Jacobian
Because the innovation $\mathbf{y}$ measures $\delta\boldsymbol{\theta}$ directly and has no direct dependence on gyro bias, the measurement Jacobian $H \in \mathbb{R}^{3 \times 6}$ is constant and exact:
$$H = \begin{bmatrix} I_3 & \mathbf{0}_{3 \times 3} \end{bmatrix}$$

---

## 6. The Joseph Update and Covariance Reset

With innovation $\mathbf{y}$, innovation covariance $S = H P^- H^T + R$, and Kalman gain $K = P^- H^T S^{-1}$, the error state correction is:
$$\Delta\hat{\mathbf{x}} = \begin{bmatrix} \Delta\hat{\boldsymbol{\theta}} \\ \Delta\hat{\mathbf{b}}_g \end{bmatrix} = K \mathbf{y}$$

The covariance is updated via the numerically stable Joseph form:
$$P^+ = (I - KH) P^- (I - KH)^T + K R K^T$$

### The Reset Step (Moving the Operating Point)
In an MEKF, the nominal quaternion is updated immediately after the measurement update by folding the estimated error angle into the nominal orientation:
$$\bar{\mathbf{q}}^+ = \left[ \bar{\mathbf{q}}^- \otimes \delta\mathbf{q}(\Delta\hat{\boldsymbol{\theta}}) \right]_{\text{normalized}}$$
$$\hat{\mathbf{b}}_g^+ = \hat{\mathbf{b}}_g^- + \Delta\hat{\mathbf{b}}_g$$

Once this correction is incorporated, the expected value of the attitude error is zeroed out:
$$\delta\hat{\boldsymbol{\theta}} \leftarrow \mathbf{0}$$

### Covariance Reset Transformation
Because the coordinate frame of $\delta\boldsymbol{\theta}$ was shifted by the rotation $\delta\mathbf{q}(\Delta\hat{\boldsymbol{\theta}})$, the covariance matrix must be rotated to account for this change of coordinates:
$$\delta\boldsymbol{\theta}_{\text{new}} \approx \left( I_3 - \frac{1}{2}[\Delta\hat{\boldsymbol{\theta}}]_\times \right) \delta\boldsymbol{\theta}_{\text{old}}$$

The complete 6x6 reset Jacobian is:
$$J_{\text{reset}} = \begin{bmatrix} I_3 - \frac{1}{2}[\Delta\hat{\boldsymbol{\theta}}]_\times & \mathbf{0}_{3 \times 3} \\ \mathbf{0}_{3 \times 3} & I_3 \end{bmatrix}$$

The reset covariance is:
$$P \leftarrow J_{\text{reset}} P^+ J_{\text{reset}}^T$$

Because $[\Delta\hat{\boldsymbol{\theta}}]_\times$ is skew-symmetric, $\text{tr}(J_{\text{reset}}) = 6$, and the trace of the covariance is invariant to first order. This reset step guarantees that the error state always operates around zero, ensuring that small-angle approximations remain valid indefinitely!

---

## 7. Nonlinear Relative Range Measurement Update

During proximity operations, the chaser spacecraft measures the scalar distance to the target vehicle:
$$\rho = \|\mathbf{r}_t - \mathbf{r}_s\| = \sqrt{(x_t - x_s)^2 + (y_t - y_s)^2 + (z_t - z_s)^2}$$

### The Analytical Jacobian $H_\rho$
The measurement is scalar ($m = 1$), updating the 6-element translational state $\mathbf{x} = [\mathbf{r}_s^T, \mathbf{v}_s^T]^T \in \mathbb{R}^6$.
Taking the gradient of $\rho$ with respect to spacecraft position $\mathbf{r}_s$:
$$\frac{\partial \rho}{\partial \mathbf{r}_s} = -\frac{(\mathbf{r}_t - \mathbf{r}_s)^T}{\|\mathbf{r}_t - \mathbf{r}_s\|} = -\hat{\mathbf{u}}_{\text{LOS}}^T$$
where $\hat{\mathbf{u}}_{\text{LOS}}$ is the unit line-of-sight vector pointing from chaser to target.
Because range does not depend on velocity, the velocity partials are zero:
$$H_\rho = \begin{bmatrix} -\hat{\mathbf{u}}_{\text{LOS}}^T & \mathbf{0}_{1 \times 3} \end{bmatrix} \in \mathbb{R}^{1 \times 6}$$

### Singularity Handling: The Coincident Geometry Defect
When $\mathbf{r}_s \to \mathbf{r}_t$, the range $\rho \to 0$, and the unit vector $\hat{\mathbf{u}}_{\text{LOS}}$ is undefined ($\frac{0}{0}$).
In AstraDock, this condition is treated as an explicit domain singularity. If separation is less than $10^{-6}$ m (1 micrometre), the estimator rejects the update loudly with `std::domain_error`, preventing numerical NaN corruption from polluting the navigation covariance.

### Observable vs Unobservable Directions
A single scalar range measurement provides information **only along the line-of-sight vector** $\hat{\mathbf{u}}_{\text{LOS}}$:
- Uncertainty in the direction of $\hat{\mathbf{u}}_{\text{LOS}}$ collapses immediately.
- Cross-line-of-sight position uncertainty is unaffected by a single range fix.
- Over time, as the orbital geometry changes (or as line-of-sight bearings rotate), successive range measurements across different geometries render all three translational axes observable.

---

## 8. Summary of Milestone M13C Deliverables

1. `cpp/estimation/attitude_ekf.hpp`: Production MEKF with nominal quaternion propagation, continuous Jacobian $F$, discrete state transition $\Phi$, process noise $Q$, star-tracker update with double-cover alignment, Joseph form, and covariance reset.
2. `cpp/estimation/range_update.hpp`: Nonlinear range prediction, analytical 1x6 Jacobian, singularity detection, scalar Joseph update.
3. `cpp/estimation/translational_ekf.hpp`: Integrated `update_range` method and diagnostics.
4. `tests/cpp/test_attitude_ekf.cpp`: Canonical cases A–E, Jacobian finite differences, covariance reset audit, dropout verification, 100-seed Monte Carlo NEES/NIS consistency, truth non-interference.
5. `tests/cpp/test_range_update.cpp`: Analytical hand-calculated test, Jacobian finite differences, translation invariance, singularity rejection, filter integration, Monte Carlo NIS consistency.
6. `tools/attitude_ekf_demo.cpp`: Produces authoritative CSV telemetry for attitude and range scenarios.
7. `python/audit/independent_attitude_ekf_reference.py` & `python/audit/independent_range_reference.py`: Independent Python oracles verifying all mathematical steps from scratch.
8. `python/analysis/plot_attitude_estimation.py`: Publication-quality visualizations (Plots A–F) in `artifacts/figures/`.
