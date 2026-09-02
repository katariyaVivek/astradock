# Lesson 010 — Integrated 6-DOF Spacecraft State

## 1. Physical Motivation

In space mission engineering, a spacecraft is not merely an orbital point mass traversing Keplerian ellipses, nor is it an isolated gyroscope spinning in a laboratory. A real spacecraft possesses **six physical degrees of freedom (6-DOF)**:

- **3 Translational Degrees of Freedom**: The motion of the spacecraft center of mass (COM) through three-dimensional space, governed by orbital dynamics.
- **3 Rotational Degrees of Freedom**: The angular orientation and rotation of the spacecraft body about its center of mass in $SO(3)$, governed by attitude dynamics.

Real-world spaceflight operations require continuous, synchronized simulation of both translation and rotation:
- **Solar Power Generation**: Solar panels fixed to the spacecraft body must point toward the Sun while the spacecraft traverses its orbit.
- **Payload Pointing**: Earth-observation cameras, radar antennas, and star trackers must track ground targets or inertial stars across the orbital arc.
- **Proximity Operations & Rendezvous**: Docking mechanisms (such as AstraDock's target interface) require precise relative position *and* relative orientation alignment.
- **Orbit Maneuvers**: Thrusters mounted to the spacecraft body fire in body-fixed directions ($\mathbf{F}_B$), meaning orbital velocity changes ($\Delta \mathbf{v}$) depend critically on body attitude.

Milestone M10 unifies the translational subsystem (developed in M01–M07) and the rotational subsystem (developed in M08–M09) into a clean, modular, deterministic 6-DOF simulation engine.

---

## 2. Degrees of Freedom vs Numerical State Representation

A common source of confusion in aerospace software is the distinction between **physical degrees of freedom** and the **number of numerical floating-point variables** in the simulation state vector.

```
+-----------------------------------------------------------------------------+
|               ASTRADOCK 13-COMPONENT COMPOSITE 6-DOF STATE                  |
+------------------------------------+----------------------------------------+
| Translational State (3 DOF, 6 nums)| Rotational State (3 DOF, 7 nums)       |
|   - r in R^3 (ECI position, m)     |   - q in S^3 (Attitude quaternion)     |
|   - v in R^3 (ECI velocity, m/s)   |   - omega in R^3 (Body rate, rad/s)    |
+------------------------------------+----------------------------------------+
| 6 Numerical Variables (r_x..v_z)   | 7 Numerical Variables (q_w..omega_z)   |
+------------------------------------+----------------------------------------+
| Total: 6 Physical Degrees of Freedom represented by 13 numerical values.   |
| Constraint: ||q|| = 1 removes the 1 redundant rotational degree of freedom. |
+-----------------------------------------------------------------------------+
```

1. **Translational Subsystem**:
   - 3 physical degrees of freedom $\to 6$ numerical components:
     $$\mathbf{r} = [r_x, r_y, r_z]^T \in \mathbb{R}^3 \quad (\text{m}), \quad \mathbf{v} = [v_x, v_y, v_z]^T \in \mathbb{R}^3 \quad (\text{m/s})$$
2. **Rotational Subsystem**:
   - 3 physical degrees of freedom ($SO(3)$) $\to 7$ numerical components:
     $$q = [q_w, q_x, q_y, q_z]^T \in \mathbb{S}^3 \quad (\text{dimensionless}), \quad \boldsymbol{\omega}_B = [\omega_x, \omega_y, \omega_z]^T \in \mathbb{R}^3 \quad (\text{rad/s})$$
   - The algebraic constraint $\|q\|^2 = q_w^2 + q_x^2 + q_y^2 + q_z^2 = 1$ ensures that the 4 quaternion parameters represent exactly 3 independent rotational degrees of freedom without the singularities (gimbal lock) inherent to 3-parameter Euler angles.

---

## 3. Coordinate Frames & Governing Equations

### 3.1 Coordinate Frame Conventions

To maintain strict physical clarity, AstraDock explicitly documents the reference frame of every quantity:

| Quantity | Symbol | Reference Frame | Units | Description |
| :--- | :---: | :--- | :--- | :--- |
| Position | $\mathbf{r}$ | **ECI** ($\mathcal{I}$) | $\text{m}$ | COM position relative to Earth center |
| Velocity | $\mathbf{v}$ | **ECI** ($\mathcal{I}$) | $\text{m/s}$ | Inertial velocity of COM |
| Attitude | $q_{\mathcal{I}\_\mathcal{B}}$ | $\mathcal{B} \to \mathcal{I}$ | Dimensionless | Active rotation mapping Body vectors to ECI |
| Angular Rate | $\boldsymbol{\omega}_B$ | **BODY** ($\mathcal{B}$) | $\text{rad/s}$ | Angular velocity of Body wrt ECI, in Body axes |
| External Force | $\mathbf{F}_{\mathcal{I}}$ | **ECI** ($\mathcal{I}$) | $\text{N}$ | Total non-gravitational force on COM |
| External Torque | $\boldsymbol{\tau}_B$ | **BODY** ($\mathcal{B}$) | $\text{N}\cdot\text{m}$ | Total external torque about COM |
| Inertia Tensor | $\mathbf{I}$ | **BODY** ($\mathcal{B}$) | $\text{kg}\cdot\text{m}^2$ | Principal diagonal inertia $\operatorname{diag}(I_{xx}, I_{yy}, I_{zz})$ |

### 3.2 Unified First-Order Equations of Motion

The composite state $\mathbf{X} = [\mathbf{r}^T, \mathbf{v}^T, q^T, \boldsymbol{\omega}_B^T]^T \in \mathbb{R}^{13}$ evolves under the unified ODE system:

$$\frac{d\mathbf{X}}{dt} = \mathbf{f}(t, \mathbf{X}) = \begin{bmatrix} \dot{\mathbf{r}} \\ \dot{\mathbf{v}} \\ \dot{q} \\ \dot{\boldsymbol{\omega}}_B \end{bmatrix} = \begin{bmatrix} \mathbf{v} \\ -\dfrac{\mu}{r^3}\mathbf{r} + \dfrac{\mathbf{F}_{\mathcal{I}}}{m} \\ \dfrac{1}{2} q \otimes [0, \boldsymbol{\omega}_B]^T \\ \mathbf{I}^{-1}\left( \boldsymbol{\tau}_B - \boldsymbol{\omega}_B \times (\mathbf{I}\boldsymbol{\omega}_B) \right) \end{bmatrix}$$

Where:
- $\mu = 3.986004418 \times 10^{14}\,\text{m}^3/\text{s}^2$ is Earth's gravitational parameter.
- In two-body central gravity, point-mass gravitational acceleration $\mathbf{a}_{\text{grav}} = -\frac{\mu}{r^3}\mathbf{r}$ is independent of spacecraft mass $m$.
- $\otimes$ denotes the Hamilton quaternion product: $q \otimes p$.

---

## 4. Physical Decoupling in Baseline Central Two-Body Physics

In the introductory baseline model (M10):
1. **Orbital Motion Does Not Depend on Attitude**: Central two-body gravity acts on the spacecraft center of mass as a point mass. Because Earth's gravity field is modeled as spherically symmetric, the spacecraft orientation $q$ has zero influence on $\ddot{\mathbf{r}}$.
2. **Attitude Motion Does Not Depend on Orbit**: A uniform or point-mass central gravity field exerts no net torque about the center of mass:
   $$\boldsymbol{\tau}_{\text{grav}} = \int_V \mathbf{r}' \times d\mathbf{F}_{\text{grav}} = \mathbf{0}$$
   (Torque arises only when higher-order gravity gradients $\nabla \mathbf{g}$ across a finite-sized body are considered, which is intentionally deferred to future milestones).

This clean decoupling provides a powerful verification invariant: **a spacecraft's translational orbit must remain identical whether it is spin-stabilized, tumbling at high rates, or subjected to internal/body torques.**

---

## 5. Numerical Integration & Quaternion Normalization Policy

The 13-component `SpacecraftState` implements linear vector-space addition (`operator+`) and scalar scaling (`operator*`). This allows generic numerical integrators (such as AstraDock's `numerics::rk4_step`) to operate seamlessly on composite spacecraft states.

### Normalization Policy:
1. **Intermediate RK4 Stages**: Slopes $k_1, k_2, k_3, k_4$ sample unconstrained linear state space $\mathbb{R}^{13}$. **No normalization** is applied during intermediate stages, preserving the exact 4th-order Runge-Kutta order of convergence $\mathcal{O}(\Delta t^4)$.
2. **Step Completion**: At the discrete boundary after the final stage combination:
   $$\mathbf{X}_{n+1} = \mathbf{X}_n + \frac{\Delta t}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$
   the attitude quaternion is reprojected onto the unit hypersphere:
   $$q_{n+1} \leftarrow \frac{q_{n+1}}{\|q_{n+1}\|}$$
   This completely prevents numerical norm drift over millions of integration steps while preserving machine-precision energy conservation.

---

## 6. Orientation Error Metric: Shortest Rotation Angle

To compare estimated vs reference quaternions or evaluate tracking performance without double-cover ambiguities ($q \equiv -q$), AstraDock implements the geodesic metric on $SO(3)$:

$$q_{\text{err}} = q_{\text{ref}}^* \otimes q_{\text{est}}$$

$$\theta_{\text{err}} = 2 \arccos\left(\operatorname{clamp}(|q_{\text{err}, w}|, 0.0, 1.0)\right) \in [0, \pi]\,\text{rad}$$

- **Double-Cover Invariant**: If $q_{\text{est}} = -q_{\text{ref}}$, $q_{\text{err}} = [-1, 0, 0, 0]^T \implies |w| = 1.0 \implies \theta_{\text{err}} = 0.0\,\text{rad}$.
- **Numerical Robustness**: Clamping $|w|$ to $[0.0, 1.0]$ prevents domain errors in `std::acos` from floating-point roundoff exceeding $1 + \epsilon$.

---

## 7. Educational Takeaways & Summary

1. **Unified State, Independent Physics**: Translating and rotating simultaneously is accomplished by structuring state variables as a composite hierarchy (`translational` + `rotational`), preserving the purity of each subsystem.
2. **Modular Architecture**: M10 reuses canonical M02 gravity (`dynamics::two_body_acceleration`) and M09 rigid-body kinematics/dynamics (`attitude::rotational_state_derivative`) with zero duplicated equations.
3. **Defensive Rigor**: State validation catches non-finite values and zero orbital radius while respecting the unconstrained intermediate states required by Runge-Kutta integration.
4. **Verified Baseline**: 141 tests verify exact regression equivalence against M04 orbital motion, exact regression against M09 attitude motion, fourth-order convergence, and independent Python reference parity.
