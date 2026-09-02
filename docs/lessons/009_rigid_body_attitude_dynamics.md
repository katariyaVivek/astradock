# Lesson 009: Rigid-Body Spacecraft Attitude Dynamics & Quaternion Kinematics

> *"Force causes linear acceleration and changes where a spacecraft goes; torque causes angular acceleration and gyroscopic precession, changing where it points."*

---

## 1. Physical Problem & Objective

In **M08**, we established how to mathematically represent a static 3D rotation using unit quaternions $q = [w, x, y, z]$ without coordinate singularities or gimbal lock.

In **M09**, we address the dynamic problem:
> **How does a spacecraft's orientation evolve over time under the influence of angular velocity, moments of inertia, and external torques?**

Understanding rotational motion is critical for aerospace engineers because:
1. Spacecraft must steer antennas toward ground stations, point solar panels at the Sun, and aim optical cameras at docking targets.
2. Spacecraft undergo complex multi-axis gyroscopic cross-coupling when rotating freely in Low Earth Orbit.
3. Unlike translational motion where mass $m$ is a simple scalar constant, rotational inertia is a tensorial property $\mathbf{I}$, giving rise to nonlinear gyroscopic phenomena.

```text
+----------------------------------------------------------------------------+
|                         Rotational State Propagation                       |
+----------------------------------------------------------------------------+
       Attitude q          Inertia Tensor I          External Torque tau_B
            |                     |                           |
            v                     v                           v
   +------------------------------------------------------------------+
   | Euler Rigid-Body Dynamics:                                       |
   | I * d(omega)/dt + omega x (I * omega) = tau_B                    |
   | => d(omega)/dt = I^(-1) * [ tau_B - omega x (I * omega) ]        |
   +------------------------------------------------------------------+
                                  |
                                  v
                    Updated Angular Velocity omega_B
                                  |
                                  v
   +------------------------------------------------------------------+
   | Quaternion Kinematics:                                           |
   | dq/dt = 1/2 * q (x) [0, omega_B]                                 |
   +------------------------------------------------------------------+
                                  |
                                  v
                        Updated Attitude q
```

---

## 2. Coordinate Frames & Conventions

To prevent frame confusion, every kinematic and dynamic variable in AstraDock has an explicit frame definition:

1. **Inertial Reference Frame ($\mathcal{I}$ / ECI):**
   - Non-rotating, unaccelerated frame (Earth-Centered Inertial).
   - Vectors expressed in $\mathcal{I}$ are denoted with subscript $I$ (e.g. $\mathbf{H}_I$).

2. **Spacecraft Body-Fixed Frame ($\mathcal{B}$):**
   - Fixed to the spacecraft structure, origin at the center of mass.
   - Aligned with the spacecraft's **principal axes of inertia**.
   - Vectors expressed in $\mathcal{B}$ are denoted with subscript $B$ (e.g. $\boldsymbol{\omega}_B$, $\boldsymbol{\tau}_B$, $\mathbf{H}_B$).

3. **Angular Velocity ($\boldsymbol{\omega}_B$):**
   - The angular velocity vector of the body frame $\mathcal{B}$ relative to the inertial frame $\mathcal{I}$, **expressed in body coordinates** ($\text{rad/s}$).

4. **External Torque ($\boldsymbol{\tau}_B$):**
   - The resultant moment acting on the spacecraft about its center of mass, **expressed in body coordinates** ($\text{N}\cdot\text{m}$).

5. **Attitude Quaternion ($q_{\mathcal{I}\_\mathcal{B}}$):**
   - Scalar-first unit quaternion $q = [w, x, y, z] = [w, \mathbf{v}]$ mapping vector coordinates from the body frame $\mathcal{B}$ into the inertial frame $\mathcal{I}$:
     $$\mathbf{v}_I = q_{\mathcal{I}\_\mathcal{B}} \otimes [0, \mathbf{v}_B] \otimes q_{\mathcal{I}\_\mathcal{B}}^* = C_{\mathcal{I}\_\mathcal{B}} \mathbf{v}_B$$

---

## 3. Moment of Inertia Tensor & Principal Axes

The mass distribution of a rigid spacecraft is described by its $3 \times 3$ symmetric inertia tensor $\mathbf{I}$:

$$\mathbf{I} = \begin{bmatrix}
I_{xx} & -I_{xy} & -I_{xz} \\
-I_{yx} & I_{yy} & -I_{yz} \\
-I_{zx} & -I_{zy} & I_{zz}
\end{bmatrix}$$

Because $\mathbf{I}$ is symmetric and positive-definite, there always exists an orthogonal coordinate system—the **principal axes of inertia**—in which all off-diagonal products of inertia vanish ($I_{xy} = I_{xz} = I_{yz} = 0$).

In AstraDock, the spacecraft body frame $\mathcal{B}$ is chosen to coincide with these principal axes:

$$\mathbf{I} = \begin{bmatrix}
I_{xx} & 0 & 0 \\
0 & I_{yy} & 0 \\
0 & 0 & I_{zz}
\end{bmatrix}, \quad I_{xx}, I_{yy}, I_{zz} > 0 \quad (\text{kg}\cdot\text{m}^2)$$

---

## 4. Why $\tau / I$ is Not the Whole Story: Gyroscopic Coupling

In high-school physics or 1D planar rotation, angular acceleration is simply:
$$\dot{\omega} = \frac{\tau}{I}$$

In 3D space, **this formula is incomplete and will fail catastrophically** for any rotating rigid body.

By Newton-Euler mechanics, angular momentum $\mathbf{H}$ obeys:
$$\left( \frac{d\mathbf{H}}{dt} \right)_{\mathcal{I}} = \boldsymbol{\tau}$$

Because the inertia tensor $\mathbf{I}$ is constant only when viewed from the *rotating body frame* $\mathcal{B}$, we must apply the **Transport Theorem** (Coriolis' kinematic identity for rotating frames):
$$\left( \frac{d\mathbf{H}}{dt} \right)_{\mathcal{I}} = \left( \frac{d\mathbf{H}}{dt} \right)_{\mathcal{B}} + \boldsymbol{\omega}_B \times \mathbf{H}_B$$

Since $\mathbf{H}_B = \mathbf{I}\boldsymbol{\omega}_B$ and $\left( \frac{d\mathbf{H}}{dt} \right)_{\mathcal{B}} = \mathbf{I}\dot{\boldsymbol{\omega}}_B$:

$$\mathbf{I}\dot{\boldsymbol{\omega}}_B + \boldsymbol{\omega}_B \times (\mathbf{I}\boldsymbol{\omega}_B) = \boldsymbol{\tau}_B$$

Rearranging for angular acceleration yields **Euler's Rigid-Body Equations**:

$$\dot{\boldsymbol{\omega}}_B = \mathbf{I}^{-1} \left[ \boldsymbol{\tau}_B - \boldsymbol{\omega}_B \times (\mathbf{I}\boldsymbol{\omega}_B) \right]$$

### Principal Axis Component Expansion:

$$\dot{\omega}_x = \frac{\tau_x - (I_{zz} - I_{yy})\omega_y\omega_z}{I_{xx}}$$
$$\dot{\omega}_y = \frac{\tau_y - (I_{xx} - I_{zz})\omega_z\omega_x}{I_{yy}}$$
$$\dot{\omega}_z = \frac{\tau_z - (I_{yy} - I_{xx})\omega_x\omega_y}{I_{zz}}$$

### Physical Example of Gyroscopic Coupling:
Suppose an asymmetric satellite has $I_{xx} = 10, I_{yy} = 20, I_{zz} = 30\,\text{kg}\cdot\text{m}^2$.
It is spinning with rates $\omega_x = 0.2, \omega_y = 0.3, \omega_z = 0.1\,\text{rad/s}$ in pure free space with **zero external torque** ($\boldsymbol{\tau} = \mathbf{0}$).

Evaluating the angular acceleration about the $Y$ axis:
$$\dot{\omega}_y = \frac{0 - (10 - 30)(0.1)(0.2)}{20} = \frac{-(-20)(0.02)}{20} = +0.02\,\text{rad/s}^2$$

> [!IMPORTANT]
> Even though applied torque $\tau_y = 0$, the satellite experiences an angular acceleration of $+0.02\,\text{rad/s}^2$ about $Y$ purely due to gyroscopic energy transfer between the $X$ and $Z$ axes!

---

## 5. Body-Frame vs Inertial-Frame Angular Momentum

A critical point of confusion for beginners is the distinction between angular momentum in the body frame versus the inertial frame:

### Body Frame:
$$\mathbf{H}_B = \mathbf{I}\boldsymbol{\omega}_B = \begin{bmatrix} I_{xx}\omega_x \\ I_{yy}\omega_y \\ I_{zz}\omega_z \end{bmatrix}$$
In torque-free asymmetric rotation, **the components $[H_{Bx}, H_{By}, H_{Bz}]$ are NOT constant**. They oscillate continuously as the body frame rotates.

### Inertial Frame:
$$\mathbf{H}_I = C_{\mathcal{I}\_\mathcal{B}} \mathbf{H}_B = q \otimes [0, \mathbf{H}_B] \otimes q^*$$
In torque-free motion ($\boldsymbol{\tau} = \mathbf{0}$):
$$\frac{d\mathbf{H}_I}{dt} = \mathbf{0} \implies \mathbf{H}_I(t) = \text{constant vector}$$

Furthermore, because rotation preserves vector magnitude, the Euclidean norm is invariant in both frames:
$$\|\mathbf{H}_I(t)\| = \|\mathbf{H}_B(t)\| = \text{constant}$$

---

## 6. Rotational Kinetic Energy

The rotational kinetic energy of a rigid spacecraft is:

$$E_{rot} = \frac{1}{2}\boldsymbol{\omega}_B^T \mathbf{I} \boldsymbol{\omega}_B = \frac{1}{2}\left( I_{xx}\omega_x^2 + I_{yy}\omega_y^2 + I_{zz}\omega_z^2 \right)$$

### Conservation Properties:
1. **Torque-Free Motion ($\boldsymbol{\tau} = \mathbf{0}$):**
   $$E_{rot}(t) = \text{constant}$$
   Combined with $\|\mathbf{H}_I\| = \text{constant}$, rotational motion is constrained to the intersection of two ellipsoids in angular velocity space: the **energy ellipsoid** (Poinsot ellipsoid) and the **momentum ellipsoid**.

2. **Driven Motion ($\boldsymbol{\tau} \ne \mathbf{0}$):**
   $$\frac{dE_{rot}}{dt} = \boldsymbol{\tau}_B \cdot \boldsymbol{\omega}_B$$
   External torque does mechanical work on the spacecraft, altering its kinetic energy at rate $P = \boldsymbol{\tau} \cdot \boldsymbol{\omega}$.

---

## 7. Mathematical Derivation of Quaternion Kinematics

How does the attitude quaternion $q(t)$ change as the spacecraft rotates with angular velocity $\boldsymbol{\omega}_B(t)$?

### Step 1: Incremental Rotation
Consider a small time interval $\Delta t$. Over $\Delta t$, the body frame rotates by an incremental angle $\Delta \theta = \|\boldsymbol{\omega}_B\|\Delta t$ about unit axis $\hat{\mathbf{u}} = \boldsymbol{\omega}_B / \|\boldsymbol{\omega}_B\|$.

The quaternion representing this incremental rotation of $\mathcal{B}(t+\Delta t)$ relative to $\mathcal{B}(t)$ is:
$$\Delta q = \begin{bmatrix} \cos(\Delta\theta/2) \\ \sin(\Delta\theta/2)\hat{\mathbf{u}} \end{bmatrix} \approx \begin{bmatrix} 1 \\ \frac{1}{2}\boldsymbol{\omega}_B \Delta t \end{bmatrix} = \begin{bmatrix} 1 \\ \mathbf{0} \end{bmatrix} + \frac{1}{2} \begin{bmatrix} 0 \\ \boldsymbol{\omega}_B \end{bmatrix} \Delta t$$

### Step 2: Rotation Composition
According to the M08 active frame composition convention ($q_{A\_C} = q_{A\_B} \otimes q_{B\_C}$):
$$q(t + \Delta t) = q(t) \otimes \Delta q$$

Substituting $\Delta q$:
$$q(t + \Delta t) = q(t) \otimes \left( \begin{bmatrix} 1 \\ \mathbf{0} \end{bmatrix} + \frac{1}{2}\begin{bmatrix} 0 \\ \boldsymbol{\omega}_B \end{bmatrix}\Delta t \right) = q(t) + \frac{1}{2} q(t) \otimes \begin{bmatrix} 0 \\ \boldsymbol{\omega}_B \end{bmatrix}\Delta t$$

### Step 3: Taking the Limit $\Delta t \to 0$
$$\dot{q} = \lim_{\Delta t \to 0} \frac{q(t + \Delta t) - q(t)}{\Delta t} = \frac{1}{2} q \otimes \begin{bmatrix} 0 \\ \boldsymbol{\omega}_B \end{bmatrix}$$

### Explicit Component Expansion:
Let $q = [w, x, y, z]$ and $\boldsymbol{\omega}_B = [\omega_x, \omega_y, \omega_z]$:

$$\dot{w} = -\frac{1}{2}(x\omega_x + y\omega_y + z\omega_z)$$
$$\dot{x} = \frac{1}{2}(w\omega_x + y\omega_z - z\omega_y)$$
$$\dot{y} = \frac{1}{2}(w\omega_y + z\omega_x - x\omega_z)$$
$$\dot{z} = \frac{1}{2}(w\omega_z + x\omega_y - y\omega_x)$$

> [!WARNING]
> If angular velocity was instead given in the *inertial frame* $\boldsymbol{\omega}_I$, the product order would reverse: $\dot{q} = \frac{1}{2} [0, \boldsymbol{\omega}_I] \otimes q$. Since flight sensors (gyroscopes) measure rates in the **body frame**, the body-rate form $\dot{q} = \frac{1}{2} q \otimes [0, \boldsymbol{\omega}_B]$ is universally used in onboard flight software.

---

## 8. Numerical Integration & Quaternion Normalization Policy

When integrating $\dot{q} = \frac{1}{2} q \otimes \boldsymbol{\omega}$ with discrete numerical integrators such as RK4:
- The four intermediate stages sample the continuous linear derivative space.
- Due to truncation error, finite-precision floating-point arithmetic causes the computed norm $\|q(t)\|$ to drift slightly away from $1.0$ (typically $\sim 10^{-14}$ to $10^{-8}$ over thousands of steps).

### AstraDock Policy:
1. **Never normalize intermediate RK4 substages:** Artificially normalizing $k_1, k_2, k_3$ distorts the integrator weights and destroys 4th-order convergence.
2. **Reproject onto $S^3$ at the step boundary:**
   $$q_{next} = \frac{q_{RK4}}{\|q_{RK4}\|}$$
   This bounds quaternion norm error to machine precision ($< 10^{-15}$) across indefinite simulation spans.

---

## 9. Verification & Analytical Benchmarks

M09 establishes three layers of verification:

1. **Analytical Constant-Torque Response:**
   For $\boldsymbol{\tau} = [\tau_x, 0, 0]$ on principal axis $X$:
   $$\omega_x(t) = \frac{\tau_x}{I_{xx}} t, \quad \theta_x(t) = \frac{1}{2}\frac{\tau_x}{I_{xx}} t^2, \quad q(t) = \left[ \cos\left(\frac{\theta_x}{2}\right), \sin\left(\frac{\theta_x}{2}\right), 0, 0 \right]$$

2. **Invariant Conservation in Torque-Free Tumbling:**
   For asymmetric inertia ($I_{xx} \ne I_{yy} \ne I_{zz}$) and multi-axis initial spin:
   $$\frac{|E_{rot}(t) - E_{rot}(0)|}{E_{rot}(0)} < 10^{-10}, \quad \frac{\|\mathbf{H}_I(t) - \mathbf{H}_I(0)\|}{\|\mathbf{H}_I(0)\|} < 10^{-10}$$

3. **Independent Python Oracle:**
   A standalone Python audit reference verifies C++ telemetry to sub-micro-residual accuracy ($\Delta \omega < 10^{-12}\,\text{rad/s}$, $\Delta q < 10^{-12}$).

---

## 10. Summary & Next Steps

With M09 complete, AstraDock possesses two fully verified 3-DOF dynamic halves:
- **Translational State (M01–M07):** $\mathbf{r}, \mathbf{v}$ in ECI under central gravity.
- **Rotational State (M08–M09):** $q, \boldsymbol{\omega}$ in Body/ECI under rigid-body dynamics.

In **M10**, we will combine translation and rotation into an **Integrated 6-DOF Spacecraft State**.
