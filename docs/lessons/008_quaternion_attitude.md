# Lesson 008: Spacecraft Attitude Representation and Quaternion Mathematics

> *"Orbital mechanics tells you where a spacecraft is in space. Attitude mathematics tells you which way it is pointing."*

---

## 1. Spacecraft Attitude vs Orbital State

In spaceflight engineering, a spacecraft's full physical state is fundamentally decoupled into two independent domains:

1. **Translational (Orbital) State:**
   - Position and velocity vectors $\mathbf{r}, \mathbf{v} \in \mathbb{R}^3$ in an Earth-Centered Inertial (ECI) frame.
   - Dictated by gravitational acceleration and orbital mechanics (M01–M07).

2. **Rotational (Attitude) State:**
   - The instantaneous spatial orientation of the spacecraft body-fixed frame ($\mathcal{B}$) relative to a reference frame (such as ECI $\mathcal{I}$ or local orbital LVLH $\mathcal{L}$).
   - Independent of position: a satellite can be at a specific point along its orbit while pointing its solar panels at the Sun, its antenna at Earth, or tumbling freely.

```text
               +-------------------------------------------+
               |             Spacecraft State              |
               +-------------------------------------------+
                        |                         |
                        v                         v
               +-----------------+       +-----------------+
               |  Orbital State  |       | Attitude State  |
               |   r, v in ECI   |       | Orientation of  |
               | (Translational) |       |  Body vs Ref    |
               +-----------------+       +-----------------+
```

---

## 2. Comparing Attitude Representations

Aerospace engineers use four primary mathematical representations for 3D rotations, each with distinct trade-offs:

| Representation | Parameters | Constraints | Advantages | Limitations |
| :--- | :---: | :---: | :--- | :--- |
| **Euler Angles** (e.g. Yaw-Pitch-Roll) | 3 | None | Highly intuitive; matches human perception of aircraft/spacecraft rotations. | **Gimbal lock:** coordinate singularity at pitch $\theta = \pm 90^\circ$; non-unique representations; computationally expensive trigonometric evaluations. |
| **Direction Cosine Matrix (DCM)** | 9 ($3 \times 3$) | 6 constraints ($C^T C = I, \det(C) = +1$) | Singularity-free; straightforward vector transformation via matrix multiplication; easy frame chaining. | Redundant (9 parameters for 3 DOF); prone to numerical drift requiring frequent matrix orthogonalization; expensive to propagate. |
| **Axis-Angle** $(\hat{\mathbf{u}}, \theta)$ | 4 ($3 + 1$) | $\|\hat{\mathbf{u}}\| = 1$ | Direct physical interpretation based on Euler's rotation theorem (single rotation about one axis). | Cannot be directly composed or multiplied without converting to another representation; singular at $\theta = 0$. |
| **Unit Quaternion** $q = [w, x, y, z]$ | 4 | 1 constraint ($\|q\| = 1$) | **Singularity-free;** minimal redundant parameter (only 1 constraint); simple normalization ($\mathbf{q} / \|\mathbf{q}\|$); fast algebraic composition via Hamilton product (no trig functions); bilinear kinematics. | Less directly intuitive; 2-to-1 double cover ($+q$ and $-q$ represent the same physical rotation). |

---

## 3. Gimbal Lock: Why Euler Angles Fail for Global Spacecraft GNC

Euler angles represent a 3D rotation as three successive rotations about moving axes (e.g., Z-Y-X yaw, pitch, roll). 

When the pitch angle $\theta$ approaches $+90^\circ$ or $-90^\circ$:
- The first rotation axis ($Z$) and the third rotation axis ($X''$) become parallel and align along the same spatial line.
- One degree of freedom is lost: changing yaw or roll produces the exact same physical rotation.
- In the inverse transformation from a rotation matrix back to Euler angles, $\psi$ and $\phi$ appear only as the sum or difference $(\phi \pm \psi)$, rendering individual angles indeterminate.

> [!WARNING]
> Gimbal lock is a **mathematical artifact of 3-parameter coordinate systems** (Euler's rotation theorem proves no 3-parameter global parameterization of $SO(3)$ can be both global and singularity-free). It is NOT a physical locking of the spacecraft. To avoid numerical crashes, flight software uses unit quaternions as the authoritative attitude state.

---

## 4. Quaternion Anatomy & Conventions

A quaternion is a four-dimensional hypercomplex number discovered by William Rowan Hamilton in 1843:
$$q = w + x\mathbf{i} + y\mathbf{j} + z\mathbf{k}$$
where the imaginary basis units satisfy Hamilton's fundamental rules:
$$\mathbf{i}^2 = \mathbf{j}^2 = \mathbf{k}^2 = \mathbf{i}\mathbf{j}\mathbf{k} = -1$$
$$\mathbf{i}\mathbf{j} = \mathbf{k} = -\mathbf{j}\mathbf{i}, \quad \mathbf{j}\mathbf{k} = \mathbf{i} = -\mathbf{k}\mathbf{j}, \quad \mathbf{k}\mathbf{i} = \mathbf{j} = -\mathbf{i}\mathbf{k}$$

### AstraDock Scalar-First Notation
AstraDock adopts the standard aerospace **scalar-first** layout:
$$q = \begin{bmatrix} w \\ x \\ y \\ z \end{bmatrix} = \begin{bmatrix} w \\ \mathbf{v} \end{bmatrix}$$
where:
- $w \in \mathbb{R}$ is the scalar (real) component.
- $\mathbf{v} = [x, y, z]^T \in \mathbb{R}^3$ is the vector (imaginary) component.

---

## 5. Foundational Quaternion Algebra

### 5.1 Euclidean Norm
$$\|q\| = \sqrt{w^2 + x^2 + y^2 + z^2}$$
For physical attitude representations, we restrict ourselves to **unit quaternions** ($\|q\| = 1$).

### 5.2 Normalization
$$\text{normalize}(q) = \frac{q}{\|q\|}$$
*Defensive Policy:* If $\|q\| \le 10^{-14}$ or non-finite, AstraDock throws `std::domain_error`.

### 5.3 Conjugate
The conjugate flips the signs of the imaginary vector components:
$$q^* = \begin{bmatrix} w \\ -x \\ -y \\ -z \end{bmatrix} = \begin{bmatrix} w \\ -\mathbf{v} \end{bmatrix}$$
Physically, $q^*$ represents the **reverse rotation** (rotation by $-\theta$ about $\hat{\mathbf{u}}$).

### 5.4 Inverse
$$q^{-1} = \frac{q^*}{\|q\|^2}$$
For any unit quaternion ($\|q\| = 1$), the inverse is identically the conjugate:
$$q^{-1} = q^*$$

---

## 6. The Hamilton Product (Quaternion Multiplication)

Given $q_1 = [w_1, \mathbf{v}_1]$ and $q_2 = [w_2, \mathbf{v}_2]$, their Hamilton product $q_1 \otimes q_2$ is:

$$q_1 \otimes q_2 = \begin{bmatrix} w_1 w_2 - \mathbf{v}_1 \cdot \mathbf{v}_2 \\ w_1 \mathbf{v}_2 + w_2 \mathbf{v}_1 + \mathbf{v}_1 \times \mathbf{v}_2 \end{bmatrix}$$

In explicit scalar-first components:
$$q_1 \otimes q_2 = \begin{bmatrix}
w_1 w_2 - x_1 x_2 - y_1 y_2 - z_1 z_2 \\
w_1 x_2 + x_1 w_2 + y_1 z_2 - z_1 y_2 \\
w_1 y_2 - x_1 z_2 + y_1 w_2 + z_1 x_2 \\
w_1 z_2 + x_1 y_2 - y_1 x_2 + z_1 w_2
\end{bmatrix}$$

### Key Properties:
1. **Non-Commutative:** $q_1 \otimes q_2 \ne q_2 \otimes q_1$ in general (because 3D cross products are anticommutative).
2. **Associative:** $(q_1 \otimes q_2) \otimes q_3 = q_1 \otimes (q_2 \otimes q_3)$.
3. **Identity Element:** $q_I = [1, 0, 0, 0]$ satisfies $q \otimes q_I = q_I \otimes q = q$.
4. **Norm Multiplication:** $\|q_1 \otimes q_2\| = \|q_1\| \cdot \|q_2\|$.

---

## 7. Vector Rotation: Active Transformation Formula

To rotate a 3D vector $\mathbf{v} \in \mathbb{R}^3$ by a unit quaternion $q$, we embed $\mathbf{v}$ into a pure quaternion $p = [0, \mathbf{v}]$ and evaluate the conjugation:
$$p' = q \otimes p \otimes q^*$$

Expanding the quaternion product yields the efficient vector rotation formula (Rodrigues' formula in quaternion form):
$$\mathbf{v}' = \mathbf{v} + 2w(\mathbf{u} \times \mathbf{v}) + 2(\mathbf{u} \times (\mathbf{u} \times \mathbf{v}))$$
where $\mathbf{u} = [x, y, z]^T$.

This formula requires only 15 multiplications and 15 additions, without evaluating any matrix elements.

---

## 8. Worked Canonical Example: 90° Rotation About +Z

Let us rotate the vector $\mathbf{v} = [1, 0, 0]^T$ by $\theta = 90^\circ$ ($\pi/2$ rad) around the $+Z$ axis ($\hat{\mathbf{u}} = [0, 0, 1]^T$).

### Step 1: Axis-Angle Construction
$$\text{Half angle: } \theta/2 = 45^\circ = \pi/4 \implies \cos(\pi/4) = \frac{1}{\sqrt{2}}, \quad \sin(\pi/4) = \frac{1}{\sqrt{2}}$$
$$q = \begin{bmatrix} \cos(\pi/4) \\ 0 \\ 0 \\ \sin(\pi/4) \end{bmatrix} = \begin{bmatrix} 1/\sqrt{2} \\ 0 \\ 0 \\ 1/\sqrt{2} \end{bmatrix} \approx \begin{bmatrix} 0.707107 \\ 0 \\ 0 \\ 0.707107 \end{bmatrix}$$

### Step 2: Vector Rotation Evaluation
- $\mathbf{u} = [0, 0, 1/\sqrt{2}]^T$, $\mathbf{v} = [1, 0, 0]^T$, $w = 1/\sqrt{2}$.
- $\mathbf{u} \times \mathbf{v} = \begin{bmatrix} 0 \\ 0 \\ 1/\sqrt{2} \end{bmatrix} \times \begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix} = \begin{bmatrix} 0 \\ 1/\sqrt{2} \\ 0 \end{bmatrix}$.
- $2w(\mathbf{u} \times \mathbf{v}) = 2(1/\sqrt{2})\begin{bmatrix} 0 \\ 1/\sqrt{2} \\ 0 \end{bmatrix} = \begin{bmatrix} 0 \\ 1 \\ 0 \end{bmatrix}$.
- $\mathbf{u} \times (\mathbf{u} \times \mathbf{v}) = \begin{bmatrix} 0 \\ 0 \\ 1/\sqrt{2} \end{bmatrix} \times \begin{bmatrix} 0 \\ 1/\sqrt{2} \\ 0 \end{bmatrix} = \begin{bmatrix} -1/2 \\ 0 \\ 0 \end{bmatrix}$.
- $2(\mathbf{u} \times (\mathbf{u} \times \mathbf{v})) = \begin{bmatrix} -1 \\ 0 \\ 0 \end{bmatrix}$.
- $\mathbf{v}' = \mathbf{v} + 2w(\mathbf{u} \times \mathbf{v}) + 2(\mathbf{u} \times (\mathbf{u} \times \mathbf{v})) = \begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix} + \begin{bmatrix} 0 \\ 1 \\ 0 \end{bmatrix} + \begin{bmatrix} -1 \\ 0 \\ 0 \end{bmatrix} = \begin{bmatrix} 0 \\ 1 \\ 0 \end{bmatrix}$.

The vector $\hat{\mathbf{x}}$ is rotated directly into $+\hat{\mathbf{y}}$, as expected for an active counter-clockwise rotation about $+Z$.

---

## 9. The Double-Cover Property of Quaternions

A fundamental geometric property of quaternions is that the unit 3-sphere $S^3$ is a **double cover of the rotation group $SO(3)$**:

$$\mathbf{R}(q) \equiv \mathbf{R}(-q)$$

Both $q = [w, x, y, z]$ and $-q = [-w, -x, -y, -z]$ represent the **exact same physical 3D orientation**.

### Why?
In the axis-angle formulation:
- $q(\hat{\mathbf{u}}, \theta) = [\cos(\theta/2), \sin(\theta/2)\hat{\mathbf{u}}]$.
- If we rotate by $(\theta + 2\pi)$ (a full additional revolution):
  $$\cos\left(\frac{\theta + 2\pi}{2}\right) = \cos(\theta/2 + \pi) = -\cos(\theta/2)$$
  $$\sin\left(\frac{\theta + 2\pi}{2}\right) = \sin(\theta/2 + \pi) = -\sin(\theta/2)$$
- Thus, rotating by an extra $360^\circ$ negates the quaternion (giving $-q$), while leaving the physical orientation completely unchanged!

> [!IMPORTANT]
> When testing quaternion equality or computing attitude errors, never write `q1 == q2`. Always use physical equivalence:
> `represents_same_rotation(q1, q2)` which checks `q1 ≈ q2 || q1 ≈ -q2`.

---

## 10. Quaternion $\leftrightarrow$ Direction Cosine Matrix (DCM)

### 10.1 Quaternion $\to$ DCM
For a unit quaternion $q = [w, x, y, z]$:
$$C(q) = \begin{bmatrix}
1 - 2(y^2 + z^2) & 2(xy - wz) & 2(xz + wy) \\
2(xy + wz) & 1 - 2(x^2 + z^2) & 2(yz - wx) \\
2(xz - wy) & 2(yz + wx) & 1 - 2(x^2 + y^2)
\end{bmatrix}$$

This matrix satisfies:
- Orthonormality: $C^T C = I$
- Proper rotation: $\det(C) = +1$
- Vector rotation: $C \mathbf{v} \equiv q \otimes [0, \mathbf{v}] \otimes q^*$

### 10.2 DCM $\to$ Quaternion (Shepperd's Stable Algorithm)
Converting an arbitrary DCM to a quaternion requires finding the maximum diagonal element to avoid division by near-zero numbers:
1. Examine the trace $\operatorname{tr}(C) = c_{00} + c_{11} + c_{22}$ and the three diagonal elements $c_{00}, c_{11}, c_{22}$.
2. Branch into four numerically conditioned cases:
   - If $\operatorname{tr}(C) > \max(c_{00}, c_{11}, c_{22})$: compute $w = 0.5\sqrt{1 + \operatorname{tr}(C)}$, then divide by $4w$.
   - If $c_{00}$ is largest: compute $x = 0.5\sqrt{1 + 2c_{00} - \operatorname{tr}(C)}$, then divide by $4x$.
   - If $c_{11}$ is largest: compute $y = 0.5\sqrt{1 + 2c_{11} - \operatorname{tr}(C)}$, then divide by $4y$.
   - If $c_{22}$ is largest: compute $z = 0.5\sqrt{1 + 2c_{22} - \operatorname{tr}(C)}$, then divide by $4z$.
3. Normalize the resulting quaternion to guarantee $\|q\| = 1.0$.

This algorithm is unconditionally robust for all rotation angles ($0^\circ, 90^\circ, 180^\circ$, near-$180^\circ$).

---

## 11. Quaternion Convention Used by AstraDock

To prevent any ambiguity across the AstraDock codebase:

1. **Storage Order:** Scalar-first: $q = [w, x, y, z] = [w, \mathbf{v}]$.
2. **Action / Transformation:** Active rotation: $v' = q \otimes v \otimes q^* = C(q) v$.
3. **Frame Interpretation:** $q_{A\_B}$ represents the orientation of frame $\mathcal{B}$ relative to frame $\mathcal{A}$, mapping vector coordinates from $\mathcal{B}$ into $\mathcal{A}$:
   $$\mathbf{v}_A = q_{A\_B} \otimes [0, \mathbf{v}_B] \otimes q_{A\_B}^* = C_{A\_B} \mathbf{v}_B$$
4. **Composition Order:** Successive rotations chain identically to DCMs:
   $$q_{A\_C} = q_{A\_B} \otimes q_{B\_C} \iff C_{A\_C} = C_{A\_B} C_{B\_C}$$

---

## 12. Summary & Key Takeaways

1. **Decoupled States:** Orbital position ($\mathbf{r}, \mathbf{v}$) and attitude orientation ($q$) are distinct physical state variables.
2. **Singularity-Free Representation:** Unit quaternions represent full 3D rotations globally without gimbal lock singularities.
3. **Double Cover:** Every physical rotation corresponds to a pair of antipodal quaternions $\pm q$.
4. **Active DCM Compatibility:** The quaternion algebra matches the M06 Direction Cosine Matrix conventions exactly.
5. **Next Step:** M08 establishes the kinematic representation. M09 will introduce rigid-body angular velocity, moments of inertia, and torque-free rotational dynamics.
