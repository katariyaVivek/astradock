# Lesson 006: Coordinate Frames and Orbital Reference Frames

## Learning objective

M05 established that AstraDock propagates orbital trajectories with characterized numerical accuracy.
M06 introduces the single most vital principle in guidance, navigation, control, and spacecraft engineering:

> **A physical vector is not just three numbers. A vector is meaningless without knowing the reference frame in which it is expressed.**

```text
same physical vector
        ↓
different coordinate frame
        ↓
different numerical components
        ↓
same physical magnitude
```

---

## 1. Physical Vectors vs Numerical Components

A **physical vector** represents a geometric quantity possessing magnitude and direction (e.g., position relative to Earth center, instantaneous orbital velocity, thrust force, sunlight vector).

A physical vector exists independently of any coordinate system.

However, to compute with a vector on a digital computer, we must choose a set of coordinate axes (a **reference frame**) and measure the projections of the physical vector onto those axes:

$$\mathbf{v} = v_1 \mathbf{e}_1 + v_2 \mathbf{e}_2 + v_3 \mathbf{e}_3$$

The three numbers $[v_1, v_2, v_3]^T$ are the **components** of $\mathbf{v}$ relative to the basis triad $\{\mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3\}$. If you change the basis triad, the three numbers change, even though the physical vector itself has not changed.

### The Golden Rule of Spacecraft Engineering

> **Never silently mix reference frames.**
> Adding an ECI position vector to an LVLH relative vector or a body-frame force vector is physically meaningless.

---

## 2. Coordinate Frame Hierarchy in AstraDock

AstraDock establishes the following clear frame hierarchy:

```text
                     ECI (Earth-Centered Inertial)
                                  │
                                  │ orbital state (r, v)
                                  ▼
                   LVLH (Local Vertical Local Horizontal)
                                  │
                                  │ planned attitude milestone
                                  ▼
                   Spacecraft Body Frame (PLANNED)
```

---

## 3. Earth-Centered Inertial (ECI) Frame

In AstraDock, **ECI** is the fundamental non-rotating Cartesian coordinate frame used for solving the orbital equations of motion:

- **Origin:** Earth's center of mass.
- **Fundamental Plane:** Earth's equatorial plane (in this idealized model).
- **Axes:**
  - $+X_{ECI}$: Fixed reference direction in the equatorial plane.
  - $+Y_{ECI}$: In the equatorial plane, $90^\circ$ east of $+X_{ECI}$.
  - $+Z_{ECI}$: Perpendicular to the equatorial plane, pointing along Earth's north celestial pole.

### Idealization Note
At this educational stage:
- No operational epoch, precession, nutation, or Earth orientation parameters (EOP) are modeled.
- No Earth rotation or ECEF (Earth-Centered Earth-Fixed) frame is implemented.
- We do not claim high-fidelity realization of ICRF/J2000. It is an idealized Newtonian inertial frame for two-body orbital dynamics.

---

## 4. Local Vertical Local Horizontal (LVLH) Frame

While orbital equations of motion are simplest to integrate in an inertial frame (ECI), spacecraft operations, rendezvous, relative navigation, sensor pointing, and docking are far easier to reason about in an **orbital local frame** attached to the spacecraft.

The **LVLH frame** moves with the spacecraft and rotates as the spacecraft orbits Earth.

### 4.1 Basis Construction

Given the spacecraft instantaneous state in ECI coordinates:
- Position: $\mathbf{r}_{ECI}$
- Velocity: $\mathbf{v}_{ECI}$

We construct three orthonormal unit vectors:

1. **Radial unit vector ($\hat{\mathbf{r}}$ / $\mathbf{e}_r$):** Points outward from Earth's center through the spacecraft:
   $$\hat{\mathbf{r}} = \frac{\mathbf{r}_{ECI}}{\|\mathbf{r}_{ECI}\|}$$

2. **Specific angular momentum ($\mathbf{h}$):** Normal to the instantaneous orbital plane:
   $$\mathbf{h} = \mathbf{r}_{ECI} \times \mathbf{v}_{ECI}$$

3. **Orbit-normal unit vector ($\hat{\mathbf{h}}$ / $\mathbf{e}_h$):**
   $$\hat{\mathbf{h}} = \frac{\mathbf{h}}{\|\mathbf{h}\|}$$

4. **Along-track unit vector ($\hat{\mathbf{t}}$ / $\mathbf{e}_t$):** Completes the right-handed triad in the forward direction of motion:
   $$\hat{\mathbf{t}} = \hat{\mathbf{h}} \times \hat{\mathbf{r}}$$

### 4.2 AstraDock LVLH Axis Convention

AstraDock explicitly adopts the following standard right-handed definition:

| Axis | Unit Vector | Physical Meaning | Direction |
| :--- | :--- | :--- | :--- |
| **LVLH X** | $\mathbf{e}_r = \hat{\mathbf{r}}$ | Radial / Local Vertical | Outward from Earth center through spacecraft |
| **LVLH Y** | $\mathbf{e}_t = \hat{\mathbf{t}}$ | Along-Track / Local Horizontal | Forward in direction of orbital motion |
| **LVLH Z** | $\mathbf{e}_h = \hat{\mathbf{h}}$ | Cross-Track / Orbit-Normal | Normal to orbital plane (aligned with $\mathbf{h}$) |

Right-handed orientation check:
$$\mathbf{e}_r \times \mathbf{e}_t = \hat{\mathbf{r}} \times (\hat{\mathbf{h}} \times \hat{\mathbf{r}}) = (\hat{\mathbf{r}} \cdot \hat{\mathbf{r}})\hat{\mathbf{h}} - (\hat{\mathbf{r}} \cdot \hat{\mathbf{h}})\hat{\mathbf{r}} = 1 \cdot \hat{\mathbf{h}} - 0 = \hat{\mathbf{h}} = \mathbf{e}_h$$

---

## 5. Orthonormal Basis Representation: `FrameBasis`

A 3D coordinate frame basis is represented in AstraDock as three unit vectors expressed in the parent frame:

```cpp
namespace astradock::frames {

struct FrameBasis {
    math::Vector3 x; // e_r (radial)
    math::Vector3 y; // e_t (along-track)
    math::Vector3 z; // e_h (orbit-normal)
};

}
```

A valid basis triad must satisfy:
1. **Unit length:** $\|\mathbf{x}\| = 1, \|\mathbf{y}\| = 1, \|\mathbf{z}\| = 1$
2. **Mutual orthogonality:** $\mathbf{x} \cdot \mathbf{y} = 0, \mathbf{y} \cdot \mathbf{z} = 0, \mathbf{z} \cdot \mathbf{x} = 0$
3. **Right-handed orientation:** $\mathbf{x} \times \mathbf{y} = \mathbf{z}$, and $(\mathbf{x} \times \mathbf{y}) \cdot \mathbf{z} = +1$

---

## 6. Direction Cosine Matrices (DCM)

A **Direction Cosine Matrix (DCM)** is a $3 \times 3$ orthogonal matrix that maps the numerical components of a vector from one coordinate frame to another.

### 6.1 Transformation Direction Convention

To avoid ambiguity, AstraDock uses the explicit notation:

$$\mathbf{v}_A = C_{A\_B} \cdot \mathbf{v}_B$$

`C_A_B` (or `dcm_a_from_b`) transforms coordinate components of a vector expressed in Frame B into coordinate components of the **same physical vector** expressed in Frame A.

### 6.2 Constructing $C_{LVLH\_ECI}$

Let $\mathbf{v}$ be a physical vector with components $\mathbf{v}_{ECI} = [v_X, v_Y, v_Z]^T$.
Its components along the LVLH axes are given by scalar projections (dot products) with the basis vectors:

$$v_{LVLH, x} = \mathbf{e}_r \cdot \mathbf{v}_{ECI}$$
$$v_{LVLH, y} = \mathbf{e}_t \cdot \mathbf{v}_{ECI}$$
$$v_{LVLH, z} = \mathbf{e}_h \cdot \mathbf{v}_{ECI}$$

In matrix form:

$$\begin{bmatrix} v_{LVLH, x} \\ v_{LVLH, y} \\ v_{LVLH, z} \end{bmatrix} = \begin{bmatrix} e_{r,x} & e_{r,y} & e_{r,z} \\ e_{t,x} & e_{t,y} & e_{t,z} \\ e_{h,x} & e_{h,y} & e_{h,z} \end{bmatrix} \begin{bmatrix} v_{ECI, x} \\ v_{ECI, y} \\ v_{ECI, z} \end{bmatrix}$$

Therefore:
> **The rows of $C_{LVLH\_ECI}$ are the LVLH basis vectors expressed in ECI coordinates.**

$$C_{LVLH\_ECI} = \begin{bmatrix} \mathbf{e}_r^T \\ \mathbf{e}_t^T \\ \mathbf{e}_h^T \end{bmatrix}$$

### 6.3 Inverse Transformation $C_{ECI\_LVLH}$

Because the basis vectors are orthonormal, the DCM is orthogonal:

$$C_{LVLH\_ECI} \cdot C_{LVLH\_ECI}^T = I$$

Hence, the inverse transform matrix is simply the matrix **transpose**:

$$C_{ECI\_LVLH} = C_{LVLH\_ECI}^{-1} = C_{LVLH\_ECI}^T = \begin{bmatrix} \mathbf{e}_r & \mathbf{e}_t & \mathbf{e}_h \end{bmatrix}$$

> **The columns of $C_{ECI\_LVLH}$ are the LVLH basis vectors expressed in ECI coordinates.**

---

## 7. Fundamental Properties of Rotation Matrices

For any proper Direction Cosine Matrix $C$:

1. **Orthogonality:** $C^T C = C C^T = I$
2. **Inverse equals Transpose:** $C^{-1} = C^T$
3. **Determinant:** $\det(C) = +1.0$ (a reflection matrix has $\det = -1.0$; a proper rotation has $\det = +1.0$)
4. **Norm Preservation:** For any vector $\mathbf{v}$,
   $$\|C \mathbf{v}\|^2 = (C \mathbf{v})^T (C \mathbf{v}) = \mathbf{v}^T C^T C \mathbf{v} = \mathbf{v}^T I \mathbf{v} = \|\mathbf{v}\|^2$$
   $$\|C \mathbf{v}\| = \|\mathbf{v}\|$$
   A coordinate rotation changes component numbers, but preserves vector length.
5. **Round-Trip Identity:** $C^T (C \mathbf{v}) = \mathbf{v}$

---

## 8. Worked Example: Same Vector, Different Numbers

Consider a spacecraft in a 500 km circular equatorial orbit ($r = 6,878,137$ m, $v = 7,612.608$ m/s).

### State A: At $\theta = 0^\circ$ ($t = 0$)
- $\mathbf{r}_{ECI} = [6878137.0, 0.0, 0.0]^T$ m
- $\mathbf{v}_{ECI} = [0.0, 7612.608, 0.0]^T$ m/s

Basis vectors in ECI:
- $\mathbf{e}_r = [1, 0, 0]^T$
- $\mathbf{e}_t = [0, 1, 0]^T$
- $\mathbf{e}_h = [0, 0, 1]^T$

Transforming position and velocity into LVLH:
$$\mathbf{r}_{LVLH} = C_{LVLH\_ECI} \mathbf{r}_{ECI} = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix} \begin{bmatrix} 6878137.0 \\ 0.0 \\ 0.0 \end{bmatrix} = \begin{bmatrix} 6878137.0 \\ 0.0 \\ 0.0 \end{bmatrix}\,\text{m}$$
$$\mathbf{v}_{LVLH} = C_{LVLH\_ECI} \mathbf{v}_{ECI} = \begin{bmatrix} 0.0 \\ 7612.608 \\ 0.0 \end{bmatrix}\,\text{m/s}$$

### State B: At Quarter Orbit ($\theta = 90^\circ$, $t \approx 1419.2$ s)
- $\mathbf{r}_{ECI} = [0.0, 6878137.0, 0.0]^T$ m
- $\mathbf{v}_{ECI} = [-7612.608, 0.0, 0.0]^T$ m/s

Basis vectors in ECI:
- $\mathbf{e}_r = [0, 1, 0]^T$ (points along $+Y_{ECI}$)
- $\mathbf{e}_t = [-1, 0, 0]^T$ (points along $-X_{ECI}$)
- $\mathbf{e}_h = [0, 0, 1]^T$ (points along $+Z_{ECI}$)

DCM:
$$C_{LVLH\_ECI} = \begin{bmatrix} 0 & 1 & 0 \\ -1 & 0 & 0 \\ 0 & 0 & 1 \end{bmatrix}$$

Transforming position and velocity into LVLH:
$$\mathbf{r}_{LVLH} = \begin{bmatrix} 0 & 1 & 0 \\ -1 & 0 & 0 \\ 0 & 0 & 1 \end{bmatrix} \begin{bmatrix} 0.0 \\ 6878137.0 \\ 0.0 \end{bmatrix} = \begin{bmatrix} 6878137.0 \\ 0.0 \\ 0.0 \end{bmatrix}\,\text{m}$$
$$\mathbf{v}_{LVLH} = \begin{bmatrix} 0 & 1 & 0 \\ -1 & 0 & 0 \\ 0 & 0 & 1 \end{bmatrix} \begin{bmatrix} -7612.608 \\ 0.0 \\ 0.0 \end{bmatrix} = \begin{bmatrix} 0.0 \\ 7612.608 \\ 0.0 \end{bmatrix}\,\text{m/s}$$

### Key Observations:
1. In **ECI**, the position and velocity vectors oscillate sinusoidally as the spacecraft orbits:
   - At $t=0$: $r_X = +R, r_Y = 0$; $v_X = 0, v_Y = +V$.
   - At quarter orbit: $r_X = 0, r_Y = +R$; $v_X = -V, v_Y = 0$.
2. In **LVLH**, the position and velocity representations are **constant**:
   - $\mathbf{r}_{LVLH} = [R, 0, 0]^T$ (purely radial).
   - $\mathbf{v}_{LVLH} = [0, V, 0]^T$ (purely along-track).
3. The physical vector norms are identical in both frames:
   $$\|\mathbf{r}_{ECI}\| = \|\mathbf{r}_{LVLH}\| = 6,878,137.0\,\text{m}$$
   $$\|\mathbf{v}_{ECI}\| = \|\mathbf{v}_{LVLH}\| = 7,612.608\,\text{m/s}$$

---

## 9. Critical Distinction: The LVLH Frame Rotates

The ECI frame is inertial (fixed in space). The LVLH frame is **non-inertial** because its axes rotate as the spacecraft moves along the orbit:

$$\boldsymbol{\omega}_{LVLH/ECI} \approx \begin{bmatrix} 0 \\ 0 \\ \dot{\theta} \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \\ n \end{bmatrix}$$

where $n = \sqrt{\mu / r^3}$ is the orbital mean motion.

### Geometric Coordinate Transformation vs Rotating-Frame Kinematics
In M06, `transform_eci_to_lvlh` performs an instantaneous **geometric projection**:
$$\mathbf{v}_{LVLH} = C_{LVLH\_ECI} \cdot \mathbf{v}_{ECI}$$
This gives the components of the inertial velocity resolved along the instantaneous LVLH axes.

It is **not** the relative velocity $\mathbf{v}_{rel}$ that would be observed by an observer fixed to the rotating LVLH frame:
$$\left(\frac{d\mathbf{r}}{dt}\right)_{LVLH} = \left(\frac{d\mathbf{r}}{dt}\right)_{ECI} - \boldsymbol{\omega}_{LVLH/ECI} \times \mathbf{r}$$
Rotating-frame kinematics and relative dynamics (Clohessy-Wiltshire equations) will be formally studied in M15.

---

## 10. Frame Rotation vs Spacecraft Attitude Rotation

Do not confuse the rotation of the LVLH frame with the rotation of the spacecraft body:
- **LVLH-pointing (Nadir pointing):** A spacecraft maintaining a fixed orientation relative to LVLH (e.g. camera pointing at Earth) rotates in inertial space at rate $\omega = n$.
- **Inertially-pointing (Sun pointing / Star pointing):** A spacecraft maintaining a fixed orientation in ECI rotates relative to LVLH at rate $\omega_{rel} = -n$.

Spacecraft body frames, attitude kinematics, and quaternions will be implemented in M07–M09.

---

## 11. Degenerate States & Safety

The LVLH frame is physically undefined under the following conditions:
1. **Zero position magnitude ($\|\mathbf{r}\| = 0$):** The radial direction $\hat{\mathbf{r}} = \mathbf{r} / \|\mathbf{r}\|$ is undefined at the central body origin.
2. **Zero angular momentum ($\|\mathbf{r} \times \mathbf{v}\| = 0$):** Occurs when velocity is zero or collinear with position (pure radial trajectory). The orbital plane normal $\hat{\mathbf{h}}$ is undefined.
3. **Non-finite numbers:** NaN or Inf in state vectors must fail immediately with `std::domain_error`.

AstraDock rejects these states with explicit exceptions rather than silently producing corrupt basis vectors.

---

## 12. Summary of M06 Capabilities

- `astradock::math::Matrix3`: Minimal 3x3 double-precision matrix arithmetic, transposition, determinants, trace, and orthonormality checks.
- `astradock::frames::FrameBasis`: Triad representation with right-handed orthonormality validation.
- `astradock::frames::compute_lvlh_basis`: Deterministic LVLH basis construction from orbital state.
- `astradock::frames::dcm_lvlh_from_eci` / `dcm_eci_from_lvlh`: Explicit DCM construction.
- `astradock::frames::transform_eci_to_lvlh` / `transform_lvlh_to_eci`: Geometric vector coordinate transformations.
