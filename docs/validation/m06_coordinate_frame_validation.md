# M06 Coordinate Frame Validation Report

> **Scope:** This report documents the mathematical and numerical validation of coordinate frames, reference frame basis triads, Direction Cosine Matrices (DCMs), and vector transformations implemented in Milestone 06.

---

## 1. Frame Conventions and Definitions

### 1.1 Earth-Centered Inertial (ECI)
- **Origin:** Earth center of mass.
- **Orientation:** Idealized non-rotating Cartesian coordinate axes with $+Z_{ECI}$ aligned with Earth's north celestial pole and $+X_{ECI}, +Y_{ECI}$ spanning the equatorial plane.
- **Idealization Boundary:** No epoch, precession, nutation, Earth orientation parameters, or Earth rotation are implemented in this milestone.

### 1.2 Local Vertical Local Horizontal (LVLH)
Derived instantaneously from spacecraft state $\mathbf{r}_{ECI}, \mathbf{v}_{ECI}$:

$$\hat{\mathbf{r}} = \frac{\mathbf{r}_{ECI}}{\|\mathbf{r}_{ECI}\|}, \quad \mathbf{h} = \mathbf{r}_{ECI} \times \mathbf{v}_{ECI}, \quad \hat{\mathbf{h}} = \frac{\mathbf{h}}{\|\mathbf{h}\|}, \quad \hat{\mathbf{t}} = \hat{\mathbf{h}} \times \hat{\mathbf{r}}$$

$$\mathbf{e}_r = \hat{\mathbf{r}} \quad (\text{Radial}), \qquad \mathbf{e}_t = \hat{\mathbf{t}} \quad (\text{Along-Track}), \qquad \mathbf{e}_h = \hat{\mathbf{h}} \quad (\text{Orbit-Normal})$$

Triad orientation: Right-handed ($\mathbf{e}_r \times \mathbf{e}_t = \mathbf{e}_h$).

### 1.3 Direction Cosine Matrix Convention
AstraDock defines $C_{A\_B}$ such that:

$$\mathbf{v}_A = C_{A\_B} \cdot \mathbf{v}_B$$

where $\mathbf{v}_B$ and $\mathbf{v}_A$ represent the coordinate components of the same physical vector in Frame B and Frame A, respectively.

- **ECI to LVLH Matrix ($C_{LVLH\_ECI}$):**
  $$C_{LVLH\_ECI} = \begin{bmatrix} \mathbf{e}_r^T \\ \mathbf{e}_t^T \\ \mathbf{e}_h^T \end{bmatrix} = \begin{bmatrix} e_{r,x} & e_{r,y} & e_{r,z} \\ e_{t,x} & e_{t,y} & e_{t,z} \\ e_{h,x} & e_{h,y} & e_{h,z} \end{bmatrix}$$
  *Basis vectors are the **rows** of $C_{LVLH\_ECI}$.*

- **LVLH to ECI Matrix ($C_{ECI\_LVLH}$):**
  $$C_{ECI\_LVLH} = C_{LVLH\_ECI}^T = \begin{bmatrix} \mathbf{e}_r & \mathbf{e}_t & \mathbf{e}_h \end{bmatrix} = \begin{bmatrix} e_{r,x} & e_{t,x} & e_{h,x} \\ e_{r,y} & e_{t,y} & e_{h,y} \\ e_{r,z} & e_{t,z} & e_{h,z} \end{bmatrix}$$
  *Basis vectors are the **columns** of $C_{ECI\_LVLH}$.*

---

## 2. Analytical Regression Baselines

For a 500 km circular equatorial orbit ($R = 6,878,137$ m, $V = 7,612.608173$ m/s), the analytical reference frames at four orbital quadrants are:

| Orbital Phase $\theta$ | ECI Position $\mathbf{r}$ (km) | ECI Velocity $\mathbf{v}$ (km/s) | LVLH $\mathbf{e}_r$ | LVLH $\mathbf{e}_t$ | LVLH $\mathbf{e}_h$ | DCM $C_{LVLH\_ECI}$ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **$0^\circ$ (Initial)** | $[6878.137, 0, 0]$ | $[0, 7.6126, 0]$ | $[+1, 0, 0]$ | $[0, +1, 0]$ | $[0, 0, +1]$ | $\begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix}$ |
| **$90^\circ$ (Quarter)** | $[0, 6878.137, 0]$ | $[-7.6126, 0, 0]$ | $[0, +1, 0]$ | $[-1, 0, 0]$ | $[0, 0, +1]$ | $\begin{bmatrix} 0 & 1 & 0 \\ -1 & 0 & 0 \\ 0 & 0 & 1 \end{bmatrix}$ |
| **$180^\circ$ (Half)** | $[-6878.137, 0, 0]$ | $[0, -7.6126, 0]$ | $[-1, 0, 0]$ | $[0, -1, 0]$ | $[0, 0, +1]$ | $\begin{bmatrix} -1 & 0 & 0 \\ 0 & -1 & 0 \\ 0 & 0 & 1 \end{bmatrix}$ |
| **$270^\circ$ (Three-Quarter)** | $[0, -6878.137, 0]$ | $[+7.6126, 0, 0]$ | $[0, -1, 0]$ | $[+1, 0, 0]$ | $[0, 0, +1]$ | $\begin{bmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 1 \end{bmatrix}$ |

All quadrant baselines are verified in `test_coordinate_frames.cpp` to within machine epsilon ($\le 1.0 \times 10^{-12}$).

---

## 3. Measured Numerical Tolerances and Invariants

Across a full 500 km RK4 orbit integration (569 samples, $T = 5,676.978$ s, $\Delta t = 10.0$ s), frame metrics were measured with the following results:

| Verification Metric | Definition | Theoretical Target | Measured Maximum Error | Status |
| :--- | :--- | :--- | :--- | :--- |
| **Basis Orthonormality Error** | $\max_{i, j} |\|\mathbf{e}_i \cdot \mathbf{e}_j - \delta_{ij}\||$ | $0.0$ | $2.22 \times 10^{-16}$ | **PASS** |
| **Determinant Error** | $|\det(C) - 1.0|$ | $0.0$ | $4.44 \times 10^{-16}$ | **PASS** |
| **Right-Handed Scalar Triple Product** | $|(\mathbf{e}_r \times \mathbf{e}_t) \cdot \mathbf{e}_h - 1.0|$ | $0.0$ | $2.22 \times 10^{-16}$ | **PASS** |
| **Position Norm Invariance** | $|\|\mathbf{r}_{LVLH}\| - \|\mathbf{r}_{ECI}\||$ | $0.0\,\text{m}$ | $1.86 \times 10^{-9}\,\text{m}$ | **PASS** |
| **Round-Trip Transform Error** | $\|C^T (C \mathbf{r}) - \mathbf{r}\|$ | $0.0\,\text{m}$ | $2.63 \times 10^{-9}\,\text{m}$ | **PASS** |
| **Velocity Along-Track Invariant** | $|v_{LVLH, y} - v_c|$ | $0.0\,\text{m/s}$ | $2.14 \times 10^{-3}\,\text{m/s}$ | **PASS** |
| **Velocity Out-of-Plane Invariant** | $|v_{LVLH, z}|$ | $0.0\,\text{m/s}$ | $0.00\,\text{m/s}$ | **PASS** |

*Note:* The position norm difference ($1.86 \times 10^{-9}$ m) on a vector magnitude of $6.88 \times 10^6$ m corresponds to a relative error of $2.7 \times 10^{-16}$, directly at double-precision floating-point precision limit.

---

## 4. Degenerate and Error Handling

The implementation explicitly guards against unphysical or mathematically ill-conditioned states:

| Degenerate Condition | Physical Cause | Implementation Action | Verified Exception |
| :--- | :--- | :--- | :--- |
| $\|\mathbf{r}\| = 0$ | Spacecraft at central body origin | Reject | `std::domain_error` |
| $\|\mathbf{v}\| = 0$ | Zero velocity (no orbital plane) | Reject | `std::domain_error` |
| $\mathbf{r} \parallel \mathbf{v}$ | Collinear radial trajectory ($\mathbf{h} = 0$) | Reject | `std::domain_error` |
| $\text{NaN} \in \mathbf{r}$ or $\mathbf{v}$ | Numerical arithmetic failure | Reject | `std::domain_error` |
| $\pm\infty \in \mathbf{r}$ or $\mathbf{v}$ | Numerical divergence / overflow | Reject | `std::domain_error` |

All degenerate cases are verified in deterministic unit tests.

---

## 5. Artifacts and Generated Evidence

- C++ Data Export: `artifacts/data/m06_frames.csv`
- Plot A (Single-Point Frame Triads): `artifacts/figures/m06_frame_single_point.png`
- Plot B (Orbital Triad Evolution): `artifacts/figures/m06_frame_orbit_evolution.png`
- Plot C (Vector Transformation Sanity): `artifacts/figures/m06_frame_transform_sanity.png`

---

## 6. Test Suite Status

Current project test baseline: **91 / 91 tests passing** (100%).
- 74 baseline tests from M01–M05 (vectors, two-body gravity, integrators, orbital propagation, validation diagnostics)
- 9 Matrix3 unit tests (factories, operators, transpose, determinant, trace, orthonormality, finite checks)
- 8 Coordinate frame unit/integration tests (basis validation, circular quadrant tests, inclined orbits, DCM transformation, norm preservation, round-trip, degenerate rejections, full-orbit propagation)
