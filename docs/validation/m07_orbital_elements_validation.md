# Milestone 07: Classical Orbital Elements Validation Report

> **Validation Objective:** Verify that classical Keplerian orbital elements $(a, e, i, \Omega, \omega, \nu)$ and Cartesian states $(\mathbf{r}, \mathbf{v})$ in ECI convert bidirectionally with sub-millimetre round-trip consistency, that analytical orbital properties match independent derivations, and that orbital elements remain invariant under multi-orbit numerical propagation.

---

## 1. Classical Orbital Element Definitions

AstraDock models Keplerian orbital geometry using six classical orbital elements defined relative to the Earth-Centered Inertial (ECI) coordinate frame:

```text
a:  Semi-major axis (m)                 — size of the orbital ellipse
e:  Eccentricity (-)                    — shape / elongation of the ellipse
i:  Inclination (rad)                   — angle between orbital normal h and ECI +Z axis
Ω:  RAAN (rad)                          — right ascension of the ascending node from ECI +X axis
ω:  Argument of periapsis (rad)         — angle in orbital plane from node line to periapsis
ν:  True anomaly (rad)                  — instantaneous angle in orbital plane from periapsis to satellite
```

---

## 2. Bidirectional Conversion Mathematics

### 2.1 State $\to$ Elements
- $\mathbf{h} = \mathbf{r} \times \mathbf{v}$, $h = \|\mathbf{h}\|$, $\hat{\mathbf{h}} = \mathbf{h} / h$
- $\epsilon = \frac{v^2}{2} - \frac{\mu}{r} \implies a = -\frac{\mu}{2\epsilon}$
- $\mathbf{e} = \frac{\mathbf{v} \times \mathbf{h}}{\mu} - \frac{\mathbf{r}}{r}, \quad e = \|\mathbf{e}\|$
- $i = \arccos\left(\operatorname{clamp}\left(\frac{h_z}{h}, -1, 1\right)\right)$
- $\mathbf{n} = \mathbf{k} \times \mathbf{h} = [-h_y, h_x, 0]^T, \quad n = \|\mathbf{n}\|$
- $\Omega = \operatorname{atan2}(n_y, n_x) \pmod{2\pi}$ (if $n > \text{tol}$, else $0$)
- $\omega = \operatorname{atan2}((\mathbf{n} \times \mathbf{e}) \cdot \hat{\mathbf{h}}, \mathbf{n} \cdot \mathbf{e}) \pmod{2\pi}$ (if $n > \text{tol}, e > \text{tol}$)
- $\nu = \operatorname{atan2}((\mathbf{e} \times \mathbf{r}) \cdot \hat{\mathbf{h}}, \mathbf{e} \cdot \mathbf{r}) \pmod{2\pi}$ (if $e > \text{tol}$)

### 2.2 Elements $\to$ State (Perifocal Frame $PQW$)
- Semi-latus rectum: $p = a(1 - e^2)$
- Radial distance: $r = \frac{p}{1 + e\cos\nu}$
- Perifocal state:
  $$\mathbf{r}_{PQW} = \begin{bmatrix} r \cos\nu \\ r \sin\nu \\ 0 \end{bmatrix}, \quad \mathbf{v}_{PQW} = \sqrt{\frac{\mu}{p}} \begin{bmatrix} -\sin\nu \\ e + \cos\nu \\ 0 \end{bmatrix}$$
- Perifocal to ECI rotation matrix: $C_{ECI\_PQW} = R_3(-\Omega) R_1(-i) R_3(-\omega)$
- Cartesian state in ECI: $\mathbf{r}_{ECI} = C_{ECI\_PQW} \mathbf{r}_{PQW}, \quad \mathbf{v}_{ECI} = C_{ECI\_PQW} \mathbf{v}_{PQW}$

---

## 3. Quantitative Round-Trip Verification Results

The bidirectional round-trip transformation ($\text{State} \to \text{Elements} \to \text{State}'$) was tested across six distinct orbital regimes spanning circular, inclined, equatorial, eccentric, and near-circular trajectories:

| Scenario | Input Elements ($a, e, i, \Omega, \omega, \nu$) | Round-Trip Position Error ($\|\mathbf{r}' - \mathbf{r}\|$) | Round-Trip Velocity Error ($\|\mathbf{v}' - \mathbf{v}\|$) | Evaluation |
| :--- | :--- | :---: | :---: | :---: |
| **Case A: Circular Equatorial** | $a=7,000\,\text{km}, e=0, i=0^\circ, \Omega=0^\circ, \omega=0^\circ, \nu=45^\circ$ | $0.000 \times 10^0\,\text{m}$ | $2.034 \times 10^{-12}\,\text{m/s}$ | Exact machine precision |
| **Case B: Circular Inclined** | $a=7,000\,\text{km}, e=0, i=45^\circ, \Omega=60^\circ, \omega=0^\circ, \nu=90^\circ$ | $2.794 \times 10^{-9}\,\text{m}$ | $2.315 \times 10^{-12}\,\text{m/s}$ | Sub-nanometre |
| **Case C: Elliptical Equatorial** | $a=10,000\,\text{km}, e=0.2, i=0^\circ, \Omega=0^\circ, \omega=30^\circ, \nu=60^\circ$ | $4.657 \times 10^{-10}\,\text{m}$ | $9.095 \times 10^{-13}\,\text{m/s}$ | Sub-nanometre |
| **Case D: Elliptical Inclined** | $a=10,000\,\text{km}, e=0.2, i=45^\circ, \Omega=120^\circ, \omega=60^\circ, \nu=30^\circ$ | $2.281 \times 10^{-9}\,\text{m}$ | $2.193 \times 10^{-12}\,\text{m/s}$ | Sub-nanometre |
| **Case E: High-Inclination (SSO)** | $a=8,000\,\text{km}, e=0.15, i=98^\circ, \Omega=45^\circ, \omega=270^\circ, \nu=150^\circ$ | $2.281 \times 10^{-9}\,\text{m}$ | $1.575 \times 10^{-12}\,\text{m/s}$ | Sub-nanometre |
| **Case F: Near-Circular** | $a=7,000\,\text{km}, e=10^{-5}, i=30^\circ, \Omega=15^\circ, \omega=45^\circ, \nu=120^\circ$ | $2.134 \times 10^{-9}\,\text{m}$ | $2.652 \times 10^{-12}\,\text{m/s}$ | Sub-nanometre |

---

## 4. Independent Python Oracle Comparison

An independent pure-Python astrodynamics oracle ([`python/audit/independent_orbital_elements.py`](file:///C:/Users/user/AstraDock/python/audit/independent_orbital_elements.py)) was executed without importing or calling AstraDock C++ code:

- **Case D (Inclined Elliptical Orbit) Cross-Check:**
  * Analytical Orbit Radius at $\nu = 30^\circ$: $8,182,712.603\,\text{m}$ (C++ matches Python to 14 digits)
  * Analytical Speed at $\nu = 30^\circ$: $7,587.158167\,\text{m/s}$ (C++ matches Python to 14 digits)
  * Analytical Orbital Period: $9,952.014050\,\text{s}$ (165.867 min)
  * Maximum cross-language element discrepancy: $< 1.0 \times 10^{-14}$ across all components.

---

## 5. Multi-Orbit Invariance Under Numerical Propagation

Under ideal two-body gravitational dynamics, orbital elements $a, e, i, \Omega, \omega$ are strictly conserved invariants of motion.

An inclined eccentric orbit (Case D: $a = 10,000\,\text{km}, e = 0.2, i = 45^\circ, \Omega = 120^\circ, \omega = 60^\circ$) was propagated for **3 full orbital periods** ($29,856\,\text{s}$, $2,987$ RK4 steps at $\Delta t = 10\,\text{s}$):

| Orbital Element | Initial Reference Value | Maximum Measured Drift Over 3 Orbits | Relative Invariant Error |
| :--- | :---: | :---: | :---: |
| **Semi-Major Axis ($a$)** | $10,000,000.000\,\text{m}$ | $4.30 \times 10^{-4}\,\text{m}$ | $4.30 \times 10^{-11}$ |
| **Eccentricity ($e$)** | $0.200000000$ | $8.22 \times 10^{-11}$ | $4.11 \times 10^{-10}$ |
| **Inclination ($i$)** | $45.000000000^\circ$ | $2.66 \times 10^{-15}\,\text{rad}$ ($1.53 \times 10^{-13\,\circ}$) | $3.39 \times 10^{-15}$ |
| **RAAN ($\Omega$)** | $120.000000000^\circ$ | $1.78 \times 10^{-15}\,\text{rad}$ ($1.02 \times 10^{-13\,\circ}$) | $8.48 \times 10^{-16}$ |
| **Arg of Periapsis ($\omega$)** | $60.000000000^\circ$ | $2.85 \times 10^{-9}\,\text{rad}$ ($1.64 \times 10^{-7\,\circ}$) | $2.72 \times 10^{-9}$ |
| **Max Round-Trip Pos Err** | -- | $1.27 \times 10^{-8}\,\text{m}$ | $1.27 \times 10^{-15}$ |
| **Max Round-Trip Vel Err** | -- | $6.81 \times 10^{-12}\,\text{m/s}$ | $8.97 \times 10^{-16}$ |

---

## 6. Generated Visualizations

1. **Plot A ([`artifacts/figures/m07_orbital_geometry.png`](file:///C:/Users/user/AstraDock/artifacts/figures/m07_orbital_geometry.png)):** 3D perspective visualization showing Earth sphere, equatorial plane, ascending node line ($\Omega$), periapsis direction ($\omega$), instantaneous spacecraft position ($\nu$), and the full 3D orbital ellipse in ECI.
2. **Plot B ([`artifacts/figures/m07_elements_time_evolution.png`](file:///C:/Users/user/AstraDock/artifacts/figures/m07_elements_time_evolution.png)):** Time-series subplots demonstrating that $a(t), e(t), i(t), \Omega(t), \omega(t)$ remain strictly flat lines under two-body propagation, while true anomaly $\nu(t)$ advances smoothly through successive orbital revolutions.

---

## 7. Conclusions & Boundary of Scope

- **Milestone 07 Complete:** AstraDock now possesses an independently validated, non-singular bidirectional conversion between 6-state Cartesian vectors and classical orbital elements, complete with analytical geometry helpers, angle normalization, and explicit singularity policies.
- **Scope Boundary Maintained:** Perturbations ($J_2$, drag), non-spherical harmonics, attitude kinematics/quaternions, equinoctial elements, and relative motion models remain strictly out of scope for M07 and are scheduled for subsequent roadmap milestones.
