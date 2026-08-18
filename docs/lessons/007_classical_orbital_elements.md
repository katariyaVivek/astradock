# Lesson 007: Classical Orbital Elements and State Conversion

> *"A Cartesian state tells you where a spacecraft is and how fast it is moving at one instant. Orbital elements tell you what orbit it is on."*

---

## 1. The Dual Representation of an Orbit

In aerospace engineering and orbital mechanics, a spacecraft's instantaneous translational motion can be fully described in two distinct mathematical representations:

1. **Cartesian State Vector ($\mathbf{r}, \mathbf{v}$):**
   - 6 components: $[x, y, z, v_x, v_y, v_z]^T$.
   - Expressed in an Earth-Centered Inertial (ECI) coordinate frame.
   - **Advantage:** Unambiguous, non-singular, and directly suitable for numerical ODE integrators (such as RK4).
   - **Limitation:** Obscures the geometric shape, orientation, and periodicity of the trajectory. Two Cartesian states that look completely different may belong to the exact same closed orbit at different times.

2. **Classical Orbital Elements ($a, e, i, \Omega, \omega, \nu$):**
   - 6 Keplerian parameters: semi-major axis, eccentricity, inclination, RAAN, argument of periapsis, and true anomaly.
   - **Advantage:** Directly describes the size, eccentricity, spatial orientation, and in-orbit position. In an ideal two-body gravitational field, **five of the six elements are constant invariants of motion**, while only the true anomaly ($\nu$) evolves over time.
   - **Limitation:** Exhibits mathematical coordinate singularities for circular ($e=0$) and equatorial ($i=0$) orbits.

```text
               +-------------------------------------------+
               |        Cartesian State (ECI Frame)        |
               |          r = [x, y, z]  (m)               |
               |          v = [vx, vy, vz] (m/s)           |
               +-------------------------------------------+
                                   |   ^
            state_to_classical_    |   |  classical_elements_
            elements()             |   |  to_state()
                                   v   |
               +-------------------------------------------+
               |        Classical Keplerian Elements       |
               |  a:  Semi-Major Axis (m)                  |
               |  e:  Eccentricity (-)                     |
               |  i:  Inclination (rad)                    |
               |  Ω:  Right Ascension of Ascending Node    |
               |  ω:  Argument of Periapsis (rad)          |
               |  ν:  True Anomaly (rad)                   |
               +-------------------------------------------+
```

---

## 2. The Six Classical Orbital Elements

```
                          Orbit Normal (h)
                                 ^
                                 |   / Apoapsis (ν = 180°)
                                 |  /
                                 | /
                     +-----------|/-----------+
                    /            O (Earth)     \
                   /            /               \
                  /            /                 \
                 +------------/-------------------+
                             / Periapsis (ν = 0°)
                            v
                     Eccentricity Vector (e)
```

| Element | Symbol | Name | Physical Meaning | Geometric Role | Valid Domain |
| :--- | :---: | :--- | :--- | :--- | :---: |
| **Semi-Major Axis** | $a$ | Size / Scale | Mean orbital radius; sets orbital energy and period $T = 2\pi\sqrt{a^3/\mu}$. | Semi-major axis of the ellipse | $a > 0$ (bound) |
| **Eccentricity** | $e$ | Shape | Deviation of the orbit from a circle ($e=0$ is circle, $0<e<1$ is ellipse). | Ellipse elongation | $0 \le e < 1$ |
| **Inclination** | $i$ | Tilt | Angle between orbital angular momentum $\mathbf{h}$ and Earth's north pole $+Z_{ECI}$. | Orbital plane tilt | $0 \le i \le \pi$ |
| **RAAN** | $\Omega$ | Swivel | Angle in equatorial plane from $+X_{ECI}$ (vernal equinox) to ascending node. | Line of nodes orientation | $[0, 2\pi)$ |
| **Arg of Periapsis**| $\omega$| Orientation| Angle in orbital plane from ascending node to periapsis (closest approach). | Major axis orientation | $[0, 2\pi)$ |
| **True Anomaly** | $\nu$ | Position | Instantaneous angle in orbital plane from periapsis to spacecraft position $\mathbf{r}$.| In-orbit progress | $[0, 2\pi)$ |

---

## 3. Fundamental Vector Quantities

To derive the orbital elements from $\mathbf{r}$ and $\mathbf{v}$, we compute four foundational astrodynamics vectors:

### 3.1 Specific Angular Momentum Vector ($\mathbf{h}$)
$$\mathbf{h} = \mathbf{r} \times \mathbf{v}$$
- **Physical Meaning:** Orbital angular momentum per unit mass ($\text{m}^2/\text{s}$).
- **Geometry:** $\mathbf{h}$ is strictly perpendicular to the instantaneous orbital plane. Because gravity is a central force ($\mathbf{r} \times \mathbf{a} = \mathbf{0}$), $\mathbf{h}$ is a conserved vector invariant in ideal two-body motion.
- Unit vector: $\hat{\mathbf{h}} = \mathbf{h} / \|\mathbf{h}\|$.

### 3.2 Ascending Node Vector ($\mathbf{n}$)
$$\mathbf{n} = \mathbf{k} \times \mathbf{h} = \begin{bmatrix} 0 \\ 0 \\ 1 \end{bmatrix} \times \begin{bmatrix} h_x \\ h_y \\ h_z \end{bmatrix} = \begin{bmatrix} -h_y \\ h_x \\ 0 \end{bmatrix}$$
- **Physical Meaning:** Points along the line of intersection between the equatorial plane ($XY_{ECI}$) and the orbital plane, directed toward the **ascending node** (where the spacecraft crosses from south to north).
- Magnitude: $n = \|\mathbf{n}\| = \sqrt{h_x^2 + h_y^2}$.

### 3.3 Eccentricity Vector ($\mathbf{e}$)
$$\mathbf{e} = \frac{\mathbf{v} \times \mathbf{h}}{\mu} - \frac{\mathbf{r}}{\|\mathbf{r}\|}$$
- **Physical Meaning:** A conserved dimensionless vector pointing along the major axis directly toward **periapsis** (point of closest approach).
- Magnitude: $e = \|\mathbf{e}\|$.

### 3.4 Specific Orbital Energy ($\epsilon$)
$$\epsilon = \frac{v^2}{2} - \frac{\mu}{r}$$
- For bound elliptical trajectories: $\epsilon < 0$.
- Relationship to semi-major axis:
  $$a = -\frac{\mu}{2\epsilon}$$

---

## 4. Mathematical Derivation: State $\to$ Elements

Given $\mathbf{r}_{ECI}, \mathbf{v}_{ECI}$, and central gravitational parameter $\mu$:

1. **Magnitudes:**
   $$r = \|\mathbf{r}\|, \quad v = \|\mathbf{v}\|$$

2. **Angular Momentum & Energy:**
   $$\mathbf{h} = \mathbf{r} \times \mathbf{v}, \quad h = \|\mathbf{h}\|, \quad \hat{\mathbf{h}} = \mathbf{h} / h$$
   $$\epsilon = \frac{v^2}{2} - \frac{\mu}{r}$$

3. **Semi-Major Axis ($a$):**
   $$a = -\frac{\mu}{2\epsilon}$$

4. **Eccentricity ($e$):**
   $$\mathbf{e} = \frac{\mathbf{v} \times \mathbf{h}}{\mu} - \frac{\mathbf{r}}{r}, \quad e = \|\mathbf{e}\|$$

5. **Inclination ($i$):**
   $$i = \arccos\left(\operatorname{clamp}\left(\frac{h_z}{h}, -1, 1\right)\right) \in [0, \pi]$$

6. **Ascending Node Vector ($\mathbf{n}$):**
   $$\mathbf{n} = [-h_y, h_x, 0]^T, \quad n = \sqrt{h_x^2 + h_y^2}$$

7. **Right Ascension of Ascending Node ($\Omega$):**
   - For inclined orbits ($n > 10^{-11}$):
     $$\Omega = \operatorname{atan2}(n_y, n_x) = \operatorname{atan2}(h_x, -h_y) \pmod{2\pi}$$
   - For equatorial orbits ($n \le 10^{-11}$): $\Omega \equiv 0.0$ rad (convention).

8. **Argument of Periapsis ($\omega$):**
   - For inclined eccentric orbits ($n > 10^{-11}, e > 10^{-11}$):
     $$\cos\omega = \mathbf{n} \cdot \mathbf{e}, \quad \sin\omega = (\mathbf{n} \times \mathbf{e}) \cdot \hat{\mathbf{h}}$$
     $$\omega = \operatorname{atan2}(\sin\omega, \cos\omega) \pmod{2\pi}$$
   - For equatorial eccentric orbits ($n \le 10^{-11}, e > 10^{-11}$):
     $$\omega = \varpi = \operatorname{atan2}(e_y, e_x) \pmod{2\pi}$$
   - For circular orbits ($e \le 10^{-11}$): $\omega \equiv 0.0$ rad (convention).

9. **True Anomaly ($\nu$):**
   - For eccentric orbits ($e > 10^{-11}$):
     $$\cos\nu = \mathbf{e} \cdot \mathbf{r}, \quad \sin\nu = (\mathbf{e} \times \mathbf{r}) \cdot \hat{\mathbf{h}}$$
     $$\nu = \operatorname{atan2}(\sin\nu, \cos\nu) \pmod{2\pi}$$
   - For circular inclined orbits ($e \le 10^{-11}, n > 10^{-11}$):
     $$\nu = u = \operatorname{atan2}((\mathbf{n} \times \mathbf{r}) \cdot \hat{\mathbf{h}}, \mathbf{n} \cdot \mathbf{r}) \pmod{2\pi}$$
   - For circular equatorial orbits:
     $$\nu = \lambda = \operatorname{atan2}(r_y, r_x) \pmod{2\pi}$$

---

## 5. Mathematical Derivation: Elements $\to$ State

To convert $(a, e, i, \Omega, \omega, \nu)$ back into Cartesian position and velocity in ECI, we construct the state in the **Perifocal Coordinate Frame ($PQW$)** and rotate into ECI.

### 5.1 The Perifocal Frame ($PQW$)
- **$\mathbf{P}$ Axis:** Unit vector pointing toward periapsis in the orbital plane.
- **$\mathbf{W}$ Axis:** Unit vector along orbital angular momentum $\hat{\mathbf{h}}$ (orbit normal).
- **$\mathbf{Q}$ Axis:** Completes the right-handed triad in the orbital plane ($\mathbf{Q} = \mathbf{W} \times \mathbf{P}$, in the direction of motion at $\nu = 90^\circ$).

### 5.2 Perifocal Coordinates
1. **Semi-Latus Rectum:**
   $$p = a(1 - e^2)$$
2. **Radial Distance:**
   $$r = \frac{p}{1 + e \cos\nu}$$
3. **Position in Perifocal Frame:**
   $$\mathbf{r}_{PQW} = \begin{bmatrix} r \cos\nu \\ r \sin\nu \\ 0 \end{bmatrix}$$
4. **Velocity in Perifocal Frame:**
   $$\mathbf{v}_{PQW} = \sqrt{\frac{\mu}{p}} \begin{bmatrix} -\sin\nu \\ e + \cos\nu \\ 0 \end{bmatrix}$$

### 5.3 Rotation from Perifocal to ECI ($C_{ECI\_PQW}$)
The orientation of the perifocal frame relative to ECI is defined by the 3-1-3 Euler rotation sequence:
$$C_{PQW\_ECI} = R_3(\omega) R_1(i) R_3(\Omega)$$
$$C_{ECI\_PQW} = C_{PQW\_ECI}^T = R_3(-\Omega) R_1(-i) R_3(-\omega)$$

The columns of $C_{ECI\_PQW}$ are the unit vectors $\mathbf{P}, \mathbf{Q}, \mathbf{W}$ expressed in ECI:
$$\mathbf{P} = \begin{bmatrix} \cos\Omega \cos\omega - \sin\Omega \sin\omega \cos i \\ \sin\Omega \cos\omega + \cos\Omega \sin\omega \cos i \\ \sin\omega \sin i \end{bmatrix}$$
$$\mathbf{Q} = \begin{bmatrix} -\cos\Omega \sin\omega - \sin\Omega \cos\omega \cos i \\ -\sin\Omega \sin\omega + \cos\Omega \cos\omega \cos i \\ \cos\omega \sin i \end{bmatrix}$$
$$\mathbf{W} = \begin{bmatrix} \sin\Omega \sin i \\ -\cos\Omega \sin i \\ \cos i \end{bmatrix}$$

Finally:
$$\mathbf{r}_{ECI} = (r \cos\nu) \mathbf{P} + (r \sin\nu) \mathbf{Q}$$
$$\mathbf{v}_{ECI} = \left(\sqrt{\frac{\mu}{p}} (-\sin\nu)\right) \mathbf{P} + \left(\sqrt{\frac{\mu}{p}} (e + \cos\nu)\right) \mathbf{Q}$$

---

## 6. Critical Concept: Why Six Numbers Do Not Always Mean Six Independent Angles

A central lesson in aerospace astrodynamics is that **classical orbital elements possess coordinate singularities**:

### 6.1 Circular Orbit Singularity ($e = 0$)
When an orbit is circular, the orbit is rotationally symmetric. There is no periapsis (closest approach point) or apoapsis (furthest approach point).
- The eccentricity vector is zero ($\mathbf{e} = \mathbf{0}$).
- The argument of periapsis ($\omega$) is **physically undefined**.
- True anomaly ($\nu$), which is measured from periapsis, is also physically undefined.
- **Convention:** We set $\omega = 0$ and measure angular position as **argument of latitude** ($u = \omega + \nu$) from the ascending node.

### 6.2 Equatorial Orbit Singularity ($i = 0$ or $i = \pi$)
When an orbit lies entirely in the equatorial plane ($XY_{ECI}$):
- The orbital plane does not intersect the equatorial plane along a distinct line; the planes are parallel or coincident.
- The node vector is zero ($\mathbf{n} = \mathbf{0}$).
- The Right Ascension of Ascending Node ($\Omega$) is **physically undefined**.
- **Convention:** We set $\Omega = 0$ and measure periapsis orientation as **longitude of periapsis** ($\varpi = \Omega + \omega$) from $+X_{ECI}$.

### 6.3 Circular Equatorial Singularity ($e = 0$ and $i = 0$)
The most degenerate case: both $\Omega$ and $\omega$ are undefined.
- **Convention:** $\Omega = 0, \omega = 0$, and position is represented by **true longitude** ($\lambda = \Omega + \omega + \nu = \operatorname{atan2}(r_y, r_x)$).

> [!IMPORTANT]
> Because of these singularities, attempting to numerically differentiate or filter classical orbital elements for near-circular or near-equatorial satellites can lead to severe numerical instabilities (divisions by zero). In later missions, alternative nonsingular representations (e.g. Equinoctial Elements or Cartesian coordinates) are used for GNC filtering.

---

## 7. Angle Normalization Policy

All angular calculations normalize angles to $[0, 2\pi)$ radians using [`normalize_angle_2pi_rad()`](file:///C:/Users/user/AstraDock/cpp/math/angle.hpp):

| Input Angle | Equivalent Radians | Normalized $[0, 2\pi)$ | Output Degrees |
| :---: | :---: | :---: | :---: |
| $-10^\circ$ | $-0.174533$ rad | $5.934119$ rad | $350^\circ$ |
| $350^\circ$ | $6.108652$ rad | $6.108652$ rad | $350^\circ$ |
| $370^\circ$ | $6.457718$ rad | $0.174533$ rad | $10^\circ$ |
| $0^\circ$ | $0.000000$ rad | $0.000000$ rad | $0^\circ$ |
| $360^\circ$ | $6.283185$ rad | $0.000000$ rad | $0^\circ$ |

---

## 8. Summary of What Was Learned

1. **State Transformation:** A 6-state Cartesian vector in ECI maps uniquely to 6 classical Keplerian elements ($a, e, i, \Omega, \omega, \nu$) for any bound non-circular non-equatorial orbit.
2. **Perifocal Frame:** The Perifocal frame ($PQW$) provides an orthogonal basis aligned with the orbit geometry, allowing simple closed-form computation of $\mathbf{r}$ and $\mathbf{v}$ before rotating into ECI via $C_{ECI\_PQW}$.
3. **Orbital Invariants:** In ideal two-body gravitational propagation, $a, e, i, \Omega, \omega$ are strictly conserved invariants of motion; only $\nu(t)$ advances along the orbit.
4. **Degeneracies:** When $e \to 0$ or $i \to 0$, classical elements lose uniqueness, requiring explicit convention-based representations.
