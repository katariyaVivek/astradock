# Lesson 011 — Spacecraft Environment & Force/Torque Models

## 1. Physical Motivation

In Milestones M01 through M10, AstraDock simulated spacecraft motion under the idealized assumptions of:
1. Spherically symmetric point-mass Earth gravity (pure Keplerian two-body motion),
2. A perfect hard vacuum with zero aerodynamic drag,
3. No gravitational influence from external celestial bodies (such as the Moon or Sun),
4. Point-mass attitude decoupling (zero gravitational torque on the spacecraft body).

While these assumptions establish foundational dynamics, a spacecraft operating in real near-Earth orbits experiences continuous environmental forces and torques:
- **Earth Oblateness ($J_2$)**: Earth is not a sphere; its equatorial diameter exceeds its polar diameter by approximately $42.8\,\text{km}$ due to rotation. This equatorial mass bulge exerts non-spherical gravitational forces that rotate the orbital plane (nodal precession) and rotate the ellipse within its plane (apsidal precession).
- **Atmospheric Drag**: In Low Earth Orbit (LEO, below $1000\,\text{km}$), neutral atmospheric gases collide with the spacecraft, opposing its motion relative to the rotating atmosphere, dissipating orbital mechanical energy, and causing altitude decay.
- **Third-Body Gravitational Tidal Forces**: In high orbits (such as Geostationary Orbit, GEO), the gravitational pull of the Moon and Sun creates differential (tidal) accelerations relative to Earth's center of mass.
- **Gravity-Gradient Torques**: Because gravity weakens with distance as $1/r^2$, the side of a finite spacecraft closer to Earth experiences slightly stronger gravity than the side farther away. This differential gravitational force exerts a restoring torque that aligns the spacecraft's axis of minimum moment of inertia with the local nadir/zenith direction.

Milestone M11 implements these four foundational environmental models in AstraDock.

---

## 2. Mathematical Models & Governing Equations

### 2.1 Earth Oblateness ($J_2$) Perturbation

Earth's gravitational potential expressed in spherical harmonics is:
$$U(r, \phi) = \frac{\mu}{r} \left[ 1 - J_2 \left(\frac{R_E}{r}\right)^2 \frac{3\sin^2\phi - 1}{2} + \dots \right]$$

Taking the gradient $\mathbf{a}_{J2} = \nabla U_{J2}$ in Earth-Centered Inertial (ECI) coordinates $\mathbf{r} = [x, y, z]^T$:

$$\mathbf{a}_{J2}(\mathbf{r}) = -\frac{3 \mu J_2 R_E^2}{2 r^5} \begin{bmatrix} x \left( 1 - 5 \dfrac{z^2}{r^2} \right) \\ y \left( 1 - 5 \dfrac{z^2}{r^2} \right) \\ z \left( 3 - 5 \dfrac{z^2}{r^2} \right) \end{bmatrix}$$

Where:
- $\mu = 3.986004418 \times 10^{14}\,\text{m}^3/\text{s}^2$ (Earth gravitational parameter),
- $R_E = 6378137.0\,\text{m}$ (Earth equatorial reference radius, WGS 84),
- $J_2 = 1.08262668 \times 10^{-3}$ (second zonal harmonic coefficient).

#### Secular Orbital Precession
First-order perturbation theory shows that $J_2$ causes zero secular change in semi-major axis $a$, eccentricity $e$, or inclination $i$. However, it drives steady secular rates in the orientation angles:
1. **Nodal Regression Rate** (Right Ascension of the Ascending Node, $\Omega$):
   $$\dot{\Omega} = -\frac{3}{2} n J_2 \left(\frac{R_E}{p}\right)^2 \cos(i)$$
   For prograde orbits ($i < 90^\circ$), $\cos(i) > 0$, so $\dot{\Omega} < 0$ (the orbital plane regresses westward).
2. **Apsidal Precession Rate** (Argument of Periapsis, $\omega$):
   $$\dot{\omega} = \frac{3}{4} n J_2 \left(\frac{R_E}{p}\right)^2 \left( 5\cos^2(i) - 1 \right)$$
   At the **critical inclination** $i = \arccos(1/\sqrt{5}) \approx 63.4349^\circ$ (and $116.565^\circ$), $5\cos^2(i) - 1 = 0$, so apsidal precession halts ($\dot{\omega} = 0$, the basis of Molniya orbits).

---

### 2.2 Atmospheric Drag & Earth Rotation

#### Atmospheric Velocity Relative to Spacecraft
Earth's atmosphere co-rotates with the planet about the inertial $+Z$ axis with angular rate vector $\boldsymbol{\omega}_E = [0, 0, \omega_E]^T$ where $\omega_E = 7.2921150 \times 10^{-5}\,\text{rad/s}$.

The atmospheric velocity in the ECI frame at position $\mathbf{r}_{\text{ECI}}$ is:
$$\mathbf{v}_{\text{atm}} = \boldsymbol{\omega}_E \times \mathbf{r}_{\text{ECI}} = \begin{bmatrix} -\omega_E y \\ \omega_E x \\ 0 \end{bmatrix}$$

The spacecraft's velocity relative to the atmosphere is:
$$\mathbf{v}_{\text{rel}} = \mathbf{v}_{\text{ECI}} - \mathbf{v}_{\text{atm}} = \begin{bmatrix} v_x + \omega_E y \\ v_y - \omega_E x \\ v_z \end{bmatrix}$$

#### Exponential Density Model
The neutral density $\rho(h)$ at geometric altitude $h = \|\mathbf{r}_{\text{ECI}}\| - R_E$ is parameterized by a reference base altitude $h_0$, reference density $\rho_0$, and scale height $H$:
$$\rho(h) = \rho_0 \exp\left( -\frac{h - h_0}{H} \right)$$

#### Drag Acceleration Vector
The aerodynamic drag force opposes $\mathbf{v}_{\text{rel}}$:
$$\mathbf{a}_{\text{drag}} = -\frac{1}{2} C_D \frac{A}{m} \rho(h) \|\mathbf{v}_{\text{rel}}\| \mathbf{v}_{\text{rel}}$$

Where:
- $C_D$ is the dimensionless drag coefficient (typically $2.0 - 2.2$ for spacecraft in free molecular flow),
- $A$ is the cross-sectional reference area ($\text{m}^2$),
- $m$ is the spacecraft mass ($\text{kg}$).

Because $\mathbf{a}_{\text{drag}} \cdot \mathbf{v}_{\text{rel}} < 0$, drag extracts mechanical energy, decreasing specific orbital energy $\varepsilon = \frac{v^2}{2} - \frac{\mu}{r}$ and causing secular semi-major axis shrinkage $\frac{da}{dt} < 0$.

---

### 2.3 Third-Body Tidal Gravitational Perturbation

When a spacecraft orbits Earth in the presence of a third body (such as the Moon with position $\mathbf{r}_3$ and parameter $\mu_3$), the net acceleration relative to Earth's center of mass consists of two parts:
1. **Direct gravitational attraction** of the third body on the spacecraft:
   $$\mathbf{a}_{\text{direct}} = \mu_3 \frac{\mathbf{r}_3 - \mathbf{r}}{\|\mathbf{r}_3 - \mathbf{r}\|^3}$$
2. **Indirect acceleration** of Earth's center of mass toward the third body:
   $$\mathbf{a}_{\text{indirect}} = \mu_3 \frac{\mathbf{r}_3}{\|\mathbf{r}_3\|^3}$$

The net tidal perturbation in ECI is the difference:
$$\mathbf{a}_{3B} = \mathbf{a}_{\text{direct}} - \mathbf{a}_{\text{indirect}} = \mu_3 \left( \frac{\mathbf{r}_3 - \mathbf{r}}{\|\mathbf{r}_3 - \mathbf{r}\|^3} - \frac{\mathbf{r}_3}{\|\mathbf{r}_3\|^3} \right)$$

Along the Earth–Moon line of centers ($\|\mathbf{r}\| \ll \|\mathbf{r}_3\|$), this produces the classic tidal dipole:
$$a_{\text{tidal}} \approx \frac{2 \mu_3 r}{r_3^3}$$

---

### 2.4 Gravity-Gradient Attitude Torque

In a non-uniform gravitational field, a rigid body of non-spherical mass distribution experiences an attitude torque about its center of mass.

Let $\hat{\mathbf{r}}_{\mathcal{I}} = \mathbf{r}_{\text{ECI}} / \|\mathbf{r}_{\text{ECI}}\|$ be the unit radial vector from Earth center to the spacecraft in ECI.

1. Transform the radial direction into the **Spacecraft Body Frame** ($\mathcal{B}$):
   $$\hat{\mathbf{r}}_{\mathcal{B}} = q^* \otimes \hat{\mathbf{r}}_{\mathcal{I}} \otimes q = \begin{bmatrix} u_x \\ u_y \\ u_z \end{bmatrix}$$
2. The gravity-gradient torque expressed in Body coordinates is:
   $$\boldsymbol{\tau}_{gg} = \frac{3\mu}{r^3} \hat{\mathbf{r}}_{\mathcal{B}} \times (\mathbf{I} \hat{\mathbf{r}}_{\mathcal{B}}) = \frac{3\mu}{r^3} \begin{bmatrix} (I_{zz} - I_{yy}) u_y u_z \\ (I_{xx} - I_{zz}) u_x u_z \\ (I_{yy} - I_{xx}) u_x u_y \end{bmatrix}$$

#### Physical Properties:
- **Spherical Symmetries**: If $I_{xx} = I_{yy} = I_{zz}$, $\boldsymbol{\tau}_{gg} \equiv \mathbf{0}$.
- **Principal Axis Alignment**: If $\hat{\mathbf{r}}_{\mathcal{B}}$ aligns with a principal axis (e.g. $[1,0,0]^T$), $\boldsymbol{\tau}_{gg} = \mathbf{0}$ (equilibrium attitude).
- **Pitch Libration**: If a satellite with $I_{zz} > I_{xx}$ (a dumbbell shape) is perturbed in pitch by angle $\theta$, gravity-gradient torque acts as a restoring spring ($\tau_y \propto -\sin(2\theta)$), causing pitch oscillations (librations) with frequency $\omega_{\text{lib}} = \sqrt{3(I_{zz} - I_{xx})/I_{yy}} \cdot \omega_{\text{orb}}$.

---

## 3. Coordinate Frames & Units

AstraDock maintains strict frame discipline across all environmental models:

| Quantity | Symbol | Coordinate Frame | Units |
| :--- | :---: | :--- | :--- |
| Earth Position | $\mathbf{r}$ | **ECI** ($\mathcal{I}$) | $\text{m}$ |
| Spacecraft Velocity | $\mathbf{v}$ | **ECI** ($\mathcal{I}$) | $\text{m/s}$ |
| $J_2$ Perturbation Acceleration | $\mathbf{a}_{J2}$ | **ECI** ($\mathcal{I}$) | $\text{m/s}^2$ |
| Atmospheric Drag Acceleration | $\mathbf{a}_{\text{drag}}$ | **ECI** ($\mathcal{I}$) | $\text{m/s}^2$ |
| Third-Body Tidal Acceleration | $\mathbf{a}_{3B}$ | **ECI** ($\mathcal{I}$) | $\text{m/s}^2$ |
| Gravity-Gradient Torque | $\boldsymbol{\tau}_{gg}$ | **BODY** ($\mathcal{B}$) | $\text{N}\cdot\text{m}$ |

---

## 4. Software Architecture & Configuration Switches

The environmental subsystem is structured as pure, stateless physics functions in `cpp/environment/`:
- `j2_gravity.hpp`: `j2_acceleration_eci()`
- `atmospheric_drag.hpp`: `relative_atmospheric_velocity_eci()`, `exponential_atmospheric_density()`, `drag_acceleration_eci()`
- `third_body_gravity.hpp`: `third_body_acceleration_eci()`
- `gravity_gradient.hpp`: `gravity_gradient_torque_body()`
- `environment_models.hpp`: `EnvironmentConfiguration`, `EnvironmentalParameters`, `EnvironmentalEffects`, `spacecraft_environmental_derivative()`, `propagate_spacecraft_environmental()`

### Zero-Overhead Regression Guarantee
`EnvironmentConfiguration` contains four booleans:
```cpp
struct EnvironmentConfiguration {
    bool enable_j2{false};
    bool enable_drag{false};
    bool enable_third_body{false};
    bool enable_gravity_gradient{false};
};
```
When all flags are `false`, `spacecraft_environmental_derivative` evaluates zero additional forces or torques and reproduces M10 6-DOF propagation bitwise identically.

---

## 5. Verification Summary

The M11 environment models have been verified through:
1. **Analytical Benchmarks**: Special cases (equatorial plane, polar axis, critical latitude, spherical inertia nulls, exact line-of-centers tidal acceleration).
2. **Orbital Precession**: Multi-orbit simulation reproducing secular RAAN regression and argument of periapsis advancement under $J_2$.
3. **Decay Dissipation**: Strict monotonicity of mechanical energy dissipation and altitude decay under atmospheric drag.
4. **Attitude Libration**: Gravity-gradient pitch restoring torque producing stable oscillatory dynamics.
5. **Deterministic Regression**: 100% bitwise match with M10 when disabled, and deterministic multi-run repeatability.
6. **Independent Python Oracle**: `python/audit/independent_environment_reference.py` verifying C++ telemetry with discrepancies $< 10^{-17}\,\text{m/s}^2$.
