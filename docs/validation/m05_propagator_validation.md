# M05 Propagator Validation Report

> **Scope:** This report characterizes the numerical accuracy of the AstraDock
> two-body propagator. It measures how integration error depends on timestep
> and integrator choice, quantifies invariant preservation, and establishes
> verification baselines. It does not evaluate physical-model fidelity against
> a higher-fidelity perturbed Earth model (reserved for M28).

---

## 1. Test Setup and Reference Quantities

### 1.1 Physical Constants (WGS 84 / Standard Earth)

| Quantity | Symbol | Value | Unit |
| :--- | :--- | ---: | :--- |
| Earth Gravitational Parameter | $\mu$ | $3.986004418 \times 10^{14}$ | $\text{m}^3/\text{s}^2$ |
| Earth Reference Radius | $R_\oplus$ | $6,378,137.0$ | $\text{m}$ |
| Ratio of Circumference to Diameter | $\pi$ | $3.141592653589793$ | -- |

### 1.2 Analytical Reference Orbit (500 km Equatorial Circular)

All analytical reference quantities are derived strictly from the physical constants above:

| Quantity | Formula | Analytical Value | Unit |
| :--- | :--- | ---: | :--- |
| Orbital Altitude | $h$ | $500,000.0$ | $\text{m}$ |
| Orbital Radius | $r_0 = R_\oplus + h$ | $6,878,137.0$ | $\text{m}$ |
| Circular Velocity | $v_c = \sqrt{\mu / r_0}$ | $7,612.608173$ | $\text{m}/\text{s}$ |
| Orbital Period | $T = 2\pi \sqrt{r_0^3 / \mu}$ | $5,676.978029$ | $\text{s}$ |
| Mean Motion | $n = \sqrt{\mu / r_0^3}$ | $0.001106783$ | $\text{rad}/\text{s}$ |
| Specific Mechanical Energy | $\epsilon_0 = -\mu / (2 r_0)$ | $-28,975,901.600$ | $\text{m}^2/\text{s}^2$ |
| Specific Angular Momentum | $h_0 = r_0 v_c$ | $52,360,561,942.754$ | $\text{m}^2/\text{s}$ |
| Gravitational Acceleration | $g(r_0) = \mu / r_0^2$ | $8.425509$ | $\text{m}/\text{s}^2$ |

Initial state vector at $t = 0$:
$$\mathbf{r}_0 = \begin{bmatrix} 6,878,137.0 \\ 0.0 \\ 0.0 \end{bmatrix}\,\text{m}, \quad \mathbf{v}_0 = \begin{bmatrix} 0.0 \\ 7,612.608173 \\ 0.0 \end{bmatrix}\,\text{m}/\text{s}$$

### 1.3 Propagation Scenarios and Integrators

- **Integration Interval:** One analytical period ($T = 5,676.978029$ s).
- **Tested Timesteps:** $\Delta t \in \{40.0, 20.0, 10.0, 5.0, 2.5\}$ s.
- **High-Accuracy Reference:** Classical RK4 at $\Delta t = 0.25$ s.
- **Integrators Evaluated:**
  1. Forward Euler (1st-order explicit)
  2. Classical Runge--Kutta 4th-Order (RK4, 4th-order explicit)

Within each sweep, initial conditions, duration, equations of motion, and compiler flags are held strictly identical; only the timestep $\Delta t$ and integrator method vary.

---

## 2. Validation Diagnostics and Error Metrics

| Metric | Symbol | Definition | Physical Significance |
| :--- | :--- | :--- | :--- |
| **Position Closure Error** | $\|\Delta \mathbf{r}\|$ | $\|\mathbf{r}(T) - \mathbf{r}(0)\|$ | Distance between final state and initial state after exactly one analytical period $T$. |
| **Velocity Closure Error** | $\|\Delta \mathbf{v}\|$ | $\|\mathbf{v}(T) - \mathbf{v}(0)\|$ | Velocity difference after one analytical period. |
| **Max Radial Deviation** | $\delta r_\text{max}$ | $\max_t |\|\mathbf{r}(t)\| - r_0|$ | Maximum departure from the circular orbit radius over the entire trajectory. |
| **Max Relative Energy Drift** | $\delta \epsilon_\text{max}$ | $\max_t \frac{|\epsilon(t) - \epsilon_0|}{|\epsilon_0|}$ | Conservation of specific mechanical energy $\epsilon = \frac{v^2}{2} - \frac{\mu}{r}$. |
| **Max Relative $h$ Drift** | $\delta h_\text{max}$ | $\max_t \frac{|\|\mathbf{h}(t)\| - h_0|}{h_0}$ | Conservation of angular momentum magnitude $\mathbf{h} = \mathbf{r} \times \mathbf{v}$. |
| **Max $h$ Direction Drift** | $\delta \phi_h$ | $\max_t \arccos\left(\frac{\mathbf{h}(0)\cdot\mathbf{h}(t)}{\|\mathbf{h}(0)\|\|\mathbf{h}(t)\|}\right)$ | Precession / tilt of the orbital plane normal. |
| **Max Phase Error** | $\delta \theta_\text{max}$ | $\max_t |\theta_\text{unwrapped}(t) - n t|$ | Along-track phase error relative to uniform mean motion. |
| **Period Estimation Error** | $\delta T$ | $T_\text{estimated} - T_\text{analytical}$ | Difference between trajectory zero-crossing time and analytical period. |

### 2.1 Empirical Convergence Order Formulation

For any error metric $E$ evaluated at two step sizes $h_1 > h_2 > 0$:
$$p = \frac{\log(E_1 / E_2)}{\log(h_1 / h_2)}$$
When halving the step size ($h_1 / h_2 = 2$):
- First-order method ($p = 1$): error halves ($E_1 / E_2 \approx 2$).
- Fourth-order method ($p = 4$): error reduces by a factor of 16 ($E_1 / E_2 \approx 16$).

---

## 3. Measured Results and Convergence Analysis

### 3.1 RK4 Timestep Convergence Sweep

| $\Delta t$ (s) | $\|\Delta\mathbf{r}\|$ (m) | $\|\Delta\mathbf{v}\|$ (m/s) | $\delta r_\text{max}$ (m) | $\delta\epsilon_\text{max}$ | $\delta h_\text{max}$ | $\delta\phi_h$ (rad) | $\delta\theta_\text{max}$ (rad) | $\delta T$ (s) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 40.0 | 4.767891 | 5.268e-03 | 0.931107 | 2.963e-08 | 1.482e-08 | 0.0 | 6.926e-07 | -6.259e-04 |
| 20.0 | 0.267956 | 2.962e-04 | 0.054847 | 9.261e-10 | 4.630e-10 | 0.0 | 3.909e-08 | -3.519e-05 |
| 10.0 | 0.015806 | 1.748e-05 | 0.003326 | 2.896e-11 | 1.448e-11 | 0.0 | 2.329e-09 | -2.076e-06 |
| 5.0 | 0.000958 | 1.060e-06 | 0.000205 | 9.041e-13 | 4.521e-13 | 0.0 | 1.422e-10 | -1.259e-07 |
| 2.5 | 5.915e-05 | 6.547e-08 | 1.272e-05 | 3.536e-14 | 1.763e-14 | 0.0 | 8.805e-12 | -7.770e-09 |
| *0.25 (ref)* | *7.058e-08* | *2.034e-10* | *7.730e-08* | *1.144e-14* | *5.683e-15* | *0.0* | *9.193e-13* | *-9.095e-12* |

**Empirical Convergence Orders (RK4):**

| Timestep Pair | Position Error Ratio | Measured Order $p_\text{pos}$ | Energy Error Ratio | Measured Order $p_\text{energy}$ |
| :--- | :---: | :---: | :---: | :---: |
| 40 s $\to$ 20 s | 17.7936 | **4.1533** | 31.9953 | **4.9998** |
| 20 s $\to$ 10 s | 16.9532 | **4.0835** | 31.9803 | **4.9991** |
| 10 s $\to$ 5 s | 16.4922 | **4.0437** | 32.0311 | **5.0014** |
| 5 s $\to$ 2.5 s | 16.2021 | **4.0181** | 25.5709 | **4.6764** |

*Analysis:*
1. The measured convergence order for position error approaches $4.0$ monotonically from above ($4.15 \to 4.08 \to 4.04 \to 4.02$), proving true asymptotic $O(\Delta t^4)$ convergence.
2. Energy conservation behaves as $O(\Delta t^5)$ in the smooth circular regime before hitting the double-precision round-off plateau ($\sim 10^{-14}$) at $\Delta t = 2.5$ s and $\Delta t = 0.25$ s.

---

### 3.2 Forward Euler Timestep Sweep

| $\Delta t$ (s) | $\|\Delta\mathbf{r}\|$ (m) | $\|\Delta\mathbf{v}\|$ (m/s) | $\delta r_\text{max}$ (m) | $\delta\epsilon_\text{max}$ | $\delta h_\text{max}$ | $\delta\phi_h$ (rad) | $\delta\theta_\text{max}$ (rad) | $\delta T$ (s) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :--- |
| 40.0 | 12,486,800.2 | 9,786.42 | 3,541,546.5 | 0.284498 | 0.177492 | 0.0 | 1.571087 | *N/A (no crossing)* |
| 20.0 | 7,502,956.0 | 6,780.83 | 1,808,622.4 | 0.183166 | 0.105900 | 0.0 | 0.980892 | *N/A (no crossing)* |
| 10.0 | 4,152,340.1 | 4,086.81 | 919,840.2 | 0.109568 | 0.059691 | 0.0 | 0.560188 | *N/A (no crossing)* |
| 5.0 | 2,187,390.5 | 2,258.49 | 466,911.6 | 0.061160 | 0.032054 | 0.0 | 0.301798 | *N/A (no crossing)* |
| 2.5 | 1,122,634.9 | 1,189.37 | 235,948.3 | 0.032525 | 0.016670 | 0.0 | 0.157065 | *N/A (no crossing)* |

**Empirical Convergence Orders (Euler):**

| Timestep Pair | Position Error Ratio | Measured Order $p_\text{pos}$ | Energy Error Ratio | Measured Order $p_\text{energy}$ |
| :--- | :---: | :---: | :---: | :---: |
| 40 s $\to$ 20 s | 1.6642 | **0.7349** | 1.5532 | **0.6353** |
| 20 s $\to$ 10 s | 1.8069 | **0.8535** | 1.6717 | **0.7413** |
| 10 s $\to$ 5 s | 1.8983 | **0.9247** | 1.7915 | **0.8412** |
| 5 s $\to$ 2.5 s | 1.9484 | **0.9623** | 1.8804 | **0.9110** |

*Analysis:*
1. Euler converges at approximately **first order** ($p \approx 1$), approaching $1.0$ from below as the timestep is refined ($0.73 \to 0.85 \to 0.92 \to 0.96$).
2. At $\Delta t = 40$ s and $\Delta t = 20$ s, the trajectory spirals outward so violently (altitude increases by $> 3,500$ km) that the orbital period cannot be estimated within one analytical period.
3. Even at $\Delta t = 2.5$ s (over 2,270 integration steps), Euler accumulates over **$1,122$ km** of position closure error and $3.25\%$ energy error.

---

### 3.3 Richardson-Style Error Analysis

Comparing solutions at $\Delta t = 20$ s ($h$), $\Delta t = 10$ s ($h/2$), and $\Delta t = 5$ s ($h/4$):

$$\frac{\|\mathbf{r}_{20} - \mathbf{r}_{10}\|}{\|\mathbf{r}_{10} - \mathbf{r}_5\|} = \frac{0.252150\,\text{m}}{0.014847\,\text{m}} \approx 16.98$$

For a 4th-order method, theoretical Richardson error scaling is $2^4 = 16$. The measured ratio $16.98$ confirms asymptotic 4th-order behavior without relying on analytical truth.

---

## 4. Addressing Core Engineering Questions

### 4.1 Does RK4 demonstrate approximately fourth-order convergence?
**Yes.** Empirical convergence order for position error is measured between **$4.15$** and **$4.02$**, strictly converging toward $4.000$ as $\Delta t \to 0$.

### 4.2 Does Euler demonstrate approximately first-order convergence?
**Yes.** Empirical convergence order for Euler is measured between **$0.73$** (coarse, distorted regime) and **$0.96$** (asymptotic regime), converging toward $1.000$.

### 4.3 How does timestep affect orbital phase error?
Phase error scales strictly with the integrator's global truncation order. For RK4, reducing $\Delta t$ from $40$ s to $2.5$ s decreases phase error from $6.93 \times 10^{-7}$ rad to $8.81 \times 10^{-12}$ rad (a factor of $\sim 80,000$). For Euler, phase error remains catastrophic ($0.56$ rad at $10$ s, equivalent to being $32^\circ$ behind the expected position along the orbit).

### 4.4 How does timestep affect energy drift?
In two-body motion, mechanical energy is an exact invariant. Euler continuously injects spurious energy ($+11\%$ per orbit at $\Delta t = 10$ s), causing an unstable outward spiral. RK4 maintains energy drift below $2.9 \times 10^{-11}$ at $\Delta t = 10$ s, dropping to the machine-precision noise floor ($\sim 3.5 \times 10^{-14}$) at $\Delta t \le 2.5$ s.

### 4.5 What timestep is currently recommended for nominal simulation?
**$\Delta t = 10$ s** is the recommended baseline for 500 km LEO educational demonstrations:
- Position closure error: $1.58$ cm per orbit ($0.0158$ m).
- Specific energy drift: $2.9 \times 10^{-11}$.
- Specific angular momentum drift: $1.45 \times 10^{-11}$.
- Computational effort: only $568$ integration steps per orbit.
For high-precision guidance and proximity operations, $\Delta t = 1.0$ s to $2.5$ s yields sub-millimetre closure ($< 0.06$ mm).

### 4.6 Is the current propagator suitable as a foundation for future GNC work?
- **Software architecture:** **Yes.** Modular, deterministic, header-only C++20, SI units, no hidden global state.
- **Numerical integrator:** **Yes.** RK4 provides sub-centimetre precision at moderate step sizes. (Symplectic methods may be considered in M28 for multi-year propagation).
- **Physical fidelity:** **No (by design).** The idealized point-mass spherical Earth model does not yet include $J_2$ oblateness, atmospheric drag, third-body gravity, solar radiation pressure, or coordinate transformations (ECI/ECEF/LVLH).

---

## 5. Scope Guard Confirmation

This milestone strictly added **numerical validation and QA diagnostics**. None of the following were added:
- Perturbations ($J_2$, atmospheric drag, lunisolar gravity)
- Coordinate frame transformations (ECEF, LVLH, Hill)
- Orbital elements subsystem (Keplerian / equinoctial elements)
- Attitude dynamics, quaternions, rigid-body mechanics
- Sensors, state estimation (EKF), guidance, control, rendezvous, computer vision, or machine learning.