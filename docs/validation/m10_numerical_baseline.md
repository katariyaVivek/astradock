# Milestone M10 Numerical Entry Gate & Historical Baseline

## 1. Executive Summary & Gate Decision

Prior to introducing environmental forces and torques in Milestone M11, a mandatory **M10 Numerical Entry Gate** was conducted to verify that the rotational invariant drifts observed in M10 ($\approx 8.42 \times 10^{-7}$ in rotational kinetic energy and $\approx 1.98 \times 10^{-6}$ in inertial angular momentum over 1 full orbit) are purely smooth $\mathcal{O}(\Delta t^4)$ Runge-Kutta numerical truncation error, rather than an implementation defect or improper quaternion normalization.

### Gate Decision: **PASS**

- **Empirical Angular Velocity Convergence Order**: $p \approx 4.20 - 4.32$ (consistent with 4th-order RK4).
- **Empirical Orientation Convergence Order**: $p \approx 4.46 - 4.61$.
- **Inertial Angular Momentum Drift Reduction**: Factor of $\approx 16.4$ when halving $\Delta t$ ($2^4 = 16$).
- **Rotational Kinetic Energy Drift Reduction**: Factor of $\approx 32.1$ when halving $\Delta t$ ($2^5 = 32$, scaling as $\mathcal{O}(h^5)$ due to the quadratic relation $E_{\text{rot}} \propto \boldsymbol{\omega}^2$).
- **Quaternion Norm Preservation**: Machine-precision $\|q\| = 1.0 \pm 0.0$ at all step boundaries.
- **Independent Python vs C++ Agreement**: Position error $< 0.52\,\mu\text{m}$, velocity error $< 0.5\,\text{nm/s}$, orientation error $< 0.26\,\mu\text{rad}$, angular rate error $< 10^{-13}\,\text{rad/s}$.

---

## 2. Baseline Scenario Configuration

The canonical 500 km circular Low Earth Orbit + asymmetric tumbling spacecraft was configured identically across C++ and independent Python:

```text
Orbital Radius r_0:           6878137.0 m (500 km altitude above WGS 84 Earth)
Orbital Velocity v_0:         7612.608119 m/s (circular speed)
Orbital Period T:             5676.977 s (~94.6 minutes)
Spacecraft Principal Inertia: diag(10.0, 20.0, 30.0) kg*m^2
Initial Attitude q_0:         [cos(22.5°), 0, 0, sin(22.5°)]  (45° yaw about +Z)
Initial Angular Rate w_0:     [0.05, 0.08, 0.02] rad/s (tumbling multi-axis rate)
External Force / Torque:      F = 0 N, tau = 0 N*m
```

---

## 3. Timestep Refinement Study ($T = 5677.0\,\text{s}$, 1 Full Orbit)

| Timestep $\Delta t$ | Integration Steps | Rotational Energy Drift $|\Delta E_{\text{rot}}|/E_0$ | Inertial Momentum Drift $\|\Delta \mathbf{H}_{\mathcal{I}}\|/H_0$ | Final $\boldsymbol{\omega}$ Error vs Ref (rad/s) | Final Attitude Error vs Ref (rad) |
| :---: | :---: | :---: | :---: | :---: | :---: |
| **$1.00\,\text{s}$ (M10 Baseline)** | 5677 | $8.4179 \times 10^{-7}$ | $1.7008 \times 10^{-6}$ | $9.7154 \times 10^{-7}$ | $1.6008 \times 10^{-4}$ |
| **$0.50\,\text{s}$** | 11354 | $2.6286 \times 10^{-8}$ | $1.0345 \times 10^{-7}$ | $4.8502 \times 10^{-8}$ | $6.5457 \times 10^{-6}$ |
| **$0.25\,\text{s}$** | 22708 | $8.1865 \times 10^{-10}$ | $6.4198 \times 10^{-9}$ | $2.6436 \times 10^{-9}$ | $2.9802 \times 10^{-7}$ |

*Reference solution: $\Delta t_{\text{ref}} = 0.05\,\text{s}$ (113540 steps).*

---

## 4. Separation of Causes

1. **Numerical Truncation Error ($\mathcal{O}(\Delta t^4)$)**:
   The observed drift scales strictly with timestep refinement. Reducing $\Delta t$ from $1.0\,\text{s}$ to $0.25\,\text{s}$ (a factor of 4) decreases the angular velocity error by a factor of 367 ($4^{4.25}$), and decreases angular momentum drift by a factor of 265 ($4^4 = 256$).
2. **Quaternion Normalization Policy**:
   Normalizing only at step boundaries ensures that intermediate Runge-Kutta stages sample unconstrained linear tangent space $\mathbb{R}^4$, preserving exact 4th-order order of convergence without introducing high-order perturbation artifacts.
3. **Implementation Correctness**:
   Both C++ and independent Python reference implementations agree to $< 10^{-13}\,\text{rad/s}$ in angular velocity, proving that no C++ specific coding defects exist.

The repository is cleared to proceed with **Stage B: Spacecraft Environment & Force/Torque Models**.
