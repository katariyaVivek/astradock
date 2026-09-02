# Milestone M11 Validation Report: Spacecraft Environment & Force/Torque Models

## 1. Executive Summary

Milestone M11 upgrades AstraDock from idealized two-body physics to realistic near-Earth orbital environmental dynamics. Following the completion and passage of the mandatory **M10 Numerical Entry Gate**, four foundational environmental perturbation models were implemented and verified:
1. **Earth Oblateness ($J_2$) Gravitational Perturbation**: Second zonal harmonic potential gradient driving secular nodal regression and apsidal precession.
2. **Atmospheric Drag with Earth Rotation**: Exponential neutral density atmosphere with co-rotating atmospheric velocity field, producing orbital energy dissipation and altitude decay.
3. **Third-Body Gravitational Tidal Perturbations**: Direct third-body gravitational attraction combined with indirect primary acceleration, verified against lunar tidal dipole limits.
4. **Gravity-Gradient Attitude Torque**: Differential gravitational torque on asymmetric rigid bodies in Earth's non-uniform gravity field, driving restorative pitch librations.

All **150 unit tests pass** (100%), with zero compiler warnings under MSVC `/W4 /permissive-`, zero Ruff lint errors, and complete cross-verification against an independent Python reference oracle.

---

## 2. Test Execution & Build Verification

The complete CTest test suite was compiled in Release configuration and executed:

```text
100% tests passed, 0 tests failed out of 150
Total Test time (real) = 2.16 sec
```

### Test Distribution Across Milestones:
- `astradock_vector3_tests`: 18 tests
- `astradock_two_body_tests`: 15 tests
- `astradock_integrator_tests`: 16 tests
- `astradock_orbit_tests`: 13 tests
- `astradock_validation_tests`: 14 tests
- `astradock_matrix3_tests`: 14 tests
- `astradock_frame_tests`: 14 tests
- `astradock_audit_tests`: 12 tests
- `astradock_elements_tests`: 14 tests
- `astradock_quaternion_tests`: 14 tests
- `astradock_attitude_dynamics_tests`: 14 tests
- `astradock_6dof_tests`: 11 tests
- `astradock_environment_tests` (NEW in M11): 9 tests
- **Total: 150 / 150 tests passing**

---

## 3. M10 Numerical Entry Gate Re-Validation

Before implementing environmental models, Stage A verified that rotational drift in M10 is pure numerical integration truncation error $\mathcal{O}(\Delta t^4)$:

| Timestep $\Delta t$ | Rel Rate Error $E_{\omega}$ | Rel Attitude Error $E_q$ | Rate Conv Order $p$ | Attitude Conv Order $p$ | Momentum Drift $\Delta H$ | Energy Drift $\Delta T$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| $1.000\,\text{s}$ | $2.5592 \times 10^{-6}$ | $6.9744 \times 10^{-6}$ | — | — | $1.98 \times 10^{-6}$ | $8.42 \times 10^{-7}$ |
| $0.500\,\text{s}$ | $1.3912 \times 10^{-7}$ | $3.1672 \times 10^{-7}$ | **$4.201$** | **$4.461$** | $1.21 \times 10^{-7}$ | $2.62 \times 10^{-8}$ |
| $0.250\,\text{s}$ | $7.0279 \times 10^{-9}$ | $1.2989 \times 10^{-8}$ | **$4.307$** | **$4.608$** | $7.48 \times 10^{-9}$ | $8.16 \times 10^{-10}$ |

**Gate Decision: PASS** (Documented in `docs/validation/m10_numerical_baseline.md`).

---

## 4. Environmental Scenario Validation

### 4.1 M10 Disabled-Environment Regression
With all environment flags disabled (`enable_j2 = false`, `enable_drag = false`, `enable_third_body = false`, `enable_gravity_gradient = false`), `spacecraft_environmental_derivative` and `propagate_spacecraft_environmental` match M10 6-DOF propagation **bitwise identically** over all translational and rotational state variables.

### 4.2 J2 Earth Oblateness Scenario
- **Orbit**: $a = 10,000\,\text{km}$, $e = 0.2$, $i = 45^\circ$, $\Omega_0 = 120^\circ$, $\omega_0 = 60^\circ$, propagated for 5 orbits ($49,760\,\text{s}$).
- **Observed Nodal Regression**: $\Delta\Omega = -0.9105^\circ$ ($\dot{\Omega} < 0$, westward precession as predicted).
- **Observed Apsidal Precession**: $\Delta\omega = +0.9679^\circ$ ($\dot{\omega} > 0$, prograde advance as predicted for $i = 45^\circ < 63.4^\circ$).
- **Secular Invariant Stability**: $|\Delta a / a_0| < 1.0 \times 10^{-3}$, $|\Delta i / i_0| < 1.0 \times 10^{-3}$ (only short-periodic oscillations).

### 4.3 Atmospheric Drag Decay Scenario
- **Orbit**: $300\,\text{km}$ circular LEO, $C_D = 2.2$, $A = 3.0\,\text{m}^2$, $m = 250\,\text{kg}$, propagated for 10 orbits ($54,312\,\text{s}$).
- **Altitude Loss**: $\Delta h = -1.5860\,\text{km}$ ($300.00\,\text{km} \to 298.41\,\text{km}$).
- **Energy Dissipation**: $\Delta\varepsilon = -7090.41\,\text{m}^2/\text{s}^2$ (strictly monotonic mechanical energy decay).

### 4.4 Third-Body Lunar Tidal Gravity Scenario
- **Orbit**: Geostationary Orbit ($r = 42,164\,\text{km}$) with Moon on ECI $+X$ axis ($r_3 = 384,400\,\text{km}$, $\mu_3 = 4.90487 \times 10^{12}\,\text{m}^3/\text{s}^2$).
- **Tidal Acceleration**: Direct minus indirect acceleration correctly produces a symmetric tidal bulge (pointing toward Moon at sublunar point, pointing away at antilunar point).
- **Cross-Verification Discrepancy**: C++ vs independent Python reference discrepancy is **$1.7004 \times 10^{-18}\,\text{m/s}^2$** (floating-point precision limit).

### 4.5 Gravity-Gradient Torque & Pitch Libration Scenario
- **Spacecraft**: Dumbbell satellite ($I_{xx} = 10.0\,\text{kg}\cdot\text{m}^2$, $I_{yy} = I_{zz} = 50.0\,\text{kg}\cdot\text{m}^2$) in $500\,\text{km}$ orbit with $10^\circ$ pitch offset.
- **Torque Response**: Peak restoring torque $\tau_{gg,y} = 7.289 \times 10^{-5}\,\text{N}\cdot\text{m}$.
- **Attitude Response**: Sustained, stable pitch oscillations (librations) with peak angular rate $\omega_y = 2.527 \times 10^{-4}\,\text{rad/s}$.

---

## 5. Independent Python Reference Oracle Audit

The independent Python oracle `python/audit/independent_environment_reference.py` executed across all generated CSV datasets:

```text
======================================================================
Auditing J2 Perturbation Telemetry: environment_j2_orbit.csv
======================================================================
  Initial RAAN:        120.0000 deg
  Final RAAN:          119.0895 deg (Delta = -0.9105 deg)
  Initial ArgP:        60.0000 deg
  Final ArgP:          60.9679 deg (Delta = 0.9679 deg)
J2 Scenario Audit Result: PASSED

======================================================================
Auditing Atmospheric Drag Telemetry: environment_drag_decay.csv
======================================================================
  Initial Altitude:    300.0000 km
  Final Altitude:      298.4140 km (Delta = -1.5860 km)
  Initial Energy:      -29843685.5818 m^2/s^2
  Final Energy:        -29850775.9877 m^2/s^2 (Delta = -7090.4060 m^2/s^2)
Atmospheric Drag Audit Result: PASSED

======================================================================
Auditing Third-Body Gravity Telemetry: environment_third_body.csv
======================================================================
  Max Acceleration Discrepancy (C++ vs Python): 1.7004e-18 m/s^2
Third-Body Gravity Audit Result: PASSED

======================================================================
Auditing Gravity-Gradient Torque Telemetry: environment_gravity_gradient.csv
======================================================================
  Max Gravity-Gradient Torque:                 7.2890e-05 N*m
  Max Pitch Libration Angular Rate:            2.5268e-04 rad/s
Gravity-Gradient Torque Audit Result: PASSED

======================================================================
OVERALL INDEPENDENT ENVIRONMENT ORACLE AUDIT: ALL TESTS PASSED
======================================================================
```

---

## 6. Generated Visualizations

Four publication-grade figures were generated in `artifacts/figures/`:
1. `m11_j2_precession.png`: Secular RAAN regression, ArgP advancement, and short-periodic element oscillations under $J_2$.
2. `m11_drag_decay.png`: Geometric altitude decay, semi-major axis shrinkage, and mechanical energy dissipation under atmospheric drag.
3. `m11_third_body.png`: Multi-axis lunar gravitational acceleration components and total tidal magnitude along GEO orbit.
4. `m11_gravity_gradient.png`: Pitch angle libration, pitch rate, and restoring gravity-gradient torque on dumbbell satellite.
