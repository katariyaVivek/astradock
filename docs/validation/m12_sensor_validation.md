# Milestone M12 Validation Report — Sensor Simulation & Measurement Models

## 1. Executive Summary

Milestone **M12** introduces realistic, educational sensor simulation models to AstraDock. All sensors generate discrete measurement packets with explicit timestamps, coordinate frames, SI units, sampling schedules, constant biases, zero-mean Gaussian noise, and deterministic failure dropouts.

Verification was performed across three independent layers:
1. **Layer 1 (C++ Catch2 Unit & Integration Tests)**: 159 / 159 tests passing (100%).
2. **Layer 2 (Statistical Property & Timing Checks)**: Sample means and standard deviations verified against theoretical distributions over $N = 10,000$ samples.
3. **Layer 3 (Independent Python Reference Oracle)**: Telemetry verified against pure Python reference computations across all four sensors.

---

## 2. Sensor Verification Scorecard

| Sensor / Metric | Expected Configuration | Sampled / Measured Value | Verification Result |
| :--- | :--- | :--- | :---: |
| **IMU Gyroscope Bias** | $[0.0005, -0.0003, 0.0002]\,\text{rad/s}$ | $[0.000500, -0.000301, 0.000202]\,\text{rad/s}$ | **PASS** |
| **IMU Gyroscope Noise $\sigma$** | $0.001000\,\text{rad/s}$ | $[0.001000, 0.001000, 0.001001]\,\text{rad/s}$ | **PASS** |
| **IMU Accelerometer Bias** | $[0.0020, -0.0010, 0.0015]\,\text{m/s}^2$ | $[0.002007, -0.000997, 0.001501]\,\text{m/s}^2$ | **PASS** |
| **IMU Accelerometer Noise $\sigma$** | $0.005000\,\text{m/s}^2$ | $[0.004992, 0.005000, 0.005001]\,\text{m/s}^2$ | **PASS** |
| **Free-Fall Specific Force** | $f_{\text{ideal}} \equiv 0$ | $f_{\text{drag}} \approx 8.90 \times 10^{-7}\,\text{m/s}^2 \ll g$ | **PASS** |
| **GNSS Position Bias** | $[2.500, -1.500, 3.000]\,\text{m}$ | $[2.442, -1.525, 3.039]\,\text{m}$ | **PASS** |
| **GNSS Position Noise $\sigma$** | $3.000\,\text{m}$ | $[2.970, 3.029, 2.984]\,\text{m}$ | **PASS** |
| **GNSS Velocity Bias** | $[0.0200, -0.0100, 0.0150]\,\text{m/s}$ | $[0.0197, -0.0108, 0.0152]\,\text{m/s}$ | **PASS** |
| **GNSS Velocity Noise $\sigma$** | $0.0300\,\text{m/s}$ | $[0.0302, 0.0297, 0.0302]\,\text{m/s}$ | **PASS** |
| **Star Tracker Quat Norm** | $\|q\| \equiv 1.0$ | $\max |\|q\| - 1| = 2.22 \times 10^{-16}$ | **PASS** |
| **Star Tracker Error Mean** | $\approx 0.0005\,\text{rad}$ ($103''$) | $0.000527\,\text{rad}$ ($108.61''$) | **PASS** |
| **Range Sensor Bias** | $+0.5000\,\text{m}$ | $+0.5011\,\text{m}$ | **PASS** |
| **Range Sensor Noise $\sigma$** | $0.1000\,\text{m}$ | $0.0999\,\text{m}$ | **PASS** |
| **Deterministic Dropout Window** | $[2000, 2100]\,\text{s} \implies \text{valid} = \text{false}$ | Exactly flagged invalid, zero fabricated data | **PASS** |
| **Truth Non-Interference** | 100% Bitwise Identical | Zero perturbation on truth state | **PASS** |

---

## 3. Specific Force Orbital Physics Validation

During orbital free fall at $500\,\text{km}$ altitude, the gravitational acceleration is $g \approx 8.43\,\text{m/s}^2$.

The accelerometer simulation was validated against two canonical test cases:
1. **Pure Gravitational Orbit**: In the absence of non-gravitational acceleration, $\mathbf{a}_{\text{inertial}} = \mathbf{g} \implies \mathbf{f}_I = \mathbf{a}_{\text{inertial}} - \mathbf{g} \equiv \mathbf{0}$. The ideal specific force evaluated to identically $0.0\,\text{m/s}^2$.
2. **Atmospheric Drag Acceleration**: When atmospheric drag was enabled in M11, the drag force produced an aerodynamic deceleration of $\sim 8.90 \times 10^{-7}\,\text{m/s}^2$. The accelerometer correctly registered this non-gravitational contact acceleration transformed into the Spacecraft Body frame.

---

## 4. Generated Figures

The demonstration tool `tools/sensor_demo.cpp` generated comprehensive telemetry datasets in `data/`, which were rendered into publication-grade figures via `python/analysis/plot_sensors.py`:

1. `artifacts/figures/m12_imu_telemetry.png`: Gyroscope body angular rates and Accelerometer specific force demonstrating near-zero free-fall contact acceleration.
2. `artifacts/figures/m12_gnss_telemetry.png`: ECI Cartesian position and velocity truth vs noisy GNSS measurements and Gaussian error histograms.
3. `artifacts/figures/m12_star_tracker_telemetry.png`: Star tracker attitude orientation error angle in arcseconds and probability density distribution.
4. `artifacts/figures/m12_range_telemetry.png`: Line-of-sight relative range to target spacecraft and measurement residuals.
5. `artifacts/figures/m12_sensor_sampling_rates.png`: Multi-rate timeline demonstrating 100 Hz IMU, 10 Hz Star Tracker, 10 Hz Range, 1 Hz GNSS, and dropout intervals.

---

## 5. Scope Guard Compliance

In accordance with `AGENTS.md` and the M12 specification, the following were **strictly excluded**:
- Extended Kalman Filters (EKF) and Unscented Kalman Filters (UKF).
- Sensor fusion and state estimation algorithms.
- Feedback guidance, attitude control, and thruster firing logic.
- Rendezvous, docking, computer vision, and machine learning.

These components are deferred to Milestone M13 (State Estimation & EKF) and subsequent milestones.
