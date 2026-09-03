# M16 Guidance & Control Foundations — Validation Report

## Scope

`cpp/control/attitude_pd.hpp` (error vector, geodesic angle, PD, rate damping,
gain suggestion), `cpp/control/relative_pd.hpp` (CW feedforward + relative PD),
`cpp/control/double_integrator_lqr.hpp` (closed-form LQR, eigenvalues, damping),
`cpp/control/control_saturation.hpp` (desired/achieved/deficit), tests
(`tests/cpp/test_control.cpp`, 11 cases / 65 assertions), demo
(`tools/control_demo.cpp`), oracle
(`python/audit/independent_control_reference.py`), telemetry (`data/m16_*.csv`),
figures (`artifacts/figures/m16_*`), lesson (`docs/lessons/016_guidance_control.md`).

## Analytical checks (all pass)

- Error vector: 10° X-rotation → `[sin5°, 0, 0]`; sign per axis; identity → 0.
- Double cover: `q` vs `−q` identical error; geodesic angle exact.
- PD linearity in both errors; restoring direction (positive error → negative
  torque); zero-error zero-torque.
- Rate damping principal-axis values; energy removal (`tau·w < 0` per axis).
- Feedforward `[3n²x + 2nvy, −2nvx, −n²z]` to 1e-12; exact setpoint cancellation.
- LQR: hand values, both ARE residuals to 1e-12, Hurwitz eigenvalues, unit-weight
  `K = [1, √3]`, damping `√3/2`; action opposes displacement.
- Saturation: per-axis clamp, deficit identity, symmetric negative side.

## Independent verification

Pure-Python oracle reimplements all five laws; seven spot cross-checks pass. No
production code imported.

## Scenario validation (metrics)

- Detumble of master-spec tumble `[0.2, −0.1, 0.15]` rad/s on a (10, 20, 30)
  bus with ζ = 0.9 / ts = 30 s tuning: settles (< 0.57° and < 0.29°/s) inside
  60 s of a 100 s window; terminal error 0.00031°, terminal rate 6.3e-05°/s;
  peak torque 0.87 Nm under the 0.5 Nm-per-axis demo cap with early saturation
  flagged (documents the M14-sizing lesson: this torque class needs thrusters
  or larger wheels — allocation care for M17).
- 90° yaw slew with the same law: terminal error 0.00016°.
- 100 m radial station keeping (CW plant, Kp = 1e-5, Kd = 6e-3): 0.0044 m and
  sub-0.05 m/s residual after 1 h; peak command 1.36 mm/s². (Initial gains
  2e-6/3e-3 closed only to 5.4 m — retuned on evidence, documented here rather
  than hidden.)

## Regression

New `cpp/control/` + tests/demo only; no existing-file behavior changed. Full
suite 256/256.

## Limitations

No integral action; no closed loop (M17); no thruster allocation (M17); CW
feedforward inherits CW validity limits; gains are bus-specific.
