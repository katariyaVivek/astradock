# M14 Actuator Dynamics & Modeling — Validation Report

## Scope

Reaction wheels (`cpp/actuators/reaction_wheel.hpp`), thrusters + propellant
(`cpp/actuators/thruster.hpp`), assembly (`cpp/actuators/actuator_assembly.hpp`),
tests (`tests/cpp/test_actuators.cpp`, 11 cases / 112 assertions), demo
(`tools/actuator_demo.cpp`), oracle (`python/audit/independent_actuator_reference.py`),
telemetry (`data/m14_*.csv`), figures (`artifacts/figures/m14_*.png`), lesson
(`docs/lessons/014_actuators.md`).

## Analytical checks (all pass)

- Wheel sign: +X motor torque -> `-X` bus torque; `H_w` along +axis.
- Exact integration: `I d(omega) = tau dt` closes to 1e-12.
- Lag: `1 - exp(-dt/tau)` to 1e-14; 200-step convergence to target within 1e-9.
- Saturation sizing: `I w_max / tau = 1e-3*100/0.05 = 2.0 s` buildup observed.
- Thruster hand case: `[0,.5,0]x[2,0,0] = [0,0,-1.0]` N*m to 1e-12.
- ECI rotation: 90-deg yaw maps body +X thrust into +Y ECI to 1e-12.
- Mass flow: `-T/(Isp g0)` to 1e-14; 10 s burn depletion exact to 1e-12.
- Assembly: wheel/thruster/momentum/mass aggregation verified component-wise.
- Momentum conservation (bus + wheel, torque-free): total inertial Z residual
  `-2.2e-16 Nms` after 5 s of 0.1 Nm wheel torque on the 6-DOF bus.

## Independent verification

Pure-Python oracle reimplements every scalar equation independently; six spot
cross-checks pass (speed, lag, saturate, force, torque, mass flow). No production
code imported.

## Regression

Zero actuator command produces identically zero force/torque/propellant and
unchanged mass; the `ForceTorqueInput{}` path through `rk4_step_spacecraft` is
untouched, so M10/M11 truth baselines are preserved by construction (no modified
files in their path).

## Limitations

No friction/cogging/thermal/power-bus/bearing models; ideal on/off thrusters
(no MIB, plume, valve lag, slosh); constant Isp; no desaturation controller
(that is GNC scope, M16+).
