# Lesson 014 — Spacecraft Actuator Dynamics & Modeling

## 1. Physical Problem & Motivation

Through M13 the simulator commands *abstract* forces and torques: `ForceTorqueInput`
appears directly in the 6-DOF derivative. A real spacecraft cannot do this. It owns
a finite set of physical devices that convert stored resources (electrical energy,
propellant, angular momentum capacity) into force and torque subject to hard limits:

- **Reaction wheels** spin a flywheel with an electric motor. Accelerating the wheel
  one way applies an equal-and-opposite reaction torque to the bus. No propellant is
  spent, but the wheel has a maximum speed (momentum capacity): once saturated, it
  cannot produce torque in the saturated direction until desaturated (classically by
  thrusters or magnetorquers).
- **Thrusters** expel mass. An offset thruster produces both force
  (translation/orbit control) and torque `tau = r x F` (attitude control, wheel
  desaturation) at the cost of propellant. Mass depletion changes every later
  force-to-acceleration computation, so thrust modeling introduces a live mass state.

M14 inserts this layer between future control commands and truth dynamics. M14 does
**not** implement guidance or control — it only models what the hardware does with
a command.

## 2. State Variables, Inputs, Outputs

| Symbol | Frame | Units | Meaning |
|---|---|---|---|
| `omega_w` | wheel spin axis (signed scalar) | rad/s | flywheel speed |
| `tau_motor` | wheel spin axis (signed scalar) | N*m | motor torque actually achieved |
| `H_w = I_w omega_w` | BODY vector along spin axis | N*m*s | stored wheel momentum |
| `tau_body = -axis * tau_motor` | BODY | N*m | reaction torque felt by the bus |
| `F_B = dir * T` | BODY | N | thruster force on the spacecraft |
| `tau_B = r_B x F_B` | BODY | N*m | thruster torque about the CM |
| `F_I = q (*) F_B` | ECI | N | thruster force in inertial axes |
| `m_dot = -T/(Isp g0)` | scalar | kg/s | propellant mass flow |
| `m_wet` | scalar | kg | spacecraft wet mass (dry floor) |

Commands (`ActuatorCommand`): per-wheel signed motor torque, per-thruster thrust
magnitude. Outputs (`ActuatorOutput`): summed body/ECI force, summed body torque,
per-source torque split, total wheel momentum, propellant used, and one explicit
flag per limiting mechanism (torque/speed/authority/thrust/depletion). Nothing is
silently clipped: every limiter sets a flag the caller must handle.

## 3. Governing Equations

Wheel spin: `I_w d(omega_w)/dt = tau_achieved`; stored energy `E = I omega^2/2`;
motor lag `tau(t) = target(1 - exp(-t/tau_lag))` integrated exactly as
`a += (target - a)(1 - exp(-dt/tau_lag))`; bus reaction `tau_bus = -axis * tau`.
Total (bus + wheels) inertial angular momentum is conserved under internal wheel
torques — verified in demo Scenario 3 to `-2.2e-16 Nms` residual.

Thruster: `F_B = dir*T`, `tau_B = r_B x F_B`, `m_dot = -T/(Isp*g0)` with
`g0 = 9.80665 m/s^2` exact. Piecewise-constant steps integrate exactly.

## 4. Assumptions and Failure Modes

Rigid balanced wheels (no friction/cogging/thermal/power/flexibility); ideal
on/off thrusters (no minimum impulse bit, plume, valve lag, slosh); constant Isp.
Known failure modes made explicit: torque saturation (demand > motor max),
**speed saturation with authority loss** (at the limit, further demand in the same
direction is refused — bus feels zero torque), propellant depletion (thrust scaled
to exactly reach the dry floor). Mass was unnecessary for ideal point-mass gravity
(`a = -mu r/r^3` is mass-independent) but is required here because `F = m a` and
the rocket equation couple thrust to a depleting mass state.

## 5. Verification

11 Catch2 test cases / 112 assertions: sign convention per axis, exact speed
integration with momentum-book closure `I d(omega) = tau dt`, saturation reporting,
lag closed form `1 - exp(-dt/tau)`, 2-second deterministic buildup to a computed
`I w_max / tau` limit with refused authority and recovery-direction authority,
thruster `r x F` hand case (`[0,.5,0]x[2,0,0] = [0,0,-1]`), 90-degree-yaw ECI
rotation, `m_dot` equation check, exact 10-second burn depletion, dry-floor clamp,
assembly aggregation, zero-command regression (identically zero input), invalid
rejection. Independent Python oracle `independent_actuator_reference.py`
cross-checks every scalar equation. Demo telemetry + two figures in
`artifacts/figures/m14_*`.

## What you should now understand

1. Why does a positive wheel motor torque produce a negative bus torque?
2. Why is total (bus + wheel) momentum conserved but bus-only momentum is not?
3. What distinguishes speed saturation from authority loss?
4. Why does `tau = r x F` vanish for a through-CM thruster?
5. Why did mass become a state variable only at M14?
6. Why must actuator limits raise flags instead of silently clipping?
