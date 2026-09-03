# Lesson 015 — Rotating-Frame Kinematics & Relative Orbital Dynamics

## 1. Physical Problem & Motivation

Rendezvous is flown relative to the target, not relative to Earth. ECI truth
states differ by kilometres with millimetre-level operational significance buried
in the difference, so the mission needs a target-centered LVLH frame carrying
relative position/velocity with direct meaning: x radial, y along-track,
z orbit-normal (AstraDock M06 convention). M06 built that frame geometrically but
froze it in time. The frame rotates at orbit rate, so forming relative velocity
demands the frame's angular velocity and the transport theorem — otherwise the
"relative velocity" is silently the projected inertial velocity, wrong by
`omega x rho` (≈ 2.2 m/s at 2 km in LEO, caught explicitly by test M15C).

For a circular reference and small separations, the linearized Clohessy-Wiltshire
(Hill) equations predict relative motion in closed form — the targeting engine of
far-field rendezvous. M15 implements the exact CW transition plus the nonlinear
ECI truth it must be judged against, and maps where the linearization fails.

## 2. State Variables, Inputs, Outputs

| Symbol | Frame | Units | Meaning |
|---|---|---|---|
| `omega_LVLH/ECI = h_vec / r^2` | ECI | rad/s | LVLH frame rate (exact, any eccentricity) |
| `alpha = dh/dt/r^2 - h(2r.v/r^4)` | ECI | rad/s^2 | LVLH angular acceleration |
| `rho_lvlh` | LVLH | m | chaser-minus-target relative position |
| `v_rel_lvlh` | LVLH | m/s | rotating-frame derivative (sensor-observed rate) |
| CW `X = [rho; v_rel]` | LVLH | m, m/s | 6-state linearized relative state |
| `n = sqrt(mu/R^3)` | scalar | rad/s | circular mean motion (CW parameter) |

Target/chaser roles are function arguments, never inferred.

## 3. Governing Equations

Frame rate: `omega = h/r^2` since `r x v = omega r^2` (radial velocity drops out
of the cross product; omega ⟂ r). Circular case: `|omega| = n`.
Transport: `(dv/dt)_I = (dv/dt)_R + omega x v`; acceleration composition with
Coriolis `2w x v`, centrifugal `w x (w x r)`, Euler `alpha x r`, plus origin term.
Relative state: `rho_lvlh = C rho_eci`, `v_rel = C dv_eci - omega_lvlh x rho_lvlh`.
CW: `x_dd - 2n y_d - 3n^2 x = ax`, `y_dd + 2n x_d = ay`, `z_dd + n^2 z = az`, with
the Schaub & Junkins exact `Phi(t)` closed form; forced arcs via RK4 on the CW
right-hand side.

## 4. Assumptions and Failure Modes

CW assumes circular reference, `separation << R`, two-body gravity both vehicles,
short-to-moderate horizons. Failure modes quantified (M15E): error scales
quadratically with separation (0.027 m at 100 m → 1713 m at 25 km after one
orbit) and linearly with horizon (3-orbit ≈ 3× 1-orbit). Eccentric reference,
differential J2/drag, and long horizons all break CW — none are modeled inside
it. Collinear/zero-momentum states reject (LVLH undefined, same policy as M06).

## 5. Verification

11 cases / 69 assertions: circular `n` along orbit normal, `h/r^2` identity with
radial-velocity independence, transport splits both directions, Coriolis/
centrifugal/Euler hand values, full composition identity, circular `alpha = 0`,
rate-consistency residual, ECI↔LVLH round trip, antisymmetry/zero cases,
projection-vs-rotating distinction, CW textbook cases (drift hold, radial growth,
cross-track oscillator, bounded ellipse closing over one period, `Phi(0) = I`),
RHS equation match, RK4-vs-closed-form to 1e-6, invalid rejection. Independent
Python oracle cross-checks rate, relative state, and CW prediction. Demo sweep +
5 km overlay telemetry with two figures (`artifacts/figures/m15_*`).

## What you should now understand

1. Why does `omega = h/r^2` hold for eccentric orbits, not just circular ones?
2. Why is the projected inertial velocity difference not the relative velocity?
3. What three terms appear when differentiating twice in a rotating frame?
4. Why does `vy0 = -2n x0` give a bounded relative ellipse?
5. Why does CW error scale as separation squared and grow with horizon?
6. When must you abandon CW for nonlinear ECI truth?
