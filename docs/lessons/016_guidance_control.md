# Lesson 016 — Guidance & Control Foundations

## 1. Physical Problem & Motivation

M14 built the muscles (wheels, thrusters, limits); M15 built the map (LVLH,
relative state, CW prediction). M16 connects intent to force: given an estimated
attitude/rate or relative state, what torque/acceleration should the actuators be
asked for? Two decoupled problems, deliberately open-loop here (no
estimate→command→actuator→truth closure until M17):

- **Attitude**: detumble the bus and hold/slew an orientation with body torque.
  The error lives on SO(3): raw quaternion subtraction is geometrically
  meaningless (breaks the unit norm, mixes the double cover). The controller
  uses the error-quaternion vector part, the same object as the M13C MEKF
  residual.
- **Translation**: hold/track an LVLH offset with thrust. The CW linearization
  says the local plant is coupled oscillators plus drift, so the law adds the
  CW natural acceleration as feedforward and lets PD fight only residual error.

## 2. State Variables, Inputs, Outputs

| Symbol | Frame | Units | Meaning |
|---|---|---|---|
| `e_att = vec(q_des* ⊗ q_cur)`, w ≥ 0 | desired-frame vector | rad-ish (half-angle) | attitude error |
| `theta = 2 acos(|q_err.w|)` | scalar | rad in [0, π] | geodesic error angle |
| `e_rate = w − w_des` | BODY | rad/s | rate error |
| `tau_des = −Kp e_att − Kd e_rate` | BODY | N*m | desired torque |
| `e_rho, e_v` | LVLH | m, m/s | relative position/velocity errors |
| `a_cmd = a_des − ff − Kp e_rho − Kd e_v` | LVLH | m/s² | commanded acceleration |
| `K = [kp, kv]` | scalar axis | 1/s², 1/s | LQR optimal gains |

Commands enter as estimates by contract; controllers never touch truth (M17
enforces the wiring). Torque/force limits live downstream in M14; M16E only
pairs desired with achieved and raises the flag.

## 3. Governing Equations

Attitude PD with the half-angle note: `q_err ≈ [1, e/2]` for small errors, so
the effective angular stiffness is `Kp/2` and one axis reads
`I θ_dd + Kd θ_d + (Kp/2) θ = 0` (stable for Kp, Kd > 0, ζ = Kd/√(2KpI)).
Rate damping is the attitude-blind special case `tau = −Kw(w − w_cmd)`.
Gain suggestion inverts the 2% settle-time relation `wn = 4/(ζ ts)`:
`Kp = 2I wn²`, `Kd = 2ζwnI` — an engineering basis, not a magic number.
Translation feedforward `ff = [3n²x + 2nvy, −2nvx, −n²z]` cancels the CW plant
exactly when the linearization holds. LQR on the double integrator
`A = [[0,1],[0,0]], B = [0,1]'` closes in radicals: `p12 = √(q_pos·r)`,
`p22 = √(r(2p12 + q_vel))`, `K = [p12, p22]/r`; unit weights give the textbook
`K = [1, √3]`, ζ = √3/2. Saturation is a per-axis clamp with deficit reporting.

## 4. Assumptions and Failure Modes

Rigid diagonal-inertia bus; constant positive gains; no integral action (bias
torques like gravity gradient become M17 disturbances); circular reference for
the translation feedforward; allocation of LVLH acceleration onto body thrusters
deferred to M17. Explicit failure accounting: over-limit commands saturate with
a flag (M16E); authority loss and desaturation remain M14/M17 concerns. Gains
tuned for a (10, 20, 30) bus do not transfer to other inertias — retune via
`suggest_pd_gains`, never by copying numbers.

## 5. Verification

11 cases / 65 assertions: per-axis hand values, double-cover invariance,
linearity in both errors, restoring direction, PD≡rate-damping reduction,
feedforward closed form + exact setpoint cancellation, Riccati residuals,
Hurwitz eigenvalues, per-axis clamp symmetry, invalid rejection. Metrics:
detumble of the master-spec tumble settles in < 60 s (budget 100 s) with bounded
transient and sub-0.005° terminal error; 100 m radial offset closes to 4 mm in
one hour with peak command 1.4 mm/s². Independent Python oracle cross-checks all
five laws. Demo telemetry (hold, 90° slew, station keeping) + three figures.

## What you should now understand

1. Why is raw quaternion subtraction not an attitude error?
2. Where does the factor of 2 between Kp and angular stiffness come from?
3. Why does PD with zero attitude error reduce to rate damping?
4. What does the CW feedforward cancel, and when is the cancellation exact?
5. What tradeoff do LQR weights Q and R express?
6. Why do controllers never clip, and where does saturation get accounted?
