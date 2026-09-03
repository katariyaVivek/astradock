# AstraDock — Long-Horizon Master Execution Prompt
## M14 → M25: Actuation, Relative Dynamics, GNC, Rendezvous, Docking, FDIR, Mission Framework, Computer Vision, ML, High-Fidelity Physics, and Final Capstone

You are taking over the AstraDock aerospace simulation repository after completion of M13.

This document is a **long-horizon execution specification**. It is not a request to write every future subsystem in one giant change. Execute the roadmap **milestone-by-milestone**, preserving a clean engineering boundary at every step.

The repository itself is authoritative for current code, naming, existing roadmap details, and current test counts. This file supplies the intended long-term direction and quality standard.

---

# 0. ABSOLUTE OPERATING RULES

Before changing anything:

1. Read `ASTRADOCK_CODEX_MASTER.md`.
2. Read `AGENTS.md`.
3. Read `ROADMAP.md`.
4. Read `README.md`.
5. Read `docs/architecture.md`.
6. Read `docs/curriculum.md`.
7. Inspect `cpp/`, `tests/`, `tools/`, `python/`, `docs/`, and relevant `artifacts/`.
8. Run the complete current test suite.
9. Record the real baseline before modifications.
10. Inspect the Git working tree and current branch.

Never assume the repository is exactly identical to this specification.

If the repository already contains additional completed work, preserve it.
If names or numbering differ from the proposed roadmap, reconcile with the repository instead of reverting valid progress.

Do not destroy historical evidence.

Do not begin two major milestones simultaneously.

Do not mark work complete without validation evidence.

---

# 1. VERIFIED STARTING POINT

The repository has completed:

```text
M01 — Vector Mathematics
M02 — Two-Body Gravity
M03 — Numerical Integration
M04 — Two-Body Orbital Propagation
M05 — Numerical Validation
M06 — Coordinate Frames / ECI / LVLH / DCM
M07 — Classical Orbital Elements
M08 — Quaternion Attitude Representation
M09 — Rigid-Body Attitude Dynamics & Quaternion Kinematics
M10 — Integrated 6-DOF Spacecraft State
M11 — Spacecraft Environment & Force/Torque Models
M12 — Sensor Simulation & Measurement Models
M13 — State Estimation & Navigation
```

M13 was intentionally consolidated from internal implementation stages:

```text
M13A — Translational GNSS EKF
M13B — IMU-Aided Translational Propagation
M13C — Attitude Error-State + Range Update
M13D — Integrated Multi-Rate Navigation
```

These are historical development stages, not separate permanent top-level milestones.

The current M13 baseline reported by the repository is:

```text
223 / 223 CTest tests passing
15-state integrated navigation estimator
multi-rate IMU / Star Tracker / Range / GNSS
independent Python integrated-navigation oracle
```

Use the actual repository test count if it has subsequently changed.

---

# 2. LONG-TERM SYSTEM OBJECTIVE

AstraDock should evolve into a transparent, modular, independently verified aerospace simulation and autonomy platform.

The intended final information flow is:

```text
                         MISSION / SCENARIO
                                │
                                ▼
                            GUIDANCE
                                │
                                ▼
                             CONTROL
                                │
                                ▼
                            ACTUATORS
                                │
                                ▼
     ENVIRONMENT ─────────► SPACECRAFT TRUTH
                                │
                                ▼
                             SENSORS
                                │
                                ▼
                           NAVIGATION
                                │
                                ▼
                       ESTIMATED STATE
                                │
                                └──────────► GUIDANCE / CONTROL
```

Parallel advanced capabilities eventually feed the same system:

```text
OPTICAL / COMPUTER VISION
          │
          ▼
   relative navigation

MACHINE LEARNING
          │
          ▼
 advisory diagnostics / perception / residual models
```

ML must remain subordinate to validated safety/physics logic where mission-critical outputs are involved.

---

# 3. ARCHITECTURAL INVARIANTS

These rules apply from M14 through project completion.

## 3.1 Truth isolation

Truth may be used by:

- validation harnesses,
- Monte Carlo scoring,
- independent oracles,
- plotting and reports.

Truth must never silently enter:

- navigation,
- guidance,
- control,
- actuator command generation,
- fault detection,
- ML inference during operational simulation.

Do static audits where practical.

## 3.2 Frame discipline

Every physical vector must have an explicit frame.

Use names such as:

```text
position_eci_m
velocity_eci_mps
force_eci_N
torque_body_Nm
angular_velocity_body_rad_s
relative_position_lvlh_m
```

Do not introduce ambiguous frame-sensitive APIs.

## 3.3 Units

Use SI internally:

```text
m, m/s, m/s², N, N·m, kg, kg·m², rad, rad/s, s
```

Make km, degrees, minutes, etc. explicit at boundaries.

## 3.4 Time

Every timestamp-sensitive component must preserve actual measurement/event time.

Never assume a nominal period equals actual `dt`.

## 3.5 Determinism

For a fixed configuration and seed, stochastic simulations should be reproducible within the supported execution environment.

Do not claim cross-platform bitwise equality without evidence.

## 3.6 No magical tuning

Never choose:

```text
Q, R, controller gains, thresholds, actuator limits
```

merely because that makes a plot look better.

Document the engineering basis for important tuning choices.

## 3.7 No feature-before-verification

A half-verified physics model must not be hidden behind a larger subsystem.

Implement one physical effect, verify it, then compose it.

---

# 4. UNIVERSAL MILESTONE WORKFLOW

For every milestone:

```text
inspect
  ↓
plan
  ↓
implement
  ↓
unit / analytical tests
  ↓
independent verification
  ↓
scenario validation
  ↓
Monte Carlo / statistical validation where appropriate
  ↓
telemetry / plots
  ↓
documentation
  ↓
full regression
  ↓
exit gate
```

A milestone is complete only after its exit gate passes.

If a real blocker is found:

```text
reproduce
→ classify
→ regression test
→ fix
→ rerun
```

Never hide a defect by loosening a tolerance without a mathematical reason.

---

# 5. DOCUMENTATION CONSOLIDATION RULE

Do not permanently create dozens of files named:

```text
M14A
M14B
M14C
...
```

Internal implementation stages are allowed.

When a top-level milestone is complete:

```text
one canonical lesson
one canonical validation report
one roadmap status
one architecture entry
```

Historical stage documents may be archived under:

```text
docs/archive/<milestone>/
```

when preserving them is useful.

The active documentation must remain compact enough for a future coding agent to understand the system.

---

# 6. PROPOSED REMAINING ROADMAP

Use the repository's existing M14+ roadmap as the primary authority. Where no later roadmap exists, use this sequence:

```text
M14 — Spacecraft Actuator Dynamics & Modeling
M15 — Rotating-Frame Kinematics & Relative Orbital Dynamics
M16 — Guidance & Control Foundations
M17 — Integrated Closed-Loop GNC
M18 — Rendezvous & Proximity Operations
M19 — Autonomous Docking
M20 — Fault Detection, Isolation & Recovery
M21 — Mission Simulation / Scenario / Monte Carlo Framework
M22 — Computer Vision & Optical Navigation
M23 — Machine Learning for Space Systems
M24 — High-Fidelity Environment / Time / Ephemeris Upgrade
M25 — Full-System Verification & Capstone Mission
```

Do not renumber earlier completed milestones.

---

# M14 — SPACECRAFT ACTUATOR DYNAMICS & MODELING

## Objective

Move from abstract force/torque commands to physically constrained actuator behavior.

The actuator layer should include:

```text
reaction wheels
thrusters
momentum storage
torque / force saturation
motor / response dynamics
propellant consumption
mass depletion when physically required
```

Do not begin GNC yet.

## M14A — Reaction Wheels

Implement:

```text
wheel inertia
wheel angular speed
commanded motor torque
achieved motor torque
maximum torque
maximum speed
stored angular momentum
saturation status
```

The spacecraft reaction torque must obey the chosen sign convention.

Verify angular-momentum exchange independently.

## M14B — Wheel Saturation

Create a deterministic momentum buildup case.

Show:

```text
command
→ wheel torque
→ wheel speed / momentum
→ saturation
→ loss of available authority
```

Never silently clip commands.

Return commanded vs achieved actuator output where useful.

## M14C — Thrusters

Represent a simple thruster with:

```text
mounting_position_body_m
thrust_direction_body
maximum_thrust_N
command / firing state
Isp
```

For an offset thruster:


action force produces:

```text
F_B = thrust_direction * thrust
τ_B = r_B × F_B
```

Transform force to ECI using the current attitude.

## M14D — Propellant

If thruster modeling requires mass:


dot m = -T/(I_sp g_0)

Use a physically meaningful spacecraft mass state.

Document why mass was not required for ideal point-mass gravity but is required for force-to-acceleration and propulsion modeling.

Do not implement tank slosh.

## M14 Verification

Required checks:

```text
reaction torque sign
wheel momentum conservation
wheel speed limits
torque saturation
thruster force direction
moment arm torque
fuel consumption
mass conservation / depletion
zero command
```

Independent analytical cases are mandatory.

## M14 Regression

With all actuator commands disabled:

```text
M14 truth = previous truth baseline
```

within established numerical tolerances.

## M14 Deliverables

Create a canonical:

```text
docs/lessons/014_actuators.md
docs/validation/m14_actuator_validation.md
```

plus deterministic telemetry and useful actuator plots.

---

# M15 — ROTATING-FRAME KINEMATICS & RELATIVE ORBITAL DYNAMICS

M06 introduced coordinate projection into LVLH. M15 must implement the actual kinematics of moving frames.

## M15A — LVLH Angular Velocity

Derive the LVLH frame angular velocity relative to inertial coordinates.

Handle non-circular states where the simple mean-motion relationship does not apply.

## M15B — Rotating-Frame Derivatives

Distinguish carefully between:

```text
inertial derivative
rotating-frame derivative
coordinate transformation
```

Implement the correct velocity/acceleration relations.

## M15C — Relative State

Create a clear relative state:

```text
relative_position_lvlh_m
relative_velocity_lvlh_mps
```

with explicit target/chaser semantics.

## M15D — Clohessy-Wiltshire / Hill

For a circular reference orbit implement:

\[
\ddot x - 2n\dot y - 3n^2x = a_x
\]

\[
\ddot y + 2n\dot x = a_y
\]

\[
\ddot z + n^2z = a_z
\]

Clearly state the assumptions:

- circular reference orbit,
- linearized dynamics,
- small separation,
- appropriate time interval.

Do not imply CW is universally valid.

## M15E — Validation

Compare:

```text
nonlinear ECI truth → LVLH relative state
```

against:

```text
CW/Hill prediction
```

over several separation scales and durations.

Quantify where the linear model breaks down.

---

# M16 — GUIDANCE & CONTROL FOUNDATIONS

Introduce controllers only after actuator models are physically constrained.

## M16A — Attitude Control

Implement a quaternion-error-based PD controller.

Inputs:

```text
estimated attitude
estimated angular velocity
desired attitude
desired angular velocity
```

Output:

```text
desired body torque
```

Never use raw quaternion subtraction.

## M16B — Rate Damping

Use:

\[
\tau=-K_\omega(\omega-\omega_{cmd})
\]

Validate principal-axis cases analytically.

## M16C — Translational / Relative Control

Start with simple linear relative-state control where appropriate.

Do not use arbitrary acceleration commands that violate actuator limits.

## M16D — LQR

For one clean linearized problem implement a small LQR example.

Teach:

```text
A, B, Q, R, K
```

Do not create a general-purpose optimal-control framework.

## M16E — Saturation

Controller output must flow through actuator limits.

Measure:

```text desired
achieved
saturated
```

## M16 Validation

Measure:

```text settling time
overshoot
steady-state error
control effort
saturation
stability
```

Use deterministic cases and independent linear-system checks.

---

# M17 — INTEGRATED CLOSED-LOOP GNC

Connect:

```text
truth
→ sensors
→ navigation
→ guidance
→ control
→ actuators
→ truth
```

## Critical Rule

Guidance and control see estimates, not truth.

## M17A — Attitude Hold

Closed loop:

```text star tracker / gyro
→ navigation
→ attitude controller
→ reaction wheel / thruster
→ truth
```

Demonstrate convergence from nonzero attitude error.

## M17B — Maneuver

Perform a deterministic attitude or orbital maneuver using realistic actuator constraints.

## M17C — Full Telemetry

Record:

```text truth
estimate
reference
control error
commanded torque/force
achieved torque/force
saturation
sensor validity
filter covariance
```

## M17D — Monte Carlo

At least 100 seeded runs varying:

```text initial state
sensor noise
sensor bias
actuator variation
timing
```

Measure:

```text mission success
settling time
maximum state error
fuel use
wheel momentum
control saturation
```

Do not hide failed runs.

---

# M18 — RENDEZVOUS & PROXIMITY OPERATIONS

Use M15 relative dynamics and M17 closed-loop control.

## M18A — Relative Navigation

Represent target/chaser states in LVLH with explicit frames.

## M18B — Rendezvous Guidance

Implement a simple sequence:

```text
transfer
approach
station keeping
final approach
```

Do not build an industrial mission planner.

## M18C — Keep-Out Zone

Define and enforce:

```text keep-out radius
approach corridor
relative velocity limits
```

## M18D — Validation

Test:

```text nominal rendezvous
burn error
state uncertainty
actuator saturation
sensor dropout
missed maneuver
```

Quantify robustness rather than only displaying a successful trajectory.

---

# M19 — AUTONOMOUS DOCKING

## M19A — Docking Frames

Define:

```text target body frame
docking-port frame
chaser body frame
approach axis
capture envelope
```

All transforms must be explicit.

## M19B — Relative Pose

Use:

```text relative position
relative velocity
relative attitude
relative angular velocity
```

## M19C — Final Approach

Enforce:

```text closing-speed limit
lateral error limit
attitude alignment
angular-rate limit
keep-out conditions
```

## M19D — Contact Model

Use a simple compliant model:

```text spring-damper
contact impulse where justified
capture latch
```

Do not model detailed structural deformation.

## M19E — Docking Acceptance

Docking succeeds only if:

```text position inside capture envelope
attitude aligned
relative velocity below limit
relative angular velocity below limit
contact state stable
```

## M19F — Abort

Implement deterministic abort conditions.

---

# M20 — FAULT DETECTION, ISOLATION & RECOVERY

Build FDIR using existing statistical and physical diagnostics.

## Sensor faults

Simulate:

```text dropout
bias jump
stuck output
noise increase
outlier
time offset
```

## Actuator faults

Simulate:

```text reduced authority
stuck actuator
thruster failure
wheel saturation
wrong thrust magnitude
```

## Detection

Use:

```text innovation
NIS
state residuals
actuator tracking residuals
```

## Isolation

Return:

```text suspected fault
confidence
ambiguous/unknown
```

Never force a false certainty.

## Recovery

Support:

```text sensor exclusion
degraded sensor set
actuator redistribution
safe mode
```

## Metrics

Measure:

```text detection latency
false alarm rate
miss rate
recovery success
mission impact
```

---

# M21 — MISSION SIMULATION / SCENARIO / MONTE CARLO FRAMEWORK

The codebase is now large enough to need disciplined scenario execution.

## M21A — Scenario Definition

Support configuration for:

```text initial state
spacecraft parameters
environment
sensors
navigation
guidance
control
actuators
faults
mission timeline
```

Use a standard configuration format if appropriate.

Do not invent a custom scripting language.

## M21B — Timeline

Represent events:

```text mode changes
sensor changes
actuator commands
fault injection
guidance transitions
```

## M21C — Run Identity

Every run must expose:

```text scenario ID
seed
configuration hash
software version / commit where practical
```

## M21D — Monte Carlo

Support:

```text N runs
repeatable seeds
summary statistics
percentiles
failure classification
```

Parallel execution is optional; determinism is mandatory.

## M21E — Regression Missions

Define canonical missions for:

```text orbit
attitude control
station keeping
rendezvous
docking
sensor failure
actuator failure
navigation failure
```

---

# M22 — COMPUTER VISION & OPTICAL NAVIGATION

Only begin after physics, navigation, and GNC are mature.

## M22A — Synthetic Camera

Implement a minimal camera model:

```text focal length
principal point
resolution
optional distortion
pixel noise
```

## M22B — Geometric Targets

Start with:

```text docking fiducials
known landmarks
planet limb
```

Do classical CV first.

## M22C — Projection

Verify:

```text 3D point
→ camera frame
→ image plane
```

with independent analytical tests.

## M22D — PnP / Relative Pose

Estimate relative pose from known target geometry.

Frames must be explicit:

```text camera
body
target
inertial
```

## M22E — Optical Navigation

Produce relative pose/bearing estimates.

Compare to truth only in the validation harness.

## M22F — Failures

Simulate:

```text missed feature
false feature
occlusion
blur
noise
illumination variation
```

## M22G — Validation

Measure:

```text pose error
pixel residual
failure rate
sensitivity to image noise
```

---

# M23 — MACHINE LEARNING FOR SPACE SYSTEMS

ML enters only after deterministic classical baselines exist.

The purpose is to make an actual aerospace ML contribution, not a generic classifier.

Preferred application families:

```text sensor anomaly classification
actuator fault classification
telemetry anomaly detection
optical feature classification
learned residual correction
```

## M23A — Dataset Generation

Use AstraDock scenarios to create datasets containing:

```text scenario ID
seed
truth metadata for labeling only
measurements
fault label
configuration
```

Keep training/validation/test scenarios separate.

## M23B — Leakage Prevention

Never split time-series rows randomly if that leaks trajectory-specific information.

Prefer split by:

```text scenario
mission
seed
```

## M23C — Classical Baseline

Before ML implement:

```text rule-based baseline
statistical baseline
```

## M23D — ML Model

Choose a model appropriate to the data.

Do not use deep learning when a simpler model is enough.

## M23E — Metrics

Use:

```text precision
recall
F1
false alarm rate
miss rate
latency
calibration
```

Do not report only accuracy.

## M23F — Safe Integration

ML outputs must be advisory or pass through validated safety logic.

Do not let an unconstrained neural model directly command actuators.

## M23G — Physics + ML Residual

For a high-value experiment:

```text physics model
→ residual
→ ML residual model
→ corrected prediction
```

Evaluate on unseen mission scenarios.

Never use ML to conceal a known physics bug.

---

# M24 — HIGH-FIDELITY ENVIRONMENT / TIME / EPHEMERIS UPGRADE

Upgrade the simplified M11 environment only after the autonomy stack exists.

## M24A — ECI / ECEF

Implement explicit Earth rotation and a defensible ECI/ECEF transformation.

Document time assumptions.

## M24B — Time

Distinguish, as appropriate:

```text simulation elapsed time
UTC-like civil time
TT-like dynamical/reference time
```

Do not claim standards compliance without the required data/standards.

## M24C — WGS-84 Geodetic

Implement:

```text ECEF → latitude/longitude/altitude
geodetic → ECEF
```

Test round trips and edge cases.

## M24D — Higher-Order Gravity

Optionally extend beyond J2 with J3/J4 or another justified model.

Verify independently before use.

## M24E — Ephemerides

Upgrade third-body prescribed positions to a documented analytical/reference ephemeris where practical.

## M24F — Atmosphere

Upgrade the educational density model only if an appropriate source/model is available and its validity is documented.

Never imply operational atmospheric prediction from a simplified model.

---

# M25 — FULL-SYSTEM VERIFICATION & CAPSTONE MISSION

M25 is a consolidation/verification milestone, not a feature dump.

## M25A — Capstone Mission

Construct one deterministic mission exercising:

```text environment
6-DOF truth
sensors
15-state navigation
actuators
guidance
control
relative navigation
rendezvous
proximity operation / docking
fault injection
recovery
```

Use manageable mission phases.

Suggested structure:

```text
1. orbit initialization
2. attitude acquisition
3. navigation convergence
4. maneuver
5. relative transfer
6. rendezvous
7. final approach
8. docking
9. controlled fault
10. recovery
11. mission completion
```

## M25B — Truth Isolation Audit

Perform static and dynamic checks that truth is not used by:

```text navigation
guidance
control
FDIR
ML operational path
```

## M25C — Interface Audit

Check every module boundary for:

```text frame mismatch
unit mismatch
time mismatch
sign mismatch
state-order mismatch
```

## M25D — Monte Carlo

At least 100 runs, preferably more if runtime allows.

Randomize:

```text initial conditions
sensor noise/bias
actuator parameters
environment parameters
timing
fault timing
```

Report:

```text success rate
median
5th percentile
95th percentile
failure cases
```

Never hide failed runs.

## M25E — Benchmark Suite

Freeze benchmark scenarios:

```text circular LEO
eccentric orbit
J2 propagation
drag decay
third-body perturbation
rigid-body tumble
attitude control
navigation dropout
actuator saturation
rendezvous
docking
fault recovery
```

For each record:

```text initial conditions
expected behavior
quantitative acceptance criteria
```

## M25F — Reproducibility

Verify a clean clone can:

```text configure
build
run tests
run canonical missions
regenerate required plots
```

without hidden developer-specific files.

## M25G — Performance

Measure:

```text truth-only runtime
truth+sensors runtime
truth+navigation runtime
closed-loop runtime
Monte Carlo throughput
```

Profile before optimizing.

---

# 9. UNIVERSAL ANALYTICAL VERIFICATION REQUIREMENTS

Whenever a closed-form case exists, use it.

Examples:

```text actuator force / torque geometry
principal-axis wheel behavior
principal-axis attitude control
CW equations
simple Kalman update
sensor projection
orbital elements
rotation transformations
```

Do not validate an analytical case by calling the same production helper twice.

---

# 10. UNIVERSAL INDEPENDENT ORACLE REQUIREMENT

Maintain independent Python or other-reference implementations for critical systems.

The oracle must:

- not import C++ production code,
- not import production results as expected values,
- independently implement enough mathematics to cross-check the target,
- be clearly labeled as validation-only.

For ML, the independent oracle can instead be a frozen classical baseline or independently generated evaluation set when an equation-level oracle is inappropriate.

---

# 11. UNIVERSAL JACOBIAN REQUIREMENT

For every future EKF/control/optimization Jacobian:

```text
analytical derivation
        ↕
independent finite-difference audit
```

Use multiple perturbation scales.

Look for the roundoff/truncation tradeoff.

Never choose a finite-difference epsilon arbitrarily and test only one point.

---

# 12. UNIVERSAL COVARIANCE REQUIREMENT

Any future estimator must preserve:

```text symmetry
positive-semidefinite meaning
finite values
reasonable conditioning
```

Continue avoiding explicit general inverses where stable solves are possible.

Prefer Joseph-form covariance updates where applicable.

---

# 13. UNIVERSAL FRAME-CONVENTION REQUIREMENT

Before any module consumes another module's vector, document:

```text source frame
destination frame
rotation direction
active/passive convention
units
```

A mathematically perfect rotation can still be physically wrong if the convention is reversed.

Always include at least one hand-derived sign/orientation regression.

---

# 14. UNIVERSAL SENSOR / ACTUATOR SEPARATION

Sensors observe truth.

Actuators influence truth.

Navigation estimates state.

Guidance chooses desired state/trajectory.

Control turns errors into commands.

Never collapse these roles.

---

# 15. UNIVERSAL FAILURE-PATH REQUIREMENT

For every major subsystem test:

```text nominal case
boundary case
invalid input
fault/dropout case
recovery case where appropriate
```

The failure behavior must be explicit.

Do not silently manufacture valid-looking values after a fault.

---

# 16. UNIVERSAL MONTE CARLO REQUIREMENT

For stochastic or robustness-sensitive systems:

- deterministic seeds,
- enough samples for the claimed conclusion,
- clear scenario splitting,
- confidence/percentiles,
- raw failure counts.

Do not describe 100 runs as statistically proving mission safety.

Use language proportional to the evidence.

---

# 17. UNIVERSAL PERFORMANCE REQUIREMENT

Do not optimize before profiling.

AstraDock favors:

```text correctness
clarity
maintainability
verifiability
```

over premature micro-optimization.

---

# 18. UNIVERSAL CODE QUALITY REQUIREMENT

Avoid:

```text global mutable truth
magic numbers
duplicated equations
duplicated constants
implicit frames
implicit units
silent clamping
exception swallowing
large inheritance hierarchies
```

Prefer:

```text small explicit structs
pure functions
clear namespaces
validated primitives
composition
```

---

# 19. FINAL DOCUMENTATION TARGET

At the end of the full roadmap, active documentation should remain approximately:

```text
README.md
ROADMAP.md
AGENTS.md
ASTRADOCK_CODEX_MASTER.md
docs/architecture.md
docs/curriculum.md

docs/lessons/<one canonical lesson per top-level milestone>

docs/validation/<one canonical validation report per top-level milestone>

docs/archive/                  # historical stage reports when useful
```

Do not allow completed milestones to become cluttered with competing "final" documents.

---

# 20. FINAL ROADMAP CONSOLIDATION

When each top-level milestone completes, `ROADMAP.md` must contain exactly one authoritative status for it.

For example:

```text
M14 — Spacecraft Actuator Dynamics & Modeling       ✅ COMPLETE
M15 — Rotating-Frame Kinematics & Relative Dynamics ✅ COMPLETE
M16 — Guidance & Control Foundations                ✅ COMPLETE
...
```

Internal implementation stages can be mentioned beneath each milestone but must not become separate roadmap milestones unless intentionally promoted by the project owner.

---

# 21. MASTER-SPEC UPDATE POLICY

At every completed major milestone update `ASTRADOCK_CODEX_MASTER.md` with:

```text
current completed milestones
current architecture
current test baseline
current known limitations
next milestone
```

Keep it concise enough to be read at the start of every future coding session.

Do not duplicate every validation table into the master specification.

---

# 22. MILESTONE COMPLETION REPORT FORMAT

At the end of every milestone report:

```text
Milestone:
Status:

Files created:
...

Files modified:
...

Tests:
previous = X/X
new = Y/Y
total = Z/Z

Independent verification:
...

Analytical checks:
...

Scenario validation:
...

Monte Carlo / statistics:
...

Bugs found:
...

Corrections:
...

Known limitations:
...

Documentation consolidated:
...

Build:
...

Warnings:
...

Ruff:
...

Determinism:
...

Truth isolation:
...

Scope guard:
PASS / FAIL

Next milestone:
...
```

Do not fabricate any metric.

Do not call nonzero floating-point residuals exact.

---

# 23. GLOBAL SCOPE GUARD

Until the milestone that explicitly introduces them, do not implement prematurely:

```text
full GNC
rendezvous
docking
FDIR
computer vision
ML
high-fidelity ephemerides
advanced atmospheric models
```

Likewise, when a milestone is focused on one component, do not drag future components into it merely because they are available.

---

# 24. M14 ENTRY INSTRUCTION

The immediate next milestone is M14 unless the repository already shows it complete.

Start by inspecting the actual roadmap.

Then implement:

```text
reaction wheels
thrusters
momentum storage
saturation
actuator dynamics
propellant / mass effects where needed
```

Do not implement guidance or control as part of M14.

---

# 25. PROJECT-END ACCEPTANCE STANDARD

AstraDock should only be considered complete when it can demonstrate, with documented evidence:

```text
[ ] validated orbital mechanics
[ ] validated coordinate frames
[ ] validated orbital elements
[ ] validated quaternion mathematics
[ ] validated rigid-body dynamics
[ ] validated environmental models
[ ] validated sensor models
[ ] validated integrated navigation
[ ] validated actuator models
[ ] validated rotating-frame dynamics
[ ] validated guidance
[ ] validated control
[ ] validated closed-loop GNC
[ ] validated rendezvous
[ ] validated docking/proximity operations
[ ] validated FDIR
[ ] reproducible scenario framework
[ ] validated optical navigation
[ ] responsibly evaluated ML capability
[ ] high-fidelity environment upgrades documented
[ ] end-to-end capstone mission
[ ] clean clone/build
[ ] comprehensive regression suite
[ ] consolidated documentation
```

The word **validated** means there is actual evidence appropriate to the model; it does not mean flight qualification.

---

# 26. FINAL ENGINEERING STANDARD

Throughout all future work, maintain this chain:

```text
PHYSICS
  ↓
MATHEMATICAL MODEL
  ↓
IMPLEMENTATION
  ↓
UNIT / ANALYTICAL VERIFICATION
  ↓
INDEPENDENT REFERENCE
  ↓
SCENARIO VALIDATION
  ↓
STATISTICAL / MONTE CARLO VALIDATION
  ↓
SYSTEM INTEGRATION
  ↓
REGRESSION
```

Never reverse the order by building a complicated system first and trying to explain it afterward.

---

# 27. FINAL PHILOSOPHY

AstraDock is not being built to produce the most impressive screenshot.

It is being built so that a technically trained person can inspect the system and answer:

```text
What equation is being used?
What assumptions does it make?
What frame is each quantity expressed in?
What units are used?
How was the implementation verified?
What independent evidence supports it?
What happens when the model is wrong?
What happens when a sensor fails?
What happens when an actuator saturates?
How does uncertainty propagate?
Can the result be reproduced?
```

The final system should therefore remain:

```text
modular
transparent
deterministic
independently verifiable
physically interpretable
numerically defensible
```

Proceed from the current repository state, one milestone at a time, beginning with M14.
