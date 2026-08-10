# AstraDock — Codex Master Project Specification

> **Purpose:** Build a serious, long-term aerospace engineering portfolio project that teaches orbital mechanics, spacecraft dynamics, GNC, estimation, control, computer vision, and machine learning while Codex acts as a senior engineering pair-programmer.
>
> **Primary constraint:** This repository is both an engineering project and a learning environment. Do not optimize for maximum code output. Optimize for correctness, physical clarity, testability, verification, and teaching value.

---

# 1. Project Overview

## Project Name

**AstraDock**

## Mission

Build an **Autonomous Orbital Rendezvous & Docking Simulator**.

A chaser spacecraft starts several kilometres away from a target spacecraft in Low Earth Orbit. The target may be tumbling.

The chaser must eventually be capable of:

1. propagating orbital motion,
2. simulating 6-DOF rigid-body spacecraft dynamics,
3. simulating realistic onboard sensors,
4. estimating position, velocity, and attitude from imperfect measurements,
5. navigating relative to the target spacecraft,
6. planning a safe rendezvous trajectory,
7. controlling translation and attitude,
8. entering a docking corridor,
9. estimating target pose using computer vision,
10. detecting sensor and actuator faults,
11. autonomously aborting unsafe approaches,
12. docking successfully,
13. validating performance through Monte Carlo testing,
14. comparing important dynamics results against an established astrodynamics framework such as Basilisk.

The final conceptual architecture is:

```text
                 ┌─────────────────────────┐
                 │   Orbital Dynamics      │
                 │ Target + Chaser 6-DOF   │
                 └────────────┬────────────┘
                              │ truth state
                 ┌────────────┴────────────┐
                 │                         │
             IMU / GNSS              Docking Camera
                 │                         │
                 ▼                         ▼
          Sensor Simulation        Vision Pose Estimator
                 │                         │
                 └────────────┬────────────┘
                              ▼
                        State Estimator
                         EKF / UKF later
                              │
                        estimated state
                              ▼
                           Guidance
                              │
                        desired state
                              ▼
                           Control
                              │
                 Thrusters / Reaction Wheels
                              │
                              ▼
                       Spacecraft Dynamics
                              │
                              └────────── LOOP
```

---

# 2. Why This Project Exists

The learner has a Machine Learning / Computer Science background and wants to become capable of working in the space sector.

This project must therefore bridge the gap between:

```text
Machine Learning / Software
            ↓
Numerical Engineering
            ↓
Orbital Mechanics
            ↓
Spacecraft Dynamics
            ↓
State Estimation
            ↓
Guidance, Navigation & Control
            ↓
Autonomous Space Systems
            ↓
Computer Vision + ML
```

Do **not** treat machine learning as the starting point.

Classical aerospace methods must be learned first.

Examples:

- analytical orbital mechanics before learned trajectory models,
- EKF before neural state estimation,
- PID/LQR before reinforcement-learning control,
- geometric pose estimation before learned end-to-end docking,
- deterministic anomaly checks before ML anomaly detection.

---

# 3. Codex Role

Codex must behave like:

- a senior aerospace software engineer,
- a numerical methods engineer,
- a GNC engineer,
- a C++ reviewer,
- a Python scientific-computing reviewer,
- a patient technical instructor,
- a test engineer.

Codex must **not** behave like a code generator whose goal is to fill the repository as quickly as possible.

---

# 4. Mandatory Learning Rules

Before implementing any unfamiliar concept involving:

- orbital mechanics,
- dynamics,
- coordinate frames,
- numerical integration,
- quaternions,
- estimation,
- filtering,
- guidance,
- control,
- spacecraft sensors,
- relative motion,
- rendezvous,
- docking,
- computer vision,
- machine learning,

Codex must first document the following:

1. **Physical problem**
2. **Why the concept is needed**
3. **Coordinate frames involved**
4. **Units**
5. **State variables**
6. **Inputs**
7. **Outputs**
8. **Governing equations**
9. **Assumptions**
10. **Numerical method**
11. **Known failure modes**
12. **Verification strategy**

Every major new concept must produce a lesson under:

```text
docs/lessons/
```

Example:

```text
docs/lessons/001_vectors.md
docs/lessons/002_two_body_dynamics.md
docs/lessons/003_rk4.md
docs/lessons/004_coordinate_frames.md
docs/lessons/005_quaternions.md
```

The lesson must explain the mathematics in understandable engineering language.

---

# 5. Fundamental Engineering Rules

## Units

Use **SI units internally** unless explicitly documented otherwise.

Examples:

```text
distance        metre
time            second
mass            kilogram
velocity        metre/second
acceleration    metre/second²
force           newton
torque          newton-metre
angle           radian
angular rate    radian/second
```

User-facing visualization may display kilometres or degrees.

Internal simulation code must remain unambiguous.

---

## Coordinate Frames

Never use position, velocity, acceleration, orientation, or force vectors without clearly specifying the frame.

Relevant frames will eventually include:

```text
ECI
ECEF
LVLH / Hill
spacecraft body frame
sensor frame
camera frame
target body frame
docking frame
```

Frame conversions must be explicit.

No hidden transformations.

---

## State Separation

Every simulation must distinguish:

```text
truth state
sensor measurement
estimated state
commanded state
```

Never overwrite one with another.

Example:

```text
truth.position_eci
measurement.gnss_position
estimate.position_eci
guidance.commanded_position
```

---

## Determinism

Every stochastic simulation must accept a random seed.

Given the same seed, results must be reproducible.

---

## Numerical Stability

Do not implement numerical algorithms without discussing:

- integration step,
- numerical error,
- convergence,
- conditioning,
- stability,
- tolerances.

---

## Assertions

Invalid physical states should fail clearly.

Examples:

- zero quaternion norm,
- negative mass,
- invalid time step,
- NaN state,
- impossible covariance matrix,
- invalid coordinate transformation.

---

# 6. Technology Stack

## Core Engineering Language

**C++20**

Use C++ for:

- mathematical primitives,
- orbital mechanics,
- dynamics,
- integrators,
- sensors,
- estimation,
- guidance,
- control,
- simulation core.

Avoid overengineering.

Prefer readable modern C++.

---

## Scientific / ML Language

**Python 3.12+**

Use Python for:

- plotting,
- analysis,
- scenario tooling,
- notebooks,
- ML,
- Monte Carlo post-processing,
- visualization,
- validation.

Useful libraries may later include:

```text
numpy
scipy
matplotlib
pandas
pyyaml
pytest
torch
opencv-python
```

Do not add dependencies until they are needed.

---

## Build

Use:

```text
CMake
CTest
```

---

## Testing

C++:

```text
Catch2 or GoogleTest
```

Choose one and document the decision.

Python:

```text
pytest
```

---

## Formatting / Static Analysis

Where practical:

```text
clang-format
clang-tidy
ruff
```

Keep configuration minimal and maintainable.

---

# 7. Suggested Repository Architecture

```text
astradock/
│
├── AGENTS.md
├── README.md
├── ROADMAP.md
├── CMakeLists.txt
├── pyproject.toml
│
├── cpp/
│   ├── math/
│   │   ├── vector3.hpp
│   │   ├── matrix.hpp
│   │   ├── quaternion.hpp
│   │   ├── rotations.hpp
│   │   └── frames.hpp
│   │
│   ├── dynamics/
│   │   ├── gravity.hpp
│   │   ├── orbit.hpp
│   │   ├── rigid_body.hpp
│   │   └── integrator.hpp
│   │
│   ├── sensors/
│   │   ├── imu.hpp
│   │   ├── gnss.hpp
│   │   ├── star_tracker.hpp
│   │   ├── range_sensor.hpp
│   │   └── camera.hpp
│   │
│   ├── navigation/
│   │   ├── kalman_filter.hpp
│   │   ├── ekf.hpp
│   │   └── relative_navigation.hpp
│   │
│   ├── guidance/
│   │   ├── rendezvous.hpp
│   │   ├── docking_corridor.hpp
│   │   └── trajectory.hpp
│   │
│   ├── control/
│   │   ├── pid.hpp
│   │   ├── lqr.hpp
│   │   └── attitude_controller.hpp
│   │
│   ├── faults/
│   │   └── fault_injector.hpp
│   │
│   └── sim/
│       ├── spacecraft.hpp
│       ├── target.hpp
│       ├── mission.hpp
│       └── simulation.hpp
│
├── python/
│   ├── analysis/
│   ├── visualization/
│   ├── scenarios/
│   └── bindings/
│
├── ml/
│   ├── pose_estimation/
│   └── anomaly_detection/
│
├── scenarios/
│   ├── basic_orbit.yaml
│   ├── attitude_control.yaml
│   ├── rendezvous.yaml
│   ├── docking.yaml
│   └── sensor_failure.yaml
│
├── tests/
│   ├── cpp/
│   └── python/
│
├── docs/
│   ├── architecture.md
│   ├── curriculum.md
│   ├── mathematics/
│   ├── gnc/
│   ├── orbital_mechanics/
│   └── lessons/
│
├── notebooks/
│
└── tools/
```

This architecture is a target, not a requirement to create every directory immediately.

Create directories only when needed.

---

# 8. Project Curriculum

The repository itself acts as the curriculum.

---

## PHASE 1 — Mathematical Foundation

Learn:

- scalars,
- vectors,
- dot product,
- cross product,
- norms,
- normalization,
- matrices,
- coordinate representation.

Implement:

```text
Vector3
```

Required tests:

- vector addition,
- subtraction,
- scalar multiplication,
- magnitude,
- normalization,
- dot product,
- cross product,
- orthogonality,
- known analytical cases.

Lesson:

```text
docs/lessons/001_vectors.md
```

---

# PHASE 2 — Two-Body Orbital Dynamics

Learn:

- Newtonian gravity,
- gravitational parameter,
- orbital state vectors,
- position,
- velocity,
- acceleration,
- two-body problem.

Implement gravitational acceleration:

```text
r_ddot = -mu * r / |r|^3
```

Do not hide this behind an astrodynamics library.

Learn:

```text
Earth gravitational parameter
Earth radius
orbital altitude
orbital velocity
specific orbital energy
specific angular momentum
```

Lesson:

```text
docs/lessons/002_two_body_dynamics.md
```

Verification:

For circular orbit:

```text
v = sqrt(mu / r)
```

Check gravitational acceleration analytically.

---

# PHASE 3 — Numerical Integration

Implement:

1. Euler integration
2. RK4

Compare them.

Demonstrate why Euler produces larger orbital energy drift.

Record:

```text
specific orbital energy vs time
angular momentum vs time
position error
```

Lesson:

```text
docs/lessons/003_numerical_integration.md
```

---

# PHASE 4 — Orbital Simulation

Create a complete propagator.

Initial scenario example:

```text
orbit altitude = 500 km
inclination    = 45 degrees
```

Expose:

```text
x
y
z
vx
vy
vz
```

Visualization should eventually show:

- Earth,
- orbit,
- spacecraft position,
- altitude,
- velocity,
- orbital period,
- orbital energy,
- angular momentum.

---

# PHASE 5 — Coordinate Frames

Implement and understand:

```text
ECI
ECEF
LVLH / Hill
body frame
```

Later:

```text
camera frame
sensor frame
target frame
docking frame
```

Every transform must have:

- documented convention,
- tests,
- inverse check,
- round-trip test.

Lesson:

```text
docs/lessons/004_coordinate_frames.md
```

---

# PHASE 6 — Quaternions and Attitude

Learn:

- Euler angles,
- rotation matrices,
- quaternions,
- quaternion normalization,
- quaternion multiplication,
- vector rotation,
- attitude propagation,
- gimbal lock.

Implement quaternion representation.

State:

```text
q = [qw, qx, qy, qz]
omega = [wx, wy, wz]
```

Tests:

- identity rotation,
- 90-degree known rotations,
- inverse rotation,
- norm preservation,
- quaternion ↔ rotation matrix consistency.

Lesson:

```text
docs/lessons/005_quaternions.md
```

---

# PHASE 7 — Rigid Body Dynamics

Simulate spacecraft rotation.

Learn:

- torque,
- angular momentum,
- inertia tensor,
- Euler rigid-body equations.

State:

```text
attitude
angular velocity
moment of inertia
external torque
```

Allow deliberate tumbling.

---

# PHASE 8 — Spacecraft Sensors

The navigation system must stop receiving truth directly.

Simulate:

## IMU

Accelerometer:

- white noise,
- bias,
- drift.

Gyroscope:

- white noise,
- bias,
- drift.

## GNSS

Simulate:

- position noise,
- velocity noise,
- update rate,
- dropout.

## Star Tracker

Simulate:

- attitude measurement noise,
- update rate,
- dropout.

## Range Sensor

Simulate relative range.

Architecture:

```text
true state
   ↓
sensor model
   ↓
bias + noise + sampling
   ↓
measurement
```

Lesson:

```text
docs/lessons/006_sensor_models.md
```

---

# PHASE 9 — Kalman Filtering

Begin with a basic linear Kalman filter.

Understand:

```text
state
measurement
prediction
innovation
covariance
process noise
measurement noise
Kalman gain
```

Then implement Extended Kalman Filter.

Dashboard/analysis must show:

```text
truth
raw measurements
estimate
estimation error
covariance
innovation
```

Never only display the final estimate.

Lesson:

```text
docs/lessons/007_kalman_filter.md
docs/lessons/008_extended_kalman_filter.md
```

---

# PHASE 10 — Attitude Estimation / Sensor Fusion

Fuse sensors such as:

```text
gyro
star tracker
other attitude measurements
```

Estimate:

```text
quaternion
angular velocity
gyro bias
```

Evaluate attitude error.

---

# PHASE 11 — Spacecraft #2

Introduce:

```text
Target spacecraft
Chaser spacecraft
```

Target may tumble.

Typical initial separation:

```text
5000 m
```

Both spacecraft propagate in orbit.

---

# PHASE 12 — Relative Orbital Mechanics

Learn:

- relative position,
- relative velocity,
- Hill frame,
- Clohessy-Wiltshire / Hill equations.

Mission distances:

```text
5000 m
1000 m
250 m
50 m
10 m
2 m
contact
```

Lesson:

```text
docs/lessons/009_relative_motion.md
```

---

# PHASE 13 — Guidance

Build an autonomous mission state machine.

Example:

```text
ACQUIRE_TARGET
      ↓
FAR_APPROACH
      ↓
RENDEZVOUS
      ↓
PROXIMITY_OPERATIONS
      ↓
FINAL_APPROACH
      ↓
ATTITUDE_MATCH
      ↓
DOCK
      ↓
CAPTURE_CONFIRMED
```

Abort conditions:

```text
SENSOR_FAILURE
EXCESSIVE_CLOSING_RATE
APPROACH_CORRIDOR_VIOLATION
NAVIGATION_DIVERGENCE
ATTITUDE_ERROR
THRUSTER_FAILURE
```

Unsafe states must trigger ABORT.

---

# PHASE 14 — Translational Control

Implement simplified thruster models.

Learn:

- force,
- acceleration,
- impulse,
- thrust limits,
- actuator saturation.

Start simple.

Do not immediately model realistic propulsion systems.

---

# PHASE 15 — Attitude Control

Start with:

```text
PID
```

Then:

```text
LQR
```

Simulate initial tumble:

```text
omega = [0.2, -0.1, 0.15] rad/s
```

Measure stabilization performance.

Plot:

```text
attitude error
angular velocity
control torque
settling time
```

Lesson:

```text
docs/lessons/010_pid_attitude_control.md
docs/lessons/011_lqr_control.md
```

---

# PHASE 16 — Docking Corridor

Implement safety geometry.

Define:

- approach axis,
- corridor radius,
- allowable closing velocity,
- hold points,
- keep-out sphere,
- abort maneuver.

Example hold points:

```text
250 m
50 m
10 m
2 m
```

The guidance system must not simply fly directly at the target.

---

# PHASE 17 — Docking Camera

Add a virtual docking camera.

Camera model must eventually include:

- intrinsic matrix,
- focal length,
- field of view,
- resolution,
- coordinate transform,
- projection,
- noise.

Start with deterministic synthetic geometry.

Do not start with ML.

---

# PHASE 18 — Classical Computer Vision

Before ML, implement or study:

- feature detection,
- feature matching,
- keypoints,
- PnP,
- relative pose,
- reprojection error.

Pipeline:

```text
image
 ↓
target detection
 ↓
keypoints
 ↓
geometric pose estimation
 ↓
relative position + attitude
```

Lesson:

```text
docs/lessons/012_visual_pose_estimation.md
```

---

# PHASE 19 — Machine Learning Pose Estimation

Only now introduce ML.

Potential approaches:

```text
CNN
Vision Transformer
keypoint regression
pose regression
segmentation + geometry
```

Synthetic training data may be generated using:

```text
Blender
procedural renderer
simulation camera
```

Model outputs may include:

```text
relative translation
relative quaternion
keypoints
confidence
```

The learned output should feed into the estimator.

---

# PHASE 20 — Fault Injection

Create a fault framework.

Faults:

```text
gyro_bias
gyro_failure
accelerometer_bias
GNSS_dropout
star_tracker_dropout
camera_occlusion
range_sensor_failure
thruster_failure
sensor_spike
sensor_drift
measurement_delay
stuck_measurement
```

Scenario example:

```yaml
fault:
  type: gyro_bias
  start_time: 120.0
  axis: z
  magnitude: 0.2
```

Navigation must react appropriately.

---

# PHASE 21 — Classical Fault Detection

Before ML anomaly detection, implement:

- range checks,
- innovation monitoring,
- residual thresholds,
- covariance consistency checks,
- rate-of-change checks.

Example:

```text
GYRO-Z ANOMALY

Expected:    0.013 rad/s
Measured:    0.216 rad/s
Residual:    7.9 sigma

NAVIGATION MODE:
DEGRADED
```

---

# PHASE 22 — ML Anomaly Detection

Then investigate:

- isolation forest,
- autoencoder,
- temporal model,
- sequence anomaly detection.

Compare ML against classical methods.

Do not assume ML is superior.

Measure:

```text
precision
recall
false positive rate
detection latency
```

---

# PHASE 23 — Monte Carlo Validation

Run many randomized missions.

Examples:

```text
100
500
1000
```

Randomize:

- initial orbital state,
- target tumble,
- sensor noise,
- bias,
- spacecraft mass,
- thruster error,
- camera noise,
- timing variation.

Compute:

```text
docking success rate
abort rate
collision rate
mean docking time
propellant usage
position error
velocity error
attitude error
failure-detection rate
```

Example summary:

```text
Docking success       94.7 %
Aborted                4.1 %
Collision              1.2 %

Mean docking time      342 s
Mean propellant        0.83 kg
95% position error     0.14 m
```

---

# PHASE 24 — External Validation

After foundational algorithms are understood, compare selected results against an established framework such as Basilisk.

Possible comparisons:

```text
orbit propagation
rigid-body attitude propagation
relative motion
spacecraft dynamics
```

Architecture:

```text
AstraDock simulation
        ↓
same initial condition
        ↑
reference simulation
        ↓
difference analysis
```

Do not merely replace AstraDock code with an external framework.

The purpose is verification.

---

# 9. Scenario System

Eventually support human-readable scenarios.

Example:

```yaml
mission:
  name: nominal_rendezvous
  seed: 42

earth:
  gravity_model: two_body

target:
  altitude_m: 500000
  initial_tumble_rad_s:
    x: 0.02
    y: -0.01
    z: 0.015

chaser:
  relative_position_m:
    x: -5000
    y: 0
    z: 0

simulation:
  dt_s: 0.1
  duration_s: 1800
```

Do not build the complete schema immediately.

Introduce fields as features are implemented.

---

# 10. Visualization Goals

Python visualization should eventually support:

## Orbit view

```text
Earth
target orbit
chaser orbit
current positions
```

## Telemetry

```text
Mission Time
Altitude
Velocity
Position ECI
Relative Distance
Closing Velocity
Attitude
Angular Rate
Mission Mode
```

## Estimation

Plot:

```text
truth vs measurement vs estimate
position error
velocity error
attitude error
covariance
innovation
```

## Control

Plot:

```text
attitude error
control torque
thruster command
settling time
```

## Rendezvous

Plot:

```text
relative trajectory in LVLH
hold points
keep-out zone
docking corridor
```

---

# 11. Testing Philosophy

Tests are mandatory.

Every mathematical concept must include analytical or invariant-based validation.

Examples:

## Vector

```text
a × b perpendicular to a
a × b perpendicular to b
a · b known result
```

## Orbit

For two-body propagation:

```text
specific orbital energy approximately conserved
specific angular momentum approximately conserved
```

## Quaternion

```text
norm remains 1
rotation preserves vector length
q * q^-1 = identity
```

## Coordinate transformation

```text
A → B → A returns original vector
```

## RK4

Compare with an ODE having a known analytical solution.

## Kalman filter

Use deterministic synthetic system with known truth.

## Controller

Verify expected convergence from known initial state.

---

# 12. Verification Levels

Each major subsystem should eventually be verified at several levels:

```text
Level 1 — unit tests
Level 2 — analytical comparison
Level 3 — conservation/invariant checks
Level 4 — scenario-level tests
Level 5 — Monte Carlo testing
Level 6 — external reference comparison
```

---

# 13. Documentation Rules

Every major module must document:

```text
What problem does it solve?
What coordinate frame does it use?
What units does it expect?
What mathematical model does it implement?
What assumptions does it make?
How was it tested?
What is not modeled?
```

Do not write meaningless comments.

Bad:

```cpp
// update velocity
velocity += acceleration * dt;
```

Good documentation should explain engineering intent or constraints.

---

# 14. README Goals

README should eventually explain:

1. what AstraDock is,
2. why it exists,
3. current capabilities,
4. architecture,
5. screenshots/plots,
6. how to build,
7. how to run scenarios,
8. verification results,
9. learning curriculum,
10. limitations,
11. roadmap.

Clearly distinguish:

```text
IMPLEMENTED
IN PROGRESS
PLANNED
```

Never claim planned functionality exists.

---

# 15. AGENTS.md Requirements

Create a root `AGENTS.md` containing the following behavioral principles.

Codex must follow them for the entire repository.

Include wording equivalent to:

```text
This repository is both an engineering project and an aerospace
learning environment.

Never implement a new mathematical, astrodynamics, GNC, estimation,
control, computer-vision, or machine-learning concept without first
explaining the physical problem and documenting the governing model.

For major concepts document:

- physical problem
- coordinate frames
- units
- state variables
- governing equations
- assumptions
- numerical method
- verification strategy

Do not hide foundational physics behind third-party aerospace libraries.

During the learning phases implement foundational algorithms ourselves.

External libraries may later be used for validation.

Never introduce machine learning where a deterministic classical method
should be understood first.

Examples:

- EKF before neural state estimation
- PID/LQR before reinforcement-learning control
- analytical orbital mechanics before learned trajectory models
- geometric pose estimation before learned end-to-end pose estimation

All simulations must be deterministic when given a random seed.

SI units must be used internally unless explicitly documented otherwise.

Every major simulation must expose truth state, estimated state,
measurement state, and command state separately.

Before implementing a major subsystem:
1. inspect existing code,
2. inspect relevant lessons and architecture,
3. write a concise implementation plan,
4. identify verification criteria,
5. implement,
6. compile,
7. run tests,
8. inspect failures,
9. fix failures,
10. summarize what changed and what was learned.

Do not optimize for maximum code output.

Optimize for correctness, testability, physical clarity, maintainability,
and educational value.
```

---

# 16. Roadmap Milestones

Create `ROADMAP.md` with sequential milestones.

Recommended order:

```text
M00 repository foundation
M01 vectors and math primitives
M02 two-body gravity
M03 numerical integration
M04 orbital propagation
M05 orbital invariants
M06 coordinate frames
M07 quaternions
M08 attitude propagation
M09 rigid-body dynamics
M10 sensor models
M11 linear Kalman filter
M12 extended Kalman filter
M13 attitude estimation
M14 second spacecraft
M15 relative motion / Hill frame
M16 rendezvous guidance
M17 translational control
M18 attitude control
M19 docking corridor
M20 mission state machine
M21 virtual camera
M22 classical visual pose estimation
M23 ML pose estimation
M24 fault injection
M25 classical fault detection
M26 ML anomaly detection
M27 Monte Carlo validation
M28 external dynamics validation
M29 final portfolio documentation
```

Every milestone must have:

```text
objective
concepts learned
deliverables
tests
completion criteria
```

---

# 17. CURRENT TASK — START HERE

Do **not** build the entire project immediately.

For the first Codex run, implement only the repository and mathematical foundation.

## Task 1 Goals

### A. Inspect Environment

Before editing files:

- inspect operating system,
- compiler availability,
- CMake version,
- Python version,
- Git state,
- directory contents.

Report blockers only if they genuinely prevent work.

Use reasonable fallbacks when possible.

---

### B. Write Implementation Plan

Before modifying files, write a concise plan covering:

- repository setup,
- C++ math library,
- tests,
- Python tooling,
- documentation.

Then execute the plan.

Do not wait for user confirmation unless destructive action would be required.

---

### C. Establish Repository

Create only what is currently useful.

Minimum expected files:

```text
AGENTS.md
README.md
ROADMAP.md
CMakeLists.txt
pyproject.toml
docs/architecture.md
docs/curriculum.md
docs/lessons/001_vectors.md
```

Plus source/test files needed for Vector3.

---

### D. C++ Foundation

Use C++20.

Implement:

```text
Vector3
```

Required functionality:

```text
constructor
x/y/z access
addition
subtraction
scalar multiplication
scalar division
dot product
cross product
squared norm
norm
normalization
approximate equality helper if appropriate
```

Avoid unnecessary templates.

Prefer readable double-precision implementation initially.

---

### E. Constants

Create a minimal constants location.

Only include constants needed now or imminently.

Potential examples:

```text
pi
```

Do not populate large aerospace constant catalogs yet.

---

### F. Unit Convention

Document internal SI convention.

No unit framework is needed yet unless justified.

---

### G. Tests

Add deterministic C++ tests for Vector3.

Must include:

1. addition,
2. subtraction,
3. scalar multiplication,
4. scalar division,
5. dot product,
6. cross product,
7. norm,
8. normalization,
9. cross-product orthogonality,
10. zero-vector normalization behavior.

Clearly choose and document expected zero-normalization behavior.

---

### H. Python Foundation

Create minimal Python project tooling for future analysis.

Do not add heavy dependencies.

At most establish:

```text
pytest
ruff
```

if appropriate.

No ML libraries yet.

---

### I. Documentation

Create:

```text
docs/architecture.md
```

Describe:

- current architecture,
- future architecture,
- separation between simulation core and analysis,
- truth / measurement / estimate / command concept.

Clearly label future functionality as planned.

Create:

```text
docs/curriculum.md
```

Map each roadmap stage to the concepts learned.

Create:

```text
docs/lessons/001_vectors.md
```

The lesson must explain:

- what a vector represents physically,
- magnitude,
- direction,
- dot product,
- cross product,
- why vectors matter in spacecraft engineering,
- examples using position, velocity, acceleration, force, torque,
- frame dependence,
- units,
- common mistakes.

---

### J. Build Verification

After implementation:

1. configure CMake,
2. compile,
3. run all C++ tests,
4. run Python checks if configured,
5. inspect warnings/errors,
6. fix project-caused failures,
7. run tests again.

---

### K. Final Report

At the end of this first task report:

```text
files created
files modified
build command
test command
tests passed
important design decisions
concepts introduced
known limitations
next recommended milestone
```

Do not implement orbital propagation in Task 1.

---

# 18. NEXT TASK AFTER TASK 1

Once Task 1 is clean and tests pass, proceed to:

## Task 2 — Two-Body Gravity

Before implementation, teach and document:

```text
Newton's law of gravitation
gravitational parameter mu
state vector
position vector
velocity vector
acceleration vector
Earth-centred inertial frame
circular orbital velocity
specific orbital energy
specific angular momentum
```

Implement:

```text
two_body_acceleration(position, mu)
```

Validation must include:

```text
known radius
analytical acceleration magnitude
direction toward central body
circular-orbit velocity
```

Create:

```text
docs/lessons/002_two_body_dynamics.md
```

Do not implement full RK4 propagation until gravity unit tests pass.

---

# 19. Codex Working Protocol for All Future Tasks

For each milestone:

## Before coding

Codex must:

1. inspect repository state,
2. inspect relevant previous implementations,
3. inspect architecture documentation,
4. state the physical objective,
5. define frames,
6. define units,
7. define state variables,
8. list equations,
9. list assumptions,
10. list tests.

---

## During coding

Codex must:

- make small coherent changes,
- avoid unrelated refactors,
- keep APIs explicit,
- use clear names,
- add tests alongside code,
- keep truth and estimates separate,
- avoid premature optimization,
- avoid premature abstraction.

---

## After coding

Codex must:

1. compile,
2. run tests,
3. inspect warnings,
4. inspect numerical outputs,
5. fix failures,
6. run tests again,
7. update docs,
8. summarize.

---

# 20. Learning Checkpoints

After each major lesson Codex should include a short section:

```text
What you should now understand
```

containing approximately 3–8 questions.

Example after orbital mechanics:

```text
1. Why is orbital velocity perpendicular to the radius vector in a circular orbit?
2. What does mu represent?
3. Why does the gravitational acceleration scale with 1/r²?
4. Why are position and velocity together required to define orbital state?
5. What quantities should stay constant in an ideal two-body orbit?
```

Do not block implementation waiting for answers.

These are review prompts for the learner.

---

# 21. Code Quality Expectations

Prefer code that another engineer can inspect.

Avoid:

- giant files,
- magic numbers,
- hidden global state,
- ambiguous frames,
- ambiguous units,
- unnecessary inheritance,
- complex metaprogramming,
- premature ECS architectures,
- premature GPU acceleration,
- unnecessary service architecture,
- web backends before simulation exists,
- ML before physics.

---

# 22. Error Handling

Use explicit errors for physically invalid operations.

Examples:

```text
normalize zero vector
invalid timestep
negative mass
invalid inertia matrix
NaN state
invalid quaternion
```

Tests should verify important error behavior.

---

# 23. Performance Philosophy

Correctness first.

Expected progression:

```text
correct
↓
verified
↓
clear
↓
profiled
↓
optimized
```

Never reverse this sequence.

---

# 24. Git Practice

Keep work logically separated.

Suggested milestone commits:

```text
feat(math): add Vector3 primitives
feat(dynamics): add two-body gravity
feat(sim): add RK4 propagation
feat(frames): add ECI and LVLH transforms
feat(attitude): add quaternion math
feat(nav): add Kalman filter
```

Codex should not commit unless the environment/user workflow allows it and committing is appropriate.

Never rewrite unrelated history.

---

# 25. Portfolio Goal

The final repository should allow a recruiter or engineer to see that the learner understands:

```text
C++
Python
numerical methods
orbital mechanics
spacecraft dynamics
coordinate frames
quaternions
sensor modeling
Kalman filtering
sensor fusion
guidance
control
relative navigation
autonomous rendezvous
computer vision
machine learning
fault detection
Monte Carlo validation
engineering verification
```

The final project should **not** merely say:

> ML graduate interested in space.

It should demonstrate:

> Software/ML engineer who has independently built and verified substantial components of an autonomous spacecraft simulation stack.

---

# 26. Non-Goals

Do not attempt to build high-fidelity versions of these early:

```text
full atmospheric model
high-order gravity
thermal model
structural FEM
real combustion
real flight hardware drivers
radiation model
complete spacecraft bus
high-fidelity CFD
full mission planning system
production flight certification
```

These may be explored much later.

---

# 27. Safety / Realism Boundary

This project is an educational simulation.

Keep work focused on:

- peaceful spacecraft simulation,
- rendezvous,
- docking,
- navigation,
- robotics,
- research,
- Earth observation,
- autonomous space systems.

Do not expand the project into weapon targeting or weapon-delivery optimization.

---

# 28. Definition of Done for First Run

The first Codex run is complete when all of the following are true:

```text
[ ] repository inspected
[ ] implementation plan written
[ ] AGENTS.md created
[ ] README.md created
[ ] ROADMAP.md created
[ ] architecture documented
[ ] curriculum documented
[ ] C++20 build configured
[ ] Vector3 implemented
[ ] Vector3 tests implemented
[ ] vectors lesson written
[ ] build succeeds
[ ] tests pass
[ ] project-caused warnings/errors addressed
[ ] final implementation summary produced
```

Do not continue to orbital mechanics automatically in the same run unless explicitly instructed.

---

# 29. First Command to Codex

Once this file exists in the repository, the user may give Codex this instruction:

```text
Read ASTRADOCK_CODEX_MASTER.md in full.

Treat it as the master specification for this repository.

Then execute ONLY the section titled:

"CURRENT TASK — START HERE"

Follow all engineering, learning, testing, documentation, and verification
rules defined in the master specification.

Do not proceed to Task 2 in this run.

Inspect the environment first, write a concise implementation plan, execute
the task, compile, test, fix project-caused failures, and report the result.
```

---

# 30. Long-Term Principle

The goal is not to let Codex build AstraDock *for* the learner.

The goal is to use Codex to build AstraDock **with** the learner while the repository itself becomes a structured aerospace engineering education.

Every milestone should leave the learner with:

```text
new mathematics understood
+
new engineering code implemented
+
tests proving it works
+
visual evidence
+
documentation explaining why
```

That principle takes priority over speed.
