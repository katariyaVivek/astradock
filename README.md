# AstraDock

AstraDock is an educational C++20 simulator for learning how autonomous
spacecraft rendezvous and docking systems are built and verified. The long-term
project will cover orbital and rigid-body dynamics, sensing, estimation,
guidance, control, docking vision, fault handling, and Monte Carlo validation.
It begins with transparent mathematical foundations rather than hiding the
physics behind an astrodynamics framework.

## Project status

### Implemented

- A C++20 header-only `Vector3` primitive using double precision.
- Addition, subtraction, scalar arithmetic, dot and cross products, norms,
  normalization, and approximate comparison.
- Explicit errors for division by zero and normalization of the zero vector.
- A minimal constants header containing `pi`, the WGS 84 Earth gravitational
  parameter, and the WGS 84 Earth reference radius.
- Central-body two-body gravitational acceleration in a documented
  central-body-centered inertial frame, using SI units.
- Explicit rejection of the gravity singularity at zero position and invalid
  physical parameters.
- Reusable fixed-step Forward Euler and classical RK4 integration primitives
  for scalar and `Vector3` states.
- Documented support for zero and negative finite time steps, with explicit
  rejection of non-finite M03 numerical inputs and results.
- A six-component Cartesian translational orbital state in a simplified
  Earth-centered inertial frame.
- A two-body state derivative that reuses the canonical M02 gravity function.
- Deterministic signed fixed-step propagation with Euler/RK4 selection,
  endpoint sampling, shortened final steps, and backward-time support.
- Circular speed and period helpers plus specific orbital energy and angular-
  momentum diagnostics.
- A C++ 500 km circular-orbit demo that exports one- and five-orbit CSV data.
- Python plots for the RK4 trajectory, Euler/RK4 altitude, energy drift, and
  angular-momentum drift.
- Deterministic Catch2 tests registered with CTest.
- Architecture, curriculum, roadmap, and lessons for vectors, two-body
  dynamics, numerical integration, orbital propagation, and propagator
  validation.
- Minimal Python 3.12+ metadata and optional `pytest`/`ruff` development tools.
- M05 validation: systematic timestep-sweep convergence study, empirical
  order measurement, phase error, period estimation, repeatability checks,
  and regression baselines.
- Reusable orbital diagnostics: phase angle, analytical phase, unwrapped
  phase, max phase error, period estimation, angular-momentum direction
  drift, and empirical convergence order.
- M06 Coordinate Frames:
  - `Matrix3` 3x3 double-precision matrix primitive (transposition, determinant, trace, and orthonormality checks).
  - `FrameBasis` triad representation with right-handed orthonormality validation.
  - Orthonormal LVLH (Local Vertical Local Horizontal) basis construction from instantaneous ECI state vectors.
  - Direction Cosine Matrices (`dcm_lvlh_from_eci`, `dcm_eci_from_lvlh`) with explicit row/column conventions.
  - Geometric coordinate transformations (`transform_eci_to_lvlh`, `transform_lvlh_to_eci`) preserving vector norm.
  - Explicit domain error rejection for degenerate states (zero position, collinear position/velocity, non-finite values).
  - C++ frame demonstration tool (`astradock_frame_demo`) exporting full-orbit frame telemetry to CSV.
  - Python analysis plots (`plot_frames.py`) showing single-point frame triads, orbital triad evolution, and vector transformation sanity.
- M07 Classical Orbital Elements:
  - `ClassicalOrbitalElements` representation ($a, e, i, \Omega, \omega, \nu$) in SI units and radians.
  - Perifocal Coordinate Frame ($PQW$) construction and $C_{ECI\_PQW}$ Direction Cosine Matrix.
  - Bidirectional conversions: `state_to_classical_elements` and `classical_elements_to_state` with sub-nanometre round-trip fidelity.
  - Analytical helpers for semi-latus rectum ($p$), periapsis/apoapsis radii ($r_p, r_a$), mean motion ($n$), orbital period ($T$), and orbit classification (`classify_orbit`).
  - Strict $[0, 2\pi)$ angle normalization and angular distance utilities (`cpp/math/angle.hpp`).
  - Singularity handling and defensive conventions for circular, equatorial, and circular-equatorial orbits.
  - Multi-orbit element invariance verification ($a, e, i, \Omega, \omega$ drift $< 10^{-9}$ over multiple orbits).
  - C++ orbital elements demo (`astradock_elements_demo`) exporting full-orbit telemetry to CSV.
- M08 Attitude Representation & Quaternion Mathematics:
  - `Quaternion` representation ($q = [w, x, y, z]$ scalar-first) in double precision.
  - Quaternion algebra: norm, unit check, safe normalization, conjugate ($q^*$), and inverse ($q^{-1}$).
  - Hamilton product ($\otimes$) with non-commutative multiplication and rotation composition ($q_{A\_C} = q_{A\_B} \otimes q_{B\_C}$).
  - Active vector rotation: $v' = q \otimes v \otimes q^* \equiv C(q) v$.
  - Bidirectional conversions between `Quaternion` and `Matrix3` (DCM) using Shepperd's numerically stable algorithm.
  - Double-cover property verification ($+q$ and $-q$ represent the same physical rotation).
  - Secondary `EulerAngles` representation (ZYX yaw-pitch-roll $\psi, \theta, \phi$) with gimbal lock detection at $\theta = \pm 90^\circ$.
  - Spacecraft `AttitudeState` representation (`cpp/attitude/attitude_state.hpp`).
  - Scale-aware vector norm preservation across 15 orders of magnitude ($10^{-9}$ to $10^6$).
  - C++ attitude demonstration tool (`astradock_attitude_demo`) exporting telemetry to CSV.
  - Python 3D body frame, quaternion vs DCM consistency, and composition plots (`plot_attitude.py`).
  - Lesson 008 on Spacecraft Attitude & Quaternions and M08 Validation Report.
- M09 Rigid-Body Attitude Dynamics & Quaternion Kinematics:
  - `PrincipalInertia` representation ($I_{xx}, I_{yy}, I_{zz} > 0\,\text{kg}\cdot\text{m}^2$) with physical positive/finite validation.
  - `RotationalState` 7-component state ($q, \boldsymbol{\omega}$) with ODE vector-space arithmetic operators.
  - Euler's rigid-body equations in principal axes: $\dot{\boldsymbol{\omega}} = \mathbf{I}^{-1} [\boldsymbol{\tau} - \boldsymbol{\omega} \times (\mathbf{I}\boldsymbol{\omega})]$.
  - Quaternion kinematics differential equation: $\dot{q} = \frac{1}{2} q \otimes [0, \boldsymbol{\omega}_B]$.
  - Rotational invariants: rotational kinetic energy ($E_{rot} = \frac{1}{2}\boldsymbol{\omega}^T \mathbf{I}\boldsymbol{\omega}$) and inertial angular momentum ($\mathbf{H}_I = C_{\mathcal{I}\_\mathcal{B}}\mathbf{H}_B$).
  - Step-boundary quaternion normalization policy bounding norm error to $< 10^{-15}$ without distorting RK4 4th-order convergence.
  - C++ attitude dynamics demo tool (`astradock_attitude_dynamics_demo`) exporting telemetry for principal spin, constant torque, and asymmetric tumble.
  - Independent Python reference oracle (`python/audit/independent_attitude_reference.py`) achieving sub-micro-residual cross-validation ($\Delta \omega < 10^{-12}\,\text{rad/s}$, $\Delta q < 10^{-12}$).
  - Python analysis plots (`plot_attitude_dynamics.py`) showing rate evolution, quaternion trajectory, conservation invariants, analytical torque comparisons, and body vs inertial momentum distinction.
  - Lesson 009 on Rigid-Body Attitude Dynamics & Quaternion Kinematics and M09 Validation Report.
- M10 Integrated 6-DOF Spacecraft State:
  - `SpacecraftState` 13-component composite state combining `CartesianState` (position, velocity in ECI) and `RotationalState` (quaternion, body rate).
  - Physical 6-DOF vs 13-component numerical representation with unit-quaternion constraint ($\|q\| = 1$).
  - `SpacecraftParameters` (gravitational parameter $\mu$, principal inertia moments $\mathbf{I}$).
  - `ForceTorqueInput` with explicit frame definitions: force in **ECI** ($\mathbf{F}_{\mathcal{I}}$), torque in **BODY** ($\boldsymbol{\tau}_B$).
  - Unified first-order derivative `spacecraft_state_derivative` reusing canonical M02 gravity and M09 rigid-body rotational dynamics.
  - Geodesic orientation error metric `quaternion_orientation_error_rad` invariant to quaternion double cover.
  - Fixed-step RK4 propagation engine `rk4_step_spacecraft` and `propagate_spacecraft_fixed_step` with post-step quaternion normalization.
  - Canonical 500 km circular orbit + asymmetric tumbling simulation verifying simultaneous conservation of all 5 orbital and rotational invariants.
  - Verified physical decoupling of translation and attitude in central two-body gravity.
  - C++ 6-DOF demo tool (`astradock_6dof_demo`) exporting telemetry to CSV.
  - Independent pure-Python reference oracle (`python/audit/independent_6dof_reference.py`) achieving sub-micrometre orbital and sub-micro-radian attitude cross-validation.
  - Python analysis plots (`plot_6dof.py`) generating 3D orbit + body triads, translational telemetry, rotational invariants, and subsystem decoupling comparisons.
  - Lesson 010 on Integrated 6-DOF Spacecraft State and M10 Validation Report.
- M11 Spacecraft Environment & Force/Torque Models:
  - M10 Numerical Entry Gate re-validation confirming rotational drift is pure $\mathcal{O}(\Delta t^4)$ truncation error.
  - Earth Oblateness ($J_2$) gravity perturbation acceleration in ECI (`j2_acceleration_eci`), verifying secular RAAN regression $\dot{\Omega}$ and apsidal precession $\dot{\omega}$.
  - Atmospheric Drag model with rotating atmosphere (`relative_atmospheric_velocity_eci`), exponential density (`exponential_atmospheric_density`), and drag acceleration (`drag_acceleration_eci`), verifying orbital decay and monotonic energy dissipation.
  - Third-Body Lunar/Solar gravitational tidal acceleration (`third_body_acceleration_eci`) with direct and indirect tidal terms.
  - Gravity-Gradient body torque (`gravity_gradient_torque_body`) on asymmetric spacecraft, verifying spherical symmetry nulls, principal alignment, and restorative pitch libration.
  - Modular `EnvironmentConfiguration` switches with guaranteed bitwise regression to M10 when disabled.
  - Unified environmental 6-DOF propagator `propagate_spacecraft_environmental`.
  - C++ environment demo tool (`astradock_environment_demo`) exporting telemetry to CSV.
  - Independent pure-Python reference oracle (`python/audit/independent_environment_reference.py`) achieving $< 10^{-17}\,\text{m/s}^2$ cross-validation.
- M12 Spacecraft Sensor Simulation & Measurement Models:
  - `DeterministicRng` wrapping 64-bit Mersenne Twister (`std::mt19937_64`) for reproducible Gaussian noise and uniform spherical vector sampling.
  - `SensorSchedule` managing multi-rate asynchronous update schedules and `DropoutWindow` handling deterministic sensor failures (`valid = false`).
  - 6-Axis IMU (`ImuSensor`, `ImuConfig`, `ImuMeasurement`) at 100 Hz:
    - Gyroscope measuring body rates with constant bias and Gaussian noise.
    - Accelerometer measuring specific force $\mathbf{f} = \mathbf{a} - \mathbf{g}$ in Spacecraft BODY frame, demonstrating near-zero free-fall contact acceleration in orbit and nonzero response to non-gravitational drag.
  - Spaceborne GNSS Receiver (`GnssSensor`, `GnssConfig`, `GnssMeasurement`) at 1 Hz measuring absolute ECI position and velocity.
  - Optical Star Tracker (`StarTrackerSensor`, `StarTrackerConfig`, `StarTrackerMeasurement`) at 10 Hz with proper physical $SO(3)$ rotation perturbations and exact unit norm preservation.
  - Line-of-Sight Relative Range Sensor (`RangeSensor`, `RangeSensorConfig`, `RangeMeasurement`) at 10 Hz with translation invariance.
  - Truth Non-Interference regression: active sensors leave simulation truth states 100% bitwise identical.
  - C++ sensor demo tool (`astradock_sensor_demo`) exporting telemetry to CSV.
  - Independent pure-Python reference oracle (`python/audit/independent_sensor_reference.py`) cross-verifying telemetry statistics and physical invariants.
  - Python analysis plots (`plot_sensors.py`) generating IMU, GNSS, Star Tracker, Range, and multi-rate timeline figures.
  - Lesson 012 on Sensor Simulation & Measurement Models and M12 Validation Report.
- M13 Spacecraft State Estimation & Integrated Multi-Rate Navigation (COMPLETE):
  - Consolidated navigation architecture uniting translational motion, attitude on $S^3$, and dynamic sensor biases.
  - Unified 15-state estimation error representation: $\delta\mathbf{x} = [\delta\mathbf{r}^T, \delta\mathbf{v}^T, \delta\boldsymbol{\theta}^T, \delta\mathbf{b}_a^T, \delta\mathbf{b}_g^T]^T \in \mathbb{R}^{15}$.
  - Analytical continuous error dynamics Jacobian $F \in \mathbb{R}^{15 \times 15}$ incorporating physical cross-coupling blocks:
    - Attitude-to-velocity coupling: $\partial \delta\dot{\mathbf{v}} / \partial \delta\boldsymbol{\theta} = -C_I^B(\hat{\mathbf{q}}) [\hat{\mathbf{f}}_B]_\times$ (attitude errors rotate specific force into false inertial acceleration).
    - Accelerometer-bias-to-velocity coupling: $\partial \delta\dot{\mathbf{v}} / \partial \delta\mathbf{b}_a = -C_I^B(\hat{\mathbf{q}})$.
    - Analytical gravity gradient tensor: $G(\hat{\mathbf{r}}) = -\frac{\mu}{r^3}(I_3 - 3\hat{\mathbf{r}}\hat{\mathbf{r}}^T)$.
  - Multiplicative Extended Kalman Filter (MEKF) on $S^3$ using body-frame error vector $\delta\boldsymbol{\theta}$, rate-gyro bias $\delta\mathbf{b}_g$, and antipodal sign alignment $\text{sign}(\mathbf{q}_m \cdot \hat{\mathbf{q}})$.
  - Asynchronous multi-rate sequential updates (100 Hz IMU, 10 Hz star tracker, 10 Hz range, 1 Hz GNSS) using actual timestamp intervals ($\Delta t$).
  - Joseph-form covariance updates on unified $15 \times 15$ covariance $P$ with first-order covariance reset $J_{\text{reset}}$.
  - Scalar nonlinear relative range updates with analytical 1x15 line-of-sight Jacobian $H_\rho$ and coincident singularity rejection ($\rho < 10^{-6}$ m).
  - C++ integrated navigation demo tool (`astradock_integrated_navigation_demo`) exporting multi-rate telemetry, sensor combination comparisons, and 100-seed Monte Carlo results to CSV.
  - Independent pure-Python reference oracle (`python/audit/independent_integrated_navigation_reference.py`) verifying all operations and finite-difference Jacobians ($< 10^{-10}$).
  - Comprehensive visualization suite (`plot_integrated_navigation.py`, Plots A–H in `artifacts/figures/`).
  - Canonical Lesson 013 (`docs/lessons/013_state_estimation_and_navigation.md`) and consolidated validation report (`docs/validation/m13_state_estimation_validation.md`).

### Verification status

M00-M13 were configured, built, and verified with CMake, MSVC, and Clang C++20.
CTest reports all **223 / 223 registered test cases passing (100%)** across 19 test suites:
- 74 from M01-M05 (math, two-body, integrators, orbits, validation)
- 9 Matrix3 tests
- 8 Coordinate Frame tests
- 9 M01-M06 audit tests
- 7 Classical Orbital Elements tests
- 9 Quaternion/Attitude tests
- 14 Attitude Dynamics tests
- 11 6-DOF tests
- 9 Environment tests
- 9 Sensor tests
- 19 State Estimation (M13A) tests
- 19 IMU-Aided Navigation (M13B) tests
- 10 Attitude Error-State EKF (M13C) tests
- 6 Nonlinear Range Update (M13C) tests
- 10 Integrated 15-State Navigation (M13D) tests

All simulations are bitwise deterministic given a seed. All independent Python oracles pass 100%. Ruff passed with zero linter errors.

### Planned

Milestone M14 (actuator dynamics & modeling: reaction wheels and thrusters) is the next scheduled milestone in AstraDock. Actuators, guidance, control, rendezvous, docking vision, fault handling, machine learning, and mission-level Monte Carlo remain planned.

## Repository layout

```text
cpp/math/             C++ mathematical primitives (Vector3, Matrix3, Quaternion, EulerAngles, angle, constants)
cpp/dynamics/         instantaneous physical dynamics models
cpp/numerics/         reusable fixed-step ODE integration primitives
cpp/orbit/            orbital Cartesian states, diagnostics, classical elements, and perifocal frame
cpp/frames/           orthonormal frame bases, LVLH frames, and Direction Cosine Matrices
cpp/attitude/         attitude quaternions, principal inertia, rotational states, and rigid-body dynamics
cpp/spacecraft/       composite 6-DOF spacecraft states, parameters, force/torque inputs, and unified propagator
cpp/environment/      orbital perturbations (J2, atmospheric drag, third-body gravity, gravity gradient)
cpp/sensors/          sensor simulation models (IMU, GNSS, Star Tracker, Range, DeterministicRng, schedules)
cpp/estimation/       Kalman filter algebra, gravity Jacobian, translational EKF (GNSS-only and IMU-aided prediction)
tools/                deterministic C++ demonstration, validation, frame, elements, attitude, dynamics, 6-DOF, environment, sensor, EKF, and IMU-EKF demo executables
python/analysis/      CSV-driven plotting, analysis, and validation plots
python/audit/         independent pure-Python astrodynamics, orbital elements, quaternion, attitude, 6-DOF, environment, sensor, EKF, and IMU-EKF audit oracles
tests/cpp/            deterministic C++ Catch2 unit, validation, coordinate frame, audit, elements, quaternion, dynamics, 6-DOF, environment, sensor, and estimation tests
docs/lessons/         engineering lessons that precede implementations
docs/validation/      numerical validation reports and verification matrices
docs/architecture.md  current and planned system architecture
docs/curriculum.md    learning outcomes by milestone
ROADMAP.md             sequential project milestones
pyproject.toml         Python analysis and linting configuration
```

## Build and test

Prerequisites:

- CMake 3.24 or newer,
- a C++20 compiler,
- Git and network access on the first configure if Catch2 3 is not already
  installed.

Catch2 3.7.1 is pinned and fetched by CMake only when an installed Catch2 3
package cannot be found. Catch2 was chosen over GoogleTest because its concise
test-case style is well suited to analytical math examples while still
integrating cleanly with CTest.

Configure, build, and run all C++ tests:

```powershell
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build --config Debug
ctest --test-dir build -C Debug --output-on-failure
```

For a single-configuration generator, omit `-C Debug` from the CTest command.

## Run the 500 km orbit demo

After building, run the executable path produced by your generator. For a
single-configuration build:

```powershell
.\build\astradock_orbit_demo.exe artifacts\data
py -3.12 python\analysis\plot_orbit.py `
    artifacts\data\orbit_500km_one_orbit.csv `
    artifacts\data\orbit_500km_five_orbits.csv `
    artifacts
```

The C++ executable is authoritative: it propagates both integrators and writes
state plus invariant diagnostics. Python consumes those CSV rows and performs
only display-unit conversion and plotting. Generated `artifacts/` are ignored
by Git.

## Python development tools

Python now provides CSV-driven orbit plots and is reserved for future analysis,
scenarios, and validation; it does not implement simulation physics. Matplotlib
is an optional analysis dependency. To prepare analysis and development tools
with Python 3.12 or newer:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python -m pip install -e ".[analysis,dev]"
.\.venv\Scripts\ruff check .
.\.venv\Scripts\pytest
```

## Engineering conventions

- SI units are used internally unless a documented interface says otherwise.
- Every physical vector must be associated with a coordinate frame at the API
  boundary; `Vector3` itself intentionally does not guess a frame.
- Future simulations will keep truth, measurement, estimate, and command state
  separate.
- Major mathematical or aerospace concepts receive a lesson and verification
  plan before implementation.
- Stochastic simulations will require an explicit seed.

Start with [the vector lesson](docs/lessons/001_vectors.md), continue to the
[two-body dynamics lesson](docs/lessons/002_two_body_dynamics.md), then read
the [numerical integration lesson](docs/lessons/003_numerical_integration.md).
Continue with the
[orbital propagation lesson](docs/lessons/004_orbital_propagation.md) and the
[propagator validation lesson](docs/lessons/005_propagator_validation.md). The
[validation report](docs/validation/m05_propagator_validation.md) quantifies
RK4/Euler convergence. The
[architecture](docs/architecture.md), [curriculum](docs/curriculum.md), and
[roadmap](ROADMAP.md) define the larger sequence.

## Current limitations

- The fixed-step integrators provide no adaptive error control, event handling,
  dense output, or stiffness detection.
- Propagation stores every integration sample in memory and provides no event,
  maneuver, ephemeris, or general scenario system.
- The Earth-centered inertial frame is idealized; no epoch, Earth rotation,
  ECI/ECEF conversion, or LVLH/Hill transform exists.
- Gravity is a spherical point-mass model without atmosphere, drag, J2,
  third-body gravity, or other perturbations.
- `Vector3` does not encode units or frames in its C++ type. This keeps the first
  primitive readable, but APIs must document both explicitly.
- Only an exactly zero vector is rejected by `normalized()`; callers working
  with poorly scaled data must decide and document a suitable near-zero policy.
- The current test dependency may require a one-time network fetch.
