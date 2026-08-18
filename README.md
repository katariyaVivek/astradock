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

### Verification status

M00-M08 were configured and built successfully with CMake, Ninja,
and Clang 21.1.0 through Zig's C++ frontend. CTest reported all 116 registered
test cases passing (74 from M01-M05, 9 Matrix3 tests, 8 Coordinate Frame tests,
9 M01-M06 audit tests, 7 Classical Orbital Elements tests, 9 Quaternion/Attitude tests).
Ruff passed with zero linter errors. All simulations, frame transformations,
orbital element conversions, and quaternion attitude mathematics are deterministic.

### Planned

Rigid body rotational dynamics, angular velocity integration ($\dot{q} = \frac{1}{2} q \otimes \omega$),
inertia tensors, torque models, sensors, estimation, guidance, control, rendezvous/docking,
vision, fault handling, machine learning, and Monte Carlo analysis remain planned.
M09 — Attitude Propagation & Quaternion Kinematics is the next recommended milestone.

## Repository layout

```text
cpp/math/             C++ mathematical primitives (Vector3, Matrix3, Quaternion, EulerAngles, angle, constants)
cpp/dynamics/         instantaneous physical dynamics models
cpp/numerics/         reusable fixed-step ODE integration primitives
cpp/orbit/            Cartesian orbital state, derivative, diagnostics, classical elements, and validation tools
cpp/frames/           coordinate frame bases, LVLH, Direction Cosine Matrices, and transformations
cpp/attitude/         spacecraft attitude state representation (pure orientation)
tools/                deterministic C++ demonstration, validation, frame, elements, and attitude demo executables
python/analysis/      CSV-driven plotting, analysis, and validation plots
python/audit/         independent pure-Python astrodynamics, orbital elements, and quaternion audit oracles
tests/cpp/            deterministic C++ unit, validation, coordinate frame, audit, elements, and quaternion tests
docs/lessons/         engineering lessons that precede implementations
docs/validation/      numerical validation reports and verification matrices
docs/architecture.md  current and planned system architecture
docs/curriculum.md    learning outcomes by milestone
ROADMAP.md             sequential project milestones
pyproject.toml         future Python analysis-tool configuration
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
