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
  dynamics, and numerical integration.
- Minimal Python 3.12+ metadata and optional `pytest`/`ruff` development tools.

### Verification status

M00-M04 were configured and built successfully with CMake 4.4.2, Ninja 1.13,
and Clang 21.1.0 through Zig's C++ frontend. CTest reported all 49 registered
test cases passing. The host did not have a persistent C++ compiler or CMake
installation, so verification used isolated portable tools under the system
temporary directory; future local builds still require the prerequisites below.

### Planned

ECI/ECEF conversion, LVLH/Hill coordinates, orbital elements, perturbations,
J2, drag, attitude dynamics, quaternions, sensors, estimation, guidance,
control, rendezvous/docking, vision, fault handling, machine learning, and
Monte Carlo analysis remain planned. The current orbit is an educational,
spherical point-mass Earth simulation rather than a high-fidelity ephemeris.

## Repository layout

```text
cpp/math/             C++ mathematical primitives
cpp/dynamics/         instantaneous physical dynamics models
cpp/numerics/         reusable fixed-step ODE integration primitives
cpp/orbit/            Cartesian orbital state, derivative, and diagnostics
tools/                deterministic C++ demonstration executables
python/analysis/      CSV-driven plotting and analysis
tests/cpp/            deterministic C++ unit tests
docs/lessons/         engineering lessons that precede implementations
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
[two-body dynamics lesson](docs/lessons/002_two_body_dynamics.md), and then read
the [numerical integration lesson](docs/lessons/003_numerical_integration.md).
Continue with the
[orbital propagation lesson](docs/lessons/004_orbital_propagation.md). The
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
