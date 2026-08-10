# AstraDock Architecture

## Purpose and scope

AstraDock separates a deterministic C++ simulation core from Python analysis
tools. The current implementation contains the `Vector3` mathematical
foundation and one instantaneous physical model: central two-body gravitational
acceleration. It also contains reusable Forward Euler and classical RK4
single-step ODE integrators. M04 connects them through a Cartesian orbital
state, two-body state derivative, fixed-step propagation driver, diagnostic
helpers, CSV demonstration, and Python plotting layer.

## Current architecture

```text
C++ simulation core

Vector3 + Earth constants
          |
          +--> two_body_acceleration
          |             |
CartesianState -------->+--> two_body_state_derivative
                                      |
Euler / RK4 --------------------------+
          |
          v
propagate_fixed_step
          |
          v
time-stamped Cartesian trajectory
          |
          +--> energy / angular-momentum diagnostics
          |
          v
astradock_orbit_demo --> deterministic CSV
                                  |
                                  v
Python analysis            plot_orbit.py --> PNG figures
```

`AstraDock::math`, `AstraDock::dynamics`, `AstraDock::numerics`, and
`AstraDock::orbit` are header-only CMake interface targets requiring C++20.
Orbit depends on dynamics and numerics; dynamics and numerics each depend on
math. The `astradock_orbit_demo` executable links only the C++ core and exports
plain CSV. CTest owns test execution, while Catch2 provides test cases and
assertions.

| Area | Current responsibility | Current non-responsibility |
| --- | --- | --- |
| `cpp/math` | General double-precision vector operations and minimal mathematical/Earth constants | Frames, unit types, matrices, dynamics |
| `cpp/dynamics` | Instantaneous point-mass central gravitational acceleration | Time integration, orbit propagation, perturbations |
| `cpp/numerics` | Generic Euler/RK4 steps plus deterministic signed fixed-step propagation and endpoint policy | Physical derivative assembly, adaptive steps, events |
| `cpp/orbit` | Cartesian orbital state, two-body derivative assembly, circular-orbit helpers, and invariants | Frame transforms, perturbations, orbital elements |
| `tools` | Nominal 500 km C++ scenario, quantitative summaries, and deterministic CSV export | General scenario or telemetry framework |
| `python/analysis` | Plotting and display-unit conversion from authoritative C++ CSV | Gravity, integration, or independent propagation |
| `tests/cpp` | Deterministic analytical, invariant, ODE convergence, endpoint, and orbit scenario tests | Higher-fidelity external reference validation |
| `docs` | Architecture, curriculum, and concept lessons | Claims of unimplemented behavior |
| `pyproject.toml` | Python 3.12+ tooling and optional Matplotlib analysis dependency | Python simulation physics |

## C++ simulation core and Python analysis

The C++ core will own physical models and time evolution because it is the
authoritative implementation being learned and verified. Planned C++ areas
include math, dynamics, frames, sensors, navigation, guidance, control, faults,
and simulation orchestration.

Python now consumes exported C++ orbit results for plotting. Future Python work
may add diagnostics, Monte Carlo post-processing, notebooks, and machine
learning. Python analysis must not silently reimplement a different physical
model. When Python contains a reference calculation, its role and tolerance
must be explicit.

The planned dependency direction is:

```text
scenario inputs
      |
      v
C++ simulation core ---> deterministic telemetry/results ---> Python analysis
      ^                                                        |
      |                                                        v
tests and analytical references                       plots and reports
```

## Planned state separation

Future simulation interfaces will keep four kinds of state distinct:

| State kind | Meaning | Example planned field |
| --- | --- | --- |
| Truth | The simulator's physical state | `truth.position_eci_m` |
| Measurement | A sensor's sampled, imperfect output | `measurement.gnss_position_eci_m` |
| Estimate | Navigation's belief based on models and measurements | `estimate.position_eci_m` |
| Command | Guidance or control's desired state/action | `command.force_body_n` |

The arrows between them are one-way responsibilities, not aliases:

```text
truth --> sensor model --> measurement --> estimator --> estimate
                                                    |
                                                    v
truth dynamics <-- actuators <-- control <-- guidance/command
```

A controller must not read truth in place of an estimate, and a sensor must not
overwrite truth. This separation remains planned. The M04 trajectory contains
only deterministic translational truth samples; no measurement, estimate, or
command state exists yet.

## Units and coordinate frames

Internal physical quantities use SI units: metres, seconds, kilograms, radians,
newtons, and derived SI units. Human-facing plots may convert units if labels
make the conversion explicit.

`Vector3` intentionally represents only three numeric components. The meaning
must be supplied by its containing type, variable name, or function contract.
For example, `position_eci_m` and `force_body_n` may both use `Vector3`, but they
cannot be added because they have different dimensions and frames. Planned
frame transforms will be explicit functions with documented conventions,
inverse tests, and round-trip tests.

The M02 gravity API names its position as central-body-centered inertial and
uses metres, m^3/s^2, and m/s^2 for position, gravitational parameter, and
acceleration. For Earth this is conceptually ECI, but no concrete ECI definition
or frame transformation is implemented yet.

The M03 integration APIs are frame-neutral and unit-neutral. The caller owns
those meanings: when time is measured in seconds, the derivative must use
state-units per second. Finite positive and negative steps are supported, while
a zero step is an identity operation. The integrators validate finite time and
the known `double`/`Vector3` values but do not choose a physical timestep or
error budget for the caller.

M04's `CartesianState` uses position in metres and velocity in m/s in an
idealized Earth-centered inertial Cartesian frame. The derivative's same-shaped
slots represent m/s and m/s^2. Propagation time is seconds. Specific energy and
angular momentum use m^2/s^2 and m^2/s. No operational ECI definition, epoch,
Earth rotation, or coordinate transformation is implemented.

## Planned subsystem architecture

The following structure is a direction, not a claim of implementation:

```text
orbital and rigid-body truth dynamics
          |
          +--> sensor models --> measurements --> state estimation
          |                                      |
          |                                      v
          +<-- actuators <-- control <-- guidance and mission state machine
                                             |
                                             v
                                  safety monitoring and aborts

virtual camera --> classical geometry --> pose measurement --> estimation
```

Machine learning is planned only after classical estimation, control, geometric
vision, and fault detection can serve as understandable baselines.

## Verification architecture

Each subsystem should progress through unit tests, analytical comparisons,
invariant checks, scenario tests, Monte Carlo tests, and finally selected
external-reference comparisons. The current `Vector3` suite covers exact known
cases, approximate floating-point cases, cross-product orthogonality, norm
behavior, and invalid zero normalization.

The two-body suite adds analytical magnitude/direction cases, axis independence,
anti-parallel vector invariants, inverse-square scaling, Earth and LEO sanity
checks, and explicit invalid-input behavior.

The numerical integration suite adds one-step hand calculations, scalar and
vector ODEs with analytical solutions, first- and fourth-order convergence
trends, equal-step accuracy comparison, constant derivatives, backward and
zero-step policies, and non-finite-value rejection. It does not include the
gravity header or propagate a spacecraft state.

The M04 suite adds Cartesian state algebra, derivative assembly against the
canonical gravity function, exact/shortened/zero/backward endpoint behavior,
invalid propagation parameters, short-duration sanity, analytical one-period
closure, radius stability, energy and angular-momentum conservation, and a
five-orbit Euler/RK4 drift comparison. Plot generation is checked separately
from core numerical unit tests.

## Current assumptions and limitations

- Arithmetic uses IEEE 754 `double` values.
- No frame or unit type system exists yet.
- `Vector3::normalized()` rejects exactly zero magnitude and does not invent a
  direction.
- The only physical model is spherical point-mass Earth gravity.
- Euler and RK4 are explicit, fixed-step, single-step primitives with no
  adaptive error control, event detection, or stiffness detection.
- The propagation driver stores all samples and provides no maneuver, event,
  general scenario, or ephemeris system.
- There is no coordinate transformation, orbital-element representation,
  perturbation model, attitude state, concurrency, or external astrodynamics
  dependency.
- Catch2 is a test-only dependency and is not part of the simulation API.
