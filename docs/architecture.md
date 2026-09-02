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

Vector3 / Matrix3 / Quaternion / EulerAngles + Constants
          |
          +--> two_body_acceleration
          |             |
CartesianState -------->+--> two_body_state_derivative
     ^                                |
     | Euler / RK4 -------------------+
     |        |
     |        v
     |   propagate_fixed_step
     |        |
     |        v
     |   time-stamped Cartesian trajectory
     |        |
     |        +--> energy / angular-momentum diagnostics
     |        +--> phase error, period estimation, convergence order
     |        +--> LVLH FrameBasis & DCMs (AstraDock::frames)
     |        +--> ClassicalOrbitalElements & Perifocal PQW (AstraDock::orbit)
     |        +--> AttitudeState & Quaternion Rotations (AstraDock::attitude)
     |        +--> PrincipalInertia, RotationalState & Euler Dynamics (AstraDock::attitude)
     |
     v
`astradock_orbit_demo -------------> deterministic trajectory CSV
astradock_validation_demo --------> convergence CSV + QA
astradock_frame_demo -------------> coordinate frames CSV + QA
astradock_elements_demo ----------> orbital elements propagation CSV + QA
astradock_attitude_demo ----------> attitude rotation telemetry CSV + QA
astradock_attitude_dynamics_demo -> rigid-body dynamics & kinematics CSV + QA
astradock_6dof_demo --------------> integrated 6-DOF telemetry CSV + QA
astradock_environment_demo -------> environmental perturbation CSV + QA
astradock_sensor_demo ------------> multi-rate sensor telemetry CSV + QA
astradock_ekf_demo ---------------> estimator truth/measurement/estimate CSV,
                                    Monte Carlo summary CSV + QA
                                        |
                                        v
Python analysis            plot_orbit.py -------------> trajectory PNG figures
                           plot_validation.py --------> convergence plots
                           plot_frames.py ------------> coordinate frame plots
                           plot_orbital_elements.py --> 3D orbital geometry & element evolution
                           plot_attitude.py ----------> 3D body frame, quaternion/DCM & composition
                           plot_attitude_dynamics.py -> rates, quaternions, energy/momentum invariants
                           plot_estimation.py ---------> truth vs GNSS vs estimate, covariance, innovation/NIS

Frame & Attitude Hierarchy:
ECI Frame (Inertial Reference)
 │
 ├── Orbital State (Cartesian r, v / Classical Elements a, e, i, Ω, ω, ν)
 │
 └── LVLH Frame (Orbital Local)
       │
       └── Spacecraft Attitude & Rotational Dynamics
              │
              ▼
          SPACECRAFT BODY FRAME (Principal axes: I = diag(Ixx, Iyy, Izz), omega_B, tau_B, q_ECI_Body)
```

`AstraDock::math`, `AstraDock::dynamics`, `AstraDock::numerics`, `AstraDock::orbit`,
`AstraDock::frames`, and `AstraDock::attitude` are header-only CMake interface targets requiring C++20.
Orbit depends on dynamics, numerics, and math; dynamics, frames, and attitude each depend on
math and numerics. Demo executables link the C++ core and export plain CSV. CTest owns test execution,
while Catch2 provides test cases and assertions.

| Area | Current responsibility | Current non-responsibility |
| --- | --- | --- |
| `cpp/math` | General double-precision `Vector3`, `Matrix3`, `Quaternion`, `EulerAngles` (ZYX), arithmetic, transpose, determinant, Shepperd DCM $\leftrightarrow$ quaternion conversion, angle normalization, and constants | Frames, unit types, dynamics |
| `cpp/dynamics` | Instantaneous point-mass central gravitational acceleration | Time integration, orbit propagation, perturbations |
| `cpp/numerics` | Generic Euler/RK4 steps plus deterministic signed fixed-step propagation and endpoint policy | Physical derivative assembly, adaptive steps, events |
| `cpp/orbit` | Cartesian orbital state, two-body derivative assembly, circular-orbit helpers, invariants, validation diagnostics, `ClassicalOrbitalElements`, perifocal coordinate frame ($PQW$), and bidirectional state $\leftrightarrow$ elements conversions | Perturbations ($J_2$, drag), equinoctial elements |
| `cpp/frames` | Orthonormal `FrameBasis`, LVLH basis construction, Direction Cosine Matrices (`dcm_lvlh_from_eci`, `dcm_eci_from_lvlh`), vector coordinate transformations, and degenerate-state rejection | ECEF, Earth rotation, attitude dynamics |
| `cpp/attitude` | `AttitudeState`, `PrincipalInertia`, `RotationalState`, Euler's rigid-body equations, quaternion kinematics, rotational kinetic energy, body/inertial angular momentum, and RK4 step with unit-norm reprojection | Actuators (reaction wheels, thrusters), sensors, control, perturbations |
| `cpp/spacecraft` | Composite 13-component `SpacecraftState` (Cartesian translational + rotational), `SpacecraftParameters`, `ForceTorqueInput`, unified `spacecraft_state_derivative`, fixed-step RK4 propagation, and geodesic `quaternion_orientation_error_rad` metric | Cross-coupled environmental perturbations, actuators, control |
| `cpp/environment` | Environmental perturbations: Earth oblateness ($J_2$), atmospheric drag with rotating atmosphere, third-body lunar/solar tidal gravity, gravity-gradient body torque, `EnvironmentConfiguration`, `EnvironmentalParameters`, and environmental 6-DOF propagation | Actuators, sensors, active guidance and control |
| `cpp/sensors` | Sensor simulation layer: `DeterministicRng`, `SensorSchedule`, `DropoutWindow`, 6-axis IMU (gyroscope body rate + accelerometer specific force $f = a - g$), GNSS receiver (ECI pos/vel), Star Tracker ($SO(3)$ attitude quaternion perturbation), and line-of-sight relative Range Sensor | State estimation, EKF/UKF, sensor fusion, active guidance |
| `cpp/estimation` | Navigation layer: fixed-size `math::Matrix<R,C>` with Cholesky SPD solves, generic Joseph-form updates with NIS/NEES diagnostics, analytical gravity Jacobian with finite-difference audits, first-order discrete transition $\Phi \approx I + F\Delta t$, coupled process noise $Q$, GNSS measurement model ($H_{\text{GNSS}}$), star tracker measurement model ($H_{\text{ST}}$ with double-cover alignment), relative range measurement model ($H_{\text{range}}$ with singularity guard), 6-state `TranslationalEkf`, 6-state `AttitudeEkf`, IMU-aided specific-force dead reckoning, and the unified 15-state `IntegratedNavigationEkf` ($\delta\mathbf{x} = [\delta\mathbf{r}^T, \delta\mathbf{v}^T, \delta\boldsymbol{\theta}^T, \delta\mathbf{b}_a^T, \delta\mathbf{b}_g^T]^T \in \mathbb{R}^{15}$) with full cross-covariance coupling (attitude/accel-bias to velocity) and first-order covariance reset $J_{\text{reset}}$ | Actuators (M14), guidance, control |
| `tools` | Nominal scenarios, quantitative summaries, deterministic CSV export, validation sweeps, frame demonstrations, elements demos, attitude demos, attitude dynamics demos, 6-DOF integrated demos, environmental perturbation demos, multi-rate sensor demos, EKF demos, IMU-EKF demos, attitude/range EKF demos, and integrated navigation filter demos | General scenario or telemetry framework |
| `python/analysis` | Plotting, display-unit conversion, convergence analysis, frame visualizations, 3D orbital geometry, 3D attitude figures, attitude dynamics invariants, 3D 6-DOF orbit + attitude triad figures, environmental perturbation figures, sensor telemetry, EKF telemetry, IMU-EKF plots, attitude/range EKF plots, and 15-state integrated navigation plots (Plots A–H) from authoritative C++ CSV | Authority over physics or integration |
| `tests/cpp` | Deterministic analytical, invariant, ODE convergence, endpoint, orbit scenario, regression, matrix algebra, coordinate frame, audit property, classical element, quaternion attitude, attitude dynamics, 6-DOF integrated, environmental perturbation, and sensor measurement tests | Higher-fidelity external reference validation |
| `docs` | Architecture, curriculum, concept lessons, and numerical validation reports | Claims of unimplemented behavior |
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

`Vector3` and `Matrix3` intentionally represent generic numeric components. The
physical meaning and coordinate system must be supplied by the containing type,
variable name, or function contract. For example, `position_eci_m` and
`position_lvlh_m` both use `Vector3`, but cannot be added directly without an
explicit frame transformation.

### Coordinate Frame Hierarchy

```text
                     ECI (Earth-Centered Inertial)
                                  │
                                  │ orbital state (r, v)
                                  ▼
                   LVLH (Local Vertical Local Horizontal)
                                  │
                                  │ PLANNED (future attitude milestones)
                                  ▼
                   Spacecraft Body Frame (PLANNED)
```

- **ECI:** Idealized non-rotating Cartesian coordinate frame with origin at Earth's
  center of mass, equatorial fundamental plane, and $+Z_{ECI}$ pointing toward the north
  celestial pole.
- **LVLH:** Spacecraft orbital-local frame where $+X$ points radially outward ($\mathbf{e}_r = \mathbf{r}/\|\mathbf{r}\|$),
  $+Z$ points normal to the orbital plane ($\mathbf{e}_h = \mathbf{h}/\|\mathbf{h}\|$), and $+Y$ completes the
  right-handed triad along-track ($\mathbf{e}_t = \mathbf{e}_h \times \mathbf{e}_r$).
- **Spacecraft Body Frame:** Planned for attitude dynamics and sensor alignment in M07–M09.

### DCM Convention

Direction Cosine Matrices follow the explicit convention $C_{A\_B}$ (or `dcm_a_from_b`):
$$\mathbf{v}_A = C_{A\_B} \cdot \mathbf{v}_B$$
- $C_{LVLH\_ECI}$: Transforms ECI coordinates to LVLH coordinates; its **rows** are the LVLH basis vectors expressed in ECI.
- $C_{ECI\_LVLH} = C_{LVLH\_ECI}^T$: Transforms LVLH coordinates to ECI coordinates; its **columns** are the LVLH basis vectors expressed in ECI.

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
five-orbit Euler/RK4 drift comparison.

The M05 suite adds phase diagnostics (angle, analytical, unwrapped), period
estimation from trajectory crossings, empirical convergence-order measurement,
angular-momentum direction drift, repeatability checks, regression baselines
against M04 values, circular-orbit radius/speed sanity, and timestep-sensitivity
tests for both Euler and RK4.

Plot generation is checked separately from core numerical unit tests.

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
- The M05 validation suite confirms RK4 demonstrates approximately fourth-order
  convergence on the two-body orbital problem, and that the propagator is
  deterministic and repeatable.
