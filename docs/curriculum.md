# AstraDock Curriculum

The repository is a staged aerospace curriculum. Each milestone introduces the
minimum theory needed for its implementation, then requires code and evidence
that make the learning inspectable. Stages remain sequential so that classical
physics and estimation become baselines before higher-level autonomy or machine
learning is attempted.

| Stage | Concepts to learn | Evidence produced |
| --- | --- | --- |
| M00 Repository foundation | Reproducible builds, target-based CMake, CTest, repository conventions, SI/frame discipline. | A documented project that configures, builds, and runs checks from clean instructions. |
| M01 Vectors and math primitives | Physical vectors, magnitude/direction, dot/cross products, normalization, units, frame dependence, floating-point tolerance. | Lesson 001, readable `Vector3`, analytical and invariant unit tests. |
| M02 Two-body gravity | Newtonian gravity, gravitational parameter, Cartesian orbital state, ECI, circular velocity, specific energy and angular momentum. | Lesson 002 and analytical acceleration tests. |
| M03 Numerical integration | ODEs, Euler and RK4, local/global error, timestep, convergence and stability. | Lesson 003, generic fixed-step primitives, and scalar/`Vector3` known-solution convergence results. |
| M04 Orbital propagation | Cartesian state derivatives, initial conditions, sampling, circular orbit geometry, period, conservation diagnostics, and Euler/RK4 drift. | Lesson 004, deterministic 500 km C++ scenarios, analytical/invariant tests, CSV output, and orbit/drift plots. |
| M05 Orbital invariants, accuracy characterization, and propagator validation | Timestep sensitivity, empirical convergence order, phase error, period estimation, energy/angular-momentum drift, repeatability, regression baselines, verification vs validation. | Lesson 005, orbital diagnostics header, timestep-sweep validation demo, 25 C++ validation tests, convergence/phase/plots, and a technical validation report. |
| M06 Coordinate frames | Physical vectors vs components, ECI, LVLH basis construction (radial, along-track, orbit-normal), Direction Cosine Matrices, row/column DCM conventions, norm preservation, moving vs inertial frames. | Lesson 006, `Matrix3`, `FrameBasis`, LVLH APIs, `astradock_frame_demo`, 17 C++ frame/matrix unit tests, 3 frame plots, and an M06 validation report. |

| M07 Classical orbital elements | Keplerian elements ($a, e, i, \Omega, \omega, \nu$), perifocal frame ($PQW$), bidirectional state $\leftrightarrow$ elements conversion, angle normalization, and singularity policies. | Lesson 007, `ClassicalOrbitalElements`, `astradock_elements_demo`, 7 C++ unit tests, 3D orbit geometry/evolution plots, and an M07 validation report. |
| M08 Quaternions & attitude representation | Attitude conventions, scalar-first quaternion, Hamilton product, vector rotation, DCM $\leftrightarrow$ quaternion conversion (Shepperd algorithm), double cover, gimbal lock. | Lesson 008, `Quaternion`, `EulerAngles`, `AttitudeState`, `astradock_attitude_demo`, 9 C++ unit tests, 3 attitude plots, and an M08 validation report. |
| M09 Attitude propagation & kinematics | Quaternion kinematics ($\dot{q} = \frac{1}{2} q \otimes \omega$), angular-rate frames, integration and renormalization. | Constant-rate analytical comparison. |
| M10 Rigid-body dynamics | Torque, inertia, Euler equations, rotational energy and angular momentum, tumbling. | Verified torque-free and forced rotation scenarios. |
| M11 Sensor models | Sampling, white noise, bias, drift, dropout, calibration, seeded randomness. | Lesson on sensors and measurement-statistics tests. |
| M12 Linear Kalman filter | State/covariance prediction, innovation, process and measurement noise, Kalman gain. | Lesson and truth/measurement/estimate plots. |
| M13 Extended Kalman filter | Nonlinear models, Jacobians, linearization, conditioning, covariance consistency. | Lesson and Jacobian/nonlinear tests. |
| M14 Attitude estimation | Quaternion error, gyro bias estimation, multi-rate fusion, attitude error metrics. | Sensor-fusion scenarios and error analysis. |
| M15 Second spacecraft & relative motion | Ownership of target/chaser states, relative initial conditions, Hill/LVLH relative state, Clohessy-Wiltshire equations. | Lesson and analytical relative-motion checks. |
| M16 Rendezvous guidance | Guidance references, hold points, approach profiles, closing-rate and path constraints. | Safe reference trajectories through staged distances. |
| M17 Translational control | Force, mass, acceleration, impulse, saturation and tracking error. | Bounded-thrust tracking tests. |
| M18 Attitude control | PID, LQR, stability, tuning, actuator bounds, settling time. | Detumble and pointing tracking plots. |
| M19 Docking corridor | Approach axis, keep-out zones, corridor geometry, safety margins and abort maneuvers. | Boundary tests and visualized safety geometry. |
| M20 Mission state machine | Modes, guarded transitions, safety interlocks, degraded operation and abort priority. | Nominal and faulted transition tests. |
| M21 Virtual camera | Camera/sensor frames, intrinsics, extrinsics, perspective projection, FOV and clipping. | Synthetic observations with analytical pixel checks. |
| M22 Classical visual pose estimation | Keypoints, correspondences, PnP, pose conventions, reprojection error and outliers. | Known-pose synthetic validation. |
| M23 ML pose estimation | Dataset generation, leakage prevention, pose representations, confidence, generalization and baseline comparison. | Reproducible training/evaluation artifacts. |
| M24 Fault injection | Fault taxonomies, timing, persistence, isolation and deterministic scenario design. | Reproducible sensor/actuator fault cases. |
| M25 Classical fault detection | Residuals, innovations, thresholds, covariance consistency, false alarms and latency. | Detection scorecard against injected faults. |
| M26 ML anomaly detection | Imbalanced evaluation, temporal splits, precision/recall, calibration and classical-baseline comparison. | Honest comparative anomaly study. |
| M27 Monte Carlo validation | Parameter distributions, seeded sampling, confidence intervals, aggregation and mission metrics. | Reproducible robustness report. |
| M28 External dynamics validation | Reference alignment, initial-condition matching, model-fidelity differences and error attribution. | Independent comparison with external tools. |
| M29 Final portfolio documentation | Traceability, reproducibility, technical communication, limitations and evidence-based claims. | A complete portfolio narrative backed by runnable results. |

## Learning workflow for every stage

Before coding, state the physical problem, why the model is needed, coordinate
frames, units, state, inputs, outputs, governing equations, assumptions,
numerical method, known failures, and verification strategy. Then implement a
small change, compile it, run deterministic tests, inspect warnings and
numerical evidence, fix project-caused failures, and update the documentation.

Machine learning stages deliberately follow analytical mechanics, classical
filtering/control, geometric vision, and deterministic fault detection. Their
results must be compared with those baselines rather than assumed superior.

## Current checkpoint

M00-M08 are implemented. The learner understands both orbital translational
state and spacecraft attitude state, the mathematics of unit quaternions
($q = [w, x, y, z]$), vector rotations ($v' = q \otimes v \otimes q^*$),
Hamilton composition ($q_{A\_C} = q_{A\_B} \otimes q_{B\_C}$), Direction
Cosine Matrix conversions (Shepperd algorithm), double-cover properties
($\pm q$), and Euler-angle singularities (gimbal lock).

The immediate boundary is **M09: Attitude Propagation & Quaternion Kinematics**.
Do not implement angular velocity propagation, inertia, or torque ahead of M09.
