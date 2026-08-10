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
| M05 Orbital invariants, accuracy characterization, and propagator validation | Timestep sensitivity, measured convergence on orbital dynamics, long-duration error budgets, and independent validation strategy. | A systematic propagator accuracy report with justified operating tolerances. |
| M06 Coordinate frames | ECI, ECEF, LVLH/Hill, bases, handedness, rotations, inverse transforms. | A coordinate-frames lesson with known-axis and round-trip tests. |
| M07 Quaternions | Attitude conventions, multiplication, inverse, vector rotation, rotation matrices, gimbal lock. | A quaternion lesson and known-rotation tests. |
| M08 Attitude propagation | Quaternion kinematics, angular-rate frames, integration and renormalization. | Constant-rate analytical comparison. |
| M09 Rigid-body dynamics | Torque, inertia, Euler equations, rotational energy and angular momentum, tumbling. | Verified torque-free and forced rotation scenarios. |
| M10 Sensor models | Sampling, white noise, bias, drift, dropout, calibration, seeded randomness. | Lesson 006 and measurement-statistics tests. |
| M11 Linear Kalman filter | State/covariance prediction, innovation, process and measurement noise, Kalman gain. | Lesson 007 and truth/measurement/estimate plots. |
| M12 Extended Kalman filter | Nonlinear models, Jacobians, linearization, conditioning, covariance consistency. | Lesson 008 and Jacobian/nonlinear tests. |
| M13 Attitude estimation | Quaternion error, gyro bias estimation, multi-rate fusion, attitude error metrics. | Sensor-fusion scenarios and error analysis. |
| M14 Second spacecraft | Ownership of target/chaser states, relative initial conditions, independent propagation and tumble. | Two-vehicle deterministic scenario. |
| M15 Relative motion / Hill frame | Relative state, Hill basis, Clohessy-Wiltshire equations, linear-model limitations. | Lesson 009 and analytical relative-motion checks. |
| M16 Rendezvous guidance | Guidance references, hold points, approach profiles, closing-rate and path constraints. | Safe reference trajectories through staged distances. |
| M17 Translational control | Force, mass, acceleration, impulse, saturation and tracking error. | Bounded-thrust tracking tests. |
| M18 Attitude control | PID, LQR, stability, tuning, actuator bounds, settling time. | Lessons 010-011 and detumble/tracking plots. |
| M19 Docking corridor | Approach axis, keep-out zones, corridor geometry, safety margins and abort maneuvers. | Boundary tests and visualized safety geometry. |
| M20 Mission state machine | Modes, guarded transitions, safety interlocks, degraded operation and abort priority. | Nominal and faulted transition tests. |
| M21 Virtual camera | Camera/sensor frames, intrinsics, extrinsics, perspective projection, FOV and clipping. | Synthetic observations with analytical pixel checks. |
| M22 Classical visual pose estimation | Keypoints, correspondences, PnP, pose conventions, reprojection error and outliers. | Lesson 012 and known-pose synthetic validation. |
| M23 ML pose estimation | Dataset generation, leakage prevention, pose representations, confidence, generalization and baseline comparison. | Reproducible training/evaluation artifacts. |
| M24 Fault injection | Fault taxonomies, timing, persistence, isolation and deterministic scenario design. | Reproducible sensor/actuator fault cases. |
| M25 Classical fault detection | Residuals, innovations, thresholds, covariance consistency, false alarms and latency. | Detection scorecard against injected faults. |
| M26 ML anomaly detection | Imbalanced evaluation, temporal splits, precision/recall, calibration and classical-baseline comparison. | Honest comparative anomaly study. |
| M27 Monte Carlo validation | Parameter distributions, seeded sampling, confidence intervals, aggregation and mission metrics. | Reproducible robustness report. |
| M28 External dynamics validation | Reference alignment, initial-condition matching, model-fidelity differences and error attribution. | Independent comparison with Basilisk or a justified equivalent. |
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

M00-M04 are implemented. The learner should now be able to connect central
gravity to a first-order position/velocity state derivative, explain how
repeated Euler/RK4 steps create a trajectory, initialize an analytical circular
orbit, and use closure, radius, specific energy, and specific angular momentum
to evaluate numerical fidelity. The 500 km demonstration also makes the
difference between physical-model error and integration error visible and
quantitative.
