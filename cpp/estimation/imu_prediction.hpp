#pragma once

#include "dynamics/two_body.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"

#include <cmath>
#include <stdexcept>

namespace astradock::estimation {

// M13B — IMU-aided translational prediction models.
//
// Physical problem
// ----------------
// M13A predicted the spacecraft state by re-integrating its own dynamics
// model (two-body gravity) between GNSS fixes. That is only legitimate when
// the filter's process model matches reality. A real spacecraft cannot do
// this: it must sense how it is actually accelerating. The onboard inertial
// sensor chain provides exactly that information:
//
//     accelerometer : f_m = f_true + b_a + n_a      (specific force, BODY)
//     gyroscope     : w_m = w_true + b_g + n_g      (body rates, BODY)
//
// where specific force is the NON-GRAVITATIONAL acceleration felt by the
// instrument:
//
//     f = a_inertial - g_gravity
//
// In free fall f is near zero even though gravity is large -- an accelerometer
// in orbit reads approximately ZERO, not one g. To recover the inertial
// acceleration that drives translation the filter must therefore undo the
// measurement physics:
//
//     a_I = C_I_B(q_hat) * (f_m - b_a) + g(r_hat)
//
// This single line is the heart of inertial navigation and carries three
// distinct responsibilities:
//   1. bias removal        (b_a subtracted; KNOWN configuration in M13B),
//   2. frame transformation(BODY -> ECI through the ATTITUDE ESTIMATE q_hat),
//   3. gravity compensation(g(r_hat) added back at the estimated position).
//
// Coordinate frames and units
// ---------------------------
//   f_m, b_a       : spacecraft BODY frame, m/s^2
//   q_hat          : unit quaternion mapping BODY coordinates into ECI
//                    (project convention: q_A_B transforms B-frame vectors
//                    into A-frame; here A = ECI, B = body)
//   r_hat, a, g    : Earth-Centered Inertial frame, m, m/s^2
//   mu             : m^3/s^2
//
// Known-attitude assumption (deliberate staging)
// ----------------------------------------------
// M13B receives the attitude used for the BODY->ECI transformation as an
// EXPLICIT navigation-state input supplied by the caller. It is treated as an
// estimate, never as truth-state access: nothing in this header or in
// TranslationalEkf accepts a SpacecraftState or any other truth container.
// The translational error-state covariance deliberately EXCLUDES attitude
// uncertainty columns (the Jacobian of a_I with respect to attitude error is
// -[C(q_hat)(f_m - b_a)]x, which vanishes identically in free fall where
// f ~= 0). M13C will estimate attitude jointly and fold this coupling into
// the filter; until then the limitation is documented rather than hidden.
//
// Gyroscope handling
// ------------------
// The gyro sample travels through the same pipeline and timestamp discipline
// as the accelerometer, but M13B does NOT use it inside the filter: attitude
// estimation belongs to M13C. Pretending otherwise would duplicate M09/M10
// attitude propagation behind the filter's back.

// Explicit navigation-layer representation of the externally supplied
// attitude estimate. Wrapping the raw quaternion keeps the estimator API
// self-documenting: this is a NAVIGATION STATE consumed by the filter, not a
// truth parameter smuggled in from the simulator.
struct NavigationAttitudeEstimate {
    // Unit quaternion q_ECI_BODY mapping body-frame vector coordinates into
    // ECI coordinates (scalar-first [w, x, y, z], dimensionless).
    math::Quaternion orientation_eci_from_body{};

    // Validates finiteness and approximate unit norm, returning a defensively
    // normalized copy. Tolerance is loose (1e-6) because callers propagate
    // attitudes with their own numerical machinery; gross violations throw.
    [[nodiscard]] math::Quaternion validated_orientation() const {
        if (!math::is_finite(orientation_eci_from_body)) {
            throw std::domain_error(
                "Navigation attitude estimate must contain only finite values");
        }
        const double norm = orientation_eci_from_body.norm();
        if (!std::isfinite(norm) || std::abs(norm - 1.0) > 1.0e-6) {
            throw std::domain_error(
                "Navigation attitude estimate must be a unit quaternion "
                "(|q| - 1 <= 1e-6)");
        }
        return orientation_eci_from_body.normalized();
    }
};

// Known IMU configuration consumed by the translational predictor.
//
// M13B treats both entries as CALIBRATION FACTS, not estimated states:
//   - accelerometer_bias_body_mps2 is SUBTRACTED from every measurement;
//     if it equals the true sensor bias the correction is exact. Section 14
//     of the M13B specification deliberately runs mismatched-bias scenarios
//     to expose the resulting drift, motivating bias states in M13C/D.
//   - accelerometer_noise_std_mps2 (per axis, isotropic white noise, m/s^2)
//     is the physically justified driver of the translational process noise:
//     between updates the dominant stochastic forcing of the estimate IS the
//     accelerometer noise. See imu_prediction_process_noise().
struct ImuPredictionConfig {
    math::Vector3 accelerometer_bias_body_mps2{};
    double accelerometer_noise_std_mps2{0.0};

    void validate() const {
        if (!math::is_finite(accelerometer_bias_body_mps2)) {
            throw std::domain_error("IMU accelerometer bias must be finite");
        }
        if (!std::isfinite(accelerometer_noise_std_mps2)
            || accelerometer_noise_std_mps2 < 0.0) {
            throw std::domain_error(
                "IMU accelerometer noise standard deviation must be finite and >= 0");
        }
    }
};

// One IMU-aided prediction request. All fields are measurements or estimates;
// no truth type can be expressed through this interface.
struct ImuPredictionInput {
    // Interval to integrate across, s. Must be finite and > 0: the multi-rate
    // loop computes dt from actual measurement timestamps, so non-monotone
    // time flow is a caller bug and fails loudly.
    double dt_s{0.0};

    // Measured specific force (bias and noise INCLUDED), BODY frame, m/s^2.
    math::Vector3 measured_specific_force_body_mps2{};

    // Measured angular velocity, BODY frame, rad/s. Carried through the
    // pipeline for timestamping/logging symmetry; NOT consumed by the M13B
    // translational predictor (attitude estimation is deferred to M13C).
    math::Vector3 measured_angular_velocity_body_rad_s{};

    // Externally supplied attitude estimate used for the BODY -> ECI
    // transformation across this interval (zero-order hold).
    NavigationAttitudeEstimate attitude{};
};

// Diagnostics returned by each IMU-aided prediction step.
struct ImuPredictionDiagnostics {
    double dt_s{0.0};
    // Reconstructed inertial acceleration actually integrated:
    //   a = C(q_hat)(f_m - b_a) + g(r_hat_before_step)
    math::Vector3 inertial_acceleration_eci_mps2{};
    // Post-prediction mean state (prior estimate at the next epoch).
    math::Vector3 predicted_position_eci_m{};
    math::Vector3 predicted_velocity_eci_mps{};
};

// Reconstructs the inertial acceleration from one IMU observation.
//
//     a_I = C_I_B(q_hat) * (f_m - b_a) + g(r_hat)
//
// This is the production implementation of the central M13B equation chain:
// specific force (BODY) -> attitude rotation -> gravity compensation ->
// inertial acceleration (ECI). It is a pure function of estimates,
// measurements, and configuration, so tests can audit it against independent
// hand computations without any truth involvement.
//
// Note on the gravity term: g is evaluated at the CURRENT estimate position
// (before integration). Across finite dt the position changes, which is why
// the predictor integrates the full nonlinear ODE with RK4 rather than
// freezing this instantaneous value.
[[nodiscard]] inline math::Vector3 imu_inertial_acceleration_eci(
    const math::Vector3& measured_specific_force_body_mps2,
    const math::Vector3& configured_accel_bias_body_mps2,
    const math::Quaternion& attitude_eci_from_body_estimate,
    const math::Vector3& position_eci_m,
    double gravitational_parameter_m3_per_s2) {
    if (!math::is_finite(measured_specific_force_body_mps2)
        || !math::is_finite(configured_accel_bias_body_mps2)
        || !math::is_finite(position_eci_m)) {
        throw std::domain_error("IMU acceleration inputs must contain only finite values");
    }

    // 1. Bias removal: the sensor ADDED bias at measurement time; the filter
    //    removes its best (M13B: configured) knowledge of it.
    const math::Vector3 corrected_specific_force_body =
        measured_specific_force_body_mps2 - configured_accel_bias_body_mps2;

    // 2. Frame transformation BODY -> ECI through the supplied ESTIMATE.
    //    rotate_vector implements v_ECI = q * v_BODY * q^*, i.e. C_I_B v_B.
    const math::Vector3 specific_force_eci =
        attitude_eci_from_body_estimate.rotate_vector(corrected_specific_force_body);

    // 3. Gravity compensation: specific force EXCLUDED gravity, so it must be
    //    added back to obtain the inertial acceleration.
    const math::Vector3 gravity_eci = dynamics::two_body_acceleration(
        position_eci_m, gravitational_parameter_m3_per_s2);

    return specific_force_eci + gravity_eci;
}

// The discrete translational process noise driven by this configuration
// (imu_prediction_process_noise) lives in translational_ekf.hpp beside the
// generic white-acceleration discretization it reuses, keeping this header
// free of filter-state types.

}  // namespace astradock::estimation
