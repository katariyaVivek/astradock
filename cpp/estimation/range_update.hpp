#pragma once

#include "estimation/kalman.hpp"
#include "math/linalg.hpp"
#include "math/vector3.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace astradock::estimation {

// M13C — Relative Range Measurement Update.
//
// Physical problem
// ----------------
// During orbital proximity operations and rendezvous, a line-of-sight range sensor
// (e.g., laser rangefinder, radar, or optical time-of-flight) provides direct scalar
// distance measurements between the chaser spacecraft and a target:
//
//     rho = || r_target - r_spacecraft ||   (m)
//
// Unlike GNSS position fixes (which provide direct linear observations of absolute ECI
// coordinates, H = [I, 0]), range is fundamentally NONLINEAR in the Cartesian position
// coordinates:
//
//     rho(r_s) = sqrt( (r_tx - r_sx)^2 + (r_ty - r_sy)^2 + (r_tz - r_sz)^2 )
//
// Measurement Jacobian
// --------------------
// Differentiating rho with respect to spacecraft position r_s:
//
//     d(rho) / d(r_s) = - (r_target - r_spacecraft)^T / || r_target - r_spacecraft ||
//                     = - u_los^T
//
// where u_los = (r_target - r_s) / ||r_target - r_s|| is the unit line-of-sight vector
// pointing FROM the chaser spacecraft TO the target.
//
// Properties:
//   1. Directional sensitivity: Range observes position errors ONLY along the line-of-sight.
//      Cross-track position errors are orthogonal to u_los and remain unobservable from
//      a single range measurement (their Jacobian projection is zero to first order).
//   2. Translation Invariance: Shifting both spacecraft and target by any vector c leaves
//      rho and H_rho strictly unchanged: ||(r_t + c) - (r_s + c)|| = ||r_t - r_s||.
//   3. Geometric Singularity: At r_s == r_target (rho = 0), the line-of-sight direction
//      is undefined (0 / 0). The estimator explicitly detects and rejects this condition.

inline constexpr double k_min_range_separation_m = 1.0e-6;

// Diagnostics from a scalar range measurement update.
struct RangeUpdateDiagnostics {
    double predicted_range_m{0.0};
    double measured_range_m{0.0};
    double innovation_m{0.0};
    double innovation_variance_m2{0.0};
    double normalized_innovation_squared{0.0};
    math::ColVector<6> kalman_gain{};
};

// Computes predicted scalar range from spacecraft position to target position:
//
//     rho = || r_target - r_spacecraft ||
//
// Throws std::domain_error if inputs are non-finite or if bodies coincide (rho < 1e-6 m).
[[nodiscard]] inline double predicted_range(
    const math::Vector3& spacecraft_position_eci_m,
    const math::Vector3& target_position_eci_m) {
    if (!math::is_finite(spacecraft_position_eci_m) || !math::is_finite(target_position_eci_m)) {
        throw std::domain_error("Positions must contain only finite values for range prediction");
    }

    const math::Vector3 delta_r = target_position_eci_m - spacecraft_position_eci_m;
    const double rho = delta_r.norm();
    if (rho < k_min_range_separation_m) {
        throw std::domain_error("Range measurement geometry is singular: target and spacecraft positions coincide");
    }
    return rho;
}

// Evaluates the 1x6 measurement Jacobian H_rho of scalar range with respect to the
// 6-component translational state [r_s, v_s]^T:
//
//     H_rho = [ - (r_target - r_s)^T / ||r_target - r_s||,  0_1x3 ]
[[nodiscard]] inline math::Matrix<1, 6> range_measurement_jacobian(
    const math::Vector3& spacecraft_position_eci_m,
    const math::Vector3& target_position_eci_m) {
    if (!math::is_finite(spacecraft_position_eci_m) || !math::is_finite(target_position_eci_m)) {
        throw std::domain_error("Positions must contain only finite values for range Jacobian");
    }

    const math::Vector3 delta_r = target_position_eci_m - spacecraft_position_eci_m;
    const double rho = delta_r.norm();
    if (rho < k_min_range_separation_m) {
        throw std::domain_error("Range measurement Jacobian is undefined at zero separation");
    }

    const double inv_rho = 1.0 / rho;
    // Note negative sign: d(rho)/d(r_s) = - (r_t - r_s) / ||r_t - r_s||
    math::Matrix<1, 6> h;
    h(0, 0) = -delta_r.x() * inv_rho;
    h(0, 1) = -delta_r.y() * inv_rho;
    h(0, 2) = -delta_r.z() * inv_rho;
    h(0, 3) = 0.0;
    h(0, 4) = 0.0;
    h(0, 5) = 0.0;
    return h;
}

// Performs a scalar range measurement update on a 6-component translational state and covariance:
//
// Inputs:
//   prior_state: 6x1 [r_s, v_s]^T in ECI (m, m/s)
//   prior_covariance: 6x6 covariance matrix P
//   measured_range_m: scalar range observation from sensor (m)
//   target_position_eci_m: known target position in ECI (m)
//   range_noise_std_m: 1-sigma sensor noise (m)
//
// Returns EkfUpdateResult<6, 1> and updates out_diagnostics.
[[nodiscard]] inline EkfUpdateResult<6, 1> execute_range_update(
    const math::ColVector<6>& prior_state,
    const math::Matrix<6, 6>& prior_covariance,
    double measured_range_m,
    const math::Vector3& target_position_eci_m,
    double range_noise_std_m,
    RangeUpdateDiagnostics* out_diagnostics = nullptr) {
    if (!prior_state.all_finite() || !std::isfinite(measured_range_m)
        || !math::is_finite(target_position_eci_m)) {
        throw std::domain_error("Range update inputs must contain only finite values");
    }
    if (!std::isfinite(range_noise_std_m) || range_noise_std_m <= 0.0) {
        throw std::domain_error("Range noise standard deviation must be positive and finite");
    }
    validate_covariance(prior_covariance, "Range update: prior covariance");

    const math::Vector3 sc_pos{prior_state(0, 0), prior_state(1, 0), prior_state(2, 0)};
    const double pred_rho = predicted_range(sc_pos, target_position_eci_m);
    const math::Matrix<1, 6> h = range_measurement_jacobian(sc_pos, target_position_eci_m);

    math::Matrix<1, 1> z_meas;
    z_meas(0, 0) = measured_range_m;

    math::Matrix<1, 1> z_pred;
    z_pred(0, 0) = pred_rho;

    math::Matrix<1, 1> r_mat;
    r_mat(0, 0) = range_noise_std_m * range_noise_std_m;

    const EkfUpdateResult<6, 1> result = ekf_update(
        prior_state,
        prior_covariance,
        z_meas,
        z_pred,
        h,
        r_mat
    );

    if (out_diagnostics != nullptr) {
        out_diagnostics->predicted_range_m = pred_rho;
        out_diagnostics->measured_range_m = measured_range_m;
        out_diagnostics->innovation_m = result.innovation(0, 0);
        out_diagnostics->innovation_variance_m2 = result.innovation_covariance(0, 0);
        out_diagnostics->normalized_innovation_squared = result.normalized_innovation_squared;
        out_diagnostics->kalman_gain = result.kalman_gain;
    }

    return result;
}

}  // namespace astradock::estimation
