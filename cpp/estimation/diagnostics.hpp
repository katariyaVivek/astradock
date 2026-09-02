#pragma once

#include "estimation/translational_ekf.hpp"
#include "math/linalg.hpp"

#include <cmath>
#include <stdexcept>

namespace astradock::estimation {

// Evaluation-only diagnostics. These functions REQUIRE the true state and are
// therefore FORBIDDEN inside the filter itself: a navigation filter that reads
// truth is not navigating. Test harnesses and validation tools call them to
// quantify how well the filter's self-assessed uncertainty matches reality.

// Normalized Estimation Error Squared for the translational state:
//
//     NEES = e^T P^{-1} e,    e = x_true - x_estimate   (6x1)
//
// For a consistent filter and 6 degrees of freedom, NEES follows approximately
// a chi-square distribution with 6 DOF: mean 6, central 95% interval
// [1.64, 14.45]. A filter can track well yet be statistically overconfident
// (NEES >> 6): low error does not imply calibrated covariance.
//
// Computed via Cholesky solve; never forms P^{-1} explicitly.
[[nodiscard]] inline double translational_nees(
    const orbit::CartesianState& truth_state,
    const orbit::CartesianState& estimated_state,
    const TranslationalCovariance& covariance) {
    const math::Vector3 position_error = truth_state.position - estimated_state.position;
    const math::Vector3 velocity_error = truth_state.velocity - estimated_state.velocity;

    const TranslationalMeanState error(std::array<double, 6>{
        position_error.x(),
        position_error.y(),
        position_error.z(),
        velocity_error.x(),
        velocity_error.y(),
        velocity_error.z(),
    });

    if (!error.all_finite() || !covariance.all_finite()) {
        throw std::domain_error("NEES inputs must contain only finite values");
    }

    const TranslationalMeanState whitened = math::solve_spd(covariance, error);
    double nees = 0.0;
    for (std::size_t i = 0; i < k_translational_state_dim; ++i) {
        nees += error(i, 0) * whitened(i, 0);
    }
    if (!std::isfinite(nees) || nees < 0.0) {
        throw std::overflow_error("NEES computation produced an invalid value");
    }
    return nees;
}

// Reference quantiles of the chi-square distribution with 6 degrees of freedom
// (standard tabulated values), used by consistency tests:
inline constexpr double k_chi2_6dof_95_lower = 1.63538;
inline constexpr double k_chi2_6dof_95_upper = 14.4494;
inline constexpr double k_chi2_6dof_mean = 6.0;

// Reference quantiles for 3 degrees of freedom (attitude only, or 3D star tracker innovation):
inline constexpr double k_chi2_3dof_95_lower = 0.215795;
inline constexpr double k_chi2_3dof_95_upper = 9.34840;
inline constexpr double k_chi2_3dof_mean = 3.0;

// Reference quantiles for 1 degree of freedom (scalar range measurement innovation):
inline constexpr double k_chi2_1dof_95_lower = 0.000982;
inline constexpr double k_chi2_1dof_95_upper = 5.02389;
inline constexpr double k_chi2_1dof_mean = 1.0;

// Reference quantiles for 15 degrees of freedom (full integrated navigation state):
inline constexpr double k_chi2_15dof_95_lower = 6.26214;
inline constexpr double k_chi2_15dof_95_upper = 27.4884;
inline constexpr double k_chi2_15dof_mean = 15.0;

// Computes the 3-element physical attitude error vector delta_theta in the spacecraft Body frame:
//
//     q_true = q_nom ⊗ delta_q(delta_theta)
//     delta_q = q_nom^* ⊗ q_true
//     delta_theta ≈ 2 * delta_q.vector_part()
[[nodiscard]] inline math::Vector3 attitude_error_vector_rad(
    const math::Quaternion& truth_orientation,
    const math::Quaternion& estimated_nominal_orientation) {
    if (!math::is_finite(truth_orientation) || !math::is_finite(estimated_nominal_orientation)) {
        throw std::domain_error("Attitude error calculation requires finite quaternions");
    }

    math::Quaternion q_t = truth_orientation;
    const double dot_prod = q_t.w() * estimated_nominal_orientation.w()
                          + q_t.x() * estimated_nominal_orientation.x()
                          + q_t.y() * estimated_nominal_orientation.y()
                          + q_t.z() * estimated_nominal_orientation.z();
    if (dot_prod < 0.0) {
        q_t = math::Quaternion(-q_t.w(), -q_t.x(), -q_t.y(), -q_t.z());
    }

    math::Quaternion q_err = estimated_nominal_orientation.conjugate() * q_t;
    if (q_err.w() < 0.0) {
        q_err = math::Quaternion(-q_err.w(), -q_err.x(), -q_err.y(), -q_err.z());
    }

    return q_err.vector_part() * 2.0;
}

// Normalized Estimation Error Squared for the 3-DOF attitude error:
//
//     NEES_att = delta_theta^T * P_theta^{-1} * delta_theta
//
// Follows approximately a chi-square distribution with 3 degrees of freedom.
[[nodiscard]] inline double attitude_nees(
    const math::Quaternion& truth_orientation,
    const math::Quaternion& estimated_nominal_orientation,
    const math::Matrix<3, 3>& attitude_covariance) {
    const math::Vector3 err = attitude_error_vector_rad(
        truth_orientation, estimated_nominal_orientation);
    const math::ColVector<3> error(std::array<double, 3>{err.x(), err.y(), err.z()});

    if (!error.all_finite() || !attitude_covariance.all_finite()) {
        throw std::domain_error("Attitude NEES inputs must contain only finite values");
    }

    const math::ColVector<3> whitened = math::solve_spd(attitude_covariance, error);
    double nees = 0.0;
    for (std::size_t i = 0; i < 3; ++i) {
        nees += error(i, 0) * whitened(i, 0);
    }
    if (!std::isfinite(nees) || nees < 0.0) {
        throw std::overflow_error("Attitude NEES computation produced an invalid value");
    }
    return nees;
}

// Composite 6-DOF NEES for attitude error and gyro bias error:
//
//     NEES_6 = [delta_theta, delta_b_g]^T * P^{-1} * [delta_theta, delta_b_g]
//
// Follows approximately a chi-square distribution with 6 degrees of freedom.
[[nodiscard]] inline double attitude_bias_nees(
    const math::Quaternion& truth_orientation,
    const math::Vector3& truth_gyro_bias_rad_s,
    const math::Quaternion& estimated_nominal_orientation,
    const math::Vector3& estimated_gyro_bias_rad_s,
    const math::Matrix<6, 6>& covariance) {
    const math::Vector3 att_err = attitude_error_vector_rad(
        truth_orientation, estimated_nominal_orientation);
    const math::Vector3 bias_err = truth_gyro_bias_rad_s - estimated_gyro_bias_rad_s;

    const math::ColVector<6> error(std::array<double, 6>{
        att_err.x(), att_err.y(), att_err.z(),
        bias_err.x(), bias_err.y(), bias_err.z()
    });

    if (!error.all_finite() || !covariance.all_finite()) {
        throw std::domain_error("Composite NEES inputs must contain only finite values");
    }

    const math::ColVector<6> whitened = math::solve_spd(covariance, error);
    double nees = 0.0;
    for (std::size_t i = 0; i < 6; ++i) {
        nees += error(i, 0) * whitened(i, 0);
    }
    if (!std::isfinite(nees) || nees < 0.0) {
        throw std::overflow_error("Composite NEES computation produced an invalid value");
    }
    return nees;
}

// Computes full 15-state NEES:
// delta_x = [delta_r, delta_v, delta_theta, delta_ba, delta_bg]^T in R^15
// NEES_15 = delta_x^T * P^{-1} * delta_x
[[nodiscard]] inline double full_state_nees_15(
    const math::Vector3& truth_pos_eci_m,
    const math::Vector3& truth_vel_eci_mps,
    const math::Quaternion& truth_att_eci_from_body,
    const math::Vector3& truth_accel_bias_body_mps2,
    const math::Vector3& truth_gyro_bias_body_rad_s,
    const math::Vector3& est_pos_eci_m,
    const math::Vector3& est_vel_eci_mps,
    const math::Quaternion& est_att_eci_from_body,
    const math::Vector3& est_accel_bias_body_mps2,
    const math::Vector3& est_gyro_bias_body_rad_s,
    const math::Matrix<15, 15>& covariance) {

    const math::Vector3 pos_err = truth_pos_eci_m - est_pos_eci_m;
    const math::Vector3 vel_err = truth_vel_eci_mps - est_vel_eci_mps;
    const math::Vector3 att_err = attitude_error_vector_rad(
        truth_att_eci_from_body, est_att_eci_from_body);
    const math::Vector3 ba_err = truth_accel_bias_body_mps2 - est_accel_bias_body_mps2;
    const math::Vector3 bg_err = truth_gyro_bias_body_rad_s - est_gyro_bias_body_rad_s;

    math::ColVector<15> error;
    error(0, 0) = pos_err.x(); error(1, 0) = pos_err.y(); error(2, 0) = pos_err.z();
    error(3, 0) = vel_err.x(); error(4, 0) = vel_err.y(); error(5, 0) = vel_err.z();
    error(6, 0) = att_err.x(); error(7, 0) = att_err.y(); error(8, 0) = att_err.z();
    error(9, 0) = ba_err.x();  error(10, 0) = ba_err.y(); error(11, 0) = ba_err.z();
    error(12, 0) = bg_err.x(); error(13, 0) = bg_err.y(); error(14, 0) = bg_err.z();

    if (!error.all_finite() || !covariance.all_finite()) {
        throw std::domain_error("full_state_nees_15 inputs must contain only finite values");
    }

    const math::ColVector<15> whitened = math::solve_spd(covariance, error);
    double nees = 0.0;
    for (std::size_t i = 0; i < 15; ++i) {
        nees += error(i, 0) * whitened(i, 0);
    }
    if (!std::isfinite(nees) || nees < 0.0) {
        throw std::overflow_error("full_state_nees_15 computation produced an invalid value");
    }
    return nees;
}

}  // namespace astradock::estimation
