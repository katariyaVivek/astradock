#pragma once

#include "estimation/kalman.hpp"
#include "estimation/imu_prediction.hpp"
#include "math/constants.hpp"
#include "math/linalg.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"
#include "sensors/star_tracker.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>

namespace astradock::estimation {

// M13C — Attitude Error-State Extended Kalman Filter (MEKF).
//
// Physical problem
// ----------------
// Spacecraft orientation must be estimated from rate-gyro and star-tracker
// measurements. While an orientation in 3D space has 3 physical degrees of
// freedom, a unit quaternion has 4 numerical components constrained by ||q|| = 1.
// Directly running an additive EKF on a 4-component quaternion:
//   1. Violates the unit-norm constraint during additive corrections (x + K y).
//   2. Suffers from a singular 4x4 covariance matrix (rank 3 on the tangent space).
//   3. Fails to account for the double-cover property (q and -q represent identical
//      physical rotations, yet their component difference is 2*q).
//
// The Multiplicative Extended Kalman Filter (MEKF) resolves this by splitting
// the attitude state into:
//   - A 4-component NOMINAL unit quaternion q_nom representing the large rotation.
//   - A 3-component LOCAL error-vector delta_theta parameterizing small rotations
//     away from the nominal attitude in the spacecraft BODY frame.
//
// State representation
// --------------------
// Nominal state:
//     q_nom       orientation of spacecraft Body frame relative to ECI (q_ECI_Body)
//     b_g         gyroscope bias in spacecraft Body frame (rad/s)
//
// Error state (6x1, Body frame, SI units):
//     delta_x = [ delta_theta_x, delta_theta_y, delta_theta_z,
//                 delta_b_gx,    delta_b_gy,    delta_b_gz    ]^T
//
// Convention:
//     q_true = q_nom ⊗ delta_q(delta_theta)
//     delta_q(delta_theta) ≈ [ 1, 0.5 * delta_theta ]^T
//
// Covariance P is 6x6 with indices:
//     0..2: attitude error covariance (rad^2)
//     3..5: gyro bias error covariance ((rad/s)^2)

inline constexpr std::size_t k_attitude_state_dim = 6;
inline constexpr std::size_t k_attitude_meas_dim = 3;

using AttitudeCovariance = math::Matrix<k_attitude_state_dim, k_attitude_state_dim>;
using AttitudeErrorState = math::ColVector<k_attitude_state_dim>;

// Estimated attitude state consisting of nominal quaternion and estimated gyro bias.
struct AttitudeEstimate {
    math::Quaternion nominal_orientation{math::Quaternion::identity()};
    math::Vector3 gyro_bias_rad_s{0.0, 0.0, 0.0};

    [[nodiscard]] NavigationAttitudeEstimate to_navigation_attitude() const noexcept {
        return NavigationAttitudeEstimate{nominal_orientation};
    }
};

// Configuration parameters for the attitude error-state EKF.
struct AttitudeEkfConfig {
    // Gyro rate noise standard deviation (rad/s)
    double gyro_noise_std_rad_s{1.0e-3};

    // Gyro bias random-walk diffusion standard deviation (rad/s^2)
    // 0.0 models strictly constant bias during propagation
    double gyro_bias_walk_std_rad_s2{0.0};

    // Star tracker 1-sigma isotropic attitude noise (radians)
    double star_tracker_noise_std_rad{1.0e-3};

    // Initial 1-sigma uncertainties
    double initial_attitude_error_std_rad{0.05};       // ~2.9 deg
    double initial_gyro_bias_std_rad_s{1.0e-3};        // 1 mrad/s

    constexpr AttitudeEkfConfig() noexcept = default;

    [[nodiscard]] AttitudeCovariance initial_covariance() const {
        AttitudeCovariance p;
        const double var_att = initial_attitude_error_std_rad * initial_attitude_error_std_rad;
        const double var_bias = initial_gyro_bias_std_rad_s * initial_gyro_bias_std_rad_s;
        p(0, 0) = var_att;
        p(1, 1) = var_att;
        p(2, 2) = var_att;
        p(3, 3) = var_bias;
        p(4, 4) = var_bias;
        p(5, 5) = var_bias;
        return p;
    }
};

// Diagnostics from a star tracker measurement update.
struct StarTrackerUpdateDiagnostics {
    math::Vector3 innovation_rad{};
    math::Matrix<3, 3> innovation_covariance{};
    double normalized_innovation_squared{0.0};
    math::Vector3 attitude_correction_rad{};
    math::Vector3 bias_correction_rad_s{};
};

// Skew-symmetric cross-product matrix: [v]x * u = v x u.
[[nodiscard]] inline math::Matrix<3, 3> skew_symmetric(const math::Vector3& v) noexcept {
    math::Matrix<3, 3> m;
    m(0, 1) = -v.z();
    m(0, 2) =  v.y();
    m(1, 0) =  v.z();
    m(1, 2) = -v.x();
    m(2, 0) = -v.y();
    m(2, 1) =  v.x();
    return m;
}

// Continuous-time attitude error dynamics Jacobian F (6x6):
//
//   delta_dot_theta = - [omega_hat]x * delta_theta - delta_b_g - n_g
//   delta_dot_b_g   = n_bg
//
//         [ -[omega_hat]x   -I_3 ]
//     F = [                      ]
//         [    0_3x3        0_3x3]
//
// Governing derivation:
//   q_true = q_nom ⊗ delta_q(delta_theta)
//   dot(q_true) = 0.5 * q_true ⊗ [0, omega]
//   dot(q_nom)  = 0.5 * q_nom ⊗ [0, omega_hat]
//   Differentiating and expanding to first order in delta_theta and delta_b_g
//   yields delta_dot_theta = -omega_hat x delta_theta - delta_b_g.
[[nodiscard]] inline AttitudeCovariance attitude_error_dynamics_matrix(
    const math::Vector3& estimated_angular_velocity_rad_s) {
    if (!math::is_finite(estimated_angular_velocity_rad_s)) {
        throw std::domain_error("Angular velocity must contain only finite values for dynamics matrix");
    }

    AttitudeCovariance f;
    const math::Matrix<3, 3> neg_cross = skew_symmetric(estimated_angular_velocity_rad_s) * (-1.0);
    for (std::size_t r = 0; r < 3; ++r) {
        for (std::size_t c = 0; c < 3; ++c) {
            f(r, c) = neg_cross(r, c);
        }
        f(r, r + 3) = -1.0;
    }
    return f;
}

// Discrete-time state transition matrix Phi (6x6) to first order in dt:
//
//     Phi ≈ I_6 + F * dt
//
//         [ I_3 - [omega_hat]x * dt   -I_3 * dt ]
//     Phi = [                                     ]
//         [          0_3x3               I_3    ]
[[nodiscard]] inline AttitudeCovariance attitude_discrete_state_transition(
    const math::Vector3& estimated_angular_velocity_rad_s,
    double dt_s) {
    if (!std::isfinite(dt_s)) {
        throw std::domain_error("Timestep must be finite for discrete state transition");
    }
    const AttitudeCovariance f = attitude_error_dynamics_matrix(estimated_angular_velocity_rad_s);
    return AttitudeCovariance::identity() + f * dt_s;
}

// Discrete process noise covariance Q (6x6):
//
//     Q_theta = sigma_g^2 * dt * I_3
//     Q_bias  = sigma_bg^2 * dt * I_3
[[nodiscard]] inline AttitudeCovariance attitude_discrete_process_noise(
    const AttitudeEkfConfig& config,
    double dt_s) {
    if (!std::isfinite(dt_s) || dt_s < 0.0) {
        throw std::domain_error("Timestep must be non-negative and finite for process noise");
    }
    AttitudeCovariance q;
    const double q_theta = config.gyro_noise_std_rad_s * config.gyro_noise_std_rad_s * dt_s;
    const double q_bias = config.gyro_bias_walk_std_rad_s2 * config.gyro_bias_walk_std_rad_s2 * dt_s;

    q(0, 0) = q_theta;
    q(1, 1) = q_theta;
    q(2, 2) = q_theta;
    q(3, 3) = q_bias;
    q(4, 4) = q_bias;
    q(5, 5) = q_bias;
    return q;
}

// Advances the nominal quaternion under constant estimated angular velocity over dt:
//
//     theta = ||omega_hat|| * dt
//     delta_q = [ cos(theta/2), (omega_hat / ||omega_hat||) * sin(theta/2) ]
//     q_nom_new = (q_nom ⊗ delta_q).normalized()
[[nodiscard]] inline math::Quaternion propagate_nominal_quaternion(
    const math::Quaternion& nominal_q,
    const math::Vector3& estimated_angular_velocity_rad_s,
    double dt_s) {
    if (!math::is_finite(nominal_q) || !math::is_finite(estimated_angular_velocity_rad_s)
        || !std::isfinite(dt_s)) {
        throw std::domain_error("Inputs must be finite for quaternion propagation");
    }
    const double rate = estimated_angular_velocity_rad_s.norm();
    const double angle = rate * dt_s;

    math::Quaternion delta_q;
    if (angle < 1.0e-8) {
        // High-order series expansion for small angles to prevent division by near-zero rate
        const double half_angle = 0.5 * angle;
        const double half_dt = 0.5 * dt_s;
        const double series_s = half_dt * (1.0 - (angle * angle) / 24.0);
        delta_q = math::Quaternion(
            1.0 - 0.5 * half_angle * half_angle,
            estimated_angular_velocity_rad_s.x() * series_s,
            estimated_angular_velocity_rad_s.y() * series_s,
            estimated_angular_velocity_rad_s.z() * series_s
        );
    } else {
        const double half_angle = 0.5 * angle;
        const double sin_coeff = std::sin(half_angle) / rate;
        delta_q = math::Quaternion(
            std::cos(half_angle),
            estimated_angular_velocity_rad_s.x() * sin_coeff,
            estimated_angular_velocity_rad_s.y() * sin_coeff,
            estimated_angular_velocity_rad_s.z() * sin_coeff
        );
    }

    return (nominal_q * delta_q).normalized();
}

// Evaluates the 3-element physical attitude residual from a star tracker measurement:
//
// Convention:
//   q_meas ≈ q_nom ⊗ delta_q(delta_theta)
//   => delta_q = q_nom^* ⊗ q_meas
//
// Double-cover resolution:
//   In SO(3), q and -q represent the identical physical rotation. If q_meas has a negative
//   inner product with q_nom, we flip the sign of q_meas before forming delta_q.
//   Then for small errors:
//   delta_theta_meas = 2 * delta_q.vector_part().
[[nodiscard]] inline math::Vector3 star_tracker_measurement_residual(
    const math::Quaternion& nominal_q,
    const math::Quaternion& measured_q) {
    if (!math::is_finite(nominal_q) || !math::is_finite(measured_q)) {
        throw std::domain_error("Quaternions must be finite for residual calculation");
    }

    math::Quaternion q_m = measured_q;
    // Sign-align measured quaternion with nominal quaternion to preserve shortest rotation arc
    const double dot_prod = q_m.w() * nominal_q.w() + q_m.x() * nominal_q.x()
                          + q_m.y() * nominal_q.y() + q_m.z() * nominal_q.z();
    if (dot_prod < 0.0) {
        q_m = math::Quaternion(-q_m.w(), -q_m.x(), -q_m.y(), -q_m.z());
    }

    math::Quaternion q_err = nominal_q.conjugate() * q_m;
    // Ensure shortest representation along error quaternion
    if (q_err.w() < 0.0) {
        q_err = math::Quaternion(-q_err.w(), -q_err.x(), -q_err.y(), -q_err.z());
    }

    return q_err.vector_part() * 2.0;
}

// Measurement matrix H (3x6) mapping the error state delta_x to the attitude innovation:
//
//     y = delta_theta_meas = H * delta_x + v
//     H = [ I_3,  0_3x3 ]
[[nodiscard]] inline math::Matrix<3, 6> star_tracker_measurement_jacobian() noexcept {
    math::Matrix<3, 6> h;
    h(0, 0) = 1.0;
    h(1, 1) = 1.0;
    h(2, 2) = 1.0;
    return h;
}

// First-order covariance reset transformation following attitude error injection:
//
// When the estimated error delta_theta is injected into the nominal quaternion:
//     q_nom^+ = q_nom^- ⊗ delta_q(delta_theta)
// the local coordinate frame changes. The new post-injection error state delta_theta^+
// relates to the pre-reset error via the Jacobian:
//     G_reset = I_3 - 0.5 * [delta_theta]x
//
// For the 6-element error state [delta_theta, delta_b_g]^T:
//     J_reset = [ I_3 - 0.5 * [delta_theta]x    0_3x3 ]
//               [           0_3x3                I_3  ]
//
// Covariance reset:
//     P_reset = J_reset * P * J_reset^T
[[nodiscard]] inline AttitudeCovariance attitude_covariance_reset(
    const AttitudeCovariance& covariance,
    const math::Vector3& delta_theta) {
    if (!covariance.all_finite() || !math::is_finite(delta_theta)) {
        throw std::domain_error("Covariance reset inputs must be finite");
    }

    AttitudeCovariance j = AttitudeCovariance::identity();
    const math::Matrix<3, 3> skew = skew_symmetric(delta_theta) * 0.5;
    for (std::size_t r = 0; r < 3; ++r) {
        for (std::size_t c = 0; c < 3; ++c) {
            j(r, c) -= skew(r, c);
        }
    }

    return (j * covariance * j.transpose()).symmetrized();
}

// Multiplicative Extended Kalman Filter for attitude and gyro bias.
class AttitudeEkf {
public:
    explicit AttitudeEkf(
        const AttitudeEkfConfig& config,
        const AttitudeEstimate& initial_estimate = AttitudeEstimate{})
        : config_(config),
          estimate_(initial_estimate),
          covariance_(config.initial_covariance()) {
        validate_config(config_);
        if (!math::is_finite(estimate_.nominal_orientation)
            || !math::is_finite(estimate_.gyro_bias_rad_s)) {
            throw std::domain_error("Initial attitude estimate must be finite");
        }
        estimate_.nominal_orientation = estimate_.nominal_orientation.normalized();
    }

    [[nodiscard]] const AttitudeEkfConfig& config() const noexcept {
        return config_;
    }

    [[nodiscard]] const AttitudeEstimate& estimate() const noexcept {
        return estimate_;
    }

    [[nodiscard]] const AttitudeCovariance& covariance() const noexcept {
        return covariance_;
    }

    [[nodiscard]] const StarTrackerUpdateDiagnostics& last_star_tracker_diagnostics() const noexcept {
        return last_st_diagnostics_;
    }

    // Resets filter state and covariance to initial conditions.
    void reset(
        const AttitudeEstimate& estimate = AttitudeEstimate{},
        const std::optional<AttitudeCovariance>& covariance = std::nullopt) {
        estimate_ = estimate;
        estimate_.nominal_orientation = estimate_.nominal_orientation.normalized();
        if (!covariance.has_value()) {
            covariance_ = config_.initial_covariance();
        } else {
            covariance_ = *covariance;
        }
        last_st_diagnostics_ = StarTrackerUpdateDiagnostics{};
    }

    // Time update (prediction) driven by body rate-gyro measurement:
    //
    // 1. Unbias measured body rate: omega_hat = omega_meas - b_g_hat.
    // 2. Propagate nominal quaternion: q_nom(t + dt) = q_nom(t) ⊗ delta_q(omega_hat * dt).
    // 3. Propagate error covariance: P = Phi * P * Phi^T + Q.
    void predict(const math::Vector3& gyro_meas_rad_s, double dt_s) {
        if (!math::is_finite(gyro_meas_rad_s)) {
            throw std::domain_error("Gyro measurement must be finite for EKF prediction");
        }
        if (!std::isfinite(dt_s) || dt_s <= 0.0) {
            throw std::domain_error("Prediction timestep must be positive and finite");
        }

        const math::Vector3 omega_hat = gyro_meas_rad_s - estimate_.gyro_bias_rad_s;

        // Nominal state propagation
        estimate_.nominal_orientation = propagate_nominal_quaternion(
            estimate_.nominal_orientation, omega_hat, dt_s);

        // Linearized error covariance propagation
        const AttitudeCovariance phi = attitude_discrete_state_transition(omega_hat, dt_s);
        const AttitudeCovariance q = attitude_discrete_process_noise(config_, dt_s);

        covariance_ = (phi * covariance_ * phi.transpose() + q).symmetrized();
        validate_covariance(covariance_, "Attitude EKF predict: covariance");
    }

    // Measurement update from an optical star tracker:
    //
    // 1. Extract physical 3-element attitude innovation y = 2 * delta_q(q_nom, q_meas).vector.
    // 2. Compute innovation covariance S = H P H^T + R = P_theta + R.
    // 3. Compute Kalman gain K = P H^T S^(-1).
    // 4. Update error state: [delta_theta, delta_b_g]^T = K y.
    // 5. Compute Joseph-form updated covariance P^+.
    // 6. Inject attitude correction: q_nom = (q_nom ⊗ delta_q(delta_theta)).normalized().
    // 7. Inject bias correction: b_g_hat += delta_b_g.
    // 8. Apply first-order covariance reset: P = J_reset * P^+ * J_reset^T.
    // 9. Reset error state: delta_x = 0.
    bool update_star_tracker(const sensors::StarTrackerMeasurement& meas) {
        if (!meas.valid) {
            return false;
        }
        if (!math::is_finite(meas.orientation_eci_from_body)) {
            throw std::domain_error("Star tracker orientation must be finite");
        }

        validate_covariance(covariance_, "Attitude EKF update: prior covariance");

        // 1. Physical 3-element attitude innovation
        const math::Vector3 innovation_vec = star_tracker_measurement_residual(
            estimate_.nominal_orientation, meas.orientation_eci_from_body);
        const math::ColVector<3> innovation(std::array<double, 3>{
            innovation_vec.x(), innovation_vec.y(), innovation_vec.z()
        });

        // 2. Innovation covariance S = H P H^T + R
        // If star_tracker_noise_std_rad is the total 1-sigma isotropic rotation angle sigma,
        // then the per-axis variance is sigma^2 / 3 because E[u_i^2] = 1/3 for uniform unit vectors on S^2.
        const math::Matrix<3, 6> h = star_tracker_measurement_jacobian();
        math::Matrix<3, 3> r_mat;
        const double var_st = (config_.star_tracker_noise_std_rad * config_.star_tracker_noise_std_rad) / 3.0;
        r_mat(0, 0) = var_st;
        r_mat(1, 1) = var_st;
        r_mat(2, 2) = var_st;

        const math::Matrix<3, 3> s_mat = (h * covariance_ * h.transpose() + r_mat).symmetrized();

        // 3. Kalman gain K = P H^T S^(-1)
        const math::Matrix<6, 3> kalman_gain =
            math::solve_spd_multi(s_mat, h * covariance_).transpose();

        // 4. Normalized innovation squared (NIS) = y^T S^(-1) y
        const math::ColVector<3> whitened = math::solve_spd(s_mat, innovation);
        double nis = 0.0;
        for (std::size_t i = 0; i < 3; ++i) {
            nis += innovation(i, 0) * whitened(i, 0);
        }

        // 5. Error state correction
        const AttitudeErrorState delta_x = kalman_gain * innovation;
        const math::Vector3 delta_theta{delta_x(0, 0), delta_x(1, 0), delta_x(2, 0)};
        const math::Vector3 delta_bias{delta_x(3, 0), delta_x(4, 0), delta_x(5, 0)};

        // 6. Joseph-form covariance update
        const AttitudeCovariance eye_minus_kh = AttitudeCovariance::identity() - kalman_gain * h;
        const AttitudeCovariance p_post =
            (eye_minus_kh * covariance_ * eye_minus_kh.transpose()
             + kalman_gain * r_mat * kalman_gain.transpose()).symmetrized();

        // 7. Inject attitude correction
        const double theta_corr_norm = delta_theta.norm();
        math::Quaternion q_corr;
        if (theta_corr_norm < 1.0e-8) {
            q_corr = math::Quaternion(1.0, 0.5 * delta_theta.x(), 0.5 * delta_theta.y(), 0.5 * delta_theta.z());
        } else {
            q_corr = math::Quaternion::from_axis_angle(delta_theta / theta_corr_norm, theta_corr_norm);
        }
        estimate_.nominal_orientation = (estimate_.nominal_orientation * q_corr).normalized();

        // 8. Inject gyro bias correction
        estimate_.gyro_bias_rad_s = estimate_.gyro_bias_rad_s + delta_bias;

        // 9. Apply covariance reset
        covariance_ = attitude_covariance_reset(p_post, delta_theta);
        validate_covariance(covariance_, "Attitude EKF update: reset covariance");

        // Record diagnostics
        last_st_diagnostics_ = StarTrackerUpdateDiagnostics{
            innovation_vec,
            s_mat,
            nis,
            delta_theta,
            delta_bias
        };

        return true;
    }

private:
    static void validate_config(const AttitudeEkfConfig& config) {
        if (!std::isfinite(config.gyro_noise_std_rad_s) || config.gyro_noise_std_rad_s <= 0.0) {
            throw std::domain_error("Attitude EKF gyro noise standard deviation must be positive and finite");
        }
        if (!std::isfinite(config.gyro_bias_walk_std_rad_s2) || config.gyro_bias_walk_std_rad_s2 < 0.0) {
            throw std::domain_error("Attitude EKF gyro bias walk standard deviation must be non-negative and finite");
        }
        if (!std::isfinite(config.star_tracker_noise_std_rad) || config.star_tracker_noise_std_rad <= 0.0) {
            throw std::domain_error("Attitude EKF star tracker noise standard deviation must be positive and finite");
        }
    }

    AttitudeEkfConfig config_;
    AttitudeEstimate estimate_;
    AttitudeCovariance covariance_;
    StarTrackerUpdateDiagnostics last_st_diagnostics_{};
};

}  // namespace astradock::estimation
