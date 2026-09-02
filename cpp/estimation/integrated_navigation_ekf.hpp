#pragma once

#include "estimation/kalman.hpp"
#include "estimation/translational_ekf.hpp"
#include "math/constants.hpp"
#include "math/linalg.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"
#include "numerics/fixed_step_propagation.hpp"
#include "orbit/cartesian_state.hpp"
#include "orbit/two_body_orbit.hpp"
#include "sensors/gnss.hpp"
#include "sensors/star_tracker.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>

namespace astradock::estimation {

// M13D — Integrated Multi-Rate 15-State Navigation Filter.
//
// Physical Problem:
// In proximity operations and orbital rendezvous, translational motion and attitude
// are inextricably coupled through specific-force measurements:
//
//     a_I = C_I_B(q)(f_m - b_a) + g(r_I)
//
// A small attitude error delta_theta rotates the large specific-force vector f_m
// into the wrong inertial direction, generating a spurious acceleration:
//
//     delta_a_I approx -C_I_B(q) [f_m - b_a]x delta_theta
//
// Similarly, an uncalibrated accelerometer bias delta_b_a produces a continuous
// translational acceleration error:
//
//     delta_a_I approx -C_I_B(q) delta_b_a
//
// Separating translation and attitude into independent filters ignores this cross-coupling,
// producing inconsistent covariance estimates, understating position uncertainty during
// thrust or non-gravitational maneuvers, and preventing translational observations (e.g. GNSS)
// from calibrating sensor biases or attitude errors.
//
// M13D unifies translation, attitude, and sensor biases into one 15-state estimation error vector:
//
//     delta_x = [ delta_r^T, delta_v^T, delta_theta^T, delta_b_a^T, delta_b_g^T ]^T in R^15
//
// State Ordering (strictly 0-indexed):
//   delta_r   : indices 0..2  (ECI position error, m)
//   delta_v   : indices 3..5  (ECI velocity error, m/s)
//   delta_theta: indices 6..8  (Spacecraft BODY attitude error vector, rad)
//   delta_b_a : indices 9..11 (Spacecraft BODY accelerometer bias error, m/s^2)
//   delta_b_g : indices 12..14(Spacecraft BODY gyroscope bias error, rad/s)
//
// Nominal Navigation State:
// Maintains the 4-component unit quaternion q in S^3 along with 3-vectors:
//     x_nom = { r_I, v_I, q_I_B, b_a_B, b_g_B }
//
// Covariance:
//     P in R^(15x15)

inline constexpr std::size_t k_nav_state_dim = 15;

inline constexpr std::size_t k_idx_pos = 0;
inline constexpr std::size_t k_idx_vel = 3;
inline constexpr std::size_t k_idx_att = 6;
inline constexpr std::size_t k_idx_acc_bias = 9;
inline constexpr std::size_t k_idx_gyro_bias = 12;

using NavigationCovariance = math::Matrix<k_nav_state_dim, k_nav_state_dim>;
using NavigationErrorState = math::Matrix<k_nav_state_dim, 1>;

// 16-component nominal navigation state representation (15 physical/bias DOF + 1 quaternion constraint).
struct NavigationState {
    math::Vector3 position_eci_m{};
    math::Vector3 velocity_eci_mps{};
    math::Quaternion attitude_body_to_eci{math::Quaternion::identity()};
    math::Vector3 accelerometer_bias_body_mps2{};
    math::Vector3 gyro_bias_body_rad_s{};

    [[nodiscard]] bool all_finite() const noexcept {
        return math::is_finite(position_eci_m) &&
               math::is_finite(velocity_eci_mps) &&
               math::is_finite(attitude_body_to_eci) &&
               math::is_finite(accelerometer_bias_body_mps2) &&
               math::is_finite(gyro_bias_body_rad_s);
    }

    [[nodiscard]] orbit::CartesianState to_cartesian_state() const noexcept {
        return orbit::CartesianState{position_eci_m, velocity_eci_mps};
    }
};

// Configuration parameters for the integrated 15-state navigation EKF.
struct IntegratedNavigationConfig {
    double gravitational_parameter_m3_per_s2{constants::earth_gravitational_parameter_m3_per_s2};

    // IMU white noise continuous spectral density parameters
    double accelerometer_noise_std_mps2{1.0e-3}; // drives velocity process noise
    double gyro_noise_std_rad_s{1.0e-3};         // drives attitude process noise

    // Bias random walk parameters (zero for constant bias model)
    double accelerometer_bias_random_walk_mps2_per_sqrt_s{0.0};
    double gyro_bias_random_walk_rad_s_per_sqrt_s{0.0};

    // Sensor measurement noise parameters
    double gnss_position_noise_std_m{10.0};
    double gnss_velocity_noise_std_mps{0.05};
    double star_tracker_noise_std_rad{1.0e-3};
    double range_noise_std_m{1.0};

    // Initial 1-sigma uncertainty parameters
    double initial_position_std_m{50.0};
    double initial_velocity_std_mps{0.5};
    double initial_attitude_std_rad{0.02}; // ~1.15 deg
    double initial_accel_bias_std_mps2{0.01};
    double initial_gyro_bias_std_rad_s{0.005};

    [[nodiscard]] NavigationCovariance initial_covariance() const {
        NavigationCovariance P;
        const double var_r = initial_position_std_m * initial_position_std_m;
        const double var_v = initial_velocity_std_mps * initial_velocity_std_mps;
        const double var_att = initial_attitude_std_rad * initial_attitude_std_rad;
        const double var_ba = initial_accel_bias_std_mps2 * initial_accel_bias_std_mps2;
        const double var_bg = initial_gyro_bias_std_rad_s * initial_gyro_bias_std_rad_s;

        for (std::size_t i = 0; i < 3; ++i) {
            P(k_idx_pos + i, k_idx_pos + i) = var_r;
            P(k_idx_vel + i, k_idx_vel + i) = var_v;
            P(k_idx_att + i, k_idx_att + i) = var_att;
            P(k_idx_acc_bias + i, k_idx_acc_bias + i) = var_ba;
            P(k_idx_gyro_bias + i, k_idx_gyro_bias + i) = var_bg;
        }
        return P;
    }
};

// Analytical 15x15 continuous-time error dynamics Jacobian F:
//
//     F = [  0       I_3      0                    0             0     ]
//         [ G(r)      0   -C(q)[f_B]x          -C(q)             0     ]
//         [  0        0     -[omega_B]x            0           -I_3    ]
//         [  0        0       0                    0             0     ]
//         [  0        0       0                    0             0     ]
//
// where:
//   G(r) = -mu/r^3 * (I_3 - 3 r r^T / r^2) is the gravity gradient tensor in ECI
//   C(q) = C_I_B(q) is the DCM mapping Body to ECI
//   f_B = f_m - b_a is the bias-corrected specific force in the Body frame
//   omega_B = omega_m - b_g is the bias-corrected angular velocity in the Body frame
[[nodiscard]] inline NavigationCovariance integrated_error_dynamics_jacobian(
    const math::Vector3& position_eci_m,
    const math::Quaternion& attitude_body_to_eci,
    const math::Vector3& corrected_specific_force_body_mps2,
    const math::Vector3& corrected_angular_velocity_body_rad_s,
    double gravitational_parameter_m3_per_s2) {

    NavigationCovariance f;

    // Block (0, 1): d(delta_r_dot) / d(delta_v) = I_3
    for (std::size_t i = 0; i < 3; ++i) {
        f(k_idx_pos + i, k_idx_vel + i) = 1.0;
    }

    // Block (1, 0): d(delta_v_dot) / d(delta_r) = G(r)
    const math::Matrix<3, 3> G = gravity_position_jacobian(
        position_eci_m, gravitational_parameter_m3_per_s2);
    for (std::size_t r = 0; r < 3; ++r) {
        for (std::size_t c = 0; c < 3; ++c) {
            f(k_idx_vel + r, k_idx_pos + c) = G(r, c);
        }
    }

    // DCM C_I_B mapping Body vectors into ECI
    const math::Matrix3 dcm = attitude_body_to_eci.to_rotation_matrix();
    math::Matrix<3, 3> C_mat;
    for (std::size_t r = 0; r < 3; ++r) {
        for (std::size_t c = 0; c < 3; ++c) {
            C_mat(r, c) = dcm(r, c);
        }
    }

    // Cross product matrix [f_B]x
    const double fx = corrected_specific_force_body_mps2.x();
    const double fy = corrected_specific_force_body_mps2.y();
    const double fz = corrected_specific_force_body_mps2.z();
    math::Matrix<3, 3> cross_f;
    cross_f(0, 1) = -fz; cross_f(0, 2) = fy;
    cross_f(1, 0) = fz;  cross_f(1, 2) = -fx;
    cross_f(2, 0) = -fy; cross_f(2, 1) = fx;

    // Block (1, 2): d(delta_v_dot) / d(delta_theta) = -C_I_B * [f_B]x
    const math::Matrix<3, 3> att_to_acc = (C_mat * cross_f) * (-1.0);
    for (std::size_t r = 0; r < 3; ++r) {
        for (std::size_t c = 0; c < 3; ++c) {
            f(k_idx_vel + r, k_idx_att + c) = att_to_acc(r, c);
        }
    }

    // Block (1, 3): d(delta_v_dot) / d(delta_b_a) = -C_I_B
    for (std::size_t r = 0; r < 3; ++r) {
        for (std::size_t c = 0; c < 3; ++c) {
            f(k_idx_vel + r, k_idx_acc_bias + c) = -C_mat(r, c);
        }
    }

    // Block (2, 2): d(delta_theta_dot) / d(delta_theta) = -[omega_B]x
    const double wx = corrected_angular_velocity_body_rad_s.x();
    const double wy = corrected_angular_velocity_body_rad_s.y();
    const double wz = corrected_angular_velocity_body_rad_s.z();
    math::Matrix<3, 3> cross_w;
    cross_w(0, 1) = -wz; cross_w(0, 2) = wy;
    cross_w(1, 0) = wz;  cross_w(1, 2) = -wx;
    cross_w(2, 0) = -wy; cross_w(2, 1) = wx;
    const math::Matrix<3, 3> neg_cross_w = cross_w * (-1.0);
    for (std::size_t r = 0; r < 3; ++r) {
        for (std::size_t c = 0; c < 3; ++c) {
            f(k_idx_att + r, k_idx_att + c) = neg_cross_w(r, c);
        }
    }

    // Block (2, 4): d(delta_theta_dot) / d(delta_b_g) = -I_3
    for (std::size_t i = 0; i < 3; ++i) {
        f(k_idx_att + i, k_idx_gyro_bias + i) = -1.0;
    }

    return f;
}

// First-order discrete state transition matrix Phi = I_15 + F * dt.
[[nodiscard]] inline NavigationCovariance integrated_discrete_state_transition(
    const NavigationCovariance& F,
    double dt_s) noexcept {
    NavigationCovariance phi = NavigationCovariance::identity();
    for (std::size_t r = 0; r < k_nav_state_dim; ++r) {
        for (std::size_t c = 0; c < k_nav_state_dim; ++c) {
            phi(r, c) += F(r, c) * dt_s;
        }
    }
    return phi;
}

// Discrete process noise covariance Q in R^(15x15) driven by IMU noise and bias random walk.
[[nodiscard]] inline NavigationCovariance integrated_discrete_process_noise(
    const IntegratedNavigationConfig& config,
    double dt_s) noexcept {
    NavigationCovariance q;
    const double dt2 = dt_s * dt_s;
    const double dt3 = dt2 * dt_s;

    const double var_a = config.accelerometer_noise_std_mps2 * config.accelerometer_noise_std_mps2;
    const double q_rr = (1.0 / 3.0) * var_a * dt3;
    const double q_rv = 0.5 * var_a * dt2;
    const double q_vv = var_a * dt_s;

    const double var_g = config.gyro_noise_std_rad_s * config.gyro_noise_std_rad_s;
    const double q_theta = var_g * dt_s;

    const double var_rw_ba = config.accelerometer_bias_random_walk_mps2_per_sqrt_s *
                             config.accelerometer_bias_random_walk_mps2_per_sqrt_s;
    const double q_ba = var_rw_ba * dt_s;

    const double var_rw_bg = config.gyro_bias_random_walk_rad_s_per_sqrt_s *
                             config.gyro_bias_random_walk_rad_s_per_sqrt_s;
    const double q_bg = var_rw_bg * dt_s;

    for (std::size_t i = 0; i < 3; ++i) {
        q(k_idx_pos + i, k_idx_pos + i) = q_rr;
        q(k_idx_pos + i, k_idx_vel + i) = q_rv;
        q(k_idx_vel + i, k_idx_pos + i) = q_rv;
        q(k_idx_vel + i, k_idx_vel + i) = q_vv;

        q(k_idx_att + i, k_idx_att + i) = q_theta;
        q(k_idx_acc_bias + i, k_idx_acc_bias + i) = q_ba;
        q(k_idx_gyro_bias + i, k_idx_gyro_bias + i) = q_bg;
    }
    return q;
}

// 15x15 covariance reset Jacobian for coordinate transformation following attitude correction delta_theta:
//
//     J_reset = diag(I_3, I_3, I_3 - 0.5 * [delta_theta]x, I_3, I_3)
[[nodiscard]] inline NavigationCovariance integrated_covariance_reset_jacobian(
    const math::Vector3& delta_theta) noexcept {
    NavigationCovariance J = NavigationCovariance::identity();

    const double tx = delta_theta.x();
    const double ty = delta_theta.y();
    const double tz = delta_theta.z();

    // -0.5 * [delta_theta]x
    J(k_idx_att + 0, k_idx_att + 1) = 0.5 * tz;
    J(k_idx_att + 0, k_idx_att + 2) = -0.5 * ty;
    J(k_idx_att + 1, k_idx_att + 0) = -0.5 * tz;
    J(k_idx_att + 1, k_idx_att + 2) = 0.5 * tx;
    J(k_idx_att + 2, k_idx_att + 0) = 0.5 * ty;
    J(k_idx_att + 2, k_idx_att + 1) = -0.5 * tx;

    return J;
}

// Applies the first-order covariance reset transformation: P <- J_reset * P * J_reset^T.
[[nodiscard]] inline NavigationCovariance integrated_covariance_reset(
    const NavigationCovariance& covariance,
    const math::Vector3& delta_theta) noexcept {
    const NavigationCovariance J = integrated_covariance_reset_jacobian(delta_theta);
    return (J * covariance * J.transpose()).symmetrized();
}

// Measurement Jacobians in the 15-state ordering:

// GNSS position and velocity measurement matrix (6 x 15):
//     H_gnss = [ I_3  0    0  0  0 ]
//              [  0  I_3   0  0  0 ]
[[nodiscard]] inline math::Matrix<6, k_nav_state_dim> gnss_measurement_jacobian_15() noexcept {
    math::Matrix<6, k_nav_state_dim> H;
    for (std::size_t i = 0; i < 3; ++i) {
        H(i, k_idx_pos + i) = 1.0;
        H(3 + i, k_idx_vel + i) = 1.0;
    }
    return H;
}

// Star tracker attitude residual measurement matrix (3 x 15):
//     H_st = [ 0  0  I_3  0  0 ]
[[nodiscard]] inline math::Matrix<3, k_nav_state_dim> star_tracker_measurement_jacobian_15() noexcept {
    math::Matrix<3, k_nav_state_dim> H;
    for (std::size_t i = 0; i < 3; ++i) {
        H(i, k_idx_att + i) = 1.0;
    }
    return H;
}

// Nonlinear scalar range measurement Jacobian (1 x 15):
//     H_range = [ -u_LOS^T  0  0  0  0 ]
[[nodiscard]] inline math::Matrix<1, k_nav_state_dim> range_measurement_jacobian_15(
    const math::Vector3& sc_position_eci_m,
    const math::Vector3& target_position_eci_m) {
    const math::Vector3 rel_pos = target_position_eci_m - sc_position_eci_m;
    const double sep = rel_pos.norm();
    if (!std::isfinite(sep) || sep < 1.0e-6) {
        throw std::domain_error(
            "range_measurement_jacobian_15: target and spacecraft positions are coincident "
            "(separation < 1e-6 m); line-of-sight direction is undefined");
    }
    const math::Vector3 u_los = rel_pos / sep;

    math::Matrix<1, k_nav_state_dim> H;
    H(0, k_idx_pos + 0) = -u_los.x();
    H(0, k_idx_pos + 1) = -u_los.y();
    H(0, k_idx_pos + 2) = -u_los.z();
    return H;
}

// Diagnostic containers for sequential measurement updates:
struct GnssUpdate15Diagnostics {
    math::ColVector<6> innovation{};
    double normalized_innovation_squared{0.0};
    bool accepted{false};
};

struct StarTrackerUpdate15Diagnostics {
    math::Vector3 innovation_rad{};
    double normalized_innovation_squared{0.0};
    bool accepted{false};
};

struct RangeUpdate15Diagnostics {
    double predicted_range_m{0.0};
    double innovation_m{0.0};
    double normalized_innovation_squared{0.0};
    bool accepted{false};
};

// ============================================================================
// Production 15-State Integrated Navigation Filter Class
// ============================================================================
class IntegratedNavigationEkf {
public:
    IntegratedNavigationEkf(
        const IntegratedNavigationConfig& config,
        const NavigationState& initial_state,
        const std::optional<NavigationCovariance>& initial_covariance = std::nullopt)
        : config_(config),
          state_(initial_state) {
        if (!state_.all_finite()) {
            throw std::domain_error("IntegratedNavigationEkf: initial state must contain only finite values");
        }
        state_.attitude_body_to_eci = state_.attitude_body_to_eci.normalized();
        if (initial_covariance.has_value()) {
            covariance_ = *initial_covariance;
        } else {
            covariance_ = config_.initial_covariance();
        }
        validate_covariance(covariance_, "IntegratedNavigationEkf: initial covariance");
    }

    [[nodiscard]] const NavigationState& state() const noexcept { return state_; }
    [[nodiscard]] const NavigationCovariance& covariance() const noexcept { return covariance_; }
    [[nodiscard]] const IntegratedNavigationConfig& config() const noexcept { return config_; }

    [[nodiscard]] const GnssUpdate15Diagnostics& last_gnss_diagnostics() const noexcept {
        return last_gnss_diag_;
    }
    [[nodiscard]] const StarTrackerUpdate15Diagnostics& last_star_tracker_diagnostics() const noexcept {
        return last_st_diag_;
    }
    [[nodiscard]] const RangeUpdate15Diagnostics& last_range_diagnostics() const noexcept {
        return last_range_diag_;
    }

    // Time update (prediction) driven by IMU measurement across actual interval dt_s:
    // Reconstructs inertial acceleration:
    //     a_I = C_I_B(q)(f_m - b_a) + g(r_I)
    // Propagates nominal state via kinematics, and advances 15x15 covariance via:
    //     P^- = Phi P^+ Phi^T + Q
    void predict(
        const math::Vector3& measured_specific_force_body_mps2,
        const math::Vector3& measured_angular_velocity_body_rad_s,
        double dt_s) {
        if (!std::isfinite(dt_s) || dt_s <= 0.0) {
            throw std::domain_error("IntegratedNavigationEkf: dt must be finite and positive");
        }
        if (!math::is_finite(measured_specific_force_body_mps2) ||
            !math::is_finite(measured_angular_velocity_body_rad_s)) {
            throw std::domain_error("IntegratedNavigationEkf: IMU measurements must be finite");
        }

        // 1. Bias correction
        const math::Vector3 f_b = measured_specific_force_body_mps2 - state_.accelerometer_bias_body_mps2;
        const math::Vector3 w_b = measured_angular_velocity_body_rad_s - state_.gyro_bias_body_rad_s;

        // 2. Inertial acceleration reconstruction
        const math::Vector3 f_eci = state_.attitude_body_to_eci.rotate_vector(f_b);
        const math::Vector3 g_eci = dynamics::two_body_acceleration(
            state_.position_eci_m, config_.gravitational_parameter_m3_per_s2);
        const math::Vector3 a_eci = f_eci + g_eci;

        // 3. Nominal state propagation across interval dt_s:
        // Position and velocity updated via midpoint/Euler with reconstructed acceleration
        const math::Vector3 next_pos = state_.position_eci_m + state_.velocity_eci_mps * dt_s + a_eci * (0.5 * dt_s * dt_s);
        const math::Vector3 next_vel = state_.velocity_eci_mps + a_eci * dt_s;

        // Attitude quaternion propagated via angular velocity w_b
        const double w_norm = w_b.norm();
        const double angle = w_norm * dt_s;
        math::Quaternion dq;
        if (angle < 1.0e-8) {
            dq = math::Quaternion{1.0, 0.5 * w_b.x() * dt_s, 0.5 * w_b.y() * dt_s, 0.5 * w_b.z() * dt_s};
        } else {
            dq = math::Quaternion::from_axis_angle(w_b / w_norm, angle);
        }
        const math::Quaternion next_q = (state_.attitude_body_to_eci * dq).normalized();

        state_.position_eci_m = next_pos;
        state_.velocity_eci_mps = next_vel;
        state_.attitude_body_to_eci = next_q;
        // Biases remain constant during prediction

        // 4. Analytical Error Dynamics Jacobian F and Transition Phi
        const NavigationCovariance F = integrated_error_dynamics_jacobian(
            state_.position_eci_m,
            state_.attitude_body_to_eci,
            f_b,
            w_b,
            config_.gravitational_parameter_m3_per_s2);

        const NavigationCovariance phi = integrated_discrete_state_transition(F, dt_s);
        const NavigationCovariance Q = integrated_discrete_process_noise(config_, dt_s);

        // 5. Covariance propagation P^- = Phi P Phi^T + Q
        covariance_ = (phi * covariance_ * phi.transpose() + Q).symmetrized();
        validate_covariance(covariance_, "IntegratedNavigationEkf: post-prediction covariance");
    }

    // GNSS measurement update (6-DOF position and velocity in ECI):
    bool update_gnss(const sensors::GnssMeasurement& meas) {
        if (!meas.valid) {
            last_gnss_diag_.accepted = false;
            return false;
        }

        // 1. Innovation vector in R^6
        math::ColVector<6> z;
        z(0, 0) = meas.position_eci_m.x();
        z(1, 0) = meas.position_eci_m.y();
        z(2, 0) = meas.position_eci_m.z();
        z(3, 0) = meas.velocity_eci_mps.x();
        z(4, 0) = meas.velocity_eci_mps.y();
        z(5, 0) = meas.velocity_eci_mps.z();

        math::ColVector<6> h_x;
        h_x(0, 0) = state_.position_eci_m.x();
        h_x(1, 0) = state_.position_eci_m.y();
        h_x(2, 0) = state_.position_eci_m.z();
        h_x(3, 0) = state_.velocity_eci_mps.x();
        h_x(4, 0) = state_.velocity_eci_mps.y();
        h_x(5, 0) = state_.velocity_eci_mps.z();

        const math::ColVector<6> y = z - h_x;

        // 2. Measurement noise R
        math::Matrix<6, 6> R;
        const double var_r = config_.gnss_position_noise_std_m * config_.gnss_position_noise_std_m;
        const double var_v = config_.gnss_velocity_noise_std_mps * config_.gnss_velocity_noise_std_mps;
        for (std::size_t i = 0; i < 3; ++i) {
            R(i, i) = var_r;
            R(3 + i, 3 + i) = var_v;
        }

        const math::Matrix<6, k_nav_state_dim> H = gnss_measurement_jacobian_15();

        // 3. S = H P H^T + R and Kalman Gain K = P H^T S^{-1}
        const math::Matrix<6, 6> S = (H * covariance_ * H.transpose() + R).symmetrized();
        const math::Matrix<k_nav_state_dim, 6> PHT = covariance_ * H.transpose();
        const math::Matrix<k_nav_state_dim, 6> K = math::solve_spd_multi(S, PHT.transpose()).transpose();

        // 4. Correction delta_x = K y
        const NavigationErrorState delta_x = K * y;

        // 5. Joseph-form covariance update
        const NavigationCovariance I15 = NavigationCovariance::identity();
        const NavigationCovariance I_minus_KH = I15 - (K * H);
        covariance_ = (I_minus_KH * covariance_ * I_minus_KH.transpose() + K * R * K.transpose()).symmetrized();

        // 6. Incorporate corrections into nominal state
        state_.position_eci_m = state_.position_eci_m + math::Vector3{
            delta_x(k_idx_pos + 0, 0),
            delta_x(k_idx_pos + 1, 0),
            delta_x(k_idx_pos + 2, 0)
        };

        state_.velocity_eci_mps = state_.velocity_eci_mps + math::Vector3{
            delta_x(k_idx_vel + 0, 0),
            delta_x(k_idx_vel + 1, 0),
            delta_x(k_idx_vel + 2, 0)
        };

        const math::Vector3 delta_theta{
            delta_x(k_idx_att + 0, 0),
            delta_x(k_idx_att + 1, 0),
            delta_x(k_idx_att + 2, 0)
        };
        const double th_norm = delta_theta.norm();
        if (th_norm > 0.0) {
            const math::Quaternion dq = math::Quaternion::from_axis_angle(delta_theta / th_norm, th_norm);
            state_.attitude_body_to_eci = (state_.attitude_body_to_eci * dq).normalized();
            covariance_ = integrated_covariance_reset(covariance_, delta_theta);
        }

        state_.accelerometer_bias_body_mps2 = state_.accelerometer_bias_body_mps2 + math::Vector3{
            delta_x(k_idx_acc_bias + 0, 0),
            delta_x(k_idx_acc_bias + 1, 0),
            delta_x(k_idx_acc_bias + 2, 0)
        };

        state_.gyro_bias_body_rad_s = state_.gyro_bias_body_rad_s + math::Vector3{
            delta_x(k_idx_gyro_bias + 0, 0),
            delta_x(k_idx_gyro_bias + 1, 0),
            delta_x(k_idx_gyro_bias + 2, 0)
        };

        // 7. Diagnostics
        const math::ColVector<6> Sinv_y = math::solve_spd(S, y);
        double nis = 0.0;
        for (std::size_t i = 0; i < 6; ++i) {
            nis += y(i, 0) * Sinv_y(i, 0);
        }

        last_gnss_diag_.innovation = y;
        last_gnss_diag_.normalized_innovation_squared = nis;
        last_gnss_diag_.accepted = true;

        validate_covariance(covariance_, "IntegratedNavigationEkf: post-GNSS covariance");
        return true;
    }

    // Star tracker measurement update (attitude error residual in Body frame with double-cover alignment):
    bool update_star_tracker(const sensors::StarTrackerMeasurement& meas) {
        if (!meas.valid) {
            last_st_diag_.accepted = false;
            return false;
        }

        // 1. Antipodal sign alignment
        math::Quaternion q_m = meas.orientation_eci_from_body.normalized();
        const double dot = state_.attitude_body_to_eci.w() * q_m.w() +
                           state_.attitude_body_to_eci.x() * q_m.x() +
                           state_.attitude_body_to_eci.y() * q_m.y() +
                           state_.attitude_body_to_eci.z() * q_m.z();
        if (dot < 0.0) {
            q_m = math::Quaternion{-q_m.w(), -q_m.x(), -q_m.y(), -q_m.z()};
        }

        // 2. Physical error quaternion q_err = q_nom^* * q_m
        const math::Quaternion q_err = state_.attitude_body_to_eci.conjugate() * q_m;
        const math::Vector3 inno_vec{2.0 * q_err.x(), 2.0 * q_err.y(), 2.0 * q_err.z()};

        math::ColVector<3> y;
        y(0, 0) = inno_vec.x();
        y(1, 0) = inno_vec.y();
        y(2, 0) = inno_vec.z();

        // 3. Measurement noise R (isotropic per-axis variance = sigma^2 / 3)
        const double var_st = (config_.star_tracker_noise_std_rad * config_.star_tracker_noise_std_rad) / 3.0;
        math::Matrix<3, 3> R;
        R(0, 0) = var_st; R(1, 1) = var_st; R(2, 2) = var_st;

        const math::Matrix<3, k_nav_state_dim> H = star_tracker_measurement_jacobian_15();

        // 4. S = H P H^T + R and Kalman Gain K = P H^T S^{-1}
        const math::Matrix<3, 3> S = (H * covariance_ * H.transpose() + R).symmetrized();
        const math::Matrix<k_nav_state_dim, 3> PHT = covariance_ * H.transpose();
        const math::Matrix<k_nav_state_dim, 3> K = math::solve_spd_multi(S, PHT.transpose()).transpose();

        // 5. Correction delta_x = K y
        const NavigationErrorState delta_x = K * y;

        // 6. Joseph-form covariance update
        const NavigationCovariance I15 = NavigationCovariance::identity();
        const NavigationCovariance I_minus_KH = I15 - (K * H);
        covariance_ = (I_minus_KH * covariance_ * I_minus_KH.transpose() + K * R * K.transpose()).symmetrized();

        // 7. Incorporate corrections into nominal state
        state_.position_eci_m = state_.position_eci_m + math::Vector3{
            delta_x(k_idx_pos + 0, 0),
            delta_x(k_idx_pos + 1, 0),
            delta_x(k_idx_pos + 2, 0)
        };

        state_.velocity_eci_mps = state_.velocity_eci_mps + math::Vector3{
            delta_x(k_idx_vel + 0, 0),
            delta_x(k_idx_vel + 1, 0),
            delta_x(k_idx_vel + 2, 0)
        };

        const math::Vector3 delta_theta{
            delta_x(k_idx_att + 0, 0),
            delta_x(k_idx_att + 1, 0),
            delta_x(k_idx_att + 2, 0)
        };
        const double th_norm = delta_theta.norm();
        if (th_norm > 0.0) {
            const math::Quaternion dq = math::Quaternion::from_axis_angle(delta_theta / th_norm, th_norm);
            state_.attitude_body_to_eci = (state_.attitude_body_to_eci * dq).normalized();
            covariance_ = integrated_covariance_reset(covariance_, delta_theta);
        }

        state_.accelerometer_bias_body_mps2 = state_.accelerometer_bias_body_mps2 + math::Vector3{
            delta_x(k_idx_acc_bias + 0, 0),
            delta_x(k_idx_acc_bias + 1, 0),
            delta_x(k_idx_acc_bias + 2, 0)
        };

        state_.gyro_bias_body_rad_s = state_.gyro_bias_body_rad_s + math::Vector3{
            delta_x(k_idx_gyro_bias + 0, 0),
            delta_x(k_idx_gyro_bias + 1, 0),
            delta_x(k_idx_gyro_bias + 2, 0)
        };

        // 8. Diagnostics
        const math::ColVector<3> Sinv_y = math::solve_spd(S, y);
        double nis = 0.0;
        for (std::size_t i = 0; i < 3; ++i) {
            nis += y(i, 0) * Sinv_y(i, 0);
        }

        last_st_diag_.innovation_rad = inno_vec;
        last_st_diag_.normalized_innovation_squared = nis;
        last_st_diag_.accepted = true;

        validate_covariance(covariance_, "IntegratedNavigationEkf: post-StarTracker covariance");
        return true;
    }

    // Scalar relative range measurement update:
    bool update_range(
        double measured_range_m,
        const math::Vector3& target_position_eci_m,
        double range_noise_std_m) {
        if (!std::isfinite(measured_range_m) || measured_range_m <= 0.0) {
            last_range_diag_.accepted = false;
            return false;
        }

        const double pred_range = (target_position_eci_m - state_.position_eci_m).norm();
        if (pred_range < 1.0e-6) {
            throw std::domain_error("update_range: coincident geometry (separation < 1e-6 m)");
        }

        const double innovation = measured_range_m - pred_range;
        const auto H = range_measurement_jacobian_15(state_.position_eci_m, target_position_eci_m);

        const double R = range_noise_std_m * range_noise_std_m;

        // Scalar innovation variance S = H P H^T + R
        double HPH = 0.0;
        for (std::size_t r = 0; r < k_nav_state_dim; ++r) {
            if (H(0, r) == 0.0) continue;
            for (std::size_t c = 0; c < k_nav_state_dim; ++c) {
                if (H(0, c) == 0.0) continue;
                HPH += H(0, r) * covariance_(r, c) * H(0, c);
            }
        }
        const double S = HPH + R;
        if (S <= 0.0 || !std::isfinite(S)) {
            throw std::domain_error("update_range: non-positive or non-finite innovation variance");
        }

        // Kalman gain K = P H^T / S (15 x 1)
        NavigationErrorState K;
        for (std::size_t i = 0; i < k_nav_state_dim; ++i) {
            double PHT_i = 0.0;
            for (std::size_t j = 0; j < 3; ++j) { // H only non-zero on position (0..2)
                PHT_i += covariance_(i, k_idx_pos + j) * H(0, k_idx_pos + j);
            }
            K(i, 0) = PHT_i / S;
        }

        // State correction delta_x = K * innovation
        const NavigationErrorState delta_x = K * innovation;

        // Joseph-form covariance update: P^+ = (I - K H) P (I - K H)^T + K R K^T
        const NavigationCovariance I15 = NavigationCovariance::identity();
        const NavigationCovariance I_minus_KH = I15 - (K * H);
        covariance_ = (I_minus_KH * covariance_ * I_minus_KH.transpose() + (K * K.transpose()) * R).symmetrized();

        // Incorporate corrections
        state_.position_eci_m = state_.position_eci_m + math::Vector3{
            delta_x(k_idx_pos + 0, 0),
            delta_x(k_idx_pos + 1, 0),
            delta_x(k_idx_pos + 2, 0)
        };

        state_.velocity_eci_mps = state_.velocity_eci_mps + math::Vector3{
            delta_x(k_idx_vel + 0, 0),
            delta_x(k_idx_vel + 1, 0),
            delta_x(k_idx_vel + 2, 0)
        };

        const math::Vector3 delta_theta{
            delta_x(k_idx_att + 0, 0),
            delta_x(k_idx_att + 1, 0),
            delta_x(k_idx_att + 2, 0)
        };
        const double th_norm = delta_theta.norm();
        if (th_norm > 0.0) {
            const math::Quaternion dq = math::Quaternion::from_axis_angle(delta_theta / th_norm, th_norm);
            state_.attitude_body_to_eci = (state_.attitude_body_to_eci * dq).normalized();
            covariance_ = integrated_covariance_reset(covariance_, delta_theta);
        }

        state_.accelerometer_bias_body_mps2 = state_.accelerometer_bias_body_mps2 + math::Vector3{
            delta_x(k_idx_acc_bias + 0, 0),
            delta_x(k_idx_acc_bias + 1, 0),
            delta_x(k_idx_acc_bias + 2, 0)
        };

        state_.gyro_bias_body_rad_s = state_.gyro_bias_body_rad_s + math::Vector3{
            delta_x(k_idx_gyro_bias + 0, 0),
            delta_x(k_idx_gyro_bias + 1, 0),
            delta_x(k_idx_gyro_bias + 2, 0)
        };

        const double nis = (innovation * innovation) / S;

        last_range_diag_.predicted_range_m = pred_range;
        last_range_diag_.innovation_m = innovation;
        last_range_diag_.normalized_innovation_squared = nis;
        last_range_diag_.accepted = true;

        validate_covariance(covariance_, "IntegratedNavigationEkf: post-Range covariance");
        return true;
    }

    // Explicit reset:
    void reset(
        const NavigationState& state,
        const std::optional<NavigationCovariance>& covariance = std::nullopt) {
        state_ = state;
        state_.attitude_body_to_eci = state_.attitude_body_to_eci.normalized();
        if (covariance.has_value()) {
            covariance_ = *covariance;
        } else {
            covariance_ = config_.initial_covariance();
        }
        last_gnss_diag_ = GnssUpdate15Diagnostics{};
        last_st_diag_ = StarTrackerUpdate15Diagnostics{};
        last_range_diag_ = RangeUpdate15Diagnostics{};
        validate_covariance(covariance_, "IntegratedNavigationEkf: reset covariance");
    }

private:
    IntegratedNavigationConfig config_;
    NavigationState state_;
    NavigationCovariance covariance_;

    GnssUpdate15Diagnostics last_gnss_diag_{};
    StarTrackerUpdate15Diagnostics last_st_diag_{};
    RangeUpdate15Diagnostics last_range_diag_{};
};

} // namespace astradock::estimation
