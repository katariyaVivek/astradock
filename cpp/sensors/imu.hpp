#pragma once

#include "math/quaternion.hpp"
#include "math/vector3.hpp"
#include "sensors/sensor_common.hpp"
#include "spacecraft/spacecraft_state.hpp"

#include <cmath>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <vector>

namespace astradock::sensors {

// Configuration parameters for the 6-axis Inertial Measurement Unit (IMU).
struct ImuConfig {
    // Gyroscope parameters (BODY frame, rad/s)
    math::Vector3 gyro_bias_rad_s{};
    double gyro_noise_std_rad_s{0.0};

    // Accelerometer parameters (BODY frame, m/s^2)
    math::Vector3 accel_bias_mps2{};
    double accel_noise_std_mps2{0.0};

    // Timing & sampling parameters (nominal 100 Hz = 0.01 s)
    double sample_period_s{0.01};
    std::vector<DropoutWindow> dropouts{};
    std::uint64_t random_seed{0};

    constexpr ImuConfig() noexcept = default;
};

// 6-axis IMU measurement packet.
// Contains specific force (accelerometer) and angular velocity (gyroscope),
// both expressed in the Spacecraft Body coordinate frame.
struct ImuMeasurement {
    double timestamp_s{0.0};
    math::Vector3 specific_force_body_mps2{};
    math::Vector3 angular_velocity_body_rad_s{};
    bool valid{true};
};

// Inertial Measurement Unit (IMU) sensor simulation.
//
// Models a 3-axis gyroscope and a 3-axis accelerometer rigidly mounted to the
// spacecraft structure at the center of mass:
//
// 1. Gyroscope: Measures the angular velocity of the spacecraft body frame
//    relative to the inertial frame, expressed in the body frame:
//      omega_meas = omega_true_body + b_g + n_g   (rad/s, BODY frame)
//
// 2. Accelerometer: Measures specific force (non-gravitational acceleration):
//      f = a_inertial - g_gravity
//    In pure gravitational free fall (two-body orbit), a_inertial = g_gravity,
//    so the ideal specific force is zero. Non-gravitational perturbations
//    (e.g., atmospheric drag, solar radiation pressure, thruster firings)
//    produce nonzero specific force:
//      f_I = a_non_grav_I + F_thrust_I / m
//      f_B = C_{B_I}(q) * f_I = q^* ⊗ f_I ⊗ q   (m/s^2, BODY frame)
//      f_meas = f_B + b_a + n_a                 (m/s^2, BODY frame)
class ImuSensor {
public:
    explicit ImuSensor(const ImuConfig& config)
        : config_(config),
          rng_(config.random_seed),
          schedule_(config.sample_period_s, 0.0) {
        validate_config(config);
    }

    [[nodiscard]] const ImuConfig& config() const noexcept {
        return config_;
    }

    // Resets random number generator and sampling schedule to initial epoch.
    void reset(double start_time_s = 0.0) noexcept {
        rng_.reset();
        schedule_.reset(start_time_s);
    }

    // Evaluates instantaneous measurement without modifying sampling schedule.
    [[nodiscard]] ImuMeasurement measure(
        double timestamp_s,
        const spacecraft::SpacecraftState& truth_state,
        const math::Vector3& non_grav_accel_eci_mps2 = math::Vector3{}) {
        if (!spacecraft::is_finite(truth_state)) {
            throw std::domain_error("IMU truth state must be finite");
        }
        if (!math::is_finite(non_grav_accel_eci_mps2)) {
            throw std::domain_error("IMU non-gravitational acceleration must be finite");
        }

        const bool in_dropout = is_in_dropout_window(timestamp_s, config_.dropouts);
        if (in_dropout) {
            return ImuMeasurement{
                timestamp_s,
                math::Vector3{},
                math::Vector3{},
                false
            };
        }

        // 1. Gyroscope measurement: omega_meas = omega_true + b_g + n_g (BODY)
        const math::Vector3 true_omega_body = truth_state.angular_velocity_rad_per_s();
        const math::Vector3 gyro_noise = rng_.gaussian_vector3(config_.gyro_noise_std_rad_s);
        const math::Vector3 measured_omega_body = true_omega_body + config_.gyro_bias_rad_s + gyro_noise;

        // 2. Accelerometer specific force measurement:
        //    Ideal specific force in ECI: f_I = non-gravitational acceleration
        //    Transform to Body: f_B = q* ⊗ f_I ⊗ q
        const math::Quaternion q_eci_from_body = truth_state.orientation().normalized();
        const math::Vector3 true_f_body = q_eci_from_body.conjugate().rotate_vector(non_grav_accel_eci_mps2);
        const math::Vector3 accel_noise = rng_.gaussian_vector3(config_.accel_noise_std_mps2);
        const math::Vector3 measured_f_body = true_f_body + config_.accel_bias_mps2 + accel_noise;

        return ImuMeasurement{
            timestamp_s,
            measured_f_body,
            measured_omega_body,
            true
        };
    }

    // Advances the sampling schedule and returns a new sample if scheduled at current_t_s.
    [[nodiscard]] std::optional<ImuMeasurement> sample(
        double current_t_s,
        const spacecraft::SpacecraftState& truth_state,
        const math::Vector3& non_grav_accel_eci_mps2 = math::Vector3{}) {
        if (!schedule_.should_sample(current_t_s)) {
            return std::nullopt;
        }

        const double sample_time_s = schedule_.next_sample_time_s();
        schedule_.advance(current_t_s);

        return measure(sample_time_s, truth_state, non_grav_accel_eci_mps2);
    }

private:
    static void validate_config(const ImuConfig& cfg) {
        if (!std::isfinite(cfg.sample_period_s) || cfg.sample_period_s <= 0.0) {
            throw std::domain_error("IMU sample period must be positive and finite");
        }
        if (!math::is_finite(cfg.gyro_bias_rad_s) || !math::is_finite(cfg.accel_bias_mps2)) {
            throw std::domain_error("IMU biases must be finite");
        }
        if (!std::isfinite(cfg.gyro_noise_std_rad_s) || cfg.gyro_noise_std_rad_s < 0.0 ||
            !std::isfinite(cfg.accel_noise_std_mps2) || cfg.accel_noise_std_mps2 < 0.0) {
            throw std::domain_error("IMU noise standard deviations must be finite and non-negative");
        }
    }

    ImuConfig config_;
    DeterministicRng rng_;
    SensorSchedule schedule_;
};

}  // namespace astradock::sensors
