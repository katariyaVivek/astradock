#pragma once

#include "math/vector3.hpp"
#include "sensors/sensor_common.hpp"
#include "spacecraft/spacecraft_state.hpp"

#include <cmath>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <vector>

namespace astradock::sensors {

// Configuration parameters for the GNSS receiver sensor simulation.
struct GnssConfig {
    // Position error parameters (ECI frame, m)
    math::Vector3 position_bias_eci_m{};
    double position_noise_std_m{0.0};

    // Velocity error parameters (ECI frame, m/s)
    math::Vector3 velocity_bias_eci_mps{};
    double velocity_noise_std_mps{0.0};

    // Timing & sampling parameters (nominal 1 Hz = 1.0 s)
    double sample_period_s{1.0};
    std::vector<DropoutWindow> dropouts{};
    std::uint64_t random_seed{0};

    constexpr GnssConfig() noexcept = default;
};

// GNSS measurement packet.
// Contains measured Cartesian position and velocity in Earth-Centered Inertial (ECI) coordinates.
struct GnssMeasurement {
    double timestamp_s{0.0};
    math::Vector3 position_eci_m{};
    math::Vector3 velocity_eci_mps{};
    bool valid{true};
};

// GNSS Receiver sensor simulation.
//
// Models an idealized spaceborne GNSS receiver providing absolute Cartesian
// position and velocity solutions in the Earth-Centered Inertial (ECI) coordinate frame:
//
//   r_meas = r_true_ECI + b_r + n_r   (m, ECI frame)
//   v_meas = v_true_ECI + b_v + n_v   (m/s, ECI frame)
//
// where:
//   b_r, b_v = constant receiver position and velocity biases
//   n_r ~ N(0, sigma_r^2 I), n_v ~ N(0, sigma_v^2 I) = independent isotropic Gaussian noise
class GnssSensor {
public:
    explicit GnssSensor(const GnssConfig& config)
        : config_(config),
          rng_(config.random_seed),
          schedule_(config.sample_period_s, 0.0) {
        validate_config(config);
    }

    [[nodiscard]] const GnssConfig& config() const noexcept {
        return config_;
    }

    // Resets random number generator and sampling schedule to initial epoch.
    void reset(double start_time_s = 0.0) noexcept {
        rng_.reset();
        schedule_.reset(start_time_s);
    }

    // Evaluates instantaneous measurement without modifying sampling schedule.
    [[nodiscard]] GnssMeasurement measure(
        double timestamp_s,
        const spacecraft::SpacecraftState& truth_state) {
        if (!spacecraft::is_finite(truth_state)) {
            throw std::domain_error("GNSS truth state must be finite");
        }

        const bool in_dropout = is_in_dropout_window(timestamp_s, config_.dropouts);
        if (in_dropout) {
            return GnssMeasurement{
                timestamp_s,
                math::Vector3{},
                math::Vector3{},
                false
            };
        }

        const math::Vector3 pos_noise = rng_.gaussian_vector3(config_.position_noise_std_m);
        const math::Vector3 vel_noise = rng_.gaussian_vector3(config_.velocity_noise_std_mps);

        const math::Vector3 measured_pos = truth_state.position() + config_.position_bias_eci_m + pos_noise;
        const math::Vector3 measured_vel = truth_state.velocity() + config_.velocity_bias_eci_mps + vel_noise;

        return GnssMeasurement{
            timestamp_s,
            measured_pos,
            measured_vel,
            true
        };
    }

    // Advances sampling schedule and returns a sample if scheduled at current_t_s.
    [[nodiscard]] std::optional<GnssMeasurement> sample(
        double current_t_s,
        const spacecraft::SpacecraftState& truth_state) {
        if (!schedule_.should_sample(current_t_s)) {
            return std::nullopt;
        }

        const double sample_time_s = schedule_.next_sample_time_s();
        schedule_.advance(current_t_s);

        return measure(sample_time_s, truth_state);
    }

private:
    static void validate_config(const GnssConfig& cfg) {
        if (!std::isfinite(cfg.sample_period_s) || cfg.sample_period_s <= 0.0) {
            throw std::domain_error("GNSS sample period must be positive and finite");
        }
        if (!math::is_finite(cfg.position_bias_eci_m) || !math::is_finite(cfg.velocity_bias_eci_mps)) {
            throw std::domain_error("GNSS biases must be finite");
        }
        if (!std::isfinite(cfg.position_noise_std_m) || cfg.position_noise_std_m < 0.0 ||
            !std::isfinite(cfg.velocity_noise_std_mps) || cfg.velocity_noise_std_mps < 0.0) {
            throw std::domain_error("GNSS noise standard deviations must be finite and non-negative");
        }
    }

    GnssConfig config_;
    DeterministicRng rng_;
    SensorSchedule schedule_;
};

}  // namespace astradock::sensors
