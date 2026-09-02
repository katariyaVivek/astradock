#pragma once

#include "math/vector3.hpp"
#include "sensors/sensor_common.hpp"
#include "spacecraft/spacecraft_state.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <vector>

namespace astradock::sensors {

// Configuration parameters for the Range Sensor simulation (e.g., laser / radar / radio rangefinder).
struct RangeSensorConfig {
    // Constant range measurement bias (m)
    double bias_m{0.0};

    // 1-sigma Gaussian range noise (m)
    double noise_std_m{0.0};

    // Timing & sampling parameters (nominal 10 Hz = 0.1 s)
    double sample_period_s{0.1};
    std::vector<DropoutWindow> dropouts{};
    std::uint64_t random_seed{0};

    constexpr RangeSensorConfig() noexcept = default;
};

// Range measurement packet.
// Contains measured scalar distance between spacecraft and target.
struct RangeMeasurement {
    double timestamp_s{0.0};
    double range_m{0.0};
    bool valid{true};
};

// Relative Range Sensor simulation.
//
// Models a line-of-sight range measurement (e.g. lidar, radar, optical time-of-flight,
// or RF ranging) between the chaser spacecraft and a target location/spacecraft:
//
//   rho_true = ||r_target_ECI - r_spacecraft_ECI||   (m)
//   rho_meas = max(0, rho_true + b_rho + n_rho)      (m)
//
// Properties:
//   - Frame Invariant: Translation of both bodies by an identical vector c leaves
//     the range measurement unchanged: ||(r_target + c) - (r_sc + c)|| = ||r_target - r_sc||.
//   - Non-Negative: Physical range is bounded below by zero.
class RangeSensor {
public:
    explicit RangeSensor(const RangeSensorConfig& config)
        : config_(config),
          rng_(config.random_seed),
          schedule_(config.sample_period_s, 0.0) {
        validate_config(config);
    }

    [[nodiscard]] const RangeSensorConfig& config() const noexcept {
        return config_;
    }

    // Resets random number generator and sampling schedule to initial epoch.
    void reset(double start_time_s = 0.0) noexcept {
        rng_.reset();
        schedule_.reset(start_time_s);
    }

    // Evaluates instantaneous measurement without modifying sampling schedule.
    [[nodiscard]] RangeMeasurement measure(
        double timestamp_s,
        const spacecraft::SpacecraftState& truth_state,
        const math::Vector3& target_position_eci_m) {
        if (!spacecraft::is_finite(truth_state)) {
            throw std::domain_error("Range sensor truth state must be finite");
        }
        if (!math::is_finite(target_position_eci_m)) {
            throw std::domain_error("Target position must be finite");
        }

        const bool in_dropout = is_in_dropout_window(timestamp_s, config_.dropouts);
        if (in_dropout) {
            return RangeMeasurement{
                timestamp_s,
                0.0,
                false
            };
        }

        const double true_range = (target_position_eci_m - truth_state.position()).norm();
        const double noise = rng_.gaussian(0.0, config_.noise_std_m);
        const double measured_range = std::max(0.0, true_range + config_.bias_m + noise);

        return RangeMeasurement{
            timestamp_s,
            measured_range,
            true
        };
    }

    // Advances sampling schedule and returns a sample if scheduled at current_t_s.
    [[nodiscard]] std::optional<RangeMeasurement> sample(
        double current_t_s,
        const spacecraft::SpacecraftState& truth_state,
        const math::Vector3& target_position_eci_m) {
        if (!schedule_.should_sample(current_t_s)) {
            return std::nullopt;
        }

        const double sample_time_s = schedule_.next_sample_time_s();
        schedule_.advance(current_t_s);

        return measure(sample_time_s, truth_state, target_position_eci_m);
    }

private:
    static void validate_config(const RangeSensorConfig& cfg) {
        if (!std::isfinite(cfg.sample_period_s) || cfg.sample_period_s <= 0.0) {
            throw std::domain_error("Range sensor sample period must be positive and finite");
        }
        if (!std::isfinite(cfg.bias_m)) {
            throw std::domain_error("Range sensor bias must be finite");
        }
        if (!std::isfinite(cfg.noise_std_m) || cfg.noise_std_m < 0.0) {
            throw std::domain_error("Range sensor noise standard deviation must be finite and non-negative");
        }
    }

    RangeSensorConfig config_;
    DeterministicRng rng_;
    SensorSchedule schedule_;
};

}  // namespace astradock::sensors
