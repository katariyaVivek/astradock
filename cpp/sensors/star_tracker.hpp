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

// Configuration parameters for the Star Tracker attitude sensor simulation.
struct StarTrackerConfig {
    // 1-sigma isotropic angular attitude error (radians)
    double noise_std_rad{0.0};

    // Constant boresight mounting bias (axis and rotation angle in radians)
    math::Vector3 bias_axis{0.0, 0.0, 1.0};
    double bias_angle_rad{0.0};

    // Timing & sampling parameters (nominal 10 Hz = 0.1 s)
    double sample_period_s{0.1};
    std::vector<DropoutWindow> dropouts{};
    std::uint64_t random_seed{0};

    constexpr StarTrackerConfig() noexcept = default;
};

// Star Tracker measurement packet.
// Contains measured unit quaternion mapping spacecraft Body frame into ECI frame (q_ECI_Body).
struct StarTrackerMeasurement {
    double timestamp_s{0.0};
    math::Quaternion orientation_eci_from_body{math::Quaternion::identity()};
    bool valid{true};
};

// Star Tracker attitude sensor simulation.
//
// Models an optical star tracker measuring spacecraft 3D orientation relative to
// inertial star catalogues, outputting a scalar-first unit quaternion q (ECI from Body).
//
// Physical Perturbation Model:
//   Quaternion errors are modeled strictly as physical 3D rotations in SO(3),
//   rather than component-wise additions:
//
//     1. Constant Bias Rotation: q_bias = from_axis_angle(u_bias, theta_bias)
//     2. Stochastic Angular Noise: delta_q = from_axis_angle(u_rand, theta_noise)
//        where u_rand is uniformly distributed on S^2 and theta_noise ~ N(0, sigma^2).
//     3. Composite Physical Attitude Measurement:
//        q_meas = (delta_q ⊗ q_bias ⊗ q_true).normalized()
//
// Properties:
//   - Guaranteed to produce valid unit quaternions (||q_meas|| = 1).
//   - Complies with M08/M10 quaternion conventions (scalar-first, active rotation).
//   - Preserves double-cover equivalence (q and -q represent identical physical attitudes).
class StarTrackerSensor {
public:
    explicit StarTrackerSensor(const StarTrackerConfig& config)
        : config_(config),
          rng_(config.random_seed),
          schedule_(config.sample_period_s, 0.0),
          bias_quaternion_(compute_bias_quaternion(config)) {
        validate_config(config);
    }

    [[nodiscard]] const StarTrackerConfig& config() const noexcept {
        return config_;
    }

    // Resets random number generator and sampling schedule to initial epoch.
    void reset(double start_time_s = 0.0) noexcept {
        rng_.reset();
        schedule_.reset(start_time_s);
    }

    // Evaluates instantaneous measurement without modifying sampling schedule.
    [[nodiscard]] StarTrackerMeasurement measure(
        double timestamp_s,
        const spacecraft::SpacecraftState& truth_state) {
        if (!spacecraft::is_finite(truth_state)) {
            throw std::domain_error("Star tracker truth state must be finite");
        }

        const bool in_dropout = is_in_dropout_window(timestamp_s, config_.dropouts);
        if (in_dropout) {
            return StarTrackerMeasurement{
                timestamp_s,
                math::Quaternion::identity(),
                false
            };
        }

        const math::Quaternion true_q = truth_state.orientation().normalized();

        // Generate small physical rotation perturbation
        math::Quaternion noise_q = math::Quaternion::identity();
        if (config_.noise_std_rad > 0.0) {
            const double theta_err = rng_.gaussian(0.0, config_.noise_std_rad);
            const math::Vector3 u_axis = rng_.uniform_unit_vector3();
            if (std::abs(theta_err) > 0.0) {
                noise_q = math::Quaternion::from_axis_angle(u_axis, theta_err);
            }
        }

        // Apply rotation perturbation: q_meas = (noise_q ⊗ bias_q ⊗ true_q).normalized()
        const math::Quaternion perturbed_q = (noise_q * bias_quaternion_ * true_q).normalized();

        return StarTrackerMeasurement{
            timestamp_s,
            perturbed_q,
            true
        };
    }

    // Advances sampling schedule and returns a sample if scheduled at current_t_s.
    [[nodiscard]] std::optional<StarTrackerMeasurement> sample(
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
    [[nodiscard]] static math::Quaternion compute_bias_quaternion(const StarTrackerConfig& cfg) {
        if (cfg.bias_angle_rad == 0.0 || cfg.bias_axis.norm() == 0.0) {
            return math::Quaternion::identity();
        }
        return math::Quaternion::from_axis_angle(cfg.bias_axis, cfg.bias_angle_rad);
    }

    static void validate_config(const StarTrackerConfig& cfg) {
        if (!std::isfinite(cfg.sample_period_s) || cfg.sample_period_s <= 0.0) {
            throw std::domain_error("Star tracker sample period must be positive and finite");
        }
        if (!std::isfinite(cfg.noise_std_rad) || cfg.noise_std_rad < 0.0) {
            throw std::domain_error("Star tracker noise standard deviation must be finite and non-negative");
        }
        if (!std::isfinite(cfg.bias_angle_rad)) {
            throw std::domain_error("Star tracker bias angle must be finite");
        }
        if (!math::is_finite(cfg.bias_axis)) {
            throw std::domain_error("Star tracker bias axis must be finite");
        }
    }

    StarTrackerConfig config_;
    DeterministicRng rng_;
    SensorSchedule schedule_;
    math::Quaternion bias_quaternion_{math::Quaternion::identity()};
};

}  // namespace astradock::sensors
