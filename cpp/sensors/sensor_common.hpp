#pragma once

#include "math/vector3.hpp"

#include <cmath>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>

namespace astradock::sensors {

// Deterministic pseudo-random number generator for stochastic sensor simulation.
// Wraps std::mt19937_64 to ensure reproducible Gaussian noise generation across
// runs with identical random seeds.
class DeterministicRng {
public:
    explicit DeterministicRng(std::uint64_t seed = 0) noexcept
        : rng_(seed), initial_seed_(seed) {}

    // Reseeds the generator with a specific seed value.
    void reset(std::uint64_t seed) noexcept {
        initial_seed_ = seed;
        rng_.seed(seed);
    }

    // Reseeds back to the initial configured seed value.
    void reset() noexcept {
        rng_.seed(initial_seed_);
    }

    [[nodiscard]] std::uint64_t initial_seed() const noexcept {
        return initial_seed_;
    }

    // Draws a scalar sample from a 1D normal (Gaussian) distribution: X ~ N(mean, std_dev^2).
    [[nodiscard]] double gaussian(double mean = 0.0, double std_dev = 1.0) {
        if (!std::isfinite(mean) || !std::isfinite(std_dev) || std_dev < 0.0) {
            throw std::domain_error("Gaussian parameters must be finite and std_dev >= 0");
        }
        if (std_dev == 0.0) {
            return mean;
        }
        std::normal_distribution<double> dist(mean, std_dev);
        return dist(rng_);
    }

    // Draws a 3D Vector3 with isotropic zero-mean Gaussian components: v_i ~ N(0, std_dev^2).
    [[nodiscard]] math::Vector3 gaussian_vector3(double std_dev) {
        if (std_dev == 0.0) {
            return math::Vector3{};
        }
        return math::Vector3{
            gaussian(0.0, std_dev),
            gaussian(0.0, std_dev),
            gaussian(0.0, std_dev)
        };
    }

    // Draws a 3D Vector3 with independent per-axis Gaussian components: v_i ~ N(0, std_dev_i^2).
    [[nodiscard]] math::Vector3 gaussian_vector3(const math::Vector3& std_dev) {
        return math::Vector3{
            gaussian(0.0, std_dev.x()),
            gaussian(0.0, std_dev.y()),
            gaussian(0.0, std_dev.z())
        };
    }

    // Draws a uniform random direction unit vector on the 2-sphere S^2.
    [[nodiscard]] math::Vector3 uniform_unit_vector3() {
        const math::Vector3 raw = gaussian_vector3(1.0);
        const double n = raw.norm();
        if (n == 0.0) {
            return math::Vector3{1.0, 0.0, 0.0};
        }
        return raw / n;
    }

private:
    std::mt19937_64 rng_;
    std::uint64_t initial_seed_{0};
};

// Common metadata envelope returned with every discrete sensor measurement.
struct SensorSampleMetadata {
    double timestamp_s{0.0};
    bool valid{true};
};

// Time interval defining a deterministic sensor dropout / failure window.
// During an active dropout, the sensor produces samples flagged with valid = false.
struct DropoutWindow {
    double start_time_s{0.0};
    double end_time_s{0.0};

    constexpr DropoutWindow() noexcept = default;
    constexpr DropoutWindow(double start_s, double end_s) noexcept
        : start_time_s(start_s), end_time_s(end_s) {}

    [[nodiscard]] constexpr bool is_active(double t_s) const noexcept {
        return t_s >= start_time_s && t_s <= end_time_s;
    }
};

// Evaluates whether a given timestamp falls within any configured dropout window.
[[nodiscard]] inline bool is_in_dropout_window(
    double t_s,
    const std::vector<DropoutWindow>& windows) noexcept {
    for (const auto& w : windows) {
        if (w.is_active(t_s)) {
            return true;
        }
    }
    return false;
}

// Deterministic sampling schedule tracker.
// Manages sensor update frequency and next scheduled sample time without floating-point modulo traps.
class SensorSchedule {
public:
    explicit SensorSchedule(
        double sample_period_s = 1.0,
        double start_time_s = 0.0)
        : sample_period_s_(sample_period_s),
          next_sample_time_s_(start_time_s) {
        if (!std::isfinite(sample_period_s) || sample_period_s <= 0.0) {
            throw std::domain_error("Sample period must be positive and finite");
        }
        if (!std::isfinite(start_time_s)) {
            throw std::domain_error("Start time must be finite");
        }
    }

    [[nodiscard]] constexpr double sample_period_s() const noexcept {
        return sample_period_s_;
    }

    [[nodiscard]] constexpr double next_sample_time_s() const noexcept {
        return next_sample_time_s_;
    }

    // Evaluates whether a measurement should be generated at current_t_s.
    [[nodiscard]] constexpr bool should_sample(
        double current_t_s,
        double tolerance_s = 1.0e-9) const noexcept {
        return current_t_s >= next_sample_time_s_ - tolerance_s;
    }

    // Advances the scheduled next sample time past current_t_s.
    void advance(double current_t_s, double tolerance_s = 1.0e-9) noexcept {
        while (next_sample_time_s_ <= current_t_s + tolerance_s) {
            next_sample_time_s_ += sample_period_s_;
        }
    }

    // Resets schedule to a new starting epoch.
    void reset(double start_time_s = 0.0) noexcept {
        next_sample_time_s_ = start_time_s;
    }

private:
    double sample_period_s_{1.0};
    double next_sample_time_s_{0.0};
};

}  // namespace astradock::sensors
