#pragma once

#include "math/constants.hpp"

#include <cmath>

namespace astradock::math {

// Normalizes an angle in radians to the half-open interval [0, 2pi).
// Handles arbitrary positive and negative values, and avoids producing -0.0.
[[nodiscard]] inline double normalize_angle_2pi_rad(double angle_rad) noexcept {
    if (!std::isfinite(angle_rad)) {
        return angle_rad;
    }

    constexpr double two_pi = 2.0 * constants::pi;
    double normalized = std::fmod(angle_rad, two_pi);
    if (normalized < 0.0) {
        normalized += two_pi;
    }

    // Handle floating-point edge case where normalized is exactly 2pi or -0.0
    if (normalized >= two_pi || std::abs(normalized) < 1.0e-15) {
        normalized = 0.0;
    }

    return normalized;
}

// Normalizes an angle in radians to the symmetric interval [-pi, pi).
[[nodiscard]] inline double normalize_angle_pi_rad(double angle_rad) noexcept {
    if (!std::isfinite(angle_rad)) {
        return angle_rad;
    }

    constexpr double two_pi = 2.0 * constants::pi;
    double normalized = normalize_angle_2pi_rad(angle_rad);
    if (normalized >= constants::pi) {
        normalized -= two_pi;
    }
    return normalized;
}

// Computes the shortest angular distance between two angles in radians,
// returning a value in [0, pi].
[[nodiscard]] inline double angular_distance_rad(double angle_a_rad, double angle_b_rad) noexcept {
    const double diff = normalize_angle_pi_rad(angle_a_rad - angle_b_rad);
    return std::abs(diff);
}

}  // namespace astradock::math
