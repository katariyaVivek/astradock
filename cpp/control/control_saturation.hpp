#pragma once

// AstraDock M16E — Controller output saturation discipline.
//
// Physical problem:
//   M16A–M16D compute DESIRED force/torque. M14 owns the hardware limits. This
//   header is the explicit boundary between them: it pairs every desired vector
//   with the achieved vector after limits and raises a flag when they differ.
//   Controllers never clip internally (their math stays linear and analyzable);
//   the assembly never clips silently (every limit sets a flag). Placing the
//   accounting here keeps both contracts intact and gives M17 the
//   desired/achieved/saturated telemetry the master spec requires.
//
// Frames/units: frame-agnostic (BODY torque or LVLH force — the caller labels);
//   units follow the input (N*m or N). Pure function on Vector3 triples.

#include "math/vector3.hpp"

#include <cmath>
#include <stdexcept>

namespace astradock::control {

// Desired-vs-achieved control output with per-axis clamp accounting.
struct SaturatedControl {
    math::Vector3 desired{};
    math::Vector3 achieved{};
    bool saturated{false};
    math::Vector3 deficit{};  // desired - achieved (zero when unsaturated)
};

[[nodiscard]] inline SaturatedControl saturate_control_vector(
    const math::Vector3& desired,
    const math::Vector3& limit) {
    if (!math::is_finite(desired) || !math::is_finite(limit)) {
        throw std::domain_error("Saturation inputs must contain only finite values");
    }
    if (limit.x() <= 0.0 || limit.y() <= 0.0 || limit.z() <= 0.0) {
        throw std::domain_error("Saturation limits must be strictly positive per axis");
    }
    const auto clamp_axis = [](double value, double max_abs) {
        if (value > max_abs) {
            return max_abs;
        }
        if (value < -max_abs) {
            return -max_abs;
        }
        return value;
    };
    const math::Vector3 achieved{
        clamp_axis(desired.x(), limit.x()),
        clamp_axis(desired.y(), limit.y()),
        clamp_axis(desired.z(), limit.z()),
    };
    const math::Vector3 deficit = desired - achieved;
    const bool saturated = deficit.squared_norm() > 0.0;
    return {desired, achieved, saturated, deficit};
}

}  // namespace astradock::control
