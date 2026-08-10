#pragma once

#include "math/vector3.hpp"

#include <cmath>
#include <stdexcept>

namespace astradock::dynamics {

// Computes point-mass central gravitational acceleration in a
// central-body-centered inertial frame.
//
// Inputs:
//   position_central_body_inertial_m: position from the central body's center,
//       expressed in metres in one inertial frame.
//   gravitational_parameter_m3_per_s2: positive central-body mu = G*M, in
//       cubic metres per square second.
// Output:
//   acceleration in the same inertial frame, in metres per square second.
//
// Model: a = -mu * r / |r|^3. The model assumes point masses and includes no
// atmosphere, oblateness, third bodies, or other perturbations.
[[nodiscard]] inline math::Vector3 two_body_acceleration(
    const math::Vector3& position_central_body_inertial_m,
    double gravitational_parameter_m3_per_s2) {
    if (!std::isfinite(gravitational_parameter_m3_per_s2)
        || gravitational_parameter_m3_per_s2 <= 0.0) {
        throw std::domain_error("Two-body gravitational parameter must be finite and positive");
    }

    if (!std::isfinite(position_central_body_inertial_m.x())
        || !std::isfinite(position_central_body_inertial_m.y())
        || !std::isfinite(position_central_body_inertial_m.z())) {
        throw std::domain_error("Two-body position must contain only finite components");
    }

    const double radius_m = position_central_body_inertial_m.norm();
    if (radius_m == 0.0) {
        throw std::domain_error("Two-body acceleration is undefined at the central-body origin");
    }
    if (!std::isfinite(radius_m)) {
        throw std::domain_error("Two-body position magnitude must be finite");
    }

    const double inverse_radius_m = 1.0 / radius_m;
    const double acceleration_scale_per_s2 =
        -gravitational_parameter_m3_per_s2 * inverse_radius_m * inverse_radius_m
        * inverse_radius_m;

    if (!std::isfinite(acceleration_scale_per_s2)) {
        throw std::overflow_error("Two-body acceleration scale is not representable");
    }

    const math::Vector3 acceleration_central_body_inertial_m_per_s2 =
        position_central_body_inertial_m * acceleration_scale_per_s2;
    if (!std::isfinite(acceleration_central_body_inertial_m_per_s2.x())
        || !std::isfinite(acceleration_central_body_inertial_m_per_s2.y())
        || !std::isfinite(acceleration_central_body_inertial_m_per_s2.z())) {
        throw std::overflow_error("Two-body acceleration is not representable");
    }

    return acceleration_central_body_inertial_m_per_s2;
}

}  // namespace astradock::dynamics
