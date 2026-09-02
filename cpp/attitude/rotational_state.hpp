#pragma once

#include "math/quaternion.hpp"
#include "math/vector3.hpp"

#include <cmath>

namespace astradock::attitude {

// Seven-component rigid-body rotational state.
//
// Represents the instantaneous rotational state of a spacecraft:
//   - orientation: scalar-first unit quaternion q = [w, x, y, z] mapping
//                  body-frame coordinates into inertial-frame coordinates
//                  (dimensionless).
//   - angular_velocity_rad_per_s: angular velocity vector of the body frame
//                                relative to the inertial frame, expressed in
//                                the body frame (rad/s).
//
// In accordance with numerical ODE integration architecture, this structure
// also serves as its own algebraic state derivative (where the orientation slot
// stores dq/dt and the angular velocity slot stores d(omega)/dt in rad/s^2).
struct RotationalState {
    math::Quaternion orientation{math::Quaternion::identity()};
    math::Vector3 angular_velocity_rad_per_s{};

    // Linear vector-space addition for Runge--Kutta / Euler stage combination.
    [[nodiscard]] constexpr RotationalState operator+(const RotationalState& other) const noexcept {
        return {
            orientation + other.orientation,
            angular_velocity_rad_per_s + other.angular_velocity_rad_per_s
        };
    }

    // Scalar scaling for weighted slope combinations.
    [[nodiscard]] constexpr RotationalState operator*(double scalar) const noexcept {
        return {
            orientation * scalar,
            angular_velocity_rad_per_s * scalar
        };
    }

    // Returns a copy of the rotational state with the quaternion normalized to unit length.
    // Used at discrete integration step boundaries to prevent numerical norm drift.
    [[nodiscard]] RotationalState normalized(double min_norm_tolerance = 1.0e-14) const {
        return {
            orientation.normalized(min_norm_tolerance),
            angular_velocity_rad_per_s
        };
    }
};

[[nodiscard]] constexpr RotationalState operator*(
    double scalar,
    const RotationalState& state) noexcept {
    return state * scalar;
}

[[nodiscard]] inline bool is_finite(const RotationalState& state) noexcept {
    return math::is_finite(state.orientation) && math::is_finite(state.angular_velocity_rad_per_s);
}

// Customization point found via Argument-Dependent Lookup (ADL) by the generic
// M03 integrators (euler_step, rk4_step).
[[nodiscard]] inline bool is_finite_state(const RotationalState& state) noexcept {
    return is_finite(state);
}

}  // namespace astradock::attitude
