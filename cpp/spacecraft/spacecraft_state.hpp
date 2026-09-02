#pragma once

#include "attitude/rotational_state.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"
#include "orbit/cartesian_state.hpp"

#include <cmath>

namespace astradock::spacecraft {

// Thirteen-component composite rigid-spacecraft simulation state.
//
// Represents the complete six-degree-of-freedom (6-DOF) physical state of a rigid spacecraft:
//
//   1. Translational State (3 physical DOF, 6 numerical components):
//      - position: center-of-mass position in Earth-Centered Inertial (ECI) coordinates (m).
//      - velocity: center-of-mass velocity in Earth-Centered Inertial (ECI) coordinates (m/s).
//
//   2. Rotational State (3 physical DOF, 7 numerical components):
//      - orientation: scalar-first unit quaternion q = [w, x, y, z] mapping spacecraft body
//                     frame coordinates into ECI coordinates (dimensionless, subject to ||q|| = 1).
//      - angular_velocity_rad_per_s: angular velocity vector of the spacecraft body frame
//                                     relative to the ECI frame, expressed in the body frame (rad/s).
//
// Physical DOF vs Numerical Representation:
//   Although this structure contains 13 numerical floating-point values (3 + 3 + 4 + 3),
//   it models exactly 6 physical degrees of freedom (3 translational + 3 rotational) because
//   the attitude quaternion carries an algebraic unit-norm constraint (||q|| = 1).
//
// The ODE integrators also use this algebraic shape for the composite state derivative:
//   - translational.position holds dr/dt (m/s).
//   - translational.velocity holds dv/dt (m/s^2).
//   - rotational.orientation holds dq/dt (1/s).
//   - rotational.angular_velocity_rad_per_s holds d(omega)/dt (rad/s^2).
struct SpacecraftState {
    orbit::CartesianState translational{};
    attitude::RotationalState rotational{};

    // Convenience accessors
    [[nodiscard]] constexpr const math::Vector3& position() const noexcept {
        return translational.position;
    }

    [[nodiscard]] constexpr const math::Vector3& velocity() const noexcept {
        return translational.velocity;
    }

    [[nodiscard]] constexpr const math::Quaternion& orientation() const noexcept {
        return rotational.orientation;
    }

    [[nodiscard]] constexpr const math::Vector3& angular_velocity_rad_per_s() const noexcept {
        return rotational.angular_velocity_rad_per_s;
    }

    // Linear vector-space addition for Runge--Kutta / Euler stage combination.
    [[nodiscard]] constexpr SpacecraftState operator+(const SpacecraftState& other) const noexcept {
        return {
            translational + other.translational,
            rotational + other.rotational
        };
    }

    // Scalar scaling for weighted slope combinations.
    [[nodiscard]] constexpr SpacecraftState operator*(double scalar) const noexcept {
        return {
            translational * scalar,
            rotational * scalar
        };
    }

    // Returns a copy of the spacecraft state with the attitude quaternion reprojected
    // onto the unit hypersphere S^3 to prevent numerical integration norm drift.
    [[nodiscard]] SpacecraftState normalized(double min_norm_tolerance = 1.0e-14) const {
        return {
            translational,
            rotational.normalized(min_norm_tolerance)
        };
    }
};

[[nodiscard]] constexpr SpacecraftState operator*(
    double scalar,
    const SpacecraftState& state) noexcept {
    return state * scalar;
}

[[nodiscard]] inline bool is_finite(const SpacecraftState& state) noexcept {
    return orbit::is_finite(state.translational) && attitude::is_finite(state.rotational);
}

// Customization point found via Argument-Dependent Lookup (ADL) by generic M03 integrators.
[[nodiscard]] inline bool is_finite_state(const SpacecraftState& state) noexcept {
    return is_finite(state);
}

}  // namespace astradock::spacecraft
