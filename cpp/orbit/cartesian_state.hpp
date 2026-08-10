#pragma once

#include "math/vector3.hpp"

#include <cmath>

namespace astradock::orbit {

// Six-component translational state in one idealized Earth-centered inertial
// Cartesian frame. For a physical state, position is in metres and velocity is
// in metres per second.
//
// The ODE integrators also use this algebraic shape for a state derivative. In
// that role, the position slot stores dr/dt in m/s and the velocity slot stores
// dv/dt in m/s^2. The shared shape keeps Euler/RK4 arithmetic small; callers
// must preserve the documented derivative semantics.
struct CartesianState {
    math::Vector3 position;
    math::Vector3 velocity;

    [[nodiscard]] constexpr CartesianState operator+(const CartesianState& other) const noexcept {
        return {position + other.position, velocity + other.velocity};
    }

    [[nodiscard]] constexpr CartesianState operator*(double scalar) const noexcept {
        return {position * scalar, velocity * scalar};
    }
};

[[nodiscard]] constexpr CartesianState operator*(
    double scalar,
    const CartesianState& state) noexcept {
    return state * scalar;
}

[[nodiscard]] inline bool is_finite(const CartesianState& state) noexcept {
    return std::isfinite(state.position.x()) && std::isfinite(state.position.y())
        && std::isfinite(state.position.z()) && std::isfinite(state.velocity.x())
        && std::isfinite(state.velocity.y()) && std::isfinite(state.velocity.z());
}

// Customization found through argument-dependent lookup by the generic M03
// integrators, allowing them to validate the M04 state without depending on it.
[[nodiscard]] inline bool is_finite_state(const CartesianState& state) noexcept {
    return is_finite(state);
}

}  // namespace astradock::orbit
