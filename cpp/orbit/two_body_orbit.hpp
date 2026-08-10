#pragma once

#include "dynamics/two_body.hpp"
#include "math/constants.hpp"
#include "orbit/cartesian_state.hpp"

#include <cmath>
#include <stdexcept>

namespace astradock::orbit {

namespace detail {

inline void require_positive_finite(double value, const char* message) {
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::domain_error(message);
    }
}

inline void require_valid_orbital_state(const CartesianState& state) {
    if (!is_finite(state)) {
        throw std::domain_error("Cartesian orbital state must contain only finite values");
    }
    if (state.position.norm() == 0.0) {
        throw std::domain_error("Cartesian orbital state is undefined at Earth's center");
    }
}

}  // namespace detail

// Returns the first-order derivative of a Cartesian two-body orbital state.
// The time argument is intentionally unnamed because point-mass gravity is
// autonomous: acceleration depends on state, not absolute time.
//
// Input state:
//   position: idealized Earth-centered inertial position, m
//   velocity: idealized Earth-centered inertial velocity, m/s
// Returned state-shaped derivative:
//   position slot: dr/dt = velocity, m/s
//   velocity slot: dv/dt = two-body acceleration, m/s^2
[[nodiscard]] inline CartesianState two_body_state_derivative(
    double,
    const CartesianState& state,
    double gravitational_parameter_m3_per_s2) {
    if (!is_finite(state)) {
        throw std::domain_error("Cartesian orbital state must contain only finite values");
    }

    return {
        state.velocity,
        dynamics::two_body_acceleration(
            state.position,
            gravitational_parameter_m3_per_s2),
    };
}

// Specific mechanical energy, epsilon = |v|^2/2 - mu/|r|, in m^2/s^2
// (equivalently J/kg). It is conserved by the ideal continuous two-body model.
[[nodiscard]] inline double specific_orbital_energy_m2_per_s2(
    const CartesianState& state,
    double gravitational_parameter_m3_per_s2) {
    detail::require_valid_orbital_state(state);
    detail::require_positive_finite(
        gravitational_parameter_m3_per_s2,
        "Orbital gravitational parameter must be finite and positive");

    const double radius_m = state.position.norm();
    const double speed_m_per_s = state.velocity.norm();
    const double energy_m2_per_s2 =
        0.5 * speed_m_per_s * speed_m_per_s
        - gravitational_parameter_m3_per_s2 / radius_m;
    if (!std::isfinite(energy_m2_per_s2)) {
        throw std::overflow_error("Specific orbital energy is not representable");
    }
    return energy_m2_per_s2;
}

// Specific angular momentum, h = r x v, in m^2/s. Its direction is normal to
// the orbital plane and is conserved by ideal central gravity.
[[nodiscard]] inline math::Vector3 specific_angular_momentum_m2_per_s(
    const CartesianState& state) {
    detail::require_valid_orbital_state(state);

    const math::Vector3 angular_momentum_m2_per_s = state.position.cross(state.velocity);
    if (!std::isfinite(angular_momentum_m2_per_s.x())
        || !std::isfinite(angular_momentum_m2_per_s.y())
        || !std::isfinite(angular_momentum_m2_per_s.z())) {
        throw std::overflow_error("Specific angular momentum is not representable");
    }
    return angular_momentum_m2_per_s;
}

[[nodiscard]] inline double circular_orbit_speed_m_per_s(
    double gravitational_parameter_m3_per_s2,
    double orbital_radius_m) {
    detail::require_positive_finite(
        gravitational_parameter_m3_per_s2,
        "Circular-orbit gravitational parameter must be finite and positive");
    detail::require_positive_finite(
        orbital_radius_m,
        "Circular-orbit radius must be finite and positive");

    const double speed_m_per_s =
        std::sqrt(gravitational_parameter_m3_per_s2 / orbital_radius_m);
    if (!std::isfinite(speed_m_per_s)) {
        throw std::overflow_error("Circular-orbit speed is not representable");
    }
    return speed_m_per_s;
}

[[nodiscard]] inline double circular_orbit_period_s(
    double gravitational_parameter_m3_per_s2,
    double orbital_radius_m) {
    detail::require_positive_finite(
        gravitational_parameter_m3_per_s2,
        "Circular-orbit gravitational parameter must be finite and positive");
    detail::require_positive_finite(
        orbital_radius_m,
        "Circular-orbit radius must be finite and positive");

    // r * sqrt(r / mu) is algebraically sqrt(r^3 / mu) while avoiding an
    // unnecessary intermediate r^3 overflow for large but representable radii.
    const double period_s = 2.0 * constants::pi * orbital_radius_m
        * std::sqrt(orbital_radius_m / gravitational_parameter_m3_per_s2);
    if (!std::isfinite(period_s)) {
        throw std::overflow_error("Circular-orbit period is not representable");
    }
    return period_s;
}

}  // namespace astradock::orbit
