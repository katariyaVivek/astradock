#pragma once

#include "math/vector3.hpp"

#include <cmath>
#include <functional>
#include <stdexcept>

namespace astradock::numerics {

namespace detail {

[[nodiscard]] inline bool is_finite_state(double value) noexcept {
    return std::isfinite(value);
}

[[nodiscard]] inline bool is_finite_state(const math::Vector3& value) noexcept {
    return std::isfinite(value.x()) && std::isfinite(value.y()) && std::isfinite(value.z());
}

// A custom state may define arithmetic without exposing a general finite-value
// query. Such a state remains usable, but its own type or derivative function
// must validate its components. The integrator directly validates the state
// types used in M03: double and Vector3.
template <typename State>
[[nodiscard]] constexpr bool is_finite_state(const State&) noexcept {
    return true;
}

inline void require_finite_time(double value, const char* message) {
    if (!std::isfinite(value)) {
        throw std::domain_error(message);
    }
}

template <typename State>
void require_finite_input_state(const State& state) {
    if (!is_finite_state(state)) {
        throw std::domain_error("Integrator state must contain only finite values");
    }
}

template <typename State>
void require_finite_derivative(const State& derivative) {
    if (!is_finite_state(derivative)) {
        throw std::domain_error("Derivative function returned a non-finite value");
    }
}

template <typename State>
void require_finite_calculation(const State& state) {
    if (!is_finite_state(state)) {
        throw std::overflow_error("Integrator calculation produced a non-finite state");
    }
}

[[nodiscard]] inline double checked_step_time(double t, double offset) {
    const double step_time = t + offset;
    if (!std::isfinite(step_time)) {
        throw std::overflow_error("Integrator step time is not representable");
    }
    return step_time;
}

}  // namespace detail

// Advances x' = f(t, x) by one fixed Forward Euler step:
//     x_next = x + dt * f(t, x)
//
// t and dt must be finite, and their sum must be representable. A zero dt is
// an identity operation and does not evaluate derivative. Negative dt is
// supported for backward integration. For double and Vector3 states, supplied
// and calculated values must be finite.
template <typename State, typename Derivative>
[[nodiscard]] State euler_step(
    double t,
    const State& state,
    double dt,
    Derivative&& derivative) {
    detail::require_finite_time(t, "Integrator time must be finite");
    detail::require_finite_time(dt, "Integrator time step must be finite");
    detail::require_finite_input_state(state);
    static_cast<void>(detail::checked_step_time(t, dt));

    if (dt == 0.0) {
        return state;
    }

    const auto slope = std::invoke(derivative, t, state);
    detail::require_finite_derivative(slope);

    const auto next_state = state + slope * dt;
    detail::require_finite_calculation(next_state);
    return next_state;
}

// Advances x' = f(t, x) by one fixed classical fourth-order Runge--Kutta
// step. The four slopes sample the start, two estimated midpoints, and the end
// of the interval before being combined with weights 1, 2, 2, 1.
//
// The time-step and finite-value policy is the same as euler_step().
template <typename State, typename Derivative>
[[nodiscard]] State rk4_step(
    double t,
    const State& state,
    double dt,
    Derivative&& derivative) {
    detail::require_finite_time(t, "Integrator time must be finite");
    detail::require_finite_time(dt, "Integrator time step must be finite");
    detail::require_finite_input_state(state);

    const double half_dt = dt * 0.5;
    const double midpoint_time = detail::checked_step_time(t, half_dt);
    const double next_time = detail::checked_step_time(t, dt);

    if (dt == 0.0) {
        return state;
    }

    const auto k1 = std::invoke(derivative, t, state);
    detail::require_finite_derivative(k1);

    const auto midpoint_state_1 = state + k1 * half_dt;
    detail::require_finite_calculation(midpoint_state_1);
    const auto k2 = std::invoke(derivative, midpoint_time, midpoint_state_1);
    detail::require_finite_derivative(k2);

    const auto midpoint_state_2 = state + k2 * half_dt;
    detail::require_finite_calculation(midpoint_state_2);
    const auto k3 = std::invoke(derivative, midpoint_time, midpoint_state_2);
    detail::require_finite_derivative(k3);

    const auto end_state = state + k3 * dt;
    detail::require_finite_calculation(end_state);
    const auto k4 = std::invoke(derivative, next_time, end_state);
    detail::require_finite_derivative(k4);

    const auto weighted_slope = k1 + k2 * 2.0 + k3 * 2.0 + k4;
    detail::require_finite_calculation(weighted_slope);
    const auto next_state = state + weighted_slope * (dt / 6.0);
    detail::require_finite_calculation(next_state);
    return next_state;
}

}  // namespace astradock::numerics
