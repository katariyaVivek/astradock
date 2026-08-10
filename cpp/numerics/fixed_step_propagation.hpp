#pragma once

#include "numerics/integrators.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace astradock::numerics {

enum class IntegrationMethod {
    forward_euler,
    classical_rk4,
};

[[nodiscard]] inline const char* integration_method_name(IntegrationMethod method) {
    switch (method) {
    case IntegrationMethod::forward_euler:
        return "euler";
    case IntegrationMethod::classical_rk4:
        return "rk4";
    }
    throw std::invalid_argument("Unknown fixed-step integration method");
}

template <typename State>
struct StateSample {
    double time_s;
    State state;
};

namespace detail {

inline void validate_propagation_times(
    double start_time_s,
    double end_time_s,
    double nominal_dt_s,
    IntegrationMethod method) {
    if (!std::isfinite(start_time_s) || !std::isfinite(end_time_s)) {
        throw std::domain_error("Propagation start and end times must be finite");
    }
    if (!std::isfinite(nominal_dt_s) || nominal_dt_s == 0.0) {
        throw std::domain_error("Propagation time step must be finite and nonzero");
    }

    static_cast<void>(integration_method_name(method));

    if (end_time_s > start_time_s && nominal_dt_s < 0.0) {
        throw std::domain_error("Forward propagation requires a positive time step");
    }
    if (end_time_s < start_time_s && nominal_dt_s > 0.0) {
        throw std::domain_error("Backward propagation requires a negative time step");
    }
}

template <typename State, typename Derivative>
[[nodiscard]] State take_selected_step(
    IntegrationMethod method,
    double t,
    const State& state,
    double dt,
    Derivative& derivative) {
    switch (method) {
    case IntegrationMethod::forward_euler:
        return euler_step(t, state, dt, derivative);
    case IntegrationMethod::classical_rk4:
        return rk4_step(t, state, dt, derivative);
    }
    throw std::invalid_argument("Unknown fixed-step integration method");
}

[[nodiscard]] inline double endpoint_tolerance_s(
    double current_time_s,
    double end_time_s,
    double nominal_dt_s) noexcept {
    const double scale = std::max(
        {1.0, std::abs(current_time_s), std::abs(end_time_s), std::abs(nominal_dt_s)});
    return 16.0 * std::numeric_limits<double>::epsilon() * scale;
}

}  // namespace detail

// Repeatedly applies an existing M03 one-step integrator and returns every
// state sample, including the initial state and the requested endpoint.
//
// nominal_dt_s is signed: positive for forward propagation and negative for
// backward propagation. Full steps are used until the endpoint is within one
// nominal step, then one shortened step lands exactly on end_time_s without
// overshoot. Equal start/end times return only the initial sample.
template <typename State, typename Derivative>
[[nodiscard]] std::vector<StateSample<State>> propagate_fixed_step(
    double start_time_s,
    double end_time_s,
    double nominal_dt_s,
    const State& initial_state,
    IntegrationMethod method,
    Derivative&& derivative) {
    detail::validate_propagation_times(start_time_s, end_time_s, nominal_dt_s, method);

    std::vector<StateSample<State>> samples;
    samples.push_back({start_time_s, initial_state});

    if (start_time_s == end_time_s) {
        return samples;
    }

    double current_time_s = start_time_s;
    State current_state = initial_state;

    while (current_time_s != end_time_s) {
        const double remaining_time_s = end_time_s - current_time_s;
        const double tolerance_s =
            detail::endpoint_tolerance_s(current_time_s, end_time_s, nominal_dt_s);
        const bool final_step =
            std::abs(remaining_time_s) <= std::abs(nominal_dt_s) + tolerance_s;
        const double step_dt_s = final_step ? remaining_time_s : nominal_dt_s;
        const double next_time_s = final_step ? end_time_s : current_time_s + step_dt_s;

        if (!std::isfinite(next_time_s) || next_time_s == current_time_s) {
            throw std::overflow_error("Propagation time step cannot advance representable time");
        }

        current_state = detail::take_selected_step(
            method,
            current_time_s,
            current_state,
            step_dt_s,
            derivative);
        current_time_s = next_time_s;
        samples.push_back({current_time_s, current_state});
    }

    return samples;
}

}  // namespace astradock::numerics
