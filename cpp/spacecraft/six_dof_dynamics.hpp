#pragma once

#include "attitude/rigid_body.hpp"
#include "dynamics/two_body.hpp"
#include "math/constants.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"
#include "numerics/fixed_step_propagation.hpp"
#include "numerics/integrators.hpp"
#include "orbit/two_body_orbit.hpp"
#include "spacecraft/force_torque.hpp"
#include "spacecraft/spacecraft_parameters.hpp"
#include "spacecraft/spacecraft_state.hpp"

#include <algorithm>
#include <cmath>
#include <concepts>
#include <functional>
#include <stdexcept>
#include <vector>

namespace astradock::spacecraft {

namespace detail {

inline void require_valid_spacecraft_state(const SpacecraftState& state) {
    if (!is_finite(state)) {
        throw std::domain_error("Spacecraft state must contain only finite components");
    }
    if (state.position().norm() == 0.0) {
        throw std::domain_error("Spacecraft position is undefined at the central-body origin");
    }
}

inline void require_valid_parameters(const SpacecraftParameters& params) {
    if (!is_valid(params)) {
        throw std::domain_error("Spacecraft parameters must contain positive and finite values");
    }
}

inline void require_finite_force_torque(const ForceTorqueInput& input) {
    if (!is_finite(input)) {
        throw std::domain_error("Force and torque inputs must contain only finite components");
    }
}

}  // namespace detail

// Computes the shortest physical rotation angle between a reference attitude
// and an estimated attitude quaternion.
//
// Mathematical derivation:
//   Relative rotation quaternion: q_err = q_reference^(-1) ⊗ q_estimate
//   Smallest rotation angle:      theta_err = 2 * acos(|q_err.w()|)
//
// Properties:
//   - Invariant to quaternion double-cover (q and -q produce identical theta_err = 0).
//   - Handles floating-point clamping robustly for |q_err.w| near 1.0.
//   - Range: [0, pi] radians.
[[nodiscard]] inline double quaternion_orientation_error_rad(
    const math::Quaternion& q_reference,
    const math::Quaternion& q_estimate) {
    if (!math::is_finite(q_reference) || !math::is_finite(q_estimate)) {
        throw std::domain_error("Quaternions must contain only finite values for orientation error");
    }
    if (!q_reference.is_unit(1.0e-3) || !q_estimate.is_unit(1.0e-3)) {
        throw std::domain_error("Quaternions must have unit norm for orientation error calculation");
    }

    const math::Quaternion q_err = q_reference.conjugate() * q_estimate;
    const double clamped_w = std::clamp(std::abs(q_err.w()), 0.0, 1.0);
    return 2.0 * std::acos(clamped_w);
}

// Evaluates the unified first-order time derivative of the 13-component 6-DOF SpacecraftState.
//
// Governing Equations:
//   1. Translation (ECI frame):
//      dr/dt = v
//      dv/dt = a_grav(r, mu)
//   2. Rotation (Body frame):
//      dq/dt = 0.5 * q ⊗ [0, omega_B]
//      d(omega_B)/dt = I^(-1) * [ tau_B - omega_B x (I * omega_B) ]
//
// Modularity:
//   Reuses canonical M02 point-mass two-body gravity (dynamics::two_body_acceleration)
//   and canonical M09 rigid-body Euler kinematics/dynamics (attitude::rotational_state_derivative).
[[nodiscard]] inline SpacecraftState spacecraft_state_derivative(
    double t,
    const SpacecraftState& state,
    const SpacecraftParameters& params,
    const ForceTorqueInput& input = ForceTorqueInput{}) {
    detail::require_valid_spacecraft_state(state);
    detail::require_valid_parameters(params);
    detail::require_finite_force_torque(input);

    // Translational branch (evaluated in ECI)
    const math::Vector3 dr_dt = state.velocity();
    const math::Vector3 dv_dt = dynamics::two_body_acceleration(
        state.position(),
        params.gravitational_parameter_m3_per_s2
    );

    // Rotational branch (evaluated in Body frame)
    const attitude::RotationalState rot_derivative = attitude::rotational_state_derivative(
        t,
        state.rotational,
        params.inertia,
        input.torque_body_Nm
    );

    return SpacecraftState{
        orbit::CartesianState{dr_dt, dv_dt},
        rot_derivative
    };
}

// Advances the 6-DOF spacecraft state by one fixed RK4 integration step under constant inputs.
//
// Normalization Policy:
//   Intermediate RK4 substages sample unconstrained linear state space.
//   At the completion of the full step, the attitude quaternion is reprojected onto S^3
//   (if renormalize_quaternion is true) to prevent norm drift without degrading 4th-order accuracy.
[[nodiscard]] inline SpacecraftState rk4_step_spacecraft(
    double t,
    const SpacecraftState& state,
    double dt,
    const SpacecraftParameters& params,
    const ForceTorqueInput& input = ForceTorqueInput{},
    bool renormalize_quaternion = true) {
    detail::require_valid_parameters(params);

    const auto derivative_evaluator = [&](double step_t, const SpacecraftState& step_state) -> SpacecraftState {
        return spacecraft_state_derivative(step_t, step_state, params, input);
    };

    const SpacecraftState raw_next_state = numerics::rk4_step(
        t,
        state,
        dt,
        derivative_evaluator
    );

    return renormalize_quaternion ? raw_next_state.normalized() : raw_next_state;
}

// Overload for time- or state-varying force/torque callables: input_func(t, state) -> ForceTorqueInput.
template <typename ForceTorqueCallable>
requires std::invocable<ForceTorqueCallable, double, const SpacecraftState&>
[[nodiscard]] inline SpacecraftState rk4_step_spacecraft(
    double t,
    const SpacecraftState& state,
    double dt,
    const SpacecraftParameters& params,
    ForceTorqueCallable&& input_callable,
    bool renormalize_quaternion = true) {
    detail::require_valid_parameters(params);

    const auto derivative_evaluator = [&](double step_t, const SpacecraftState& step_state) -> SpacecraftState {
        const ForceTorqueInput input = std::invoke(input_callable, step_t, step_state);
        return spacecraft_state_derivative(step_t, step_state, params, input);
    };

    const SpacecraftState raw_next_state = numerics::rk4_step(
        t,
        state,
        dt,
        derivative_evaluator
    );

    return renormalize_quaternion ? raw_next_state.normalized() : raw_next_state;
}

// Propagates a 6-DOF spacecraft state across a time interval using fixed-step RK4 with
// post-step quaternion normalization.
[[nodiscard]] inline std::vector<numerics::StateSample<SpacecraftState>> propagate_spacecraft_fixed_step(
    double start_time_s,
    double end_time_s,
    double nominal_dt_s,
    const SpacecraftState& initial_state,
    const SpacecraftParameters& params,
    const ForceTorqueInput& input = ForceTorqueInput{},
    bool renormalize_quaternion = true) {
    detail::require_valid_parameters(params);

    if (!std::isfinite(start_time_s) || !std::isfinite(end_time_s) || !std::isfinite(nominal_dt_s) || nominal_dt_s == 0.0) {
        throw std::domain_error("Propagation times must be finite and nominal dt nonzero");
    }

    std::vector<numerics::StateSample<SpacecraftState>> samples;
    samples.push_back({start_time_s, initial_state});

    if (start_time_s == end_time_s) {
        return samples;
    }

    double current_t = start_time_s;
    SpacecraftState current_state = initial_state;

    while (current_t != end_time_s) {
        const double remaining_t = end_time_s - current_t;
        const double tol_s = numerics::detail::endpoint_tolerance_s(current_t, end_time_s, nominal_dt_s);
        const bool is_final_step = std::abs(remaining_t) <= std::abs(nominal_dt_s) + tol_s;
        const double step_dt = is_final_step ? remaining_t : nominal_dt_s;
        const double next_t = is_final_step ? end_time_s : current_t + step_dt;

        current_state = rk4_step_spacecraft(
            current_t,
            current_state,
            step_dt,
            params,
            input,
            renormalize_quaternion
        );
        current_t = next_t;
        samples.push_back({current_t, current_state});
    }

    return samples;
}

}  // namespace astradock::spacecraft
