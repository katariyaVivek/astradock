#pragma once

#include "atmospheric_drag.hpp"
#include "gravity_gradient.hpp"
#include "j2_gravity.hpp"
#include "math/constants.hpp"
#include "math/vector3.hpp"
#include "numerics/fixed_step_propagation.hpp"
#include "numerics/integrators.hpp"
#include "spacecraft/force_torque.hpp"
#include "spacecraft/six_dof_dynamics.hpp"
#include "spacecraft/spacecraft_parameters.hpp"
#include "spacecraft/spacecraft_state.hpp"
#include "third_body_gravity.hpp"

#include <cmath>
#include <vector>

namespace astradock::environment {

// Configuration switches to enable/disable specific environmental perturbation models.
struct EnvironmentConfiguration {
    bool enable_j2{false};
    bool enable_drag{false};
    bool enable_third_body{false};
    bool enable_gravity_gradient{false};

    constexpr EnvironmentConfiguration() noexcept = default;

    constexpr EnvironmentConfiguration(
        bool j2,
        bool drag,
        bool third_body,
        bool gravity_gradient) noexcept
        : enable_j2(j2),
          enable_drag(drag),
          enable_third_body(third_body),
          enable_gravity_gradient(gravity_gradient) {}

    [[nodiscard]] constexpr bool any_enabled() const noexcept {
        return enable_j2 || enable_drag || enable_third_body || enable_gravity_gradient;
    }
};

// Physical and atmospheric parameters for environmental perturbation models.
struct EnvironmentalParameters {
    // Spacecraft aerodynamic properties
    double mass_kg{500.0};
    double drag_coefficient_cd{2.2};
    double drag_reference_area_m2{2.0};

    // Planetary constants
    double earth_reference_radius_m{constants::earth_reference_radius_m};
    double earth_j2{constants::earth_j2};
    double earth_rotation_rate_rad_per_s{constants::earth_rotation_rate_rad_per_s};

    // Atmospheric exponential model parameters
    double atmosphere_ref_alt_m{500.0e3};
    double atmosphere_ref_density_kg_per_m3{6.967e-13};
    double atmosphere_scale_height_m{63.8e3};

    // Third-body parameters (default: Moon at canonical mean orbital radius)
    math::Vector3 third_body_position_eci_m{384'400'000.0, 0.0, 0.0};
    double third_body_mu_m3_per_s2{constants::moon_gravitational_parameter_m3_per_s2};

    constexpr EnvironmentalParameters() noexcept = default;
};

// Aggregated environmental forces and torques.
struct EnvironmentalEffects {
    math::Vector3 acceleration_eci_mps2{0.0, 0.0, 0.0};  // Translational perturbation acceleration in ECI (m/s^2)
    math::Vector3 torque_body_Nm{0.0, 0.0, 0.0};         // Environmental torque in Spacecraft BODY frame (N*m)
};

// Evaluates all active environmental effects for a given spacecraft state.
[[nodiscard]] inline EnvironmentalEffects evaluate_environmental_effects(
    const spacecraft::SpacecraftState& state,
    const spacecraft::SpacecraftParameters& params,
    const EnvironmentalParameters& env_params,
    const EnvironmentConfiguration& config) {
    EnvironmentalEffects effects;

    if (!config.any_enabled()) {
        return effects;
    }

    // 1. Earth J2 Gravitational Perturbation (ECI)
    if (config.enable_j2) {
        effects.acceleration_eci_mps2 = effects.acceleration_eci_mps2 + j2_acceleration_eci(
            state.position(),
            params.gravitational_parameter_m3_per_s2,
            env_params.earth_reference_radius_m,
            env_params.earth_j2
        );
    }

    // 2. Atmospheric Drag Perturbation (ECI)
    if (config.enable_drag) {
        effects.acceleration_eci_mps2 = effects.acceleration_eci_mps2 + drag_acceleration_eci(
            state.velocity(),
            state.position(),
            env_params.mass_kg,
            env_params.drag_coefficient_cd,
            env_params.drag_reference_area_m2,
            env_params.earth_reference_radius_m,
            env_params.earth_rotation_rate_rad_per_s,
            env_params.atmosphere_ref_alt_m,
            env_params.atmosphere_ref_density_kg_per_m3,
            env_params.atmosphere_scale_height_m
        );
    }

    // 3. Third-Body Tidal Gravitational Perturbation (ECI)
    if (config.enable_third_body) {
        effects.acceleration_eci_mps2 = effects.acceleration_eci_mps2 + third_body_acceleration_eci(
            state.position(),
            env_params.third_body_position_eci_m,
            env_params.third_body_mu_m3_per_s2
        );
    }

    // 4. Gravity-Gradient Environmental Torque (BODY)
    if (config.enable_gravity_gradient) {
        // Safe check on quaternion norm for intermediate stages
        const double q_norm = state.rotational.orientation.norm();
        if (q_norm > 0.0) {
            const math::Quaternion q_unit = state.rotational.orientation.normalized();
            effects.torque_body_Nm = effects.torque_body_Nm + gravity_gradient_torque_body(
                state.position(),
                q_unit,
                params.inertia,
                params.gravitational_parameter_m3_per_s2
            );
        }
    }

    return effects;
}

// Evaluates the unified first-order time derivative of the 13-component SpacecraftState
// incorporating environmental force and torque models.
[[nodiscard]] inline spacecraft::SpacecraftState spacecraft_environmental_derivative(
    double t,
    const spacecraft::SpacecraftState& state,
    const spacecraft::SpacecraftParameters& params,
    const EnvironmentalParameters& env_params,
    const EnvironmentConfiguration& config,
    const spacecraft::ForceTorqueInput& extra_input = {}) {
    // 1. Evaluate environmental perturbations
    const EnvironmentalEffects env_effects = evaluate_environmental_effects(
        state,
        params,
        env_params,
        config
    );

    // 2. Translational kinematics and dynamics
    const math::Vector3 dr_dt = state.velocity();

    // Central two-body gravity acceleration
    const math::Vector3 a_twobody = dynamics::two_body_acceleration(
        state.position(),
        params.gravitational_parameter_m3_per_s2
    );

    // Total translational acceleration = Two-body gravity + Environmental accelerations + Non-grav forces / mass
    math::Vector3 dv_dt = a_twobody + env_effects.acceleration_eci_mps2;
    if (extra_input.force_eci_N.squared_norm() > 0.0) {
        dv_dt = dv_dt + (extra_input.force_eci_N / env_params.mass_kg);
    }

    // 3. Rotational kinematics and dynamics
    const math::Vector3 total_torque_body = extra_input.torque_body_Nm + env_effects.torque_body_Nm;
    const attitude::RotationalState rot_derivative = attitude::rotational_state_derivative(
        t,
        state.rotational,
        params.inertia,
        total_torque_body
    );

    return spacecraft::SpacecraftState{
        orbit::CartesianState{dr_dt, dv_dt},
        rot_derivative
    };
}

// Single-step classical RK4 integrator for integrated 6-DOF spacecraft with environmental models.
// Enforces post-step unit quaternion reprojection.
[[nodiscard]] inline spacecraft::SpacecraftState rk4_step_spacecraft_environmental(
    double t,
    const spacecraft::SpacecraftState& state,
    double dt,
    const spacecraft::SpacecraftParameters& params,
    const EnvironmentalParameters& env_params,
    const EnvironmentConfiguration& config,
    const spacecraft::ForceTorqueInput& input = {},
    bool renormalize_quaternion = true) {
    const auto derivative_evaluator = [&](double step_t, const spacecraft::SpacecraftState& step_state) {
        return spacecraft_environmental_derivative(step_t, step_state, params, env_params, config, input);
    };

    const spacecraft::SpacecraftState raw_next_state = numerics::rk4_step(
        t,
        state,
        dt,
        derivative_evaluator
    );

    return renormalize_quaternion ? raw_next_state.normalized() : raw_next_state;
}

// Fixed-step propagation loop for integrated 6-DOF spacecraft under environmental models.
[[nodiscard]] inline std::vector<numerics::StateSample<spacecraft::SpacecraftState>> propagate_spacecraft_environmental(
    const spacecraft::SpacecraftState& initial_state,
    double duration_s,
    double dt_s,
    const spacecraft::SpacecraftParameters& params,
    const EnvironmentalParameters& env_params,
    const EnvironmentConfiguration& config,
    const spacecraft::ForceTorqueInput& input = {},
    bool renormalize_quaternion = true) {
    if (!spacecraft::is_finite(initial_state)) {
        throw std::domain_error("Initial spacecraft state must be finite");
    }
    if (!spacecraft::is_valid(params)) {
        throw std::domain_error("Spacecraft parameters must be valid and finite");
    }
    if (!std::isfinite(duration_s) || !std::isfinite(dt_s) || dt_s <= 0.0) {
        throw std::domain_error("Propagation duration and timestep must be finite and positive");
    }

    std::vector<numerics::StateSample<spacecraft::SpacecraftState>> samples;
    samples.push_back({0.0, initial_state});

    if (duration_s == 0.0) {
        return samples;
    }

    double current_t = 0.0;
    spacecraft::SpacecraftState current_state = initial_state;

    while (current_t < duration_s) {
        const double remaining_t = duration_s - current_t;
        const double tol_s = numerics::detail::endpoint_tolerance_s(current_t, duration_s, dt_s);
        const bool is_final_step = std::abs(remaining_t) <= std::abs(dt_s) + tol_s;
        const double step_dt = is_final_step ? remaining_t : dt_s;
        const double next_t = is_final_step ? duration_s : current_t + step_dt;

        current_state = rk4_step_spacecraft_environmental(
            current_t,
            current_state,
            step_dt,
            params,
            env_params,
            config,
            input,
            renormalize_quaternion
        );
        current_t = next_t;
        samples.push_back({current_t, current_state});
    }

    return samples;
}

}  // namespace astradock::environment
