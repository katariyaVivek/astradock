#pragma once

// AstraDock M14 — Spacecraft assembly: bus mass + reaction-wheel array + thrusters.
//
// Physical problem:
//   M14A models one wheel and M14C one thruster. A spacecraft carries an array of
//   both: N wheels on (ideally orthogonal) spin axes provide 3-axis torque without
//   propellant, and M thrusters provide force plus offset torque at the cost of
//   mass. This header composes them into a single deterministic actuator assembly
//   whose outputs plug directly into the 6-DOF truth dynamics:
//
//     sum over wheels:    tau_wheels_body = -sum(applied_i * axis_i)
//     sum over thrusters: F_body = sum(F_i), tau_thr_body = sum(r_i x F_i)
//     mass:               m -= sum(propellant_i)
//
// Frames: all wheel axes, thruster mounts/directions, summed force/torque in BODY;
//   summed force additionally provided in ECI via the attitude quaternion.
// Units: SI throughout (N, N*m, kg, s, rad/s, N*m*s).
//
// Momentum bookkeeping (M14B result):
//   total_wheel_momentum_body_Nms = sum(I_w,i * omega_w,i * axis_i).
//   The bus-plus-wheels system conserves total angular momentum under internal
//   wheel torques: any change in wheel momentum appears with opposite sign in the
//   bus reaction torque. External thruster/environment torques change the total.

#include "actuators/reaction_wheel.hpp"
#include "actuators/thruster.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace astradock::actuators {

// Wet/dry mass state. Propellant mass = wet - dry; thruster steps clamp at dry.
struct SpacecraftMassState {
    double wet_mass_kg{500.0};
    double dry_mass_kg{400.0};
};

[[nodiscard]] inline bool is_valid_mass(const SpacecraftMassState& mass) noexcept {
    return std::isfinite(mass.wet_mass_kg) && std::isfinite(mass.dry_mass_kg)
        && mass.wet_mass_kg > 0.0 && mass.dry_mass_kg > 0.0
        && mass.dry_mass_kg <= mass.wet_mass_kg;
}

// Complete actuator command for one assembly step.
struct ActuatorCommand {
    std::vector<double> wheel_torques_Nm{};  // signed motor torque per wheel
    std::vector<double> thruster_forces_N{};  // thrust magnitude per thruster (>= 0)
};

// Aggregated actuator output for one assembly step. Commanded vs achieved are
// both reported: any limiting (torque/speed/propellant) sets the matching flag.
struct ActuatorOutput {
    math::Vector3 total_force_body_N{};
    math::Vector3 total_force_eci_N{};
    math::Vector3 total_torque_body_Nm{};
    math::Vector3 wheel_reaction_torque_body_Nm{};
    math::Vector3 thruster_torque_body_Nm{};
    double total_wheel_momentum_body_x_Nms{0.0};
    double total_wheel_momentum_body_y_Nms{0.0};
    double total_wheel_momentum_body_z_Nms{0.0};
    double propellant_used_kg{0.0};
    bool any_torque_saturated{false};
    bool any_speed_saturated{false};
    bool any_authority_lost{false};
    bool any_thrust_saturated{false};
    bool any_propellant_depleted{false};
};

[[nodiscard]] inline math::Vector3 total_wheel_momentum_body_Nms(
    const std::vector<ReactionWheelParameters>& wheel_params,
    const std::vector<ReactionWheelState>& wheel_states) {
    if (wheel_params.size() != wheel_states.size()) {
        throw std::domain_error("Wheel parameter and state counts must match");
    }
    math::Vector3 total{};
    for (std::size_t i = 0; i < wheel_params.size(); ++i) {
        if (!is_valid(wheel_params[i])) {
            throw std::domain_error("Reaction wheel parameters must be valid");
        }
        if (!std::isfinite(wheel_states[i].wheel_speed_rad_s)
            || !std::isfinite(wheel_states[i].achieved_motor_torque_Nm)) {
            throw std::domain_error("Reaction wheel state must be finite");
        }
        total = total + wheel_momentum_body_Nms(wheel_params[i], wheel_states[i].wheel_speed_rad_s);
    }
    return total;
}

// Steps the full assembly: wheels first (reaction torque + momentum), then
// thrusters sequentially against the shared depleting mass state, then frame
// aggregation. Wheel/thruster command vectors must match array sizes exactly.
[[nodiscard]] inline ActuatorOutput step_actuator_assembly(
    const std::vector<ReactionWheelParameters>& wheel_params,
    std::vector<ReactionWheelState>& wheel_states,
    const std::vector<ThrusterParameters>& thruster_params,
    const math::Quaternion& attitude_body_to_eci,
    const ActuatorCommand& command,
    SpacecraftMassState& mass,
    double dt_s) {
    if (wheel_params.size() != wheel_states.size()) {
        throw std::domain_error("Wheel parameter and state counts must match");
    }
    if (command.wheel_torques_Nm.size() != wheel_params.size()) {
        throw std::domain_error("Wheel command count must match wheel count");
    }
    if (command.thruster_forces_N.size() != thruster_params.size()) {
        throw std::domain_error("Thruster command count must match thruster count");
    }
    if (!math::is_finite(attitude_body_to_eci)) {
        throw std::domain_error("Attitude quaternion must contain only finite values");
    }
    if (!attitude_body_to_eci.is_unit(1.0e-3)) {
        throw std::domain_error("Attitude quaternion must be approximately a unit quaternion");
    }
    if (!is_valid_mass(mass)) {
        throw std::domain_error("Spacecraft mass state must be valid (positive, dry <= wet)");
    }
    if (!std::isfinite(dt_s) || dt_s <= 0.0) {
        throw std::domain_error("Actuator step size dt must be finite and strictly positive");
    }

    ActuatorOutput output;
    math::Vector3 wheel_torque_sum{};
    for (std::size_t i = 0; i < wheel_params.size(); ++i) {
        const ReactionWheelStep step =
            step_reaction_wheel(wheel_params[i], wheel_states[i], command.wheel_torques_Nm[i], dt_s);
        wheel_states[i] = step.state;
        wheel_torque_sum = wheel_torque_sum + step.reaction_torque_body_Nm;
        output.any_torque_saturated = output.any_torque_saturated || step.torque_saturated;
        output.any_speed_saturated = output.any_speed_saturated || step.speed_saturated;
        output.any_authority_lost = output.any_authority_lost || step.authority_lost;
    }
    output.wheel_reaction_torque_body_Nm = wheel_torque_sum;

    math::Vector3 force_body_sum{};
    math::Vector3 thr_torque_sum{};
    for (std::size_t i = 0; i < thruster_params.size(); ++i) {
        const ThrusterStep step = step_thruster(
            thruster_params[i], attitude_body_to_eci, command.thruster_forces_N[i],
            mass.wet_mass_kg, mass.dry_mass_kg, dt_s);
        force_body_sum = force_body_sum + step.force_body_N;
        thr_torque_sum = thr_torque_sum + step.torque_body_Nm;
        mass.wet_mass_kg = step.spacecraft_mass_after_kg;
        output.propellant_used_kg += step.propellant_used_kg;
        output.any_thrust_saturated = output.any_thrust_saturated || step.thrust_saturated;
        output.any_propellant_depleted = output.any_propellant_depleted || step.propellant_depleted;
    }
    output.thruster_torque_body_Nm = thr_torque_sum;
    output.total_force_body_N = force_body_sum;
    output.total_force_eci_N = attitude_body_to_eci.rotate_vector(force_body_sum);
    output.total_torque_body_Nm = wheel_torque_sum + thr_torque_sum;

    const math::Vector3 wheel_momentum = total_wheel_momentum_body_Nms(wheel_params, wheel_states);
    output.total_wheel_momentum_body_x_Nms = wheel_momentum.x();
    output.total_wheel_momentum_body_y_Nms = wheel_momentum.y();
    output.total_wheel_momentum_body_z_Nms = wheel_momentum.z();
    return output;
}

}  // namespace astradock::actuators
