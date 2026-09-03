#pragma once

// AstraDock M14A — Reaction wheel actuator dynamics.
//
// Physical problem:
//   Spacecraft attitude must be controlled without expending propellant. A reaction
//   wheel is an electric-motor-driven flywheel rigidly mounted to the spacecraft bus.
//   Accelerating the wheel in one direction applies an equal-and-opposite reaction
//   torque to the bus (Newton's third law), exchanging angular momentum between the
//   wheel and the bus while total system momentum is conserved.
//
// Frames:
//   The wheel spins about a fixed unit spin axis expressed in the spacecraft BODY
//   frame (spin_axis_body). Scalar quantities (speed, motor torque, stored momentum
//   magnitude) are signed about +spin_axis_body. Vector outputs
//   (reaction_torque_body_Nm, stored_momentum_body_Nms) are BODY-frame vectors.
//
// Units (SI):
//   wheel_inertia_kg_m2 (kg*m^2), torques (N*m), wheel speed (rad/s),
//   angular momentum (N*m*s), kinetic energy (J), time (s).
//
// Governing equations:
//   Wheel spin:            I_w * d(omega_w)/dt = tau_motor_achieved
//   Stored momentum:       H_w = I_w * omega_w            (scalar about +spin axis)
//   Stored energy:         E_w = 0.5 * I_w * omega_w^2
//   Spacecraft reaction:   tau_body = -spin_axis_body * tau_motor_achieved
//   Motor lag (1st order): tau_achieved(t) = tau_target * (1 - exp(-t / tau_lag))
//   System conservation:   d/dt (H_bus + H_wheels)_inertial = tau_external_only
//
// Assumptions and non-goals:
//   - Rigid, balanced wheel with a single spin degree of freedom.
//   - Motor torque and speed limits are hard; exceeding them saturates (flagged)
//     rather than damaging the model. No friction, cogging, thermal, power-bus,
//     bearing-noise, or structural-flexibility modeling.
//   - Speed-limit saturation models authority loss: a torque command that would
//     drive the wheel further past its speed limit produces no acceleration and
//     raises authority_lost. The caller (later: momentum management / thrusters)
//     must handle desaturation. Actuator commands are NEVER silently clipped:
//     every limiting event sets an explicit flag.
//
// Numerical method:
//   Constant-torque steps are integrated exactly (omega += tau * dt / I_w) before
//   the hard speed clamp, so the only "error" is the modeled limiter itself.
//   Motor lag uses the exact discrete update of a first-order low-pass filter.

#include "math/vector3.hpp"

#include <cmath>
#include <stdexcept>

namespace astradock::actuators {

namespace detail {

inline void require_valid_wheel_parameters(
    double inertia_kg_m2,
    double max_torque_Nm,
    double max_speed_rad_s,
    double time_constant_s,
    const math::Vector3& spin_axis) {
    if (!std::isfinite(inertia_kg_m2) || inertia_kg_m2 <= 0.0) {
        throw std::domain_error("Reaction wheel inertia must be finite and strictly positive");
    }
    if (!std::isfinite(max_torque_Nm) || max_torque_Nm <= 0.0) {
        throw std::domain_error("Reaction wheel maximum motor torque must be finite and strictly positive");
    }
    if (!std::isfinite(max_speed_rad_s) || max_speed_rad_s <= 0.0) {
        throw std::domain_error("Reaction wheel maximum speed must be finite and strictly positive");
    }
    if (!std::isfinite(time_constant_s) || time_constant_s < 0.0) {
        throw std::domain_error("Reaction wheel motor time constant must be finite and non-negative");
    }
    if (!math::is_finite(spin_axis)) {
        throw std::domain_error("Reaction wheel spin axis must contain only finite components");
    }
    if (std::abs(spin_axis.norm() - 1.0) > 1.0e-9) {
        throw std::domain_error("Reaction wheel spin axis must be a unit vector in the BODY frame");
    }
}

inline void require_valid_command(double commanded_torque_Nm, double dt_s) {
    if (!std::isfinite(commanded_torque_Nm)) {
        throw std::domain_error("Commanded wheel torque must be finite");
    }
    if (!std::isfinite(dt_s) || dt_s <= 0.0) {
        throw std::domain_error("Actuator step size dt must be finite and strictly positive");
    }
}

}  // namespace detail

// Static physical parameters of one reaction wheel. A time constant of 0 selects
// an ideal instantaneous motor response; any positive value models first-order lag.
struct ReactionWheelParameters {
    double wheel_inertia_kg_m2{1.0e-3};
    double max_motor_torque_Nm{0.05};
    double max_wheel_speed_rad_s{500.0};
    double motor_time_constant_s{0.0};
    math::Vector3 spin_axis_body{0.0, 0.0, 1.0};

    constexpr ReactionWheelParameters() noexcept = default;

    ReactionWheelParameters(
        double inertia_kg_m2,
        double max_torque_Nm,
        double max_speed_rad_s,
        double time_constant_s,
        const math::Vector3& spin_axis_body_in)
        : wheel_inertia_kg_m2(inertia_kg_m2),
          max_motor_torque_Nm(max_torque_Nm),
          max_wheel_speed_rad_s(max_speed_rad_s),
          motor_time_constant_s(time_constant_s),
          spin_axis_body(spin_axis_body_in) {
        detail::require_valid_wheel_parameters(
            inertia_kg_m2, max_torque_Nm, max_speed_rad_s, time_constant_s, spin_axis_body_in);
    }
};

[[nodiscard]] inline bool is_valid(const ReactionWheelParameters& params) noexcept {
    return std::isfinite(params.wheel_inertia_kg_m2) && params.wheel_inertia_kg_m2 > 0.0
        && std::isfinite(params.max_motor_torque_Nm) && params.max_motor_torque_Nm > 0.0
        && std::isfinite(params.max_wheel_speed_rad_s) && params.max_wheel_speed_rad_s > 0.0
        && std::isfinite(params.motor_time_constant_s) && params.motor_time_constant_s >= 0.0
        && math::is_finite(params.spin_axis_body)
        && std::abs(params.spin_axis_body.norm() - 1.0) <= 1.0e-9;
}

// Dynamic state of one reaction wheel: signed speed about +spin_axis_body and the
// motor torque actually achieved at the end of the previous step (lag memory).
struct ReactionWheelState {
    double wheel_speed_rad_s{0.0};
    double achieved_motor_torque_Nm{0.0};
};

// Result of hard torque limiting. `saturated == true` iff the command magnitude
// exceeded the motor capability; the limited value is what the motor was asked for.
struct TorqueSaturation {
    double limited_torque_Nm{0.0};
    bool saturated{false};
};

// Applies the hard motor torque limit. Never throws away the fact of limiting:
// the caller must inspect `saturated`.
[[nodiscard]] inline TorqueSaturation saturate_motor_torque(
    const ReactionWheelParameters& params,
    double commanded_torque_Nm) {
    if (!is_valid(params)) {
        throw std::domain_error("Reaction wheel parameters must be valid");
    }
    if (!std::isfinite(commanded_torque_Nm)) {
        throw std::domain_error("Commanded wheel torque must be finite");
    }
    const double limit = params.max_motor_torque_Nm;
    if (commanded_torque_Nm > limit) {
        return {limit, true};
    }
    if (commanded_torque_Nm < -limit) {
        return {-limit, true};
    }
    return {commanded_torque_Nm, false};
}

// Exact discrete update of first-order motor lag under a constant torque target:
//   achieved += (target - achieved) * (1 - exp(-dt / tau_lag)).
// A zero time constant gives the ideal instantaneous response.
[[nodiscard]] inline double apply_motor_lag(
    double previous_achieved_torque_Nm,
    double target_torque_Nm,
    double dt_s,
    double time_constant_s) {
    if (!std::isfinite(previous_achieved_torque_Nm) || !std::isfinite(target_torque_Nm)) {
        throw std::domain_error("Motor torques must be finite");
    }
    detail::require_valid_command(0.0, dt_s);
    if (!std::isfinite(time_constant_s) || time_constant_s < 0.0) {
        throw std::domain_error("Motor time constant must be finite and non-negative");
    }
    if (time_constant_s == 0.0) {
        return target_torque_Nm;
    }
    const double blend = 1.0 - std::exp(-dt_s / time_constant_s);
    return previous_achieved_torque_Nm + (target_torque_Nm - previous_achieved_torque_Nm) * blend;
}

// Result of advancing the wheel speed under a constant achieved motor torque.
struct WheelSpeedStep {
    double wheel_speed_rad_s{0.0};  // updated speed, clamped to +-max_wheel_speed_rad_s
    double applied_torque_Nm{0.0};  // torque that actually changed wheel momentum (0 when
                                    // driven into the speed limit with no room to accelerate)
    bool speed_saturated{false};    // wheel is now resting on its speed limit
    bool authority_lost{false};     // at the limit AND the torque pushes further outward,
                                    // so the requested momentum change was refused
};

// Exact constant-torque speed update followed by the hard speed clamp. The applied
// torque is back-computed from the achieved speed change so that callers can close
// the momentum books: I_w * d(omega_w) == applied_torque * dt exactly.
[[nodiscard]] inline WheelSpeedStep step_wheel_speed(
    const ReactionWheelParameters& params,
    double wheel_speed_rad_s,
    double achieved_motor_torque_Nm,
    double dt_s) {
    if (!is_valid(params)) {
        throw std::domain_error("Reaction wheel parameters must be valid");
    }
    if (!std::isfinite(wheel_speed_rad_s) || !std::isfinite(achieved_motor_torque_Nm)) {
        throw std::domain_error("Wheel speed and motor torque must be finite");
    }
    detail::require_valid_command(0.0, dt_s);

    const double unconstrained =
        wheel_speed_rad_s + achieved_motor_torque_Nm * dt_s / params.wheel_inertia_kg_m2;
    const double limit = params.max_wheel_speed_rad_s;

    double clamped = unconstrained;
    if (clamped > limit) {
        clamped = limit;
    } else if (clamped < -limit) {
        clamped = -limit;
    }

    // speed_saturated reports whether the wheel now rests on its speed limit.
    // authority_lost reports whether the requested momentum change was refused:
    // the wheel started the step already at the limit and the torque demanded
    // motion further outward, so no speed change was possible. Landing exactly on
    // the limit from the interior still applied the full requested change.
    const bool at_limit = (clamped == limit) || (clamped == -limit);
    const bool pushing_outward = (clamped == limit && achieved_motor_torque_Nm > 0.0)
        || (clamped == -limit && achieved_motor_torque_Nm < 0.0);
    const bool already_at_limit = (wheel_speed_rad_s == limit) || (wheel_speed_rad_s == -limit);
    const bool refused = already_at_limit && pushing_outward && (clamped == wheel_speed_rad_s);

    const double applied =
        (clamped - wheel_speed_rad_s) * params.wheel_inertia_kg_m2 / dt_s;

    return {clamped, applied, at_limit, refused};}

// Scalar stored angular momentum about +spin_axis_body: H_w = I_w * omega_w (N*m*s).
[[nodiscard]] inline double wheel_momentum_scalar_Nms(
    const ReactionWheelParameters& params,
    double wheel_speed_rad_s) {
    if (!is_valid(params)) {
        throw std::domain_error("Reaction wheel parameters must be valid");
    }
    if (!std::isfinite(wheel_speed_rad_s)) {
        throw std::domain_error("Wheel speed must be finite");
    }
    return params.wheel_inertia_kg_m2 * wheel_speed_rad_s;
}

// Stored angular momentum as a BODY-frame vector: H_w_vec = spin_axis * I_w * omega_w.
[[nodiscard]] inline math::Vector3 wheel_momentum_body_Nms(
    const ReactionWheelParameters& params,
    double wheel_speed_rad_s) {
    return params.spin_axis_body * wheel_momentum_scalar_Nms(params, wheel_speed_rad_s);
}

// Spacecraft reaction torque (BODY frame) for an achieved motor torque.
// Sign convention: a positive motor torque accelerates the wheel toward
// +spin_axis_body, so the bus feels tau_body = -spin_axis_body * tau_motor.
[[nodiscard]] inline math::Vector3 reaction_torque_body_Nm(
    const ReactionWheelParameters& params,
    double achieved_motor_torque_Nm) {
    if (!is_valid(params)) {
        throw std::domain_error("Reaction wheel parameters must be valid");
    }
    if (!std::isfinite(achieved_motor_torque_Nm)) {
        throw std::domain_error("Achieved motor torque must be finite");
    }
    return params.spin_axis_body * (-achieved_motor_torque_Nm);
}

// Stored kinetic energy: E_w = 0.5 * I_w * omega_w^2 (J). Needed for M14B momentum
// buildup analysis (energy cost of saturation approach).
[[nodiscard]] inline double wheel_kinetic_energy_J(
    const ReactionWheelParameters& params,
    double wheel_speed_rad_s) {
    if (!is_valid(params)) {
        throw std::domain_error("Reaction wheel parameters must be valid");
    }
    if (!std::isfinite(wheel_speed_rad_s)) {
        throw std::domain_error("Wheel speed must be finite");
    }
    return 0.5 * params.wheel_inertia_kg_m2 * wheel_speed_rad_s * wheel_speed_rad_s;
}

// Complete single-wheel step: saturate command -> motor lag -> speed integration,
// plus the BODY reaction torque and saturation bookkeeping for this step.
struct ReactionWheelStep {
    ReactionWheelState state{};
    math::Vector3 reaction_torque_body_Nm{};
    double stored_momentum_Nms{0.0};
    bool torque_saturated{false};
    bool speed_saturated{false};
    bool authority_lost{false};
};

[[nodiscard]] inline ReactionWheelStep step_reaction_wheel(
    const ReactionWheelParameters& params,
    const ReactionWheelState& state,
    double commanded_torque_Nm,
    double dt_s) {
    if (!is_valid(params)) {
        throw std::domain_error("Reaction wheel parameters must be valid");
    }
    if (!std::isfinite(state.wheel_speed_rad_s) || !std::isfinite(state.achieved_motor_torque_Nm)) {
        throw std::domain_error("Reaction wheel state must be finite");
    }
    detail::require_valid_command(commanded_torque_Nm, dt_s);

    const TorqueSaturation limited = saturate_motor_torque(params, commanded_torque_Nm);
    const double achieved = apply_motor_lag(
        state.achieved_motor_torque_Nm, limited.limited_torque_Nm, dt_s, params.motor_time_constant_s);
    const WheelSpeedStep speed = step_wheel_speed(params, state.wheel_speed_rad_s, achieved, dt_s);

    ReactionWheelStep result;
    result.state = {speed.wheel_speed_rad_s, achieved};
    // The torque the bus actually feels is the reaction to the momentum change
    // imparted to the wheel: tau_bus = -applied * spin_axis. When authority is
    // lost (driven into the speed limit), applied == 0 and the bus feels nothing.
    result.reaction_torque_body_Nm = params.spin_axis_body * (-speed.applied_torque_Nm);
    result.stored_momentum_Nms = wheel_momentum_scalar_Nms(params, speed.wheel_speed_rad_s);
    result.torque_saturated = limited.saturated;
    result.speed_saturated = speed.speed_saturated;
    result.authority_lost = speed.authority_lost;
    return result;
}

}  // namespace astradock::actuators
