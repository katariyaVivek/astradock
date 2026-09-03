#pragma once

// AstraDock M14C/M14D — Thruster actuator model with propellant consumption.
//
// Physical problem:
//   Reaction wheels exchange momentum but cannot change total system momentum or
//   orbit. Thrusters expel mass to produce force (orbit/translation control) and,
//   when offset from the center of mass, torque (attitude control + wheel
//   desaturation). They consume propellant: mass depletion changes subsequent
//   force-to-acceleration mapping, so thrust modeling requires a live mass state.
//
// Frames:
//   mounting_position_body_m: thruster location in the BODY frame (m).
//   thrust_direction_body: unit vector of the exhaust-plume force direction in the
//     BODY frame (the force ON the spacecraft, not the plume direction).
//   Body outputs: force_body_N, torque_body_Nm = r_body x F_body.
//   ECI output: force_eci_N = C(q)_eci_from_body * force_body_N via the
//     spacecraft attitude quaternion.
//
// Units (SI): thrust (N), Isp (s), g0 = 9.80665 m/s^2 (exact definitional
//   constant), mass flow (kg/s), mass (kg), position (m), torque (N*m).
//
// Governing equations:
//   F_B = thrust_direction_body * thrust_achieved_N
//   tau_B = r_B x F_B
//   F_I = q_eci_from_body (*) F_B          (active rotation via quaternion)
//   m_dot = -T / (Isp * g0)                (T = thrust magnitude actually produced)
//   m(t+dt) = m(t) + m_dot * dt            (exact for piecewise-constant thrust)
//
// Assumptions and non-goals:
//   - Ideal on/off thruster: instantaneous rise, commanded thrust achieved exactly
//     up to the hard max; no minimum-impulse-bit, plume, cosine loss (fold canting
//     into thrust_direction_body), valve lag, or tank slosh/ullage modeling.
//   - Mass depletion is tracked per total thrust magnitude; Isp constant.
//   - Achieved-vs-commanded reporting: when propellant is exhausted (mass at dry
//     floor) or the command exceeds max thrust, achieved < commanded and the
//     corresponding flag is set. Commands are NEVER silently clipped.

#include "math/quaternion.hpp"
#include "math/vector3.hpp"

#include <cmath>
#include <stdexcept>

namespace astradock::actuators {

namespace detail {

inline void require_valid_thruster_parameters(
    double max_thrust_N,
    double specific_impulse_s,
    const math::Vector3& mounting_position_body_m,
    const math::Vector3& thrust_direction_body) {
    if (!std::isfinite(max_thrust_N) || max_thrust_N <= 0.0) {
        throw std::domain_error("Thruster maximum thrust must be finite and strictly positive");
    }
    if (!std::isfinite(specific_impulse_s) || specific_impulse_s <= 0.0) {
        throw std::domain_error("Thruster specific impulse must be finite and strictly positive");
    }
    if (!math::is_finite(mounting_position_body_m)) {
        throw std::domain_error("Thruster mounting position must contain only finite components");
    }
    if (!math::is_finite(thrust_direction_body)) {
        throw std::domain_error("Thruster direction must contain only finite components");
    }
    if (std::abs(thrust_direction_body.norm() - 1.0) > 1.0e-9) {
        throw std::domain_error("Thruster direction must be a unit vector in the BODY frame");
    }
}

}  // namespace detail

// Exact definitional standard gravity for Isp <-> exhaust-velocity conversion.
inline constexpr double k_standard_gravity_mps2 = 9.80665;

struct ThrusterParameters {
    double max_thrust_N{1.0};
    double specific_impulse_s{220.0};
    math::Vector3 mounting_position_body_m{};
    math::Vector3 thrust_direction_body{1.0, 0.0, 0.0};

    constexpr ThrusterParameters() noexcept = default;

    ThrusterParameters(
        double max_thrust_N_in,
        double specific_impulse_s_in,
        const math::Vector3& mounting_position_body_m_in,
        const math::Vector3& thrust_direction_body_in)
        : max_thrust_N(max_thrust_N_in),
          specific_impulse_s(specific_impulse_s_in),
          mounting_position_body_m(mounting_position_body_m_in),
          thrust_direction_body(thrust_direction_body_in) {
        detail::require_valid_thruster_parameters(
            max_thrust_N_in, specific_impulse_s_in,
            mounting_position_body_m_in, thrust_direction_body_in);
    }
};

[[nodiscard]] inline bool is_valid(const ThrusterParameters& params) noexcept {
    return std::isfinite(params.max_thrust_N) && params.max_thrust_N > 0.0
        && std::isfinite(params.specific_impulse_s) && params.specific_impulse_s > 0.0
        && math::is_finite(params.mounting_position_body_m)
        && math::is_finite(params.thrust_direction_body)
        && std::abs(params.thrust_direction_body.norm() - 1.0) <= 1.0e-9;
}

// Mass flow rate for a produced thrust magnitude: m_dot = -T / (Isp * g0) (kg/s).
[[nodiscard]] inline double thruster_mass_flow_rate_kg_per_s(
    const ThrusterParameters& params,
    double thrust_achieved_N) {
    if (!is_valid(params)) {
        throw std::domain_error("Thruster parameters must be valid");
    }
    if (!std::isfinite(thrust_achieved_N) || thrust_achieved_N < 0.0) {
        throw std::domain_error("Achieved thrust must be finite and non-negative");
    }
    return -thrust_achieved_N / (params.specific_impulse_s * k_standard_gravity_mps2);
}

// Body-frame force for an achieved thrust magnitude: F_B = direction * T.
[[nodiscard]] inline math::Vector3 thruster_force_body_N(
    const ThrusterParameters& params,
    double thrust_achieved_N) {
    if (!is_valid(params)) {
        throw std::domain_error("Thruster parameters must be valid");
    }
    if (!std::isfinite(thrust_achieved_N) || thrust_achieved_N < 0.0) {
        throw std::domain_error("Achieved thrust must be finite and non-negative");
    }
    return params.thrust_direction_body * thrust_achieved_N;
}

// Body-frame torque from mounting offset: tau_B = r_B x F_B.
[[nodiscard]] inline math::Vector3 thruster_torque_body_Nm(
    const ThrusterParameters& params,
    double thrust_achieved_N) {
    return params.mounting_position_body_m.cross(thruster_force_body_N(params, thrust_achieved_N));
}

// ECI-frame force via attitude: F_I = q_eci_from_body (*) F_B.
[[nodiscard]] inline math::Vector3 thruster_force_eci_N(
    const ThrusterParameters& params,
    const math::Quaternion& attitude_body_to_eci,
    double thrust_achieved_N) {
    if (!math::is_finite(attitude_body_to_eci)) {
        throw std::domain_error("Attitude quaternion must contain only finite values");
    }
    if (!attitude_body_to_eci.is_unit(1.0e-3)) {
        throw std::domain_error("Attitude quaternion must be approximately a unit quaternion");
    }
    return attitude_body_to_eci.rotate_vector(thruster_force_body_N(params, thrust_achieved_N));
}

// Result of firing one thruster for one step with propellant bookkeeping.
struct ThrusterStep {
    double thrust_achieved_N{0.0};      // thrust actually produced (N)
    math::Vector3 force_body_N{};       // achieved force in BODY (N)
    math::Vector3 torque_body_Nm{};     // achieved torque in BODY (N*m)
    math::Vector3 force_eci_N{};        // achieved force in ECI (N)
    double propellant_used_kg{0.0};     // mass consumed this step (kg, >= 0)
    double spacecraft_mass_after_kg{0.0};
    bool thrust_saturated{false};       // command exceeded max_thrust_N
    bool propellant_depleted{false};    // hit the dry-mass floor this step
};

// Fires a thruster for dt_s at commanded thrust with hard limits and mass update.
//   commanded_thrust_N: requested thrust magnitude (>= 0; on/off via 0 vs max).
//   spacecraft_mass_kg: wet mass at step start; dry_mass_kg is the hard floor.
// Thrust that would consume more propellant than (mass - dry) allows is scaled
// down to exactly reach the dry floor, and propellant_depleted is set.
[[nodiscard]] inline ThrusterStep step_thruster(
    const ThrusterParameters& params,
    const math::Quaternion& attitude_body_to_eci,
    double commanded_thrust_N,
    double spacecraft_mass_kg,
    double dry_mass_kg,
    double dt_s) {
    if (!is_valid(params)) {
        throw std::domain_error("Thruster parameters must be valid");
    }
    if (!std::isfinite(commanded_thrust_N) || commanded_thrust_N < 0.0) {
        throw std::domain_error("Commanded thrust must be finite and non-negative");
    }
    if (!std::isfinite(spacecraft_mass_kg) || !std::isfinite(dry_mass_kg)
        || spacecraft_mass_kg <= 0.0 || dry_mass_kg <= 0.0
        || dry_mass_kg > spacecraft_mass_kg) {
        throw std::domain_error("Spacecraft and dry masses must be finite, positive, with dry <= wet");
    }
    if (!std::isfinite(dt_s) || dt_s <= 0.0) {
        throw std::domain_error("Actuator step size dt must be finite and strictly positive");
    }

    double thrust = commanded_thrust_N;
    bool saturated = false;
    if (thrust > params.max_thrust_N) {
        thrust = params.max_thrust_N;
        saturated = true;
    }

    // Propellant-limited thrust: available mass sets the achievable impulse.
    const double mass_flow_full = thrust / (params.specific_impulse_s * k_standard_gravity_mps2);
    const double available_propellant_kg = spacecraft_mass_kg - dry_mass_kg;
    double consumed_kg = mass_flow_full * dt_s;
    bool depleted = false;
    if (consumed_kg > available_propellant_kg) {
        consumed_kg = available_propellant_kg;
        thrust = (dt_s > 0.0) ? consumed_kg * params.specific_impulse_s * k_standard_gravity_mps2 / dt_s : 0.0;
        depleted = true;
    }

    ThrusterStep result;
    result.thrust_achieved_N = thrust;
    result.force_body_N = thruster_force_body_N(params, thrust);
    result.torque_body_Nm = thruster_torque_body_Nm(params, thrust);
    result.force_eci_N = thruster_force_eci_N(params, attitude_body_to_eci, thrust);
    result.propellant_used_kg = consumed_kg;
    result.spacecraft_mass_after_kg = spacecraft_mass_kg - consumed_kg;
    result.thrust_saturated = saturated;
    result.propellant_depleted = depleted;
    return result;
}

}  // namespace astradock::actuators
