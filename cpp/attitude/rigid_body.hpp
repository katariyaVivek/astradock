#pragma once

#include "attitude/principal_inertia.hpp"
#include "attitude/rotational_state.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"
#include "numerics/integrators.hpp"

#include <cmath>
#include <concepts>
#include <functional>
#include <stdexcept>

namespace astradock::attitude {

namespace detail {

inline void require_finite_angular_velocity(const math::Vector3& omega) {
    if (!math::is_finite(omega)) {
        throw std::domain_error("Angular velocity must contain only finite components");
    }
}

inline void require_finite_torque(const math::Vector3& torque) {
    if (!math::is_finite(torque)) {
        throw std::domain_error("Applied torque must contain only finite components");
    }
}

inline void require_valid_inertia(const PrincipalInertia& inertia) {
    if (!is_valid(inertia)) {
        throw std::domain_error("Principal inertia moments must be finite and strictly positive");
    }
}

inline void require_valid_rotational_state(const RotationalState& state) {
    if (!is_finite(state)) {
        throw std::domain_error("Rotational state must contain only finite values");
    }
    if (!state.orientation.is_unit(1.0e-3)) {
        throw std::domain_error("Attitude quaternion must be approximately a unit quaternion");
    }
}

}  // namespace detail

// Computes the time derivative of the attitude quaternion from body angular velocity.
//
// Governing equation:
//   dq/dt = 0.5 * q ⊗ [0, omega_B]
//
// Frame & Convention:
//   q is scalar-first [w, x, y, z] mapping coordinates from the spacecraft body
//   frame B into the inertial reference frame I: v_I = q * [0, v_B] * q*.
//   omega_B is the angular velocity of frame B relative to frame I, expressed
//   in frame B coordinates (rad/s).
[[nodiscard]] inline math::Quaternion quaternion_derivative(
    const math::Quaternion& orientation,
    const math::Vector3& angular_velocity_body_rad_per_s) {
    if (!math::is_finite(orientation)) {
        throw std::domain_error("Attitude quaternion must contain only finite values");
    }
    detail::require_finite_angular_velocity(angular_velocity_body_rad_per_s);

    const math::Quaternion pure_omega(
        0.0,
        angular_velocity_body_rad_per_s.x(),
        angular_velocity_body_rad_per_s.y(),
        angular_velocity_body_rad_per_s.z()
    );

    return (orientation * pure_omega) * 0.5;
}

// Computes rigid-body angular acceleration using Euler's equations of rotational motion
// in the spacecraft principal-axis body-fixed frame.
//
// Governing equation:
//   I * d(omega)/dt + omega x (I * omega) = tau
//   => d(omega)/dt = I^(-1) * [ tau - omega x (I * omega) ]
//
// Principal axes expansion:
//   d(omega_x)/dt = [ tau_x - (Izz - Iyy) * omega_y * omega_z ] / Ixx
//   d(omega_y)/dt = [ tau_y - (Ixx - Izz) * omega_z * omega_x ] / Iyy
//   d(omega_z)/dt = [ tau_z - (Iyy - Ixx) * omega_x * omega_y ] / Izz
//
// Inputs:
//   angular_velocity_body_rad_per_s: omega expressed in body frame (rad/s).
//   inertia: principal moments of inertia Ixx, Iyy, Izz (kg*m^2).
//   torque_body_Nm: external torque expressed in body frame (N*m).
// Output:
//   angular acceleration d(omega)/dt in body frame (rad/s^2).
[[nodiscard]] inline math::Vector3 euler_rotational_acceleration(
    const math::Vector3& angular_velocity_body_rad_per_s,
    const PrincipalInertia& inertia,
    const math::Vector3& torque_body_Nm) {
    detail::require_finite_angular_velocity(angular_velocity_body_rad_per_s);
    detail::require_valid_inertia(inertia);
    detail::require_finite_torque(torque_body_Nm);

    const double wx = angular_velocity_body_rad_per_s.x();
    const double wy = angular_velocity_body_rad_per_s.y();
    const double wz = angular_velocity_body_rad_per_s.z();

    // Body angular momentum H_B = I * omega
    const math::Vector3 h_body(
        inertia.Ixx_kg_m2 * wx,
        inertia.Iyy_kg_m2 * wy,
        inertia.Izz_kg_m2 * wz
    );

    // Gyroscopic coupling term: omega x (I * omega)
    const math::Vector3 gyroscopic_torque = angular_velocity_body_rad_per_s.cross(h_body);

    // Net torque available for angular acceleration: tau_ext - omega x (I * omega)
    const math::Vector3 net_torque = torque_body_Nm - gyroscopic_torque;

    const math::Vector3 alpha(
        net_torque.x() / inertia.Ixx_kg_m2,
        net_torque.y() / inertia.Iyy_kg_m2,
        net_torque.z() / inertia.Izz_kg_m2
    );

    if (!math::is_finite(alpha)) {
        throw std::overflow_error("Euler rotational acceleration calculation produced non-finite values");
    }

    return alpha;
}

// Computes the combined first-order time derivative of the full rotational state.
//
// Returned state slots:
//   orientation slot: dq/dt (quaternion kinematics derivative)
//   angular_velocity slot: d(omega)/dt in rad/s^2 (Euler rotational dynamics)
[[nodiscard]] inline RotationalState rotational_state_derivative(
    double /*time_s*/,
    const RotationalState& state,
    const PrincipalInertia& inertia,
    const math::Vector3& torque_body_Nm) {
    if (!is_finite(state)) {
        throw std::domain_error("Rotational state must contain only finite values");
    }

    return {
        quaternion_derivative(state.orientation, state.angular_velocity_rad_per_s),
        euler_rotational_acceleration(state.angular_velocity_rad_per_s, inertia, torque_body_Nm)
    };
}

// Computes the rotational kinetic energy of the rigid body:
//   E_rot = 0.5 * omega^T * I * omega = 0.5 * (Ixx * wx^2 + Iyy * wy^2 + Izz * wz^2)
//
// In torque-free motion (tau = 0), E_rot is an exact analytical invariant.
// Units: Joules (J = kg * m^2 / s^2).
[[nodiscard]] inline double rotational_kinetic_energy_J(
    const math::Vector3& angular_velocity_body_rad_per_s,
    const PrincipalInertia& inertia) {
    detail::require_finite_angular_velocity(angular_velocity_body_rad_per_s);
    detail::require_valid_inertia(inertia);

    const double wx = angular_velocity_body_rad_per_s.x();
    const double wy = angular_velocity_body_rad_per_s.y();
    const double wz = angular_velocity_body_rad_per_s.z();

    const double energy_J = 0.5 * (
        inertia.Ixx_kg_m2 * wx * wx
        + inertia.Iyy_kg_m2 * wy * wy
        + inertia.Izz_kg_m2 * wz * wz
    );

    if (!std::isfinite(energy_J)) {
        throw std::overflow_error("Rotational kinetic energy calculation is not representable");
    }

    return energy_J;
}

// Computes angular momentum expressed in spacecraft BODY coordinates:
//   H_body = I * omega = [Ixx * wx, Iyy * wy, Izz * wz]^T
//
// Note: In asymmetric torque-free tumbling, H_body components vary over time
// because the body frame itself rotates relative to inertial space.
// Units: N * m * s (or kg * m^2 / s).
[[nodiscard]] inline math::Vector3 body_angular_momentum_kg_m2_per_s(
    const math::Vector3& angular_velocity_body_rad_per_s,
    const PrincipalInertia& inertia) {
    detail::require_finite_angular_velocity(angular_velocity_body_rad_per_s);
    detail::require_valid_inertia(inertia);

    const math::Vector3 h_body(
        inertia.Ixx_kg_m2 * angular_velocity_body_rad_per_s.x(),
        inertia.Iyy_kg_m2 * angular_velocity_body_rad_per_s.y(),
        inertia.Izz_kg_m2 * angular_velocity_body_rad_per_s.z()
    );

    if (!math::is_finite(h_body)) {
        throw std::overflow_error("Body angular momentum calculation produced non-finite values");
    }

    return h_body;
}

// Computes angular momentum transformed into the INERTIAL reference frame:
//   H_inertial = C_{I_B} * H_body = q ⊗ [0, H_body] ⊗ q*
//
// In torque-free motion (tau = 0), H_inertial is a conserved constant vector.
// Units: N * m * s (or kg * m^2 / s).
[[nodiscard]] inline math::Vector3 inertial_angular_momentum_kg_m2_per_s(
    const RotationalState& state,
    const PrincipalInertia& inertia) {
    detail::require_valid_rotational_state(state);
    detail::require_valid_inertia(inertia);

    const math::Vector3 h_body = body_angular_momentum_kg_m2_per_s(
        state.angular_velocity_rad_per_s, inertia
    );

    return state.orientation.rotate_vector(h_body);
}

// Advances the rotational state by one fixed RK4 step under constant external torque.
//
// Normalization Policy:
//   Intermediate RK4 substages operate in unconstrained linear state space.
//   At the completion of the full RK4 step, the attitude quaternion is
//   reprojected onto the unit hypersphere S^3 (if renormalize_quaternion is true)
//   to prevent numerical integration drift over long durations.
[[nodiscard]] inline RotationalState rk4_step_rotational(
    double t,
    const RotationalState& state,
    double dt,
    const PrincipalInertia& inertia,
    const math::Vector3& constant_torque_body_Nm,
    bool renormalize_quaternion = true) {
    detail::require_valid_inertia(inertia);

    const auto derivative_evaluator = [&](double step_t, const RotationalState& step_state) -> RotationalState {
        return rotational_state_derivative(step_t, step_state, inertia, constant_torque_body_Nm);
    };

    const RotationalState raw_next_state = numerics::rk4_step(
        t,
        state,
        dt,
        derivative_evaluator
    );

    return renormalize_quaternion ? raw_next_state.normalized() : raw_next_state;
}

// Overload for time- or state-varying torque callables: torque_func(t, state) -> Vector3.
template <typename TorqueCallable>
requires std::invocable<TorqueCallable, double, const RotationalState&>
[[nodiscard]] inline RotationalState rk4_step_rotational(
    double t,
    const RotationalState& state,
    double dt,
    const PrincipalInertia& inertia,
    TorqueCallable&& torque_callable,
    bool renormalize_quaternion = true) {
    detail::require_valid_inertia(inertia);

    const auto derivative_evaluator = [&](double step_t, const RotationalState& step_state) -> RotationalState {
        const math::Vector3 torque = std::invoke(torque_callable, step_t, step_state);
        return rotational_state_derivative(step_t, step_state, inertia, torque);
    };

    const RotationalState raw_next_state = numerics::rk4_step(
        t,
        state,
        dt,
        derivative_evaluator
    );

    return renormalize_quaternion ? raw_next_state.normalized() : raw_next_state;
}

}  // namespace astradock::attitude
