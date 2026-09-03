#pragma once

// AstraDock M15A/M15B — LVLH angular velocity and rotating-frame kinematics.
//
// Physical problem:
//   M06 builds the LVLH basis and projects vectors between ECI and LVLH, but only
//   geometrically: it treats the frame as frozen at one instant. The LVLH frame
//   actually rotates as the reference spacecraft orbits. Expressing chaser motion
//   relative to a target (M15C) therefore requires the frame's angular velocity
//   omega_LVLH/ECI and the transport theorem that separates the inertial
//   derivative seen in ECI from the derivative seen by an observer rotating
//   with LVLH.
//
// Frames:
//   All vectors are 3-component ECI-coordinate triples unless the name says
//   `_lvlh`. `lvlh_angular_velocity_rad_s` returns omega of the LVLH frame
//   relative to ECI, expressed in ECI coordinates. The rotating-frame transport
//   helpers take explicit omega (ECI components) so the frame pairing can never
//   be implicit.
//
// Units (SI): positions (m), velocities (m/s), accelerations (m/s^2),
//   angular velocity (rad/s), angular acceleration (rad/s^2).
//
// Governing equations:
//   Orbital angular velocity of the reference point (exact, any eccentricity):
//     h_vec = r x v,  h = |h_vec|
//     omega = h_vec / |r|^2
//   because v_r = (r.v/r) r_hat is parallel to r and drops out of r x v, so
//   r x v = r x (omega x r) = omega |r|^2 for omega perpendicular to r.
//   Circular special case: |omega| = n = sqrt(mu / R^3).
//   Transport theorem (first order):
//     (dv/dt)_inertial = (dv/dt)_rotating + omega x v
//   Transport theorem (second order, acceleration composition):
//     a_inertial = a_rel + alpha x r_rel + 2 omega x v_rel
//                         + omega x (omega x r_rel) + a_frame_origin
//   with Coriolis 2 omega x v_rel and centrifugal omega x (omega x r_rel).
//
// Assumptions and non-goals:
//   - The reference ECI state defines an osculating orbital frame; no
//     perturbation model is assumed. omega follows purely from (r, v).
//   - Collinear/zero-momentum states are rejected: the orbital plane (and hence
//     LVLH) is undefined there, exactly as in M06.
//   - No guidance, control, or CW linearization in this header (M15C/M15D).

#include "frames/lvlh.hpp"
#include "math/matrix3.hpp"
#include "math/vector3.hpp"

#include <cmath>
#include <stdexcept>

namespace astradock::frames {

namespace detail {

inline void require_valid_reference_state(
    const math::Vector3& position_eci_m,
    const math::Vector3& velocity_eci_mps) {
    if (!math::is_finite(position_eci_m) || !math::is_finite(velocity_eci_mps)) {
        throw std::domain_error("Reference orbital state must contain only finite values");
    }
    if (position_eci_m.norm() == 0.0) {
        throw std::domain_error("Reference position is undefined at the central-body origin");
    }
    if (position_eci_m.cross(velocity_eci_mps).norm() == 0.0) {
        throw std::domain_error(
            "Reference angular momentum is zero (collinear position and velocity); "
            "LVLH angular velocity is undefined");
    }
}

inline void require_finite_kinematics(
    const math::Vector3& omega_rad_s,
    const math::Vector3& a,
    const math::Vector3& b) {
    if (!math::is_finite(omega_rad_s) || !math::is_finite(a) || !math::is_finite(b)) {
        throw std::domain_error("Rotating-frame kinematic inputs must contain only finite values");
    }
}

}  // namespace detail

// Angular velocity of the LVLH frame relative to ECI, expressed in ECI
// coordinates: omega_LVLH/ECI = h_vec / |r|^2, with h_vec = r x v.
//
// Derivation: split v = v_radial + omega x r with v_radial parallel to r.
// Then r x v = r x (omega x r) = omega (r.r) - r (r.omega) = omega |r|^2
// since omega (along h) is perpendicular to r. Valid for eccentric orbits;
// for circular orbits |omega| reduces to the mean motion n = sqrt(mu/R^3).
[[nodiscard]] inline math::Vector3 lvlh_angular_velocity_rad_s(
    const math::Vector3& position_eci_m,
    const math::Vector3& velocity_eci_mps) {
    detail::require_valid_reference_state(position_eci_m, velocity_eci_mps);
    const double r_squared = position_eci_m.squared_norm();
    return position_eci_m.cross(velocity_eci_mps) / r_squared;
}

// Magnitude of the LVLH angular velocity (rad/s). For a circular orbit this is
// the mean motion n.
[[nodiscard]] inline double lvlh_rate_rad_per_s(
    const math::Vector3& position_eci_m,
    const math::Vector3& velocity_eci_mps) {
    return lvlh_angular_velocity_rad_s(position_eci_m, velocity_eci_mps).norm();
}

// First-order transport theorem: splits an inertial time derivative into the
// derivative seen in the rotating frame plus the frame-rotation term.
//
//   (dv/dt)_inertial = (dv/dt)_rotating + omega x v
//
// Inputs: v (any physical vector, ECI components), omega (frame rate, ECI
// components), dv_rotating (time derivative observed in the rotating frame,
// ECI components). Output: inertial derivative, ECI components.
[[nodiscard]] inline math::Vector3 transport_first_derivative(
    const math::Vector3& vector_eci,
    const math::Vector3& omega_rad_s,
    const math::Vector3& derivative_in_rotating_frame) {
    detail::require_finite_kinematics(omega_rad_s, vector_eci, derivative_in_rotating_frame);
    return derivative_in_rotating_frame + omega_rad_s.cross(vector_eci);
}

// Inverse split: rotating-frame derivative from a known inertial derivative.
[[nodiscard]] inline math::Vector3 rotating_frame_derivative(
    const math::Vector3& vector_eci,
    const math::Vector3& omega_rad_s,
    const math::Vector3& derivative_in_inertial_frame) {
    detail::require_finite_kinematics(omega_rad_s, vector_eci, derivative_in_inertial_frame);
    return derivative_in_inertial_frame - omega_rad_s.cross(vector_eci);
}

// Coriolis acceleration for relative velocity v_rel in a frame rotating at
// omega: a_cor = 2 omega x v_rel (ECI components in/out).
[[nodiscard]] inline math::Vector3 coriolis_acceleration(
    const math::Vector3& omega_rad_s,
    const math::Vector3& relative_velocity) {
    detail::require_finite_kinematics(omega_rad_s, relative_velocity, relative_velocity);
    return omega_rad_s.cross(relative_velocity) * 2.0;
}

// Centrifugal acceleration at relative position r_rel: a_cen = omega x (omega x r).
[[nodiscard]] inline math::Vector3 centrifugal_acceleration(
    const math::Vector3& omega_rad_s,
    const math::Vector3& relative_position) {
    detail::require_finite_kinematics(omega_rad_s, relative_position, relative_position);
    return omega_rad_s.cross(omega_rad_s.cross(relative_position));
}

// Euler acceleration from frame angular acceleration: a_eul = alpha x r_rel.
[[nodiscard]] inline math::Vector3 euler_acceleration(
    const math::Vector3& angular_acceleration_rad_per_s2,
    const math::Vector3& relative_position) {
    if (!math::is_finite(angular_acceleration_rad_per_s2) || !math::is_finite(relative_position)) {
        throw std::domain_error("Euler acceleration inputs must contain only finite values");
    }
    return angular_acceleration_rad_per_s2.cross(relative_position);
}

// Full second-order composition: inertial acceleration of a point tracked in a
// rotating frame whose origin itself accelerates at origin_accel_inertial:
//
//   a_inertial = a_rel + alpha x r_rel + 2 omega x v_rel
//                       + omega x (omega x r_rel) + a_origin
[[nodiscard]] inline math::Vector3 compose_inertial_acceleration(
    const math::Vector3& relative_position,
    const math::Vector3& relative_velocity,
    const math::Vector3& relative_acceleration,
    const math::Vector3& omega_rad_s,
    const math::Vector3& angular_acceleration_rad_per_s2,
    const math::Vector3& origin_acceleration_inertial) {
    if (!math::is_finite(relative_position) || !math::is_finite(relative_velocity)
        || !math::is_finite(relative_acceleration) || !math::is_finite(omega_rad_s)
        || !math::is_finite(angular_acceleration_rad_per_s2)
        || !math::is_finite(origin_acceleration_inertial)) {
        throw std::domain_error("Acceleration composition inputs must contain only finite values");
    }
    return relative_acceleration + euler_acceleration(angular_acceleration_rad_per_s2, relative_position)
        + coriolis_acceleration(omega_rad_s, relative_velocity)
        + centrifugal_acceleration(omega_rad_s, relative_position) + origin_acceleration_inertial;
}

// Time derivative of the LVLH angular velocity (frame angular acceleration)
// from reference-state acceleration, via the quotient rule on omega = h/r^2:
//
//   d(omega)/dt = (dh/dt)/r^2 - h (2 r.v / r^4),
//   dh/dt = r x a_ref (central + perturbative acceleration of the reference).
//
// For Keplerian circular motion dh/dt = 0 and r.v = 0, so alpha = 0 exactly.
[[nodiscard]] inline math::Vector3 lvlh_angular_acceleration_rad_per_s2(
    const math::Vector3& position_eci_m,
    const math::Vector3& velocity_eci_mps,
    const math::Vector3& reference_acceleration_eci_mps2) {
    detail::require_valid_reference_state(position_eci_m, velocity_eci_mps);
    if (!math::is_finite(reference_acceleration_eci_mps2)) {
        throw std::domain_error("Reference acceleration must contain only finite values");
    }
    const double r_squared = position_eci_m.squared_norm();
    const math::Vector3 h_vec = position_eci_m.cross(velocity_eci_mps);
    const math::Vector3 dh_dt = position_eci_m.cross(reference_acceleration_eci_mps2);
    const double r_dot_v = position_eci_m.dot(velocity_eci_mps);
    return dh_dt / r_squared - h_vec * (2.0 * r_dot_v / (r_squared * r_squared));
}

// Kinematic consistency audit: reconstructs the ECI velocity of the reference
// point from the LVLH frame motion, v_check = omega x r. For the reference
// point itself the rotating-frame derivative of r vanishes in the radial/tangential
// sense only for circular orbits; the exact identity is:
//
//   v = (r.v/r) r_hat + omega x r,
//
// i.e. radial rate plus frame transport. Returns v minus that reconstruction
// (zero iff the omega definition is consistent with the state).
[[nodiscard]] inline math::Vector3 lvlh_rate_consistency_residual(
    const math::Vector3& position_eci_m,
    const math::Vector3& velocity_eci_mps) {
    detail::require_valid_reference_state(position_eci_m, velocity_eci_mps);
    const math::Vector3 omega = lvlh_angular_velocity_rad_s(position_eci_m, velocity_eci_mps);
    const math::Vector3 r_hat = position_eci_m / position_eci_m.norm();
    const math::Vector3 radial_velocity = r_hat * position_eci_m.dot(velocity_eci_mps) / position_eci_m.norm();
    return velocity_eci_mps - radial_velocity - omega.cross(position_eci_m);
}

}  // namespace astradock::frames
