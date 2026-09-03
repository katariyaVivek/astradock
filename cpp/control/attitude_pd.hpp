#pragma once

// AstraDock M16A/M16B — Quaternion-error PD attitude control + rate damping.
//
// Physical problem:
//   Detumble a tumbling bus and hold/track an attitude using body torque. The
//   error lives on SO(3), not in a vector space: subtracting quaternion
//   components is geometrically meaningless (it breaks the unit norm and mixes
//   the double cover), so the controller extracts the 3-parameter vector part
//   of the error quaternion q_err = q_des* ⊗ q_cur, exactly as the M13C MEKF
//   residual does, and applies PD feedback on that rotation vector plus rate.
//
// Frames:
//   q_current / q_desired: body-to-ECI attitudes (v_eci = q (*) v_body).
//   omega_current/desired: body rates, BODY components (rad/s).
//   Output torque: BODY components (N*m). Gains are 3-vectors for per-axis
//   tuning on a diagonal-inertia bus.
//
// Units (SI): angle (rad, via q_err vector part ≈ half-angle axis for small
//   errors), rate (rad/s), torque (N*m). Gains: Kp (N*m/rad), Kd (N*m*s/rad).
//
// Governing equations:
//   q_err = q_des.conjugate() * q_cur  (rotation taking desired body axes to
//           current body axes, expressed in the desired frame)
//   e_att = shortest-arc vector part of q_err (double-cover sign aligned)
//   e_rate = omega_cur - omega_des      (both in BODY components)
//   tau_des = -Kp .* e_att - Kd .* e_rate
// Rate damping (M16B) is the special case q_des = q_cur (e_att = 0):
//   tau = -K_w (omega - omega_cmd).
//
// Small-angle reading: q_err ≈ [1, e/2], so e_att ≈ (angle/2) axis. The 1/2
// scale is absorbed into Kp (documented, not hidden): quoting Kp without noting
// the half-angle convention would mislead gain comparisons with literature that
// feeds back the full angle. Stability intuition: for one principal axis with
// inertia I, closed loop ≈ I θ_dd + Kd θ_d + (Kp/2) θ = 0, stable for Kp,Kd > 0
// with damping ratio ζ = Kd / sqrt(2 Kp I).
//
// Assumptions and non-goals:
//   - Rigid bus, diagonal inertia, torque exerted exactly (saturation is applied
//     downstream by M16E through the M14 actuator models, never here).
//   - Gains constant and positive; no integral action (bias torques like gravity
//     gradient are an M17+ disturbance-rejection concern), no gain scheduling.
//   - No state estimation here: inputs are estimates by contract (M17 wires them).

#include "math/quaternion.hpp"
#include "math/vector3.hpp"

#include <cmath>
#include <stdexcept>

namespace astradock::control {

namespace detail {

inline void require_finite_attitude_inputs(
    const math::Quaternion& q_current,
    const math::Quaternion& q_desired,
    const math::Vector3& omega_current,
    const math::Vector3& omega_desired) {
    if (!math::is_finite(q_current) || !math::is_finite(q_desired)) {
        throw std::domain_error("Attitude quaternions must contain only finite values");
    }
    if (!q_current.is_unit(1.0e-3) || !q_desired.is_unit(1.0e-3)) {
        throw std::domain_error("Attitude quaternions must be approximately unit quaternions");
    }
    if (!math::is_finite(omega_current) || !math::is_finite(omega_desired)) {
        throw std::domain_error("Angular velocities must contain only finite values");
    }
}

inline void require_positive_gains(const math::Vector3& kp, const math::Vector3& kd) {
    if (!math::is_finite(kp) || !math::is_finite(kd)) {
        throw std::domain_error("Controller gains must contain only finite values");
    }
    if (kp.x() <= 0.0 || kp.y() <= 0.0 || kp.z() <= 0.0 || kd.x() <= 0.0 || kd.y() <= 0.0
        || kd.z() <= 0.0) {
        throw std::domain_error("PD controller gains must be strictly positive per axis");
    }
}

}  // namespace detail

// Gains for a 3-axis decoupled PD law. Per-axis so a (10, 20, 30) inertia bus
// can carry axis-matched damping without one shared compromise gain.
struct AttitudePdGains {
    math::Vector3 kp_Nm_per_rad{1.0, 1.0, 1.0};
    math::Vector3 kd_Nm_s_per_rad{1.0, 1.0, 1.0};
};

// Shortest-arc attitude error vector: vector part of q_err = q_des* ⊗ q_cur
// with double-cover sign alignment (q_err.w >= 0). Matches the M13C star-tracker
// residual convention. Output ≈ (angle/2) * axis for small errors, BODY-frame
// relevant components.
[[nodiscard]] inline math::Vector3 attitude_error_vector(
    const math::Quaternion& q_current_body_to_eci,
    const math::Quaternion& q_desired_body_to_eci) {
    if (!math::is_finite(q_current_body_to_eci) || !math::is_finite(q_desired_body_to_eci)) {
        throw std::domain_error("Attitude quaternions must contain only finite values");
    }
    if (!q_current_body_to_eci.is_unit(1.0e-3) || !q_desired_body_to_eci.is_unit(1.0e-3)) {
        throw std::domain_error("Attitude quaternions must be approximately unit quaternions");
    }
    math::Quaternion q_err = q_desired_body_to_eci.conjugate() * q_current_body_to_eci;
    if (q_err.w() < 0.0) {
        q_err = math::Quaternion(-q_err.w(), -q_err.x(), -q_err.y(), -q_err.z());
    }
    return q_err.vector_part();
}

// Geodesic attitude error angle: theta = 2 acos(|q_err.w|), range [0, pi].
[[nodiscard]] inline double attitude_error_angle_rad(
    const math::Quaternion& q_current_body_to_eci,
    const math::Quaternion& q_desired_body_to_eci) {
    if (!math::is_finite(q_current_body_to_eci) || !math::is_finite(q_desired_body_to_eci)) {
        throw std::domain_error("Attitude quaternions must contain only finite values");
    }
    if (!q_current_body_to_eci.is_unit(1.0e-3) || !q_desired_body_to_eci.is_unit(1.0e-3)) {
        throw std::domain_error("Attitude quaternions must be approximately unit quaternions");
    }
    const math::Quaternion q_err = q_desired_body_to_eci.conjugate() * q_current_body_to_eci;
    const double clamped = std::clamp(std::abs(q_err.w()), 0.0, 1.0);
    return 2.0 * std::acos(clamped);
}

// Full quaternion-error PD law: tau_des_body = -Kp .* e_att - Kd .* e_rate.
[[nodiscard]] inline math::Vector3 attitude_pd_torque_body_Nm(
    const math::Quaternion& q_current_body_to_eci,
    const math::Quaternion& q_desired_body_to_eci,
    const math::Vector3& omega_current_body_rad_s,
    const math::Vector3& omega_desired_body_rad_s,
    const AttitudePdGains& gains) {
    detail::require_finite_attitude_inputs(
        q_current_body_to_eci, q_desired_body_to_eci, omega_current_body_rad_s, omega_desired_body_rad_s);
    detail::require_positive_gains(gains.kp_Nm_per_rad, gains.kd_Nm_s_per_rad);
    const math::Vector3 e_att =
        attitude_error_vector(q_current_body_to_eci, q_desired_body_to_eci);
    const math::Vector3 e_rate = omega_current_body_rad_s - omega_desired_body_rad_s;
    return {
        -gains.kp_Nm_per_rad.x() * e_att.x() - gains.kd_Nm_s_per_rad.x() * e_rate.x(),
        -gains.kp_Nm_per_rad.y() * e_att.y() - gains.kd_Nm_s_per_rad.y() * e_rate.y(),
        -gains.kp_Nm_per_rad.z() * e_att.z() - gains.kd_Nm_s_per_rad.z() * e_rate.z(),
    };
}

// Rate damping law (M16B): tau = -K_w .* (omega - omega_cmd). Pure derivative
// action on body rates; the attitude-independent special case of the PD law.
[[nodiscard]] inline math::Vector3 rate_damping_torque_body_Nm(
    const math::Vector3& omega_current_body_rad_s,
    const math::Vector3& omega_command_body_rad_s,
    const math::Vector3& k_w_Nm_s_per_rad) {
    if (!math::is_finite(omega_current_body_rad_s) || !math::is_finite(omega_command_body_rad_s)
        || !math::is_finite(k_w_Nm_s_per_rad)) {
        throw std::domain_error("Rate damping inputs must contain only finite values");
    }
    if (k_w_Nm_s_per_rad.x() <= 0.0 || k_w_Nm_s_per_rad.y() <= 0.0 || k_w_Nm_s_per_rad.z() <= 0.0) {
        throw std::domain_error("Rate damping gains must be strictly positive per axis");
    }
    const math::Vector3 e_rate = omega_current_body_rad_s - omega_command_body_rad_s;
    return {
        -k_w_Nm_s_per_rad.x() * e_rate.x(),
        -k_w_Nm_s_per_rad.y() * e_rate.y(),
        -k_w_Nm_s_per_rad.z() * e_rate.z(),
    };
}

// Suggests per-axis PD gains from inertia + desired damping ratio and settle
// time, so tuning has an engineering basis instead of a magic number. Uses the
// small-angle closed loop I θ_dd + Kd θ_d + (Kp/2) θ = 0 with wn = 4/(ζ ts):
//   Kp = 2 I wn^2,  Kd = 2 ζ wn I.
struct PdTuningRequest {
    math::Vector3 inertia_kg_m2{1.0, 1.0, 1.0};
    double damping_ratio{0.9};
    double settle_time_s{30.0};
};

[[nodiscard]] inline AttitudePdGains suggest_pd_gains(const PdTuningRequest& request) {
    if (!math::is_finite(request.inertia_kg_m2)
        || request.inertia_kg_m2.x() <= 0.0 || request.inertia_kg_m2.y() <= 0.0
        || request.inertia_kg_m2.z() <= 0.0) {
        throw std::domain_error("Tuning inertias must be finite and strictly positive");
    }
    if (!std::isfinite(request.damping_ratio) || request.damping_ratio <= 0.0) {
        throw std::domain_error("Damping ratio must be finite and strictly positive");
    }
    if (!std::isfinite(request.settle_time_s) || request.settle_time_s <= 0.0) {
        throw std::domain_error("Settle time must be finite and strictly positive");
    }
    const double wn = 4.0 / (request.damping_ratio * request.settle_time_s);
    const auto axis_gain = [&](double inertia) {
        return 2.0 * inertia * wn * wn;
    };
    const auto axis_damp = [&](double inertia) {
        return 2.0 * request.damping_ratio * wn * inertia;
    };
    return {
        {axis_gain(request.inertia_kg_m2.x()),
         axis_gain(request.inertia_kg_m2.y()),
         axis_gain(request.inertia_kg_m2.z())},
        {axis_damp(request.inertia_kg_m2.x()),
         axis_damp(request.inertia_kg_m2.y()),
         axis_damp(request.inertia_kg_m2.z())},
    };
}

}  // namespace astradock::control
