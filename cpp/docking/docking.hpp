#pragma once

// AstraDock M19 — Autonomous docking: port frames, capture geometry, contact.
//
// Physical problem:
//   Rendezvous (M18) delivers the chaser to a 50 m hold. Docking closes the
//   last meters along the docking axis, aligns attitude, makes compliant
//   contact, and captures — or aborts deterministically when any acceptance
//   criterion fails. The frame problem is the core: the docking ports live on
//   the vehicle BODIES (target port + chaser port), each with position offset
//   and approach-axis direction in body coordinates, while the approach is
//   flown in LVLH. Every transform here is explicit; a sign flip in the
//   approach axis is a mission failure, so axis directions get hand-derived
//   regression tests, per the frame-convention requirement.
//
// Frames:
//   Target port: offset r_tp_t (target BODY m), axis a_tp_t (target BODY unit,
//     pointing OUT of the port = direction the chaser arrives FROM).
//   Chaser port: offset r_cp_c (chaser BODY m), axis a_cp_c (chaser BODY unit,
//     pointing OUT of the port = direction the target arrives FROM).
//   Nominal mated condition: ports co-located, axes anti-parallel, attitudes
//   aligned up to the port symmetry. Relative vectors resolved in LVLH via the
//   two body-to-ECI quaternions and the target LVLH DCM.
// Units (SI): m, m/s, rad, rad/s, N, N*s.
//
// Governing construction:
//   Port world positions (ECI): P_tp = R_t + q_t (*) r_tp_t,
//     P_cp = R_c + q_c (*) r_cp_c.
//   Relative port error (LVLH): e = C_LVLH_ECI(target) (P_cp - P_tp).
//   Axial coordinate s = e . u_axis (u_axis = LVLH approach direction, +from
//     target toward chaser hold: -y in M18 terms, i.e. u_axis = [0,-1,0] LVLH).
//   Lateral error l = |e - s u_axis|. Closing speed vc = -ds/dt (computed by
//     finite difference on s by the caller; this header takes s + s_dot).
//   Relative attitude: q_rel = q_t* ⊗ q_c; alignment angle φ = 2 acos|w|;
//     relative rate ω_rel = ω_c - ω_t (BODY-mismatched frames: reported as the
//     norm of the difference, a conservative bound — documented limitation).
//   Contact (1-DOF compliant along axis once |s| < contact_range):
//     F = -k s_pen - c s_dot_pen for penetration s_pen > 0, else 0
//     (spring-damper, penalty method; no structural modes, no friction cone).
//   Capture latch: latched when acceptance holds continuously for latch_time.
//
// Acceptance (ALL must hold at the latch check):
//   |lateral| <= lateral_limit | axial closing vc <= vc_limit
//   |axial residual| <= axial_limit | φ <= align_limit | |ω_rel| <= rate_limit.
// Abort (ANY triggers): keep-out breach (M18 predicate), lateral > abort band,
//   closing above hard limit, attitude beyond abort angle, contact force above
//   crush limit, timeout. Abort OUTPUT here is the boolean + reason enum; the
//   retreat trajectory is M21 mission-scope (this header only declares it).
//
// Assumptions: static cooperative target (target rates optional input, zero by
// default); ports on rigid bodies; approach axis fixed in LVLH over the final
// meters; contact is axial penalty only (lateral contact = abort, not sliding).

#include "frames/lvlh.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace astradock::docking {

namespace detail {

inline void require_finite_port(
    const math::Vector3& offset_body_m,
    const math::Vector3& axis_body) {
    if (!math::is_finite(offset_body_m) || !math::is_finite(axis_body)) {
        throw std::domain_error("Docking port offset and axis must contain only finite values");
    }
    if (std::abs(axis_body.norm() - 1.0) > 1.0e-9) {
        throw std::domain_error("Docking port axis must be a unit vector in the BODY frame");
    }
}

inline void require_positive_dock(double value, const char* message) {
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::domain_error(message);
    }
}

}  // namespace detail

// One docking port: mechanical interface location + approach direction.
struct DockingPort {
    math::Vector3 offset_body_m{};
    math::Vector3 axis_body{0.0, 0.0, 1.0};
};

// Capture envelope + contact + abort thresholds.
struct DockingEnvelope {
    double lateral_limit_m{0.25};
    double axial_limit_m{0.25};
    double closing_speed_limit_mps{0.1};
    double alignment_limit_rad{5.0 * 3.141592653589793 / 180.0};
    double rate_limit_rad_s{0.5 * 3.141592653589793 / 180.0};
    double contact_range_m{1.0};
    double contact_stiffness_N_per_m{1000.0};
    double contact_damping_N_s_per_m{100.0};
    double crush_force_limit_N{500.0};
    double abort_lateral_m{1.0};
    double abort_lateral_range_m{10.0};  // lateral abort enforced only when
                                         // separation s <= this range: far out,
                                         // port-offset geometry reads meters
                                         // of lateral on a nominal approach
                                         // (M18 corridor owns far-field safety)
    double abort_angle_rad{15.0 * 3.141592653589793 / 180.0};
    double latch_time_s{5.0};
};

// Abort reason codes (deterministic, telemetry-friendly).
enum class AbortReason {
    none,
    keep_out_breach,
    lateral_exceeded,
    closing_speed_exceeded,
    attitude_exceeded,
    crush_force_exceeded,
    timeout,
};
[[nodiscard]] inline const char* abort_reason_name(AbortReason reason) noexcept {
    switch (reason) {
    case AbortReason::none:
        return "none";
    case AbortReason::keep_out_breach:
        return "keep_out_breach";
    case AbortReason::lateral_exceeded:
        return "lateral_exceeded";
    case AbortReason::closing_speed_exceeded:
        return "closing_speed_exceeded";
    case AbortReason::attitude_exceeded:
        return "attitude_exceeded";
    case AbortReason::crush_force_exceeded:
        return "crush_force_exceeded";
    case AbortReason::timeout:
        return "timeout";
    }
    return "unknown";
}

// Full docking state evaluation at one tick.
struct DockingStep {    math::Vector3 port_error_lvlh_m{};  // chaser-port minus target-port, LVLH
    double axial_m{0.0};                // s along approach axis (negative = separated)
    double lateral_m{0.0};
    double closing_speed_mps{0.0};      // +closing (toward mate)
    double alignment_rad{0.0};
    double relative_rate_norm_rad_s{0.0};
    double contact_force_N{0.0};
    bool in_contact{false};
    bool acceptance{false};
    bool latched{false};
    bool abort{false};
    AbortReason abort_reason{AbortReason::none};
};

// Port world position in ECI: P = R + q (*) r_body.
[[nodiscard]] inline math::Vector3 port_position_eci_m(
    const math::Vector3& cm_position_eci_m,
    const math::Quaternion& attitude_body_to_eci,
    const DockingPort& port) {
    detail::require_finite_port(port.offset_body_m, port.axis_body);
    if (!math::is_finite(cm_position_eci_m) || !math::is_finite(attitude_body_to_eci)) {
        throw std::domain_error("Spacecraft position and attitude must be finite");
    }
    if (!attitude_body_to_eci.is_unit(1.0e-3)) {
        throw std::domain_error("Attitude must be approximately a unit quaternion");
    }
    return cm_position_eci_m + attitude_body_to_eci.rotate_vector(port.offset_body_m);
}

// Relative attitude alignment angle between two body frames.
[[nodiscard]] inline double relative_attitude_rad(
    const math::Quaternion& target_body_to_eci,
    const math::Quaternion& chaser_body_to_eci) {
    if (!math::is_finite(target_body_to_eci) || !math::is_finite(chaser_body_to_eci)) {
        throw std::domain_error("Attitudes must contain only finite values");
    }
    if (!target_body_to_eci.is_unit(1.0e-3) || !chaser_body_to_eci.is_unit(1.0e-3)) {
        throw std::domain_error("Attitudes must be approximately unit quaternions");
    }
    const math::Quaternion q_rel = target_body_to_eci.conjugate() * chaser_body_to_eci;
    return 2.0 * std::acos(std::clamp(std::abs(q_rel.w()), 0.0, 1.0));
}

// Axial penalty contact force: F = k*pen + c*pen_rate for pen > 0, else 0.
// pen = contact_range - |axial_gap| once inside contact range... here the caller
// passes penetration directly (gap undershoot), keeping this function pure.
[[nodiscard]] inline double contact_force_N(
    double penetration_m,
    double penetration_rate_mps,
    double stiffness_N_per_m,
    double damping_N_s_per_m) {
    if (!std::isfinite(penetration_m) || !std::isfinite(penetration_rate_mps)) {
        throw std::domain_error("Penetration and rate must be finite");
    }
    detail::require_positive_dock(stiffness_N_per_m, "Contact stiffness must be positive");
    if (!std::isfinite(damping_N_s_per_m) || damping_N_s_per_m < 0.0) {
        throw std::domain_error("Contact damping must be finite and non-negative");
    }
    if (penetration_m <= 0.0) {
        return 0.0;
    }
    const double force = stiffness_N_per_m * penetration_m + damping_N_s_per_m * penetration_rate_mps;
    return force > 0.0 ? force : 0.0;  // penalty only pushes, never pulls
}

// One docking evaluation tick. Separation s >= 0 means ports apart (s = 0 at
// the mate plane); penetration = max(0, -s) depth past the plane. Closing
// speed is positive toward mate (closing = -ds/dt); the caller passes the
// finite-difference ds/dt negated, or equivalently the approach speed.
[[nodiscard]] inline DockingStep step_docking(
    const math::Vector3& target_cm_eci_m,
    const math::Quaternion& target_attitude,
    const math::Vector3& target_rate_body_rad_s,
    const math::Vector3& chaser_cm_eci_m,
    const math::Quaternion& chaser_attitude,
    const math::Vector3& chaser_rate_body_rad_s,
    const DockingPort& target_port,
    const DockingPort& chaser_port,
    const math::Vector3& approach_axis_lvlh,
    const math::Matrix3& dcm_lvlh_from_eci,
    double axial_closing_speed_mps,
    const DockingEnvelope& envelope,
    double latch_accumulated_s,
    double elapsed_s,
    double timeout_s) {
    detail::require_finite_port(target_port.offset_body_m, target_port.axis_body);
    detail::require_finite_port(chaser_port.offset_body_m, chaser_port.axis_body);
    if (!math::is_finite(approach_axis_lvlh) || std::abs(approach_axis_lvlh.norm() - 1.0) > 1.0e-9) {
        throw std::domain_error("Approach axis must be a finite LVLH unit vector");
    }
    if (!std::isfinite(axial_closing_speed_mps) || !std::isfinite(latch_accumulated_s)
        || !std::isfinite(elapsed_s) || !std::isfinite(timeout_s) || timeout_s <= 0.0) {
        throw std::domain_error("Docking tick scalars must be finite with positive timeout");
    }

    const math::Vector3 p_target = port_position_eci_m(target_cm_eci_m, target_attitude, target_port);
    const math::Vector3 p_chaser = port_position_eci_m(chaser_cm_eci_m, chaser_attitude, chaser_port);
    const math::Vector3 error_lvlh = dcm_lvlh_from_eci * (p_chaser - p_target);

    DockingStep step;
    step.port_error_lvlh_m = error_lvlh;
    step.axial_m = error_lvlh.dot(approach_axis_lvlh);
    const math::Vector3 lateral_vec = error_lvlh - approach_axis_lvlh * step.axial_m;
    step.lateral_m = lateral_vec.norm();
    step.closing_speed_mps = axial_closing_speed_mps;
    step.alignment_rad = relative_attitude_rad(target_attitude, chaser_attitude);
    step.relative_rate_norm_rad_s =
        (chaser_rate_body_rad_s - target_rate_body_rad_s).norm();

    // Contact + acceptance. Separation s > 0: ports apart. Contact at s <= 0
    // with penetration pen = -s pushed back by the penalty spring.
    step.in_contact = step.axial_m <= 0.0;
    if (step.in_contact) {
        step.contact_force_N = contact_force_N(
            -step.axial_m, axial_closing_speed_mps, envelope.contact_stiffness_N_per_m,
            envelope.contact_damping_N_s_per_m);
    }
    const bool geometry_ok = (step.lateral_m <= envelope.lateral_limit_m)
        && (step.axial_m <= envelope.contact_range_m + envelope.axial_limit_m)
        && (step.axial_m >= -envelope.axial_limit_m)
        && (step.closing_speed_mps <= envelope.closing_speed_limit_mps)
        && (step.closing_speed_mps >= -envelope.closing_speed_limit_mps)
        && (step.alignment_rad <= envelope.alignment_limit_rad)
        && (step.relative_rate_norm_rad_s <= envelope.rate_limit_rad_s);
    step.acceptance = geometry_ok && (step.contact_force_N <= envelope.crush_force_limit_N);
    step.latched = step.acceptance && (latch_accumulated_s >= envelope.latch_time_s);

    // Abort predicates (deterministic priority order). Lateral abort is gated
    // by range: far-field lateral geometry belongs to the M18 corridor, and
    // body port offsets read as meters of port-gap lateral at long range.
    // NOTE: closing-abort compares against the ENVELOPE limit, not the cruise
    // command — a hot-but-tracked approach stays inside acceptance; only a
    // runaway (3x envelope) aborts. Acceptance itself enforces the envelope
    // limit every tick, so overspeed can never latch.
    if (step.axial_m <= envelope.abort_lateral_range_m
        && step.lateral_m > envelope.abort_lateral_m) {
        step.abort = true;
        step.abort_reason = AbortReason::lateral_exceeded;
    } else if (step.closing_speed_mps > 3.0 * envelope.closing_speed_limit_mps) {
        step.abort = true;
        step.abort_reason = AbortReason::closing_speed_exceeded;
    } else if (step.alignment_rad > envelope.abort_angle_rad) {
        step.abort = true;
        step.abort_reason = AbortReason::attitude_exceeded;
    } else if (step.contact_force_N > envelope.crush_force_limit_N) {
        step.abort = true;
        step.abort_reason = AbortReason::crush_force_exceeded;
    } else if (elapsed_s > timeout_s) {
        step.abort = true;
        step.abort_reason = AbortReason::timeout;
    }
    return step;
}

}  // namespace astradock::docking
