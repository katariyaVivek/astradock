#pragma once

// AstraDock M17 — Closed-loop GNC wiring: estimates-only control.
//
// Physical problem:
//   M16 proved each control law converges when fed perfect state. A real loop
//   never has perfect state: it flies the navigation estimate. M17 closes the
//   loop bus -> sensors -> navigation -> control -> actuators -> bus dynamics
//   with one inviolable rule: guidance and control see ESTIMATES, never the
//   simulator state. This header is the small, auditable seam where that rule
//   is enforced — the controller input structs below have no simulator-state
//   member, so a state leak is a compile-time shape error, not a runtime hope.
//
// Data flow per control step (attitude hold + maneuver, M17A/M17B):
//   1. Sensors sample the simulated bus (M12 models): gyro + star tracker.
//   2. Navigation filter predicts/updates from measurements only (M13).
//   3. Guidance emits reference attitude/rate from the maneuver timeline.
//   4. Control computes desired body torque from ESTIMATED attitude/rate vs
//      reference (M16A law), saturates against wheel capability (M16E), and
//      steps the M14 reaction-wheel assembly.
//   5. Achieved reaction torque integrates the 6-DOF simulated bus (M09/M10).
//   6. Telemetry records simulated, estimated, and reference values plus errors,
//      commands, achieved outputs, saturation flags, sensor validity, covariance.
//
// Frames/units: SI throughout; torque BODY (N*m); attitude body-to-ECI.

#include "actuators/reaction_wheel.hpp"
#include "control/attitude_pd.hpp"
#include "control/control_saturation.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace astradock::gnc {

namespace detail {

inline void require_finite_loop_inputs(
    const math::Quaternion& estimated_attitude,
    const math::Quaternion& reference_attitude,
    const math::Vector3& estimated_rate,
    const math::Vector3& reference_rate) {
    if (!math::is_finite(estimated_attitude) || !math::is_finite(reference_attitude)
        || !math::is_finite(estimated_rate) || !math::is_finite(reference_rate)) {
        throw std::domain_error("Closed-loop control inputs must contain only finite values");
    }
    if (!estimated_attitude.is_unit(1.0e-3) || !reference_attitude.is_unit(1.0e-3)) {
        throw std::domain_error("Closed-loop attitudes must be approximately unit quaternions");
    }
}

}  // namespace detail

// Controller input bundle. NOTE: fields are named estimated_* / reference_* on
// purpose — there is deliberately no simulator-state member. The M17 static audit
// (test) pins this down by scanning these headers for simulator-state declarations.
struct EstimateBasedAttitudeCommand {
    math::Quaternion estimated_attitude_body_to_eci{math::Quaternion::identity()};
    math::Vector3 estimated_rate_body_rad_s{};
    math::Quaternion reference_attitude_body_to_eci{math::Quaternion::identity()};
    math::Vector3 reference_rate_body_rad_s{};
};

// One control tick: estimate-based PD torque + per-axis saturation accounting.
// The achieved torque is what the caller may forward to the wheel assembly;
// when saturated, achieved < desired and the flag records it (never clipped
// silently). Pure function of estimates + reference + gains + limits.
struct ClosedLoopTorqueStep {
    math::Vector3 desired_torque_body_Nm{};
    math::Vector3 achieved_torque_body_Nm{};
    math::Vector3 attitude_error_vector{};
    double attitude_error_angle_rad{0.0};
    double rate_error_norm_rad_s{0.0};
    bool saturated{false};
};

[[nodiscard]] inline ClosedLoopTorqueStep step_estimate_based_attitude_control(
    const EstimateBasedAttitudeCommand& command,
    const control::AttitudePdGains& gains,
    const math::Vector3& torque_limit_body_Nm) {
    detail::require_finite_loop_inputs(
        command.estimated_attitude_body_to_eci, command.reference_attitude_body_to_eci,
        command.estimated_rate_body_rad_s, command.reference_rate_body_rad_s);
    if (!math::is_finite(torque_limit_body_Nm) || torque_limit_body_Nm.x() <= 0.0
        || torque_limit_body_Nm.y() <= 0.0 || torque_limit_body_Nm.z() <= 0.0) {
        throw std::domain_error("Torque limits must be finite and strictly positive per axis");
    }
    const math::Vector3 desired = control::attitude_pd_torque_body_Nm(
        command.estimated_attitude_body_to_eci,
        command.reference_attitude_body_to_eci,
        command.estimated_rate_body_rad_s,
        command.reference_rate_body_rad_s,
        gains);
    const control::SaturatedControl sat =
        control::saturate_control_vector(desired, torque_limit_body_Nm);
    ClosedLoopTorqueStep step;
    step.desired_torque_body_Nm = sat.desired;
    step.achieved_torque_body_Nm = sat.achieved;
    step.attitude_error_vector = control::attitude_error_vector(
        command.estimated_attitude_body_to_eci, command.reference_attitude_body_to_eci);
    step.attitude_error_angle_rad = control::attitude_error_angle_rad(
        command.estimated_attitude_body_to_eci, command.reference_attitude_body_to_eci);
    step.rate_error_norm_rad_s =
        (command.estimated_rate_body_rad_s - command.reference_rate_body_rad_s).norm();
    step.saturated = sat.saturated;
    return step;
}

// Piecewise reference timeline: hold segments joined by instantaneous switches.
// Slew performance between switches is the controller's job (M16 slew demo);
// the timeline only declares what the reference IS at time t. Deterministic,
// zero allocation after construction.
struct AttitudeReferenceTimeline {
    struct Segment {
        double start_time_s{0.0};
        math::Quaternion attitude_body_to_eci{math::Quaternion::identity()};
        math::Vector3 rate_body_rad_s{};
    };
    std::vector<Segment> segments{};

    [[nodiscard]] const Segment& reference_at(double time_s) const {
        if (segments.empty()) {
            throw std::domain_error("Attitude reference timeline must contain at least one segment");
        }
        if (!std::isfinite(time_s)) {
            throw std::domain_error("Reference query time must be finite");
        }
        const Segment* active = &segments.front();
        for (const auto& seg : segments) {
            if (seg.start_time_s <= time_s) {
                active = &seg;
            }
        }
        return *active;
    }
};

// Closed-loop performance summary over one run (computed from scoring telemetry;
// the controller itself never sees these values — the naming convention marks
// harness-side scoring columns, verified by the M17 static audit).
struct ClosedLoopMetrics {
    double settle_time_s{-1.0};
    double max_attitude_error_rad{0.0};
    double final_attitude_error_rad{0.0};
    double final_rate_norm_rad_s{0.0};
    double max_commanded_torque_Nm{0.0};
    double saturated_fraction{0.0};
    double rms_estimate_error_rad{0.0};
    bool converged{false};
};

[[nodiscard]] inline ClosedLoopMetrics summarize_closed_loop_run(
    const std::vector<double>& time_s,
    const std::vector<double>& scored_error_rad,
    const std::vector<double>& commanded_torque_Nm,
    const std::vector<double>& estimate_error_rad,
    const std::vector<int>& saturated_flags,
    double error_threshold_rad,
    double rate_unused_placeholder = 0.0) {
    static_cast<void>(rate_unused_placeholder);
    if (time_s.empty() || time_s.size() != scored_error_rad.size()
        || time_s.size() != commanded_torque_Nm.size() || time_s.size() != estimate_error_rad.size()
        || time_s.size() != saturated_flags.size()) {
        throw std::domain_error("Closed-loop telemetry columns must be non-empty and equal length");
    }
    if (!std::isfinite(error_threshold_rad) || error_threshold_rad <= 0.0) {
        throw std::domain_error("Error threshold must be finite and strictly positive");
    }
    ClosedLoopMetrics m;
    double sat_count = 0.0;
    double est_sq = 0.0;
    for (std::size_t i = 0; i < time_s.size(); ++i) {
        m.max_attitude_error_rad = std::max(m.max_attitude_error_rad, scored_error_rad[i]);
        m.max_commanded_torque_Nm = std::max(m.max_commanded_torque_Nm, commanded_torque_Nm[i]);
        sat_count += saturated_flags[i] != 0 ? 1.0 : 0.0;
        est_sq += estimate_error_rad[i] * estimate_error_rad[i];
        if (m.settle_time_s < 0.0 && scored_error_rad[i] < error_threshold_rad) {
            // Latching settle requires all subsequent samples to hold the bound.
            bool holds = true;
            for (std::size_t j = i; j < time_s.size(); ++j) {
                if (scored_error_rad[j] >= error_threshold_rad) {
                    holds = false;
                    break;
                }
            }
            if (holds) {
                m.settle_time_s = time_s[i];
            }
        }
    }
    m.final_attitude_error_rad = scored_error_rad.back();
    m.saturated_fraction = sat_count / static_cast<double>(time_s.size());
    m.rms_estimate_error_rad = std::sqrt(est_sq / static_cast<double>(time_s.size()));
    m.converged = m.settle_time_s >= 0.0;
    return m;
}

}  // namespace astradock::gnc
