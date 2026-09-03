#pragma once

// AstraDock M18 — Rendezvous & proximity-operations guidance.
//
// Physical problem:
//   Bring the chaser from a far hold point to a close hold point along the
//   approach axis without violating safety geometry. The reference is a
//   straight LVLH segment from waypoint to waypoint flown at a bounded closing
//   speed; the M16 translation law tracks it and the M14 thrusters execute it.
//   Safety is geometric and checked every tick: keep-out sphere (hard abort),
//   approach corridor (lateral bound around the axis), per-leg closing-speed
//   limit. No trajectory optimization, no Lambert targeting — a transparent
//   waypoint sequencer whose every decision can be hand-checked.
//
// Frames:
//   Waypoints are LVLH positions in the TARGET-centered LVLH frame (M15
//   convention: x radial, y along-track, z orbit-normal). The approach axis is
//   -y (from behind/below toward the target at the origin): waypoints march
//   from far (-5000 m y) to near (-50 m y). State input is the estimated LVLH
//   relative state; control output is an LVLH acceleration command.
//
// Units (SI): m, m/s, m/s^2, s.
//
// Governing construction (per active leg, waypoints A -> B):
//   e = B - rho (position error), d = |e|, u = e/d (line-of-sight unit).
//   Closing speed profile: v_close = min(v_max_leg, sqrt(2 a_brake d))
//     (constant-speed cruise far out, sqrt braking near the waypoint with
//     brake deceleration a_brake — the standard constant-deceleration profile).
//   v_des = v_close * u; a_des = (v_des - v_rel) / tau_track (first-order
//     velocity tracking with time constant tau_track), then the M16 relative
//     PD law adds CW feedforward + error feedback around it.
// Leg switches when d < capture_radius_m. Hysteresis: none needed — legs only
// advance (monotonic mission progress; retreat is an ABORT, M19 scope).
//
// Safety predicates (evaluated on the ESTIMATED state, flagged in telemetry):
//   keep-out: |rho| < keep_out_radius_m -> VIOLATION (abort).
//   corridor: lateral distance from the approach axis > corridor_radius_m while
//     inside corridor_range_m of the target -> VIOLATION (abort).
//   closing speed: radial-inward rate toward the target > v_limit -> WARNING
//     (telemetry flag; the profile itself commands within limits, so a warning
//     means tracking error or a fault — M20 investigates).
//
// Assumptions and non-goals:
//   - Waypoints static in LVLH (valid for short rendezvous arcs; CW drift is
//     what the tracker fights — measured, not assumed away).
//   - Single chaser, cooperative static target at the LVLH origin.
//   - No plume, no obstacle other than the keep-out sphere, no abort
//     trajectory generation (M19). Estimator dynamics are the caller's job.

#include "control/relative_pd.hpp"
#include "math/vector3.hpp"
#include "relative/relative_state.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace astradock::rendezvous {

namespace detail {

inline void require_finite_waypoint(const math::Vector3& waypoint_lvlh_m) {
    if (!math::is_finite(waypoint_lvlh_m)) {
        throw std::domain_error("Rendezvous waypoint must contain only finite values");
    }
}

inline void require_positive_guidance(double value, const char* message) {
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::domain_error(message);
    }
}

}  // namespace detail

// One hold point: a commanded LVLH offset. The leg INTO the hold point carries
// the closing-speed limit; capture radius declares arrival.
struct HoldPoint {
    math::Vector3 position_lvlh_m{};
    double closing_speed_limit_mps{1.0};
    double capture_radius_m{10.0};
};

// Safety geometry around the target (LVLH origin).
struct SafetyCorridor {
    double keep_out_radius_m{25.0};
    double corridor_radius_m{50.0};
    double corridor_range_m{1000.0};
    double closing_speed_hard_limit_mps{3.0};
};

// Default far-to-near approach sequence along -y: 5000 m -> 1000 m ->
// 250 m -> 50 m, each leg slower than the last (engineering basis: braking
// capability ~1 mm/s^2 from M16 station-keeping demo => v^2/(2a) stopping
// distance sets leg speeds: 2 m/s needs 2 km, 1 m/s needs 500 m, etc.).
[[nodiscard]] inline std::vector<HoldPoint> default_approach_sequence() {
    return {
        {math::Vector3{0.0, -5000.0, 0.0}, 2.0, 50.0},
        {math::Vector3{0.0, -1000.0, 0.0}, 1.0, 25.0},
        {math::Vector3{0.0, -250.0, 0.0}, 0.5, 10.0},
        {math::Vector3{0.0, -50.0, 0.0}, 0.2, 5.0},
    };
}

// Guidance tick output: reference + command + safety flags. Commanded
// acceleration is LVLH (caller allocates to body thrusters via live attitude).
struct GuidanceStep {
    std::size_t active_leg{0};
    math::Vector3 reference_position_lvlh_m{};
    math::Vector3 reference_velocity_lvlh_mps{};
    math::Vector3 commanded_accel_lvlh_mps2{};
    double distance_to_waypoint_m{0.0};
    double range_to_target_m{0.0};
    double lateral_error_m{0.0};
    double closing_speed_mps{0.0};
    bool leg_complete{false};
    bool keep_out_violation{false};
    bool corridor_violation{false};
    bool closing_speed_warning{false};
    bool mission_complete{false};
};

// Lateral distance from the approach (-y) axis: sqrt(x^2 + z^2).
[[nodiscard]] inline double lateral_axis_error_m(const math::Vector3& rho_lvlh_m) {
    if (!math::is_finite(rho_lvlh_m)) {
        throw std::domain_error("Relative position must be finite for corridor check");
    }
    return std::sqrt(rho_lvlh_m.x() * rho_lvlh_m.x() + rho_lvlh_m.z() * rho_lvlh_m.z());
}

// Inward closing speed toward the target: -d|r|/dt = -(rho.v)/|rho|.
[[nodiscard]] inline double closing_speed_mps(
    const math::Vector3& rho_lvlh_m,
    const math::Vector3& v_rel_lvlh_mps) {
    if (!math::is_finite(rho_lvlh_m) || !math::is_finite(v_rel_lvlh_mps)) {
        throw std::domain_error("Relative state must be finite for closing-speed check");
    }
    const double range = rho_lvlh_m.norm();
    if (range == 0.0) {
        return 0.0;
    }
    return -(rho_lvlh_m.dot(v_rel_lvlh_mps)) / range;
}

// One guidance tick: track the active leg, advance on capture, check safety.
[[nodiscard]] inline GuidanceStep step_rendezvous_guidance(
    const relative::RelativeStateLvlh& estimated_state,
    const std::vector<HoldPoint>& sequence,
    std::size_t& active_leg,
    const SafetyCorridor& safety,
    double mean_motion_rad_per_s,
    const control::RelativePdGains& gains,
    double brake_decel_mps2 = 1.0e-3,
    double track_time_constant_s = 10.0) {
    if (sequence.empty()) {
        throw std::domain_error("Rendezvous sequence must contain at least one hold point");
    }
    if (active_leg >= sequence.size()) {
        throw std::domain_error("Active leg index is past the sequence end");
    }
    if (!relative::is_finite(estimated_state)) {
        throw std::domain_error("Estimated relative state must be finite");
    }
    detail::require_positive_guidance(brake_decel_mps2, "Brake deceleration must be positive and finite");
    detail::require_positive_guidance(track_time_constant_s, "Track time constant must be positive and finite");
    if (!std::isfinite(mean_motion_rad_per_s) || mean_motion_rad_per_s <= 0.0) {
        throw std::domain_error("Mean motion must be finite and strictly positive");
    }

    const HoldPoint& leg = sequence[active_leg];
    detail::require_finite_waypoint(leg.position_lvlh_m);
    const math::Vector3 rho = estimated_state.relative_position_lvlh_m;
    const math::Vector3 vel = estimated_state.relative_velocity_lvlh_mps;

    const math::Vector3 error = leg.position_lvlh_m - rho;
    const double distance = error.norm();
    const double range = rho.norm();
    const double lateral = lateral_axis_error_m(rho);
    const double closing = closing_speed_mps(rho, vel);

    GuidanceStep step;
    step.active_leg = active_leg;
    step.reference_position_lvlh_m = leg.position_lvlh_m;
    step.distance_to_waypoint_m = distance;
    step.range_to_target_m = range;
    step.lateral_error_m = lateral;
    step.closing_speed_mps = closing;

    // Safety predicates on the estimated state.
    step.keep_out_violation = range < safety.keep_out_radius_m;
    step.corridor_violation =
        (range < safety.corridor_range_m) && (lateral > safety.corridor_radius_m);
    step.closing_speed_warning = closing > safety.closing_speed_hard_limit_mps;

    // Leg capture: advance (monotonic) or declare mission complete.
    if (distance < leg.capture_radius_m) {
        step.leg_complete = true;
        if (active_leg + 1 < sequence.size()) {
            ++active_leg;
            step.active_leg = active_leg;
            step.reference_position_lvlh_m = sequence[active_leg].position_lvlh_m;
        } else {
            step.mission_complete = true;
        }
    }

    // Velocity reference: closing-speed profile toward the (possibly new) target.
    const math::Vector3 fresh_error = step.reference_position_lvlh_m - rho;
    const double fresh_distance = fresh_error.norm();
    math::Vector3 v_des{};
    if (fresh_distance > 1.0e-9) {
        const math::Vector3 los = fresh_error / fresh_distance;
        const double cruise = sequence[step.active_leg].closing_speed_limit_mps;
        const double braking = std::sqrt(2.0 * brake_decel_mps2 * fresh_distance);
        v_des = los * std::min(cruise, braking);
    }
    step.reference_velocity_lvlh_mps = v_des;

    // Acceleration command: first-order velocity tracking + M16 PD/CW structure.
    const math::Vector3 a_track = (v_des - vel) / track_time_constant_s;
    step.commanded_accel_lvlh_mps2 = control::relative_pd_accel_lvlh_mps2(
        rho, vel, step.reference_position_lvlh_m, v_des, a_track, mean_motion_rad_per_s, gains);
    return step;
}

}  // namespace astradock::rendezvous
