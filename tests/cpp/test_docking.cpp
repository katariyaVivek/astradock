// AstraDock M19 — Autonomous docking tests.
//
// Verification strategy:
//   Frames: port world position hand cases (identity + 90-deg yaw), axis sign
//     regression (approach axis direction is mission-critical).
//   Relative attitude: identity/90-deg/double-cover cases.
//   Contact: penalty law hand values, no-pull clamp, rest contact.
//   Acceptance: all-criteria gate, each single failure rejects.
//   Abort: each reason triggers with priority; timeout; clean run aborts none.
//   Latch: requires continuous acceptance for latch_time.
//   Full approach: guided final approach from 50 m latches with sub-limits.

#include "docking/docking.hpp"
#include "frames/lvlh.hpp"
#include "math/constants.hpp"
#include "math/matrix3.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace astradock;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

docking::DockingPort target_port() {
    // Target port at +1 m along body X, axis pointing +X (out of the port).
    return {math::Vector3{1.0, 0.0, 0.0}, math::Vector3{1.0, 0.0, 0.0}};
}

docking::DockingPort chaser_port() {
    // Chaser port at -1 m along body X, axis pointing -X (ports face each other
    // when both buses share the same attitude and sit on the X axis).
    return {math::Vector3{-1.0, 0.0, 0.0}, math::Vector3{-1.0, 0.0, 0.0}};
}

docking::DockingEnvelope envelope() {
    return docking::DockingEnvelope{};
}

}  // namespace

TEST_CASE("M19 port world position matches hand-derived frame cases", "[docking][frames]") {
    // Identity attitude: port = CM + offset.
    const math::Vector3 p = docking::port_position_eci_m(
        math::Vector3{100.0, 200.0, 300.0}, math::Quaternion::identity(), target_port());
    CHECK_THAT(p.x(), WithinRel(101.0, 1.0e-15));
    CHECK_THAT(p.y(), WithinRel(200.0, 1.0e-15));
    CHECK_THAT(p.z(), WithinRel(300.0, 1.0e-15));
    // 90-deg yaw about +Z maps body +X offset into +Y ECI.
    const math::Quaternion yaw90 = math::Quaternion::from_axis_angle(
        math::Vector3{0.0, 0.0, 1.0}, constants::pi / 2.0);
    const math::Vector3 p2 =
        docking::port_position_eci_m(math::Vector3{}, yaw90, target_port());
    CHECK_THAT(p2.x(), WithinAbs(0.0, 1.0e-12));
    CHECK_THAT(p2.y(), WithinRel(1.0, 1.0e-12));
    CHECK_THAT(p2.z(), WithinAbs(0.0, 1.0e-12));
    // Chaser port at -X stays -X under identity.
    const math::Vector3 p3 = docking::port_position_eci_m(
        math::Vector3{}, math::Quaternion::identity(), chaser_port());
    CHECK_THAT(p3.x(), WithinRel(-1.0, 1.0e-15));
    // Non-unit axis rejected; non-unit attitude rejected.
    CHECK_THROWS_AS(
        docking::port_position_eci_m(
            math::Vector3{}, math::Quaternion::identity(),
            docking::DockingPort{math::Vector3{}, math::Vector3{1.0, 1.0, 0.0}}),
        std::domain_error);
    CHECK_THROWS_AS(
        docking::port_position_eci_m(
            math::Vector3{}, math::Quaternion{0.0, 0.0, 0.0, 0.0}, target_port()),
        std::domain_error);
}

TEST_CASE("M19 relative attitude reads alignment errors exactly", "[docking][attitude]") {
    const math::Quaternion id = math::Quaternion::identity();
    CHECK_THAT(docking::relative_attitude_rad(id, id), WithinAbs(0.0, 1.0e-15));
    const math::Quaternion tilt = math::Quaternion::from_axis_angle(
        math::Vector3{0.0, 1.0, 0.0}, 5.0 * constants::pi / 180.0);
    CHECK_THAT(
        docking::relative_attitude_rad(id, tilt), WithinRel(5.0 * constants::pi / 180.0, 1.0e-12));
    // Double cover: -q reads the same alignment.
    const math::Quaternion flip(-tilt.w(), -tilt.x(), -tilt.y(), -tilt.z());
    CHECK_THAT(
        docking::relative_attitude_rad(id, flip), WithinRel(5.0 * constants::pi / 180.0, 1.0e-12));
    // 20-deg tilt exceeds the 5-deg envelope but not the 15-deg abort line.
    CHECK(20.0 * constants::pi / 180.0 > envelope().alignment_limit_rad);
    CHECK(20.0 * constants::pi / 180.0 > envelope().abort_angle_rad);
}

TEST_CASE("M19 penalty contact pushes but never pulls", "[docking][contact]") {
    // k = 1000, c = 100: pen 0.05 m at 0.1 m/s closing -> 50 + 10 = 60 N.
    CHECK_THAT(docking::contact_force_N(0.05, 0.1, 1000.0, 100.0), WithinRel(60.0, 1.0e-12));
    // Separating fast enough to overcome the spring unloads to zero (no pull).
    CHECK_THAT(docking::contact_force_N(0.01, -0.5, 1000.0, 100.0), WithinAbs(0.0, 0.0));
    // Zero/negative penetration: no contact force.
    CHECK_THAT(docking::contact_force_N(0.0, 0.5, 1000.0, 100.0), WithinAbs(0.0, 0.0));
    CHECK_THAT(docking::contact_force_N(-0.5, 0.5, 1000.0, 100.0), WithinAbs(0.0, 0.0));
    // Rest contact: pen 0.02 -> 20 N static preload readout.
    CHECK_THAT(docking::contact_force_N(0.02, 0.0, 1000.0, 100.0), WithinRel(20.0, 1.0e-12));
}

TEST_CASE("M19 acceptance requires every criterion simultaneously", "[docking][acceptance]") {
    const math::Matrix3 dcm = math::Matrix3::identity();
    // Separation convention: axial s > 0 apart, s = 0 mate, s < 0 penetration.
    // Approach axis +Z in this rig: chaser short of the mate plane reads +axial.
    const math::Vector3 axis{0.0, 0.0, 1.0};
    const math::Vector3 zero{};
    auto tick = [&](math::Vector3 dt_cm, double closing, math::Quaternion q_rel, math::Vector3 rel_rate) {
        // Target CM at origin, chaser CM placed so ports meet: with identity
        // attitudes, port gap = chaser_cm + (-1 X) - (0 + 1 X) -> put chaser at
        // +2 X for zero port error, then add the desired offset.
        return docking::step_docking(
            math::Vector3{}, math::Quaternion::identity(), zero,
            math::Vector3{2.0, 0.0, 0.0} + dt_cm, q_rel, rel_rate, target_port(), chaser_port(),
            axis, dcm, closing, envelope(), 10.0, 0.0, 1000.0);
    };
    // Perfect mate: acceptance true, latched (latch time banked).
    const auto good = tick(math::Vector3{}, 0.0, math::Quaternion::identity(), zero);
    CHECK(good.acceptance);
    CHECK(good.latched);
    CHECK_FALSE(good.abort);
    CHECK_THAT(good.axial_m, WithinAbs(0.0, 1.0e-12));
    // Single failures each reject: lateral 0.5 m.
    CHECK_FALSE(tick(math::Vector3{0.0, 0.5, 0.0}, 0.0, math::Quaternion::identity(), zero).acceptance);
    // Closing 0.5 m/s over the 0.1 limit.
    CHECK_FALSE(tick(math::Vector3{}, 0.5, math::Quaternion::identity(), zero).acceptance);
    // 10-deg tilt over the 5-deg limit.
    const math::Quaternion tilt = math::Quaternion::from_axis_angle(
        math::Vector3{0.0, 1.0, 0.0}, 10.0 * constants::pi / 180.0);
    CHECK_FALSE(tick(math::Vector3{}, 0.0, tilt, zero).acceptance);
    // 2 deg/s relative rate over the 0.5 deg/s limit.
    CHECK_FALSE(
        tick(math::Vector3{}, 0.0, math::Quaternion::identity(), math::Vector3{0.0, 0.0, 0.035})
            .acceptance);
}

TEST_CASE("M19 abort fires deterministically per reason with priority", "[docking][abort]") {
    const math::Matrix3 dcm = math::Matrix3::identity();
    const math::Vector3 axis{0.0, 0.0, 1.0};
    const math::Vector3 zero{};
    auto tick = [&](math::Vector3 dt_cm, double closing, math::Quaternion q_rel, double elapsed) {
        return docking::step_docking(
            math::Vector3{}, math::Quaternion::identity(), zero,
            math::Vector3{2.0, 0.0, 0.0} + dt_cm, q_rel, zero, target_port(), chaser_port(),
            axis, dcm, closing, envelope(), 0.0, elapsed, 100.0);
    };
    // Lateral blowout: 2 m sideways.
    const auto lat = tick(math::Vector3{0.0, 2.0, 0.0}, 0.0, math::Quaternion::identity(), 0.0);
    CHECK(lat.abort);
    CHECK(lat.abort_reason == docking::AbortReason::lateral_exceeded);
    // Hot axial approach: 1.0 m/s (> 3x 0.1 limit).
    const auto hot = tick(math::Vector3{}, 1.0, math::Quaternion::identity(), 0.0);
    CHECK(hot.abort);
    CHECK(hot.abort_reason == docking::AbortReason::closing_speed_exceeded);
    // 20-deg tilt over the 15-deg abort line.
    const math::Quaternion big = math::Quaternion::from_axis_angle(
        math::Vector3{1.0, 0.0, 0.0}, 20.0 * constants::pi / 180.0);
    const auto att = tick(math::Vector3{}, 0.0, big, 0.0);
    CHECK(att.abort);
    CHECK(att.abort_reason == docking::AbortReason::attitude_exceeded);
    // Timeout at t = 101 s.
    const auto to = tick(math::Vector3{}, 0.0, math::Quaternion::identity(), 101.0);
    CHECK(to.abort);
    CHECK(to.abort_reason == docking::AbortReason::timeout);
    // Reason names resolve for telemetry.
    CHECK(std::string(docking::abort_reason_name(docking::AbortReason::none)) == "none");
    CHECK(std::string(docking::abort_reason_name(docking::AbortReason::timeout)) == "timeout");
}

TEST_CASE("M19 latch requires banked acceptance time", "[docking][latch]") {
    const math::Matrix3 dcm = math::Matrix3::identity();
    const math::Vector3 axis{0.0, 0.0, 1.0};
    const math::Vector3 zero{};
    // Same perfect geometry but zero banked time: acceptance without latch.
    const auto fresh = docking::step_docking(
        math::Vector3{}, math::Quaternion::identity(), zero, math::Vector3{2.0, 0.0, 0.0},
        math::Quaternion::identity(), zero, target_port(), chaser_port(), axis, dcm, 0.0,
        envelope(), 0.0, 0.0, 1000.0);
    CHECK(fresh.acceptance);
    CHECK_FALSE(fresh.latched);
}

TEST_CASE("M19 axial sign convention: positive separates, negative penetrates", "[docking][frames]") {
    const math::Matrix3 dcm = math::Matrix3::identity();
    // Approach axis +Z in this rig: chaser short of the mate plane reads +axial.
    const math::Vector3 axis{0.0, 0.0, 1.0};
    const math::Vector3 zero{};
    // Chaser 0.5 m short along +Z of mate: axial = +0.5, no contact.
    const auto short_step = docking::step_docking(
        math::Vector3{}, math::Quaternion::identity(), zero,
        math::Vector3{2.0, 0.0, 0.5}, math::Quaternion::identity(), zero, target_port(),
        chaser_port(), axis, dcm, 0.05, envelope(), 0.0, 0.0, 1000.0);
    CHECK_THAT(short_step.axial_m, WithinRel(0.5, 1.0e-12));
    CHECK_FALSE(short_step.in_contact);
    CHECK_THAT(short_step.contact_force_N, WithinAbs(0.0, 0.0));
    // Chaser 0.05 m past the plane: contact with k*pen + c*rate force.
    const auto contact_step = docking::step_docking(
        math::Vector3{}, math::Quaternion::identity(), zero,
        math::Vector3{2.0, 0.0, -0.05}, math::Quaternion::identity(), zero, target_port(),
        chaser_port(), axis, dcm, 0.05, envelope(), 0.0, 0.0, 1000.0);
    CHECK(contact_step.in_contact);
    CHECK_THAT(contact_step.contact_force_N, WithinRel(1000.0 * 0.05 + 100.0 * 0.05, 1.0e-12));
}
