// AstraDock M18 — Rendezvous & proximity operations tests.
//
// Verification strategy:
//   Safety geometry: lateral error hand cases, closing-speed sign convention,
//     keep-out/corridor/warning predicates incl. boundaries.
//   Guidance: capture advance + monotonic legs, mission completion, profile
//     cruise-vs-braking branches, command structure (tracking + CW ff).
//   Default sequence: approach-axis layout, decreasing speed limits.
//   Closed trajectory: waypoint following on the CW plant reaches each leg.
//   Invalid: empty sequence, overrun leg, non-finite state, bad parameters.

#include "control/relative_pd.hpp"
#include "math/constants.hpp"
#include "math/vector3.hpp"
#include "relative/relative_state.hpp"
#include "rendezvous/rendezvous_guidance.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>
#include <vector>

using namespace astradock;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

control::RelativePdGains test_gains() {
    return control::RelativePdGains{
        math::Vector3{1.0e-5, 1.0e-5, 1.0e-5}, math::Vector3{6.0e-3, 6.0e-3, 6.0e-3}};
}

rendezvous::SafetyCorridor test_safety() {
    rendezvous::SafetyCorridor s;
    s.keep_out_radius_m = 25.0;
    s.corridor_radius_m = 50.0;
    s.corridor_range_m = 1000.0;
    s.closing_speed_hard_limit_mps = 3.0;
    return s;
}

}  // namespace

TEST_CASE("M18 lateral error and closing speed match hand calculations", "[rendezvous][geometry]") {
    // Lateral: x = 30, z = 40 -> 50 m off-axis regardless of y.
    CHECK_THAT(
        rendezvous::lateral_axis_error_m(math::Vector3{30.0, -500.0, 40.0}),
        WithinRel(50.0, 1.0e-12));
    CHECK_THAT(
        rendezvous::lateral_axis_error_m(math::Vector3{0.0, -5000.0, 0.0}),
        WithinAbs(0.0, 1.0e-15));
    // Closing: rho = [0,-1000,0], v = [0,+2,0] (flying +y toward origin from -y)
    // -> closing = -((-1000)(2))/1000 = +2 m/s inward.
    CHECK_THAT(
        rendezvous::closing_speed_mps(math::Vector3{0.0, -1000.0, 0.0}, math::Vector3{0.0, 2.0, 0.0}),
        WithinRel(2.0, 1.0e-12));
    // Outbound motion reads negative.
    CHECK_THAT(
        rendezvous::closing_speed_mps(math::Vector3{0.0, -1000.0, 0.0}, math::Vector3{0.0, -1.0, 0.0}),
        WithinRel(-1.0, 1.0e-12));
    // Pure lateral velocity contributes nothing to closing speed.
    CHECK_THAT(
        rendezvous::closing_speed_mps(math::Vector3{0.0, -1000.0, 0.0}, math::Vector3{5.0, 0.0, 0.0}),
        WithinAbs(0.0, 1.0e-12));
}

TEST_CASE("M18 safety predicates flag keep-out, corridor, and speed violations", "[rendezvous][safety]") {
    const double n = 1.1e-3;
    const std::vector<rendezvous::HoldPoint> seq{
        {math::Vector3{0.0, -250.0, 0.0}, 0.5, 10.0},
    };
    const rendezvous::SafetyCorridor safety = test_safety();
    // Inside keep-out (range 10 < 25): violation.
    {
        std::size_t leg = 0;
        const relative::RelativeStateLvlh est{math::Vector3{0.0, -10.0, 0.0}, math::Vector3{}};
        const auto step = rendezvous::step_rendezvous_guidance(est, seq, leg, safety, n, test_gains());
        CHECK(step.keep_out_violation);
        CHECK(step.range_to_target_m < 25.0);
    }
    // Corridor: inside range (500 < 1000) but 60 m off-axis (> 50): violation.
    {
        std::size_t leg = 0;
        const relative::RelativeStateLvlh est{math::Vector3{60.0, -500.0, 0.0}, math::Vector3{}};
        const auto step = rendezvous::step_rendezvous_guidance(est, seq, leg, safety, n, test_gains());
        CHECK_FALSE(step.keep_out_violation);
        CHECK(step.corridor_violation);
    }
    // Same lateral offset beyond corridor range: no corridor flag.
    {
        std::size_t leg = 0;
        const relative::RelativeStateLvlh est{math::Vector3{60.0, -5000.0, 0.0}, math::Vector3{}};
        const auto step = rendezvous::step_rendezvous_guidance(est, seq, leg, safety, n, test_gains());
        CHECK_FALSE(step.corridor_violation);
    }
    // Excess closing speed: 5 m/s inward vs 3 m/s hard limit -> warning.
    {
        std::size_t leg = 0;
        const relative::RelativeStateLvlh est{
            math::Vector3{0.0, -500.0, 0.0}, math::Vector3{0.0, 5.0, 0.0}};
        const auto step = rendezvous::step_rendezvous_guidance(est, seq, leg, safety, n, test_gains());
        CHECK(step.closing_speed_warning);
        CHECK_THAT(step.closing_speed_mps, WithinRel(5.0, 1.0e-12));
    }
    // Boundary: exactly AT the keep-out radius is not a violation (strict <).
    {
        std::size_t leg = 0;
        const relative::RelativeStateLvlh est{math::Vector3{0.0, -25.0, 0.0}, math::Vector3{}};
        const auto step = rendezvous::step_rendezvous_guidance(est, seq, leg, safety, n, test_gains());
        CHECK_FALSE(step.keep_out_violation);
    }
}

TEST_CASE("M18 guidance advances legs on capture and completes the mission", "[rendezvous][sequence]") {
    const double n = 1.1e-3;
    const std::vector<rendezvous::HoldPoint> seq{
        {math::Vector3{0.0, -1000.0, 0.0}, 1.0, 25.0},
        {math::Vector3{0.0, -250.0, 0.0}, 0.5, 10.0},
    };
    const rendezvous::SafetyCorridor safety = test_safety();
    std::size_t leg = 0;
    // Far from leg 0: no advance, reference is leg 0.
    {
        const relative::RelativeStateLvlh est{math::Vector3{0.0, -5000.0, 0.0}, math::Vector3{}};
        const auto step = rendezvous::step_rendezvous_guidance(est, seq, leg, safety, n, test_gains());
        CHECK_FALSE(step.leg_complete);
        CHECK(leg == 0);
        CHECK_FALSE(step.mission_complete);
        CHECK_THAT(step.distance_to_waypoint_m, WithinRel(4000.0, 1.0e-12));
    }
    // Inside leg-0 capture radius: advance to leg 1.
    {
        const relative::RelativeStateLvlh est{math::Vector3{0.0, -1010.0, 0.0}, math::Vector3{}};
        const auto step = rendezvous::step_rendezvous_guidance(est, seq, leg, safety, n, test_gains());
        CHECK(step.leg_complete);
        CHECK(leg == 1);
        CHECK_THAT(step.reference_position_lvlh_m.y(), WithinRel(-250.0, 1.0e-15));
    }
    // Inside final capture radius: mission complete.
    {
        const relative::RelativeStateLvlh est{math::Vector3{0.0, -255.0, 0.0}, math::Vector3{}};
        const auto step = rendezvous::step_rendezvous_guidance(est, seq, leg, safety, n, test_gains());
        CHECK(step.leg_complete);
        CHECK(step.mission_complete);
    }
}

TEST_CASE("M18 closing-speed profile cruises far out and brakes near the waypoint", "[rendezvous][profile]") {
    const double n = 1.1e-3;
    const std::vector<rendezvous::HoldPoint> seq{
        {math::Vector3{0.0, -1000.0, 0.0}, 1.0, 25.0},
    };
    const rendezvous::SafetyCorridor safety = test_safety();
    // Far (4000 m): cruise branch, |v_des| == 1.0 m/s along +y.
    {
        std::size_t leg = 0;
        const relative::RelativeStateLvlh est{math::Vector3{0.0, -5000.0, 0.0}, math::Vector3{}};
        const auto step = rendezvous::step_rendezvous_guidance(
            est, seq, leg, safety, n, test_gains(), 1.0e-3, 10.0);
        CHECK_THAT(step.reference_velocity_lvlh_mps.norm(), WithinRel(1.0, 1.0e-9));
        CHECK_THAT(step.reference_velocity_lvlh_mps.y(), WithinRel(1.0, 1.0e-9));
    }
    // Near (4 m): braking branch sqrt(2*a*d) = sqrt(2*1e-3*4) = 0.0894 m/s.
    {
        std::size_t leg = 0;
        const relative::RelativeStateLvlh est{math::Vector3{0.0, -1004.0, 0.0}, math::Vector3{}};
        const auto step = rendezvous::step_rendezvous_guidance(
            est, seq, leg, safety, n, test_gains(), 1.0e-3, 10.0);
        CHECK_THAT(
            step.reference_velocity_lvlh_mps.norm(),
            WithinRel(std::sqrt(2.0 * 1.0e-3 * 4.0), 1.0e-9));
    }
}

TEST_CASE("M18 default approach sequence marches down the -y axis, slowing inward", "[rendezvous][sequence]") {
    const auto seq = rendezvous::default_approach_sequence();
    REQUIRE(seq.size() == 4);
    CHECK_THAT(seq[0].position_lvlh_m.y(), WithinRel(-5000.0, 1.0e-15));
    CHECK_THAT(seq[3].position_lvlh_m.y(), WithinRel(-50.0, 1.0e-15));
    for (const auto& hp : seq) {
        CHECK_THAT(hp.position_lvlh_m.x(), WithinAbs(0.0, 0.0));
        CHECK_THAT(hp.position_lvlh_m.z(), WithinAbs(0.0, 0.0));
    }
    CHECK(seq[0].closing_speed_limit_mps > seq[1].closing_speed_limit_mps);
    CHECK(seq[1].closing_speed_limit_mps > seq[2].closing_speed_limit_mps);
    CHECK(seq[2].closing_speed_limit_mps > seq[3].closing_speed_limit_mps);
    // Final hold parks outside the keep-out sphere.
    CHECK(seq[3].position_lvlh_m.norm() > 25.0);
}

TEST_CASE("M18 waypoint following on the CW plant reaches each leg", "[rendezvous][tracking]") {
    const double n = 1.1e-3;
    const std::vector<rendezvous::HoldPoint> seq{
        {math::Vector3{0.0, -1000.0, 0.0}, 1.0, 25.0},
        {math::Vector3{0.0, -250.0, 0.0}, 0.5, 10.0},
    };
    const rendezvous::SafetyCorridor safety = test_safety();
    relative::RelativeStateLvlh state{
        math::Vector3{0.0, -2000.0, 0.0}, math::Vector3{0.0, 0.5, 0.0}};
    std::size_t leg = 0;
    const double dt = 1.0;
    bool complete = false;
    for (int i = 0; i < 7200 && !complete; ++i) {
        const auto step = rendezvous::step_rendezvous_guidance(
            state, seq, leg, safety, n, test_gains());
        complete = step.mission_complete;
        // Exact CW drift + command integration (controller-response check).
        const auto drifted = relative::cw_predict(state, n, dt);
        math::Vector3 vel = drifted.relative_velocity_lvlh_mps + step.commanded_accel_lvlh_mps2 * dt;
        math::Vector3 rho =
            drifted.relative_position_lvlh_m + step.commanded_accel_lvlh_mps2 * (0.5 * dt * dt);
        state = {rho, vel};
    }
    CHECK(complete);
    CHECK(leg == 1);
    CHECK_THAT(state.relative_position_lvlh_m.y(), WithinAbs(-250.0, 15.0));
}

TEST_CASE("M18 guidance rejects invalid sequences and states", "[rendezvous][invalid]") {
    const double n = 1.1e-3;
    const std::vector<rendezvous::HoldPoint> seq{
        {math::Vector3{0.0, -250.0, 0.0}, 0.5, 10.0},
    };
    const rendezvous::SafetyCorridor safety = test_safety();
    const relative::RelativeStateLvlh good{math::Vector3{0.0, -500.0, 0.0}, math::Vector3{}};
    const std::vector<rendezvous::HoldPoint> empty;
    std::size_t leg = 0;
    CHECK_THROWS_AS(
        rendezvous::step_rendezvous_guidance(good, empty, leg, safety, n, test_gains()),
        std::domain_error);
    std::size_t overrun = 7;
    CHECK_THROWS_AS(
        rendezvous::step_rendezvous_guidance(good, seq, overrun, safety, n, test_gains()),
        std::domain_error);
    const relative::RelativeStateLvlh bad{
        math::Vector3{std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0}, math::Vector3{}};
    CHECK_THROWS_AS(
        rendezvous::step_rendezvous_guidance(bad, seq, leg, safety, n, test_gains()),
        std::domain_error);
    CHECK_THROWS_AS(
        rendezvous::step_rendezvous_guidance(good, seq, leg, safety, -1.0e-3, test_gains()),
        std::domain_error);
}
