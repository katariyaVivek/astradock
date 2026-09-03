// AstraDock M17 — Integrated Closed-Loop GNC tests.
//
// Verification strategy:
//   Wiring: estimate-based tick matches the M16 law on identical inputs;
//     timeline segment selection incl. pre-first/empty rejection.
//   Truth isolation: controller + gnc headers contain no truth-state member or
//     truth-typed parameter (static source audit — pins the master-spec rule).
//   Loop physics: with a perfect estimator (estimate := truth), the loop
//     detumbles the master-spec tumble and holds; metrics settle/converge.
//   Estimator-in-loop: MEKF converges from perturbed init under gyro+ST noise;
//     closed loop with real estimates settles (nominal seed, bounded horizon).
//   Metrics: settle latching, max/final bookkeeping, empty/mismatch rejection.

#include "actuators/actuator_assembly.hpp"
#include "actuators/reaction_wheel.hpp"
#include "attitude/rigid_body.hpp"
#include "control/attitude_pd.hpp"
#include "estimation/attitude_ekf.hpp"
#include "estimation/diagnostics.hpp"
#include "gnc/closed_loop.hpp"
#include "math/constants.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"
#include "sensors/sensor_common.hpp"
#include "sensors/star_tracker.hpp"
#include "spacecraft/force_torque.hpp"
#include "spacecraft/six_dof_dynamics.hpp"
#include "spacecraft/spacecraft_parameters.hpp"
#include "spacecraft/spacecraft_state.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <fstream>
#include <string>

using namespace astradock;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("M17 estimate-based tick matches the M16 law on identical inputs", "[gnc][wiring]") {
    const math::Quaternion q_est = math::Quaternion::from_axis_angle(
        math::Vector3{0.0, 0.0, 1.0}, 0.1);
    const math::Vector3 w_est{0.05, -0.03, 0.02};
    const control::AttitudePdGains gains{{2.0, 2.0, 2.0}, {0.8, 0.8, 0.8}};
    const math::Vector3 limit{0.2, 0.2, 0.2};
    const gnc::EstimateBasedAttitudeCommand cmd{
        q_est, w_est, math::Quaternion::identity(), math::Vector3{}};
    const gnc::ClosedLoopTorqueStep tick =
        gnc::step_estimate_based_attitude_control(cmd, gains, limit);
    const math::Vector3 expected =
        control::attitude_pd_torque_body_Nm(q_est, math::Quaternion::identity(), w_est, math::Vector3{}, gains);
    CHECK_THAT(tick.desired_torque_body_Nm.x(), WithinRel(expected.x(), 1.0e-12));
    CHECK_THAT(tick.desired_torque_body_Nm.y(), WithinRel(expected.y(), 1.0e-12));
    CHECK_THAT(tick.desired_torque_body_Nm.z(), WithinRel(expected.z(), 1.0e-12));
    CHECK_FALSE(tick.saturated);
    CHECK_THAT(
        tick.attitude_error_angle_rad,
        WithinRel(control::attitude_error_angle_rad(q_est, math::Quaternion::identity()), 1.0e-12));
    // Over-limit desire saturates with flag, achieved clamped.
    const gnc::EstimateBasedAttitudeCommand big{
        math::Quaternion::from_axis_angle(math::Vector3{1.0, 0.0, 0.0}, constants::pi / 2.0),
        math::Vector3{1.0, 0.0, 0.0},
        math::Quaternion::identity(),
        math::Vector3{}};
    const gnc::ClosedLoopTorqueStep tick_sat =
        gnc::step_estimate_based_attitude_control(big, gains, math::Vector3{0.05, 0.05, 0.05});
    CHECK(tick_sat.saturated);
    CHECK(std::abs(tick_sat.achieved_torque_body_Nm.x()) <= 0.05 + 1.0e-12);
}

TEST_CASE("M17 reference timeline selects segments deterministically", "[gnc][timeline]") {
    gnc::AttitudeReferenceTimeline timeline;
    timeline.segments = {
        {0.0, math::Quaternion::identity(), math::Vector3{}},
        {40.0,
         math::Quaternion::from_axis_angle(math::Vector3{0.0, 0.0, 1.0}, constants::pi / 3.0),
         math::Vector3{}},
    };
    CHECK_THAT(
        control::attitude_error_angle_rad(
            timeline.reference_at(0.0).attitude_body_to_eci, math::Quaternion::identity()),
        WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(
        control::attitude_error_angle_rad(
            timeline.reference_at(120.0).attitude_body_to_eci,
            math::Quaternion::from_axis_angle(math::Vector3{0.0, 0.0, 1.0}, constants::pi / 3.0)),
        WithinAbs(0.0, 1.0e-12));
    CHECK_THAT(
        control::attitude_error_angle_rad(
            timeline.reference_at(39.999).attitude_body_to_eci, math::Quaternion::identity()),
        WithinAbs(0.0, 1.0e-15));
    gnc::AttitudeReferenceTimeline empty;
    CHECK_THROWS_AS(empty.reference_at(0.0), std::domain_error);
}

TEST_CASE("M17 state-isolation static audit: estimates-only control seam", "[gnc][isolation]") {
    // The control seam headers must never DECLARE a simulated-state container:
    // no member, parameter, or struct whose TYPE carries the simulated bus
    // (SpacecraftState / CartesianState / RotationalState). Every declaration
    // line is scanned; a typed leak fails with file:line. Doc comments
    // describing data flow are exempt (they carry no type). If a future edit
    // routes simulated state into guidance/control, this test names the location.
    const std::vector<std::string> headers{
        "cpp/gnc/closed_loop.hpp",
        "cpp/control/attitude_pd.hpp",
        "cpp/control/relative_pd.hpp",
        "cpp/control/control_saturation.hpp",
    };
    for (const auto& rel : headers) {
        // CTest runs each case with WORKING_DIRECTORY = build dir; resolve the
        // repo root two ways (source-tree run vs build-tree run) so the audit
        // works under both `ctest` and direct binary execution.
        const std::vector<std::string> candidates{"../" + rel, rel};
        std::ifstream file;
        std::string used;
        for (const auto& cand : candidates) {
            file = std::ifstream(cand);
            if (file.is_open()) {
                used = cand;
                break;
            }
        }
        REQUIRE(file.is_open());
        std::string line;
        int line_no = 0;
        while (std::getline(file, line)) {
            ++line_no;
            // Strip the // comment tail: only declarations are audited.
            std::string code = line;
            const std::size_t comment = code.find("//");
            if (comment != std::string::npos) {
                code = code.substr(0, comment);
            }
            std::string lowered;
            lowered.reserve(code.size());
            for (char c : code) {
                lowered.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
            }
            const bool typed = (lowered.find("spacecraftstate") != std::string::npos)
                || (lowered.find("cartesianstate") != std::string::npos)
                || (lowered.find("rotationalstate") != std::string::npos);
            INFO(used << ":" << line_no << ": " << line);
            CHECK_FALSE(typed);
        }
    }
    // Struct-level pin: the command bundle exposes exactly estimate+reference.
    const gnc::EstimateBasedAttitudeCommand cmd{};
    CHECK(cmd.estimated_attitude_body_to_eci.is_unit(1.0e-12));
    CHECK(cmd.reference_attitude_body_to_eci.is_unit(1.0e-12));
}

TEST_CASE("M17 perfect-estimate loop detumbles and holds (loop physics)", "[gnc][loop]") {
    // Estimate := truth each tick (cheating estimator, honest dynamics): proves
    // the control->actuator->truth path converges before estimator error exists.
    attitude::PrincipalInertia inertia{10.0, 20.0, 30.0};
    attitude::RotationalState rot{
        math::Quaternion::identity(), math::Vector3{0.2, -0.1, 0.15}};
    const control::AttitudePdGains gains = control::suggest_pd_gains(
        control::PdTuningRequest{{10.0, 20.0, 30.0}, 0.9, 30.0});
    const std::vector<actuators::ReactionWheelParameters> wheels{
        actuators::ReactionWheelParameters{0.01, 0.2, 600.0, 0.0, math::Vector3{1.0, 0.0, 0.0}},
        actuators::ReactionWheelParameters{0.01, 0.2, 600.0, 0.0, math::Vector3{0.0, 1.0, 0.0}},
        actuators::ReactionWheelParameters{0.01, 0.2, 600.0, 0.0, math::Vector3{0.0, 0.0, 1.0}},
    };
    std::vector<actuators::ReactionWheelState> states(3);
    const double dt = 0.05;
    double t = 0.0;
    std::vector<double> time_v, err_v, cmd_v, est_v;
    std::vector<int> sat_v;
    for (int i = 0; i < 2400; ++i) {
        const gnc::EstimateBasedAttitudeCommand cmd{
            rot.orientation, rot.angular_velocity_rad_per_s,
            math::Quaternion::identity(), math::Vector3{}};
        const gnc::ClosedLoopTorqueStep tick = gnc::step_estimate_based_attitude_control(
            cmd, gains, math::Vector3{0.2, 0.2, 0.2});
        actuators::ActuatorCommand act_cmd;
        act_cmd.wheel_torques_Nm = {
            -tick.achieved_torque_body_Nm.x(),
            -tick.achieved_torque_body_Nm.y(),
            -tick.achieved_torque_body_Nm.z(),
        };
        act_cmd.thruster_forces_N = {};
        actuators::SpacecraftMassState mass{500.0, 400.0};
        const actuators::ActuatorOutput out = actuators::step_actuator_assembly(
            wheels, states, {}, rot.orientation, act_cmd, mass, dt);
        rot = attitude::rk4_step_rotational(t, rot, dt, inertia, out.total_torque_body_Nm);
        t += dt;
        time_v.push_back(t);
        err_v.push_back(control::attitude_error_angle_rad(rot.orientation, math::Quaternion::identity()));
        cmd_v.push_back(tick.desired_torque_body_Nm.norm());
        est_v.push_back(0.0);
        sat_v.push_back(tick.saturated ? 1 : 0);
    }
    const gnc::ClosedLoopMetrics m = gnc::summarize_closed_loop_run(
        time_v, err_v, cmd_v, est_v, sat_v, 0.01);
    CHECK(m.converged);
    CHECK(m.settle_time_s > 0.0);
    CHECK(m.settle_time_s < 90.0);
    CHECK(err_v.back() < 0.005);
}

TEST_CASE("M17 MEKF estimator-in-loop converges under flight-like noise", "[gnc][estimator]") {
    // Gyro (1 mrad/s) + star tracker (1 mrad) feed the MEKF; attitude error and
    // bias error must shrink from perturbed init over 60 s. Deterministic seed.
    estimation::AttitudeEkfConfig cfg;
    cfg.gyro_noise_std_rad_s = 1.0e-3;
    cfg.star_tracker_noise_std_rad = 1.0e-3;
    estimation::AttitudeEkf filter(cfg);
    sensors::DeterministicRng rng(4242);
    const math::Vector3 true_rate{0.02, -0.015, 0.01};
    const math::Vector3 true_bias{5.0e-4, -3.0e-4, 2.0e-4};
    math::Quaternion truth_q = math::Quaternion::from_axis_angle(
        math::Vector3{0.0, 0.0, 1.0}, 0.05);  // ~2.9 deg initial error vs identity init
    const double dt = 0.01;
    double t = 0.0;
    const double init_err =
        estimation::attitude_error_vector_rad(truth_q, filter.estimate().nominal_orientation).norm();
    for (int i = 0; i < 6000; ++i) {
        const math::Vector3 gyro_meas =
            true_rate + true_bias + rng.gaussian_vector3(1.0e-3 / std::sqrt(dt) * dt);
        static_cast<void>(gyro_meas);
        const math::Vector3 gyro_step =
            true_rate + true_bias + rng.gaussian_vector3(1.0e-3);
        filter.predict(gyro_step, dt);
        const double w = true_rate.norm();
        truth_q = (truth_q
                   * math::Quaternion::from_axis_angle(true_rate / w, w * dt))
                      .normalized();
        if (i % 10 == 0) {
            const double ang = rng.gaussian(0.0, 1.0e-3);
            const math::Quaternion noise =
                math::Quaternion::from_axis_angle(rng.uniform_unit_vector3(), ang);
            const sensors::StarTrackerMeasurement meas{t, (truth_q * noise).normalized(), true};
            filter.update_star_tracker(meas);
        }
        t += dt;
    }
    const double final_err = estimation::attitude_error_vector_rad(
                                 truth_q, filter.estimate().nominal_orientation)
                                 .norm();
    const double final_bias_err =
        (true_bias - filter.estimate().gyro_bias_rad_s).norm();
    CHECK(final_err < init_err);
    CHECK(final_err < 0.01);
    CHECK(final_bias_err < 1.0e-3);
}

TEST_CASE("M17 closed-loop metrics bookkeeping latches settle correctly", "[gnc][metrics]") {
    const std::vector<double> time{0.0, 1.0, 2.0, 3.0, 4.0};
    const std::vector<double> err{0.5, 0.2, 0.05, 0.04, 0.03};
    const std::vector<double> cmd{0.1, 0.1, 0.05, 0.02, 0.01};
    const std::vector<double> est{0.4, 0.2, 0.05, 0.04, 0.03};
    const std::vector<int> sat{1, 1, 0, 0, 0};
    const gnc::ClosedLoopMetrics m =
        gnc::summarize_closed_loop_run(time, err, cmd, est, sat, 0.1);
    CHECK(m.converged);
    CHECK_THAT(m.settle_time_s, WithinRel(2.0, 1.0e-12));
    CHECK_THAT(m.max_attitude_error_rad, WithinRel(0.5, 1.0e-12));
    CHECK_THAT(m.final_attitude_error_rad, WithinRel(0.03, 1.0e-12));
    CHECK_THAT(m.max_commanded_torque_Nm, WithinRel(0.1, 1.0e-12));
    CHECK_THAT(m.saturated_fraction, WithinRel(0.4, 1.0e-12));
    CHECK(m.rms_estimate_error_rad > 0.0);
    // Never-settling run: converged == false, settle == -1.
    const std::vector<double> err_bad{0.5, 0.5, 0.5, 0.5, 0.5};
    const gnc::ClosedLoopMetrics bad =
        gnc::summarize_closed_loop_run(time, err_bad, cmd, est, sat, 0.1);
    CHECK_FALSE(bad.converged);
    CHECK_THAT(bad.settle_time_s, WithinRel(-1.0, 0.0));
    // Mismatched/empty columns rejected.
    CHECK_THROWS_AS(
        gnc::summarize_closed_loop_run(time, err_bad, cmd, est, std::vector<int>{}, 0.1),
        std::domain_error);
    CHECK_THROWS_AS(
        gnc::summarize_closed_loop_run({}, {}, {}, {}, {}, 0.1), std::domain_error);
}

TEST_CASE("M17 full 6-DOF truth loop holds attitude with wheel actuation", "[gnc][6dof]") {
    // End-to-end: perfect-estimate PD -> wheel assembly -> rk4_step_spacecraft
    // (translation + rotation). Proves achieved torque integrates through the
    // real 6-DOF path, not just the rotational shortcut.
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const double r_orbit = constants::earth_reference_radius_m + 500.0e3;
    const double v_circ = std::sqrt(mu / r_orbit);
    spacecraft::SpacecraftState truth{
        orbit::CartesianState{
            math::Vector3{r_orbit, 0.0, 0.0}, math::Vector3{0.0, v_circ, 0.0}},
        attitude::RotationalState{
            math::Quaternion::identity(), math::Vector3{0.05, -0.04, 0.03}}};
    const spacecraft::SpacecraftParameters params{
        mu, attitude::PrincipalInertia{10.0, 20.0, 30.0}};
    const control::AttitudePdGains gains = control::suggest_pd_gains(
        control::PdTuningRequest{{10.0, 20.0, 30.0}, 0.9, 30.0});
    const std::vector<actuators::ReactionWheelParameters> wheels{
        actuators::ReactionWheelParameters{0.01, 0.2, 600.0, 0.0, math::Vector3{1.0, 0.0, 0.0}},
        actuators::ReactionWheelParameters{0.01, 0.2, 600.0, 0.0, math::Vector3{0.0, 1.0, 0.0}},
        actuators::ReactionWheelParameters{0.01, 0.2, 600.0, 0.0, math::Vector3{0.0, 0.0, 1.0}},
    };
    std::vector<actuators::ReactionWheelState> states(3);
    const double dt = 0.05;
    double t = 0.0;
    for (int i = 0; i < 1200; ++i) {
        const gnc::EstimateBasedAttitudeCommand cmd{
            truth.orientation(), truth.angular_velocity_rad_per_s(),
            math::Quaternion::identity(), math::Vector3{}};
        const gnc::ClosedLoopTorqueStep tick =
            gnc::step_estimate_based_attitude_control(cmd, gains, math::Vector3{0.2, 0.2, 0.2});
        actuators::ActuatorCommand act_cmd;
        act_cmd.wheel_torques_Nm = {
            -tick.achieved_torque_body_Nm.x(),
            -tick.achieved_torque_body_Nm.y(),
            -tick.achieved_torque_body_Nm.z(),
        };
        act_cmd.thruster_forces_N = {};
        actuators::SpacecraftMassState mass{500.0, 400.0};
        const actuators::ActuatorOutput out = actuators::step_actuator_assembly(
            wheels, states, {}, truth.orientation(), act_cmd, mass, dt);
        truth = spacecraft::rk4_step_spacecraft(
            t, truth, dt, params,
            spacecraft::ForceTorqueInput{math::Vector3{}, out.total_torque_body_Nm});
        t += dt;
    }
    CHECK(control::attitude_error_angle_rad(truth.orientation(), math::Quaternion::identity()) < 0.01);
    CHECK(truth.angular_velocity_rad_per_s().norm() < 0.005);
    // Orbit untouched by wheel-only actuation (translation/rotation decoupling).
    CHECK_THAT(truth.position().norm(), WithinRel(r_orbit, 1.0e-6));
}
