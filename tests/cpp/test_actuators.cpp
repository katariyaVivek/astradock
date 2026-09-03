// AstraDock M14 — Spacecraft Actuator Dynamics & Modeling tests.
//
// Verification strategy (mirrors the M14 spec):
//   M14A reaction wheels: momentum-exchange sign, exact speed integration,
//     torque saturation reporting, motor lag closed form, momentum/energy helpers.
//   M14B wheel saturation: deterministic buildup to the speed limit, authority
//     loss on further demand, and continued authority in the recovery direction.
//   M14C thrusters: body force/torque geometry (F = dir*T, tau = r x F), ECI
//     rotation through attitude, saturation reporting, zero command.
//   M14D propellant: m_dot = -T/(Isp*g0), exact depletion over a burn, dry-floor
//     clamp with achieved-thrust scale-down, mass conservation.
//   Assembly: multi-wheel momentum sum, thruster aggregation, shared mass state.
//   Regression: zero actuator command reproduces the no-input truth path.

#include "actuators/actuator_assembly.hpp"
#include "actuators/reaction_wheel.hpp"
#include "actuators/thruster.hpp"
#include "math/constants.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

using namespace astradock;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

actuators::ReactionWheelParameters test_wheel(
    double inertia = 1.0e-3,
    double max_torque = 0.05,
    double max_speed = 500.0,
    double lag = 0.0,
    math::Vector3 axis = math::Vector3{0.0, 0.0, 1.0}) {
    return actuators::ReactionWheelParameters{inertia, max_torque, max_speed, lag, axis};
}

actuators::ThrusterParameters test_thruster(
    double max_thrust = 1.0,
    double isp = 220.0,
    math::Vector3 mount = math::Vector3{0.0, 0.0, 0.0},
    math::Vector3 dir = math::Vector3{1.0, 0.0, 0.0}) {
    return actuators::ThrusterParameters{max_thrust, isp, mount, dir};
}

}  // namespace

TEST_CASE("M14A reaction torque sign opposes positive motor torque about the spin axis", "[actuators][m14a]") {
    const auto params = test_wheel();
    // +X wheel: positive motor torque must produce -X bus reaction torque.
    const auto px = test_wheel(1.0e-3, 0.05, 500.0, 0.0, math::Vector3{1.0, 0.0, 0.0});
    const math::Vector3 tau = actuators::reaction_torque_body_Nm(px, 0.02);
    CHECK_THAT(tau.x(), WithinRel(-0.02, 1.0e-15));
    CHECK_THAT(tau.y(), WithinAbs(0.0, 1.0e-18));
    CHECK_THAT(tau.z(), WithinAbs(0.0, 1.0e-18));
    // Stored wheel momentum points along +spin axis for positive speed.
    const math::Vector3 hw = actuators::wheel_momentum_body_Nms(px, 100.0);
    CHECK_THAT(hw.x(), WithinRel(1.0e-3 * 100.0, 1.0e-15));
    CHECK_THAT(hw.y(), WithinAbs(0.0, 1.0e-18));
    CHECK_THAT(hw.z(), WithinAbs(0.0, 1.0e-18));
    static_cast<void>(params);
}

TEST_CASE("M14A wheel speed integrates motor torque exactly under constant command", "[actuators][m14a]") {
    const auto params = test_wheel(2.0e-3, 1.0, 1.0e6, 0.0);
    actuators::ReactionWheelState state{10.0, 0.0};
    const double cmd = 0.04;
    const double dt = 0.25;
    const auto step = actuators::step_reaction_wheel(params, state, cmd, dt);
    const double expected_speed = 10.0 + cmd * dt / 2.0e-3;
    CHECK_THAT(step.state.wheel_speed_rad_s, WithinRel(expected_speed, 1.0e-12));
    CHECK_THAT(step.stored_momentum_Nms, WithinRel(2.0e-3 * expected_speed, 1.0e-12));
    // Bus feels the reaction: -axis * applied torque.
    CHECK_THAT(step.reaction_torque_body_Nm.z(), WithinRel(-cmd, 1.0e-12));
    CHECK_FALSE(step.torque_saturated);
    CHECK_FALSE(step.speed_saturated);
    CHECK_FALSE(step.authority_lost);
    // Momentum books close: applied torque * dt == I * d(omega).
    CHECK_THAT(2.0e-3 * (step.state.wheel_speed_rad_s - 10.0), WithinRel(cmd * dt, 1.0e-12));
}

TEST_CASE("M14A torque saturation limits the command and reports it explicitly", "[actuators][m14a]") {
    const auto params = test_wheel(1.0e-3, 0.05, 1.0e6, 0.0);
    const actuators::ReactionWheelState state{};
    const auto over = actuators::step_reaction_wheel(params, state, 0.5, 0.1);
    CHECK(over.torque_saturated);
    CHECK_THAT(over.state.achieved_motor_torque_Nm, WithinRel(0.05, 1.0e-15));
    CHECK_THAT(over.reaction_torque_body_Nm.z(), WithinRel(-0.05, 1.0e-12));
    const auto under = actuators::step_reaction_wheel(params, state, -0.5, 0.1);
    CHECK(under.torque_saturated);
    CHECK_THAT(under.state.achieved_motor_torque_Nm, WithinRel(-0.05, 1.0e-15));
    const auto exact = actuators::step_reaction_wheel(params, state, 0.05, 0.1);
    CHECK_FALSE(exact.torque_saturated);
    // Zero command: wheel coasts, bus feels nothing.
    const actuators::ReactionWheelState spinning{123.0, 0.0};
    const auto coast = actuators::step_reaction_wheel(params, spinning, 0.0, 0.1);
    CHECK_THAT(coast.state.wheel_speed_rad_s, WithinRel(123.0, 1.0e-15));
    CHECK_THAT(coast.reaction_torque_body_Nm.norm(), WithinAbs(0.0, 1.0e-18));
}

TEST_CASE("M14A first-order motor lag matches the analytical step response", "[actuators][m14a]") {
    const double tau_lag = 0.2;
    const double dt = 0.05;
    // Step target 1.0 from rest: achieved = 1 - exp(-dt/tau).
    const double achieved = actuators::apply_motor_lag(0.0, 1.0, dt, tau_lag);
    CHECK_THAT(achieved, WithinRel(1.0 - std::exp(-dt / tau_lag), 1.0e-14));
    // Cascaded steps compose: steady state approaches the target, never overshoots.
    double a = 0.0;
    for (int i = 0; i < 200; ++i) {
        a = actuators::apply_motor_lag(a, 1.0, dt, tau_lag);
    }
    CHECK_THAT(a, WithinRel(1.0, 1.0e-9));
    // Zero time constant is the ideal instantaneous motor.
    CHECK_THAT(actuators::apply_motor_lag(0.3, 0.9, dt, 0.0), WithinRel(0.9, 1.0e-15));
    // Through the full step: achieved torque feeds the speed integrator.
    const auto params = test_wheel(1.0e-3, 10.0, 1.0e6, tau_lag);
    const actuators::ReactionWheelState rest{0.0, 0.0};
    const auto step = actuators::step_reaction_wheel(params, rest, 1.0, dt);
    CHECK_THAT(step.state.achieved_motor_torque_Nm, WithinRel(achieved, 1.0e-14));
    CHECK_THAT(step.state.wheel_speed_rad_s, WithinRel(achieved * dt / 1.0e-3, 1.0e-12));
}

TEST_CASE("M14B deterministic momentum buildup saturates the wheel and loses authority", "[actuators][m14b]") {
    // Sizing: I = 1e-3, max speed 100 rad/s -> capacity 0.1 Nms. Drive at 0.05 Nm:
    // time to saturate from rest = I*max_speed/tau = 2 s.
    const auto params = test_wheel(1.0e-3, 0.05, 100.0, 0.0);
    actuators::ReactionWheelState state{};
    const double dt = 0.1;
    bool saw_speed_saturation = false;
    for (int i = 0; i < 20; ++i) {
        const auto step = actuators::step_reaction_wheel(params, state, 0.05, dt);
        state = step.state;
        if (i < 19) {
            CHECK_FALSE(step.authority_lost);
        }
    }
    CHECK_THAT(state.wheel_speed_rad_s, WithinRel(100.0, 1.0e-12));
    const auto at_limit = actuators::step_reaction_wheel(params, state, 0.05, dt);
    CHECK(at_limit.speed_saturated);
    CHECK(at_limit.authority_lost);
    CHECK_THAT(at_limit.state.wheel_speed_rad_s, WithinRel(100.0, 1.0e-15));
    CHECK_THAT(at_limit.reaction_torque_body_Nm.norm(), WithinAbs(0.0, 1.0e-18));
    static_cast<void>(saw_speed_saturation);
    // Recovery direction retains authority: negative torque spins back down.
    const auto recover = actuators::step_reaction_wheel(params, state, -0.05, dt);
    CHECK_FALSE(recover.authority_lost);
    CHECK(recover.state.wheel_speed_rad_s < 100.0);
    CHECK(recover.reaction_torque_body_Nm.z() > 0.0);
    // Stored energy at the limit matches 0.5*I*w^2.
    CHECK_THAT(
        actuators::wheel_kinetic_energy_J(params, 100.0), WithinRel(0.5 * 1.0e-3 * 100.0 * 100.0, 1.0e-15));
}

TEST_CASE("M14C thruster force direction and moment-arm torque in the body frame", "[actuators][m14c]") {
    // Offset thruster: mount (0, 0.5, 0), firing +X with T = 2 N.
    // F_B = [2, 0, 0]; tau_B = r x F = [0,0.5,0] x [2,0,0] = [0, 0, -1.0] N*m.
    const auto params = test_thruster(5.0, 220.0, math::Vector3{0.0, 0.5, 0.0}, math::Vector3{1.0, 0.0, 0.0});
    const math::Vector3 force = actuators::thruster_force_body_N(params, 2.0);
    CHECK_THAT(force.x(), WithinRel(2.0, 1.0e-15));
    CHECK_THAT(force.y(), WithinAbs(0.0, 1.0e-18));
    CHECK_THAT(force.z(), WithinAbs(0.0, 1.0e-18));
    const math::Vector3 torque = actuators::thruster_torque_body_Nm(params, 2.0);
    CHECK_THAT(torque.x(), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(torque.y(), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(torque.z(), WithinRel(-1.0, 1.0e-12));
    // Centered thruster produces pure force, zero torque.
    const auto centered = test_thruster();
    CHECK_THAT(actuators::thruster_torque_body_Nm(centered, 1.0).norm(), WithinAbs(0.0, 1.0e-18));
    // Zero thrust: zero force and torque.
    CHECK_THAT(actuators::thruster_force_body_N(params, 0.0).norm(), WithinAbs(0.0, 1.0e-18));
    CHECK_THAT(actuators::thruster_torque_body_Nm(params, 0.0).norm(), WithinAbs(0.0, 1.0e-18));
}

TEST_CASE("M14C thruster force rotates into ECI through the attitude quaternion", "[actuators][m14c]") {
    // 90 deg yaw about +Z maps body +X force into inertial +Y.
    const auto params = test_thruster();
    const math::Quaternion yaw90 =
        math::Quaternion::from_axis_angle(math::Vector3{0.0, 0.0, 1.0}, constants::pi / 2.0);
    const math::Vector3 f_eci = actuators::thruster_force_eci_N(params, yaw90, 3.0);
    CHECK_THAT(f_eci.x(), WithinAbs(0.0, 1.0e-12));
    CHECK_THAT(f_eci.y(), WithinRel(3.0, 1.0e-12));
    CHECK_THAT(f_eci.z(), WithinAbs(0.0, 1.0e-12));
    // Identity attitude: ECI force equals body force.
    const math::Vector3 f_id =
        actuators::thruster_force_eci_N(params, math::Quaternion::identity(), 3.0);
    CHECK_THAT(f_id.x(), WithinRel(3.0, 1.0e-15));
    // Command above max thrust saturates and reports.
    actuators::SpacecraftMassState mass{500.0, 400.0};
    const double wet_before = mass.wet_mass_kg;
    const auto step = actuators::step_thruster(
        params, math::Quaternion::identity(), 10.0, mass.wet_mass_kg, mass.dry_mass_kg, 1.0);
    CHECK(step.thrust_saturated);
    CHECK_THAT(step.thrust_achieved_N, WithinRel(1.0, 1.0e-15));
    CHECK_THAT(step.propellant_used_kg, WithinRel(1.0 / (220.0 * 9.80665), 1.0e-12));
    static_cast<void>(wet_before);
}

TEST_CASE("M14D propellant mass flow and depletion follow m_dot = -T/(Isp*g0)", "[actuators][m14d]") {
    const auto params = test_thruster(10.0, 300.0);
    // Instantaneous rate check against the governing equation.
    CHECK_THAT(
        actuators::thruster_mass_flow_rate_kg_per_s(params, 5.0),
        WithinRel(-5.0 / (300.0 * 9.80665), 1.0e-14));
    CHECK_THAT(actuators::thruster_mass_flow_rate_kg_per_s(params, 0.0), WithinAbs(0.0, 0.0));
    // 10 s burn at 5 N: exact depletion delta_m = -T*dt/(Isp*g0).
    const double wet = 500.0;
    const double dry = 400.0;
    const auto step = actuators::step_thruster(
        params, math::Quaternion::identity(), 5.0, wet, dry, 10.0);
    const double expected_used = 5.0 * 10.0 / (300.0 * 9.80665);
    CHECK_THAT(step.propellant_used_kg, WithinRel(expected_used, 1.0e-12));
    CHECK_THAT(step.spacecraft_mass_after_kg, WithinRel(wet - expected_used, 1.0e-12));
    CHECK_FALSE(step.propellant_depleted);
    CHECK_FALSE(step.thrust_saturated);
    // Propellant floor: nearly-dry spacecraft clamps at dry mass exactly.
    const auto floor_step = actuators::step_thruster(
        params, math::Quaternion::identity(), 10.0, dry + 1.0e-6, dry, 100.0);
    CHECK(floor_step.propellant_depleted);
    CHECK_THAT(floor_step.spacecraft_mass_after_kg, WithinRel(dry, 1.0e-15));
    CHECK_THAT(floor_step.propellant_used_kg, WithinAbs(1.0e-6, 1.0e-12));
    CHECK(floor_step.thrust_achieved_N < 10.0);
    CHECK(floor_step.thrust_achieved_N > 0.0);
    // Mass conservation across a burn: wet_after + exhaust == wet_before.
    CHECK_THAT(step.spacecraft_mass_after_kg + step.propellant_used_kg, WithinRel(wet, 1.0e-15));
}

TEST_CASE("M14 assembly aggregates wheels, thrusters, momentum, and shared mass", "[actuators][assembly]") {
    // Three orthogonal wheels (X, Y, Z axes).
    const std::vector<actuators::ReactionWheelParameters> wheels{
        test_wheel(1.0e-3, 0.05, 500.0, 0.0, math::Vector3{1.0, 0.0, 0.0}),
        test_wheel(1.0e-3, 0.05, 500.0, 0.0, math::Vector3{0.0, 1.0, 0.0}),
        test_wheel(1.0e-3, 0.05, 500.0, 0.0, math::Vector3{0.0, 0.0, 1.0}),
    };
    std::vector<actuators::ReactionWheelState> states(3);
    const std::vector<actuators::ThrusterParameters> thrusters{
        test_thruster(2.0, 220.0, math::Vector3{0.0, 0.0, 0.0}, math::Vector3{1.0, 0.0, 0.0}),
        test_thruster(2.0, 220.0, math::Vector3{0.5, 0.0, 0.0}, math::Vector3{0.0, 1.0, 0.0}),
    };
    actuators::SpacecraftMassState mass{500.0, 400.0};
    actuators::ActuatorCommand cmd;
    cmd.wheel_torques_Nm = {0.01, -0.02, 0.0};
    cmd.thruster_forces_N = {1.5, 0.5};

    const auto out = actuators::step_actuator_assembly(
        wheels, states, thrusters, math::Quaternion::identity(), cmd, mass, 0.1);

    // Wheel reaction torque: -[0.01, -0.02, 0] = [-0.01, +0.02, 0].
    CHECK_THAT(out.wheel_reaction_torque_body_Nm.x(), WithinRel(-0.01, 1.0e-12));
    CHECK_THAT(out.wheel_reaction_torque_body_Nm.y(), WithinRel(0.02, 1.0e-12));
    // Thruster 2 mount (0.5,0,0) firing +Y at 0.5 N: tau = r x F = [0, 0, 0.25].
    CHECK_THAT(out.thruster_torque_body_Nm.z(), WithinRel(0.25, 1.0e-12));
    CHECK_THAT(out.total_torque_body_Nm.z(), WithinRel(0.25, 1.0e-12));
    // Total body force [1.5, 0.5, 0]; identity attitude so ECI matches.
    CHECK_THAT(out.total_force_body_N.x(), WithinRel(1.5, 1.0e-15));
    CHECK_THAT(out.total_force_body_N.y(), WithinRel(0.5, 1.0e-15));
    CHECK_THAT(out.total_force_eci_N.x(), WithinRel(1.5, 1.0e-15));
    // Wheel momentum: each wheel gained tau*dt/I in speed; H = I*w.
    CHECK_THAT(out.total_wheel_momentum_body_x_Nms, WithinRel(0.01 * 0.1, 1.0e-12));
    CHECK_THAT(out.total_wheel_momentum_body_y_Nms, WithinRel(-0.02 * 0.1, 1.0e-12));
    // Shared mass depleted by both thrusters.
    const double expected_use = (1.5 + 0.5) * 0.1 / (220.0 * 9.80665);
    CHECK_THAT(out.propellant_used_kg, WithinRel(expected_use, 1.0e-12));
    CHECK_THAT(mass.wet_mass_kg, WithinRel(500.0 - expected_use, 1.0e-12));
    CHECK_FALSE(out.any_torque_saturated);
    CHECK_FALSE(out.any_authority_lost);
}

TEST_CASE("M14 zero actuator command leaves the spacecraft input identically zero", "[actuators][regression]") {
    const std::vector<actuators::ReactionWheelParameters> wheels{
        test_wheel(), test_wheel(1.0e-3, 0.05, 500.0, 0.0, math::Vector3{1.0, 0.0, 0.0}),
    };
    std::vector<actuators::ReactionWheelState> states(2);
    const std::vector<actuators::ThrusterParameters> thrusters{test_thruster()};
    actuators::SpacecraftMassState mass{500.0, 400.0};
    actuators::ActuatorCommand cmd;
    cmd.wheel_torques_Nm = {0.0, 0.0};
    cmd.thruster_forces_N = {0.0};
    const auto out = actuators::step_actuator_assembly(
        wheels, states, thrusters, math::Quaternion::identity(), cmd, mass, 0.1);
    CHECK_THAT(out.total_force_body_N.norm(), WithinAbs(0.0, 0.0));
    CHECK_THAT(out.total_force_eci_N.norm(), WithinAbs(0.0, 0.0));
    CHECK_THAT(out.total_torque_body_Nm.norm(), WithinAbs(0.0, 0.0));
    CHECK_THAT(out.propellant_used_kg, WithinAbs(0.0, 0.0));
    CHECK_THAT(mass.wet_mass_kg, WithinRel(500.0, 1.0e-15));
    CHECK_FALSE(out.any_torque_saturated);
    CHECK_FALSE(out.any_speed_saturated);
    CHECK_FALSE(out.any_authority_lost);
    CHECK_FALSE(out.any_thrust_saturated);
    CHECK_FALSE(out.any_propellant_depleted);
}

TEST_CASE("M14 actuator inputs reject invalid commands and configurations", "[actuators][invalid]") {
    const auto wheel = test_wheel();
    const actuators::ReactionWheelState state{};
    CHECK_THROWS_AS(actuators::step_reaction_wheel(wheel, state, 0.01, 0.0), std::domain_error);
    CHECK_THROWS_AS(actuators::step_reaction_wheel(wheel, state, 0.01, -0.1), std::domain_error);
    CHECK_THROWS_AS(actuators::ReactionWheelParameters(0.0, 0.05, 500.0, 0.0, math::Vector3{0.0, 0.0, 1.0}),
        std::domain_error);
    CHECK_THROWS_AS(actuators::ReactionWheelParameters(1.0e-3, 0.05, 500.0, 0.0, math::Vector3{1.0, 1.0, 0.0}),
        std::domain_error);
    const auto thr = test_thruster();
    CHECK_THROWS_AS(
        actuators::step_thruster(thr, math::Quaternion::identity(), -1.0, 500.0, 400.0, 0.1),
        std::domain_error);
    CHECK_THROWS_AS(
        actuators::step_thruster(thr, math::Quaternion::identity(), 1.0, 300.0, 400.0, 0.1),
        std::domain_error);
    CHECK_THROWS_AS(actuators::ThrusterParameters(1.0, 220.0, math::Vector3{}, math::Vector3{1.0, 1.0, 0.0}),
        std::domain_error);
}
