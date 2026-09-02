#include "attitude/rigid_body.hpp"
#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"
#include "dynamics/two_body.hpp"
#include "math/constants.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"
#include "orbit/cartesian_state.hpp"
#include "orbit/two_body_orbit.hpp"
#include "spacecraft/force_torque.hpp"
#include "spacecraft/six_dof_dynamics.hpp"
#include "spacecraft/spacecraft_parameters.hpp"
#include "spacecraft/spacecraft_state.hpp"

#include <cmath>
#include <vector>

using namespace astradock;
using namespace astradock::spacecraft;
using namespace astradock::math;
using namespace astradock::orbit;
using namespace astradock::attitude;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// ===========================================================================
// Test 1: SpacecraftState Vector Space Arithmetic & ADL Finite Check
// ===========================================================================
TEST_CASE("SpacecraftState arithmetic operations satisfy linear vector space axioms", "[spacecraft][6dof]") {
    const SpacecraftState s1{
        CartesianState{Vector3(1.0, 2.0, 3.0), Vector3(4.0, 5.0, 6.0)},
        RotationalState{Quaternion(1.0, 0.0, 0.0, 0.0), Vector3(0.1, 0.2, 0.3)}
    };
    const SpacecraftState s2{
        CartesianState{Vector3(10.0, 20.0, 30.0), Vector3(40.0, 50.0, 60.0)},
        RotationalState{Quaternion(0.5, 0.5, 0.5, 0.5), Vector3(1.0, 2.0, 3.0)}
    };

    const SpacecraftState sum = s1 + s2;
    CHECK_THAT(sum.position().x(), WithinAbs(11.0, 1.0e-14));
    CHECK_THAT(sum.velocity().y(), WithinAbs(55.0, 1.0e-14));
    CHECK_THAT(sum.orientation().w(), WithinAbs(1.5, 1.0e-14));
    CHECK_THAT(sum.angular_velocity_rad_per_s().z(), WithinAbs(3.3, 1.0e-14));

    const SpacecraftState scaled = s1 * 2.5;
    CHECK_THAT(scaled.position().x(), WithinAbs(2.5, 1.0e-14));
    CHECK_THAT(scaled.velocity().z(), WithinAbs(15.0, 1.0e-14));
    CHECK_THAT(scaled.orientation().w(), WithinAbs(2.5, 1.0e-14));
    CHECK_THAT(scaled.angular_velocity_rad_per_s().x(), WithinAbs(0.25, 1.0e-14));

    const SpacecraftState comm_scaled = 2.5 * s1;
    CHECK_THAT(comm_scaled.position().x(), WithinAbs(2.5, 1.0e-14));

    CHECK(is_finite_state(s1));
    CHECK(is_finite_state(s2));
}

// ===========================================================================
// Test 2: Unified Derivative Equivalence with Standalone Subsystems
// ===========================================================================
TEST_CASE("Unified spacecraft derivative exactly evaluates translational and rotational branches", "[spacecraft][6dof]") {
    const SpacecraftParameters params(
        constants::earth_gravitational_parameter_m3_per_s2,
        PrincipalInertia(12.0, 18.0, 24.0)
    );
    const ForceTorqueInput input(Vector3(0.0, 0.0, 0.0), Vector3(0.5, -0.2, 0.1));

    const SpacecraftState state{
        CartesianState{Vector3(7000.0e3, 0.0, 0.0), Vector3(0.0, 7500.0, 0.0)},
        RotationalState{
            Quaternion::from_axis_angle(Vector3(0.0, 1.0, 0.0), 0.4),
            Vector3(0.1, 0.2, -0.1)
        }
    };

    const SpacecraftState unified_deriv = spacecraft_state_derivative(0.0, state, params, input);

    // Standalone M02/M04 translational derivative
    const CartesianState standalone_trans_deriv = two_body_state_derivative(
        0.0, state.translational, params.gravitational_parameter_m3_per_s2);

    // Standalone M09 rotational derivative
    const RotationalState standalone_rot_deriv = rotational_state_derivative(
        0.0, state.rotational, params.inertia, input.torque_body_Nm);

    // Assert exact branch-by-branch agreement
    CHECK_THAT(unified_deriv.translational.position.x(), WithinAbs(standalone_trans_deriv.position.x(), 1.0e-15));
    CHECK_THAT(unified_deriv.translational.position.y(), WithinAbs(standalone_trans_deriv.position.y(), 1.0e-15));
    CHECK_THAT(unified_deriv.translational.position.z(), WithinAbs(standalone_trans_deriv.position.z(), 1.0e-15));
    CHECK_THAT(unified_deriv.translational.velocity.x(), WithinAbs(standalone_trans_deriv.velocity.x(), 1.0e-15));
    CHECK_THAT(unified_deriv.translational.velocity.y(), WithinAbs(standalone_trans_deriv.velocity.y(), 1.0e-15));
    CHECK_THAT(unified_deriv.translational.velocity.z(), WithinAbs(standalone_trans_deriv.velocity.z(), 1.0e-15));

    CHECK_THAT(unified_deriv.rotational.orientation.w(), WithinAbs(standalone_rot_deriv.orientation.w(), 1.0e-15));
    CHECK_THAT(unified_deriv.rotational.orientation.x(), WithinAbs(standalone_rot_deriv.orientation.x(), 1.0e-15));
    CHECK_THAT(unified_deriv.rotational.orientation.y(), WithinAbs(standalone_rot_deriv.orientation.y(), 1.0e-15));
    CHECK_THAT(unified_deriv.rotational.orientation.z(), WithinAbs(standalone_rot_deriv.orientation.z(), 1.0e-15));
    CHECK_THAT(unified_deriv.rotational.angular_velocity_rad_per_s.x(), WithinAbs(standalone_rot_deriv.angular_velocity_rad_per_s.x(), 1.0e-15));
    CHECK_THAT(unified_deriv.rotational.angular_velocity_rad_per_s.y(), WithinAbs(standalone_rot_deriv.angular_velocity_rad_per_s.y(), 1.0e-15));
    CHECK_THAT(unified_deriv.rotational.angular_velocity_rad_per_s.z(), WithinAbs(standalone_rot_deriv.angular_velocity_rad_per_s.z(), 1.0e-15));
}

// ===========================================================================
// Test 3: M04 Standalone Translation Regression Equivalence
// ===========================================================================
TEST_CASE("M10 6-DOF propagation matches standalone M04 orbital trajectory bitwise", "[spacecraft][regression]") {
    const double r_orbit_m = constants::earth_reference_radius_m + 500.0e3;
    const auto ref = compute_circular_orbit_reference(
        constants::earth_gravitational_parameter_m3_per_s2, r_orbit_m);

    const CartesianState m04_init{
        Vector3(ref.orbital_radius_m, 0.0, 0.0),
        Vector3(0.0, ref.speed_m_per_s, 0.0)
    };

    const SpacecraftParameters params(
        ref.gravitational_parameter_m3_per_s2,
        PrincipalInertia(15.0, 20.0, 25.0)
    );

    const SpacecraftState m10_init{
        m04_init,
        RotationalState{Quaternion::identity(), Vector3(0.01, 0.02, 0.03)}
    };

    const double dt = 5.0; // 5s step
    const double duration = 1000.0; // 1000s

    // Standalone M04 propagation using generic fixed step
    const auto m04_samples = numerics::propagate_fixed_step(
        0.0, duration, dt, m04_init, numerics::IntegrationMethod::classical_rk4,
        [&](double t, const CartesianState& s) {
            return two_body_state_derivative(t, s, params.gravitational_parameter_m3_per_s2);
        }
    );

    // M10 6-DOF propagation
    const auto m10_samples = propagate_spacecraft_fixed_step(
        0.0, duration, dt, m10_init, params, ForceTorqueInput{});

    REQUIRE(m04_samples.size() == m10_samples.size());

    double max_pos_diff = 0.0;
    double max_vel_diff = 0.0;

    for (size_t i = 0; i < m04_samples.size(); ++i) {
        const double pos_diff = (m04_samples[i].state.position - m10_samples[i].state.position()).norm();
        const double vel_diff = (m04_samples[i].state.velocity - m10_samples[i].state.velocity()).norm();
        if (pos_diff > max_pos_diff) max_pos_diff = pos_diff;
        if (vel_diff > max_vel_diff) max_vel_diff = vel_diff;
    }

    // Zero / machine-precision difference proves integration does not distort M04 translation
    CHECK_THAT(max_pos_diff, WithinAbs(0.0, 1.0e-12));
    CHECK_THAT(max_vel_diff, WithinAbs(0.0, 1.0e-12));
}

// ===========================================================================
// Test 4: M09 Standalone Attitude Dynamics Regression Equivalence
// ===========================================================================
TEST_CASE("M10 6-DOF propagation matches standalone M09 rotational trajectory bitwise", "[spacecraft][regression]") {
    const PrincipalInertia inertia(10.0, 20.0, 30.0);
    const RotationalState m09_init{
        Quaternion::from_axis_angle(Vector3(1.0, 1.0, 0.0), 0.5),
        Vector3(0.2, 0.3, 0.1)
    };

    const SpacecraftParameters params(
        constants::earth_gravitational_parameter_m3_per_s2,
        inertia
    );

    const SpacecraftState m10_init{
        CartesianState{Vector3(7000.0e3, 0.0, 0.0), Vector3(0.0, 7500.0, 0.0)},
        m09_init
    };

    const double dt = 0.01;
    const int steps = 500;

    RotationalState m09_state = m09_init;
    SpacecraftState m10_state = m10_init;

    double max_rate_diff = 0.0;
    double max_q_err = 0.0;

    for (int i = 0; i < steps; ++i) {
        const double t = i * dt;
        m09_state = rk4_step_rotational(t, m09_state, dt, inertia, Vector3{});
        m10_state = rk4_step_spacecraft(t, m10_state, dt, params, ForceTorqueInput{});

        const double rate_diff = (m09_state.angular_velocity_rad_per_s - m10_state.angular_velocity_rad_per_s()).norm();
        const double q_err = quaternion_orientation_error_rad(m09_state.orientation, m10_state.orientation());

        if (rate_diff > max_rate_diff) max_rate_diff = rate_diff;
        if (q_err > max_q_err) max_q_err = q_err;
    }

    CHECK_THAT(max_rate_diff, WithinAbs(0.0, 1.0e-12));
    CHECK_THAT(max_q_err, WithinAbs(0.0, 1.0e-6));
}

// ===========================================================================
// Test 5: Full 1-Orbit Decoupled Canonical Scenario Invariants
// ===========================================================================
TEST_CASE("Canonical 500 km orbit with tumbling attitude conserves all orbital and rotational invariants", "[spacecraft][invariants]") {
    const double r_orbit_m = constants::earth_reference_radius_m + 500.0e3;
    const auto ref = compute_circular_orbit_reference(
        constants::earth_gravitational_parameter_m3_per_s2, r_orbit_m);

    const SpacecraftParameters params(
        ref.gravitational_parameter_m3_per_s2,
        PrincipalInertia(10.0, 20.0, 30.0)
    );

    const SpacecraftState state0{
        CartesianState{Vector3(ref.orbital_radius_m, 0.0, 0.0), Vector3(0.0, ref.speed_m_per_s, 0.0)},
        RotationalState{
            Quaternion::from_axis_angle(Vector3(0.0, 0.0, 1.0), constants::pi / 4.0),
            Vector3(0.05, 0.08, 0.02)
        }
    };

    const double dt = 1.0;
    const auto samples = propagate_spacecraft_fixed_step(
        0.0, ref.period_s, dt, state0, params, ForceTorqueInput{});

    const double init_orb_energy = specific_orbital_energy_m2_per_s2(
        state0.translational, params.gravitational_parameter_m3_per_s2);
    const Vector3 init_orb_h = specific_angular_momentum_m2_per_s(state0.translational);

    const double init_rot_energy = rotational_kinetic_energy_J(
        state0.angular_velocity_rad_per_s(), params.inertia);
    const Vector3 init_rot_hi = inertial_angular_momentum_kg_m2_per_s(
        state0.rotational, params.inertia);

    double max_orb_energy_drift = 0.0;
    double max_orb_h_drift = 0.0;
    double max_rot_energy_drift = 0.0;
    double max_rot_hi_drift = 0.0;
    double max_q_norm_err = 0.0;

    for (const auto& sample : samples) {
        const double orb_e = specific_orbital_energy_m2_per_s2(
            sample.state.translational, params.gravitational_parameter_m3_per_s2);
        const Vector3 orb_h = specific_angular_momentum_m2_per_s(sample.state.translational);
        const double rot_e = rotational_kinetic_energy_J(
            sample.state.angular_velocity_rad_per_s(), params.inertia);
        const Vector3 rot_hi = inertial_angular_momentum_kg_m2_per_s(
            sample.state.rotational, params.inertia);
        const double q_norm = sample.state.orientation().norm();

        const double orb_e_rel = std::abs(orb_e - init_orb_energy) / std::abs(init_orb_energy);
        const double orb_h_rel = (orb_h - init_orb_h).norm() / init_orb_h.norm();
        const double rot_e_rel = std::abs(rot_e - init_rot_energy) / init_rot_energy;
        const double rot_hi_rel = (rot_hi - init_rot_hi).norm() / init_rot_hi.norm();
        const double q_norm_err = std::abs(q_norm - 1.0);

        if (orb_e_rel > max_orb_energy_drift) max_orb_energy_drift = orb_e_rel;
        if (orb_h_rel > max_orb_h_drift) max_orb_h_drift = orb_h_rel;
        if (rot_e_rel > max_rot_energy_drift) max_rot_energy_drift = rot_e_rel;
        if (rot_hi_rel > max_rot_hi_drift) max_rot_hi_drift = rot_hi_rel;
        if (q_norm_err > max_q_norm_err) max_q_norm_err = q_norm_err;
    }

    // Invariant conservation bounds
    CHECK(max_orb_energy_drift < 1.0e-10);
    CHECK(max_orb_h_drift < 1.0e-12);
    CHECK(max_rot_energy_drift < 1.0e-5);
    CHECK(max_rot_hi_drift < 1.0e-5);
    CHECK(max_q_norm_err < 1.0e-14);
}

// ===========================================================================
// Test 6: Cross-Independence Test (Orbit vs Attitude Decoupling)
// ===========================================================================
TEST_CASE("Translation and attitude dynamics are physically decoupled in baseline model", "[spacecraft][independence]") {
    const SpacecraftParameters params(
        constants::earth_gravitational_parameter_m3_per_s2,
        PrincipalInertia(10.0, 20.0, 30.0)
    );

    const CartesianState orbit0{
        Vector3(6878137.0, 0.0, 0.0),
        Vector3(0.0, 7612.608, 0.0)
    };

    // Spacecraft A: Identity attitude, spinning about X
    const SpacecraftState sc_A{
        orbit0,
        RotationalState{Quaternion::identity(), Vector3(0.1, 0.0, 0.0)}
    };

    // Spacecraft B: 90 deg rotated attitude, tumbling multi-axis
    const SpacecraftState sc_B{
        orbit0,
        RotationalState{
            Quaternion::from_axis_angle(Vector3(0.0, 1.0, 0.0), constants::pi / 2.0),
            Vector3(-0.05, 0.12, 0.08)
        }
    };

    const double dt = 1.0;
    const double duration = 500.0;

    const auto samples_A = propagate_spacecraft_fixed_step(
        0.0, duration, dt, sc_A, params, ForceTorqueInput{});
    const auto samples_B = propagate_spacecraft_fixed_step(
        0.0, duration, dt, sc_B, params, ForceTorqueInput{});

    // Spacecraft A and B must follow the EXACT same orbital path
    for (size_t i = 0; i < samples_A.size(); ++i) {
        const double pos_diff = (samples_A[i].state.position() - samples_B[i].state.position()).norm();
        const double vel_diff = (samples_A[i].state.velocity() - samples_B[i].state.velocity()).norm();
        CHECK_THAT(pos_diff, WithinAbs(0.0, 1.0e-14));
        CHECK_THAT(vel_diff, WithinAbs(0.0, 1.0e-14));
    }
}

// ===========================================================================
// Test 7: External Torque Acceleration with Undisturbed Orbit
// ===========================================================================
TEST_CASE("External torque produces exact analytical rotational acceleration without altering orbit", "[spacecraft][torque]") {
    const SpacecraftParameters params(
        constants::earth_gravitational_parameter_m3_per_s2,
        PrincipalInertia(10.0, 25.0, 40.0)
    );
    const ForceTorqueInput torque_input(Vector3{}, Vector3(0.0, 0.0, 0.8)); // 0.8 N*m on +Z

    const SpacecraftState sc0{
        CartesianState{Vector3(6878137.0, 0.0, 0.0), Vector3(0.0, 7612.608, 0.0)},
        RotationalState{Quaternion::identity(), Vector3(0.0, 0.0, 0.0)}
    };

    const double duration = 100.0;
    const double dt = 0.1;
    const auto samples = propagate_spacecraft_fixed_step(
        0.0, duration, dt, sc0, params, torque_input);

    const auto& final_state = samples.back().state;

    // Analytical angular velocity: alpha_z = 0.8 / 40.0 = 0.02 rad/s^2 => omega_z(100) = 2.0 rad/s
    const double expected_omega_z = 0.02 * duration;
    CHECK_THAT(final_state.angular_velocity_rad_per_s().z(), WithinAbs(expected_omega_z, 1.0e-12));

    // Analytical angle: theta_z = 0.5 * 0.02 * 100^2 = 100 rad
    const double expected_theta_z = 0.5 * 0.02 * duration * duration;
    const Quaternion expected_q = Quaternion::from_axis_angle(Vector3(0.0, 0.0, 1.0), expected_theta_z);
    CHECK(quaternion_orientation_error_rad(final_state.orientation(), expected_q) < 1.0e-4);
}

// ===========================================================================
// Test 8: Quaternion Orientation Error Metric Behavior
// ===========================================================================
TEST_CASE("quaternion_orientation_error_rad computes correct physical rotation angles", "[spacecraft][metrics]") {
    const Quaternion q_id = Quaternion::identity();
    const Quaternion q_neg_id(-1.0, 0.0, 0.0, 0.0);

    // Identity and double cover equivalence
    CHECK_THAT(quaternion_orientation_error_rad(q_id, q_id), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(quaternion_orientation_error_rad(q_id, q_neg_id), WithinAbs(0.0, 1.0e-15));

    // Known 90 deg rotation about X
    const Quaternion q_90_x = Quaternion::from_axis_angle(Vector3(1.0, 0.0, 0.0), constants::pi / 2.0);
    CHECK_THAT(quaternion_orientation_error_rad(q_id, q_90_x), WithinAbs(constants::pi / 2.0, 1.0e-14));

    // Known 180 deg rotation about Y
    const Quaternion q_180_y = Quaternion::from_axis_angle(Vector3(0.0, 1.0, 0.0), constants::pi);
    CHECK_THAT(quaternion_orientation_error_rad(q_id, q_180_y), WithinAbs(constants::pi, 1.0e-14));

    // Commutative symmetry
    const Quaternion q_rand = Quaternion::from_axis_angle(Vector3(1.0, 2.0, 3.0), 0.73);
    CHECK_THAT(
        quaternion_orientation_error_rad(q_90_x, q_rand),
        WithinAbs(quaternion_orientation_error_rad(q_rand, q_90_x), 1.0e-15)
    );
}

// ===========================================================================
// Test 9: Integrated State Timestep Convergence (4th-Order Study)
// ===========================================================================
TEST_CASE("Unified 6-DOF RK4 integrator exhibits fourth-order convergence on composite state", "[spacecraft][convergence]") {
    const SpacecraftParameters params(
        constants::earth_gravitational_parameter_m3_per_s2,
        PrincipalInertia(10.0, 20.0, 30.0)
    );
    const SpacecraftState state0{
        CartesianState{Vector3(6878137.0, 0.0, 0.0), Vector3(0.0, 7612.608, 0.0)},
        RotationalState{
            Quaternion::from_axis_angle(Vector3(1.0, 0.5, 0.2), 0.3),
            Vector3(0.1, 0.15, -0.05)
        }
    };

    // 1. Orbital position convergence over 1000s
    const double orb_time = 1000.0;
    const auto orb_10 = propagate_spacecraft_fixed_step(0.0, orb_time, 10.0, state0, params);
    const auto orb_5 = propagate_spacecraft_fixed_step(0.0, orb_time, 5.0, state0, params);
    const auto orb_2p5 = propagate_spacecraft_fixed_step(0.0, orb_time, 2.5, state0, params);
    const auto orb_ref = propagate_spacecraft_fixed_step(0.0, orb_time, 0.1, state0, params);

    const double err_pos1 = (orb_10.back().state.position() - orb_ref.back().state.position()).norm();
    const double err_pos2 = (orb_5.back().state.position() - orb_ref.back().state.position()).norm();
    const double err_pos3 = (orb_2p5.back().state.position() - orb_ref.back().state.position()).norm();

    const double ratio_pos12 = err_pos1 / err_pos2;
    const double ratio_pos23 = err_pos2 / err_pos3;

    CHECK(ratio_pos12 > 14.5);
    CHECK(ratio_pos12 < 17.5);
    CHECK(ratio_pos23 > 14.5);
    CHECK(ratio_pos23 < 17.5);

    // 2. Rotational angular rate convergence over 10s (unnormalized for pure ODE order verification)
    const double rot_time = 10.0;
    const auto rot_1 = propagate_spacecraft_fixed_step(0.0, rot_time, 0.1, state0, params, ForceTorqueInput{}, false);
    const auto rot_2 = propagate_spacecraft_fixed_step(0.0, rot_time, 0.05, state0, params, ForceTorqueInput{}, false);
    const auto rot_3 = propagate_spacecraft_fixed_step(0.0, rot_time, 0.025, state0, params, ForceTorqueInput{}, false);
    const auto rot_ref = propagate_spacecraft_fixed_step(0.0, rot_time, 0.001, state0, params, ForceTorqueInput{}, false);

    const double err_rate1 = (rot_1.back().state.angular_velocity_rad_per_s() - rot_ref.back().state.angular_velocity_rad_per_s()).norm();
    const double err_rate2 = (rot_2.back().state.angular_velocity_rad_per_s() - rot_ref.back().state.angular_velocity_rad_per_s()).norm();
    const double err_rate3 = (rot_3.back().state.angular_velocity_rad_per_s() - rot_ref.back().state.angular_velocity_rad_per_s()).norm();

    const double ratio_rate12 = err_rate1 / err_rate2;
    const double ratio_rate23 = err_rate2 / err_rate3;

    CHECK(ratio_rate12 > 14.0);
    CHECK(ratio_rate12 < 17.5);
    CHECK(ratio_rate23 > 14.0);
    CHECK(ratio_rate23 < 17.5);
}

// ===========================================================================
// Test 10: Deterministic Repeatability
// ===========================================================================
TEST_CASE("Identical 6-DOF simulations produce bitwise identical outputs", "[spacecraft][determinism]") {
    const SpacecraftParameters params(
        constants::earth_gravitational_parameter_m3_per_s2,
        PrincipalInertia(10.0, 20.0, 30.0)
    );
    const SpacecraftState state0{
        CartesianState{Vector3(6878137.0, 0.0, 0.0), Vector3(0.0, 7612.608, 0.0)},
        RotationalState{
            Quaternion::from_axis_angle(Vector3(1.0, 1.0, 1.0), 0.5),
            Vector3(0.02, 0.04, -0.01)
        }
    };

    const auto run1 = propagate_spacecraft_fixed_step(0.0, 100.0, 0.1, state0, params);
    const auto run2 = propagate_spacecraft_fixed_step(0.0, 100.0, 0.1, state0, params);

    REQUIRE(run1.size() == run2.size());
    for (size_t i = 0; i < run1.size(); ++i) {
        CHECK(run1[i].time_s == run2[i].time_s);
        CHECK(run1[i].state.position().x() == run2[i].state.position().x());
        CHECK(run1[i].state.position().y() == run2[i].state.position().y());
        CHECK(run1[i].state.position().z() == run2[i].state.position().z());
        CHECK(run1[i].state.velocity().x() == run2[i].state.velocity().x());
        CHECK(run1[i].state.velocity().y() == run2[i].state.velocity().y());
        CHECK(run1[i].state.velocity().z() == run2[i].state.velocity().z());
        CHECK(run1[i].state.orientation().w() == run2[i].state.orientation().w());
        CHECK(run1[i].state.orientation().x() == run2[i].state.orientation().x());
        CHECK(run1[i].state.orientation().y() == run2[i].state.orientation().y());
        CHECK(run1[i].state.orientation().z() == run2[i].state.orientation().z());
        CHECK(run1[i].state.angular_velocity_rad_per_s().x() == run2[i].state.angular_velocity_rad_per_s().x());
        CHECK(run1[i].state.angular_velocity_rad_per_s().y() == run2[i].state.angular_velocity_rad_per_s().y());
        CHECK(run1[i].state.angular_velocity_rad_per_s().z() == run2[i].state.angular_velocity_rad_per_s().z());
    }
}

// ===========================================================================
// Test 11: Defensive Parameter & State Validation
// ===========================================================================
TEST_CASE("Defensive validation rejects invalid spacecraft parameters and degenerate states", "[spacecraft][validation]") {
    const SpacecraftParameters valid_params(
        constants::earth_gravitational_parameter_m3_per_s2,
        PrincipalInertia(10.0, 20.0, 30.0)
    );
    const SpacecraftState valid_state{
        CartesianState{Vector3(6878137.0, 0.0, 0.0), Vector3(0.0, 7612.608, 0.0)},
        RotationalState{Quaternion::identity(), Vector3(0.0, 0.0, 0.0)}
    };

    // Invalid parameters
    CHECK_THROWS_AS(
        SpacecraftParameters(-1.0, valid_params.inertia),
        std::domain_error
    );
    CHECK_THROWS_AS(
        SpacecraftParameters(constants::earth_gravitational_parameter_m3_per_s2, PrincipalInertia(-10.0, 20.0, 30.0)),
        std::domain_error
    );

    // Degenerate position at origin
    const SpacecraftState origin_state{
        CartesianState{Vector3(0.0, 0.0, 0.0), Vector3(0.0, 7612.608, 0.0)},
        valid_state.rotational
    };
    CHECK_THROWS_AS(
        spacecraft_state_derivative(0.0, origin_state, valid_params),
        std::domain_error
    );

    // Non-finite position
    const SpacecraftState nan_pos_state{
        CartesianState{Vector3(std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0), Vector3(0.0, 7612.608, 0.0)},
        valid_state.rotational
    };
    CHECK_THROWS_AS(
        spacecraft_state_derivative(0.0, nan_pos_state, valid_params),
        std::domain_error
    );

    // Non-unit quaternion in orientation error calculation
    CHECK_THROWS_AS(
        quaternion_orientation_error_rad(Quaternion(5.0, 0.0, 0.0, 0.0), Quaternion::identity()),
        std::domain_error
    );
}
