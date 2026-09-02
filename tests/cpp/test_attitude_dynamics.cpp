#include "attitude/attitude_state.hpp"
#include "attitude/principal_inertia.hpp"
#include "attitude/rigid_body.hpp"
#include "attitude/rotational_state.hpp"
#include "math/constants.hpp"
#include "math/matrix3.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"
#include "numerics/integrators.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

using astradock::attitude::body_angular_momentum_kg_m2_per_s;
using astradock::attitude::euler_rotational_acceleration;
using astradock::attitude::inertial_angular_momentum_kg_m2_per_s;
using astradock::attitude::is_finite;
using astradock::attitude::is_valid;
using astradock::attitude::PrincipalInertia;
using astradock::attitude::quaternion_derivative;
using astradock::attitude::rk4_step_rotational;
using astradock::attitude::rotational_kinetic_energy_J;
using astradock::attitude::rotational_state_derivative;
using astradock::attitude::RotationalState;
using astradock::constants::pi;
using astradock::math::approximately_equal;
using astradock::math::Matrix3;
using astradock::math::Quaternion;
using astradock::math::represents_same_rotation;
using astradock::math::Vector3;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// ===========================================================================
// Test 1: Zero Angular Velocity and Zero Applied Torque
// ===========================================================================
TEST_CASE("Zero angular velocity and zero torque maintain stationary state", "[attitude][dynamics]") {
    const PrincipalInertia inertia(10.0, 20.0, 30.0);
    const RotationalState initial_state{
        Quaternion::identity(),
        Vector3{}
    };
    const Vector3 zero_torque = Vector3{};

    // Derivative evaluation
    const RotationalState deriv = rotational_state_derivative(0.0, initial_state, inertia, zero_torque);
    CHECK_THAT(deriv.orientation.w(), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(deriv.orientation.x(), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(deriv.orientation.y(), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(deriv.orientation.z(), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(deriv.angular_velocity_rad_per_s.norm(), WithinAbs(0.0, 1.0e-15));

    // Propagation over 10 seconds
    RotationalState state = initial_state;
    const double dt = 0.1;
    for (int step = 0; step < 100; ++step) {
        state = rk4_step_rotational(step * dt, state, dt, inertia, zero_torque);
    }

    CHECK(approximately_equal(state.orientation, Quaternion::identity(), 1.0e-15, 1.0e-15));
    CHECK_THAT(state.angular_velocity_rad_per_s.norm(), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(rotational_kinetic_energy_J(state.angular_velocity_rad_per_s, inertia), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(inertial_angular_momentum_kg_m2_per_s(state, inertia).norm(), WithinAbs(0.0, 1.0e-15));
}

// ===========================================================================
// Test 2: Principal-Axis Constant Torque Analytical Solution
// ===========================================================================
TEST_CASE("Principal-axis constant torque matches exact analytical kinematics and energy", "[attitude][dynamics]") {
    const PrincipalInertia inertia(10.0, 25.0, 40.0);
    const double tau_x = 2.0; // N*m
    const Vector3 torque(tau_x, 0.0, 0.0);
    const RotationalState initial_state{
        Quaternion::identity(),
        Vector3{}
    };

    const double alpha_x = tau_x / inertia.Ixx_kg_m2; // 0.2 rad/s^2
    const double duration_s = 5.0;
    const double dt = 0.01;
    const int num_steps = static_cast<int>(duration_s / dt);

    RotationalState state = initial_state;
    for (int step = 0; step < num_steps; ++step) {
        state = rk4_step_rotational(step * dt, state, dt, inertia, torque);
    }

    // Analytical expectations at t = 5.0 s:
    // omega_x(t) = alpha_x * t = 0.2 * 5.0 = 1.0 rad/s
    // theta_x(t) = 0.5 * alpha_x * t^2 = 0.5 * 0.2 * 25.0 = 2.5 rad
    // q(t) = [cos(theta/2), sin(theta/2), 0, 0] = [cos(1.25), sin(1.25), 0, 0]
    const double expected_omega_x = alpha_x * duration_s;
    const double expected_theta_x = 0.5 * alpha_x * duration_s * duration_s;
    const Quaternion expected_quat(std::cos(0.5 * expected_theta_x), std::sin(0.5 * expected_theta_x), 0.0, 0.0);

    CHECK_THAT(state.angular_velocity_rad_per_s.x(), WithinRel(expected_omega_x, 1.0e-11));
    CHECK_THAT(state.angular_velocity_rad_per_s.y(), WithinAbs(0.0, 1.0e-12));
    CHECK_THAT(state.angular_velocity_rad_per_s.z(), WithinAbs(0.0, 1.0e-12));

    CHECK(represents_same_rotation(state.orientation, expected_quat, 1.0e-7, 1.0e-7));

    // Rotational kinetic energy: E_rot = 0.5 * Ixx * omega_x^2 = 0.5 * 10 * 1.0 = 5.0 J
    const double energy = rotational_kinetic_energy_J(state.angular_velocity_rad_per_s, inertia);
    const double expected_energy = 0.5 * inertia.Ixx_kg_m2 * expected_omega_x * expected_omega_x;
    CHECK_THAT(energy, WithinRel(expected_energy, 1.0e-11));
}

// ===========================================================================
// Test 3: Torque-Free Principal-Axis Spin
// ===========================================================================
TEST_CASE("Torque-free principal-axis spin maintains constant rate and steady rotation", "[attitude][dynamics]") {
    const PrincipalInertia inertia(12.0, 18.0, 24.0);
    const double omega_z0 = 0.5; // rad/s
    const RotationalState initial_state{
        Quaternion::identity(),
        Vector3(0.0, 0.0, omega_z0)
    };
    const Vector3 zero_torque = Vector3{};

    // Period of one full rotation: T = 2*pi / 0.5 = 4*pi s ≈ 12.56637 s
    const double period_s = 2.0 * pi / omega_z0;
    const int num_steps = 1000;
    const double dt = period_s / num_steps;

    RotationalState state = initial_state;
    for (int step = 0; step < num_steps; ++step) {
        state = rk4_step_rotational(step * dt, state, dt, inertia, zero_torque);
    }

    // After exactly one full 360 deg rotation, attitude should return to identity (or -identity)
    CHECK_THAT(state.angular_velocity_rad_per_s.x(), WithinAbs(0.0, 1.0e-14));
    CHECK_THAT(state.angular_velocity_rad_per_s.y(), WithinAbs(0.0, 1.0e-14));
    CHECK_THAT(state.angular_velocity_rad_per_s.z(), WithinRel(omega_z0, 1.0e-12));

    CHECK(represents_same_rotation(state.orientation, Quaternion::identity(), 1.0e-10, 1.0e-10));

    // Check invariants
    const double initial_energy = rotational_kinetic_energy_J(initial_state.angular_velocity_rad_per_s, inertia);
    const double final_energy = rotational_kinetic_energy_J(state.angular_velocity_rad_per_s, inertia);
    CHECK_THAT(final_energy, WithinRel(initial_energy, 1.0e-12));
}

// ===========================================================================
// Test 4: Quaternion Kinematics Convention Guard (Multiplication Order Test)
// ===========================================================================
TEST_CASE("Quaternion derivative multiplication order matches active M08 frame convention", "[attitude][kinematics]") {
    // Under active convention v_I = q * v_B * q*, a body spin omega_z > 0 rotates
    // the spacecraft body frame counter-clockwise in inertial space.
    // Body +X axis [1, 0, 0] should rotate towards inertial +Y axis [0, 1, 0].
    const PrincipalInertia inertia(10.0, 10.0, 10.0);
    const double omega_z = 1.0; // rad/s
    const RotationalState state0{
        Quaternion::identity(),
        Vector3(0.0, 0.0, omega_z)
    };

    // Integrate for t = pi/2 s (90 deg rotation)
    const double target_time = 0.5 * pi;
    const int steps = 1000;
    const double dt = target_time / steps;

    RotationalState state = state0;
    for (int i = 0; i < steps; ++i) {
        state = rk4_step_rotational(i * dt, state, dt, inertia, Vector3{});
    }

    // Expected quaternion: 90 deg about +Z => [cos(pi/4), 0, 0, sin(pi/4)]
    const double half_sqrt2 = 0.5 * std::sqrt(2.0);
    const Quaternion expected_q(half_sqrt2, 0.0, 0.0, half_sqrt2);
    CHECK(represents_same_rotation(state.orientation, expected_q, 1.0e-10, 1.0e-10));

    // Vector transformation of body +X axis:
    const Vector3 body_x(1.0, 0.0, 0.0);
    const Vector3 inertial_vec = state.orientation.rotate_vector(body_x);
    CHECK(approximately_equal(inertial_vec, Vector3(0.0, 1.0, 0.0), 1.0e-10, 1.0e-10));

    // Negative guard: if the derivative order was reversed (dq = 0.5 * omega * q),
    // the resulting quaternion would have had opposite imaginary sign on cross-terms!
    const Quaternion reversed_derivative = (Quaternion(0.0, 0.0, 0.0, omega_z) * state0.orientation) * 0.5;
    const Quaternion nominal_derivative = quaternion_derivative(state0.orientation, state0.angular_velocity_rad_per_s);
    // For identity q, omega * q == q * omega, but for general q:
    const Quaternion q_gen = Quaternion::from_axis_angle(Vector3(1.0, 1.0, 0.0), 0.7);
    const Quaternion q_dot_nominal = quaternion_derivative(q_gen, Vector3(0.0, 0.0, omega_z));
    const Quaternion q_dot_reversed = (Quaternion(0.0, 0.0, 0.0, omega_z) * q_gen) * 0.5;
    CHECK_FALSE(approximately_equal(q_dot_nominal, q_dot_reversed, 1.0e-3, 1.0e-3));
}

// ===========================================================================
// Test 5: 90-Degree Analytical Attitude Evolution About +Y
// ===========================================================================
TEST_CASE("90-degree analytical attitude evolution about +Y axis", "[attitude][kinematics]") {
    const PrincipalInertia inertia(15.0, 15.0, 15.0);
    const double omega_y = 1.0;
    const RotationalState state0{
        Quaternion::identity(),
        Vector3(0.0, omega_y, 0.0)
    };

    const double target_time = 0.5 * pi; // 90 deg
    const int steps = 1000;
    const double dt = target_time / steps;

    RotationalState state = state0;
    for (int i = 0; i < steps; ++i) {
        state = rk4_step_rotational(i * dt, state, dt, inertia, Vector3{});
    }

    const double half_sqrt2 = 0.5 * std::sqrt(2.0);
    const Quaternion expected_q(half_sqrt2, 0.0, half_sqrt2, 0.0);
    CHECK(represents_same_rotation(state.orientation, expected_q, 1.0e-10, 1.0e-10));

    // Body +Z vector should rotate into inertial +X
    const Vector3 body_z(0.0, 0.0, 1.0);
    const Vector3 inertial_vec = state.orientation.rotate_vector(body_z);
    CHECK(approximately_equal(inertial_vec, Vector3(1.0, 0.0, 0.0), 1.0e-10, 1.0e-10));
}

// ===========================================================================
// Test 6: Quaternion Norm Behavior and Normalization Policy
// ===========================================================================
TEST_CASE("Quaternion norm drift without normalization and exact reprojection with normalization", "[attitude][numerics]") {
    const PrincipalInertia inertia(10.0, 15.0, 25.0);
    const Vector3 omega0(0.3, 0.4, 0.5);
    const RotationalState initial_state{
        Quaternion::identity(),
        omega0
    };

    const double dt = 0.05;
    const int total_steps = 1000;

    // Run 1: Without normalization (renormalize_quaternion = false)
    RotationalState unnormalized_state = initial_state;
    double max_unnormalized_drift = 0.0;
    for (int step = 0; step < total_steps; ++step) {
        unnormalized_state = rk4_step_rotational(
            step * dt, unnormalized_state, dt, inertia, Vector3{}, false
        );
        const double norm_err = std::abs(unnormalized_state.orientation.norm() - 1.0);
        if (norm_err > max_unnormalized_drift) {
            max_unnormalized_drift = norm_err;
        }
    }

    // Run 2: With normalization (renormalize_quaternion = true, default)
    RotationalState normalized_state = initial_state;
    double max_normalized_drift = 0.0;
    for (int step = 0; step < total_steps; ++step) {
        normalized_state = rk4_step_rotational(
            step * dt, normalized_state, dt, inertia, Vector3{}, true
        );
        const double norm_err = std::abs(normalized_state.orientation.norm() - 1.0);
        if (norm_err > max_normalized_drift) {
            max_normalized_drift = norm_err;
        }
    }

    // Unnormalized integration exhibits observable non-zero norm drift
    CHECK(max_unnormalized_drift > 1.0e-12);
    // Normalized integration maintains unit norm to machine precision
    CHECK(max_normalized_drift < 1.0e-15);
}

// ===========================================================================
// Test 7: Asymmetric Torque-Free Motion and Gyroscopic Coupling
// ===========================================================================
TEST_CASE("Asymmetric torque-free motion exhibits gyroscopic coupling and non-trivial rate variation", "[attitude][dynamics]") {
    const PrincipalInertia inertia(10.0, 20.0, 30.0); // Ixx != Iyy != Izz
    const Vector3 omega0(0.2, 0.3, 0.1);
    const RotationalState state0{
        Quaternion::identity(),
        omega0
    };

    // Calculate instantaneous acceleration at t=0
    // d(omega_x)/dt = -(30 - 20) * 0.3 * 0.1 / 10 = -10 * 0.03 / 10 = -0.03 rad/s^2
    // d(omega_y)/dt = -(10 - 30) * 0.1 * 0.2 / 20 = -(-20) * 0.02 / 20 = +0.02 rad/s^2
    // d(omega_z)/dt = -(20 - 10) * 0.2 * 0.3 / 30 = -10 * 0.06 / 30 = -0.02 rad/s^2
    const Vector3 alpha0 = euler_rotational_acceleration(omega0, inertia, Vector3{});
    CHECK_THAT(alpha0.x(), WithinRel(-0.03, 1.0e-14));
    CHECK_THAT(alpha0.y(), WithinRel(0.02, 1.0e-14));
    CHECK_THAT(alpha0.z(), WithinRel(-0.02, 1.0e-14));

    // Propagate for 50 seconds
    const double dt = 0.01;
    const int steps = 5000;
    RotationalState state = state0;
    for (int i = 0; i < steps; ++i) {
        state = rk4_step_rotational(i * dt, state, dt, inertia, Vector3{});
    }

    // Rates should have evolved substantially from initial values
    CHECK_FALSE(approximately_equal(state.angular_velocity_rad_per_s, omega0, 1.0e-3, 1.0e-3));
}

// ===========================================================================
// Test 8: Rotational Kinetic Energy Conservation
// ===========================================================================
TEST_CASE("Torque-free asymmetric tumbling conserves rotational kinetic energy", "[attitude][invariants]") {
    const PrincipalInertia inertia(14.0, 22.0, 35.0);
    const Vector3 omega0(0.4, -0.2, 0.3);
    const RotationalState initial_state{
        Quaternion::from_axis_angle(Vector3(1.0, 2.0, 3.0), 0.5),
        omega0
    };

    const double initial_energy = rotational_kinetic_energy_J(omega0, inertia);
    const double dt = 0.01;
    const int steps = 5000; // 50 seconds

    RotationalState state = initial_state;
    double max_rel_energy_err = 0.0;
    for (int i = 0; i < steps; ++i) {
        state = rk4_step_rotational(i * dt, state, dt, inertia, Vector3{});
        const double energy = rotational_kinetic_energy_J(state.angular_velocity_rad_per_s, inertia);
        const double rel_err = std::abs(energy - initial_energy) / initial_energy;
        if (rel_err > max_rel_energy_err) {
            max_rel_energy_err = rel_err;
        }
    }

    // RK4 conservation residual across 50 s
    CHECK(max_rel_energy_err < 1.0e-10);
}

// ===========================================================================
// Test 9: Inertial Angular Momentum Conservation
// ===========================================================================
TEST_CASE("Torque-free motion conserves the inertial angular momentum vector", "[attitude][invariants]") {
    const PrincipalInertia inertia(15.0, 25.0, 32.0);
    const Vector3 omega0(0.3, 0.5, -0.4);
    const RotationalState initial_state{
        Quaternion::from_axis_angle(Vector3(-1.0, 1.0, 2.0), 1.1),
        omega0
    };

    const Vector3 h_inertial_0 = inertial_angular_momentum_kg_m2_per_s(initial_state, inertia);
    const double h_norm_0 = h_inertial_0.norm();
    const double dt = 0.01;
    const int steps = 5000;

    RotationalState state = initial_state;
    double max_h_rel_error = 0.0;
    for (int i = 0; i < steps; ++i) {
        state = rk4_step_rotational(i * dt, state, dt, inertia, Vector3{});
        const Vector3 h_inertial = inertial_angular_momentum_kg_m2_per_s(state, inertia);
        const double vec_err = (h_inertial - h_inertial_0).norm() / h_norm_0;
        if (vec_err > max_h_rel_error) {
            max_h_rel_error = vec_err;
        }
    }

    // Inertial angular momentum conserved to high precision
    CHECK(max_h_rel_error < 1.0e-10);
}

// ===========================================================================
// Test 10: Body vs Inertial Angular Momentum Distinction
// ===========================================================================
TEST_CASE("Body angular momentum components oscillate while inertial angular momentum is stationary", "[attitude][invariants]") {
    const PrincipalInertia inertia(10.0, 20.0, 35.0);
    const Vector3 omega0(0.5, 0.4, 0.3);
    const RotationalState initial_state{
        Quaternion::identity(),
        omega0
    };

    const Vector3 h_body_0 = body_angular_momentum_kg_m2_per_s(omega0, inertia);
    const Vector3 h_inertial_0 = inertial_angular_momentum_kg_m2_per_s(initial_state, inertia);

    const double dt = 0.01;
    const int steps = 2000; // 20 s

    RotationalState state = initial_state;
    double max_body_h_deviation = 0.0;
    double max_inertial_h_deviation = 0.0;

    for (int i = 0; i < steps; ++i) {
        state = rk4_step_rotational(i * dt, state, dt, inertia, Vector3{});

        const Vector3 h_body = body_angular_momentum_kg_m2_per_s(state.angular_velocity_rad_per_s, inertia);
        const Vector3 h_inertial = inertial_angular_momentum_kg_m2_per_s(state, inertia);

        const double body_dev = (h_body - h_body_0).norm();
        const double inertial_dev = (h_inertial - h_inertial_0).norm();

        if (body_dev > max_body_h_deviation) max_body_h_deviation = body_dev;
        if (inertial_dev > max_inertial_h_deviation) max_inertial_h_deviation = inertial_dev;
    }

    // Body-frame components change significantly (due to tumbling frame)
    CHECK(max_body_h_deviation > 2.0); // Substantial variation in body components
    // Inertial-frame components remain constant
    CHECK(max_inertial_h_deviation < 1.0e-9);
}

// ===========================================================================
// Test 11: Constant-Torque Work-Energy Relation
// ===========================================================================
TEST_CASE("External torque changes rotational kinetic energy matching mechanical work", "[attitude][dynamics]") {
    const PrincipalInertia inertia(10.0, 15.0, 20.0);
    const Vector3 torque(1.5, -0.8, 0.4);
    const RotationalState state0{
        Quaternion::identity(),
        Vector3(0.1, 0.2, -0.1)
    };

    const double initial_energy = rotational_kinetic_energy_J(state0.angular_velocity_rad_per_s, inertia);
    const double dt = 0.005;
    const int steps = 1000; // 5.0 s

    RotationalState state = state0;
    double accumulated_work_J = 0.0;

    for (int i = 0; i < steps; ++i) {
        const double t = i * dt;
        // Trapezoidal integration of power P = tau . omega
        const double power_start = torque.dot(state.angular_velocity_rad_per_s);
        state = rk4_step_rotational(t, state, dt, inertia, torque);
        const double power_end = torque.dot(state.angular_velocity_rad_per_s);

        accumulated_work_J += 0.5 * (power_start + power_end) * dt;
    }

    const double final_energy = rotational_kinetic_energy_J(state.angular_velocity_rad_per_s, inertia);
    const double delta_energy = final_energy - initial_energy;

    CHECK_THAT(delta_energy, WithinRel(accumulated_work_J, 1.0e-4));
}

// ===========================================================================
// Test 12: Timestep Convergence Characterization (RK4 4th-Order)
// ===========================================================================
TEST_CASE("RK4 rotational dynamics demonstrates fourth-order timestep convergence", "[attitude][numerics]") {
    const PrincipalInertia inertia(12.0, 18.0, 28.0);
    const Vector3 omega0(0.2, 0.3, 0.4);
    const RotationalState state0{
        Quaternion::identity(),
        omega0
    };

    const double total_time = 2.0;

    // High-resolution reference solution (dt = 0.0001 s)
    const double dt_ref = 0.0001;
    const int steps_ref = static_cast<int>(total_time / dt_ref);
    RotationalState state_ref = state0;
    for (int i = 0; i < steps_ref; ++i) {
        state_ref = rk4_step_rotational(i * dt_ref, state_ref, dt_ref, inertia, Vector3{}, false);
    }

    // Coarse timesteps: dt1 = 0.1 s, dt2 = 0.05 s, dt3 = 0.025 s
    auto run_sim = [&](double dt) -> RotationalState {
        const int n = static_cast<int>(total_time / dt);
        RotationalState s = state0;
        for (int i = 0; i < n; ++i) {
            s = rk4_step_rotational(i * dt, s, dt, inertia, Vector3{}, false);
        }
        return s;
    };

    const RotationalState s1 = run_sim(0.1);
    const RotationalState s2 = run_sim(0.05);
    const RotationalState s3 = run_sim(0.025);

    const double err1 = (s1.angular_velocity_rad_per_s - state_ref.angular_velocity_rad_per_s).norm();
    const double err2 = (s2.angular_velocity_rad_per_s - state_ref.angular_velocity_rad_per_s).norm();
    const double err3 = (s3.angular_velocity_rad_per_s - state_ref.angular_velocity_rad_per_s).norm();

    const double ratio12 = err1 / err2;
    const double ratio23 = err2 / err3;

    // For 4th order, halving dt reduces error by factor of ~16 (2^4)
    CHECK(ratio12 > 14.0);
    CHECK(ratio12 < 17.5);
    CHECK(ratio23 > 14.0);
    CHECK(ratio23 < 17.5);
}

// ===========================================================================
// Test 13: Deterministic Repeatability
// ===========================================================================
TEST_CASE("Simulations with identical parameters produce deterministic identical results", "[attitude][determinism]") {
    const PrincipalInertia inertia(11.0, 22.0, 33.0);
    const RotationalState initial_state{
        Quaternion::from_axis_angle(Vector3(1.0, 1.0, 1.0), 0.8),
        Vector3(0.5, -0.3, 0.2)
    };
    const Vector3 torque(0.1, -0.2, 0.05);

    const double dt = 0.01;
    const int steps = 500;

    RotationalState run1 = initial_state;
    for (int i = 0; i < steps; ++i) {
        run1 = rk4_step_rotational(i * dt, run1, dt, inertia, torque);
    }

    RotationalState run2 = initial_state;
    for (int i = 0; i < steps; ++i) {
        run2 = rk4_step_rotational(i * dt, run2, dt, inertia, torque);
    }

    CHECK(run1.orientation.w() == run2.orientation.w());
    CHECK(run1.orientation.x() == run2.orientation.x());
    CHECK(run1.orientation.y() == run2.orientation.y());
    CHECK(run1.orientation.z() == run2.orientation.z());
    CHECK(run1.angular_velocity_rad_per_s.x() == run2.angular_velocity_rad_per_s.x());
    CHECK(run1.angular_velocity_rad_per_s.y() == run2.angular_velocity_rad_per_s.y());
    CHECK(run1.angular_velocity_rad_per_s.z() == run2.angular_velocity_rad_per_s.z());
}

// ===========================================================================
// Test 14: Defensive Validation Rejects Invalid Inputs
// ===========================================================================
TEST_CASE("Defensive validation rejects non-finite, zero, and negative parameters", "[attitude][defensive]") {
    // 1. Zero and negative inertia
    CHECK_THROWS_AS(PrincipalInertia(0.0, 10.0, 10.0), std::domain_error);
    CHECK_THROWS_AS(PrincipalInertia(10.0, -5.0, 10.0), std::domain_error);
    CHECK_THROWS_AS(PrincipalInertia(10.0, 10.0, std::numeric_limits<double>::quiet_NaN()), std::domain_error);
    CHECK_THROWS_AS(PrincipalInertia(10.0, 10.0, std::numeric_limits<double>::infinity()), std::domain_error);

    const PrincipalInertia valid_inertia(10.0, 20.0, 30.0);

    // 2. Non-finite angular velocity in derivative / acceleration
    const Vector3 nan_vec(std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0);
    const Vector3 inf_vec(0.0, std::numeric_limits<double>::infinity(), 0.0);
    CHECK_THROWS_AS(quaternion_derivative(Quaternion::identity(), nan_vec), std::domain_error);
    CHECK_THROWS_AS(euler_rotational_acceleration(nan_vec, valid_inertia, Vector3{}), std::domain_error);
    CHECK_THROWS_AS(euler_rotational_acceleration(Vector3{}, valid_inertia, inf_vec), std::domain_error);
    CHECK_THROWS_AS(rotational_kinetic_energy_J(nan_vec, valid_inertia), std::domain_error);
    CHECK_THROWS_AS(body_angular_momentum_kg_m2_per_s(inf_vec, valid_inertia), std::domain_error);

    // 3. Non-finite quaternion in state derivative
    const RotationalState nan_state{
        Quaternion(std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0, 0.0),
        Vector3{}
    };
    CHECK_THROWS_AS(rotational_state_derivative(0.0, nan_state, valid_inertia, Vector3{}), std::domain_error);
}
