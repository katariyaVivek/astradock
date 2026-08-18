#include "dynamics/two_body.hpp"
#include "frames/frame_basis.hpp"
#include "frames/lvlh.hpp"
#include "math/constants.hpp"
#include "math/matrix3.hpp"
#include "math/vector3.hpp"
#include "numerics/fixed_step_propagation.hpp"
#include "numerics/integrators.hpp"
#include "orbit/cartesian_state.hpp"
#include "orbit/orbital_diagnostics.hpp"
#include "orbit/two_body_orbit.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

using astradock::dynamics::two_body_acceleration;
using astradock::frames::compute_lvlh_basis;
using astradock::frames::dcm_eci_from_lvlh;
using astradock::frames::dcm_lvlh_from_eci;
using astradock::frames::FrameBasis;
using astradock::frames::is_orthonormal;
using astradock::frames::transform_eci_to_lvlh;
using astradock::frames::transform_lvlh_to_eci;
using astradock::math::approximately_equal;
using astradock::math::is_finite;
using astradock::math::Matrix3;
using astradock::math::Vector3;
using astradock::numerics::euler_step;
using astradock::numerics::IntegrationMethod;
using astradock::numerics::propagate_fixed_step;
using astradock::numerics::rk4_step;
using astradock::orbit::CartesianState;
using astradock::orbit::compute_circular_orbit_reference;
using astradock::orbit::specific_angular_momentum_m2_per_s;
using astradock::orbit::specific_orbital_energy_m2_per_s2;
using astradock::orbit::two_body_state_derivative;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// ===========================================================================
// M01 — Vector Mathematics Independent & Property Tests
// ===========================================================================

TEST_CASE("M01 Audit: Vector3 deterministic randomized property testing", "[audit][m01]") {
    // Fixed seed PRNG for strictly reproducible property testing
    std::mt19937_64 rng(424242ULL);
    std::uniform_real_distribution<double> dist(-1000.0, 1000.0);

    for (std::size_t sample = 0; sample < 100; ++sample) {
        const Vector3 a(dist(rng), dist(rng), dist(rng));
        const Vector3 b(dist(rng), dist(rng), dist(rng));
        const Vector3 c(dist(rng), dist(rng), dist(rng));

        // 1. Commutativity of addition: a + b == b + a
        CHECK(approximately_equal(a + b, b + a, 1.0e-13, 1.0e-13));

        // 2. Associativity of addition: (a + b) + c == a + (b + c)
        CHECK(approximately_equal((a + b) + c, a + (b + c), 1.0e-12, 1.0e-12));

        // 3. Commutativity of dot product: a · b == b · a
        CHECK_THAT(a.dot(b), WithinAbs(b.dot(a), 1.0e-11));

        // 4. Anti-commutativity of cross product: a x b == -(b x a)
        const Vector3 axb = a.cross(b);
        const Vector3 bxa = b.cross(a);
        CHECK(approximately_equal(axb, bxa * -1.0, 1.0e-11, 1.0e-11));

        // 5. Self cross product is zero: a x a == 0
        const Vector3 axa = a.cross(a);
        const double a_mag = a.norm();
        CHECK(axa.norm() <= 1.0e-13 * a_mag * a_mag + 1.0e-15);

        // 6. Cross product orthogonality: (a x b) · a == 0 and (a x b) · b == 0
        // Scaled by ||a x b|| * ||a||
        const double scale_a = axb.norm() * a.norm();
        const double scale_b = axb.norm() * b.norm();
        if (scale_a > 1.0e-6) {
            CHECK_THAT(std::abs(axb.dot(a)) / scale_a, WithinAbs(0.0, 1.0e-13));
        }
        if (scale_b > 1.0e-6) {
            CHECK_THAT(std::abs(axb.dot(b)) / scale_b, WithinAbs(0.0, 1.0e-13));
        }

        // 7. Lagrange identity: ||a x b||^2 + (a · b)^2 == ||a||^2 ||b||^2
        const double lhs = axb.squared_norm() + (a.dot(b) * a.dot(b));
        const double rhs = a.squared_norm() * b.squared_norm();
        CHECK_THAT(lhs, WithinRel(rhs, 1.0e-11));

        // 8. Norm squared identity: ||a||^2 == a · a
        CHECK_THAT(a.norm() * a.norm(), WithinRel(a.squared_norm(), 1.0e-12));

        // 9. Cauchy-Schwarz inequality: |a · b| <= ||a|| ||b||
        CHECK(std::abs(a.dot(b)) <= a.norm() * b.norm() + 1.0e-11);

        // 10. Triangle inequality: ||a + b|| <= ||a|| + ||b||
        CHECK((a + b).norm() <= a.norm() + b.norm() + 1.0e-11);

        // 11. Normalization preserves direction and achieves unit norm
        if (a.norm() > 1.0e-6) {
            const Vector3 a_hat = a.normalized();
            CHECK_THAT(a_hat.norm(), WithinAbs(1.0, 1.0e-12));
            CHECK_THAT(a_hat.dot(a), WithinRel(a.norm(), 1.0e-12));
        }
    }
}

// ===========================================================================
// M02 — Two-Body Gravity Independent Verification
// ===========================================================================

TEST_CASE("M02 Audit: Two-body gravity independent scaling and rotational invariance", "[audit][m02]") {
    const double mu = 3.986004418e14;
    const double r0_m = 6878137.0;

    // 1. Independent magnitude check: |a| = mu / r^2
    const Vector3 r_vec(r0_m, 0.0, 0.0);
    const Vector3 a_vec = two_body_acceleration(r_vec, mu);
    const double expected_mag = mu / (r0_m * r0_m);

    CHECK_THAT(a_vec.norm(), WithinRel(expected_mag, 1.0e-14));
    CHECK_THAT(a_vec.x(), WithinRel(-expected_mag, 1.0e-14));
    CHECK_THAT(a_vec.y(), WithinAbs(0.0, 1.0e-14));
    CHECK_THAT(a_vec.z(), WithinAbs(0.0, 1.0e-14));

    // 2. Inverse-square scaling: a(alpha * r) == a(r) / alpha^2
    for (double alpha : {1.5, 2.0, 3.0, 10.0}) {
        const Vector3 r_scaled = r_vec * alpha;
        const Vector3 a_scaled = two_body_acceleration(r_scaled, mu);
        const double expected_scaled_mag = expected_mag / (alpha * alpha);
        CHECK_THAT(a_scaled.norm(), WithinRel(expected_scaled_mag, 1.0e-13));
    }

    // 3. Rotational symmetry: rotating r by 90 degrees rotates a identically
    // Rotation about Z by 90 deg: [x, y, z] -> [-y, x, z]
    const Vector3 r_rot(-r_vec.y(), r_vec.x(), r_vec.z());
    const Vector3 a_rot = two_body_acceleration(r_rot, mu);
    CHECK_THAT(a_rot.x(), WithinAbs(0.0, 1.0e-14));
    CHECK_THAT(a_rot.y(), WithinRel(-expected_mag, 1.0e-14));
    CHECK_THAT(a_rot.z(), WithinAbs(0.0, 1.0e-14));

    // 4. Direction exactness: a is strictly anti-parallel to r
    const Vector3 arbitrary_r(3000000.0, -4000000.0, 5000000.0);
    const Vector3 arbitrary_a = two_body_acceleration(arbitrary_r, mu);
    // a · r == -||a|| ||r||
    CHECK_THAT(arbitrary_a.dot(arbitrary_r), WithinRel(-arbitrary_a.norm() * arbitrary_r.norm(), 1.0e-14));
    // a x r == 0
    CHECK_THAT(arbitrary_a.cross(arbitrary_r).norm(), WithinAbs(0.0, 1.0e-8));
}

// ===========================================================================
// M03 — Numerical Integration Independent Verification
// ===========================================================================

TEST_CASE("M03 Audit: Euler on independent linear ODEs", "[audit][m03]") {
    // ODE 1: y' = 1, y(0) = 5 -> y(t) = 5 + t
    double y = 5.0;
    const double dt = 0.1;
    for (int step = 0; step < 10; ++step) {
        y = euler_step(step * dt, y, dt, [](double, double) { return 1.0; });
    }
    CHECK_THAT(y, WithinAbs(6.0, 1.0e-14));

    // ODE 2: y' = -y, y(0) = 1 -> y(t) = e^(-t)
    // 1-step Euler at dt = 0.1: y_1 = 1 + 0.1 * (-1) = 0.9
    const double y1 = euler_step(0.0, 1.0, 0.1, [](double, double v) { return -v; });
    CHECK_THAT(y1, WithinAbs(0.9, 1.0e-14));
}

TEST_CASE("M03 Audit: RK4 hand-derived cubic polynomial exactness", "[audit][m03]") {
    // For y' = t^3, y(0) = 0, the analytical integral is y(h) = h^4 / 4.
    // RK4 integrates any polynomial up to degree 3 with zero truncation error.
    const double h = 2.0;
    const double y_next = rk4_step(
        0.0, 0.0, h,
        [](double t, double) { return t * t * t; });

    const double expected_hand_calc = (h * h * h * h) / 4.0; // 16 / 4 = 4.0
    CHECK_THAT(y_next, WithinAbs(expected_hand_calc, 1.0e-14));
}

TEST_CASE("M03 Audit: RK4 convergence order on exponential decay", "[audit][m03]") {
    // y' = -y, y(0) = 1, exact solution y(1) = 1/e
    const double t_end = 1.0;
    const double exact_final = std::exp(-1.0);

    auto run_rk4 = [&](double dt) {
        double y = 1.0;
        double t = 0.0;
        while (t < t_end - 1.0e-12) {
            const double step = std::min(dt, t_end - t);
            y = rk4_step(t, y, step, [](double, double v) { return -v; });
            t += step;
        }
        return std::abs(y - exact_final);
    };

    const double e_h1 = run_rk4(0.05);
    const double e_h2 = run_rk4(0.025);
    const double empirical_order = std::log(e_h1 / e_h2) / std::log(0.05 / 0.025);

    // Theoretical asymptotic order is 4.0
    CHECK_THAT(empirical_order, WithinAbs(4.0, 0.05));
}

// ===========================================================================
// M04 & M05 — Orbital Dynamics, Invariants, Multi-Orbit Stability & Determinism
// ===========================================================================

TEST_CASE("M04 & M05 Audit: Hand-derived circular orbit analytical reference values", "[audit][m04][m05]") {
    const double mu = 3.986004418e14;
    const double r_orbit = 6378137.0 + 500000.0; // 6,878,137.0 m

    // Hand-derived values computed strictly from fundamental definitions:
    const double expected_vc = std::sqrt(mu / r_orbit);               // 7612.6081729738...
    const double expected_period = 2.0 * astradock::constants::pi * std::sqrt((r_orbit * r_orbit * r_orbit) / mu); // 5676.97802927...
    const double expected_energy = -mu / (2.0 * r_orbit);              // -28975901.599517...
    const double expected_h = r_orbit * expected_vc;                  // 52360561942.7535...

    const auto ref = compute_circular_orbit_reference(mu, r_orbit);

    CHECK_THAT(ref.speed_m_per_s, WithinRel(expected_vc, 1.0e-14));
    CHECK_THAT(ref.period_s, WithinRel(expected_period, 1.0e-14));
    CHECK_THAT(ref.specific_energy_m2_per_s2, WithinRel(expected_energy, 1.0e-14));
    CHECK_THAT(ref.specific_angular_momentum_m2_per_s, WithinRel(expected_h, 1.0e-14));
}

TEST_CASE("M05 Audit: Multi-orbit stability (1, 5, 10 orbits) and deterministic repeatability", "[audit][m05]") {
    const double mu = astradock::constants::earth_gravitational_parameter_m3_per_s2;
    const double r_orbit = astradock::constants::earth_reference_radius_m + 500000.0;
    const auto ref = compute_circular_orbit_reference(mu, r_orbit);

    const CartesianState state0{
        Vector3(r_orbit, 0.0, 0.0),
        Vector3(0.0, ref.speed_m_per_s, 0.0)
    };

    const double dt = 10.0;
    const double duration_10_orbits = 10.0 * ref.period_s;

    // Run 1: 10 orbits
    const auto traj1 = propagate_fixed_step(
        0.0, duration_10_orbits, dt, state0,
        IntegrationMethod::classical_rk4,
        [mu](double t, const CartesianState& s) {
            return two_body_state_derivative(t, s, mu);
        });

    // Run 2: Repeat exact same run for determinism check
    const auto traj2 = propagate_fixed_step(
        0.0, duration_10_orbits, dt, state0,
        IntegrationMethod::classical_rk4,
        [mu](double t, const CartesianState& s) {
            return two_body_state_derivative(t, s, mu);
        });

    REQUIRE(traj1.size() == traj2.size());
    REQUIRE(traj1.size() > 5000);

    // Verify bitwise determinism between repeated runs
    for (std::size_t i = 0; i < traj1.size(); ++i) {
        CHECK(traj1[i].time_s == traj2[i].time_s);
        CHECK(traj1[i].state.position.x() == traj2[i].state.position.x());
        CHECK(traj1[i].state.position.y() == traj2[i].state.position.y());
        CHECK(traj1[i].state.position.z() == traj2[i].state.position.z());
        CHECK(traj1[i].state.velocity.x() == traj2[i].state.velocity.x());
        CHECK(traj1[i].state.velocity.y() == traj2[i].state.velocity.y());
        CHECK(traj1[i].state.velocity.z() == traj2[i].state.velocity.z());
    }

    // Inspect 10-orbit long-term behavior
    double max_energy_drift = 0.0;
    double max_h_drift = 0.0;
    double max_radius_dev = 0.0;

    for (const auto& sample : traj1) {
        const double energy = specific_orbital_energy_m2_per_s2(sample.state, mu);
        const double rel_energy_err = std::abs(energy - ref.specific_energy_m2_per_s2) / std::abs(ref.specific_energy_m2_per_s2);
        max_energy_drift = std::max(max_energy_drift, rel_energy_err);

        const double h_mag = specific_angular_momentum_m2_per_s(sample.state).norm();
        const double rel_h_err = std::abs(h_mag - ref.specific_angular_momentum_m2_per_s) / ref.specific_angular_momentum_m2_per_s;
        max_h_drift = std::max(max_h_drift, rel_h_err);

        const double r_dev = std::abs(sample.state.position.norm() - r_orbit);
        max_radius_dev = std::max(max_radius_dev, r_dev);
    }

    // Over 10 orbits (~15.7 hours, ~56,770 seconds):
    // Relative energy drift should remain < 5.0e-10
    CHECK(max_energy_drift < 5.0e-10);
    // Relative angular momentum drift should remain < 5.0e-10
    CHECK(max_h_drift < 5.0e-10);
    // Radial deviation should remain bounded within 0.05 m
    CHECK(max_radius_dev < 0.05);
}

// ===========================================================================
// M06 — Matrix3, DCM & Coordinate Frames Independent Verification
// ===========================================================================

TEST_CASE("M06 Audit: Matrix3 algebraic properties and independent determinant", "[audit][m06]") {
    // 1. Matrix multiplication associativity: (A * B) * C == A * (B * C)
    const Matrix3 A(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0);
    const Matrix3 B(2.0, 0.0, 1.0, 3.0, 1.0, 2.0, 1.0, 4.0, 1.0);
    const Matrix3 C(1.0, 1.0, 0.0, 0.0, 2.0, 1.0, 3.0, 0.0, 1.0);

    const Matrix3 AB_C = (A * B) * C;
    const Matrix3 A_BC = A * (B * C);
    CHECK(approximately_equal(AB_C, A_BC, 1.0e-11, 1.0e-11));

    // 2. Determinant check for upper-triangular matrix: det(T) = t00 * t11 * t22
    const Matrix3 T(3.0, 4.0, 5.0, 0.0, 2.0, 6.0, 0.0, 0.0, 7.0);
    CHECK_THAT(T.determinant(), WithinAbs(3.0 * 2.0 * 7.0, 1.0e-14)); // 42.0

    // 3. Matrix inversion via transpose for proper rotation: C * C^T == I
    const double angle = 0.5235987755982988; // 30 deg in rad
    const Matrix3 rot_x(
        1.0, 0.0, 0.0,
        0.0, std::cos(angle), -std::sin(angle),
        0.0, std::sin(angle), std::cos(angle)
    );
    CHECK_THAT(rot_x.determinant(), WithinAbs(1.0, 1.0e-14));
    CHECK(approximately_equal(rot_x * rot_x.transpose(), Matrix3::identity(), 1.0e-14, 1.0e-14));
}

TEST_CASE("M06 Audit: DCM explicit row/column and vector transformation manual check", "[audit][m06]") {
    // Hand-derived state: Spacecraft at 45 degrees in XY plane
    // r = [R/sqrt(2), R/sqrt(2), 0]
    // v = [-V/sqrt(2), V/sqrt(2), 0]
    const double R = 7000000.0;
    const double V = 7500.0;
    const double inv_sqrt2 = 1.0 / std::sqrt(2.0);

    const Vector3 r_eci(R * inv_sqrt2, R * inv_sqrt2, 0.0);
    const Vector3 v_eci(-V * inv_sqrt2, V * inv_sqrt2, 0.0);

    const FrameBasis basis = compute_lvlh_basis(r_eci, v_eci);

    // Hand-derived basis vectors:
    // e_r = [1/sqrt(2), 1/sqrt(2), 0]
    // e_h = [0, 0, 1] (r x v is along +Z)
    // e_t = e_h x e_r = [-1/sqrt(2), 1/sqrt(2), 0]
    CHECK_THAT(basis.x.x(), WithinAbs(inv_sqrt2, 1.0e-14));
    CHECK_THAT(basis.x.y(), WithinAbs(inv_sqrt2, 1.0e-14));
    CHECK_THAT(basis.x.z(), WithinAbs(0.0, 1.0e-14));

    CHECK_THAT(basis.y.x(), WithinAbs(-inv_sqrt2, 1.0e-14));
    CHECK_THAT(basis.y.y(), WithinAbs(inv_sqrt2, 1.0e-14));
    CHECK_THAT(basis.y.z(), WithinAbs(0.0, 1.0e-14));

    CHECK_THAT(basis.z.x(), WithinAbs(0.0, 1.0e-14));
    CHECK_THAT(basis.z.y(), WithinAbs(0.0, 1.0e-14));
    CHECK_THAT(basis.z.z(), WithinAbs(1.0, 1.0e-14));

    // Verify DCM transformation:
    // C_LVLH_ECI * r_eci should equal [R, 0, 0] exactly
    const Matrix3 c_lvlh_eci = dcm_lvlh_from_eci(basis);
    const Vector3 r_lvlh = transform_eci_to_lvlh(r_eci, c_lvlh_eci);
    CHECK_THAT(r_lvlh.x(), WithinAbs(R, 1.0e-8));
    CHECK_THAT(r_lvlh.y(), WithinAbs(0.0, 1.0e-8));
    CHECK_THAT(r_lvlh.z(), WithinAbs(0.0, 1.0e-8));

    // C_LVLH_ECI * v_eci should equal [0, V, 0] exactly
    const Vector3 v_lvlh = transform_eci_to_lvlh(v_eci, c_lvlh_eci);
    CHECK_THAT(v_lvlh.x(), WithinAbs(0.0, 1.0e-8));
    CHECK_THAT(v_lvlh.y(), WithinAbs(V, 1.0e-8));
    CHECK_THAT(v_lvlh.z(), WithinAbs(0.0, 1.0e-8));

    // Inverse transform: C_ECI_LVLH * r_lvlh == r_eci
    const Vector3 r_reconstructed = transform_lvlh_to_eci(r_lvlh, c_lvlh_eci);
    CHECK(approximately_equal(r_reconstructed, r_eci, 1.0e-8, 1.0e-8));
}
