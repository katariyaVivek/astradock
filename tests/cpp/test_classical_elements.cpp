#include "dynamics/two_body.hpp"
#include "math/angle.hpp"
#include "math/constants.hpp"
#include "math/vector3.hpp"
#include "numerics/fixed_step_propagation.hpp"
#include "orbit/cartesian_state.hpp"
#include "orbit/classical_elements.hpp"
#include "orbit/two_body_orbit.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

using astradock::math::angular_distance_rad;
using astradock::math::approximately_equal;
using astradock::constants::pi;
using astradock::math::normalize_angle_2pi_rad;
using astradock::math::normalize_angle_pi_rad;
using astradock::math::Vector3;
using astradock::numerics::IntegrationMethod;
using astradock::numerics::propagate_fixed_step;
using astradock::orbit::apoapsis_radius_m;
using astradock::orbit::CartesianState;
using astradock::orbit::classical_elements_to_state;
using astradock::orbit::ClassicalOrbitalElements;
using astradock::orbit::classify_orbit;
using astradock::orbit::mean_motion_rad_per_s;
using astradock::orbit::OrbitInclination;
using astradock::orbit::OrbitShape;
using astradock::orbit::periapsis_radius_m;
using astradock::orbit::semi_latus_rectum_m;
using astradock::orbit::state_to_classical_elements;
using astradock::orbit::two_body_state_derivative;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// ===========================================================================
// Angle Normalization and Distance Tests
// ===========================================================================

TEST_CASE("Angle normalization and wrapping handle arbitrary ranges", "[elements][angles]") {
    // 0 deg and 360 deg
    CHECK_THAT(normalize_angle_2pi_rad(0.0), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(normalize_angle_2pi_rad(2.0 * pi), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(normalize_angle_2pi_rad(4.0 * pi), WithinAbs(0.0, 1.0e-15));

    // Negative angles (-10 deg -> 350 deg)
    const double minus_10_deg = -10.0 * pi / 180.0;
    const double expected_350_deg = 350.0 * pi / 180.0;
    CHECK_THAT(normalize_angle_2pi_rad(minus_10_deg), WithinRel(expected_350_deg, 1.0e-14));

    // Over-wrapped angles (370 deg -> 10 deg)
    const double plus_370_deg = 370.0 * pi / 180.0;
    const double expected_10_deg = 10.0 * pi / 180.0;
    CHECK_THAT(normalize_angle_2pi_rad(plus_370_deg), WithinRel(expected_10_deg, 1.0e-14));

    // Symmetric interval [-pi, pi)
    CHECK_THAT(normalize_angle_pi_rad(0.0), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(normalize_angle_pi_rad(pi * 0.5), WithinAbs(pi * 0.5, 1.0e-15));
    CHECK_THAT(normalize_angle_pi_rad(350.0 * pi / 180.0), WithinRel(-10.0 * pi / 180.0, 1.0e-14));

    // Angular distance (359.999 deg and 0.001 deg are 0.002 deg apart)
    const double a1 = 359.999 * pi / 180.0;
    const double a2 = 0.001 * pi / 180.0;
    CHECK_THAT(angular_distance_rad(a1, a2), WithinRel(0.002 * pi / 180.0, 1.0e-9));
}

// ===========================================================================
// Analytical Auxiliary Functions Tests
// ===========================================================================

TEST_CASE("Analytical orbital geometry helpers satisfy mathematical definitions", "[elements][geometry]") {
    const double a = 10000000.0; // 10,000 km
    const double e = 0.2;
    const double mu = 3.986004418e14;

    // Semi-latus rectum: p = a * (1 - e^2)
    const double p = semi_latus_rectum_m(a, e);
    CHECK_THAT(p, WithinAbs(10000000.0 * (1.0 - 0.04), 1.0e-6)); // 9,600,000 m

    // Periapsis & Apoapsis: rp = a*(1-e), ra = a*(1+e)
    const double rp = periapsis_radius_m(a, e);
    const double ra = apoapsis_radius_m(a, e);
    CHECK_THAT(rp, WithinAbs(8000000.0, 1.0e-6));
    CHECK_THAT(ra, WithinAbs(12000000.0, 1.0e-6));
    CHECK_THAT((rp + ra) * 0.5, WithinAbs(a, 1.0e-6));

    // Mean motion: n = sqrt(mu / a^3)
    const double n = mean_motion_rad_per_s(a, mu);
    const double expected_n = std::sqrt(mu / (a * a * a));
    CHECK_THAT(n, WithinRel(expected_n, 1.0e-14));

    // Orbital period: T = 2*pi / n
    const double T = astradock::orbit::orbital_period_s(a, mu);
    const double expected_T = 2.0 * pi / expected_n;
    CHECK_THAT(T, WithinRel(expected_T, 1.0e-14));

    // Consistency check: p == h^2 / mu for state at periapsis
    // At periapsis, r = rp, v = sqrt(mu * (2/rp - 1/a)) = sqrt(mu/a * (1+e)/(1-e))
    const double vp = std::sqrt(mu / a * (1.0 + e) / (1.0 - e));
    const double h = rp * vp;
    const double p_from_h = (h * h) / mu;
    CHECK_THAT(p, WithinRel(p_from_h, 1.0e-13));
}

// ===========================================================================
// Orbit Classification Tests
// ===========================================================================

TEST_CASE("Orbit classification accurately detects shapes and inclinations", "[elements][classification]") {
    ClassicalOrbitalElements circ_eq{7000000.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    auto c1 = classify_orbit(circ_eq);
    CHECK(c1.is_circular);
    CHECK(c1.is_equatorial);
    CHECK(c1.shape == OrbitShape::circular);
    CHECK(c1.inclination == OrbitInclination::equatorial_prograde);

    ClassicalOrbitalElements inc_ell{10000000.0, 0.2, pi * 0.25, 0.0, 0.0, 0.0};
    auto c2 = classify_orbit(inc_ell);
    CHECK_FALSE(c2.is_circular);
    CHECK_FALSE(c2.is_equatorial);
    CHECK(c2.shape == OrbitShape::elliptic);
    CHECK(c2.inclination == OrbitInclination::inclined);

    ClassicalOrbitalElements polar_circ{7000000.0, 0.0, pi * 0.5, 0.0, 0.0, 0.0};
    auto c3 = classify_orbit(polar_circ);
    CHECK(c3.is_circular);
    CHECK_FALSE(c3.is_equatorial);
    CHECK(c3.inclination == OrbitInclination::polar);
}

// ===========================================================================
// Independent Analytical State -> Elements Verification
// ===========================================================================

TEST_CASE("Independent analytical state to classical elements calculation", "[elements][state_to_elements]") {
    const double mu = 3.986004418e14;

    // Hand-crafted non-singular inclined elliptical state:
    // a = 10,000 km, e = 0.2, i = 45 deg, RAAN = 120 deg, argp = 60 deg, nu = 30 deg
    const ClassicalOrbitalElements known_elem{
        10000000.0,
        0.2,
        45.0 * pi / 180.0,
        120.0 * pi / 180.0,
        60.0 * pi / 180.0,
        30.0 * pi / 180.0,
    };

    // Calculate Cartesian state independently:
    const CartesianState state = classical_elements_to_state(known_elem, mu);

    // Now convert state back to elements:
    const ClassicalOrbitalElements derived_elem = state_to_classical_elements(state, mu);

    // Verify each element independently:
    CHECK_THAT(derived_elem.semi_major_axis_m, WithinRel(known_elem.semi_major_axis_m, 1.0e-12));
    CHECK_THAT(derived_elem.eccentricity, WithinRel(known_elem.eccentricity, 1.0e-12));
    CHECK_THAT(derived_elem.inclination_rad, WithinRel(known_elem.inclination_rad, 1.0e-12));
    CHECK_THAT(derived_elem.raan_rad, WithinRel(known_elem.raan_rad, 1.0e-12));
    CHECK_THAT(derived_elem.argument_of_periapsis_rad, WithinRel(known_elem.argument_of_periapsis_rad, 1.0e-12));
    CHECK_THAT(derived_elem.true_anomaly_rad, WithinRel(known_elem.true_anomaly_rad, 1.0e-12));
}

// ===========================================================================
// Comprehensive Round-Trip Verification (Cases A through F)
// ===========================================================================

TEST_CASE("Comprehensive round-trip verification across all orbital regimes", "[elements][round_trip]") {
    const double mu = 3.986004418e14;

    struct TestCase {
        const char* name;
        ClassicalOrbitalElements elements;
    };

    const std::vector<TestCase> test_cases = {
        {"Case A: Circular Equatorial", {7000000.0, 0.0, 0.0, 0.0, 0.0, 45.0 * pi / 180.0}},
        {"Case B: Circular Inclined (i = 45 deg)", {7000000.0, 0.0, 45.0 * pi / 180.0, 60.0 * pi / 180.0, 0.0, 90.0 * pi / 180.0}},
        {"Case C: Elliptical Equatorial (e = 0.2)", {10000000.0, 0.2, 0.0, 0.0, 30.0 * pi / 180.0, 60.0 * pi / 180.0}},
        {"Case D: Elliptical Inclined (a = 10,000 km, e = 0.2, i = 45 deg)", {10000000.0, 0.2, 45.0 * pi / 180.0, 120.0 * pi / 180.0, 60.0 * pi / 180.0, 30.0 * pi / 180.0}},
        {"Case E: High-Inclination Elliptical (i = 98 deg Sun-Sync)", {8000000.0, 0.15, 98.0 * pi / 180.0, 45.0 * pi / 180.0, 270.0 * pi / 180.0, 150.0 * pi / 180.0}},
        {"Case F: Near-Circular (e = 1e-5)", {7000000.0, 1.0e-5, 30.0 * pi / 180.0, 15.0 * pi / 180.0, 45.0 * pi / 180.0, 120.0 * pi / 180.0}},
    };

    for (const auto& tc : test_cases) {
        DYNAMIC_SECTION(tc.name) {
            // 1. Elements -> State
            const CartesianState state = classical_elements_to_state(tc.elements, mu);

            // 2. State -> Elements
            const ClassicalOrbitalElements converted_elem = state_to_classical_elements(state, mu);

            // 3. Elements -> State'
            const CartesianState reconstructed_state = classical_elements_to_state(converted_elem, mu);

            // Verify position and velocity round-trip error
            const double pos_error = (reconstructed_state.position - state.position).norm();
            const double vel_error = (reconstructed_state.velocity - state.velocity).norm();

            // Position error should be < 1.0e-8 m (sub-millimetre)
            CHECK(pos_error < 1.0e-8);
            // Velocity error should be < 1.0e-11 m/s
            CHECK(vel_error < 1.0e-11);
        }
    }
}

// ===========================================================================
// Multi-Orbit Orbital Element Invariance Verification
// ===========================================================================

TEST_CASE("Orbital element invariants under ideal two-body propagation", "[elements][invariants]") {
    const double mu = astradock::constants::earth_gravitational_parameter_m3_per_s2;

    // Inclined elliptical orbit (Case D)
    const ClassicalOrbitalElements initial_elem{
        10000000.0,
        0.2,
        45.0 * pi / 180.0,
        120.0 * pi / 180.0,
        60.0 * pi / 180.0,
        30.0 * pi / 180.0,
    };

    const CartesianState state0 = classical_elements_to_state(initial_elem, mu);
    const double period_s = astradock::orbit::orbital_period_s(initial_elem.semi_major_axis_m, mu);

    const double dt = 10.0;
    const double duration = 3.0 * period_s; // 3 full orbits

    const auto trajectory = propagate_fixed_step(
        0.0, duration, dt, state0,
        IntegrationMethod::classical_rk4,
        [mu](double t, const CartesianState& s) {
            return two_body_state_derivative(t, s, mu);
        });

    REQUIRE(trajectory.size() > 2000);

    double max_a_err = 0.0;
    double max_e_err = 0.0;
    double max_i_err = 0.0;
    double max_raan_err = 0.0;
    double max_argp_err = 0.0;

    for (const auto& sample : trajectory) {
        const ClassicalOrbitalElements elem = state_to_classical_elements(sample.state, mu);

        const double a_err = std::abs(elem.semi_major_axis_m - initial_elem.semi_major_axis_m);
        const double e_err = std::abs(elem.eccentricity - initial_elem.eccentricity);
        const double i_err = angular_distance_rad(elem.inclination_rad, initial_elem.inclination_rad);
        const double raan_err = angular_distance_rad(elem.raan_rad, initial_elem.raan_rad);
        const double argp_err = angular_distance_rad(elem.argument_of_periapsis_rad, initial_elem.argument_of_periapsis_rad);

        max_a_err = std::max(max_a_err, a_err);
        max_e_err = std::max(max_e_err, e_err);
        max_i_err = std::max(max_i_err, i_err);
        max_raan_err = std::max(max_raan_err, raan_err);
        max_argp_err = std::max(max_argp_err, argp_err);
    }

    // In ideal two-body motion with RK4 (dt = 10 s):
    // Semi-major axis drift should be < 0.01 m
    CHECK(max_a_err < 0.01);
    // Eccentricity drift should be < 1.0e-9
    CHECK(max_e_err < 1.0e-9);
    // Inclination drift should be < 1.0e-11 rad
    CHECK(max_i_err < 1.0e-11);
    // RAAN drift should be < 1.0e-11 rad
    CHECK(max_raan_err < 1.0e-11);
    // Argument of periapsis drift should be < 1.0e-8 rad
    CHECK(max_argp_err < 1.0e-8);
}

// ===========================================================================
// Defensive / Singularity and Invalid State Tests
// ===========================================================================

TEST_CASE("Defensive validation rejects invalid and unbound states", "[elements][defensive]") {
    const double mu = 3.986004418e14;

    // Zero position
    CHECK_THROWS_AS(
        state_to_classical_elements({{0.0, 0.0, 0.0}, {0.0, 7500.0, 0.0}}, mu),
        std::domain_error);

    // Zero angular momentum (radial velocity collinear with position)
    CHECK_THROWS_AS(
        state_to_classical_elements({{7000000.0, 0.0, 0.0}, {7500.0, 0.0, 0.0}}, mu),
        std::domain_error);

    // Unbound hyperbolic trajectory (v > escape velocity: v_esc = sqrt(2*mu/r) ≈ 10672 m/s at 7000 km)
    CHECK_THROWS_AS(
        state_to_classical_elements({{7000000.0, 0.0, 0.0}, {0.0, 12000.0, 0.0}}, mu),
        std::domain_error);

    // Invalid parameters to elements_to_state
    CHECK_THROWS_AS(
        classical_elements_to_state({-10000.0, 0.1, 0.0, 0.0, 0.0, 0.0}, mu),
        std::domain_error);

    CHECK_THROWS_AS(
        classical_elements_to_state({10000000.0, 1.2, 0.0, 0.0, 0.0, 0.0}, mu),
        std::domain_error);

    CHECK_THROWS_AS(
        classical_elements_to_state({10000000.0, 0.1, 0.0, 0.0, 0.0, 0.0}, -mu),
        std::domain_error);
}
