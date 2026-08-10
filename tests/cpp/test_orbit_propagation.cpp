#include "dynamics/two_body.hpp"
#include "math/constants.hpp"
#include "math/vector3.hpp"
#include "numerics/fixed_step_propagation.hpp"
#include "orbit/cartesian_state.hpp"
#include "orbit/two_body_orbit.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

using astradock::math::Vector3;
using astradock::numerics::IntegrationMethod;
using astradock::numerics::StateSample;
using astradock::orbit::CartesianState;
using Catch::Approx;

constexpr double altitude_m = 500'000.0;
constexpr double nominal_dt_s = 10.0;
constexpr double mu_m3_per_s2 =
    astradock::constants::earth_gravitational_parameter_m3_per_s2;
constexpr double earth_radius_m = astradock::constants::earth_reference_radius_m;
constexpr double orbital_radius_m = earth_radius_m + altitude_m;

struct CircularOrbitScenario {
    CartesianState initial_state;
    double period_s;
};

struct OrbitMetrics {
    double final_position_error_m;
    double final_velocity_error_m_per_s;
    double maximum_radius_deviation_m;
    double maximum_relative_energy_drift;
    double maximum_relative_angular_momentum_drift;
};

[[nodiscard]] CircularOrbitScenario circular_orbit_scenario() {
    const double circular_speed_m_per_s =
        astradock::orbit::circular_orbit_speed_m_per_s(mu_m3_per_s2, orbital_radius_m);
    return {
        {
            {orbital_radius_m, 0.0, 0.0},
            {0.0, circular_speed_m_per_s, 0.0},
        },
        astradock::orbit::circular_orbit_period_s(mu_m3_per_s2, orbital_radius_m),
    };
}

[[nodiscard]] std::vector<StateSample<CartesianState>> propagate_orbits(
    IntegrationMethod method,
    double orbit_count) {
    const CircularOrbitScenario scenario = circular_orbit_scenario();
    const auto derivative = [](double time_s, const CartesianState& state) {
        return astradock::orbit::two_body_state_derivative(time_s, state, mu_m3_per_s2);
    };
    return astradock::numerics::propagate_fixed_step(
        0.0,
        orbit_count * scenario.period_s,
        nominal_dt_s,
        scenario.initial_state,
        method,
        derivative);
}

[[nodiscard]] OrbitMetrics calculate_metrics(
    const std::vector<StateSample<CartesianState>>& samples,
    const CartesianState& initial_state) {
    const double initial_radius_m = initial_state.position.norm();
    const double initial_energy_m2_per_s2 =
        astradock::orbit::specific_orbital_energy_m2_per_s2(initial_state, mu_m3_per_s2);
    const double initial_angular_momentum_m2_per_s =
        astradock::orbit::specific_angular_momentum_m2_per_s(initial_state).norm();

    double maximum_radius_deviation_m = 0.0;
    double maximum_relative_energy_drift = 0.0;
    double maximum_relative_angular_momentum_drift = 0.0;
    for (const auto& sample : samples) {
        const double energy_m2_per_s2 =
            astradock::orbit::specific_orbital_energy_m2_per_s2(
                sample.state,
                mu_m3_per_s2);
        const double angular_momentum_m2_per_s =
            astradock::orbit::specific_angular_momentum_m2_per_s(sample.state).norm();

        maximum_radius_deviation_m = std::max(
            maximum_radius_deviation_m,
            std::abs(sample.state.position.norm() - initial_radius_m));
        maximum_relative_energy_drift = std::max(
            maximum_relative_energy_drift,
            std::abs(
                (energy_m2_per_s2 - initial_energy_m2_per_s2)
                / std::abs(initial_energy_m2_per_s2)));
        maximum_relative_angular_momentum_drift = std::max(
            maximum_relative_angular_momentum_drift,
            std::abs(
                (angular_momentum_m2_per_s - initial_angular_momentum_m2_per_s)
                / initial_angular_momentum_m2_per_s));
    }

    return {
        (samples.back().state.position - initial_state.position).norm(),
        (samples.back().state.velocity - initial_state.velocity).norm(),
        maximum_radius_deviation_m,
        maximum_relative_energy_drift,
        maximum_relative_angular_momentum_drift,
    };
}

TEST_CASE("Cartesian state arithmetic supports integration stages", "[orbit][state]") {
    const CartesianState lhs{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
    const CartesianState rhs{{-1.0, 0.5, 2.0}, {3.0, -2.0, 1.0}};

    const CartesianState sum = lhs + rhs;
    const CartesianState right_scaled = lhs * 2.0;
    const CartesianState left_scaled = 0.5 * lhs;

    REQUIRE(astradock::math::approximately_equal(sum.position, Vector3{0.0, 2.5, 5.0}));
    REQUIRE(astradock::math::approximately_equal(sum.velocity, Vector3{7.0, 3.0, 7.0}));
    REQUIRE(astradock::math::approximately_equal(right_scaled.position, Vector3{2.0, 4.0, 6.0}));
    REQUIRE(astradock::math::approximately_equal(right_scaled.velocity, Vector3{8.0, 10.0, 12.0}));
    REQUIRE(astradock::math::approximately_equal(left_scaled.position, Vector3{0.5, 1.0, 1.5}));
    REQUIRE(astradock::math::approximately_equal(left_scaled.velocity, Vector3{2.0, 2.5, 3.0}));
}

TEST_CASE("Two-body state derivative sets position rate to velocity",
          "[orbit][derivative]") {
    const CartesianState state{{7.0e6, -1.0e6, 2.0e6}, {1.0, 2.0, 3.0}};

    const CartesianState derivative =
        astradock::orbit::two_body_state_derivative(123.0, state, mu_m3_per_s2);

    REQUIRE(astradock::math::approximately_equal(derivative.position, state.velocity));
}

TEST_CASE("Two-body state derivative reuses canonical gravitational acceleration",
          "[orbit][derivative]") {
    const CartesianState state{{7.0e6, -1.0e6, 2.0e6}, {1.0, 2.0, 3.0}};
    const Vector3 expected_acceleration =
        astradock::dynamics::two_body_acceleration(state.position, mu_m3_per_s2);

    const CartesianState derivative =
        astradock::orbit::two_body_state_derivative(123.0, state, mu_m3_per_s2);

    REQUIRE(astradock::math::approximately_equal(derivative.velocity, expected_acceleration));
}

TEST_CASE("Two-body state derivative rejects invalid orbital states",
          "[orbit][derivative]") {
    const double infinity = std::numeric_limits<double>::infinity();

    REQUIRE_THROWS_AS(
        astradock::orbit::two_body_state_derivative(
            0.0,
            CartesianState{{0.0, 0.0, 0.0}, {1.0, 2.0, 3.0}},
            mu_m3_per_s2),
        std::domain_error);
    REQUIRE_THROWS_AS(
        astradock::orbit::two_body_state_derivative(
            0.0,
            CartesianState{{orbital_radius_m, 0.0, 0.0}, {0.0, infinity, 0.0}},
            mu_m3_per_s2),
        std::domain_error);
}

TEST_CASE("Circular orbit helpers satisfy analytical mechanics identities",
          "[orbit][diagnostics]") {
    const CircularOrbitScenario scenario = circular_orbit_scenario();
    const double speed_m_per_s = scenario.initial_state.velocity.norm();
    const double energy_m2_per_s2 =
        astradock::orbit::specific_orbital_energy_m2_per_s2(
            scenario.initial_state,
            mu_m3_per_s2);
    const Vector3 angular_momentum_m2_per_s =
        astradock::orbit::specific_angular_momentum_m2_per_s(scenario.initial_state);

    REQUIRE(speed_m_per_s * speed_m_per_s
            == Approx(mu_m3_per_s2 / orbital_radius_m).epsilon(1.0e-14));
    REQUIRE(scenario.period_s
            == Approx(2.0 * astradock::constants::pi * orbital_radius_m / speed_m_per_s)
                   .epsilon(1.0e-14));
    REQUIRE(energy_m2_per_s2
            == Approx(-mu_m3_per_s2 / (2.0 * orbital_radius_m)).epsilon(1.0e-14));
    REQUIRE(angular_momentum_m2_per_s.x() == Approx(0.0));
    REQUIRE(angular_momentum_m2_per_s.y() == Approx(0.0));
    REQUIRE(angular_momentum_m2_per_s.z()
            == Approx(orbital_radius_m * speed_m_per_s).epsilon(1.0e-14));
}

TEST_CASE("Fixed-step propagation includes an exactly divisible endpoint",
          "[numerics][propagation][endpoint]") {
    const auto derivative = [](double, double) { return 1.0; };
    const auto samples = astradock::numerics::propagate_fixed_step(
        0.0,
        1.0,
        0.25,
        0.0,
        IntegrationMethod::forward_euler,
        derivative);

    REQUIRE(samples.size() == 5);
    REQUIRE(samples.front().time_s == 0.0);
    REQUIRE(samples.back().time_s == 1.0);
    REQUIRE(samples.back().state == Approx(1.0));
}

TEST_CASE("Fixed-step propagation shortens its final step without overshoot",
          "[numerics][propagation][endpoint]") {
    const auto derivative = [](double, double) { return 1.0; };
    const auto samples = astradock::numerics::propagate_fixed_step(
        0.0,
        1.0,
        0.3,
        0.0,
        IntegrationMethod::classical_rk4,
        derivative);

    REQUIRE(samples.size() == 5);
    REQUIRE(samples[samples.size() - 2].time_s == Approx(0.9));
    REQUIRE(samples.back().time_s == 1.0);
    REQUIRE(samples.back().state == Approx(1.0));
}

TEST_CASE("Zero-duration propagation returns only the initial sample",
          "[numerics][propagation][endpoint]") {
    bool derivative_was_called = false;
    const auto derivative = [&derivative_was_called](double, double) {
        derivative_was_called = true;
        return 1.0;
    };
    const auto samples = astradock::numerics::propagate_fixed_step(
        5.0,
        5.0,
        0.25,
        7.0,
        IntegrationMethod::classical_rk4,
        derivative);

    REQUIRE(samples.size() == 1);
    REQUIRE(samples.front().time_s == 5.0);
    REQUIRE(samples.front().state == Approx(7.0));
    REQUIRE_FALSE(derivative_was_called);
}

TEST_CASE("Fixed-step propagation supports signed backward steps",
          "[numerics][propagation][endpoint]") {
    const auto derivative = [](double, double) { return 2.0; };
    const auto samples = astradock::numerics::propagate_fixed_step(
        1.0,
        0.0,
        -0.3,
        0.0,
        IntegrationMethod::classical_rk4,
        derivative);

    REQUIRE(samples.size() == 5);
    REQUIRE(samples.back().time_s == 0.0);
    REQUIRE(samples.back().state == Approx(-2.0));
}

TEST_CASE("Fixed-step propagation rejects invalid time policies",
          "[numerics][propagation][policy]") {
    const double infinity = std::numeric_limits<double>::infinity();
    const double not_a_number = std::numeric_limits<double>::quiet_NaN();
    const auto derivative = [](double, double state) { return state; };

    REQUIRE_THROWS_AS(
        astradock::numerics::propagate_fixed_step(
            not_a_number,
            1.0,
            0.1,
            1.0,
            IntegrationMethod::forward_euler,
            derivative),
        std::domain_error);
    REQUIRE_THROWS_AS(
        astradock::numerics::propagate_fixed_step(
            0.0,
            infinity,
            0.1,
            1.0,
            IntegrationMethod::forward_euler,
            derivative),
        std::domain_error);
    REQUIRE_THROWS_AS(
        astradock::numerics::propagate_fixed_step(
            0.0,
            1.0,
            0.0,
            1.0,
            IntegrationMethod::forward_euler,
            derivative),
        std::domain_error);
    REQUIRE_THROWS_AS(
        astradock::numerics::propagate_fixed_step(
            0.0,
            1.0,
            -0.1,
            1.0,
            IntegrationMethod::forward_euler,
            derivative),
        std::domain_error);
    REQUIRE_THROWS_AS(
        astradock::numerics::propagate_fixed_step(
            1.0,
            0.0,
            0.1,
            1.0,
            IntegrationMethod::classical_rk4,
            derivative),
        std::domain_error);
}

TEST_CASE("Short RK4 orbit propagation changes state while remaining physical",
          "[orbit][propagation][sanity]") {
    const CircularOrbitScenario scenario = circular_orbit_scenario();
    const auto derivative = [](double time_s, const CartesianState& state) {
        return astradock::orbit::two_body_state_derivative(time_s, state, mu_m3_per_s2);
    };
    const auto samples = astradock::numerics::propagate_fixed_step(
        0.0,
        60.0,
        nominal_dt_s,
        scenario.initial_state,
        IntegrationMethod::classical_rk4,
        derivative);
    const CartesianState& final_state = samples.back().state;

    REQUIRE((final_state.position - scenario.initial_state.position).norm() > 1.0);
    REQUIRE((final_state.velocity - scenario.initial_state.velocity).norm() > 1.0e-6);
    REQUIRE(astradock::orbit::is_finite(final_state));
    REQUIRE(final_state.position.norm() > earth_radius_m);
    REQUIRE(final_state.position.norm() < orbital_radius_m + 10'000.0);
}

TEST_CASE("RK4 closes a 500 kilometre circular orbit after one analytical period",
          "[orbit][propagation][rk4]") {
    const CircularOrbitScenario scenario = circular_orbit_scenario();
    const OrbitMetrics metrics = calculate_metrics(
        propagate_orbits(IntegrationMethod::classical_rk4, 1.0),
        scenario.initial_state);

    REQUIRE(metrics.final_position_error_m < 1.0);
    REQUIRE(metrics.final_velocity_error_m_per_s < 0.002);
}

TEST_CASE("RK4 keeps circular-orbit radius stable", "[orbit][propagation][rk4]") {
    const CircularOrbitScenario scenario = circular_orbit_scenario();
    const OrbitMetrics metrics = calculate_metrics(
        propagate_orbits(IntegrationMethod::classical_rk4, 1.0),
        scenario.initial_state);

    REQUIRE(metrics.maximum_radius_deviation_m < 5.0);
}

TEST_CASE("RK4 approximately conserves specific orbital energy",
          "[orbit][propagation][rk4][invariant]") {
    const CircularOrbitScenario scenario = circular_orbit_scenario();
    const OrbitMetrics metrics = calculate_metrics(
        propagate_orbits(IntegrationMethod::classical_rk4, 1.0),
        scenario.initial_state);

    REQUIRE(metrics.maximum_relative_energy_drift < 1.0e-8);
}

TEST_CASE("RK4 approximately conserves specific angular momentum",
          "[orbit][propagation][rk4][invariant]") {
    const CircularOrbitScenario scenario = circular_orbit_scenario();
    const auto samples = propagate_orbits(IntegrationMethod::classical_rk4, 1.0);
    const OrbitMetrics metrics = calculate_metrics(samples, scenario.initial_state);
    const Vector3 initial_angular_momentum_m2_per_s =
        astradock::orbit::specific_angular_momentum_m2_per_s(scenario.initial_state);

    double maximum_relative_vector_drift = 0.0;
    for (const auto& sample : samples) {
        const Vector3 angular_momentum_m2_per_s =
            astradock::orbit::specific_angular_momentum_m2_per_s(sample.state);
        maximum_relative_vector_drift = std::max(
            maximum_relative_vector_drift,
            (angular_momentum_m2_per_s - initial_angular_momentum_m2_per_s).norm()
                / initial_angular_momentum_m2_per_s.norm());
    }

    REQUIRE(metrics.maximum_relative_angular_momentum_drift < 1.0e-8);
    REQUIRE(maximum_relative_vector_drift < 1.0e-8);
}

TEST_CASE("RK4 materially outperforms Euler for one circular orbit",
          "[orbit][propagation][comparison]") {
    const CircularOrbitScenario scenario = circular_orbit_scenario();
    const OrbitMetrics euler = calculate_metrics(
        propagate_orbits(IntegrationMethod::forward_euler, 1.0),
        scenario.initial_state);
    const OrbitMetrics rk4 = calculate_metrics(
        propagate_orbits(IntegrationMethod::classical_rk4, 1.0),
        scenario.initial_state);

    REQUIRE(euler.final_position_error_m > rk4.final_position_error_m * 1'000.0);
    REQUIRE(euler.final_velocity_error_m_per_s > rk4.final_velocity_error_m_per_s * 1'000.0);
    REQUIRE(euler.maximum_radius_deviation_m > rk4.maximum_radius_deviation_m * 1'000.0);
    REQUIRE(euler.maximum_relative_energy_drift > rk4.maximum_relative_energy_drift * 1'000.0);
    REQUIRE(euler.maximum_relative_angular_momentum_drift
            > rk4.maximum_relative_angular_momentum_drift * 1'000.0);
}

TEST_CASE("Five-orbit comparison exposes accumulated Euler drift",
          "[orbit][propagation][comparison][multi-orbit]") {
    const CircularOrbitScenario scenario = circular_orbit_scenario();
    const OrbitMetrics euler = calculate_metrics(
        propagate_orbits(IntegrationMethod::forward_euler, 5.0),
        scenario.initial_state);
    const OrbitMetrics rk4 = calculate_metrics(
        propagate_orbits(IntegrationMethod::classical_rk4, 5.0),
        scenario.initial_state);

    REQUIRE(euler.maximum_radius_deviation_m > 100'000.0);
    REQUIRE(euler.maximum_relative_energy_drift > 1.0e-3);
    REQUIRE(rk4.maximum_radius_deviation_m < 10.0);
    REQUIRE(rk4.maximum_relative_energy_drift < 1.0e-7);
    REQUIRE(rk4.maximum_relative_angular_momentum_drift < 1.0e-7);
}

}  // namespace
