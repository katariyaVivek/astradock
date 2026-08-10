#include "math/vector3.hpp"
#include "numerics/integrators.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace {

using astradock::math::Vector3;
using astradock::numerics::euler_step;
using astradock::numerics::rk4_step;
using Catch::Approx;

struct EulerStepper {
    template <typename State, typename Derivative>
    [[nodiscard]] State operator()(
        double t,
        const State& state,
        double dt,
        Derivative& derivative) const {
        return euler_step(t, state, dt, derivative);
    }
};

struct Rk4Stepper {
    template <typename State, typename Derivative>
    [[nodiscard]] State operator()(
        double t,
        const State& state,
        double dt,
        Derivative& derivative) const {
        return rk4_step(t, state, dt, derivative);
    }
};

template <typename State, typename Stepper, typename Derivative>
[[nodiscard]] State integrate_fixed_steps(
    State state,
    double initial_time,
    double dt,
    std::size_t step_count,
    Stepper stepper,
    Derivative derivative) {
    double t = initial_time;
    for (std::size_t step = 0; step < step_count; ++step) {
        state = stepper(t, state, dt, derivative);
        t += dt;
    }
    return state;
}

template <typename Stepper>
[[nodiscard]] double exponential_decay_error(double dt, Stepper stepper) {
    const auto decay_derivative = [](double, double state) { return -state; };
    const auto step_count = static_cast<std::size_t>(std::llround(1.0 / dt));
    const double numerical =
        integrate_fixed_steps(1.0, 0.0, dt, step_count, stepper, decay_derivative);
    return std::abs(numerical - std::exp(-1.0));
}

TEST_CASE("Euler one-step result matches a calculation by hand", "[numerics][euler]") {
    const auto derivative = [](double t, double state) { return t + state; };

    // y_next = 2 + 0.25 * (0 + 2) = 2.5
    const double result = euler_step(0.0, 2.0, 0.25, derivative);

    REQUIRE(result == Approx(2.5));
}

TEST_CASE("Euler approaches exponential decay as the step decreases",
          "[numerics][euler][convergence]") {
    const double error_0_1 = exponential_decay_error(0.1, EulerStepper{});
    const double error_0_05 = exponential_decay_error(0.05, EulerStepper{});
    const double error_0_025 = exponential_decay_error(0.025, EulerStepper{});

    REQUIRE(error_0_05 < error_0_1);
    REQUIRE(error_0_025 < error_0_05);
    REQUIRE(error_0_025 < 0.005);
}

TEST_CASE("Euler demonstrates first-order convergence", "[numerics][euler][convergence]") {
    const double coarse_error = exponential_decay_error(0.1, EulerStepper{});
    const double medium_error = exponential_decay_error(0.05, EulerStepper{});
    const double fine_error = exponential_decay_error(0.025, EulerStepper{});
    const double coarse_to_medium_ratio = coarse_error / medium_error;
    const double medium_to_fine_ratio = medium_error / fine_error;

    // First-order convergence predicts a ratio approaching two when dt halves.
    REQUIRE(coarse_to_medium_ratio > 1.8);
    REQUIRE(coarse_to_medium_ratio < 2.3);
    REQUIRE(medium_to_fine_ratio > 1.8);
    REQUIRE(medium_to_fine_ratio < 2.3);
}

TEST_CASE("Euler reproduces a constant derivative", "[numerics][euler]") {
    const auto constant_derivative = [](double, double) { return 2.5; };
    const double result =
        integrate_fixed_steps(-1.0, 0.0, 0.125, 8, EulerStepper{}, constant_derivative);

    REQUIRE(result == Approx(1.5).margin(1.0e-14));
}

TEST_CASE("RK4 one-step result matches the fourth-degree exponential polynomial",
          "[numerics][rk4]") {
    constexpr double dt = 0.1;
    const auto derivative = [](double, double state) { return state; };
    const double expected =
        1.0 + dt + dt * dt / 2.0 + dt * dt * dt / 6.0 + dt * dt * dt * dt / 24.0;

    const double result = rk4_step(0.0, 1.0, dt, derivative);

    REQUIRE(result == Approx(expected).margin(1.0e-14));
}

TEST_CASE("RK4 tracks exponential decay accurately", "[numerics][rk4]") {
    const auto decay_derivative = [](double, double state) { return -state; };
    const double result =
        integrate_fixed_steps(1.0, 0.0, 0.1, 10, Rk4Stepper{}, decay_derivative);

    REQUIRE(result == Approx(std::exp(-1.0)).margin(4.0e-7));
}

TEST_CASE("RK4 is materially more accurate than Euler at equal step size",
          "[numerics][euler][rk4]") {
    const double euler_error = exponential_decay_error(0.1, EulerStepper{});
    const double rk4_error = exponential_decay_error(0.1, Rk4Stepper{});

    REQUIRE(rk4_error < euler_error / 1'000.0);
}

TEST_CASE("RK4 demonstrates fourth-order convergence", "[numerics][rk4][convergence]") {
    const double coarse_error = exponential_decay_error(0.2, Rk4Stepper{});
    const double medium_error = exponential_decay_error(0.1, Rk4Stepper{});
    const double fine_error = exponential_decay_error(0.05, Rk4Stepper{});
    const double coarse_to_medium_ratio = coarse_error / medium_error;
    const double medium_to_fine_ratio = medium_error / fine_error;

    // Fourth-order convergence predicts a ratio approaching 16 when dt halves.
    REQUIRE(coarse_to_medium_ratio > 14.0);
    REQUIRE(coarse_to_medium_ratio < 20.0);
    REQUIRE(medium_to_fine_ratio > 14.0);
    REQUIRE(medium_to_fine_ratio < 20.0);
}

TEST_CASE("Vector3 state integrates exponential decay componentwise",
          "[numerics][euler][rk4][vector3]") {
    const Vector3 initial{1.0, 2.0, -3.0};
    const auto decay_derivative = [](double, const Vector3& state) { return -1.0 * state; };
    const Vector3 analytical = initial * std::exp(-1.0);

    const Vector3 euler_result =
        integrate_fixed_steps(initial, 0.0, 0.1, 10, EulerStepper{}, decay_derivative);
    const Vector3 rk4_result =
        integrate_fixed_steps(initial, 0.0, 0.1, 10, Rk4Stepper{}, decay_derivative);

    REQUIRE((euler_result - analytical).norm() < 0.08);
    REQUIRE((rk4_result - analytical).norm() < 1.5e-6);
    REQUIRE((rk4_result - analytical).norm() < (euler_result - analytical).norm());
}

TEST_CASE("A zero time step is an identity and does not evaluate the derivative",
          "[numerics][euler][rk4][policy]") {
    bool derivative_was_called = false;
    const auto derivative = [&derivative_was_called](double, double) {
        derivative_was_called = true;
        return 1.0;
    };

    REQUIRE(euler_step(2.0, 3.0, 0.0, derivative) == Approx(3.0));
    REQUIRE(rk4_step(2.0, 3.0, 0.0, derivative) == Approx(3.0));
    REQUIRE_FALSE(derivative_was_called);
}

TEST_CASE("Negative time steps support backward integration",
          "[numerics][euler][rk4][policy]") {
    const auto constant_derivative = [](double, double) { return 2.0; };

    REQUIRE(euler_step(1.0, 3.0, -0.25, constant_derivative) == Approx(2.5));
    REQUIRE(rk4_step(1.0, 3.0, -0.25, constant_derivative) == Approx(2.5));
}

TEST_CASE("Integrators reject non-finite supplied numerical values",
          "[numerics][euler][rk4][policy]") {
    const double infinity = std::numeric_limits<double>::infinity();
    const double not_a_number = std::numeric_limits<double>::quiet_NaN();
    const auto finite_derivative = [](double, double state) { return state; };

    REQUIRE_THROWS_AS(euler_step(not_a_number, 1.0, 0.1, finite_derivative), std::domain_error);
    REQUIRE_THROWS_AS(rk4_step(0.0, 1.0, infinity, finite_derivative), std::domain_error);
    REQUIRE_THROWS_AS(euler_step(0.0, infinity, 0.1, finite_derivative), std::domain_error);
    REQUIRE_THROWS_AS(
        rk4_step(0.0, Vector3{1.0, not_a_number, 3.0}, 0.1,
                 [](double, const Vector3& state) { return state; }),
        std::domain_error);
}

TEST_CASE("Integrators reject non-finite derivatives and calculations",
          "[numerics][euler][rk4][policy]") {
    const double infinity = std::numeric_limits<double>::infinity();
    const double largest = std::numeric_limits<double>::max();
    const auto non_finite_derivative = [infinity](double, double) { return infinity; };
    const auto largest_derivative = [largest](double, double) { return largest; };

    REQUIRE_THROWS_AS(euler_step(0.0, 1.0, 0.1, non_finite_derivative), std::domain_error);
    REQUIRE_THROWS_AS(rk4_step(0.0, 1.0, 0.1, non_finite_derivative), std::domain_error);
    REQUIRE_THROWS_AS(euler_step(0.0, largest, 1.0, largest_derivative), std::overflow_error);
    REQUIRE_THROWS_AS(rk4_step(0.0, largest, 1.0, largest_derivative), std::overflow_error);
    REQUIRE_THROWS_AS(euler_step(largest, 1.0, largest, largest_derivative), std::overflow_error);
}

}  // namespace
