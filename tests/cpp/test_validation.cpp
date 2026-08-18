#include "math/constants.hpp"
#include "numerics/fixed_step_propagation.hpp"
#include "orbit/cartesian_state.hpp"
#include "orbit/orbital_diagnostics.hpp"
#include "orbit/two_body_orbit.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

namespace {

using astradock::math::Vector3;
using astradock::numerics::IntegrationMethod;
using astradock::numerics::StateSample;
using astradock::orbit::CartesianState;
using astradock::orbit::CircularOrbitReference;
using Catch::Approx;

constexpr double mu_m3_per_s2 =
    astradock::constants::earth_gravitational_parameter_m3_per_s2;
constexpr double earth_radius_m =
    astradock::constants::earth_reference_radius_m;
constexpr double altitude_m = 500'000.0;
constexpr double orbital_radius_m = earth_radius_m + altitude_m;

CircularOrbitReference nominal_reference() {
    return astradock::orbit::compute_circular_orbit_reference(
        mu_m3_per_s2, orbital_radius_m);
}

CartesianState nominal_initial_state() {
    const CircularOrbitReference ref = nominal_reference();
    return {
        {ref.orbital_radius_m, 0.0, 0.0},
        {0.0, ref.speed_m_per_s, 0.0},
    };
}

std::vector<StateSample<CartesianState>> propagate(
    IntegrationMethod method, double dt_s, double duration_s) {
    const CartesianState initial_state = nominal_initial_state();
    const auto derivative = [](double t, const CartesianState& state) {
        return astradock::orbit::two_body_state_derivative(t, state, mu_m3_per_s2);
    };
    return astradock::numerics::propagate_fixed_step(
        0.0, duration_s, dt_s, initial_state, method, derivative);
}

}  // namespace

// ============================================================================
// Analytical Reference Helpers
// ============================================================================

TEST_CASE("Circular orbit reference derives exact analytical quantities",
          "[m05][reference]") {
    const CircularOrbitReference ref = nominal_reference();

    // v_c = sqrt(mu / r)
    const double expected_speed = std::sqrt(mu_m3_per_s2 / orbital_radius_m);
    REQUIRE(ref.speed_m_per_s == Approx(expected_speed).epsilon(1.0e-14));

    // T = 2*pi*sqrt(r^3 / mu)
    const double expected_period =
        2.0 * astradock::constants::pi * std::sqrt(
            orbital_radius_m * orbital_radius_m * orbital_radius_m / mu_m3_per_s2);
    REQUIRE(ref.period_s == Approx(expected_period).epsilon(1.0e-14));

    // epsilon = -mu / (2*r)
    const double expected_energy = -mu_m3_per_s2 / (2.0 * orbital_radius_m);
    REQUIRE(ref.specific_energy_m2_per_s2 == Approx(expected_energy).epsilon(1.0e-14));

    // h = r * v_c
    const double expected_h = orbital_radius_m * expected_speed;
    REQUIRE(ref.specific_angular_momentum_m2_per_s == Approx(expected_h).epsilon(1.0e-14));

    // g = mu / r^2
    const double expected_g = mu_m3_per_s2 / (orbital_radius_m * orbital_radius_m);
    REQUIRE(ref.gravitational_acceleration_m_per_s2 == Approx(expected_g).epsilon(1.0e-14));

    // mean motion n = sqrt(mu / r^3)
    const double expected_n = std::sqrt(mu_m3_per_s2 / (orbital_radius_m * orbital_radius_m * orbital_radius_m));
    REQUIRE(ref.mean_motion_rad_per_s == Approx(expected_n).epsilon(1.0e-14));
}

TEST_CASE("Circular orbit reference rejects invalid parameters",
          "[m05][reference]") {
    REQUIRE_THROWS_AS(
        astradock::orbit::compute_circular_orbit_reference(-1.0, orbital_radius_m),
        std::domain_error);
    REQUIRE_THROWS_AS(
        astradock::orbit::compute_circular_orbit_reference(mu_m3_per_s2, 0.0),
        std::domain_error);
    REQUIRE_THROWS_AS(
        astradock::orbit::compute_circular_orbit_reference(mu_m3_per_s2, -1000.0),
        std::domain_error);
}

// ============================================================================
// Phase Diagnostics
// ============================================================================

TEST_CASE("Phase angle at initial state is zero", "[m05][diagnostics][phase]") {
    const CartesianState state = nominal_initial_state();
    const double phase = astradock::orbit::diagnostics::phase_angle_rad(state);
    REQUIRE(phase == Approx(0.0).margin(1.0e-15));
}

TEST_CASE("Analytical phase grows linearly with time", "[m05][diagnostics][phase]") {
    const CircularOrbitReference ref = nominal_reference();
    const double mean_motion = ref.mean_motion_rad_per_s;

    const double phase_at_10s = astradock::orbit::diagnostics::analytical_phase_rad(
        10.0, mu_m3_per_s2, orbital_radius_m);
    const double phase_at_20s = astradock::orbit::diagnostics::analytical_phase_rad(
        20.0, mu_m3_per_s2, orbital_radius_m);

    REQUIRE(phase_at_10s == Approx(mean_motion * 10.0).epsilon(1.0e-14));
    REQUIRE(phase_at_20s == Approx(2.0 * phase_at_10s).epsilon(1.0e-14));
}

TEST_CASE("Unwrapped phase handles empty and multi-revolution trajectories",
          "[m05][diagnostics][phase]") {
    // Empty trajectory returns empty
    const std::vector<StateSample<CartesianState>> empty_samples;
    REQUIRE(astradock::orbit::diagnostics::unwrapped_phase_rad(empty_samples).empty());

    // Single sample returns initial phase
    const CartesianState state{{orbital_radius_m, 0.0, 0.0}, {0.0, 7600.0, 0.0}};
    const std::vector<StateSample<CartesianState>> single_sample{{0.0, state}};
    const auto single_res = astradock::orbit::diagnostics::unwrapped_phase_rad(single_sample);
    REQUIRE(single_res.size() == 1);
    REQUIRE(single_res.front() == Approx(0.0).margin(1.0e-15));

    // Two full orbits should unwrap smoothly to ~4*pi without discontinuous jumps
    const CircularOrbitReference ref = nominal_reference();
    const auto two_orbits = propagate(IntegrationMethod::classical_rk4, 10.0, 2.0 * ref.period_s);
    const auto unwrapped = astradock::orbit::diagnostics::unwrapped_phase_rad(two_orbits);

    REQUIRE(unwrapped.size() == two_orbits.size());
    REQUIRE(unwrapped.front() == Approx(0.0).margin(1.0e-12));
    REQUIRE(unwrapped.back() == Approx(4.0 * astradock::constants::pi).margin(1.0e-4));

    // Monotonically increasing for prograde circular orbit
    for (std::size_t i = 1; i < unwrapped.size(); ++i) {
        REQUIRE(unwrapped[i] > unwrapped[i - 1]);
    }
}

TEST_CASE("RK4 phase error is small for fine timestep", "[m05][diagnostics][phase]") {
    const CircularOrbitReference ref = nominal_reference();
    const auto samples = propagate(IntegrationMethod::classical_rk4, 1.0, ref.period_s);

    const double phase_err = astradock::orbit::diagnostics::max_phase_error_rad(
        samples, mu_m3_per_s2, orbital_radius_m);

    // Fine timestep RK4 should have phase error < 1e-4 rad
    REQUIRE(phase_err < 1.0e-4);
}

TEST_CASE("Euler phase error is large even at moderate timestep", "[m05][diagnostics][phase]") {
    const CircularOrbitReference ref = nominal_reference();
    const auto samples = propagate(IntegrationMethod::forward_euler, 10.0, ref.period_s);

    const double phase_err = astradock::orbit::diagnostics::max_phase_error_rad(
        samples, mu_m3_per_s2, orbital_radius_m);

    // Euler accumulates significant phase error (> 0.5 rad)
    REQUIRE(phase_err > 0.5);
}

// ============================================================================
// Period Estimation
// ============================================================================

TEST_CASE("Period is estimated from one RK4 orbit", "[m05][diagnostics][period]") {
    const CircularOrbitReference ref = nominal_reference();
    const auto samples = propagate(IntegrationMethod::classical_rk4, 10.0, ref.period_s);

    const double estimated = astradock::orbit::diagnostics::estimate_period_from_trajectory(
        samples, ref.period_s);

    REQUIRE(estimated > 0.0);
    // RK4 at 10s should estimate the period to within 1 ms
    REQUIRE(estimated == Approx(ref.period_s).margin(0.01));
}

TEST_CASE("Period estimation error decreases with finer timestep", "[m05][diagnostics][period]") {
    const CircularOrbitReference ref = nominal_reference();

    const auto coarse = propagate(IntegrationMethod::classical_rk4, 20.0, ref.period_s);
    const auto fine = propagate(IntegrationMethod::classical_rk4, 2.5, ref.period_s);

    const double coarse_est = astradock::orbit::diagnostics::estimate_period_from_trajectory(
        coarse, ref.period_s);
    const double fine_est = astradock::orbit::diagnostics::estimate_period_from_trajectory(
        fine, ref.period_s);

    const double coarse_err = std::abs(coarse_est - ref.period_s);
    const double fine_err = std::abs(fine_est - ref.period_s);

    REQUIRE(fine_err <= coarse_err);
    REQUIRE(fine_err < 1.0e-4);
}

TEST_CASE("Period estimation rejects invalid inputs and too-short trajectories",
          "[m05][diagnostics][period]") {
    const CircularOrbitReference ref = nominal_reference();

    // Short trajectory
    const auto samples = propagate(IntegrationMethod::classical_rk4, 10.0, 100.0);
    REQUIRE_THROWS_AS(
        astradock::orbit::diagnostics::estimate_period_from_trajectory(
            samples, ref.period_s),
        std::domain_error);

    // Invalid analytical hint
    REQUIRE_THROWS_AS(
        astradock::orbit::diagnostics::estimate_period_from_trajectory(
            samples, -100.0),
        std::domain_error);
    REQUIRE_THROWS_AS(
        astradock::orbit::diagnostics::estimate_period_from_trajectory(
            samples, 0.0),
        std::domain_error);
}

// ============================================================================
// Empirical Convergence Order
// ============================================================================

TEST_CASE("Empirical convergence order throws for invalid inputs",
          "[m05][diagnostics][convergence]") {
    using astradock::orbit::diagnostics::empirical_convergence_order;

    REQUIRE_THROWS_AS(empirical_convergence_order(10.0, 0.001, 20.0, 0.01), std::domain_error);
    REQUIRE_THROWS_AS(empirical_convergence_order(10.0, -1.0, 5.0, 0.5), std::domain_error);
    REQUIRE_THROWS_AS(empirical_convergence_order(10.0, 0.0, 5.0, 0.5), std::domain_error);
    REQUIRE_THROWS_AS(empirical_convergence_order(0.0, 1.0, 0.0, 2.0), std::domain_error);
    REQUIRE_THROWS_AS(empirical_convergence_order(10.0, 1.0, -5.0, 0.5), std::domain_error);
}

TEST_CASE("Empirical order for synthetic first-order data is exactly 1",
          "[m05][diagnostics][convergence]") {
    // If E ~ h^1, halving h halves E.
    const double p = astradock::orbit::diagnostics::empirical_convergence_order(
        0.1, 0.1, 0.05, 0.05);
    REQUIRE(p == Approx(1.0).margin(1.0e-10));
}

TEST_CASE("Empirical order for synthetic fourth-order data is exactly 4",
          "[m05][diagnostics][convergence]") {
    // If E ~ h^4, halving h divides E by 16.
    const double p = astradock::orbit::diagnostics::empirical_convergence_order(
        0.1, 0.0001, 0.05, 0.0001 / 16.0);
    REQUIRE(p == Approx(4.0).margin(1.0e-10));
}

TEST_CASE("RK4 demonstrates approximately fourth-order convergence on orbital dynamics",
          "[m05][diagnostics][convergence][rk4]") {
    const CircularOrbitReference ref = nominal_reference();

    const auto samples_10 = propagate(IntegrationMethod::classical_rk4, 10.0, ref.period_s);
    const auto samples_5 = propagate(IntegrationMethod::classical_rk4, 5.0, ref.period_s);
    const auto samples_2p5 = propagate(IntegrationMethod::classical_rk4, 2.5, ref.period_s);

    const double e10 = astradock::orbit::diagnostics::position_closure_error_m(
        samples_10.front().state, samples_10.back().state);
    const double e5 = astradock::orbit::diagnostics::position_closure_error_m(
        samples_5.front().state, samples_5.back().state);
    const double e2p5 = astradock::orbit::diagnostics::position_closure_error_m(
        samples_2p5.front().state, samples_2p5.back().state);

    REQUIRE(e10 > e5);
    REQUIRE(e5 > e2p5);

    const double p_10_5 = astradock::orbit::diagnostics::empirical_convergence_order(
        10.0, e10, 5.0, e5);
    const double p_5_2p5 = astradock::orbit::diagnostics::empirical_convergence_order(
        5.0, e5, 2.5, e2p5);

    // RK4 approaches fourth order: check it is between 3.8 and 4.2
    REQUIRE(p_10_5 > 3.8);
    REQUIRE(p_10_5 < 4.2);
    REQUIRE(p_5_2p5 > 3.8);
    REQUIRE(p_5_2p5 < 4.2);
}

TEST_CASE("Euler demonstrates approximately first-order convergence on orbital dynamics",
          "[m05][diagnostics][convergence][euler]") {
    const CircularOrbitReference ref = nominal_reference();

    const auto samples_10 = propagate(IntegrationMethod::forward_euler, 10.0, ref.period_s);
    const auto samples_5 = propagate(IntegrationMethod::forward_euler, 5.0, ref.period_s);
    const auto samples_2p5 = propagate(IntegrationMethod::forward_euler, 2.5, ref.period_s);

    const double e10 = astradock::orbit::diagnostics::position_closure_error_m(
        samples_10.front().state, samples_10.back().state);
    const double e5 = astradock::orbit::diagnostics::position_closure_error_m(
        samples_5.front().state, samples_5.back().state);
    const double e2p5 = astradock::orbit::diagnostics::position_closure_error_m(
        samples_2p5.front().state, samples_2p5.back().state);

    REQUIRE(e10 > e5);
    REQUIRE(e5 > e2p5);

    const double p_10_5 = astradock::orbit::diagnostics::empirical_convergence_order(
        10.0, e10, 5.0, e5);
    const double p_5_2p5 = astradock::orbit::diagnostics::empirical_convergence_order(
        5.0, e5, 2.5, e2p5);

    // Euler is first-order: between 0.8 and 1.2
    REQUIRE(p_10_5 > 0.8);
    REQUIRE(p_10_5 < 1.2);
    REQUIRE(p_5_2p5 > 0.8);
    REQUIRE(p_5_2p5 < 1.2);
}

// ============================================================================
// Invariants & Angular-Momentum Direction Drift
// ============================================================================

TEST_CASE("Angular-momentum direction stays nearly fixed for RK4",
          "[m05][diagnostics][angular_momentum]") {
    const CircularOrbitReference ref = nominal_reference();
    const auto samples = propagate(IntegrationMethod::classical_rk4, 10.0, ref.period_s);

    const double drift_rad =
        astradock::orbit::diagnostics::max_angular_momentum_direction_drift_rad(samples);

    // For equatorial orbit, angular momentum vector remains in +Z
    REQUIRE(drift_rad < 1.0e-8);
}

TEST_CASE("Reusable invariant helpers compute accurate errors",
          "[m05][diagnostics][invariants]") {
    const CircularOrbitReference ref = nominal_reference();
    const CartesianState state = nominal_initial_state();

    const double e_err = astradock::orbit::diagnostics::specific_energy_error_m2_per_s2(
        state, ref.specific_energy_m2_per_s2, mu_m3_per_s2);
    REQUIRE(e_err == Approx(0.0).margin(1.0e-12));

    const double rel_e_err = astradock::orbit::diagnostics::relative_energy_error(
        state, ref.specific_energy_m2_per_s2, mu_m3_per_s2);
    REQUIRE(rel_e_err == Approx(0.0).margin(1.0e-15));

    const double h_err = astradock::orbit::diagnostics::specific_angular_momentum_magnitude_error_m2_per_s(
        state, ref.specific_angular_momentum_m2_per_s);
    REQUIRE(h_err == Approx(0.0).margin(1.0e-12));

    const double rel_h_err = astradock::orbit::diagnostics::relative_angular_momentum_magnitude_error(
        state, ref.specific_angular_momentum_m2_per_s);
    REQUIRE(rel_h_err == Approx(0.0).margin(1.0e-15));
}

// ============================================================================
// Repeatability & Determinism
// ============================================================================

TEST_CASE("Same RK4 run is deterministic (byte-identical results)",
          "[m05][diagnostics][repeatability]") {
    const CircularOrbitReference ref = nominal_reference();
    const auto samples_a = propagate(IntegrationMethod::classical_rk4, 10.0, ref.period_s);
    const auto samples_b = propagate(IntegrationMethod::classical_rk4, 10.0, ref.period_s);

    REQUIRE(samples_a.size() == samples_b.size());
    for (std::size_t i = 0; i < samples_a.size(); ++i) {
        REQUIRE(samples_a[i].time_s == samples_b[i].time_s);
        REQUIRE(astradock::math::approximately_equal(samples_a[i].state.position,
                                                      samples_b[i].state.position));
        REQUIRE(astradock::math::approximately_equal(samples_a[i].state.velocity,
                                                      samples_b[i].state.velocity));
    }
}

// ============================================================================
// Regression Baselines
// ============================================================================

TEST_CASE("RK4 one-orbit regression baseline values are maintained",
          "[m05][regression][rk4]") {
    const CircularOrbitReference ref = nominal_reference();
    const auto samples = propagate(IntegrationMethod::classical_rk4, 10.0, ref.period_s);

    const double pos_err = astradock::orbit::diagnostics::position_closure_error_m(
        samples.front().state, samples.back().state);
    const double vel_err = astradock::orbit::diagnostics::velocity_closure_error_m_per_s(
        samples.front().state, samples.back().state);
    const double energy_drift = astradock::orbit::diagnostics::max_relative_energy_drift(
        samples, mu_m3_per_s2);
    const double h_drift = astradock::orbit::diagnostics::max_relative_angular_momentum_drift(samples);

    // M04 baseline: pos_err ~ 0.016 m, vel_err ~ 1.7e-5 m/s, energy_drift ~ 2.9e-11, h_drift ~ 1.4e-11
    REQUIRE(pos_err < 0.1);
    REQUIRE(vel_err < 1.0e-4);
    REQUIRE(energy_drift < 1.0e-9);
    REQUIRE(h_drift < 1.0e-9);
}

TEST_CASE("RK4 five-orbit regression baseline values are maintained",
          "[m05][regression][rk4]") {
    const CircularOrbitReference ref = nominal_reference();
    const auto samples = propagate(IntegrationMethod::classical_rk4, 10.0, 5.0 * ref.period_s);

    const double pos_err = astradock::orbit::diagnostics::position_closure_error_m(
        samples.front().state, samples.back().state);
    const double energy_drift = astradock::orbit::diagnostics::max_relative_energy_drift(
        samples, mu_m3_per_s2);

    // M04 baseline: pos_err ~ 0.098 m, energy_drift ~ 1.4e-10
    REQUIRE(pos_err < 1.0);
    REQUIRE(energy_drift < 1.0e-8);
}

TEST_CASE("Euler still degrades catastrophically (sanity regression)",
          "[m05][regression][euler]") {
    const CircularOrbitReference ref = nominal_reference();
    const auto samples = propagate(IntegrationMethod::forward_euler, 10.0, ref.period_s);
    const double pos_err = astradock::orbit::diagnostics::position_closure_error_m(
        samples.front().state, samples.back().state);
    const double energy_drift = astradock::orbit::diagnostics::max_relative_energy_drift(
        samples, mu_m3_per_s2);

    // Euler should have huge error: millions of metres and > 5% energy drift
    REQUIRE(pos_err > 100'000.0);
    REQUIRE(energy_drift > 0.05);
}

// ============================================================================
// Circular-Orbit Physical Sanity Checks
// ============================================================================

TEST_CASE("RK4 circular orbit radius is approximately constant",
          "[m05][sanity][circular]") {
    const CircularOrbitReference ref = nominal_reference();
    const auto samples = propagate(IntegrationMethod::classical_rk4, 10.0, ref.period_s);
    const double max_rad_dev = astradock::orbit::diagnostics::max_radius_deviation_m(
        samples, orbital_radius_m);

    // Maximum radial deviation should be less than 5 mm for RK4 at 10s
    REQUIRE(max_rad_dev < 0.01);
}

TEST_CASE("RK4 circular orbit speed is approximately constant",
          "[m05][sanity][circular]") {
    const CircularOrbitReference ref = nominal_reference();
    const double expected_speed = ref.speed_m_per_s;
    const auto samples = propagate(IntegrationMethod::classical_rk4, 10.0, ref.period_s);

    for (const auto& s : samples) {
        const double speed = s.state.velocity.norm();
        REQUIRE(speed == Approx(expected_speed).margin(0.01));
    }
}

// ============================================================================
// Timestep Sensitivity: finer dt strictly improves accuracy
// ============================================================================

TEST_CASE("Finer timestep reduces RK4 position closure error",
          "[m05][sensitivity]") {
    const CircularOrbitReference ref = nominal_reference();

    const auto coarse = propagate(IntegrationMethod::classical_rk4, 20.0, ref.period_s);
    const auto fine = propagate(IntegrationMethod::classical_rk4, 2.5, ref.period_s);

    const double coarse_err = astradock::orbit::diagnostics::position_closure_error_m(
        coarse.front().state, coarse.back().state);
    const double fine_err = astradock::orbit::diagnostics::position_closure_error_m(
        fine.front().state, fine.back().state);

    REQUIRE(fine_err < coarse_err);
}

TEST_CASE("Finer timestep reduces Euler position closure error",
          "[m05][sensitivity]") {
    const CircularOrbitReference ref = nominal_reference();

    const auto coarse = propagate(IntegrationMethod::forward_euler, 20.0, ref.period_s);
    const auto fine = propagate(IntegrationMethod::forward_euler, 5.0, ref.period_s);

    const double coarse_err = astradock::orbit::diagnostics::position_closure_error_m(
        coarse.front().state, coarse.back().state);
    const double fine_err = astradock::orbit::diagnostics::position_closure_error_m(
        fine.front().state, fine.back().state);

    REQUIRE(fine_err < coarse_err);
}