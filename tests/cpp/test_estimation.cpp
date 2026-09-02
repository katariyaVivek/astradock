// AstraDock M13A — State Estimation & Extended Kalman Filter test suite.
//
// Layers verified here (per the M13 independent verification standard):
//   Layer 1: analytical linear Kalman cases computed by hand
//   Layer 2: gravity-Jacobian audits against central finite differences,
//            simulation-level convergence/dropout/determinism properties,
//            lightweight Monte Carlo consistency
//   Layer 3 lives in python/audit/independent_ekf_reference.py

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "dynamics/two_body.hpp"
#include "estimation/diagnostics.hpp"
#include "estimation/kalman.hpp"
#include "estimation/translational_ekf.hpp"
#include "math/constants.hpp"
#include "math/linalg.hpp"
#include "numerics/fixed_step_propagation.hpp"
#include "orbit/two_body_orbit.hpp"
#include "sensors/sensor_common.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace astradock;
using Catch::Matchers::WithinAbs;

namespace {

constexpr double kMu = constants::earth_gravitational_parameter_m3_per_s2;

double frobenius_norm(const math::Matrix<3, 3>& m) {
    double sum = 0.0;
    for (std::size_t r = 0; r < 3; ++r) {
        for (std::size_t c = 0; c < 3; ++c) {
            sum += m(r, c) * m(r, c);
        }
    }
    return std::sqrt(sum);
}

math::Vector3 offset_along(const math::Vector3& v, std::size_t axis, double offset) {
    if (axis == 0) {
        return {v.x() + offset, v.y(), v.z()};
    }
    if (axis == 1) {
        return {v.x(), v.y() + offset, v.z()};
    }
    return {v.x(), v.y(), v.z() + offset};
}

// Shared deterministic GNSS scenario used by several simulation-level tests.
struct ScenarioResult {
    std::vector<orbit::CartesianState> truth;
    double initial_position_error_m{0.0};
    std::vector<double> position_error_m;
    std::vector<double> velocity_error_mps;
    std::vector<double> raw_position_error_m;
    std::vector<double> raw_velocity_error_mps;
    std::vector<double> nees;
    std::vector<double> nis;
    std::vector<double> covariance_trace;
    orbit::CartesianState final_estimate{};
};

// Runs one deterministic truth/measurements/EKF pass.
//
// Truth: circular 500 km orbit propagated with RK4, dt = 1 s (identical model
// to the filter's process model -> zero-mismatch baseline).
// Measurements: isotropic Gaussian GNSS noise (sigma_r = 10 m,
// sigma_v = 0.05 m/s), one sample per second from a dedicated seeded RNG.
// Estimator initialization is deliberately offset from truth; the filter never
// receives the truth state directly.
ScenarioResult run_scenario(
    std::uint64_t seed,
    double duration_s,
    bool apply_dropout,
    double dropout_start_s = 3000.0,
    double dropout_end_s = 3600.0) {
    using constants::earth_reference_radius_m;

    const double radius_m = earth_reference_radius_m + 500.0e3;
    const double speed_m_per_s = orbit::circular_orbit_speed_m_per_s(kMu, radius_m);
    const orbit::CartesianState initial_truth{
        {radius_m, 0.0, 0.0},
        {0.0, speed_m_per_s, 0.0},
    };

    // Deliberate estimation initialization error (NOT truth).
    const math::Vector3 initial_position_error_m{50.0e3, -30.0e3, 20.0e3};
    const math::Vector3 initial_velocity_error_mps{5.0, -3.0, 2.0};
    const orbit::CartesianState initial_estimate{
        initial_truth.position + initial_position_error_m,
        initial_truth.velocity + initial_velocity_error_mps,
    };

    estimation::TranslationalCovariance p0 = estimation::TranslationalCovariance::zero();
    for (std::size_t i = 0; i < 3; ++i) {
        p0(i, i) = (25.0e3) * (25.0e3);
        p0(i + 3, i + 3) = 5.0 * 5.0;
    }

    estimation::TranslationalEkfConfig config;
    config.gravitational_parameter_m3_per_s2 = kMu;
    config.acceleration_noise_std_mps2 = 1.0e-3;
    estimation::TranslationalEkf filter(config, initial_estimate, p0);

    sensors::DeterministicRng rng(seed);

    ScenarioResult result;
    result.initial_position_error_m =
        (initial_truth.position - initial_estimate.position).norm();

    const auto truth_samples = numerics::propagate_fixed_step(
        0.0,
        duration_s,
        1.0,
        initial_truth,
        numerics::IntegrationMethod::classical_rk4,
        [](double, const orbit::CartesianState& s) {
            return orbit::two_body_state_derivative(0.0, s, kMu);
        });
    result.truth.reserve(truth_samples.size());
    for (const auto& sample : truth_samples) {
        result.truth.push_back(sample.state);
    }

    for (std::size_t k = 0; k < result.truth.size(); ++k) {
        const double t_s = static_cast<double>(k);
        if (k > 0) {
            filter.predict(1.0);
        }

        const orbit::CartesianState truth_k = result.truth[k];
        const bool in_dropout =
            apply_dropout && t_s >= dropout_start_s && t_s <= dropout_end_s;

        if (!in_dropout) {
            const math::Vector3 position_noise = rng.gaussian_vector3(10.0);
            const math::Vector3 velocity_noise = rng.gaussian_vector3(0.05);
            filter.update_gnss(
                truth_k.position + position_noise,
                truth_k.velocity + velocity_noise);

            const auto& diagnostics = filter.last_update_diagnostics();
            result.raw_position_error_m.push_back(position_noise.norm());
            result.raw_velocity_error_mps.push_back(velocity_noise.norm());
            result.nis.push_back(diagnostics.normalized_innovation_squared);
        }
        // During dropout no RNG draws are consumed: measurements simply vanish.

        const orbit::CartesianState estimate = filter.estimated_state();
        result.position_error_m.push_back(
            (truth_k.position - estimate.position).norm());
        result.velocity_error_mps.push_back(
            (truth_k.velocity - estimate.velocity).norm());
        result.nees.push_back(
            estimation::translational_nees(truth_k, estimate, filter.covariance()));

        const auto diagonal = filter.covariance().diagonal();
        double trace = 0.0;
        for (std::size_t i = 0; i < 6; ++i) {
            trace += diagonal(i, 0);
        }
        result.covariance_trace.push_back(trace);

        result.final_estimate = estimate;
    }
    return result;
}

double root_mean_square(const std::vector<double>& values, std::size_t start_index) {
    if (values.size() <= start_index + 1) {
        throw std::invalid_argument("RMSE window exceeds sample count");
    }
    double sum = 0.0;
    for (std::size_t i = start_index; i < values.size(); ++i) {
        sum += values[i] * values[i];
    }
    return std::sqrt(sum / static_cast<double>(values.size() - start_index));
}

double median_of(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    const std::size_t n = values.size();
    REQUIRE(n > 0);
    if (n % 2 == 1) {
        return values[n / 2];
    }
    return 0.5 * (values[n / 2 - 1] + values[n / 2]);
}

}  // namespace

// ---------------------------------------------------------------------------
// Layer 0: small dense linear algebra
// ---------------------------------------------------------------------------

TEST_CASE("Fixed-size matrix arithmetic matches hand calculations", "[estimation][linalg]") {
    const math::Matrix<2, 2> a(std::array<double, 4>{1.0, 2.0, 3.0, 4.0});
    const math::Matrix<2, 2> b(std::array<double, 4>{5.0, 6.0, 7.0, 8.0});

    const math::Matrix<2, 2> product = a * b;
    CHECK(product(0, 0) == 19.0);
    CHECK(product(0, 1) == 22.0);
    CHECK(product(1, 0) == 43.0);
    CHECK(product(1, 1) == 50.0);

    const math::Matrix<2, 2> transposed = a.transpose();
    CHECK(transposed(0, 1) == 3.0);
    CHECK(transposed(1, 0) == 2.0);

    const math::Matrix<2, 2> identity_check = math::Matrix<2, 2>::identity() * a;
    CHECK(math::max_abs_difference(identity_check, a) == 0.0);

    const math::Matrix<2, 2> asym(std::array<double, 4>{0.0, 1.0, 0.0, 0.0});
    const math::Matrix<2, 2> sym = asym.symmetrized();
    CHECK(sym(0, 1) == 0.5);
    CHECK(sym(1, 0) == 0.5);
    CHECK(sym(0, 0) == 0.0);
    CHECK(asym.symmetry_error() == 1.0);

    const auto diagonal = a.diagonal();
    CHECK(diagonal(0, 0) == 1.0);
    CHECK(diagonal(1, 0) == 4.0);
}

TEST_CASE("Cholesky SPD solves reproduce hand-computed solutions", "[estimation][linalg]") {
    // A = [[4, 2], [2, 3]] (SPD), b = [3, 5]. Hand solution: x = [-0.125, 1.75].
    const math::Matrix<2, 2> a(std::array<double, 4>{4.0, 2.0, 2.0, 3.0});
    const math::Matrix<2, 1> b(std::array<double, 2>{3.0, 5.0});

    const math::Matrix<2, 1> x = math::solve_spd(a, b);
    CHECK_THAT(x(0, 0), WithinAbs(-0.125, 1.0e-15));
    CHECK_THAT(x(1, 0), WithinAbs(1.75, 1.0e-15));

    // Cholesky factor L = [[2, 0], [1, sqrt(2)]], LL^T = A exactly.
    const math::Matrix<2, 2> lower = math::cholesky_lower(a);
    CHECK(lower(0, 0) == 2.0);
    CHECK_THAT(lower(1, 0), WithinAbs(1.0, 1.0e-15));
    CHECK_THAT(lower(1, 1), WithinAbs(std::sqrt(2.0), 1.0e-15));
    CHECK(lower(0, 1) == 0.0);

    // Multi-RHS solve equals repeated single solves.
    const math::Matrix<2, 2> b_multi(std::array<double, 4>{3.0, 1.0, 5.0, -2.0});
    const math::Matrix<2, 2> x_multi = math::solve_spd_multi(a, b_multi);
    const math::Matrix<2, 1> col0 =
        math::solve_spd(a, math::extract_column(b_multi, 0));
    CHECK(x_multi(0, 0) == col0(0, 0));
    CHECK(x_multi(1, 0) == col0(1, 0));
}

TEST_CASE("Non-positive-definite systems fail loudly in Cholesky", "[estimation][linalg]") {
    const math::Matrix<2, 2> indefinite(std::array<double, 4>{1.0, 0.0, 0.0, -1.0});
    const math::Matrix<2, 2> singular(std::array<double, 4>{1.0, 0.0, 0.0, 0.0});
    const math::Matrix<2, 1> rhs(std::array<double, 2>{1.0, 1.0});

    REQUIRE_THROWS_AS(math::solve_spd(indefinite, rhs), std::runtime_error);
    REQUIRE_THROWS_AS(math::solve_spd(singular, rhs), std::runtime_error);

    const math::Matrix<2, 2> nan_matrix(
        std::array<double, 4>{1.0, std::nan(""), 0.0, 1.0});
    REQUIRE_THROWS_AS(math::solve_spd(nan_matrix, rhs), std::domain_error);
}

// ---------------------------------------------------------------------------
// Layer 1: analytical Kalman algebra
// ---------------------------------------------------------------------------

TEST_CASE("Scalar analytical Kalman update matches hand computation",
          "[estimation][kalman]") {
    // State x = 0, P = 1; measurement z = 2, H = 1, R = 1.
    // Hand: S = 2, K = 0.5, x+ = 1, Joseph P+ = (1-K)^2 P + K^2 R = 0.5,
    //       NIS = y^2 / S = 2.
    const math::Matrix<1, 1> prior_state(std::array<double, 1>{0.0});
    const math::Matrix<1, 1> prior_covariance(std::array<double, 1>{1.0});
    const math::Matrix<1, 1> measurement(std::array<double, 1>{2.0});
    const math::Matrix<1, 1> predicted_measurement(std::array<double, 1>{0.0});
    const math::Matrix<1, 1> h(std::array<double, 1>{1.0});
    const math::Matrix<1, 1> r(std::array<double, 1>{1.0});

    const auto result = estimation::ekf_update(
        prior_state, prior_covariance, measurement, predicted_measurement, h, r);

    CHECK_THAT(result.kalman_gain(0, 0), WithinAbs(0.5, 1.0e-15));
    CHECK_THAT(result.posterior_state(0, 0), WithinAbs(1.0, 1.0e-15));
    CHECK_THAT(result.innovation_covariance(0, 0), WithinAbs(2.0, 1.0e-15));
    CHECK_THAT(result.normalized_innovation_squared, WithinAbs(2.0, 1.0e-15));

    // Joseph form agrees with the short-form (I-KH)P in exact arithmetic.
    CHECK_THAT(result.posterior_covariance(0, 0), WithinAbs(0.5, 1.0e-15));
}

TEST_CASE("Two-state scalar-measurement update matches hand computation",
          "[estimation][kalman]") {
    // N=2, M=1, H = [1 0]; P = [[1, 0.5], [0.5, 1]]; R = 1;
    // x_prior = [0; 1]; z = 2 with h(x) = 0.
    // Hand: S = 2; K = [0.5; 0.25]; x+ = [1; 1.5];
    //       Joseph P+ = [[0.5, 0.25], [0.25, 0.875]] (verified on paper).
    const math::Matrix<2, 1> prior_state(std::array<double, 2>{0.0, 1.0});
    const math::Matrix<2, 2> prior_covariance(
        std::array<double, 4>{1.0, 0.5, 0.5, 1.0});
    const math::Matrix<1, 1> measurement(std::array<double, 1>{2.0});
    const math::Matrix<1, 1> predicted_measurement(std::array<double, 1>{0.0});
    const math::Matrix<1, 2> h(std::array<double, 2>{1.0, 0.0});
    const math::Matrix<1, 1> r(std::array<double, 1>{1.0});

    const auto result = estimation::ekf_update(
        prior_state, prior_covariance, measurement, predicted_measurement, h, r);

    CHECK_THAT(result.kalman_gain(0, 0), WithinAbs(0.5, 1.0e-15));
    CHECK_THAT(result.kalman_gain(1, 0), WithinAbs(0.25, 1.0e-15));
    CHECK_THAT(result.posterior_state(0, 0), WithinAbs(1.0, 1.0e-15));
    CHECK_THAT(result.posterior_state(1, 0), WithinAbs(1.5, 1.0e-15));
    CHECK_THAT(result.posterior_covariance(0, 0), WithinAbs(0.5, 1.0e-15));
    CHECK_THAT(result.posterior_covariance(0, 1), WithinAbs(0.25, 1.0e-15));
    CHECK_THAT(result.posterior_covariance(1, 0), WithinAbs(0.25, 1.0e-15));
    CHECK_THAT(result.posterior_covariance(1, 1), WithinAbs(0.875, 1.0e-15));
    CHECK(result.posterior_covariance.symmetry_error() == 0.0);
}

TEST_CASE("Linear covariance prediction reproduces exact discrete propagation",
          "[estimation][kalman]") {
    // Constant-velocity system: Phi = [[1, 0.5], [0, 1]] (dt = 0.5),
    // white-acceleration Q = sigma_a^2 [[t^3/3, t^2/2], [t^2/2, t]],
    // sigma_a = 2, P0 = diag(4, 9), x0 = [1; 2].
    // Hand: Phi P0 Phi^T = [[6.25, 4.5], [4.5, 9]];
    //       Q = [[1/6, 1/2], [1/2, 2]];
    //       P- = [[6.25 + 1/6, 5], [5, 11]]; mean Phi x = [2; 2].
    const math::Matrix<2, 2> phi(std::array<double, 4>{1.0, 0.5, 0.0, 1.0});
    const double sigma_a = 2.0;
    const double dt = 0.5;
    const math::Matrix<2, 2> q(std::array<double, 4>{
        sigma_a * sigma_a * dt * dt * dt / 3.0,
        sigma_a * sigma_a * dt * dt / 2.0,
        sigma_a * sigma_a * dt * dt / 2.0,
        sigma_a * sigma_a * dt,
    });
    const math::Matrix<2, 2> p0(std::array<double, 4>{4.0, 0.0, 0.0, 9.0});
    const math::Matrix<2, 1> x0(std::array<double, 2>{1.0, 2.0});

    const auto prediction = estimation::ekf_predict(
        x0,
        p0,
        [&phi](const math::Matrix<2, 1>& state) { return phi * state; },
        phi,
        q);

    CHECK_THAT(prediction.prior_state(0, 0), WithinAbs(2.0, 1.0e-15));
    CHECK_THAT(prediction.prior_state(1, 0), WithinAbs(2.0, 1.0e-15));
    CHECK_THAT(
        prediction.prior_covariance(0, 0), WithinAbs(6.25 + 1.0 / 6.0, 1.0e-15));
    CHECK_THAT(prediction.prior_covariance(0, 1), WithinAbs(5.0, 1.0e-15));
    CHECK_THAT(prediction.prior_covariance(1, 0), WithinAbs(5.0, 1.0e-15));
    CHECK_THAT(prediction.prior_covariance(1, 1), WithinAbs(11.0, 1.0e-15));
}

TEST_CASE("Measurement-noise limiting behavior trusts or rejects sensors",
          "[estimation][kalman]") {
    // Huge R: the measurement barely moves the estimate.
    // Tiny R: the estimate snaps to the measurement.
    const math::Matrix<1, 1> prior_state(std::array<double, 1>{0.0});
    const math::Matrix<1, 1> prior_covariance(std::array<double, 1>{1.0});
    const math::Matrix<1, 1> measurement(std::array<double, 1>{10.0});
    const math::Matrix<1, 1> predicted_measurement(std::array<double, 1>{0.0});
    const math::Matrix<1, 1> h(std::array<double, 1>{1.0});
    const math::Matrix<1, 1> r_huge(std::array<double, 1>{1.0e12});
    const math::Matrix<1, 1> r_tiny(std::array<double, 1>{1.0e-12});

    const auto weak = estimation::ekf_update(
        prior_state, prior_covariance, measurement, predicted_measurement, h, r_huge);
    CHECK_THAT(weak.posterior_state(0, 0), WithinAbs(1.0e-11, 1.0e-12));

    const auto strong = estimation::ekf_update(
        prior_state, prior_covariance, measurement, predicted_measurement, h, r_tiny);
    CHECK_THAT(strong.posterior_state(0, 0), WithinAbs(10.0, 1.0e-6));
}

TEST_CASE("Process-noise limiting behavior controls prediction trust",
          "[estimation][kalman]") {
    // Scalar system Phi = 0.5. Zero Q => covariance shrinks through Phi P Phi^T
    // only; large Q => large predicted covariance => larger later Kalman gain.
    const math::Matrix<1, 1> phi(std::array<double, 1>{0.5});
    const auto identity_model =
        [&phi](const math::Matrix<1, 1>& state) { return phi * state; };
    const math::Matrix<1, 1> x(std::array<double, 1>{2.0});
    const math::Matrix<1, 1> p(std::array<double, 1>{1.0});
    const math::Matrix<1, 1> q_zero(std::array<double, 1>{0.0});
    const math::Matrix<1, 1> q_large(std::array<double, 1>{10.0});

    const auto no_noise = estimation::ekf_predict(x, p, identity_model, phi, q_zero);
    // Exact hand value: P- = 0.25 * 1 = 0.25.
    CHECK(no_noise.prior_covariance(0, 0) == 0.25);

    const auto noisy = estimation::ekf_predict(x, p, identity_model, phi, q_large);
    CHECK(noisy.prior_covariance(0, 0) == 10.25);

    // Same subsequent measurement: the high-Q filter weighs it more strongly.
    const auto gain_from_low_q = estimation::ekf_update(
        no_noise.prior_state,
        no_noise.prior_covariance,
        math::Matrix<1, 1>(std::array<double, 1>{1.0}),
        no_noise.prior_state,
        phi,
        math::Matrix<1, 1>(std::array<double, 1>{1.0}));
    const auto gain_from_high_q = estimation::ekf_update(
        noisy.prior_state,
        noisy.prior_covariance,
        math::Matrix<1, 1>(std::array<double, 1>{1.0}),
        noisy.prior_state,
        phi,
        math::Matrix<1, 1>(std::array<double, 1>{1.0}));
    CHECK(gain_from_high_q.kalman_gain(0, 0) > gain_from_low_q.kalman_gain(0, 0));
    CHECK(gain_from_high_q.kalman_gain(0, 0) > 0.9);
    CHECK(gain_from_low_q.kalman_gain(0, 0) < 0.25);
}

TEST_CASE("Invalid numerical filter states are rejected, never repaired silently",
          "[estimation][kalman]") {
    const math::Matrix<1, 1> zero(std::array<double, 1>{0.0});
    const math::Matrix<1, 1> nan_value(std::array<double, 1>{std::nan("")});
    const math::Matrix<1, 1> negative_variance(std::array<double, 1>{-1.0});

    REQUIRE_THROWS_AS(
        estimation::ekf_update(zero, zero, nan_value, zero, zero, zero),
        std::domain_error);
    REQUIRE_THROWS_AS(
        estimation::ekf_update(zero, negative_variance, zero, zero, zero, zero),
        std::domain_error);
    // Singular innovation covariance (P = 0, R = 0) fails inside Cholesky.
    REQUIRE_THROWS_AS(
        estimation::ekf_update(zero, zero, zero, zero, zero, zero),
        std::runtime_error);
    REQUIRE_THROWS_AS(
        estimation::ekf_predict(
            zero, zero, [](auto& s) { return s; }, nan_value, zero),
        std::domain_error);
}

// ---------------------------------------------------------------------------
// Gravity Jacobian: analytical structure and independent audit
// ---------------------------------------------------------------------------

TEST_CASE("Analytical gravity Jacobian matches known directional structures",
          "[estimation][jacobian]") {
    // Axis-aligned position: radial eigenvalue +2 mu/r^3, transverse -mu/r^3.
    const double r = 7000.0e3;
    const math::Vector3 axis_aligned{r, 0.0, 0.0};
    const auto g_axis = estimation::gravity_position_jacobian(axis_aligned, kMu);

    const double radial_value = 2.0 * kMu / (r * r * r);
    const double transverse_value = -kMu / (r * r * r);
    CHECK_THAT(g_axis(0, 0), WithinAbs(radial_value, 1.0e-21));
    CHECK(g_axis(1, 1) == transverse_value);
    CHECK(g_axis(2, 2) == transverse_value);
    CHECK(g_axis(0, 1) == 0.0);
    CHECK(g_axis(0, 2) == 0.0);
    CHECK(g_axis(1, 2) == 0.0);

    // Body-diagonal position: G = -mu/r^3 (I - J) => zero diagonal,
    // +mu/r^3 off-diagonals.
    const double s = r / std::sqrt(3.0);
    const auto g_diagonal_case = estimation::gravity_position_jacobian(
        math::Vector3{s, s, s}, kMu);
    const double expected_off_diagonal = kMu / (r * r * r);
    for (std::size_t i = 0; i < 3; ++i) {
        CHECK(g_diagonal_case(i, i) == 0.0);
        for (std::size_t j = 0; j < 3; ++j) {
            if (i != j) {
                CHECK_THAT(
                    g_diagonal_case(i, j),
                    WithinAbs(expected_off_diagonal, 1.0e-21));
            }
        }
    }

    // General invariants: exactly symmetric, traceless (harmonic potential).
    const auto g_general = estimation::gravity_position_jacobian(
        math::Vector3{2000.0e3, -4500.0e3, 6100.0e3}, kMu);
    CHECK(g_general.symmetry_error() == 0.0);
    const auto general_diagonal = g_general.diagonal();
    const double trace =
        general_diagonal(0, 0) + general_diagonal(1, 0) + general_diagonal(2, 0);
    // Floating-point unit-vector norm is not exactly 1, so trace departs
    // from zero by O(eps * mu / r^3) ~ 1e-22 for LEO-scale positions.
    CHECK(std::abs(trace) < 1.0e-18);

    REQUIRE_THROWS_AS(
        estimation::gravity_position_jacobian(math::Vector3{0.0, 0.0, 0.0}, kMu),
        std::domain_error);
    REQUIRE_THROWS_AS(
        estimation::gravity_position_jacobian(axis_aligned, -kMu),
        std::domain_error);
    REQUIRE_THROWS_AS(
        estimation::discrete_process_noise(-1.0, 1.0),
        std::domain_error);
}

TEST_CASE("Analytical gravity Jacobian survives an independent finite-difference "
          "audit across positions and step scales",
          "[estimation][jacobian][audit]") {
    // Central-difference truncation error falls as O(h^2) until floating-point
    // cancellation error (~ eps/h) dominates. Sweeping the step documents this
    // tradeoff instead of trusting one arbitrary epsilon.
    const std::vector<math::Vector3> positions{
        math::Vector3{7000.0e3, 0.0, 0.0},
        math::Vector3{-3000.0e3, 5000.0e3, 8000.0e3},
        math::Vector3{1234.0e3, -6789.0e3, 2345.0e3},
        math::Vector3{6543.0e3, 2109.0e3, -4444.0e3},
    };
    const std::vector<double> steps_m{
        1.0e6, 1.0e5, 1.0e4, 1.0e3, 1.0e2, 10.0, 1.0, 0.1};

    double worst_best_relative_error = 0.0;

    for (const auto& position : positions) {
        const auto analytical = estimation::gravity_position_jacobian(position, kMu);
        const double analytical_norm = frobenius_norm(analytical);

        double best_relative_error = std::numeric_limits<double>::infinity();
        double best_step_m = 0.0;
        for (const double h : steps_m) {
            math::Matrix<3, 3> numerical = math::Matrix<3, 3>::zero();
            for (std::size_t j = 0; j < 3; ++j) {
                const math::Vector3 acceleration_plus =
                    dynamics::two_body_acceleration(offset_along(position, j, h), kMu);
                const math::Vector3 acceleration_minus =
                    dynamics::two_body_acceleration(offset_along(position, j, -h), kMu);
                const math::Vector3 column =
                    (acceleration_plus - acceleration_minus) / (2.0 * h);
                numerical(0, j) = column.x();
                numerical(1, j) = column.y();
                numerical(2, j) = column.z();
            }
            const double relative_error =
                frobenius_norm(analytical - numerical) / analytical_norm;
            if (relative_error < best_relative_error) {
                best_relative_error = relative_error;
                best_step_m = h;
            }
        }
        INFO("Position " << position.x() << ", " << position.y() << ", " << position.z()
             << ": best step " << best_step_m << " m gives rel err "
             << best_relative_error);
        worst_best_relative_error =
            std::max(worst_best_relative_error, best_relative_error);
    }

    CHECK(worst_best_relative_error < 1.0e-10);
}

// ---------------------------------------------------------------------------
// Translational EKF domain layer
// ---------------------------------------------------------------------------

TEST_CASE("GNSS measurement model is the identity on the documented ordering",
          "[estimation][translational]") {
    const auto h = estimation::gnss_measurement_matrix();
    const auto identity = estimation::TranslationalCovariance::identity();
    CHECK(math::max_abs_difference(h, identity) == 0.0);

    // Process-noise builder produces the exact analytic block structure.
    const double dt = 2.0;
    const double sigma = 0.01;
    const auto q = estimation::discrete_process_noise(sigma, dt);
    const double variance = sigma * sigma;
    CHECK_THAT(q(0, 0), WithinAbs(variance * dt * dt * dt / 3.0, 1.0e-30));
    CHECK_THAT(q(3, 3), WithinAbs(variance * dt, 1.0e-30));
    CHECK_THAT(q(0, 3), WithinAbs(variance * dt * dt / 2.0, 1.0e-30));
    CHECK(q(0, 4) == 0.0);
    CHECK(q(1, 5) == 0.0);
    CHECK(q(0, 1) == 0.0);
    CHECK(q(3, 4) == 0.0);

    // First-order discrete transition of the two-body dynamics at an
    // axis-aligned state: top-right I*dt, bottom-left G*dt.
    const double r = 7000.0e3;
    const auto phi = estimation::discrete_state_transition(
        estimation::continuous_dynamics_matrix(math::Vector3{r, 0.0, 0.0}, kMu), 5.0);
    CHECK(phi(0, 3) == 5.0);
    CHECK(phi(1, 4) == 5.0);
    CHECK(phi(2, 5) == 5.0);
    CHECK_THAT(
        phi(3, 0),
        WithinAbs((2.0 * kMu / (r * r * r)) * 5.0, 1.0e-24));
    CHECK_THAT(
        phi(4, 1),
        WithinAbs((-kMu / (r * r * r)) * 5.0, 1.0e-24));
    CHECK(phi(3, 1) == 0.0);
}

TEST_CASE("Predict-only filter mean propagation is bitwise identical to open-loop "
          "RK4 orbit propagation",
          "[estimation][translational]") {
    using constants::earth_reference_radius_m;
    const double radius_m = earth_reference_radius_m + 500.0e3;
    const orbit::CartesianState initial{
        {radius_m, 0.0, 0.0},
        {0.0, orbit::circular_orbit_speed_m_per_s(kMu, radius_m), 0.0},
    };

    estimation::TranslationalEkfConfig config;
    config.gravitational_parameter_m3_per_s2 = kMu;
    estimation::TranslationalEkf filter(
        config, initial, estimation::TranslationalCovariance::identity());

    const auto reference = numerics::propagate_fixed_step(
        0.0, 600.0, 1.0, initial,
        numerics::IntegrationMethod::classical_rk4,
        [](double, const orbit::CartesianState& s) {
            return orbit::two_body_state_derivative(0.0, s, kMu);
        });

    // Same integrator, same model, same operation order -> bitwise equality
    // at every step: the estimator reuses M04 machinery, it does not fork it.
    for (std::size_t step = 1; step < reference.size(); ++step) {
        filter.predict(1.0);
        const auto& expected = reference[step].state;
        const orbit::CartesianState estimated = filter.estimated_state();
        CHECK(estimated.position.x() == expected.position.x());
        CHECK(estimated.position.y() == expected.position.y());
        CHECK(estimated.position.z() == expected.position.z());
        CHECK(estimated.velocity.x() == expected.velocity.x());
        CHECK(estimated.velocity.y() == expected.velocity.y());
        CHECK(estimated.velocity.z() == expected.velocity.z());
    }
}

// ---------------------------------------------------------------------------
// Layer 2: simulation-level convergence, consistency, robustness
// ---------------------------------------------------------------------------

TEST_CASE("GNSS EKF converges from a deliberately wrong initialization and beats "
          "raw measurements",
          "[estimation][translational][convergence]") {
    // ~1.06 orbits. Initial error ~62 km / 6.2 m/s vs sigma_r = 10 m GNSS.
    const ScenarioResult run = run_scenario(42u, 6000.0, false);
    const std::size_t n = run.position_error_m.size();

    // The estimator starts deliberately wrong and ends far closer to truth.
    CHECK(run.initial_position_error_m > 60.0e3);
    const double final_position_error = run.position_error_m[n - 1];
    const double final_velocity_error = run.velocity_error_mps[n - 1];

    // Raw GNSS statistics over the converged window (last half).
    const std::size_t half = n / 2;
    const double ekf_position_rmse =
        root_mean_square(run.position_error_m, half);
    const double raw_position_rmse =
        root_mean_square(run.raw_position_error_m, half);
    const double ekf_velocity_rmse =
        root_mean_square(run.velocity_error_mps, half);
    const double raw_velocity_rmse =
        root_mean_square(run.raw_velocity_error_mps, half);

    INFO("final pos err " << final_position_error << " m, final vel err "
         << final_velocity_error << " m/s");
    INFO("EKF pos RMSE " << ekf_position_rmse << " m vs raw " << raw_position_rmse
         << " m; EKF vel RMSE " << ekf_velocity_rmse << " vs raw "
         << raw_velocity_rmse);

    CHECK(final_position_error < 100.0);
    CHECK(final_velocity_error < 0.20);
    CHECK(ekf_position_rmse < 0.75 * raw_position_rmse);
    CHECK(ekf_velocity_rmse < 0.75 * raw_velocity_rmse);

    // Single-run NEES bound (generous; statistical rigor comes from the
    // Monte Carlo study).
    CHECK(run.nees[n - 1] < 60.0);

    // Innovation consistency over the converged window: NIS ~ chi-square(6).
    double mean_nis = 0.0;
    const std::size_t nis_window_start =
        run.nis.size() > 4000 ? run.nis.size() - 4000 : 0;
    for (std::size_t i = nis_window_start; i < run.nis.size(); ++i) {
        mean_nis += run.nis[i];
    }
    mean_nis /= static_cast<double>(run.nis.size() - nis_window_start);
    INFO("mean NIS (converged window): " << mean_nis);
    CHECK(mean_nis > 4.0);
    CHECK(mean_nis < 8.0);
}

TEST_CASE("NEES computation matches a hand-built diagonal case",
          "[estimation][diagnostics]") {
    const orbit::CartesianState truth{
        {7000.0e3 + 10.0, 0.0, 0.0},
        {0.0, 7500.0, 0.0},
    };
    const orbit::CartesianState estimate{
        {7000.0e3, 0.0, 0.0},
        {0.0, 7500.0, 0.0},
    };
    estimation::TranslationalCovariance p = estimation::TranslationalCovariance::zero();
    p(0, 0) = 100.0;  // position sigma_x = 10 m
    p(1, 1) = 400.0;
    p(2, 2) = 900.0;
    p(3, 3) = 4.0;
    p(4, 4) = 9.0;
    p(5, 5) = 16.0;

    // Error is 10 m along x only: NEES = 100 / 100 = 1.
    CHECK_THAT(estimation::translational_nees(truth, estimate, p), WithinAbs(1.0, 1.0e-12));

    // Doubling the error quadruples NEES.
    const orbit::CartesianState doubled_truth{
        {7000.0e3 + 20.0, 0.0, 0.0},
        {0.0, 7500.0, 0.0},
    };
    CHECK_THAT(
        estimation::translational_nees(doubled_truth, estimate, p),
        WithinAbs(4.0, 1.0e-12));

    REQUIRE_THROWS_AS(
        estimation::translational_nees(truth, estimate,
                                       estimation::TranslationalCovariance::identity()
                                           * -1.0),
        std::runtime_error);
}

TEST_CASE("GNSS dropout grows covariance and recovery restores information",
          "[estimation][translational][dropout]") {
    // 900 s outage inside a 6000 s run.
    const ScenarioResult run = run_scenario(1234u, 6000.0, true, 3000.0, 3900.0);

    const auto index_at = [&](double t) {
        return static_cast<std::size_t>(t);
    };

    const double trace_before = run.covariance_trace[index_at(2999)];
    const double trace_during_early = run.covariance_trace[index_at(3300)];
    const double trace_during_late = run.covariance_trace[index_at(3899)];

    // Covariance grows monotonically during pure prediction.
    CHECK(trace_during_early > trace_before);
    CHECK(trace_during_late > trace_during_early);
    bool monotone = true;
    for (std::size_t i = index_at(3001); i <= index_at(3899); ++i) {
        if (run.covariance_trace[i] <= run.covariance_trace[i - 1]) {
            monotone = false;
        }
    }
    CHECK(monotone);
    CHECK(trace_during_late > 1.15 * trace_before);

    // After recovery the filter re-absorbs information quickly.
    const double trace_recovered = run.covariance_trace[index_at(4200)];
    CHECK(trace_recovered < trace_during_late);
    CHECK(trace_recovered < 1.10 * trace_before);

    // The estimate does not diverge during the outage.
    double max_dropout_position_error = 0.0;
    for (std::size_t i = index_at(3000); i <= index_at(3900); ++i) {
        max_dropout_position_error =
            std::max(max_dropout_position_error, run.position_error_m[i]);
    }
    INFO("Max position error during dropout: " << max_dropout_position_error);
    CHECK(max_dropout_position_error < 500.0);

    // Post-recovery performance returns to nominal levels.
    const double post_recovery_rmse =
        root_mean_square(run.position_error_m, index_at(4200));
    CHECK(post_recovery_rmse < 25.0);
}

TEST_CASE("Truth non-interference: estimator consumes copies and never mutates "
          "truth or measurement history",
          "[estimation][noninterference]") {
    // Generate the truth trajectory and measurement stream once.
    using constants::earth_reference_radius_m;
    const double radius_m = earth_reference_radius_m + 500.0e3;
    const orbit::CartesianState initial_truth{
        {radius_m, 0.0, 0.0},
        {0.0, orbit::circular_orbit_speed_m_per_s(kMu, radius_m), 0.0},
    };
    const auto samples = numerics::propagate_fixed_step(
        0.0, 300.0, 1.0, initial_truth,
        numerics::IntegrationMethod::classical_rk4,
        [](double, const orbit::CartesianState& s) {
            return orbit::two_body_state_derivative(0.0, s, kMu);
        });

    sensors::DeterministicRng rng(555u);
    std::vector<std::pair<math::Vector3, math::Vector3>> measurements;
    std::vector<orbit::CartesianState> truth_copy;
    for (const auto& sample : samples) {
        measurements.emplace_back(
            sample.state.position + rng.gaussian_vector3(10.0),
            sample.state.velocity + rng.gaussian_vector3(0.05));
        truth_copy.push_back(sample.state);
    }

    // Run the estimator over copies of the data.
    estimation::TranslationalEkfConfig config;
    config.gravitational_parameter_m3_per_s2 = kMu;
    estimation::TranslationalEkf filter(
        config,
        orbit::CartesianState{
            initial_truth.position + math::Vector3{1000.0, -2000.0, 500.0},
            initial_truth.velocity},
        estimation::TranslationalCovariance::identity() * 1.0e6);
    for (std::size_t k = 0; k < truth_copy.size(); ++k) {
        if (k > 0) {
            filter.predict(1.0);
        }
        static_cast<void>(filter.update_gnss(measurements[k].first, measurements[k].second));
    }

    // The original truth telemetry must be bitwise untouched by estimation.
    for (std::size_t k = 0; k < truth_copy.size(); ++k) {
        CHECK(samples[k].state.position.x() == truth_copy[k].position.x());
        CHECK(samples[k].state.position.y() == truth_copy[k].position.y());
        CHECK(samples[k].state.position.z() == truth_copy[k].position.z());
        CHECK(samples[k].state.velocity.x() == truth_copy[k].velocity.x());
        CHECK(samples[k].state.velocity.y() == truth_copy[k].velocity.y());
        CHECK(samples[k].state.velocity.z() == truth_copy[k].velocity.z());
    }
}

TEST_CASE("Estimator runs are bitwise deterministic for identical seeds",
          "[estimation][determinism]") {
    const ScenarioResult first = run_scenario(777u, 1200.0, false);
    const ScenarioResult second = run_scenario(777u, 1200.0, false);

    REQUIRE(first.truth.size() == second.truth.size());
    for (std::size_t i = 0; i < first.truth.size(); ++i) {
        CHECK(first.truth[i].position.x() == second.truth[i].position.x());
        CHECK(first.truth[i].velocity.z() == second.truth[i].velocity.z());
    }
    for (std::size_t i = 0; i < first.position_error_m.size(); ++i) {
        CHECK(first.position_error_m[i] == second.position_error_m[i]);
        CHECK(first.nis[i] == second.nis[i]);
        CHECK(first.nees[i] == second.nees[i]);
    }
}

TEST_CASE("Monte Carlo estimation study shows consistent, beating-the-sensor "
          "performance across seeds",
          "[estimation][montecarlo]") {
    constexpr int kRuns = 40;  // CI-friendly count; the demo tool runs 100.
    std::vector<double> final_errors;
    std::vector<double> nees_finals;
    int runs_beating_raw = 0;

    for (int run = 0; run < kRuns; ++run) {
        const auto seed = static_cast<std::uint64_t>(10000u + run);
        const ScenarioResult result = run_scenario(seed, 1800.0, false);
        const std::size_t n = result.position_error_m.size();
        const std::size_t half = n / 2;

        final_errors.push_back(result.position_error_m[n - 1]);
        nees_finals.push_back(result.nees[n - 1]);

        const double ekf_rmse = root_mean_square(result.position_error_m, half);
        const double raw_rmse = root_mean_square(result.raw_position_error_m, half);
        if (ekf_rmse < raw_rmse) {
            ++runs_beating_raw;
        }
    }

    const double median_final_error = median_of(final_errors);
    const double median_nees = median_of(nees_finals);
    const double beat_fraction =
        static_cast<double>(runs_beating_raw) / static_cast<double>(kRuns);

    INFO("median final error " << median_final_error << " m; median NEES "
         << median_nees << "; beat fraction " << beat_fraction);

    // The filter should reliably outperform the raw sensor.
    CHECK(median_final_error < 17.32);  // < sqrt(3) * sigma_r
    CHECK(beat_fraction >= 0.95);

    // Median NEES across seeds must sit inside the central chi-square(6)
    // band: neither overconfident nor wildly conservative.
    CHECK(median_nees > estimation::k_chi2_6dof_95_lower);
    CHECK(median_nees < estimation::k_chi2_6dof_95_upper);
}
