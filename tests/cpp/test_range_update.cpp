#include "estimation/diagnostics.hpp"
#include "estimation/range_update.hpp"
#include "estimation/translational_ekf.hpp"
#include "math/constants.hpp"
#include "math/vector3.hpp"
#include "sensors/sensor_common.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using namespace astradock;

TEST_CASE("Range update analytical test matches hand-calculated Kalman solution to machine precision", "[estimation][range][analytical]") {
    // Problem setup:
    // Target position: r_t = [1000, 0, 0] m
    // Spacecraft prior estimate: r_s = [0, 0, 0] m, v_s = [0, 0, 0] m/s
    // True spacecraft position: r_true = [10, 0, 0] m (10 m offset along +X)
    // Predicted range: rho_pred = 1000 m
    // Measured range: rho_meas = 990 m (zero noise)
    // Innovation: y = rho_meas - rho_pred = -10 m

    const math::Vector3 target_pos{1000.0, 0.0, 0.0};
    math::ColVector<6> prior_state{}; // all zeros

    // Diagonal prior covariance: sigma_r = 20 m (variance = 400 m^2), sigma_v = 1 m/s (variance = 1 m^2/s^2)
    math::Matrix<6, 6> prior_cov;
    prior_cov(0, 0) = 400.0;
    prior_cov(1, 1) = 400.0;
    prior_cov(2, 2) = 400.0;
    prior_cov(3, 3) = 1.0;
    prior_cov(4, 4) = 1.0;
    prior_cov(5, 5) = 1.0;

    const double range_noise_std = 2.0; // R = 4.0 m^2
    const double measured_range = 990.0;

    estimation::RangeUpdateDiagnostics diag;
    const auto result = estimation::execute_range_update(
        prior_state,
        prior_cov,
        measured_range,
        target_pos,
        range_noise_std,
        &diag
    );

    // Hand calculations:
    // H = [ -1, 0, 0, 0, 0, 0 ]
    // S = H P H^T + R = (-1)^2 * 400 + 4 = 404.0 m^2
    const double expected_S = 404.0;
    CHECK_THAT(diag.innovation_variance_m2, WithinAbs(expected_S, 1.0e-12));

    // y = -10.0 m
    CHECK_THAT(diag.innovation_m, WithinAbs(-10.0, 1.0e-12));

    // NIS = y^2 / S = 100 / 404
    const double expected_NIS = 100.0 / 404.0;
    CHECK_THAT(diag.normalized_innovation_squared, WithinAbs(expected_NIS, 1.0e-12));

    // K = P H^T / S = [ -400 / 404, 0, 0, 0, 0, 0 ]^T
    const double expected_K0 = -400.0 / 404.0;
    CHECK_THAT(result.kalman_gain(0, 0), WithinAbs(expected_K0, 1.0e-14));
    CHECK_THAT(result.kalman_gain(1, 0), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(result.kalman_gain(2, 0), WithinAbs(0.0, 1.0e-15));

    // Posterior state: x^+ = x^- + K y = 0 + (-400/404) * (-10) = 4000 / 404 ≈ 9.90099 m
    const double expected_x_post = 4000.0 / 404.0;
    CHECK_THAT(result.posterior_state(0, 0), WithinAbs(expected_x_post, 1.0e-14));
    CHECK_THAT(result.posterior_state(1, 0), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(result.posterior_state(2, 0), WithinAbs(0.0, 1.0e-15));

    // Posterior variance: P_00^+ = 400 / 101 ≈ 3.9603960396 m^2
    const double expected_P00 = 400.0 / 101.0;
    CHECK_THAT(result.posterior_covariance(0, 0), WithinAbs(expected_P00, 1.0e-12));
    // Cross-track variances must be completely unchanged by this scalar line-of-sight update
    CHECK_THAT(result.posterior_covariance(1, 1), WithinAbs(400.0, 1.0e-12));
    CHECK_THAT(result.posterior_covariance(2, 2), WithinAbs(400.0, 1.0e-12));
}

TEST_CASE("Range measurement Jacobian survives finite-difference audit across arbitrary non-axis-aligned geometry", "[estimation][range][jacobian]") {
    const std::vector<std::pair<math::Vector3, math::Vector3>> test_pairs = {
        {{6878137.0, 1000.0, -500.0}, {6878137.0 + 300.0, 1000.0 - 400.0, -500.0 + 1200.0}},
        {{7000.0e3, 5000.0e3, 2000.0e3}, {7000.0e3 + 50.0, 5000.0e3 - 30.0, 2000.0e3 + 80.0}},
        {{-12345.0, 67890.0, 45678.0}, {-12300.0, 67900.0, 45650.0}}
    };

    const double eps = 0.01; // 1 cm perturbation for robust finite difference on orbital scales

    for (const auto& [sc_pos, target_pos] : test_pairs) {
        const auto H_analytical = estimation::range_measurement_jacobian(sc_pos, target_pos);

        // Audit position derivatives against central finite differences
        for (std::size_t j = 0; j < 3; ++j) {
            const math::Vector3 pert_pos{
                sc_pos.x() + (j == 0 ? eps : 0.0),
                sc_pos.y() + (j == 1 ? eps : 0.0),
                sc_pos.z() + (j == 2 ? eps : 0.0)
            };
            const math::Vector3 pert_neg{
                sc_pos.x() - (j == 0 ? eps : 0.0),
                sc_pos.y() - (j == 1 ? eps : 0.0),
                sc_pos.z() - (j == 2 ? eps : 0.0)
            };

            const double rho_pos = (target_pos - pert_pos).norm();
            const double rho_neg = (target_pos - pert_neg).norm();
            const double fd_deriv = (rho_pos - rho_neg) / (2.0 * eps);

            // Note: H is 1x6
            CHECK_THAT(H_analytical(0, j), WithinAbs(fd_deriv, 1.0e-6));
        }

        // Velocity derivatives MUST be identically zero
        CHECK(H_analytical(0, 3) == 0.0);
        CHECK(H_analytical(0, 4) == 0.0);
        CHECK(H_analytical(0, 5) == 0.0);
    }
}

TEST_CASE("Range measurement and Jacobian satisfy exact translation invariance", "[estimation][range][invariance]") {
    const math::Vector3 sc_pos{6878137.0, 1000.0, 2000.0};
    const math::Vector3 target_pos{6878137.0 + 500.0, 1000.0 - 200.0, 2000.0 + 350.0};

    const double rho_base = estimation::predicted_range(sc_pos, target_pos);
    const auto H_base = estimation::range_measurement_jacobian(sc_pos, target_pos);

    // Shift both bodies by large arbitrary translation vector c
    const math::Vector3 shift{-9.876e6, 5.432e6, -1.234e5};
    const double rho_shifted = estimation::predicted_range(sc_pos + shift, target_pos + shift);
    const auto H_shifted = estimation::range_measurement_jacobian(sc_pos + shift, target_pos + shift);

    CHECK_THAT(rho_shifted, WithinAbs(rho_base, 1.0e-12));
    for (std::size_t c = 0; c < 6; ++c) {
        CHECK_THAT(H_shifted(0, c), WithinAbs(H_base(0, c), 1.0e-15));
    }
}

TEST_CASE("Range measurement rejects singular geometry when spacecraft and target coincide", "[estimation][range][singularity]") {
    const math::Vector3 pos{7000.0e3, 0.0, 0.0};

    // Exactly identical positions
    CHECK_THROWS_AS(estimation::predicted_range(pos, pos), std::domain_error);
    CHECK_THROWS_AS(estimation::range_measurement_jacobian(pos, pos), std::domain_error);

    // Coinciding within threshold (< 1e-6 m)
    const math::Vector3 tiny_offset{1.0e-7, 0.0, 0.0};
    CHECK_THROWS_AS(estimation::predicted_range(pos, pos + tiny_offset), std::domain_error);
    CHECK_THROWS_AS(estimation::range_measurement_jacobian(pos, pos + tiny_offset), std::domain_error);
}

TEST_CASE("TranslationalEkf update_range integrates cleanly with existing filter state", "[estimation][range][integration]") {
    estimation::TranslationalEkfConfig config;
    config.gravitational_parameter_m3_per_s2 = constants::earth_gravitational_parameter_m3_per_s2;
    const math::Vector3 sc_pos{6878137.0, 0.0, 0.0};
    const math::Vector3 sc_vel{0.0, 7612.6, 0.0};
    const orbit::CartesianState init_state{sc_pos, sc_vel};

    estimation::TranslationalCovariance cov;
    for (std::size_t i = 0; i < 6; ++i) {
        cov(i, i) = (i < 3) ? 100.0 : 1.0;
    }
    estimation::TranslationalEkf filter(config, init_state, cov);

    const math::Vector3 target_pos{6878137.0 + 500.0, 0.0, 0.0}; // 500 m along +X
    const double measured_range = 495.0; // 5 m shorter -> chaser should move toward target (+X)

    const bool updated = filter.update_range(measured_range, target_pos, 2.0);
    CHECK(updated);

    const auto post_state = filter.estimated_state();
    // Position X should have shifted positively towards target
    CHECK(post_state.position.x() > sc_pos.x());

    const auto diag = filter.last_range_diagnostics();
    CHECK_THAT(diag.predicted_range_m, WithinAbs(500.0, 1.0e-12));
    CHECK_THAT(diag.innovation_m, WithinAbs(-5.0, 1.0e-12));
    CHECK(diag.normalized_innovation_squared > 0.0);
}

TEST_CASE("Range measurement Monte Carlo study maintains NIS chi-square(1) consistency", "[estimation][range][monte_carlo]") {
    const std::size_t N_seeds = 300;
    double sum_nis = 0.0;

    const math::Vector3 target_pos{7000.0e3, 0.0, 0.0};
    const math::Vector3 nominal_sc_pos{7000.0e3 - 1000.0, 0.0, 0.0}; // 1000 m nominal range along -X

    const double sigma_r = 10.0;    // 10 m prior position uncertainty
    const double sigma_range = 2.0; // 2 m range sensor noise

    for (std::uint64_t seed = 1; seed <= N_seeds; ++seed) {
        sensors::DeterministicRng rng(seed);

        // Perturb truth position by prior uncertainty
        const double pos_err = rng.gaussian(0.0, sigma_r);
        const math::Vector3 truth_sc_pos = nominal_sc_pos + math::Vector3{pos_err, 0.0, 0.0};

        // Generate range measurement from truth position + noise
        const double true_range = (target_pos - truth_sc_pos).norm();
        const double noise = rng.gaussian(0.0, sigma_range);
        const double meas_range = true_range + noise;

        math::ColVector<6> prior_state{};
        prior_state(0, 0) = nominal_sc_pos.x();
        prior_state(1, 0) = nominal_sc_pos.y();
        prior_state(2, 0) = nominal_sc_pos.z();

        math::Matrix<6, 6> prior_cov;
        for (std::size_t i = 0; i < 6; ++i) {
            prior_cov(i, i) = (i < 3) ? sigma_r * sigma_r : 1.0;
        }

        estimation::RangeUpdateDiagnostics diag;
        const auto update_res = estimation::execute_range_update(
            prior_state,
            prior_cov,
            meas_range,
            target_pos,
            sigma_range,
            &diag
        );
        (void)update_res;

        sum_nis += diag.normalized_innovation_squared;
    }

    const double mean_nis = sum_nis / static_cast<double>(N_seeds);

    // NIS for df=1 has theoretical mean 1.0. Across 300 seeds, should be in [0.75, 1.25]
    CHECK_THAT(mean_nis, WithinAbs(1.0, 0.25));
}
