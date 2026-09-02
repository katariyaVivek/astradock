#include "attitude/principal_inertia.hpp"
#include "attitude/rotational_state.hpp"
#include "estimation/attitude_ekf.hpp"
#include "estimation/diagnostics.hpp"
#include "math/constants.hpp"
#include "sensors/sensor_common.hpp"
#include "sensors/star_tracker.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using namespace astradock;

TEST_CASE("Attitude EKF Canonical Case A: Perfect star tracker measurement produces zero residual", "[estimation][attitude][canonical]") {
    const math::Quaternion q_nominal = math::Quaternion::from_axis_angle(
        math::Vector3{1.0, 2.0, 3.0}.normalized(), 0.5);

    const math::Vector3 residual = estimation::star_tracker_measurement_residual(
        q_nominal, q_nominal);

    CHECK_THAT(residual.norm(), WithinAbs(0.0, 1.0e-15));
}

TEST_CASE("Attitude EKF Canonical Case B: Known small rotation matches expected physical residual", "[estimation][attitude][canonical]") {
    // Nominal attitude
    const math::Quaternion q_nominal = math::Quaternion::from_axis_angle(
        math::Vector3{0.0, 1.0, 0.0}, 0.7);

    // Apply known small body-frame rotation: delta_theta = [0.003, -0.004, 0.005] rad
    const math::Vector3 true_delta_theta{0.003, -0.004, 0.005};
    const double theta_norm = true_delta_theta.norm();
    const math::Quaternion delta_q = math::Quaternion::from_axis_angle(
        true_delta_theta / theta_norm, theta_norm);

    // q_true = q_nom ⊗ delta_q
    const math::Quaternion q_true = q_nominal * delta_q;

    const math::Vector3 residual = estimation::star_tracker_measurement_residual(
        q_nominal, q_true);

    // Residual should match true_delta_theta to second-order accuracy (< 1e-6 relative error)
    CHECK_THAT(residual.x(), WithinAbs(true_delta_theta.x(), 1.0e-6));
    CHECK_THAT(residual.y(), WithinAbs(true_delta_theta.y(), 1.0e-6));
    CHECK_THAT(residual.z(), WithinAbs(true_delta_theta.z(), 1.0e-6));
}

TEST_CASE("Attitude EKF Canonical Case C: Double-cover q and -q produce identical physical residuals", "[estimation][attitude][canonical]") {
    const math::Quaternion q_nominal = math::Quaternion::from_axis_angle(
        math::Vector3{1.0, 1.0, 1.0}.normalized(), 1.2);

    const math::Vector3 delta_theta{0.002, 0.001, -0.003};
    const double theta_norm = delta_theta.norm();
    const math::Quaternion delta_q = math::Quaternion::from_axis_angle(
        delta_theta / theta_norm, theta_norm);
    const math::Quaternion q_meas_pos = (q_nominal * delta_q).normalized();
    const math::Quaternion q_meas_neg(-q_meas_pos.w(), -q_meas_pos.x(), -q_meas_pos.y(), -q_meas_pos.z());

    const math::Vector3 res_pos = estimation::star_tracker_measurement_residual(q_nominal, q_meas_pos);
    const math::Vector3 res_neg = estimation::star_tracker_measurement_residual(q_nominal, q_meas_neg);

    CHECK_THAT(res_pos.x(), WithinAbs(res_neg.x(), 1.0e-15));
    CHECK_THAT(res_pos.y(), WithinAbs(res_neg.y(), 1.0e-15));
    CHECK_THAT(res_pos.z(), WithinAbs(res_neg.z(), 1.0e-15));

    // When measurement is -q_nominal, residual MUST be virtually zero
    const math::Quaternion neg_nominal(-q_nominal.w(), -q_nominal.x(), -q_nominal.y(), -q_nominal.z());
    const math::Vector3 res_identity = estimation::star_tracker_measurement_residual(q_nominal, neg_nominal);
    CHECK_THAT(res_identity.norm(), WithinAbs(0.0, 1.0e-15));
}

TEST_CASE("Attitude EKF Canonical Case D: Constant gyro bias is accurately estimated from star tracker fixes", "[estimation][attitude][bias]") {
    estimation::AttitudeEkfConfig config;
    config.gyro_noise_std_rad_s = 1.0e-4;
    config.star_tracker_noise_std_rad = 1.0e-4;
    config.initial_attitude_error_std_rad = 0.01;
    config.initial_gyro_bias_std_rad_s = 0.01; // 10 mrad/s initial 1-sigma

    estimation::AttitudeEkf filter(config);

    const math::Vector3 true_bias{0.005, -0.003, 0.002}; // 5, -3, 2 mrad/s
    const math::Vector3 true_omega{0.02, 0.01, -0.015};   // rad/s

    math::Quaternion truth_q = math::Quaternion::identity();
    const double dt_gyro = 0.01;   // 100 Hz
    const double dt_st = 0.1;      // 10 Hz
    const double t_total = 40.0;   // 40 seconds

    sensors::DeterministicRng rng(42);

    for (double t = 0.0; t < t_total; t += dt_gyro) {
        // True angular rate propagation
        const double angle = true_omega.norm() * dt_gyro;
        const math::Quaternion dq = math::Quaternion::from_axis_angle(true_omega.normalized(), angle);
        truth_q = (truth_q * dq).normalized();

        // Noisy, biased gyro measurement
        const math::Vector3 gyro_noise{
            rng.gaussian(0.0, config.gyro_noise_std_rad_s),
            rng.gaussian(0.0, config.gyro_noise_std_rad_s),
            rng.gaussian(0.0, config.gyro_noise_std_rad_s)
        };
        const math::Vector3 gyro_meas = true_omega + true_bias + gyro_noise;

        filter.predict(gyro_meas, dt_gyro);

        // Star tracker measurement update at 10 Hz
        if (std::fmod(t + 1.0e-6, dt_st) < dt_gyro) {
            const double st_err_angle = rng.gaussian(0.0, config.star_tracker_noise_std_rad);
            const math::Vector3 st_axis = rng.uniform_unit_vector3();
            const math::Quaternion st_noise_q = math::Quaternion::from_axis_angle(st_axis, st_err_angle);
            const math::Quaternion meas_q = (truth_q * st_noise_q).normalized();

            sensors::StarTrackerMeasurement meas{t, meas_q, true};
            filter.update_star_tracker(meas);
        }
    }

    // Filter should estimate gyro bias to within ~0.2 mrad/s (2e-4 rad/s)
    const math::Vector3 est_bias = filter.estimate().gyro_bias_rad_s;
    CHECK_THAT(est_bias.x(), WithinAbs(true_bias.x(), 2.0e-4));
    CHECK_THAT(est_bias.y(), WithinAbs(true_bias.y(), 2.0e-4));
    CHECK_THAT(est_bias.z(), WithinAbs(true_bias.z(), 2.0e-4));

    // Attitude error should be tightly bounded (< 0.5 mrad)
    const math::Vector3 att_err = estimation::attitude_error_vector_rad(
        truth_q, filter.estimate().nominal_orientation);
    CHECK_THAT(att_err.norm(), WithinAbs(0.0, 5.0e-4));
}

TEST_CASE("Attitude EKF Canonical Case E: Known angular-rate propagation matches analytical rotation", "[estimation][attitude][propagation]") {
    const math::Vector3 omega_axis{0.0, 0.0, 1.0}; // pure Z spin
    const double spin_rate = 0.5; // rad/s
    const math::Vector3 omega = omega_axis * spin_rate;

    const math::Quaternion q0 = math::Quaternion::identity();
    const double dt = 0.1;
    math::Quaternion q_num = q0;

    for (int step = 0; step < 100; ++step) {
        q_num = estimation::propagate_nominal_quaternion(q_num, omega, dt);
    }

    const double total_time = 100 * dt; // 10 s
    const double expected_angle = spin_rate * total_time; // 5.0 rad
    const math::Quaternion q_analytical = math::Quaternion::from_axis_angle(omega_axis, expected_angle);

    const math::Vector3 diff = estimation::star_tracker_measurement_residual(q_analytical, q_num);
    CHECK_THAT(diff.norm(), WithinAbs(0.0, 1.0e-12));
}

TEST_CASE("Analytical attitude error dynamics Jacobian F survives finite-difference audit across rates and attitudes", "[estimation][attitude][jacobian]") {
    const std::vector<math::Vector3> test_rates = {
        {0.01, 0.02, -0.03},
        {0.1, -0.05, 0.2},
        {-0.3, 0.2, 0.1},
        {0.0, 0.0, 0.05}
    };

    const double eps = 1.0e-6; // finite difference perturbation

    for (const auto& omega : test_rates) {
        const auto F_analytical = estimation::attitude_error_dynamics_matrix(omega);

        // Audit attitude error block: d(delta_dot_theta) / d(delta_theta) = - [omega]x
        // We evaluate rate of change of error delta_dot_theta = -omega x delta_theta
        for (std::size_t j = 0; j < 3; ++j) {
            const math::Vector3 pert_pos{
                (j == 0 ? eps : 0.0),
                (j == 1 ? eps : 0.0),
                (j == 2 ? eps : 0.0)
            };
            const math::Vector3 pert_neg{
                (j == 0 ? -eps : 0.0),
                (j == 1 ? -eps : 0.0),
                (j == 2 ? -eps : 0.0)
            };

            const math::Vector3 f_pos = (omega.cross(pert_pos)) * (-1.0);
            const math::Vector3 f_neg = (omega.cross(pert_neg)) * (-1.0);
            const math::Vector3 fd_deriv = (f_pos - f_neg) / (2.0 * eps);

            CHECK_THAT(F_analytical(0, j), WithinAbs(fd_deriv.x(), 1.0e-9));
            CHECK_THAT(F_analytical(1, j), WithinAbs(fd_deriv.y(), 1.0e-9));
            CHECK_THAT(F_analytical(2, j), WithinAbs(fd_deriv.z(), 1.0e-9));
        }

        // Audit bias block: d(delta_dot_theta) / d(delta_b_g) = -I_3
        for (std::size_t j = 0; j < 3; ++j) {
            CHECK_THAT(F_analytical(j, j + 3), WithinAbs(-1.0, 1.0e-15));
        }
    }
}

TEST_CASE("Covariance reset transformation matches analytical first-order coordinate shift", "[estimation][attitude][reset]") {
    estimation::AttitudeCovariance P;
    for (std::size_t i = 0; i < 6; ++i) {
        P(i, i) = 0.01 * (static_cast<double>(i) + 1.0);
    }
    // Add small off-diagonal correlation
    P(0, 1) = 0.001; P(1, 0) = 0.001;
    P(1, 2) = 0.0015; P(2, 1) = 0.0015;

    const math::Vector3 delta_theta{0.02, -0.015, 0.01}; // 20 mrad correction

    const estimation::AttitudeCovariance P_reset =
        estimation::attitude_covariance_reset(P, delta_theta);

    // Verify J_reset = I - 0.5 * [delta_theta]x
    // To first order, the trace of the attitude block should remain invariant because [delta_theta]x is traceless and skew-symmetric
    const double trace_before = P(0, 0) + P(1, 1) + P(2, 2);
    const double trace_after = P_reset(0, 0) + P_reset(1, 1) + P_reset(2, 2);
    CHECK_THAT(trace_after, WithinAbs(trace_before, 1.0e-5));

    // Verify symmetry of reset covariance
    CHECK(P_reset.symmetry_error() < 1.0e-12);

    // Zero correction must produce identically equal covariance
    const estimation::AttitudeCovariance P_zero_reset =
        estimation::attitude_covariance_reset(P, math::Vector3{0.0, 0.0, 0.0});
    for (std::size_t r = 0; r < 6; ++r) {
        for (std::size_t c = 0; c < 6; ++c) {
            CHECK_THAT(P_zero_reset(r, c), WithinAbs(P(r, c), 1.0e-15));
        }
    }
}

TEST_CASE("Star tracker dropout causes attitude covariance growth and recovery restores confidence", "[estimation][attitude][dropout]") {
    estimation::AttitudeEkfConfig config;
    config.gyro_noise_std_rad_s = 1.0e-3;
    config.star_tracker_noise_std_rad = 1.0e-3;
    estimation::AttitudeEkf filter(config);

    const math::Vector3 omega{0.01, -0.02, 0.015};
    const double dt = 0.1;

    // Track for 5 seconds with valid star tracker updates
    for (int i = 0; i < 50; ++i) {
        filter.predict(omega, dt);
        sensors::StarTrackerMeasurement meas{static_cast<double>(i) * dt, math::Quaternion::identity(), true};
        filter.update_star_tracker(meas);
    }
    const double sigma_attitude_tracking = std::sqrt(filter.covariance()(0, 0));

    // Dropout for 10 seconds (no star tracker updates)
    for (int i = 50; i < 150; ++i) {
        filter.predict(omega, dt);
        // Invalid measurement during dropout
        sensors::StarTrackerMeasurement meas_dropout{static_cast<double>(i) * dt, math::Quaternion::identity(), false};
        const bool updated = filter.update_star_tracker(meas_dropout);
        CHECK_FALSE(updated);
    }
    const double sigma_attitude_dropout = std::sqrt(filter.covariance()(0, 0));

    // Uncertainty MUST have grown substantially during dropout
    CHECK(sigma_attitude_dropout > sigma_attitude_tracking * 3.0);

    // Star tracker measurement returns
    sensors::StarTrackerMeasurement meas_recovery{15.0, math::Quaternion::identity(), true};
    const bool recovered = filter.update_star_tracker(meas_recovery);
    CHECK(recovered);
    const double sigma_attitude_recovered = std::sqrt(filter.covariance()(0, 0));

    // Uncertainty contracts sharply on measurement return
    CHECK(sigma_attitude_recovered < sigma_attitude_dropout);
}

TEST_CASE("Attitude EKF Monte Carlo validation maintains NEES and NIS consistency across 100 seeds", "[estimation][attitude][monte_carlo]") {
    const std::size_t N_seeds = 100;
    double sum_nis = 0.0;
    double sum_nees = 0.0;
    std::size_t total_updates = 0;

    for (std::uint64_t seed = 1; seed <= N_seeds; ++seed) {
        sensors::DeterministicRng rng(seed);

        estimation::AttitudeEkfConfig config;
        config.gyro_noise_std_rad_s = 5.0e-4;
        config.star_tracker_noise_std_rad = 5.0e-4;
        config.initial_attitude_error_std_rad = 0.005;

        // Randomize initial attitude
        const math::Vector3 init_err{
            rng.gaussian(0.0, config.initial_attitude_error_std_rad),
            rng.gaussian(0.0, config.initial_attitude_error_std_rad),
            rng.gaussian(0.0, config.initial_attitude_error_std_rad)
        };
        const double init_err_norm = init_err.norm();
        const math::Quaternion q_init_err = (init_err_norm > 0.0)
            ? math::Quaternion::from_axis_angle(init_err / init_err_norm, init_err_norm)
            : math::Quaternion::identity();

        estimation::AttitudeEstimate est;
        est.nominal_orientation = q_init_err;
        estimation::AttitudeEkf filter(config, est);

        math::Quaternion truth_q = math::Quaternion::identity();
        const math::Vector3 true_rate{0.01, -0.01, 0.02};

        // Run 5 seconds (50 steps of 0.1 s)
        for (int step = 0; step < 50; ++step) {
            const double dt = 0.1;
            const double angle = true_rate.norm() * dt;
            const math::Quaternion dq = math::Quaternion::from_axis_angle(true_rate.normalized(), angle);
            truth_q = (truth_q * dq).normalized();

            const double discrete_gyro_sigma = config.gyro_noise_std_rad_s / std::sqrt(dt);
            const math::Vector3 gyro_noise{
                rng.gaussian(0.0, discrete_gyro_sigma),
                rng.gaussian(0.0, discrete_gyro_sigma),
                rng.gaussian(0.0, discrete_gyro_sigma)
            };
            filter.predict(true_rate + gyro_noise, dt);

            // Update at each step
            const double st_err = rng.gaussian(0.0, config.star_tracker_noise_std_rad);
            const math::Vector3 st_axis = rng.uniform_unit_vector3();
            const math::Quaternion st_noise_q = math::Quaternion::from_axis_angle(st_axis, st_err);
            const math::Quaternion meas_q = (truth_q * st_noise_q).normalized();

            sensors::StarTrackerMeasurement meas{step * dt, meas_q, true};
            filter.update_star_tracker(meas);

            if (step > 20) { // steady state
                sum_nis += filter.last_star_tracker_diagnostics().normalized_innovation_squared;
                math::Matrix<3, 3> p_att;
                for (std::size_t r = 0; r < 3; ++r) {
                    for (std::size_t c = 0; c < 3; ++c) {
                        p_att(r, c) = filter.covariance()(r, c);
                    }
                }
                sum_nees += estimation::attitude_nees(truth_q, filter.estimate().nominal_orientation, p_att);
                total_updates++;
            }
        }
    }

    const double mean_nis = sum_nis / static_cast<double>(total_updates);
    const double mean_nees = sum_nees / static_cast<double>(total_updates);

    // NIS for df=3 has theoretical mean 3.0. Across 100 seeds, should be in [2.5, 3.5]
    CHECK_THAT(mean_nis, WithinAbs(3.0, 0.5));

    // NEES for df=3 has theoretical mean 3.0. Across 100 seeds, should be in [2.0, 4.0]
    CHECK_THAT(mean_nees, WithinAbs(3.0, 1.0));
}

TEST_CASE("Attitude EKF obeys truth non-interference contract", "[estimation][attitude][truth_regression]") {
    const math::Quaternion initial_q = math::Quaternion::from_axis_angle(
        math::Vector3{0.0, 1.0, 0.0}, 0.5);

    estimation::AttitudeEkfConfig config;
    estimation::AttitudeEkf filter(config);

    // Pass a copy of a measurement; verify original values are completely unchanged
    sensors::StarTrackerMeasurement meas{1.0, initial_q, true};
    const math::Quaternion meas_q_copy = meas.orientation_eci_from_body;

    filter.update_star_tracker(meas);

    CHECK(meas.orientation_eci_from_body.w() == meas_q_copy.w());
    CHECK(meas.orientation_eci_from_body.x() == meas_q_copy.x());
    CHECK(meas.orientation_eci_from_body.y() == meas_q_copy.y());
    CHECK(meas.orientation_eci_from_body.z() == meas_q_copy.z());
}
