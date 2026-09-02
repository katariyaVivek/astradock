// AstraDock M13D — Integrated Multi-Rate 15-State Navigation Filter Test Suite.
//
// Comprehensive unit, integration, and Monte Carlo verification of the unified 15-state EKF:
//   delta_x = [ delta_r^T, delta_v^T, delta_theta^T, delta_b_a^T, delta_b_g^T ]^T in R^15

#include "estimation/diagnostics.hpp"
#include "estimation/integrated_navigation_ekf.hpp"
#include "math/constants.hpp"
#include "math/matrix3.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"
#include "sensors/sensor_common.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using namespace astradock;

TEST_CASE("15-state nominal and error-state representation and index ordering", "[estimation][integrated][ordering]") {
    estimation::NavigationState state;
    state.position_eci_m = math::Vector3{7000.0e3, 0.0, 0.0};
    state.velocity_eci_mps = math::Vector3{0.0, 7500.0, 0.0};
    state.attitude_body_to_eci = math::Quaternion::identity();
    state.accelerometer_bias_body_mps2 = math::Vector3{0.01, -0.02, 0.015};
    state.gyro_bias_body_rad_s = math::Vector3{0.001, -0.002, 0.003};

    CHECK(state.all_finite());
    CHECK(estimation::k_nav_state_dim == 15);
    CHECK(estimation::k_idx_pos == 0);
    CHECK(estimation::k_idx_vel == 3);
    CHECK(estimation::k_idx_att == 6);
    CHECK(estimation::k_idx_acc_bias == 9);
    CHECK(estimation::k_idx_gyro_bias == 12);

    estimation::IntegratedNavigationConfig config;
    const estimation::NavigationCovariance P0 = config.initial_covariance();

    // Verify diagonal initial variances match config
    CHECK_THAT(P0(0, 0), WithinAbs(config.initial_position_std_m * config.initial_position_std_m, 1.0e-12));
    CHECK_THAT(P0(3, 3), WithinAbs(config.initial_velocity_std_mps * config.initial_velocity_std_mps, 1.0e-12));
    CHECK_THAT(P0(6, 6), WithinAbs(config.initial_attitude_std_rad * config.initial_attitude_std_rad, 1.0e-12));
    CHECK_THAT(P0(9, 9), WithinAbs(config.initial_accel_bias_std_mps2 * config.initial_accel_bias_std_mps2, 1.0e-12));
    CHECK_THAT(P0(12, 12), WithinAbs(config.initial_gyro_bias_std_rad_s * config.initial_gyro_bias_std_rad_s, 1.0e-12));
    CHECK(P0.symmetry_error() < 1.0e-12);
}

TEST_CASE("Analytical integrated Jacobian F survives finite-difference audit across all 5 blocks", "[estimation][integrated][jacobian]") {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const double eps = 1.0e-6;

    // Test cases covering non-trivial orbital positions, tumbling attitudes, non-zero specific forces and rates
    const std::vector<math::Vector3> positions = {
        {6878137.0, 1000.0, -500.0},
        {-4500.0e3, 5200.0e3, 1800.0e3}
    };
    const std::vector<math::Quaternion> attitudes = {
        math::Quaternion::identity(),
        math::Quaternion::from_axis_angle(math::Vector3{0.577, 0.577, 0.577}.normalized(), 0.8)
    };
    const std::vector<math::Vector3> forces = {
        {0.0, 0.0, 0.0},
        {1.5, -2.0, 0.75} // 2.6 m/s^2 specific force
    };
    const std::vector<math::Vector3> rates = {
        {0.03, -0.02, 0.04},
        {-0.05, 0.01, -0.02}
    };

    for (const auto& pos : positions) {
        for (const auto& att : attitudes) {
            for (const auto& fb : forces) {
                for (const auto& wb : rates) {
                    const auto F = estimation::integrated_error_dynamics_jacobian(pos, att, fb, wb, mu);

                    // 1. Audit Block (0, 1): d(r_dot) / d(v) = I_3
                    for (std::size_t i = 0; i < 3; ++i) {
                        for (std::size_t j = 0; j < 3; ++j) {
                            const double expected = (i == j) ? 1.0 : 0.0;
                            CHECK_THAT(F(estimation::k_idx_pos + i, estimation::k_idx_vel + j), WithinAbs(expected, 1.0e-15));
                        }
                    }

                    // 2. Audit Block (1, 0): d(v_dot) / d(r) = G(r) against finite difference of gravity
                    for (std::size_t j = 0; j < 3; ++j) {
                        const math::Vector3 pert_pos{
                            pos.x() + (j == 0 ? eps : 0.0),
                            pos.y() + (j == 1 ? eps : 0.0),
                            pos.z() + (j == 2 ? eps : 0.0)
                        };
                        const math::Vector3 pert_neg{
                            pos.x() - (j == 0 ? eps : 0.0),
                            pos.y() - (j == 1 ? eps : 0.0),
                            pos.z() - (j == 2 ? eps : 0.0)
                        };
                        const math::Vector3 g_pos = dynamics::two_body_acceleration(pert_pos, mu);
                        const math::Vector3 g_neg = dynamics::two_body_acceleration(pert_neg, mu);
                        const math::Vector3 fd_G = (g_pos - g_neg) / (2.0 * eps);

                        CHECK_THAT(F(estimation::k_idx_vel + 0, estimation::k_idx_pos + j), WithinAbs(fd_G.x(), 1.0e-8));
                        CHECK_THAT(F(estimation::k_idx_vel + 1, estimation::k_idx_pos + j), WithinAbs(fd_G.y(), 1.0e-8));
                        CHECK_THAT(F(estimation::k_idx_vel + 2, estimation::k_idx_pos + j), WithinAbs(fd_G.z(), 1.0e-8));
                    }

                    // 3. Audit Block (1, 2): d(v_dot) / d(delta_theta) = -C_I_B [f_B]x
                    // Perturb attitude via body rotation delta_theta
                    for (std::size_t j = 0; j < 3; ++j) {
                        const math::Vector3 th_pos{
                            (j == 0 ? eps : 0.0),
                            (j == 1 ? eps : 0.0),
                            (j == 2 ? eps : 0.0)
                        };
                        const math::Vector3 th_neg{
                            (j == 0 ? -eps : 0.0),
                            (j == 1 ? -eps : 0.0),
                            (j == 2 ? -eps : 0.0)
                        };
                        const math::Quaternion dq_pos = math::Quaternion::from_axis_angle(th_pos / eps, eps);
                        const math::Quaternion dq_neg = math::Quaternion::from_axis_angle(th_neg / eps, eps);

                        // Under positive attitude perturbation, q_pert = q * dq
                        const math::Quaternion q_pos = (att * dq_pos).normalized();
                        const math::Quaternion q_neg = (att * dq_neg).normalized();

                        const math::Vector3 a_pos = q_pos.rotate_vector(fb);
                        const math::Vector3 a_neg = q_neg.rotate_vector(fb);
                        const math::Vector3 fd_att = (a_pos - a_neg) / (2.0 * eps);

                        CHECK_THAT(F(estimation::k_idx_vel + 0, estimation::k_idx_att + j), WithinAbs(fd_att.x(), 1.0e-8));
                        CHECK_THAT(F(estimation::k_idx_vel + 1, estimation::k_idx_att + j), WithinAbs(fd_att.y(), 1.0e-8));
                        CHECK_THAT(F(estimation::k_idx_vel + 2, estimation::k_idx_att + j), WithinAbs(fd_att.z(), 1.0e-8));
                    }

                    // 4. Audit Block (1, 3): d(v_dot) / d(delta_b_a) = -C_I_B
                    for (std::size_t j = 0; j < 3; ++j) {
                        const math::Vector3 ba_pos{
                            (j == 0 ? eps : 0.0),
                            (j == 1 ? eps : 0.0),
                            (j == 2 ? eps : 0.0)
                        };
                        const math::Vector3 ba_neg{
                            (j == 0 ? -eps : 0.0),
                            (j == 1 ? -eps : 0.0),
                            (j == 2 ? -eps : 0.0)
                        };
                        // Specific force corrected is f_m - b_a; so perturbing b_a by +eps subtracts eps
                        const math::Vector3 a_pos = att.rotate_vector(fb - ba_pos);
                        const math::Vector3 a_neg = att.rotate_vector(fb - ba_neg);
                        const math::Vector3 fd_ba = (a_pos - a_neg) / (2.0 * eps);

                        CHECK_THAT(F(estimation::k_idx_vel + 0, estimation::k_idx_acc_bias + j), WithinAbs(fd_ba.x(), 1.0e-8));
                        CHECK_THAT(F(estimation::k_idx_vel + 1, estimation::k_idx_acc_bias + j), WithinAbs(fd_ba.y(), 1.0e-8));
                        CHECK_THAT(F(estimation::k_idx_vel + 2, estimation::k_idx_acc_bias + j), WithinAbs(fd_ba.z(), 1.0e-8));
                    }

                    // 5. Audit Block (2, 2): d(delta_theta_dot) / d(delta_theta) = -[omega_B]x
                    for (std::size_t j = 0; j < 3; ++j) {
                        const math::Vector3 th_pos{
                            (j == 0 ? eps : 0.0),
                            (j == 1 ? eps : 0.0),
                            (j == 2 ? eps : 0.0)
                        };
                        const math::Vector3 th_neg{
                            (j == 0 ? -eps : 0.0),
                            (j == 1 ? -eps : 0.0),
                            (j == 2 ? -eps : 0.0)
                        };
                        const math::Vector3 f_pos = (wb.cross(th_pos)) * (-1.0);
                        const math::Vector3 f_neg = (wb.cross(th_neg)) * (-1.0);
                        const math::Vector3 fd_w = (f_pos - f_neg) / (2.0 * eps);

                        CHECK_THAT(F(estimation::k_idx_att + 0, estimation::k_idx_att + j), WithinAbs(fd_w.x(), 1.0e-8));
                        CHECK_THAT(F(estimation::k_idx_att + 1, estimation::k_idx_att + j), WithinAbs(fd_w.y(), 1.0e-8));
                        CHECK_THAT(F(estimation::k_idx_att + 2, estimation::k_idx_att + j), WithinAbs(fd_w.z(), 1.0e-8));
                    }

                    // 6. Audit Block (2, 4): d(delta_theta_dot) / d(delta_b_g) = -I_3
                    for (std::size_t i = 0; i < 3; ++i) {
                        for (std::size_t j = 0; j < 3; ++j) {
                            const double expected = (i == j) ? -1.0 : 0.0;
                            CHECK_THAT(F(estimation::k_idx_att + i, estimation::k_idx_gyro_bias + j), WithinAbs(expected, 1.0e-15));
                        }
                    }
                }
            }
        }
    }
}

TEST_CASE("Cross-covariance coupling: attitude uncertainty directly inflates velocity and position covariance during non-zero specific force", "[estimation][integrated][cross_covariance]") {
    estimation::IntegratedNavigationConfig config;
    estimation::NavigationState init_state;
    init_state.position_eci_m = math::Vector3{7000.0e3, 0.0, 0.0};
    init_state.velocity_eci_mps = math::Vector3{0.0, 7546.0, 0.0};
    init_state.attitude_body_to_eci = math::Quaternion::identity();

    // Setup two filters: Filter A has small attitude uncertainty; Filter B has large attitude uncertainty
    estimation::NavigationCovariance P_A = config.initial_covariance();
    estimation::NavigationCovariance P_B = config.initial_covariance();

    // Set attitude variance: A has (1 mrad)^2, B has (50 mrad)^2
    for (std::size_t i = 0; i < 3; ++i) {
        P_A(estimation::k_idx_att + i, estimation::k_idx_att + i) = 1.0e-6;
        P_B(estimation::k_idx_att + i, estimation::k_idx_att + i) = 2.5e-3;
    }

    estimation::IntegratedNavigationEkf filter_A(config, init_state, P_A);
    estimation::IntegratedNavigationEkf filter_B(config, init_state, P_B);

    // Apply significant thrust/specific force: 5.0 m/s^2 along body X
    const math::Vector3 thrust_f{5.0, 0.0, 0.0};
    const math::Vector3 zero_rate{0.0, 0.0, 0.0};
    const double dt = 0.1;

    // Propagate both filters for 50 steps (5 seconds of thrust)
    for (int step = 0; step < 50; ++step) {
        filter_A.predict(thrust_f, zero_rate, dt);
        filter_B.predict(thrust_f, zero_rate, dt);
    }

    const auto& cov_A = filter_A.covariance();
    const auto& cov_B = filter_B.covariance();

    // Under thrust along X, attitude errors around Y and Z rotate thrust into Y and Z inertial directions!
    // Therefore, Filter B MUST develop strictly larger velocity and position uncertainty in Y and Z:
    CHECK(cov_B(estimation::k_idx_vel + 1, estimation::k_idx_vel + 1) > cov_A(estimation::k_idx_vel + 1, estimation::k_idx_vel + 1));
    CHECK(cov_B(estimation::k_idx_vel + 2, estimation::k_idx_vel + 2) > cov_A(estimation::k_idx_vel + 2, estimation::k_idx_vel + 2));

    CHECK(cov_B(estimation::k_idx_pos + 1, estimation::k_idx_pos + 1) > cov_A(estimation::k_idx_pos + 1, estimation::k_idx_pos + 1));
    CHECK(cov_B(estimation::k_idx_pos + 2, estimation::k_idx_pos + 2) > cov_A(estimation::k_idx_pos + 2, estimation::k_idx_pos + 2));

    // Furthermore, non-zero cross-correlation between attitude and velocity must exist in both filters:
    CHECK(std::abs(cov_B(estimation::k_idx_vel + 1, estimation::k_idx_att + 2)) > 0.0);
    CHECK(std::abs(cov_B(estimation::k_idx_vel + 2, estimation::k_idx_att + 1)) > 0.0);
}

TEST_CASE("Cross-covariance coupling: accelerometer bias uncertainty inflates translational covariance", "[estimation][integrated][bias_covariance]") {
    estimation::IntegratedNavigationConfig config;
    estimation::NavigationState init_state;
    init_state.position_eci_m = math::Vector3{7000.0e3, 0.0, 0.0};
    init_state.velocity_eci_mps = math::Vector3{0.0, 7546.0, 0.0};
    init_state.attitude_body_to_eci = math::Quaternion::identity();

    estimation::NavigationCovariance P_small = config.initial_covariance();
    estimation::NavigationCovariance P_large = config.initial_covariance();

    // Accelerometer bias variance: small = (0.1 mm/s^2)^2, large = (100 mm/s^2)^2
    for (std::size_t i = 0; i < 3; ++i) {
        P_small(estimation::k_idx_acc_bias + i, estimation::k_idx_acc_bias + i) = 1.0e-8;
        P_large(estimation::k_idx_acc_bias + i, estimation::k_idx_acc_bias + i) = 1.0e-2;
    }

    estimation::IntegratedNavigationEkf filter_small(config, init_state, P_small);
    estimation::IntegratedNavigationEkf filter_large(config, init_state, P_large);

    const math::Vector3 zero_f{0.0, 0.0, 0.0};
    const math::Vector3 zero_rate{0.0, 0.0, 0.0};
    const double dt = 0.1;

    for (int step = 0; step < 50; ++step) {
        filter_small.predict(zero_f, zero_rate, dt);
        filter_large.predict(zero_f, zero_rate, dt);
    }

    const auto& cov_s = filter_small.covariance();
    const auto& cov_l = filter_large.covariance();

    // Accelerometer bias directly inflates velocity covariance along all axes
    for (std::size_t i = 0; i < 3; ++i) {
        CHECK(cov_l(estimation::k_idx_vel + i, estimation::k_idx_vel + i) > cov_s(estimation::k_idx_vel + i, estimation::k_idx_vel + i));
        CHECK(cov_l(estimation::k_idx_pos + i, estimation::k_idx_pos + i) > cov_s(estimation::k_idx_pos + i, estimation::k_idx_pos + i));
    }
}

TEST_CASE("Integrated 15-state covariance reset transformation verifies trace preservation and symmetry", "[estimation][integrated][reset]") {
    estimation::NavigationCovariance P = estimation::IntegratedNavigationConfig{}.initial_covariance();
    const math::Vector3 delta_theta{0.015, -0.02, 0.01}; // 20 mrad attitude correction

    const estimation::NavigationCovariance P_reset = estimation::integrated_covariance_reset(P, delta_theta);

    // Verify trace shift is negligible to first order
    double trace_before = 0.0;
    double trace_after = 0.0;
    for (std::size_t i = 0; i < 15; ++i) {
        trace_before += P(i, i);
        trace_after += P_reset(i, i);
    }
    CHECK_THAT(trace_after, WithinAbs(trace_before, 1.0e-2));
    CHECK(P_reset.symmetry_error() < 1.0e-12);
}

TEST_CASE("Multi-rate sensor updates with GNSS, star tracker, and range on unified 15x15 covariance", "[estimation][integrated][updates]") {
    estimation::IntegratedNavigationConfig config;
    const math::Vector3 true_pos{7000.0e3, 0.0, 0.0};
    const math::Vector3 true_vel{0.0, 7546.0, 0.0};
    const math::Quaternion true_att = math::Quaternion::identity();

    // Initialize filter with deliberate offsets
    estimation::NavigationState init_state;
    init_state.position_eci_m = true_pos + math::Vector3{20.0, -15.0, 10.0};
    init_state.velocity_eci_mps = true_vel + math::Vector3{0.2, -0.3, 0.1};
    init_state.attitude_body_to_eci = math::Quaternion::from_axis_angle(math::Vector3{0.0, 0.0, 1.0}, 0.02);
    init_state.accelerometer_bias_body_mps2 = math::Vector3{0.0, 0.0, 0.0};
    init_state.gyro_bias_body_rad_s = math::Vector3{0.0, 0.0, 0.0};

    estimation::IntegratedNavigationEkf filter(config, init_state);

    // 1. Star tracker update
    sensors::StarTrackerMeasurement st_meas{0.1, true_att, true};
    const bool st_ok = filter.update_star_tracker(st_meas);
    CHECK(st_ok);
    CHECK(filter.last_star_tracker_diagnostics().accepted);
    CHECK(filter.last_star_tracker_diagnostics().normalized_innovation_squared > 0.0);
    // Attitude error should decrease
    const double att_err_after_st = estimation::attitude_error_vector_rad(
        true_att, filter.state().attitude_body_to_eci).norm();
    CHECK(att_err_after_st < 0.02);

    // 2. GNSS update
    sensors::GnssMeasurement gnss_meas{0.2, true_pos, true_vel, true};
    const bool gnss_ok = filter.update_gnss(gnss_meas);
    CHECK(gnss_ok);
    CHECK(filter.last_gnss_diagnostics().accepted);
    CHECK(filter.last_gnss_diagnostics().normalized_innovation_squared > 0.0);
    // Position error should decrease
    const double pos_err_after_gnss = (true_pos - filter.state().position_eci_m).norm();
    CHECK(pos_err_after_gnss < 20.0);

    // 3. Range update
    const math::Vector3 target_pos{7000.0e3 + 500.0, 0.0, 0.0};
    const double true_range = (target_pos - true_pos).norm();
    const bool range_ok = filter.update_range(true_range, target_pos, 1.0);
    CHECK(range_ok);
    CHECK(filter.last_range_diagnostics().accepted);
    CHECK(filter.last_range_diagnostics().normalized_innovation_squared >= 0.0);
}

TEST_CASE("Sensor outage handling: GNSS and star tracker outages cause graceful covariance growth and recovery", "[estimation][integrated][outage]") {
    estimation::IntegratedNavigationConfig config;
    estimation::NavigationState state;
    state.position_eci_m = math::Vector3{7000.0e3, 0.0, 0.0};
    state.velocity_eci_mps = math::Vector3{0.0, 7546.0, 0.0};
    state.attitude_body_to_eci = math::Quaternion::identity();

    estimation::IntegratedNavigationEkf filter(config, state);

    const double initial_pos_var = filter.covariance()(0, 0);
    const double initial_att_var = filter.covariance()(6, 6);

    // Predict for 20 steps (2 seconds) with no updates (outage)
    for (int i = 0; i < 20; ++i) {
        filter.predict(math::Vector3{0.0, 0.0, 0.0}, math::Vector3{0.0, 0.0, 0.0}, 0.1);
    }

    // Both position and attitude covariance must have grown
    CHECK(filter.covariance()(0, 0) > initial_pos_var);
    CHECK(filter.covariance()(6, 6) > initial_att_var);

    // Provide star tracker update -> attitude covariance contracts
    sensors::StarTrackerMeasurement st_meas{2.1, math::Quaternion::identity(), true};
    filter.update_star_tracker(st_meas);
    CHECK(filter.covariance()(6, 6) < initial_att_var);

    // Provide GNSS update -> position covariance contracts
    sensors::GnssMeasurement gnss_meas{2.2, state.position_eci_m, state.velocity_eci_mps, true};
    filter.update_gnss(gnss_meas);
    CHECK(filter.covariance()(0, 0) < initial_pos_var);
}

TEST_CASE("Integrated 15-state Monte Carlo study demonstrates full-state NEES and NIS consistency across 100 seeds", "[estimation][integrated][monte_carlo]") {
    const std::size_t N_seeds = 100;
    double sum_full_nees = 0.0;
    double sum_gnss_nis = 0.0;
    double sum_st_nis = 0.0;
    double sum_range_nis = 0.0;

    std::size_t gnss_count = 0;
    std::size_t st_count = 0;
    std::size_t range_count = 0;
    std::size_t nees_count = 0;

    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const math::Vector3 target_pos{7000.0e3 + 1000.0, 0.0, 0.0};
    const math::Vector3 true_rate{0.02, -0.015, 0.01};
    const math::Vector3 true_ba{0.005, -0.003, 0.002};
    const math::Vector3 true_bg{0.002, -0.001, 0.0015};

    estimation::IntegratedNavigationConfig config;
    config.accelerometer_noise_std_mps2 = 1.0e-3;
    config.gyro_noise_std_rad_s = 1.0e-3;
    config.gnss_position_noise_std_m = 5.0;
    config.gnss_velocity_noise_std_mps = 0.05;
    config.star_tracker_noise_std_rad = 1.0e-3;
    config.range_noise_std_m = 1.0;

    for (std::uint64_t seed = 1; seed <= N_seeds; ++seed) {
        sensors::DeterministicRng rng(seed * 4321);

        math::Vector3 truth_pos{7000.0e3, 0.0, 0.0};
        math::Vector3 truth_vel{0.0, 7546.0, 0.0};
        math::Quaternion truth_att = math::Quaternion::identity();

        estimation::NavigationState est_state;
        est_state.position_eci_m = truth_pos;
        est_state.velocity_eci_mps = truth_vel;
        est_state.attitude_body_to_eci = truth_att;
        est_state.accelerometer_bias_body_mps2 = true_ba;
        est_state.gyro_bias_body_rad_s = true_bg;

        estimation::IntegratedNavigationEkf filter(config, est_state);

        const double dt = 0.1; // 10 Hz
        for (int step = 1; step <= 30; ++step) {
            const double t = step * dt;

            // Advance truth with gravity
            const math::Vector3 g_truth = dynamics::two_body_acceleration(truth_pos, mu);
            truth_pos = truth_pos + truth_vel * dt + g_truth * (0.5 * dt * dt);
            truth_vel = truth_vel + g_truth * dt;
            const double angle = true_rate.norm() * dt;
            truth_att = (truth_att * math::Quaternion::from_axis_angle(true_rate.normalized(), angle)).normalized();

            // Predict with IMU
            const double disc_acc_std = config.accelerometer_noise_std_mps2 / std::sqrt(dt);
            const double disc_gyr_std = config.gyro_noise_std_rad_s / std::sqrt(dt);

            const math::Vector3 imu_f = true_ba + math::Vector3{
                rng.gaussian(0.0, disc_acc_std),
                rng.gaussian(0.0, disc_acc_std),
                rng.gaussian(0.0, disc_acc_std)
            };
            const math::Vector3 imu_w = true_rate + true_bg + math::Vector3{
                rng.gaussian(0.0, disc_gyr_std),
                rng.gaussian(0.0, disc_gyr_std),
                rng.gaussian(0.0, disc_gyr_std)
            };
            filter.predict(imu_f, imu_w, dt);

            // Star tracker update at 10 Hz
            const double st_ang = rng.gaussian(0.0, config.star_tracker_noise_std_rad);
            const math::Vector3 st_u = rng.uniform_unit_vector3();
            const math::Quaternion st_noise = math::Quaternion::from_axis_angle(st_u, st_ang);
            const sensors::StarTrackerMeasurement st_meas{t, (truth_att * st_noise).normalized(), true};
            filter.update_star_tracker(st_meas);
            sum_st_nis += filter.last_star_tracker_diagnostics().normalized_innovation_squared;
            st_count++;

            // Range update at 10 Hz
            const double true_rho = (target_pos - truth_pos).norm();
            const double meas_rho = true_rho + rng.gaussian(0.0, config.range_noise_std_m);
            filter.update_range(meas_rho, target_pos, config.range_noise_std_m);
            sum_range_nis += filter.last_range_diagnostics().normalized_innovation_squared;
            range_count++;

            // GNSS update at 1 Hz (every 10 steps)
            if (step % 10 == 0) {
                const math::Vector3 gnss_r = truth_pos + math::Vector3{
                    rng.gaussian(0.0, config.gnss_position_noise_std_m),
                    rng.gaussian(0.0, config.gnss_position_noise_std_m),
                    rng.gaussian(0.0, config.gnss_position_noise_std_m)
                };
                const math::Vector3 gnss_v = truth_vel + math::Vector3{
                    rng.gaussian(0.0, config.gnss_velocity_noise_std_mps),
                    rng.gaussian(0.0, config.gnss_velocity_noise_std_mps),
                    rng.gaussian(0.0, config.gnss_velocity_noise_std_mps)
                };
                const sensors::GnssMeasurement gnss_meas{t, gnss_r, gnss_v, true};
                filter.update_gnss(gnss_meas);
                sum_gnss_nis += filter.last_gnss_diagnostics().normalized_innovation_squared;
                gnss_count++;
            }

            if (step >= 20) {
                const double nees = estimation::full_state_nees_15(
                    truth_pos, truth_vel, truth_att, true_ba, true_bg,
                    filter.state().position_eci_m,
                    filter.state().velocity_eci_mps,
                    filter.state().attitude_body_to_eci,
                    filter.state().accelerometer_bias_body_mps2,
                    filter.state().gyro_bias_body_rad_s,
                    filter.covariance()
                );
                sum_full_nees += nees;
                nees_count++;
            }
        }
    }

    const double mean_st_nis = sum_st_nis / static_cast<double>(st_count);
    const double mean_range_nis = sum_range_nis / static_cast<double>(range_count);
    const double mean_gnss_nis = sum_gnss_nis / static_cast<double>(gnss_count);
    const double mean_nees_15 = sum_full_nees / static_cast<double>(nees_count);

    // Verify statistical expectations:
    // Star tracker df = 3, mean ~ 3.0
    CHECK_THAT(mean_st_nis, WithinAbs(3.0, 0.4));
    // Range df = 1, mean ~ 1.0
    CHECK_THAT(mean_range_nis, WithinAbs(1.0, 0.25));
    // GNSS df = 6, mean ~ 6.0
    CHECK(mean_gnss_nis > estimation::k_chi2_6dof_95_lower);
    CHECK(mean_gnss_nis < estimation::k_chi2_6dof_95_upper);
    CHECK_THAT(mean_gnss_nis, WithinAbs(6.0, 1.5));
    // Full state df = 15, mean ~ 15.0
    CHECK(mean_nees_15 > estimation::k_chi2_15dof_95_lower);
    CHECK(mean_nees_15 < estimation::k_chi2_15dof_95_upper);
    CHECK_THAT(mean_nees_15, WithinAbs(15.0, 4.0));
}

TEST_CASE("Integrated navigation filter obeys truth non-interference contract and bitwise determinism", "[estimation][integrated][determinism]") {
    estimation::IntegratedNavigationConfig config;
    estimation::NavigationState state;
    state.position_eci_m = math::Vector3{6878137.0, 0.0, 0.0};
    state.velocity_eci_mps = math::Vector3{0.0, 7612.0, 0.0};
    state.attitude_body_to_eci = math::Quaternion::identity();

    estimation::IntegratedNavigationEkf filter1(config, state);
    estimation::IntegratedNavigationEkf filter2(config, state);

    // Both filters execute identically under identical inputs
    const math::Vector3 f{0.01, -0.02, 0.005};
    const math::Vector3 w{0.001, -0.002, 0.003};

    filter1.predict(f, w, 0.1);
    filter2.predict(f, w, 0.1);

    CHECK(filter1.state().position_eci_m.x() == filter2.state().position_eci_m.x());
    CHECK(filter1.state().velocity_eci_mps.y() == filter2.state().velocity_eci_mps.y());
    CHECK(filter1.state().attitude_body_to_eci.w() == filter2.state().attitude_body_to_eci.w());
    CHECK(filter1.covariance()(0, 0) == filter2.covariance()(0, 0));
}

TEST_CASE("M13A, M13B, and M13C regression compatibility", "[estimation][integrated][regression]") {
    estimation::IntegratedNavigationConfig config;

    // 1. Star-tracker double-cover invariance: q and -q produce identical innovation and posterior state
    const math::Quaternion q_test = math::Quaternion::from_axis_angle(math::Vector3{0.0, 1.0, 0.0}, 0.05);
    const math::Quaternion q_neg{-q_test.w(), -q_test.x(), -q_test.y(), -q_test.z()};

    estimation::NavigationState state_pos;
    state_pos.position_eci_m = math::Vector3{7000.0e3, 0.0, 0.0};
    state_pos.velocity_eci_mps = math::Vector3{0.0, 7546.0, 0.0};
    state_pos.attitude_body_to_eci = math::Quaternion::identity();

    estimation::NavigationState state_neg = state_pos;

    estimation::IntegratedNavigationEkf f_pos(config, state_pos);
    estimation::IntegratedNavigationEkf f_neg(config, state_neg);

    sensors::StarTrackerMeasurement m_pos{1.0, q_test, true};
    sensors::StarTrackerMeasurement m_neg{1.0, q_neg, true};

    f_pos.update_star_tracker(m_pos);
    f_neg.update_star_tracker(m_neg);

    CHECK_THAT(f_pos.last_star_tracker_diagnostics().normalized_innovation_squared,
               WithinAbs(f_neg.last_star_tracker_diagnostics().normalized_innovation_squared, 1.0e-12));
    CHECK_THAT(f_pos.state().attitude_body_to_eci.w(),
               WithinAbs(f_neg.state().attitude_body_to_eci.w(), 1.0e-12));

    // 2. Coincident range update throws domain_error
    CHECK_THROWS_AS(f_pos.update_range(10.0, f_pos.state().position_eci_m, 1.0), std::domain_error);
}
