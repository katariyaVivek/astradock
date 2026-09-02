// AstraDock M13D — Integrated Multi-Rate 15-State Navigation Demonstration Tool.
//
// Simulates a full orbital scenario with multi-rate sensors (IMU at 100 Hz, Star Tracker
// at 10 Hz, Range at 10 Hz, GNSS at 1 Hz), deterministic sensor outages, accelerometer
// and gyro biases, and exports telemetry and Monte Carlo summaries for verification.

#include "estimation/diagnostics.hpp"
#include "estimation/integrated_navigation_ekf.hpp"
#include "math/constants.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"
#include "orbit/cartesian_state.hpp"
#include "orbit/two_body_orbit.hpp"
#include "sensors/gnss.hpp"
#include "sensors/imu.hpp"
#include "sensors/range_sensor.hpp"
#include "sensors/sensor_common.hpp"
#include "sensors/star_tracker.hpp"

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace astradock;

int main() {
    std::cout << "=================================================================\n";
    std::cout << "AstraDock M13D — Integrated 15-State Navigation Demonstration\n";
    std::cout << "=================================================================\n";

    std::filesystem::create_directories("data");

    // Nominal scenario parameters
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const double r_orbit = 7000.0e3; // 7000 km circular orbit
    const double v_circ = std::sqrt(mu / r_orbit); // ~7546 m/s

    math::Vector3 truth_pos{r_orbit, 0.0, 0.0};
    math::Vector3 truth_vel{0.0, v_circ, 0.0};
    math::Quaternion truth_att = math::Quaternion::identity();

    const math::Vector3 truth_rate{0.02, -0.015, 0.01}; // ~27 mrad/s tumble
    const math::Vector3 truth_ba{0.005, -0.003, 0.002}; // 5 mm/s^2 accelerometer bias
    const math::Vector3 truth_bg{0.002, -0.001, 0.0015}; // 2 mrad/s gyro bias

    // Cooperative target station in orbit ahead
    const math::Vector3 target_pos = truth_pos + math::Vector3{1000.0, 0.0, 0.0};

    // Filter configuration
    estimation::IntegratedNavigationConfig config;
    config.accelerometer_noise_std_mps2 = 1.0e-3;
    config.gyro_noise_std_rad_s = 1.0e-3;
    config.gnss_position_noise_std_m = 5.0;
    config.gnss_velocity_noise_std_mps = 0.05;
    config.star_tracker_noise_std_rad = 1.0e-3;
    config.range_noise_std_m = 1.0;

    // Deliberately perturbed initial estimate (truth non-interference)
    estimation::NavigationState est_state;
    est_state.position_eci_m = truth_pos + math::Vector3{35.0, -25.0, 15.0};
    est_state.velocity_eci_mps = truth_vel + math::Vector3{0.3, -0.4, 0.2};
    est_state.attitude_body_to_eci = (truth_att * math::Quaternion::from_axis_angle(
        math::Vector3{0.577, 0.577, 0.577}.normalized(), 0.03)).normalized(); // ~1.7 deg
    est_state.accelerometer_bias_body_mps2 = math::Vector3{0.0, 0.0, 0.0};
    est_state.gyro_bias_body_rad_s = math::Vector3{0.0, 0.0, 0.0};

    estimation::IntegratedNavigationEkf filter(config, est_state);

    sensors::DeterministicRng rng(424242);

    // Timeline: 120 seconds, dt_sim = 0.01 s (100 Hz IMU)
    const double t_end = 120.0;
    const double dt_imu = 0.01;
    const std::size_t n_steps = static_cast<std::size_t>(t_end / dt_imu);

    // Dropouts:
    // Star tracker dropout: 30 s to 45 s (15 s outage)
    // GNSS dropout: 60 s to 90 s (30 s outage)
    const double st_dropout_start = 30.0;
    const double st_dropout_end = 45.0;
    const double gnss_dropout_start = 60.0;
    const double gnss_dropout_end = 90.0;

    std::ofstream tel_file("data/m13d_integrated_telemetry.csv");
    tel_file << "time_s,truth_position_eci_x_m,truth_position_eci_y_m,truth_position_eci_z_m,"
             << "estimated_position_eci_x_m,estimated_position_eci_y_m,estimated_position_eci_z_m,"
             << "truth_velocity_eci_x_mps,truth_velocity_eci_y_mps,truth_velocity_eci_z_mps,"
             << "estimated_velocity_eci_x_mps,estimated_velocity_eci_y_mps,estimated_velocity_eci_z_mps,"
             << "attitude_error_rad,"
             << "truth_ba_x_mps2,truth_ba_y_mps2,truth_ba_z_mps2,"
             << "estimated_ba_x_mps2,estimated_ba_y_mps2,estimated_ba_z_mps2,"
             << "truth_bg_x_rad_s,truth_bg_y_rad_s,truth_bg_z_rad_s,"
             << "estimated_bg_x_rad_s,estimated_bg_y_rad_s,estimated_bg_z_rad_s,"
             << "sigma_position,sigma_velocity,sigma_attitude,sigma_ba,sigma_bg,"
             << "gnss_valid,star_tracker_valid,range_valid,"
             << "gnss_NIS,star_tracker_NIS,range_NIS,full_state_NEES\n";

    double last_gnss_nis = 0.0;
    double last_st_nis = 0.0;
    double last_range_nis = 0.0;

    for (std::size_t step = 0; step <= n_steps; ++step) {
        const double t = step * dt_imu;

        // Advance orbital truth
        if (step > 0) {
            const math::Vector3 g_truth = dynamics::two_body_acceleration(truth_pos, mu);
            truth_pos = truth_pos + truth_vel * dt_imu + g_truth * (0.5 * dt_imu * dt_imu);
            truth_vel = truth_vel + g_truth * dt_imu;

            const double w_mag = truth_rate.norm();
            truth_att = (truth_att * math::Quaternion::from_axis_angle(truth_rate / w_mag, w_mag * dt_imu)).normalized();
        }

        // Discrete IMU measurement (free-fall in orbit -> true specific force = 0)
        const double disc_acc_std = config.accelerometer_noise_std_mps2 / std::sqrt(dt_imu);
        const double disc_gyr_std = config.gyro_noise_std_rad_s / std::sqrt(dt_imu);

        const math::Vector3 meas_f = truth_ba + math::Vector3{
            rng.gaussian(0.0, disc_acc_std),
            rng.gaussian(0.0, disc_acc_std),
            rng.gaussian(0.0, disc_acc_std)
        };
        const math::Vector3 meas_w = truth_rate + truth_bg + math::Vector3{
            rng.gaussian(0.0, disc_gyr_std),
            rng.gaussian(0.0, disc_gyr_std),
            rng.gaussian(0.0, disc_gyr_std)
        };

        if (step > 0) {
            filter.predict(meas_f, meas_w, dt_imu);
        }

        bool gnss_valid = false;
        bool st_valid = false;
        bool range_valid = false;

        // Star Tracker at 10 Hz (every 10 steps)
        if (step % 10 == 0) {
            if (t < st_dropout_start || t > st_dropout_end) {
                const double st_ang = rng.gaussian(0.0, config.star_tracker_noise_std_rad);
                const math::Vector3 st_u = rng.uniform_unit_vector3();
                const math::Quaternion st_noise = math::Quaternion::from_axis_angle(st_u, st_ang);
                const sensors::StarTrackerMeasurement st_meas{t, (truth_att * st_noise).normalized(), true};
                filter.update_star_tracker(st_meas);
                last_st_nis = filter.last_star_tracker_diagnostics().normalized_innovation_squared;
                st_valid = true;
            }
        }

        // Range at 10 Hz (every 10 steps)
        if (step % 10 == 0) {
            const double true_rho = (target_pos - truth_pos).norm();
            const double meas_rho = true_rho + rng.gaussian(0.0, config.range_noise_std_m);
            filter.update_range(meas_rho, target_pos, config.range_noise_std_m);
            last_range_nis = filter.last_range_diagnostics().normalized_innovation_squared;
            range_valid = true;
        }

        // GNSS at 1 Hz (every 100 steps)
        if (step % 100 == 0) {
            if (t < gnss_dropout_start || t > gnss_dropout_end) {
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
                last_gnss_nis = filter.last_gnss_diagnostics().normalized_innovation_squared;
                gnss_valid = true;
            }
        }

        // Compute diagnostics and record telemetry every 10 steps (10 Hz export)
        if (step % 10 == 0) {
            const auto& cov = filter.covariance();
            const double sig_pos = std::sqrt(cov(0, 0) + cov(1, 1) + cov(2, 2));
            const double sig_vel = std::sqrt(cov(3, 3) + cov(4, 4) + cov(5, 5));
            const double sig_att = std::sqrt(cov(6, 6) + cov(7, 7) + cov(8, 8));
            const double sig_ba = std::sqrt(cov(9, 9) + cov(10, 10) + cov(11, 11));
            const double sig_bg = std::sqrt(cov(12, 12) + cov(13, 13) + cov(14, 14));

            const double att_err = estimation::attitude_error_vector_rad(
                truth_att, filter.state().attitude_body_to_eci).norm();

            const double full_nees = estimation::full_state_nees_15(
                truth_pos, truth_vel, truth_att, truth_ba, truth_bg,
                filter.state().position_eci_m,
                filter.state().velocity_eci_mps,
                filter.state().attitude_body_to_eci,
                filter.state().accelerometer_bias_body_mps2,
                filter.state().gyro_bias_body_rad_s,
                cov
            );

            tel_file << std::fixed << std::setprecision(4) << t << ","
                     << std::setprecision(3)
                     << truth_pos.x() << "," << truth_pos.y() << "," << truth_pos.z() << ","
                     << filter.state().position_eci_m.x() << "," << filter.state().position_eci_m.y() << "," << filter.state().position_eci_m.z() << ","
                     << std::setprecision(5)
                     << truth_vel.x() << "," << truth_vel.y() << "," << truth_vel.z() << ","
                     << filter.state().velocity_eci_mps.x() << "," << filter.state().velocity_eci_mps.y() << "," << filter.state().velocity_eci_mps.z() << ","
                     << std::setprecision(6) << att_err << ","
                     << truth_ba.x() << "," << truth_ba.y() << "," << truth_ba.z() << ","
                     << filter.state().accelerometer_bias_body_mps2.x() << "," << filter.state().accelerometer_bias_body_mps2.y() << "," << filter.state().accelerometer_bias_body_mps2.z() << ","
                     << truth_bg.x() << "," << truth_bg.y() << "," << truth_bg.z() << ","
                     << filter.state().gyro_bias_body_rad_s.x() << "," << filter.state().gyro_bias_body_rad_s.y() << "," << filter.state().gyro_bias_body_rad_s.z() << ","
                     << sig_pos << "," << sig_vel << "," << sig_att << "," << sig_ba << "," << sig_bg << ","
                     << (gnss_valid ? 1 : 0) << "," << (st_valid ? 1 : 0) << "," << (range_valid ? 1 : 0) << ","
                     << last_gnss_nis << "," << last_st_nis << "," << last_range_nis << "," << full_nees << "\n";
        }
    }
    tel_file.close();
    std::cout << "Exported nominal telemetry to data/m13d_integrated_telemetry.csv\n";

    // 2. Sensor combinations comparison study
    std::cout << "Running sensor combinations study...\n";
    std::ofstream comb_file("data/m13d_sensor_combinations.csv");
    comb_file << "configuration,position_rmse_m,velocity_rmse_mps,attitude_rmse_rad,gyro_bias_rmse_rad_s\n";

    struct ConfigSpec {
        std::string name;
        bool use_imu;
        bool use_gnss;
        bool use_st;
        bool use_range;
    };

    const std::vector<ConfigSpec> specs = {
        {"GNSS Only", false, true, false, false},
        {"IMU + GNSS", true, true, false, false},
        {"IMU + GNSS + Star Tracker", true, true, true, false},
        {"Full Suite (IMU+GNSS+ST+Range)", true, true, true, true}
    };

    for (const auto& sp : specs) {
        sensors::DeterministicRng test_rng(9999);
        math::Vector3 s_truth_pos{r_orbit, 0.0, 0.0};
        math::Vector3 s_truth_vel{0.0, v_circ, 0.0};
        math::Quaternion s_truth_att = math::Quaternion::identity();

        estimation::NavigationState s_est;
        s_est.position_eci_m = s_truth_pos + math::Vector3{20.0, -15.0, 10.0};
        s_est.velocity_eci_mps = s_truth_vel + math::Vector3{0.2, -0.2, 0.1};
        s_est.attitude_body_to_eci = (s_truth_att * math::Quaternion::from_axis_angle(math::Vector3{0.0, 1.0, 0.0}, 0.02)).normalized();

        estimation::IntegratedNavigationEkf s_filter(config, s_est);

        double sum_pos_err2 = 0.0;
        double sum_vel_err2 = 0.0;
        double sum_att_err2 = 0.0;
        double sum_bg_err2 = 0.0;
        std::size_t sample_count = 0;

        for (std::size_t step = 1; step <= 5000; ++step) { // 50 seconds
            const double t = step * dt_imu;
            const math::Vector3 g_truth = dynamics::two_body_acceleration(s_truth_pos, mu);
            s_truth_pos = s_truth_pos + s_truth_vel * dt_imu + g_truth * (0.5 * dt_imu * dt_imu);
            s_truth_vel = s_truth_vel + g_truth * dt_imu;
            const double w_mag = truth_rate.norm();
            s_truth_att = (s_truth_att * math::Quaternion::from_axis_angle(truth_rate / w_mag, w_mag * dt_imu)).normalized();

            const double disc_acc_std = config.accelerometer_noise_std_mps2 / std::sqrt(dt_imu);
            const double disc_gyr_std = config.gyro_noise_std_rad_s / std::sqrt(dt_imu);

            const math::Vector3 meas_f = truth_ba + math::Vector3{
                test_rng.gaussian(0.0, disc_acc_std),
                test_rng.gaussian(0.0, disc_acc_std),
                test_rng.gaussian(0.0, disc_acc_std)
            };
            const math::Vector3 meas_w = truth_rate + truth_bg + math::Vector3{
                test_rng.gaussian(0.0, disc_gyr_std),
                test_rng.gaussian(0.0, disc_gyr_std),
                test_rng.gaussian(0.0, disc_gyr_std)
            };

            if (sp.use_imu) {
                s_filter.predict(meas_f, meas_w, dt_imu);
            } else {
                // Dynamics only prediction
                s_filter.predict(math::Vector3{0.0, 0.0, 0.0}, math::Vector3{0.0, 0.0, 0.0}, dt_imu);
            }

            if (sp.use_st && step % 10 == 0) {
                const double st_ang = test_rng.gaussian(0.0, config.star_tracker_noise_std_rad);
                const math::Vector3 st_u = test_rng.uniform_unit_vector3();
                const math::Quaternion st_noise = math::Quaternion::from_axis_angle(st_u, st_ang);
                const sensors::StarTrackerMeasurement st_meas{t, (s_truth_att * st_noise).normalized(), true};
                s_filter.update_star_tracker(st_meas);
            }

            if (sp.use_range && step % 10 == 0) {
                const double true_rho = (target_pos - s_truth_pos).norm();
                const double meas_rho = true_rho + test_rng.gaussian(0.0, config.range_noise_std_m);
                s_filter.update_range(meas_rho, target_pos, config.range_noise_std_m);
            }

            if (sp.use_gnss && step % 100 == 0) {
                const math::Vector3 gnss_r = s_truth_pos + math::Vector3{
                    test_rng.gaussian(0.0, config.gnss_position_noise_std_m),
                    test_rng.gaussian(0.0, config.gnss_position_noise_std_m),
                    test_rng.gaussian(0.0, config.gnss_position_noise_std_m)
                };
                const math::Vector3 gnss_v = s_truth_vel + math::Vector3{
                    test_rng.gaussian(0.0, config.gnss_velocity_noise_std_mps),
                    test_rng.gaussian(0.0, config.gnss_velocity_noise_std_mps),
                    test_rng.gaussian(0.0, config.gnss_velocity_noise_std_mps)
                };
                const sensors::GnssMeasurement gnss_meas{t, gnss_r, gnss_v, true};
                s_filter.update_gnss(gnss_meas);
            }

            if (t >= 20.0) { // steady state
                const double pos_err = (s_truth_pos - s_filter.state().position_eci_m).norm();
                const double vel_err = (s_truth_vel - s_filter.state().velocity_eci_mps).norm();
                const double att_err = estimation::attitude_error_vector_rad(
                    s_truth_att, s_filter.state().attitude_body_to_eci).norm();
                const double bg_err = (truth_bg - s_filter.state().gyro_bias_body_rad_s).norm();

                sum_pos_err2 += pos_err * pos_err;
                sum_vel_err2 += vel_err * vel_err;
                sum_att_err2 += att_err * att_err;
                sum_bg_err2 += bg_err * bg_err;
                sample_count++;
            }
        }

        const double pos_rmse = std::sqrt(sum_pos_err2 / static_cast<double>(sample_count));
        const double vel_rmse = std::sqrt(sum_vel_err2 / static_cast<double>(sample_count));
        const double att_rmse = std::sqrt(sum_att_err2 / static_cast<double>(sample_count));
        const double bg_rmse = std::sqrt(sum_bg_err2 / static_cast<double>(sample_count));

        comb_file << sp.name << "," << pos_rmse << "," << vel_rmse << "," << att_rmse << "," << bg_rmse << "\n";
        std::cout << "  " << std::left << std::setw(32) << sp.name
                  << " | Pos RMSE: " << std::setw(7) << std::fixed << std::setprecision(3) << pos_rmse << " m"
                  << " | Vel RMSE: " << std::setw(7) << std::fixed << std::setprecision(4) << vel_rmse << " m/s"
                  << " | Att RMSE: " << std::setw(7) << std::fixed << std::setprecision(5) << att_rmse << " rad\n";
    }
    comb_file.close();
    std::cout << "Exported sensor combinations to data/m13d_sensor_combinations.csv\n";

    // 3. 100-seed Monte Carlo study
    std::cout << "Running 100-seed Monte Carlo study...\n";
    std::ofstream mc_file("data/m13d_integrated_monte_carlo.csv");
    mc_file << "seed,pos_rmse_m,vel_rmse_mps,att_rmse_rad,ba_rmse_mps2,bg_rmse_rad_s,"
            << "mean_full_nees,mean_gnss_nis,mean_st_nis,mean_range_nis\n";

    for (std::uint64_t seed = 1; seed <= 100; ++seed) {
        sensors::DeterministicRng mc_rng(seed * 7771);

        math::Vector3 mc_truth_pos{r_orbit, 0.0, 0.0};
        math::Vector3 mc_truth_vel{0.0, v_circ, 0.0};
        math::Quaternion mc_truth_att = math::Quaternion::identity();

        estimation::NavigationState mc_est;
        mc_est.position_eci_m = mc_truth_pos + math::Vector3{
            mc_rng.gaussian(0.0, 20.0),
            mc_rng.gaussian(0.0, 20.0),
            mc_rng.gaussian(0.0, 20.0)
        };
        mc_est.velocity_eci_mps = mc_truth_vel + math::Vector3{
            mc_rng.gaussian(0.0, 0.2),
            mc_rng.gaussian(0.0, 0.2),
            mc_rng.gaussian(0.0, 0.2)
        };
        const double init_att_err = mc_rng.gaussian(0.0, 0.02);
        mc_est.attitude_body_to_eci = (mc_truth_att * math::Quaternion::from_axis_angle(
            mc_rng.uniform_unit_vector3(), init_att_err)).normalized();
        mc_est.accelerometer_bias_body_mps2 = math::Vector3{0.0, 0.0, 0.0};
        mc_est.gyro_bias_body_rad_s = math::Vector3{0.0, 0.0, 0.0};

        estimation::IntegratedNavigationEkf mc_filter(config, mc_est);

        double mc_pos_err2 = 0.0;
        double mc_vel_err2 = 0.0;
        double mc_att_err2 = 0.0;
        double mc_ba_err2 = 0.0;
        double mc_bg_err2 = 0.0;

        double mc_sum_nees = 0.0;
        double mc_sum_gnss_nis = 0.0;
        double mc_sum_st_nis = 0.0;
        double mc_sum_range_nis = 0.0;

        std::size_t mc_gnss_n = 0;
        std::size_t mc_st_n = 0;
        std::size_t mc_range_n = 0;
        std::size_t mc_sample_n = 0;

        for (std::size_t step = 1; step <= 3000; ++step) { // 30 seconds
            const double t = step * dt_imu;
            const math::Vector3 g_truth = dynamics::two_body_acceleration(mc_truth_pos, mu);
            mc_truth_pos = mc_truth_pos + mc_truth_vel * dt_imu + g_truth * (0.5 * dt_imu * dt_imu);
            mc_truth_vel = mc_truth_vel + g_truth * dt_imu;
            const double w_mag = truth_rate.norm();
            mc_truth_att = (mc_truth_att * math::Quaternion::from_axis_angle(truth_rate / w_mag, w_mag * dt_imu)).normalized();

            const double disc_acc_std = config.accelerometer_noise_std_mps2 / std::sqrt(dt_imu);
            const double disc_gyr_std = config.gyro_noise_std_rad_s / std::sqrt(dt_imu);

            const math::Vector3 meas_f = truth_ba + math::Vector3{
                mc_rng.gaussian(0.0, disc_acc_std),
                mc_rng.gaussian(0.0, disc_acc_std),
                mc_rng.gaussian(0.0, disc_acc_std)
            };
            const math::Vector3 meas_w = truth_rate + truth_bg + math::Vector3{
                mc_rng.gaussian(0.0, disc_gyr_std),
                mc_rng.gaussian(0.0, disc_gyr_std),
                mc_rng.gaussian(0.0, disc_gyr_std)
            };

            mc_filter.predict(meas_f, meas_w, dt_imu);

            if (step % 10 == 0) {
                const double st_ang = mc_rng.gaussian(0.0, config.star_tracker_noise_std_rad);
                const math::Vector3 st_u = mc_rng.uniform_unit_vector3();
                const math::Quaternion st_noise = math::Quaternion::from_axis_angle(st_u, st_ang);
                const sensors::StarTrackerMeasurement st_meas{t, (mc_truth_att * st_noise).normalized(), true};
                mc_filter.update_star_tracker(st_meas);
                mc_sum_st_nis += mc_filter.last_star_tracker_diagnostics().normalized_innovation_squared;
                mc_st_n++;

                const double true_rho = (target_pos - mc_truth_pos).norm();
                const double meas_rho = true_rho + mc_rng.gaussian(0.0, config.range_noise_std_m);
                mc_filter.update_range(meas_rho, target_pos, config.range_noise_std_m);
                mc_sum_range_nis += mc_filter.last_range_diagnostics().normalized_innovation_squared;
                mc_range_n++;
            }

            if (step % 100 == 0) {
                const math::Vector3 gnss_r = mc_truth_pos + math::Vector3{
                    mc_rng.gaussian(0.0, config.gnss_position_noise_std_m),
                    mc_rng.gaussian(0.0, config.gnss_position_noise_std_m),
                    mc_rng.gaussian(0.0, config.gnss_position_noise_std_m)
                };
                const math::Vector3 gnss_v = mc_truth_vel + math::Vector3{
                    mc_rng.gaussian(0.0, config.gnss_velocity_noise_std_mps),
                    mc_rng.gaussian(0.0, config.gnss_velocity_noise_std_mps),
                    mc_rng.gaussian(0.0, config.gnss_velocity_noise_std_mps)
                };
                const sensors::GnssMeasurement gnss_meas{t, gnss_r, gnss_v, true};
                mc_filter.update_gnss(gnss_meas);
                mc_sum_gnss_nis += mc_filter.last_gnss_diagnostics().normalized_innovation_squared;
                mc_gnss_n++;
            }

            if (t >= 15.0) {
                const double pos_err = (mc_truth_pos - mc_filter.state().position_eci_m).norm();
                const double vel_err = (mc_truth_vel - mc_filter.state().velocity_eci_mps).norm();
                const double att_err = estimation::attitude_error_vector_rad(
                    mc_truth_att, mc_filter.state().attitude_body_to_eci).norm();
                const double ba_err = (truth_ba - mc_filter.state().accelerometer_bias_body_mps2).norm();
                const double bg_err = (truth_bg - mc_filter.state().gyro_bias_body_rad_s).norm();

                mc_pos_err2 += pos_err * pos_err;
                mc_vel_err2 += vel_err * vel_err;
                mc_att_err2 += att_err * att_err;
                mc_ba_err2 += ba_err * ba_err;
                mc_bg_err2 += bg_err * bg_err;

                const double nees = estimation::full_state_nees_15(
                    mc_truth_pos, mc_truth_vel, mc_truth_att, truth_ba, truth_bg,
                    mc_filter.state().position_eci_m,
                    mc_filter.state().velocity_eci_mps,
                    mc_filter.state().attitude_body_to_eci,
                    mc_filter.state().accelerometer_bias_body_mps2,
                    mc_filter.state().gyro_bias_body_rad_s,
                    mc_filter.covariance()
                );
                mc_sum_nees += nees;
                mc_sample_n++;
            }
        }

        mc_file << seed << ","
                << std::sqrt(mc_pos_err2 / static_cast<double>(mc_sample_n)) << ","
                << std::sqrt(mc_vel_err2 / static_cast<double>(mc_sample_n)) << ","
                << std::sqrt(mc_att_err2 / static_cast<double>(mc_sample_n)) << ","
                << std::sqrt(mc_ba_err2 / static_cast<double>(mc_sample_n)) << ","
                << std::sqrt(mc_bg_err2 / static_cast<double>(mc_sample_n)) << ","
                << (mc_sum_nees / static_cast<double>(mc_sample_n)) << ","
                << (mc_sum_gnss_nis / static_cast<double>(mc_gnss_n)) << ","
                << (mc_sum_st_nis / static_cast<double>(mc_st_n)) << ","
                << (mc_sum_range_nis / static_cast<double>(mc_range_n)) << "\n";
    }
    mc_file.close();
    std::cout << "Exported 100-seed Monte Carlo results to data/m13d_integrated_monte_carlo.csv\n";
    std::cout << "=================================================================\n";
    std::cout << "M13D demonstration completed successfully.\n";
    std::cout << "=================================================================\n";

    return 0;
}
