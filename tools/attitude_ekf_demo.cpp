// AstraDock M13C — Attitude Error-State EKF & Range Update demonstration and validation driver.
//
// Produces deterministic CSV telemetry and Monte Carlo datasets consumed by the Python
// analysis layer (python/analysis/plot_attitude_estimation.py) and independent audit oracles:
//
//   data/m13c_attitude_telemetry.csv     100 Hz attitude EKF telemetry with star tracker dropout
//   data/m13c_attitude_monte_carlo.csv   100-seed attitude EKF Monte Carlo summary statistics
//   data/m13c_range_telemetry.csv        Range measurement updates on translational state
//   data/m13c_range_monte_carlo.csv      100-seed range update Monte Carlo summary statistics

#include "attitude/principal_inertia.hpp"
#include "attitude/rigid_body.hpp"
#include "attitude/rotational_state.hpp"
#include "estimation/attitude_ekf.hpp"
#include "estimation/diagnostics.hpp"
#include "estimation/range_update.hpp"
#include "estimation/translational_ekf.hpp"
#include "math/constants.hpp"
#include "orbit/cartesian_state.hpp"
#include "orbit/two_body_orbit.hpp"
#include "sensors/sensor_common.hpp"
#include "sensors/star_tracker.hpp"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace astradock;

namespace {

void run_attitude_telemetry(const std::string& output_path) {
    std::ofstream csv(output_path);
    if (!csv.is_open()) {
        std::cerr << "Failed to open " << output_path << std::endl;
        return;
    }

    csv << "time_s,truth_attitude_error_rad,estimated_attitude_error_rad,"
        << "gyro_x,gyro_y,gyro_z,estimated_bias_x,estimated_bias_y,estimated_bias_z,"
        << "attitude_innovation_x,attitude_innovation_y,attitude_innovation_z,"
        << "attitude_NIS,sigma_attitude_x,sigma_attitude_y,sigma_attitude_z,"
        << "sigma_bias_x,sigma_bias_y,sigma_bias_z\n";

    // Spacecraft inertia and initial truth rotational state
    const attitude::PrincipalInertia inertia{10.0, 20.0, 30.0};
    math::Quaternion truth_q = math::Quaternion::identity();
    math::Vector3 truth_omega{0.04, -0.03, 0.02}; // rad/s tumbling
    const math::Vector3 true_bias{0.005, -0.003, 0.002}; // 5, -3, 2 mrad/s

    // Attitude EKF configuration
    estimation::AttitudeEkfConfig config;
    config.gyro_noise_std_rad_s = 1.0e-3;
    config.star_tracker_noise_std_rad = 1.0e-3;
    config.initial_attitude_error_std_rad = 0.05; // ~2.9 deg
    config.initial_gyro_bias_std_rad_s = 0.005;

    // Initial estimate offset
    const math::Vector3 init_err{0.03, -0.025, 0.035};
    const double init_err_norm = init_err.norm();
    const math::Quaternion q_init_err = math::Quaternion::from_axis_angle(
        init_err / init_err_norm, init_err_norm);

    estimation::AttitudeEstimate est;
    est.nominal_orientation = q_init_err;
    est.gyro_bias_rad_s = math::Vector3{0.0, 0.0, 0.0};
    estimation::AttitudeEkf filter(config, est);

    // Star tracker simulation
    sensors::StarTrackerConfig st_cfg;
    st_cfg.noise_std_rad = config.star_tracker_noise_std_rad;
    st_cfg.sample_period_s = 0.1;
    st_cfg.dropouts = {{30.0, 45.0}}; // 15-second star tracker outage
    st_cfg.random_seed = 1001;
    sensors::StarTrackerSensor star_tracker(st_cfg);

    sensors::DeterministicRng gyro_rng(2002);

    const double dt_gyro = 0.01; // 100 Hz
    const double total_duration = 75.0; // 75 seconds

    // Raw gyro-propagated attitude (open-loop baseline without updates)
    math::Quaternion raw_gyro_q = q_init_err;

    csv << std::setprecision(8) << std::fixed;

    for (double t = 0.0; t <= total_duration; t += dt_gyro) {
        // True dynamics: Euler's rigid-body equations + kinematics via RK4
        const auto rot_deriv = [&inertia](double, const attitude::RotationalState& s) {
            return attitude::rotational_state_derivative(0.0, s, inertia, math::Vector3{});
        };
        const attitude::RotationalState current_rot{truth_q, truth_omega};
        const attitude::RotationalState next_rot = numerics::rk4_step(
            t, current_rot, dt_gyro, rot_deriv);
        truth_q = next_rot.orientation.normalized();
        truth_omega = next_rot.angular_velocity_rad_per_s;

        // Noisy, biased gyro measurement
        const double discrete_gyro_sigma = config.gyro_noise_std_rad_s / std::sqrt(dt_gyro);
        const math::Vector3 gyro_noise{
            gyro_rng.gaussian(0.0, discrete_gyro_sigma),
            gyro_rng.gaussian(0.0, discrete_gyro_sigma),
            gyro_rng.gaussian(0.0, discrete_gyro_sigma)
        };
        const math::Vector3 gyro_meas = truth_omega + true_bias + gyro_noise;

        // Propagate raw uncorrected attitude
        raw_gyro_q = estimation::propagate_nominal_quaternion(raw_gyro_q, gyro_meas, dt_gyro);

        // Filter prediction
        filter.predict(gyro_meas, dt_gyro);

        // Star tracker update at 10 Hz
        math::Vector3 inno{};
        double nis = 0.0;
        if (std::fmod(t + 1.0e-6, 0.1) < dt_gyro) {
            // Truth spacecraft state wrapper for sensor
            spacecraft::SpacecraftState truth_state{
                orbit::CartesianState{},
                attitude::RotationalState{truth_q, truth_omega}
            };
            const auto st_meas = star_tracker.measure(t, truth_state);
            if (st_meas.valid) {
                filter.update_star_tracker(st_meas);
                const auto& diag = filter.last_star_tracker_diagnostics();
                inno = diag.innovation_rad;
                nis = diag.normalized_innovation_squared;
            }
        }

        // Output decimation: export every 0.05 s (20 Hz)
        if (std::fmod(t + 1.0e-6, 0.05) < dt_gyro) {
            const double raw_att_err = estimation::attitude_error_vector_rad(truth_q, raw_gyro_q).norm();
            const double ekf_att_err = estimation::attitude_error_vector_rad(truth_q, filter.estimate().nominal_orientation).norm();

            const auto& cov = filter.covariance();
            const double sig_att_x = std::sqrt(std::max(0.0, cov(0, 0)));
            const double sig_att_y = std::sqrt(std::max(0.0, cov(1, 1)));
            const double sig_att_z = std::sqrt(std::max(0.0, cov(2, 2)));
            const double sig_b_x = std::sqrt(std::max(0.0, cov(3, 3)));
            const double sig_b_y = std::sqrt(std::max(0.0, cov(4, 4)));
            const double sig_b_z = std::sqrt(std::max(0.0, cov(5, 5)));

            const auto& est_b = filter.estimate().gyro_bias_rad_s;

            csv << t << ","
                << raw_att_err << ","
                << ekf_att_err << ","
                << gyro_meas.x() << "," << gyro_meas.y() << "," << gyro_meas.z() << ","
                << est_b.x() << "," << est_b.y() << "," << est_b.z() << ","
                << inno.x() << "," << inno.y() << "," << inno.z() << ","
                << nis << ","
                << sig_att_x << "," << sig_att_y << "," << sig_att_z << ","
                << sig_b_x << "," << sig_b_y << "," << sig_b_z << "\n";
        }
    }

    std::cout << "Attitude telemetry written to " << output_path << std::endl;
}

void run_attitude_monte_carlo(const std::string& output_path) {
    std::ofstream csv(output_path);
    if (!csv.is_open()) {
        std::cerr << "Failed to open " << output_path << std::endl;
        return;
    }

    csv << "seed,attitude_rmse_rad,gyro_bias_rmse_rad_s,mean_nis,mean_nees\n";
    csv << std::setprecision(8) << std::fixed;

    const std::size_t N_seeds = 100;
    const double dt_gyro = 0.02; // 50 Hz
    const double duration = 25.0; // 25 s

    const attitude::PrincipalInertia inertia{10.0, 20.0, 30.0};
    const math::Vector3 true_rate{0.03, -0.02, 0.025};
    const math::Vector3 true_bias{0.004, -0.003, 0.002};

    estimation::AttitudeEkfConfig config;
    config.gyro_noise_std_rad_s = 5.0e-4;
    config.star_tracker_noise_std_rad = 5.0e-4;
    config.initial_attitude_error_std_rad = 0.01;
    config.initial_gyro_bias_std_rad_s = 0.005;

    for (std::uint64_t seed = 1; seed <= N_seeds; ++seed) {
        sensors::DeterministicRng rng(seed * 100);

        // Perturb initial estimate
        const math::Vector3 init_err{
            rng.gaussian(0.0, config.initial_attitude_error_std_rad),
            rng.gaussian(0.0, config.initial_attitude_error_std_rad),
            rng.gaussian(0.0, config.initial_attitude_error_std_rad)
        };
        const double err_norm = init_err.norm();
        const math::Quaternion q_init = (err_norm > 0.0)
            ? math::Quaternion::from_axis_angle(init_err / err_norm, err_norm)
            : math::Quaternion::identity();

        estimation::AttitudeEstimate est{q_init, math::Vector3{0.0, 0.0, 0.0}};
        estimation::AttitudeEkf filter(config, est);

        math::Quaternion truth_q = math::Quaternion::identity();

        double sum_sq_att_err = 0.0;
        double sum_sq_bias_err = 0.0;
        double sum_nis = 0.0;
        double sum_nees = 0.0;
        std::size_t count_updates = 0;
        std::size_t count_steps = 0;

        for (double t = 0.0; t <= duration; t += dt_gyro) {
            // Advance truth
            const double ang = true_rate.norm() * dt_gyro;
            truth_q = (truth_q * math::Quaternion::from_axis_angle(true_rate.normalized(), ang)).normalized();

            const double discrete_gyro_sigma = config.gyro_noise_std_rad_s / std::sqrt(dt_gyro);
            const math::Vector3 gyro_noise{
                rng.gaussian(0.0, discrete_gyro_sigma),
                rng.gaussian(0.0, discrete_gyro_sigma),
                rng.gaussian(0.0, discrete_gyro_sigma)
            };
            filter.predict(true_rate + true_bias + gyro_noise, dt_gyro);

            // Update at 10 Hz
            if (std::fmod(t + 1.0e-6, 0.1) < dt_gyro) {
                const double st_err_ang = rng.gaussian(0.0, config.star_tracker_noise_std_rad);
                const math::Vector3 st_ax = rng.uniform_unit_vector3();
                const math::Quaternion st_noise_q = math::Quaternion::from_axis_angle(st_ax, st_err_ang);
                const math::Quaternion meas_q = (truth_q * st_noise_q).normalized();

                sensors::StarTrackerMeasurement meas{t, meas_q, true};
                filter.update_star_tracker(meas);

                if (t > 10.0) { // steady-state
                    sum_nis += filter.last_star_tracker_diagnostics().normalized_innovation_squared;
                    math::Matrix<3, 3> p_att;
                    for (std::size_t r = 0; r < 3; ++r) {
                        for (std::size_t c = 0; c < 3; ++c) {
                            p_att(r, c) = filter.covariance()(r, c);
                        }
                    }
                    sum_nees += estimation::attitude_nees(truth_q, filter.estimate().nominal_orientation, p_att);
                    count_updates++;
                }
            }

            if (t > 10.0) {
                const math::Vector3 att_err = estimation::attitude_error_vector_rad(
                    truth_q, filter.estimate().nominal_orientation);
                sum_sq_att_err += att_err.squared_norm();
                const math::Vector3 bias_err = true_bias - filter.estimate().gyro_bias_rad_s;
                sum_sq_bias_err += bias_err.squared_norm();
                count_steps++;
            }
        }

        const double att_rmse = std::sqrt(sum_sq_att_err / static_cast<double>(count_steps));
        const double bias_rmse = std::sqrt(sum_sq_bias_err / static_cast<double>(count_steps));
        const double mean_nis = sum_nis / static_cast<double>(count_updates);
        const double mean_nees = sum_nees / static_cast<double>(count_updates);

        csv << seed << "," << att_rmse << "," << bias_rmse << "," << mean_nis << "," << mean_nees << "\n";
    }

    std::cout << "Attitude Monte Carlo written to " << output_path << std::endl;
}

void run_range_telemetry(const std::string& output_path) {
    std::ofstream csv(output_path);
    if (!csv.is_open()) {
        std::cerr << "Failed to open " << output_path << std::endl;
        return;
    }

    csv << "time_s,truth_range_m,predicted_range_m,measured_range_m,range_innovation_m,range_NIS,position_error_m\n";
    csv << std::setprecision(8) << std::fixed;

    // Target in circular 500 km LEO
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const double r_orb = constants::earth_reference_radius_m + 500.0e3;
    const double v_circ = std::sqrt(mu / r_orb);

    // Target at nominal circular state
    orbit::CartesianState target_state{
        math::Vector3{r_orb, 0.0, 0.0},
        math::Vector3{0.0, v_circ, 0.0}
    };

    // Chaser spacecraft trailing target by 500 m along along-track direction (+Y)
    orbit::CartesianState chaser_truth{
        math::Vector3{r_orb, -500.0, 0.0},
        math::Vector3{0.0, v_circ + 0.1, 0.0} // slow closing velocity 0.1 m/s
    };

    // Filter setup with initial position error
    estimation::TranslationalEkfConfig config;
    config.gravitational_parameter_m3_per_s2 = mu;
    config.acceleration_noise_std_mps2 = 1.0e-3;

    // Initial estimate has 15 m error along X and Y
    orbit::CartesianState chaser_est{
        chaser_truth.position + math::Vector3{15.0, -10.0, 5.0},
        chaser_truth.velocity
    };

    estimation::TranslationalCovariance cov;
    for (std::size_t i = 0; i < 6; ++i) {
        cov(i, i) = (i < 3) ? 400.0 : 0.01; // (20 m)^2 position uncertainty
    }
    estimation::TranslationalEkf filter(config, chaser_est, cov);

    sensors::DeterministicRng rng(3003);
    const double sigma_range = 1.0; // 1 m range noise
    const double dt = 0.1; // 10 Hz
    const double total_duration = 60.0; // 60 s

    for (double t = 0.0; t <= total_duration; t += dt) {
        // Propagate truth target and chaser via RK4 two-body
        target_state = numerics::rk4_step(t, target_state, dt, [mu](double, const orbit::CartesianState& s) {
            return orbit::two_body_state_derivative(0.0, s, mu);
        });
        chaser_truth = numerics::rk4_step(t, chaser_truth, dt, [mu](double, const orbit::CartesianState& s) {
            return orbit::two_body_state_derivative(0.0, s, mu);
        });

        // Filter predict
        filter.predict(dt);

        // Range measurement with dropout between t=20s and t=35s
        const bool in_dropout = (t >= 20.0 && t <= 35.0);
        const double true_rho = (target_state.position - chaser_truth.position).norm();
        const double noise = rng.gaussian(0.0, sigma_range);
        const double meas_rho = true_rho + noise;

        double pred_rho = (target_state.position - filter.estimated_state().position).norm();
        double inno = 0.0;
        double nis = 0.0;

        if (!in_dropout) {
            filter.update_range(meas_rho, target_state.position, sigma_range);
            const auto& diag = filter.last_range_diagnostics();
            pred_rho = diag.predicted_range_m;
            inno = diag.innovation_m;
            nis = diag.normalized_innovation_squared;
        }

        const double pos_err = (chaser_truth.position - filter.estimated_state().position).norm();

        csv << t << ","
            << true_rho << ","
            << pred_rho << ","
            << (in_dropout ? 0.0 : meas_rho) << ","
            << inno << ","
            << nis << ","
            << pos_err << "\n";
    }

    std::cout << "Range telemetry written to " << output_path << std::endl;
}

void run_range_monte_carlo(const std::string& output_path) {
    std::ofstream csv(output_path);
    if (!csv.is_open()) {
        std::cerr << "Failed to open " << output_path << std::endl;
        return;
    }

    csv << "seed,position_error_before_m,position_error_after_m,range_innovation_m,range_NIS\n";
    csv << std::setprecision(8) << std::fixed;

    const std::size_t N_seeds = 100;
    const math::Vector3 target_pos{7000.0e3, 0.0, 0.0};
    const math::Vector3 nominal_sc{7000.0e3 - 500.0, 0.0, 0.0}; // 500 m separation
    const double sigma_r = 10.0;
    const double sigma_range = 1.5;

    for (std::uint64_t seed = 1; seed <= N_seeds; ++seed) {
        sensors::DeterministicRng rng(seed * 77);

        // Perturb truth position
        const double offset_x = rng.gaussian(0.0, sigma_r);
        const double offset_y = rng.gaussian(0.0, sigma_r);
        const double offset_z = rng.gaussian(0.0, sigma_r);
        const math::Vector3 truth_sc = nominal_sc + math::Vector3{offset_x, offset_y, offset_z};

        const double true_rho = (target_pos - truth_sc).norm();
        const double meas_rho = true_rho + rng.gaussian(0.0, sigma_range);

        math::ColVector<6> prior_state{};
        prior_state(0, 0) = nominal_sc.x();
        prior_state(1, 0) = nominal_sc.y();
        prior_state(2, 0) = nominal_sc.z();

        math::Matrix<6, 6> prior_cov;
        for (std::size_t i = 0; i < 6; ++i) {
            prior_cov(i, i) = (i < 3) ? sigma_r * sigma_r : 0.1;
        }

        const double err_before = (truth_sc - nominal_sc).norm();

        estimation::RangeUpdateDiagnostics diag;
        const auto res = estimation::execute_range_update(
            prior_state,
            prior_cov,
            meas_rho,
            target_pos,
            sigma_range,
            &diag
        );

        const math::Vector3 post_pos{res.posterior_state(0, 0), res.posterior_state(1, 0), res.posterior_state(2, 0)};
        const double err_after = (truth_sc - post_pos).norm();

        csv << seed << ","
            << err_before << ","
            << err_after << ","
            << diag.innovation_m << ","
            << diag.normalized_innovation_squared << "\n";
    }

    std::cout << "Range Monte Carlo written to " << output_path << std::endl;
}

} // namespace

int main() {
    std::cout << "Starting AstraDock M13C demonstration..." << std::endl;
    run_attitude_telemetry("data/m13c_attitude_telemetry.csv");
    run_attitude_monte_carlo("data/m13c_attitude_monte_carlo.csv");
    run_range_telemetry("data/m13c_range_telemetry.csv");
    run_range_monte_carlo("data/m13c_range_monte_carlo.csv");
    std::cout << "M13C demonstration finished successfully." << std::endl;
    return 0;
}
