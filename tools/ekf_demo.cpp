// AstraDock M13A — Extended Kalman Filter demonstration and validation driver.
//
// Produces the authoritative deterministic datasets consumed by the Python
// analysis layer (python/analysis/plot_estimation.py) and the independent
// audit oracle (python/audit/independent_ekf_reference.py):
//
//   data/m13_ekf_telemetry.csv     nominal run: truth vs measurements vs
//                                  estimate vs covariance vs innovations
//   data/m13_ekf_monte_carlo.csv   per-seed Monte Carlo summary statistics
//
// Scenario (identical to the unit-test baseline):
//   truth    : circular 500 km two-body orbit, RK4, dt = 1 s
//   sensor   : GNSS 1 Hz, sigma_r = 10 m, sigma_v = 0.05 m/s, isotropic,
//              dropout window [3000, 3900] s in the telemetry run
//   estimator: 6-state translational EKF initialized deliberately ~62 km /
//              6.2 m/s away from truth; the filter NEVER receives truth.
//              Truth-dependent statistics (NEES, estimation errors) are
//              computed in this evaluation harness only.
//
// Monte Carlo note: all runs share the same deterministic initialization
// offset; seed variation therefore covers sensor-noise realization only.
// Initialization-error dispersion is studied separately in the unit suite.

#include "estimation/diagnostics.hpp"
#include "estimation/translational_ekf.hpp"
#include "math/constants.hpp"
#include "numerics/fixed_step_propagation.hpp"
#include "orbit/two_body_orbit.hpp"
#include "sensors/sensor_common.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

using namespace astradock;

namespace {

constexpr double kMu = constants::earth_gravitational_parameter_m3_per_s2;
constexpr double kOrbitRadiusM =
    constants::earth_reference_radius_m + 500.0e3;
constexpr double kSigmaPositionM = 10.0;    // GNSS position noise, m
constexpr double kSigmaVelocityMps = 0.05;  // GNSS velocity noise, m/s
constexpr double kProcessNoiseStdMps2 = 1.0e-3;
constexpr double kDtS = 1.0;

constexpr double kTelemetryDurationS = 6000.0;
constexpr double kDropoutStartS = 3000.0;
constexpr double kDropoutEndS = 3900.0;

constexpr std::uint64_t kTelemetrySeed = 20250001u;
constexpr std::uint64_t kMonteCarloSeedBase = 90000000u;
constexpr int kMonteCarloRuns = 100;
constexpr double kMonteCarloDurationS = 2000.0;

struct RunRecord {
    orbit::CartesianState truth{};
    orbit::CartesianState estimate{};
    bool measured{false};
    math::Vector3 measured_position{};
    math::Vector3 measured_velocity{};
    math::Vector3 innovation_position{};
    math::Vector3 innovation_velocity{};
    double nis{std::numeric_limits<double>::quiet_NaN()};
    double nees{0.0};
    double sigma_pos_m{0.0};
    double sigma_vel_mps{0.0};
};

orbit::CartesianState make_initial_truth() {
    const double speed_m_per_s =
        orbit::circular_orbit_speed_m_per_s(kMu, kOrbitRadiusM);
    return {{kOrbitRadiusM, 0.0, 0.0}, {0.0, speed_m_per_s, 0.0}};
}

std::vector<RunRecord> run_filter(
    std::uint64_t seed,
    double duration_s,
    bool apply_dropout) {
    const orbit::CartesianState initial_truth = make_initial_truth();

    // Deliberate initialization error; the filter starts wrong on purpose.
    const math::Vector3 initial_position_error{50.0e3, -30.0e3, 20.0e3};
    const math::Vector3 initial_velocity_error{5.0, -3.0, 2.0};
    const orbit::CartesianState initial_estimate{
        initial_truth.position + initial_position_error,
        initial_truth.velocity + initial_velocity_error,
    };

    estimation::TranslationalCovariance p0 =
        estimation::TranslationalCovariance::zero();
    for (std::size_t i = 0; i < 3; ++i) {
        p0(i, i) = (25.0e3) * (25.0e3);
        p0(i + 3, i + 3) = 5.0 * 5.0;
    }

    estimation::TranslationalEkfConfig config;
    config.gravitational_parameter_m3_per_s2 = kMu;
    config.acceleration_noise_std_mps2 = kProcessNoiseStdMps2;
    config.gnss_noise.position_variance_m2 = kSigmaPositionM * kSigmaPositionM;
    config.gnss_noise.velocity_variance_m2_per_s2 =
        kSigmaVelocityMps * kSigmaVelocityMps;
    estimation::TranslationalEkf filter(config, initial_estimate, p0);

    sensors::DeterministicRng rng(seed);

    const auto truth_samples = numerics::propagate_fixed_step(
        0.0,
        duration_s,
        kDtS,
        initial_truth,
        numerics::IntegrationMethod::classical_rk4,
        [](double, const orbit::CartesianState& s) {
            return orbit::two_body_state_derivative(0.0, s, kMu);
        });

    std::vector<RunRecord> records;
    records.reserve(truth_samples.size());
    for (std::size_t k = 0; k < truth_samples.size(); ++k) {
        const double t_s = static_cast<double>(k);
        if (k > 0) {
            filter.predict(kDtS);
        }
        const orbit::CartesianState truth_k = truth_samples[k].state;
        const bool in_dropout =
            apply_dropout && t_s >= kDropoutStartS && t_s <= kDropoutEndS;

        RunRecord record;
        record.truth = truth_k;
        record.estimate = filter.estimated_state();

        if (!in_dropout) {
            record.measured = true;
            record.measured_position =
                truth_k.position + rng.gaussian_vector3(kSigmaPositionM);
            record.measured_velocity =
                truth_k.velocity + rng.gaussian_vector3(kSigmaVelocityMps);
            filter.update_gnss(record.measured_position, record.measured_velocity);
            record.estimate = filter.estimated_state();
            const auto& diagnostics = filter.last_update_diagnostics();
            record.innovation_position = diagnostics.innovation_position_eci_m;
            record.innovation_velocity = diagnostics.innovation_velocity_eci_mps;
            record.nis = diagnostics.normalized_innovation_squared;
        }

        record.nees = estimation::translational_nees(
            truth_k, record.estimate, filter.covariance());
        const auto diagonal = filter.covariance().diagonal();
        record.sigma_pos_m = std::sqrt(
            std::max({diagonal(0, 0), diagonal(1, 0), diagonal(2, 0)}));
        record.sigma_vel_mps = std::sqrt(
            std::max({diagonal(3, 0), diagonal(4, 0), diagonal(5, 0)}));
        records.push_back(record);
    }
    return records;
}

void write_telemetry_csv(const std::vector<RunRecord>& records) {
    std::ofstream file("data/m13_ekf_telemetry.csv");
    file << "time_s,"
         << "truth_pos_x_m,truth_pos_y_m,truth_pos_z_m,"
         << "truth_vel_x_mps,truth_vel_y_mps,truth_vel_z_mps,"
         << "meas_valid,meas_pos_x_m,meas_pos_y_m,meas_pos_z_m,"
         << "meas_vel_x_mps,meas_vel_y_mps,meas_vel_z_mps,"
         << "est_pos_x_m,est_pos_y_m,est_pos_z_m,"
         << "est_vel_x_mps,est_vel_y_mps,est_vel_z_mps,"
         << "position_error_m,velocity_error_mps,"
         << "sigma_pos_max_m,sigma_vel_max_mps,"
         << "innov_pos_x_m,innov_pos_y_m,innov_pos_z_m,"
         << "innov_vel_x_mps,innov_vel_y_mps,innov_vel_z_mps,"
         << "nis,nees\n";
    file << std::setprecision(17);

    for (std::size_t k = 0; k < records.size(); ++k) {
        const RunRecord& r = records[k];
        file << static_cast<double>(k) << ','
             << r.truth.position.x() << ',' << r.truth.position.y() << ','
             << r.truth.position.z() << ','
             << r.truth.velocity.x() << ',' << r.truth.velocity.y() << ','
             << r.truth.velocity.z() << ',';
        if (r.measured) {
            file << 1 << ','
                 << r.measured_position.x() << ',' << r.measured_position.y()
                 << ',' << r.measured_position.z() << ','
                 << r.measured_velocity.x() << ',' << r.measured_velocity.y()
                 << ',' << r.measured_velocity.z() << ',';
        } else {
            file << 0 << ",,,,,,,";
        }
        file << r.estimate.position.x() << ',' << r.estimate.position.y() << ','
             << r.estimate.position.z() << ','
             << r.estimate.velocity.x() << ',' << r.estimate.velocity.y() << ','
             << r.estimate.velocity.z() << ','
             << (r.truth.position - r.estimate.position).norm() << ','
             << (r.truth.velocity - r.estimate.velocity).norm() << ','
             << r.sigma_pos_m << ',' << r.sigma_vel_mps << ',';
        if (r.measured) {
            file << r.innovation_position.x() << ',' << r.innovation_position.y()
                 << ',' << r.innovation_position.z() << ','
                 << r.innovation_velocity.x() << ',' << r.innovation_velocity.y()
                 << ',' << r.innovation_velocity.z() << ','
                 << r.nis << ',';
        } else {
            file << ",,,,,,,";
        }
        file << r.nees << '\n';
    }
}

double median_of(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    const std::size_t n = values.size();
    if (n % 2 == 1) {
        return values[n / 2];
    }
    return 0.5 * (values[n / 2 - 1] + values[n / 2]);
}

void run_monte_carlo() {
    std::ofstream mc_file("data/m13_ekf_monte_carlo.csv");
    mc_file << "seed,"
            << "final_pos_err_m,final_vel_err_mps,final_nees,"
            << "ekf_pos_rmse_m,raw_pos_rmse_m,"
            << "ekf_vel_rmse_mps,raw_vel_rmse_mps,"
            << "mean_nis_last_half,mean_nees_last_half\n";
    mc_file << std::setprecision(17);

    std::vector<double> final_errors;
    std::vector<double> final_nees_values;
    int beats_raw = 0;
    double sum_ekf_pos_rmse = 0.0;
    double sum_raw_pos_rmse = 0.0;
    double sum_ekf_vel_rmse = 0.0;
    double sum_raw_vel_rmse = 0.0;
    double sum_mean_nis = 0.0;

    for (int run = 0; run < kMonteCarloRuns; ++run) {
        const auto seed = kMonteCarloSeedBase + static_cast<std::uint64_t>(run);
        const std::vector<RunRecord> records =
            run_filter(seed, kMonteCarloDurationS, false);
        const std::size_t n = records.size();
        const std::size_t half = n / 2;

        double ekf_pos_sq = 0.0;
        double raw_pos_sq = 0.0;
        double ekf_vel_sq = 0.0;
        double raw_vel_sq = 0.0;
        double nis_sum = 0.0;
        double nees_sum = 0.0;
        std::size_t updates = 0;

        for (std::size_t i = half; i < n; ++i) {
            const double pos_err =
                (records[i].truth.position - records[i].estimate.position).norm();
            const double vel_err =
                (records[i].truth.velocity - records[i].estimate.velocity).norm();
            ekf_pos_sq += pos_err * pos_err;
            ekf_vel_sq += vel_err * vel_err;
            nees_sum += records[i].nees;
            if (records[i].measured) {
                raw_pos_sq += (records[i].truth.position
                               - records[i].measured_position)
                    .norm();
                raw_vel_sq += (records[i].truth.velocity
                               - records[i].measured_velocity)
                    .norm();
                nis_sum += records[i].nis;
                ++updates;
            }
        }

        const double samples = static_cast<double>(n - half);
        const double ekf_pos_rmse = std::sqrt(ekf_pos_sq / samples);
        const double raw_pos_rmse = std::sqrt(raw_pos_sq / samples);
        const double ekf_vel_rmse = std::sqrt(ekf_vel_sq / samples);
        const double raw_vel_rmse = std::sqrt(raw_vel_sq / samples);
        const double mean_nis = nis_sum / static_cast<double>(updates);
        const double mean_nees = nees_sum / samples;

        const double final_pos_err =
            (records[n - 1].truth.position - records[n - 1].estimate.position)
                .norm();
        const double final_vel_err =
            (records[n - 1].truth.velocity - records[n - 1].estimate.velocity)
                .norm();
        const double final_nees = records[n - 1].nees;

        mc_file << seed << ',' << final_pos_err << ',' << final_vel_err << ','
                << final_nees << ',' << ekf_pos_rmse << ',' << raw_pos_rmse
                << ',' << ekf_vel_rmse << ',' << raw_vel_rmse << ',' << mean_nis
                << ',' << mean_nees << '\n';

        final_errors.push_back(final_pos_err);
        final_nees_values.push_back(final_nees);
        sum_ekf_pos_rmse += ekf_pos_rmse;
        sum_raw_pos_rmse += raw_pos_rmse;
        sum_ekf_vel_rmse += ekf_vel_rmse;
        sum_raw_vel_rmse += raw_vel_rmse;
        sum_mean_nis += mean_nis;
        if (ekf_pos_rmse < raw_pos_rmse) {
            ++beats_raw;
        }
    }

    std::cout << "Monte Carlo (" << kMonteCarloRuns << " seeds x "
              << kMonteCarloDurationS << " s):\n";
    std::cout << "  Median final position error:   "
              << median_of(final_errors) << " m\n";
    std::cout << "  Median final NEES:             "
              << median_of(final_nees_values) << "  (chi-square(6): mean 6, "
                                                 "95% band [1.64, 14.45])\n";
    std::cout << "  Mean converged NIS:            "
              << sum_mean_nis / static_cast<double>(kMonteCarloRuns) << "\n";
    std::cout << "  Mean EKF position RMSE:        "
              << sum_ekf_pos_rmse / static_cast<double>(kMonteCarloRuns)
              << " m\n";
    std::cout << "  Mean raw position RMSE:        "
              << sum_raw_pos_rmse / static_cast<double>(kMonteCarloRuns)
              << " m\n";
    std::cout << "  Mean EKF velocity RMSE:        "
              << sum_ekf_vel_rmse / static_cast<double>(kMonteCarloRuns)
              << " m/s\n";
    std::cout << "  Mean raw velocity RMSE:        "
              << sum_raw_vel_rmse / static_cast<double>(kMonteCarloRuns)
              << " m/s\n";
    std::cout << "  Runs beating raw GNSS (pos):   " << beats_raw << "/"
              << kMonteCarloRuns << "\n";
}

}  // namespace

int main() {
    std::cout << "AstraDock M13A - Translational EKF demonstration\n";
    std::cout << "================================================\n";
    std::cout << "Truth : circular 500 km two-body orbit (RK4, dt = " << kDtS
              << " s)\n";
    std::cout << "Sensor: GNSS 1 Hz, sigma_r = " << kSigmaPositionM
              << " m, sigma_v = " << kSigmaVelocityMps << " m/s\n";
    std::cout << "Filter: 6-state EKF, sigma_a = " << kProcessNoiseStdMps2
              << " m/s^2, deliberate ~62 km init error\n\n";

    // ---- Phase 1: nominal telemetry run with dropout window ---------------
    const std::vector<RunRecord> records =
        run_filter(kTelemetrySeed, kTelemetryDurationS, true);
    write_telemetry_csv(records);

    const std::size_t n = records.size();
    const double initial_error = math::Vector3{50.0e3, -30.0e3, 20.0e3}.norm();

    double ekf_pos_sq = 0.0;
    double raw_pos_sq = 0.0;
    double ekf_vel_sq = 0.0;
    double raw_vel_sq = 0.0;
    double nis_sum = 0.0;
    double nees_sum = 0.0;
    std::size_t updates = 0;
    double max_dropout_pos_err = 0.0;

    for (std::size_t i = n / 2; i < n; ++i) {
        const double t_s = static_cast<double>(i);
        if (t_s > kDropoutStartS && t_s < kDropoutEndS) {
            max_dropout_pos_err =
                std::max(max_dropout_pos_err,
                         (records[i].truth.position - records[i].estimate.position)
                             .norm());
        }
        if (!records[i].measured) {
            continue;
        }
        const double pos_err =
            (records[i].truth.position - records[i].estimate.position).norm();
        const double vel_err =
            (records[i].truth.velocity - records[i].estimate.velocity).norm();
        ekf_pos_sq += pos_err * pos_err;
        ekf_vel_sq += vel_err * vel_err;
        raw_pos_sq +=
            (records[i].truth.position - records[i].measured_position).norm();
        raw_vel_sq +=
            (records[i].truth.velocity - records[i].measured_velocity).norm();
        nis_sum += records[i].nis;
        nees_sum += records[i].nees;
        ++updates;
    }

    const double denom = static_cast<double>(updates);

    std::cout << "Nominal run (" << kTelemetryDurationS
              << " s, dropout [" << kDropoutStartS << ", " << kDropoutEndS
              << "] s):\n";
    std::cout << "  Initial position error:        " << initial_error << " m\n";
    std::cout << "  Final position error:          "
              << (records[n - 1].truth.position - records[n - 1].estimate.position)
                     .norm()
              << " m\n";
    std::cout << "  Final velocity error:          "
              << (records[n - 1].truth.velocity - records[n - 1].estimate.velocity)
                     .norm()
              << " m/s\n";
    std::cout << "  EKF position RMSE (2nd half):  "
              << std::sqrt(ekf_pos_sq / denom) << " m\n";
    std::cout << "  Raw GNSS position RMSE:        "
              << std::sqrt(raw_pos_sq / denom) << " m\n";
    std::cout << "  EKF velocity RMSE (2nd half):  "
              << std::sqrt(ekf_vel_sq / denom) << " m/s\n";
    std::cout << "  Raw GNSS velocity RMSE:        "
              << std::sqrt(raw_vel_sq / denom) << " m/s\n";
    std::cout << "  Mean NIS (valid updates):      " << nis_sum / denom << "\n";
    std::cout << "  Mean NEES (2nd half):          " << nees_sum / denom << "\n";
    std::cout << "  Max pos error during dropout:  " << max_dropout_pos_err
              << " m\n";
    std::cout << "  Position sigma before dropout: "
              << records[static_cast<std::size_t>(kDropoutStartS) - 1]
                     .sigma_pos_m
              << " m -> late dropout: "
              << records[static_cast<std::size_t>(kDropoutEndS)].sigma_pos_m
              << " m\n";

    // ---- Phase 2: Monte Carlo study ---------------------------------------
    std::cout << '\n';
    run_monte_carlo();

    std::cout << "\nWrote data/m13_ekf_telemetry.csv and "
                 "data/m13_ekf_monte_carlo.csv\n";
    return 0;
}
