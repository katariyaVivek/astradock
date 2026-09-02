// AstraDock M13B — IMU-aided translational EKF demonstration & validation driver.
//
// Produces the authoritative deterministic datasets consumed by the Python
// analysis layer (python/analysis/plot_imu_ekf.py) and the independent audit
// oracle (python/audit/independent_imu_ekf_reference.py):
//
//   data/m13b_imu_ekf_telemetry.csv     nominal multi-rate run (IMU 100 Hz,
//                                       GNSS 1 Hz) with a GNSS dropout window;
//                                       logged at 10 Hz plus every GNSS epoch
//   data/m13b_oracle_segment.csv        full-rate 60 s segment (all IMU
//                                       samples + all GNSS fixes) for the
//                                       independent Python EKF oracle
//   data/m13b_bias_sensitivity.csv      accelerometer-bias dead-reckoning study
//   data/m13b_imu_ekf_monte_carlo.csv   100-seed Monte Carlo comparing
//                                       GNSS-only vs IMU+GNSS filtering
//
// Scenario (identical to the unit-test baseline):
//   truth    : circular 500 km two-body orbit (RK4) with a slowly tumbling
//              body frame (pure quaternion kinematics)
//   sensors  : M12 models -- IMU 100 Hz (accel noise 1e-3 m/s^2, configured
//              bias knowledge), GNSS 1 Hz (sigma_r = 10 m, sigma_v = 0.05 m/s),
//              dropout window [600, 700] s in the nominal telemetry run
//   estimator: TranslationalEkf::predict_with_imu() dead reckoning between
//              GNSS fixes using measured specific force rotated BODY -> ECI
//              through an EXTERNALLY SUPPLIED reference attitude, gravity
//              added back at the estimated position.
//
// Known-attitude idealization (deliberate M13B staging): the navigation layer
// hands the filter an attitude propagated independently by pure quaternion
// kinematics from the same initial condition as truth. No SpacecraftState is
// ever passed to the filter; truth-dependent statistics (errors, NEES) are
// computed in this evaluation harness only.
//
// Monte Carlo note: each seed randomizes the sensor-noise realization AND the
// initial state error drawn from the initial covariance; attitude stays known
// per the M13B specification.

#include "attitude/rigid_body.hpp"
#include "dynamics/two_body.hpp"
#include "estimation/diagnostics.hpp"
#include "estimation/translational_ekf.hpp"
#include "math/constants.hpp"
#include "numerics/fixed_step_propagation.hpp"
#include "orbit/two_body_orbit.hpp"
#include "sensors/gnss.hpp"
#include "sensors/imu.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace astradock;

namespace {

constexpr double kMu = constants::earth_gravitational_parameter_m3_per_s2;
constexpr double kOrbitRadiusM = constants::earth_reference_radius_m + 500.0e3;
constexpr math::Vector3 kTumbleRadS{0.01, -0.02, 0.015};

constexpr double kImuDtS = 0.01;          // 100 Hz
constexpr double kGnssEverySteps = 100;   // 1 Hz
constexpr double kSigmaPositionM = 10.0;
constexpr double kSigmaVelocityMps = 0.05;
constexpr double kAccelNoiseStdMps2 = 1.0e-3;

constexpr double kTelemetryDurationS = 1200.0;
constexpr double kDropoutStartS = 600.0;
constexpr double kDropoutEndS = 700.0;

constexpr double kOracleSegmentDurationS = 60.0;

constexpr std::uint64_t kTelemetrySeed = 20250201u;
constexpr std::uint64_t kMonteCarloSeedBase = 31'000'000u;
constexpr int kMonteCarloRuns = 100;
constexpr double kMonteCarloDurationS = 300.0;

struct TelemetryRecord {
    double time_s{0.0};
    orbit::CartesianState truth{};
    orbit::CartesianState estimate{};         // post-update posterior
    orbit::CartesianState predicted{};        // post-predict pre-update prior
    bool gnss_valid{false};
    bool gnss_processed{false};
    math::Vector3 measured_position{};        // raw fix (oracle replay)
    math::Vector3 measured_velocity{};
    math::Vector3 innovation_position{};
    math::Vector3 innovation_velocity{};
    double nis{std::numeric_limits<double>::quiet_NaN()};
    math::Vector3 specific_force_body{};
    math::Vector3 angular_velocity_body{};
    math::Quaternion attitude_supplied{};     // navigation-layer input
    math::Vector3 sigma_position_m{};
    math::Vector3 sigma_velocity_mps{};
    double nees{0.0};
};

orbit::CartesianState make_initial_truth() {
    const double speed =
        orbit::circular_orbit_speed_m_per_s(kMu, kOrbitRadiusM);
    return {{kOrbitRadiusM, 0.0, 0.0}, {0.0, speed, 0.0}};
}

math::Quaternion advance_attitude_kinematics(
    const math::Quaternion& orientation,
    const math::Vector3& omega_body_rad_s,
    double dt_s) {
    const attitude::RotationalState next = numerics::rk4_step(
        0.0,
        attitude::RotationalState{orientation, omega_body_rad_s},
        dt_s,
        [](double, const attitude::RotationalState& s) {
            return attitude::RotationalState{
                attitude::quaternion_derivative(
                    s.orientation, s.angular_velocity_rad_per_s),
                math::Vector3{},
            };
        });
    return next.orientation.normalized();
}

estimation::TranslationalCovariance make_initial_covariance() {
    estimation::TranslationalCovariance p0 =
        estimation::TranslationalCovariance::zero();
    for (std::size_t i = 0; i < 3; ++i) {
        p0(i, i) = (25.0e3) * (25.0e3);
        p0(i + 3, i + 3) = 5.0 * 5.0;
    }
    return p0;
}

estimation::TranslationalEkfConfig make_filter_config(
    math::Vector3 configured_accel_bias) {
    estimation::TranslationalEkfConfig config;
    config.gravitational_parameter_m3_per_s2 = kMu;
    // The M13A dynamics-only path keeps its own baseline sigma; the IMU path
    // uses the accelerometer-driven sigma below (see Lesson 014).
    config.acceleration_noise_std_mps2 = 1.0e-3;
    config.imu.accelerometer_bias_body_mps2 = configured_accel_bias;
    config.imu.accelerometer_noise_std_mps2 = kAccelNoiseStdMps2;
    config.gnss_noise.position_variance_m2 = kSigmaPositionM * kSigmaPositionM;
    config.gnss_noise.velocity_variance_m2_per_s2 =
        kSigmaVelocityMps * kSigmaVelocityMps;
    return config;
}

orbit::CartesianState apply_initial_error(
    const orbit::CartesianState& truth,
    const math::Vector3& position_error,
    const math::Vector3& velocity_error) {
    return {truth.position + position_error, truth.velocity + velocity_error};
}

void capture_sigmas(const estimation::TranslationalEkf& filter,
                    math::Vector3& sigma_position,
                    math::Vector3& sigma_velocity) {
    const auto diagonal = filter.covariance().diagonal();
    sigma_position = math::Vector3{
        std::sqrt(diagonal(0, 0)),
        std::sqrt(diagonal(1, 0)),
        std::sqrt(diagonal(2, 0))};
    sigma_velocity = math::Vector3{
        std::sqrt(diagonal(3, 0)),
        std::sqrt(diagonal(4, 0)),
        std::sqrt(diagonal(5, 0))};
}

// ---------------------------------------------------------------------------
// Phase 1: nominal multi-rate telemetry run (with GNSS dropout window)
// ---------------------------------------------------------------------------

std::vector<TelemetryRecord> run_nominal(bool log_full_rate_segment_only) {
    const orbit::CartesianState initial_truth = make_initial_truth();

    // Deliberate initialization error; identical to M13A for comparability.
    const math::Vector3 initial_position_error{50.0e3, -30.0e3, 20.0e3};
    const math::Vector3 initial_velocity_error{5.0, -3.0, 2.0};

    sensors::ImuConfig imu_cfg;
    imu_cfg.sample_period_s = kImuDtS;
    imu_cfg.gyro_bias_rad_s = math::Vector3{5.0e-4, -3.0e-4, 2.0e-4};  // carried, unused
    imu_cfg.gyro_noise_std_rad_s = 1.0e-3;                             // carried, unused
    imu_cfg.accel_noise_std_mps2 = kAccelNoiseStdMps2;
    imu_cfg.random_seed = kTelemetrySeed;
    sensors::ImuSensor imu(imu_cfg);

    sensors::GnssConfig gnss_cfg;
    gnss_cfg.sample_period_s = 1.0;
    gnss_cfg.position_noise_std_m = kSigmaPositionM;
    gnss_cfg.velocity_noise_std_mps = kSigmaVelocityMps;
    if (!log_full_rate_segment_only) {
        gnss_cfg.dropouts = {{kDropoutStartS, kDropoutEndS}};
    }
    gnss_cfg.random_seed = kTelemetrySeed + 1'000'000u;
    sensors::GnssSensor gnss(gnss_cfg);

    estimation::TranslationalEkf filter(
        make_filter_config(math::Vector3{}),
        apply_initial_error(initial_truth, initial_position_error,
                            initial_velocity_error),
        make_initial_covariance());

    const auto truth_samples = numerics::propagate_fixed_step(
        0.0,
        log_full_rate_segment_only ? kOracleSegmentDurationS : kTelemetryDurationS,
        kImuDtS,
        initial_truth,
        numerics::IntegrationMethod::classical_rk4,
        [](double, const orbit::CartesianState& s) {
            return orbit::two_body_state_derivative(0.0, s, kMu);
        });

    math::Quaternion truth_attitude = math::Quaternion::identity();
    math::Quaternion reference_attitude = math::Quaternion::identity();

    std::vector<TelemetryRecord> records;
    records.reserve(truth_samples.size());

    spacecraft::SpacecraftState truth_sc{};
    truth_sc.translational = initial_truth;
    truth_sc.rotational = attitude::RotationalState{
        math::Quaternion::identity(), kTumbleRadS};

    for (std::size_t k = 0; k < truth_samples.size(); ++k) {
        const double t_s = truth_samples[k].time_s;
        if (k > 0) {
            truth_attitude = advance_attitude_kinematics(
                truth_attitude, kTumbleRadS, kImuDtS);
            reference_attitude = advance_attitude_kinematics(
                reference_attitude, kTumbleRadS, kImuDtS);
        }

        truth_sc.translational = truth_samples[k].state;
        truth_sc.rotational =
            attitude::RotationalState{truth_attitude, kTumbleRadS};

        TelemetryRecord record;
        record.time_s = t_s;
        record.truth = truth_samples[k].state;

        // --- IMU-aided prediction (dead reckoning across the interval) ----
        const sensors::ImuMeasurement imu_sample = imu.measure(t_s, truth_sc);
        record.specific_force_body = imu_sample.specific_force_body_mps2;
        record.angular_velocity_body = imu_sample.angular_velocity_body_rad_s;
        record.attitude_supplied = reference_attitude;

        if (k > 0) {
            estimation::ImuPredictionInput input;
            input.dt_s = t_s - truth_samples[k - 1].time_s;  // actual interval
            input.measured_specific_force_body_mps2 =
                imu_sample.specific_force_body_mps2;
            input.measured_angular_velocity_body_rad_s =
                imu_sample.angular_velocity_body_rad_s;
            input.attitude = estimation::NavigationAttitudeEstimate{
                reference_attitude};
            static_cast<void>(filter.predict_with_imu(input));
        }
        record.predicted = filter.estimated_state();  // post-predict pre-update

        // --- GNSS update ---------------------------------------------------
        if (static_cast<std::uint64_t>(k) % 100u == 0u) {
            const sensors::GnssMeasurement fix = gnss.measure(t_s, truth_sc);
            record.gnss_valid = fix.valid;
            if (fix.valid) {
                static_cast<void>(filter.update_gnss(
                    fix.position_eci_m, fix.velocity_eci_mps));
                record.gnss_processed = true;
                record.measured_position = fix.position_eci_m;
                record.measured_velocity = fix.velocity_eci_mps;
                const auto& diagnostics = filter.last_update_diagnostics();
                record.innovation_position =
                    diagnostics.innovation_position_eci_m;
                record.innovation_velocity =
                    diagnostics.innovation_velocity_eci_mps;
                record.nis = diagnostics.normalized_innovation_squared;
            }
        }

        record.estimate = filter.estimated_state();
        capture_sigmas(filter, record.sigma_position_m, record.sigma_velocity_mps);
        record.nees = estimation::translational_nees(
            record.truth, record.estimate, filter.covariance());
        records.push_back(record);
    }
    return records;
}

void write_telemetry_row(std::ofstream& file, const TelemetryRecord& r) {
    const double position_error =
        (r.truth.position - r.estimate.position).norm();
    const double velocity_error =
        (r.truth.velocity - r.estimate.velocity).norm();
    file << r.time_s << ','
         << r.truth.position.x() << ',' << r.truth.position.y() << ','
         << r.truth.position.z() << ','
         << r.estimate.position.x() << ',' << r.estimate.position.y() << ','
         << r.estimate.position.z() << ',' << position_error << ','
         << r.truth.velocity.x() << ',' << r.truth.velocity.y() << ','
         << r.truth.velocity.z() << ','
         << r.estimate.velocity.x() << ',' << r.estimate.velocity.y() << ','
         << r.estimate.velocity.z() << ',' << velocity_error << ','
         << r.specific_force_body.x() << ',' << r.specific_force_body.y()
         << ',' << r.specific_force_body.z() << ','
         << r.angular_velocity_body.x() << ',' << r.angular_velocity_body.y()
         << ',' << r.angular_velocity_body.z() << ','
         << r.predicted.position.x() << ',' << r.predicted.position.y()
         << ',' << r.predicted.position.z() << ','
         << r.predicted.velocity.x() << ',' << r.predicted.velocity.y()
         << ',' << r.predicted.velocity.z() << ',';
        if (r.gnss_processed) {
            file << r.measured_position.x() << ',' << r.measured_position.y()
                 << ',' << r.measured_position.z() << ',';
        } else {
            file << ",,,";
        }
        file << (r.gnss_valid ? 1 : 0) << ',';
    if (r.gnss_processed) {
        file << (r.innovation_position.norm()) << ',' << r.nis << ',';
    } else {
        file << ",,";
    }
    file << r.sigma_position_m.x() << ',' << r.sigma_position_m.y() << ','
         << r.sigma_position_m.z() << ','
         << r.sigma_velocity_mps.x() << ',' << r.sigma_velocity_mps.y()
         << ',' << r.sigma_velocity_mps.z() << ',' << r.nees << '\n';
}

// Nominal telemetry: decimate the 100 Hz stream to 10 Hz (every 10th IMU
// step); all GNSS epochs are multiples of 10 steps, so no fix is lost.
void write_decimated_telemetry(const std::vector<TelemetryRecord>& records,
                               const char* path) {
    std::ofstream file(path);
    file << std::setprecision(17);
    file << "time_s,"
         << "truth_position_eci_x_m,truth_position_eci_y_m,truth_position_eci_z_m,"
         << "estimated_position_eci_x_m,estimated_position_eci_y_m,estimated_position_eci_z_m,"
         << "position_error_m,"
         << "truth_velocity_eci_x_mps,truth_velocity_eci_y_mps,truth_velocity_eci_z_mps,"
         << "estimated_velocity_eci_x_mps,estimated_velocity_eci_y_mps,estimated_velocity_eci_z_mps,"
         << "velocity_error_mps,"
         << "imu_specific_force_body_x_mps2,imu_specific_force_body_y_mps2,imu_specific_force_body_z_mps2,"
         << "imu_angular_velocity_body_x_rad_s,imu_angular_velocity_body_y_rad_s,imu_angular_velocity_body_z_rad_s,"
         << "predicted_position_eci_x_m,predicted_position_eci_y_m,predicted_position_eci_z_m,"
         << "predicted_velocity_eci_x_mps,predicted_velocity_eci_y_mps,predicted_velocity_eci_z_mps,"
         << "meas_pos_x_m,meas_pos_y_m,meas_pos_z_m,"
         << "gnss_valid,innovation_norm_m,nis,"
         << "sigma_position_x_m,sigma_position_y_m,sigma_position_z_m,"
         << "sigma_velocity_x_mps,sigma_velocity_y_mps,sigma_velocity_z_mps,"
         << "nees\n";

    for (std::size_t k = 0; k < records.size(); ++k) {
        if (k % 10u != 0u) {
            continue;
        }
        write_telemetry_row(file, records[k]);
    }
}

// Full-rate oracle segment: EVERY IMU sample with the exact inputs the filter
// consumed (specific force, supplied attitude quaternion, actual dt implied
// by timestamps) and every GNSS fix, so the independent Python oracle can
// replay the identical estimation problem from scratch.
void write_full_rate_segment(const std::vector<TelemetryRecord>& records,
                             const char* path) {
    std::ofstream file(path);
    file << std::setprecision(17);
    file << "time_s,"
         << "imu_specific_force_body_x_mps2,imu_specific_force_body_y_mps2,imu_specific_force_body_z_mps2,"
         << "attitude_w,attitude_x,attitude_y,attitude_z,"
         << "gnss_valid,gnss_pos_x_m,gnss_pos_y_m,gnss_pos_z_m,"
         << "gnss_vel_x_mps,gnss_vel_y_mps,gnss_vel_z_mps,"
         << "est_pos_x_m,est_pos_y_m,est_pos_z_m,"
         << "est_vel_x_mps,est_vel_y_mps,est_vel_z_mps\n";

    for (const TelemetryRecord& r : records) {
        file << r.time_s << ','
             << r.specific_force_body.x() << ',' << r.specific_force_body.y()
             << ',' << r.specific_force_body.z() << ','
             << r.attitude_supplied.w() << ',' << r.attitude_supplied.x()
             << ',' << r.attitude_supplied.y() << ','
             << r.attitude_supplied.z() << ','
             << (r.gnss_processed ? 1 : 0) << ',';
        if (r.gnss_processed) {
            file << r.measured_position.x() << ',' << r.measured_position.y()
                 << ',' << r.measured_position.z() << ','
                 << r.measured_velocity.x() << ',' << r.measured_velocity.y()
                 << ',' << r.measured_velocity.z() << ',';
        } else {
            file << ",,,,,,";
        }
        file << r.estimate.position.x() << ',' << r.estimate.position.y()
             << ',' << r.estimate.position.z() << ','
             << r.estimate.velocity.x() << ',' << r.estimate.velocity.y()
             << ',' << r.estimate.velocity.z() << '\n';
    }
}

void report_nominal(const std::vector<TelemetryRecord>& records) {
    const std::size_t n = records.size();
    const double initial_error =
        math::Vector3{50.0e3, -30.0e3, 20.0e3}.norm();

    double ekf_pos_sq = 0.0;
    double raw_pos_sq = 0.0;
    double ekf_vel_sq = 0.0;
    double raw_vel_sq = 0.0;
    double nis_sum = 0.0;
    double nees_sum = 0.0;
    std::size_t updates = 0;
    double max_dropout_pos_err = 0.0;
    double max_dropout_vel_err = 0.0;
    double sigma_pos_before_dropout = 0.0;
    double sigma_pos_late_dropout = 0.0;
    double sigma_pos_after_recovery = 0.0;

    for (std::size_t i = n / 2; i < n; ++i) {
        const TelemetryRecord& r = records[i];
        if (r.time_s > kDropoutStartS && r.time_s < kDropoutEndS) {
            const double err = (r.truth.position - r.estimate.position).norm();
            max_dropout_pos_err = std::max(max_dropout_pos_err, err);
            max_dropout_vel_err =
                std::max(max_dropout_vel_err,
                         (r.truth.velocity - r.estimate.velocity).norm());
        }
        if (!r.gnss_processed) {
            continue;
        }
        const double pos_err = (r.truth.position - r.estimate.position).norm();
        const double vel_err =
            (r.truth.velocity - r.estimate.velocity).norm();
        ekf_pos_sq += pos_err * pos_err;
        ekf_vel_sq += vel_err * vel_err;
        // Raw-sensor error: measurement minus truth at the same epoch.
        raw_pos_sq += (r.measured_position - r.truth.position).norm();
        raw_vel_sq += (r.measured_velocity - r.truth.velocity).norm();
        nis_sum += r.nis;
        nees_sum += r.nees;
        ++updates;
    }

    for (const TelemetryRecord& r : records) {
        if (std::abs(r.time_s - (kDropoutStartS - kImuDtS)) < 1.0e-9) {
            sigma_pos_before_dropout = r.sigma_position_m.norm();
        }
        if (std::abs(r.time_s - kDropoutEndS) < 1.0e-9) {
            sigma_pos_late_dropout = r.sigma_position_m.norm();
        }
        if (std::abs(r.time_s - (kDropoutEndS + 200.0)) < 1.0e-9) {
            sigma_pos_after_recovery = r.sigma_position_m.norm();
        }
    }

    const double denom = static_cast<double>(updates);
    std::cout << "Nominal multi-rate run (" << kTelemetryDurationS
              << " s, IMU 100 Hz / GNSS 1 Hz, dropout [" << kDropoutStartS
              << ", " << kDropoutEndS << "] s):\n";
    std::cout << "  Initial position error:         " << initial_error
              << " m\n";
    std::cout << "  Final position error:           "
              << (records[n - 1].truth.position
                  - records[n - 1].estimate.position).norm()
              << " m\n";
    std::cout << "  Final velocity error:           "
              << (records[n - 1].truth.velocity
                  - records[n - 1].estimate.velocity).norm()
              << " m/s\n";
    std::cout << "  EKF position RMSE (2nd half):   "
              << std::sqrt(ekf_pos_sq / denom) << " m\n";
    std::cout << "  Raw GNSS position RMSE:         "
              << std::sqrt(raw_pos_sq / denom) << " m\n";
    std::cout << "  EKF velocity RMSE (2nd half):   "
              << std::sqrt(ekf_vel_sq / denom) << " m/s\n";
    std::cout << "  Raw GNSS velocity RMSE:         "
              << std::sqrt(raw_vel_sq / denom) << " m/s\n";
    std::cout << "  Mean NIS (valid updates):       " << nis_sum / denom
              << "\n";
    std::cout << "  Mean NEES (2nd half):           " << nees_sum / denom
              << "\n";
    std::cout << "  Max pos error during dropout:   "
              << max_dropout_pos_err << " m\n";
    std::cout << "  Max vel error during dropout:   "
              << max_dropout_vel_err << " m/s\n";
    std::cout << "  Position |sigma| before/late/recovery: "
              << sigma_pos_before_dropout << " -> " << sigma_pos_late_dropout
              << " -> " << sigma_pos_after_recovery << " m\n";
}

// ---------------------------------------------------------------------------
// Phase 2: accelerometer-bias sensitivity study (prediction only)
// ---------------------------------------------------------------------------

struct BiasCaseResult {
    double bias_mps2{0.0};
    double final_position_drift_m{0.0};
    double final_velocity_drift_mps{0.0};
};

BiasCaseResult run_bias_case(double true_bias_x_mps2) {
    const orbit::CartesianState initial_truth = make_initial_truth();

    // Identity-attitude configuration: keeps the bias vector fixed in ECI so
    // the measured drift can be compared against the classical flat-space
    // dead-reckoning law 0.5*b*t^2 without attitude-frame rotation sweeping
    // the force direction. Horizon 300 s keeps gravity-gradient curvature
    // modest.
    constexpr double kBiasDurationS = 300.0;
    sensors::ImuConfig imu_cfg;
    imu_cfg.sample_period_s = kImuDtS;
    imu_cfg.accel_bias_mps2 = math::Vector3{true_bias_x_mps2, 0.0, 0.0};
    imu_cfg.accel_noise_std_mps2 = 0.0;  // isolate the deterministic drift
    imu_cfg.random_seed = kTelemetrySeed + 2'000'000u;
    sensors::ImuSensor imu(imu_cfg);

    // Bias UNCALIBRATED in the filter: configuration says zero. This is the
    // deliberate mismatch motivating bias states in M13C/D.
    estimation::TranslationalEkf filter(
        make_filter_config(math::Vector3{}), initial_truth,
        estimation::TranslationalCovariance::zero());

    const auto truth_samples = numerics::propagate_fixed_step(
        0.0, kBiasDurationS, kImuDtS, initial_truth,
        numerics::IntegrationMethod::classical_rk4,
        [](double, const orbit::CartesianState& s) {
            return orbit::two_body_state_derivative(0.0, s, kMu);
        });

    spacecraft::SpacecraftState truth_sc{};
    truth_sc.translational = initial_truth;
    truth_sc.rotational = attitude::RotationalState{
        math::Quaternion::identity(), math::Vector3{}};

    for (std::size_t k = 0; k < truth_samples.size(); ++k) {
        const double t_s = truth_samples[k].time_s;
        truth_sc.translational = truth_samples[k].state;

        const sensors::ImuMeasurement sample = imu.measure(t_s, truth_sc);
        if (k > 0) {
            static_cast<void>(filter.predict_with_imu(
                {t_s - truth_samples[k - 1].time_s,
                 sample.specific_force_body_mps2,
                 sample.angular_velocity_body_rad_s,
                 estimation::NavigationAttitudeEstimate{math::Quaternion::identity()}}));
        }
    }

    const orbit::CartesianState final_estimate = filter.estimated_state();
    const orbit::CartesianState final_truth = truth_samples.back().state;
    return {true_bias_x_mps2,
            (final_truth.position - final_estimate.position).norm(),
            (final_truth.velocity - final_estimate.velocity).norm()};
}

void run_bias_study() {
    std::ofstream file("data/m13b_bias_sensitivity.csv");
    file << std::setprecision(17);
    file << "bias_x_body_mps2,duration_s,"
         << "final_position_drift_m,final_velocity_drift_mps,"
         << "analytic_position_half_b_t2_m,analytic_velocity_b_t_mps\n";

    constexpr double kBiasDurationS = 300.0;
    const std::vector<double> biases{0.0, 1.0e-4, 5.0e-4};
    std::cout << "\nAccelerometer-bias sensitivity (prediction-only, identity"
                 " attitude, " << kBiasDurationS << " s):\n";
    std::cout << "  bias [m/s^2]   pos drift [m]    vel drift [m/s]   "
                 "0.5*b*t^2 [m]\n";

    for (const double bias : biases) {
        const BiasCaseResult result = run_bias_case(bias);
        const double analytic_pos = 0.5 * bias * kBiasDurationS
                                    * kBiasDurationS;
        const double analytic_vel = bias * kBiasDurationS;
        file << bias << ',' << kBiasDurationS << ','
             << result.final_position_drift_m << ','
             << result.final_velocity_drift_mps << ',' << analytic_pos << ','
             << analytic_vel << '\n';
        std::cout << "  " << std::scientific << std::setprecision(3) << bias
                  << "      " << std::fixed << std::setprecision(3)
                  << result.final_position_drift_m << "          "
                  << result.final_velocity_drift_mps << "        "
                  << analytic_pos << '\n';
    }
    std::cout << "  Position drift grows linearly with bias and quadratically"
                 " with time:\n  an unestimated constant bias integrates into"
                 " velocity (b*t) then position (b*t^2/2).\n";
}

// ---------------------------------------------------------------------------
// Phase 3: Monte Carlo — GNSS-only vs IMU+GNSS
// ---------------------------------------------------------------------------

struct MonteCarloSummary {
    double ekf_pos_rmse_m{0.0};
    double raw_pos_rmse_m{0.0};
    double ekf_vel_rmse_mps{0.0};
    double raw_vel_rmse_mps{0.0};
    double mean_nis_last_half{0.0};
    double mean_nees_last_half{0.0};
    double final_nees{0.0};
    double final_pos_err_m{0.0};
    double final_vel_err_mps{0.0};
};

MonteCarloSummary summarize(const std::vector<TelemetryRecord>& records,
                            bool imu_case) {
    const std::size_t n = records.size();
    const std::size_t half = n / 2;
    MonteCarloSummary summary;

    double ekf_pos_sq = 0.0;
    double raw_pos_sq = 0.0;
    double ekf_vel_sq = 0.0;
    double raw_vel_sq = 0.0;
    double nis_sum = 0.0;
    double nees_sum = 0.0;
    std::size_t samples = 0;
    std::size_t updates = 0;

    for (std::size_t i = half; i < n; ++i) {
        const TelemetryRecord& r = records[i];
        const double pos_err = (r.truth.position - r.estimate.position).norm();
        const double vel_err =
            (r.truth.velocity - r.estimate.velocity).norm();
        ekf_pos_sq += pos_err * pos_err;
        ekf_vel_sq += vel_err * vel_err;
        nees_sum += r.nees;
        ++samples;
        if (r.gnss_processed) {
            raw_pos_sq += r.innovation_position.norm();
            raw_vel_sq += r.innovation_velocity.norm();
            nis_sum += r.nis;
            ++updates;
        }
    }

    summary.ekf_pos_rmse_m = std::sqrt(ekf_pos_sq / static_cast<double>(samples));
    summary.raw_pos_rmse_m =
        updates > 0 ? std::sqrt(raw_pos_sq / static_cast<double>(updates)) : 0.0;
    summary.ekf_vel_rmse_mps =
        std::sqrt(ekf_vel_sq / static_cast<double>(samples));
    summary.raw_vel_rmse_mps =
        updates > 0 ? std::sqrt(raw_vel_sq / static_cast<double>(updates)) : 0.0;
    summary.mean_nis_last_half =
        updates > 0 ? nis_sum / static_cast<double>(updates) : 0.0;
    summary.mean_nees_last_half =
        nees_sum / static_cast<double>(samples);
    summary.final_nees = records[n - 1].nees;
    summary.final_pos_err_m =
        (records[n - 1].truth.position - records[n - 1].estimate.position)
            .norm();
    summary.final_vel_err_mps =
        (records[n - 1].truth.velocity - records[n - 1].estimate.velocity)
            .norm();
    static_cast<void>(imu_case);
    return summary;
}

// Runs one Monte Carlo pass producing BOTH cases over the SAME truth and the
// SAME GNSS realization (both filters consume identical fixes).
std::pair<MonteCarloSummary, MonteCarloSummary> run_monte_carlo_pass(
    std::uint64_t seed) {
    const orbit::CartesianState initial_truth = make_initial_truth();

    // Randomized initialization error (drawn from P0-scale Gaussians).
    sensors::DeterministicRng init_rng(seed + 2'000'000u);
    const math::Vector3 pos_error{
        init_rng.gaussian(0.0, 25.0e3),
        init_rng.gaussian(0.0, 25.0e3),
        init_rng.gaussian(0.0, 25.0e3)};
    const math::Vector3 vel_error{
        init_rng.gaussian(0.0, 5.0),
        init_rng.gaussian(0.0, 5.0),
        init_rng.gaussian(0.0, 5.0)};

    sensors::ImuConfig imu_cfg;
    imu_cfg.sample_period_s = kImuDtS;
    imu_cfg.accel_noise_std_mps2 = kAccelNoiseStdMps2;
    imu_cfg.random_seed = seed;
    sensors::ImuSensor imu(imu_cfg);

    sensors::GnssConfig gnss_cfg;
    gnss_cfg.sample_period_s = 1.0;
    gnss_cfg.position_noise_std_m = kSigmaPositionM;
    gnss_cfg.velocity_noise_std_mps = kSigmaVelocityMps;
    gnss_cfg.random_seed = seed + 1'000'000u;
    sensors::GnssSensor gnss(gnss_cfg);  // one stream shared by both cases

    estimation::TranslationalEkf imu_filter(
        make_filter_config(math::Vector3{}),
        apply_initial_error(initial_truth, pos_error, vel_error),
        make_initial_covariance());

    // GNSS-only case: dynamics-only M13A prediction between fixes.
    estimation::TranslationalEkfConfig gnss_only_config = make_filter_config({});
    estimation::TranslationalEkf gnss_filter(
        gnss_only_config,
        apply_initial_error(initial_truth, pos_error, vel_error),
        make_initial_covariance());

    const auto truth_samples = numerics::propagate_fixed_step(
        0.0, kMonteCarloDurationS, kImuDtS, initial_truth,
        numerics::IntegrationMethod::classical_rk4,
        [](double, const orbit::CartesianState& s) {
            return orbit::two_body_state_derivative(0.0, s, kMu);
        });

    std::vector<TelemetryRecord> imu_records;
    std::vector<TelemetryRecord> gnss_records;
    imu_records.reserve(truth_samples.size());
    gnss_records.reserve(truth_samples.size());

    math::Quaternion reference_attitude = math::Quaternion::identity();
    spacecraft::SpacecraftState truth_sc{};
    truth_sc.translational = initial_truth;
    truth_sc.rotational = attitude::RotationalState{
        math::Quaternion::identity(), kTumbleRadS};

    for (std::size_t k = 0; k < truth_samples.size(); ++k) {
        const double t_s = truth_samples[k].time_s;
        if (k > 0) {
            reference_attitude = advance_attitude_kinematics(
                reference_attitude, kTumbleRadS, kImuDtS);
            truth_sc.rotational =
                attitude::RotationalState{advance_attitude_kinematics(
                                              truth_sc.rotational.orientation,
                                              kTumbleRadS, kImuDtS),
                                          kTumbleRadS};
        }
        truth_sc.translational = truth_samples[k].state;

        const sensors::ImuMeasurement imu_sample = imu.measure(t_s, truth_sc);

        TelemetryRecord imu_record;
        TelemetryRecord gnss_record;
        imu_record.time_s = t_s;
        gnss_record.time_s = t_s;
        imu_record.truth = truth_samples[k].state;
        gnss_record.truth = truth_samples[k].state;
        imu_record.specific_force_body = imu_sample.specific_force_body_mps2;
        imu_record.angular_velocity_body =
            imu_sample.angular_velocity_body_rad_s;

        if (k > 0) {
            const double dt_s = t_s - truth_samples[k - 1].time_s;
            static_cast<void>(imu_filter.predict_with_imu(
                {dt_s,
                 imu_sample.specific_force_body_mps2,
                 imu_sample.angular_velocity_body_rad_s,
                 estimation::NavigationAttitudeEstimate{reference_attitude}}));
            gnss_filter.predict(dt_s);
        }
        imu_record.predicted = imu_filter.estimated_state();
        gnss_record.predicted = gnss_filter.estimated_state();

        if (static_cast<std::uint64_t>(k) % 100u == 0u) {
            // ONE measurement realization feeds both filters.
            const sensors::GnssMeasurement fix = gnss.measure(t_s, truth_sc);
            imu_record.gnss_valid = fix.valid;
            gnss_record.gnss_valid = fix.valid;
            if (fix.valid) {
                static_cast<void>(imu_filter.update_gnss(
                    fix.position_eci_m, fix.velocity_eci_mps));
                static_cast<void>(gnss_filter.update_gnss(
                    fix.position_eci_m, fix.velocity_eci_mps));
                imu_record.gnss_processed = true;
                gnss_record.gnss_processed = true;
                imu_record.innovation_position =
                    imu_filter.last_update_diagnostics().innovation_position_eci_m;
                imu_record.innovation_velocity =
                    imu_filter.last_update_diagnostics().innovation_velocity_eci_mps;
                imu_record.nis =
                    imu_filter.last_update_diagnostics().normalized_innovation_squared;
                gnss_record.innovation_position =
                    gnss_filter.last_update_diagnostics().innovation_position_eci_m;
                gnss_record.innovation_velocity =
                    gnss_filter.last_update_diagnostics().innovation_velocity_eci_mps;
                gnss_record.nis =
                    gnss_filter.last_update_diagnostics().normalized_innovation_squared;
            }
        }

        imu_record.estimate = imu_filter.estimated_state();
        gnss_record.estimate = gnss_filter.estimated_state();
        imu_record.nees = estimation::translational_nees(
            imu_record.truth, imu_record.estimate, imu_filter.covariance());
        gnss_record.nees = estimation::translational_nees(
            gnss_record.truth, gnss_record.estimate, gnss_filter.covariance());
        imu_records.push_back(imu_record);
        gnss_records.push_back(gnss_record);
    }

    return {summarize(gnss_records, false), summarize(imu_records, true)};
}

void run_monte_carlo() {
    std::ofstream file("data/m13b_imu_ekf_monte_carlo.csv");
    file << std::setprecision(17);
    file << "seed,case,"
         << "ekf_pos_rmse_m,raw_pos_rmse_m,"
         << "ekf_vel_rmse_mps,raw_vel_rmse_mps,"
         << "mean_nis_last_half,mean_nees_last_half,"
         << "final_nees,final_pos_err_m,final_vel_err_mps\n";

    double sum[2][6] = {};
    int beats_raw[2] = {0, 0};
    std::vector<double> final_nees_values;

    for (int run = 0; run < kMonteCarloRuns; ++run) {
        const auto seed =
            kMonteCarloSeedBase + static_cast<std::uint64_t>(run);
        const auto [gnss_summary, imu_summary] = run_monte_carlo_pass(seed);

        file << seed << ",gnss_only," << gnss_summary.ekf_pos_rmse_m << ','
             << gnss_summary.raw_pos_rmse_m << ',' << gnss_summary.ekf_vel_rmse_mps
             << ',' << gnss_summary.raw_vel_rmse_mps << ','
             << gnss_summary.mean_nis_last_half << ','
             << gnss_summary.mean_nees_last_half << ',' << gnss_summary.final_nees
             << ',' << gnss_summary.final_pos_err_m << ','
             << gnss_summary.final_vel_err_mps << '\n';
        file << seed << ",imu_gnss," << imu_summary.ekf_pos_rmse_m << ','
             << imu_summary.raw_pos_rmse_m << ',' << imu_summary.ekf_vel_rmse_mps
             << ',' << imu_summary.raw_vel_rmse_mps << ','
             << imu_summary.mean_nis_last_half << ','
             << imu_summary.mean_nees_last_half << ',' << imu_summary.final_nees
             << ',' << imu_summary.final_pos_err_m << ','
             << imu_summary.final_vel_err_mps << '\n';

        const MonteCarloSummary* summaries[2] = {&gnss_summary, &imu_summary};
        for (int c = 0; c < 2; ++c) {
            sum[c][0] += summaries[c]->ekf_pos_rmse_m;
            sum[c][1] += summaries[c]->raw_pos_rmse_m;
            sum[c][2] += summaries[c]->ekf_vel_rmse_mps;
            sum[c][3] += summaries[c]->raw_vel_rmse_mps;
            sum[c][4] += summaries[c]->mean_nis_last_half;
            sum[c][5] += summaries[c]->mean_nees_last_half;
            if (summaries[c]->ekf_pos_rmse_m < summaries[c]->raw_pos_rmse_m) {
                ++beats_raw[c];
            }
        }
        final_nees_values.push_back(imu_summary.final_nees);
    }

    const double runs = static_cast<double>(kMonteCarloRuns);
    std::sort(final_nees_values.begin(), final_nees_values.end());
    const double median_final_nees =
        final_nees_values.size() % 2 == 1
            ? final_nees_values[final_nees_values.size() / 2]
            : 0.5 * (final_nees_values[final_nees_values.size() / 2 - 1]
                     + final_nees_values[final_nees_values.size() / 2]);

    std::cout << "\nMonte Carlo (" << kMonteCarloRuns << " seeds x "
              << kMonteCarloDurationS << " s, identical fixes per seed):\n";
    std::cout << "                          GNSS-only      IMU+GNSS\n";
    std::cout << "  Mean EKF pos RMSE [m]:  " << std::setw(10) << sum[0][0] / runs
              << std::setw(14) << sum[1][0] / runs << '\n';
    std::cout << "  Mean raw pos RMSE [m]:  " << std::setw(10) << sum[0][1] / runs
              << std::setw(14) << sum[1][1] / runs << '\n';
    std::cout << "  Mean EKF vel RMSE[m/s]: " << std::setw(10) << sum[0][2] / runs
              << std::setw(14) << sum[1][2] / runs << '\n';
    std::cout << "  Mean raw vel RMSE[m/s]: " << std::setw(10) << sum[0][3] / runs
              << std::setw(14) << sum[1][3] / runs << '\n';
    std::cout << "  Mean converged NIS:     " << std::setw(10) << sum[0][4] / runs
              << std::setw(14) << sum[1][4] / runs << '\n';
    std::cout << "  Mean NEES (2nd half):   " << std::setw(10) << sum[0][5] / runs
              << std::setw(14) << sum[1][5] / runs << '\n';
    std::cout << "  Runs beating raw GNSS:  " << std::setw(8) << beats_raw[0]
              << "/" << kMonteCarloRuns << std::setw(8) << beats_raw[1] << '/'
              << kMonteCarloRuns << '\n';
    std::cout << "  Median final NEES (IMU+GNSS): " << median_final_nees
              << "  (chi-square(6): mean 6, 95% band [1.64, 14.45])\n";
}

}  // namespace

int main() {
    std::cout << "AstraDock M13B - IMU-Aided Translational EKF demonstration\n";
    std::cout << "=============================================================\n";
    std::cout << "Truth : circular 500 km two-body orbit, tumbling body frame\n";
    std::cout << "Sensors: IMU 100 Hz (sigma_f = " << kAccelNoiseStdMps2
              << " m/s^2), GNSS 1 Hz (sigma_r = " << kSigmaPositionM
              << " m, sigma_v = " << kSigmaVelocityMps << " m/s)\n";
    std::cout << "Filter: TranslationalEkf::predict_with_imu, known reference"
                 " attitude,\n        gravity compensation, configured bias"
                 " removal\n\n";

    // ---- Phase 1a: nominal telemetry run with dropout ----------------------
    std::cout << "Phase 1a: nominal telemetry run...\n";
    const std::vector<TelemetryRecord> nominal_records = run_nominal(false);
    write_decimated_telemetry(nominal_records, "data/m13b_imu_ekf_telemetry.csv");
    report_nominal(nominal_records);

    // ---- Phase 1b: full-rate oracle segment ---------------------------------
    std::cout << "\nPhase 1b: full-rate oracle segment (" << kOracleSegmentDurationS
              << " s, no dropout)...\n";
    const std::vector<TelemetryRecord> segment_records = run_nominal(true);
    write_full_rate_segment(segment_records, "data/m13b_oracle_segment.csv");
    std::cout << "  Wrote " << segment_records.size()
              << " full-rate records to data/m13b_oracle_segment.csv\n";

    // ---- Phase 2: bias sensitivity ------------------------------------------
    run_bias_study();

    // ---- Phase 3: Monte Carlo ------------------------------------------------
    run_monte_carlo();

    std::cout << "\nWrote data/m13b_imu_ekf_telemetry.csv, "
                 "data/m13b_oracle_segment.csv,\n"
                 "data/m13b_bias_sensitivity.csv, and "
                 "data/m13b_imu_ekf_monte_carlo.csv\n";
    return 0;
}
