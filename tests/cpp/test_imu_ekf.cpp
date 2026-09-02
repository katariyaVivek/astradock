// AstraDock M13B — IMU-aided translational EKF propagation test suite.
//
// Layers verified here (per the M13 independent verification standard):
//   Layer 1: analytical specific-force frame transformations, known-attitude
//            acceleration reconstruction, first-order prediction limits,
//            hand-built covariance propagation, bias drift scaling laws,
//            finite-difference audits of the IMU-aided dynamics Jacobian
//   Layer 2: free-fall and drag multi-rate simulation loops, GNSS dropout
//            growth/recovery, attitude-dependency diagnostics, timestamp-
//            jitter handling, determinism, truth non-interference
//   Layer 3 lives in python/audit/independent_imu_ekf_reference.py
//
// Known-attitude idealization: every scenario supplies the estimator with an
// externally propagated reference attitude (quaternion kinematics only). In
// the zero-error cases this equals truth BY CONSTRUCTION of the study; the
// attitude-dependency tests deliberately break that equality to show why the
// assumption matters. The filter itself never receives a SpacecraftState.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "attitude/rigid_body.hpp"
#include "dynamics/two_body.hpp"
#include "environment/environment_models.hpp"
#include "estimation/diagnostics.hpp"
#include "estimation/kalman.hpp"
#include "estimation/translational_ekf.hpp"
#include "math/constants.hpp"
#include "math/linalg.hpp"
#include "numerics/fixed_step_propagation.hpp"
#include "orbit/two_body_orbit.hpp"
#include "sensors/gnss.hpp"
#include "sensors/imu.hpp"
#include "spacecraft/spacecraft_parameters.hpp"
#include "spacecraft/spacecraft_state.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace astradock;
using Catch::Matchers::WithinAbs;

namespace {

constexpr double kMu = constants::earth_gravitational_parameter_m3_per_s2;

// Tumbling rates used across scenarios (rad/s, BODY frame): slow enough that
// per-sample attitude change is tiny, nontrivial enough to exercise C_I_B.
constexpr math::Vector3 kTumbleRadS{0.01, -0.02, 0.015};

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

// Independent DCM built from quaternion components using the standard
// textbook formula for v_ECI = C(q) v_BODY = q [0,v] q*. Deliberately NOT
// calling Quaternion::rotate_vector / to_rotation_matrix so the production
// path is audited by an independently typed-out expression.
math::Matrix<3, 3> independent_dcm_eci_from_body(const math::Quaternion& q) {
    const double w = q.w();
    const double x = q.x();
    const double y = q.y();
    const double z = q.z();
    math::Matrix<3, 3> c;
    c(0, 0) = 1.0 - 2.0 * (y * y + z * z);
    c(0, 1) = 2.0 * (x * y - w * z);
    c(0, 2) = 2.0 * (x * z + w * y);
    c(1, 0) = 2.0 * (x * y + w * z);
    c(1, 1) = 1.0 - 2.0 * (x * x + z * z);
    c(1, 2) = 2.0 * (y * z - w * x);
    c(2, 0) = 2.0 * (x * z - w * y);
    c(2, 1) = 2.0 * (y * z + w * x);
    c(2, 2) = 1.0 - 2.0 * (x * x + y * y);
    return c;
}

math::Matrix<3, 3> matrix_from_columns(const math::Vector3& c0,
                                       const math::Vector3& c1,
                                       const math::Vector3& c2) {
    math::Matrix<3, 3> m;
    m(0, 0) = c0.x();
    m(1, 0) = c0.y();
    m(2, 0) = c0.z();
    m(0, 1) = c1.x();
    m(1, 1) = c1.y();
    m(2, 1) = c1.z();
    m(0, 2) = c2.x();
    m(1, 2) = c2.y();
    m(2, 2) = c2.z();
    return m;
}

// External reference attitude propagator: constant body rates under pure
// quaternion kinematics, RK4-integrated with step-boundary normalization.
// This is the "navigation layer" supplying the known/assumed attitude to the
// estimator; in the zero-error idealization it reproduces the truth attitude.
class ReferenceAttitudePropagator {
public:
    ReferenceAttitudePropagator(
        math::Quaternion initial_orientation,
        math::Vector3 angular_velocity_body_rad_s)
        : orientation_(std::move(initial_orientation)),
          omega_(angular_velocity_body_rad_s) {}

    void advance(double dt_s) {
        const attitude::RotationalState next = numerics::rk4_step(
            0.0,
            attitude::RotationalState{orientation_, omega_},
            dt_s,
            [](double, const attitude::RotationalState& s) {
                return attitude::RotationalState{
                    attitude::quaternion_derivative(
                        s.orientation, s.angular_velocity_rad_per_s),
                    math::Vector3{},
                };
            });
        orientation_ = next.orientation.normalized();
    }

    [[nodiscard]] const math::Quaternion& orientation() const noexcept {
        return orientation_;
    }

private:
    math::Quaternion orientation_;
    math::Vector3 omega_;
};

// Advances an attitude quaternion one IMU step under the scenario tumble
// rates (pure kinematics, step-boundary normalization) -- identical to
// ReferenceAttitudePropagator::advance so truth and reference stay in
// lock-step by construction.
math::Quaternion advance_truth_attitude(math::Quaternion orientation,
                                        double dt_s) {
    const attitude::RotationalState next = numerics::rk4_step(
        0.0,
        attitude::RotationalState{orientation, kTumbleRadS},
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

struct MultiRateConfig {
    std::uint64_t seed{1};
    double duration_s{600.0};
    double imu_dt_s{0.01};

    // TRUE sensor errors (what the M12 sensor model injects).
    math::Vector3 true_accel_bias_body_mps2{};
    double accel_noise_std_mps2{0.0};

    // Filter's knowledge of the bias (M13B: configuration, not estimated).
    math::Vector3 configured_accel_bias_body_mps2{};
    double configured_accel_noise_std_mps2{0.0};

    // GNSS.
    bool gnss_enabled{true};
    double gnss_sigma_r_m{10.0};
    double gnss_sigma_v_mps{0.05};
    bool apply_dropout{false};
    double dropout_start_s{300.0};
    double dropout_end_s{400.0};

    // Deliberate initialization error of the ESTIMATE (never truth access).
    math::Vector3 initial_position_error_m{50.0e3, -30.0e3, 20.0e3};
    math::Vector3 initial_velocity_error_mps{5.0, -3.0, 2.0};

    // Identity-attitude mode (bias study): keeps the bias vector fixed in
    // ECI so the classical dead-reckoning drift law is isolated from
    // attitude-frame rotation of the force direction.
    bool zero_tumble{false};

    // Attitude error injected into the supplied estimate (diagnostics only).
    double assumed_attitude_error_rad{0.0};
};

struct MultiRateResult {
    double initial_position_error_m{0.0};
    std::vector<double> time_s;
    std::vector<double> position_error_m;
    std::vector<double> velocity_error_m;
    std::vector<double> covariance_trace;
    std::vector<double> nees;
    std::vector<bool> measured;
    std::vector<double> nis;
    std::vector<double> raw_position_error_m;
    std::vector<double> raw_velocity_error_m;
    orbit::CartesianState final_estimate{};
};

// Runs the deterministic multi-rate experiment:
//   IMU at config.imu_dt_s (nominal 100 Hz), GNSS every 100th IMU step
//   (nominal 1 Hz), IMU prediction between fixes, GNSS update on arrival.
// During IMU dropouts (valid == false samples are not generated here; the
// sensor-level dropout policy is exercised separately) prediction falls back
// to documented dynamics-only propagation via predict(dt).
MultiRateResult run_multi_rate_scenario(const MultiRateConfig& cfg) {
    using constants::earth_reference_radius_m;

    // Truth trajectory: circular 500 km two-body orbit at IMU resolution.
    const double radius_m = earth_reference_radius_m + 500.0e3;
    const math::Vector3 scenario_omega =
        cfg.zero_tumble ? math::Vector3{} : kTumbleRadS;
    const orbit::CartesianState initial_truth{
        {radius_m, 0.0, 0.0},
        {0.0, orbit::circular_orbit_speed_m_per_s(kMu, radius_m), 0.0},
    };
    const auto truth_samples = numerics::propagate_fixed_step(
        0.0, cfg.duration_s, cfg.imu_dt_s, initial_truth,
        numerics::IntegrationMethod::classical_rk4,
        [](double, const orbit::CartesianState& s) {
            return orbit::two_body_state_derivative(0.0, s, kMu);
        });

    // Truth attitude: tumbling under pure quaternion kinematics, advanced in
    // lock-step with the external reference propagator below (free fall has
    // f_true ~ 0, so attitude does not perturb translation; it exercises
    // the sensor pipeline and C_I_B bookkeeping). CRITICAL INVARIANT: truth
    // and reference attitudes integrate IDENTICAL kinematics from IDENTICAL
    // initial conditions, so the known-attitude idealization holds exactly.
    spacecraft::SpacecraftState truth_sc{};
    truth_sc.translational = initial_truth;
    truth_sc.rotational = attitude::RotationalState{
        math::Quaternion::identity(), scenario_omega};

    // Sensors (M12 models; biases/noise per scenario).
    sensors::ImuConfig imu_cfg;
    imu_cfg.sample_period_s = cfg.imu_dt_s;
    imu_cfg.accel_bias_mps2 = cfg.true_accel_bias_body_mps2;
    imu_cfg.accel_noise_std_mps2 = cfg.accel_noise_std_mps2;
    imu_cfg.gyro_noise_std_rad_s = 0.0;
    imu_cfg.random_seed = cfg.seed;
    sensors::ImuSensor imu(imu_cfg);

    sensors::GnssConfig gnss_cfg;
    gnss_cfg.sample_period_s = 1.0;
    gnss_cfg.position_noise_std_m = cfg.gnss_sigma_r_m;
    gnss_cfg.velocity_noise_std_mps = cfg.gnss_sigma_v_mps;
    if (cfg.apply_dropout) {
        gnss_cfg.dropouts = {{cfg.dropout_start_s, cfg.dropout_end_s}};
    }
    gnss_cfg.random_seed = cfg.seed + 1'000'000u;
    sensors::GnssSensor gnss(gnss_cfg);

    // Estimator: deliberately wrong initialization; NEVER receives truth.
    estimation::TranslationalEkfConfig filter_config;
    filter_config.gravitational_parameter_m3_per_s2 = kMu;
    filter_config.acceleration_noise_std_mps2 = 1.0e-3;  // M13A path (unused here)
    filter_config.imu.accelerometer_bias_body_mps2 =
        cfg.configured_accel_bias_body_mps2;
    filter_config.imu.accelerometer_noise_std_mps2 =
        cfg.configured_accel_noise_std_mps2;
    filter_config.gnss_noise.position_variance_m2 =
        cfg.gnss_sigma_r_m * cfg.gnss_sigma_r_m;
    filter_config.gnss_noise.velocity_variance_m2_per_s2 =
        cfg.gnss_sigma_v_mps * cfg.gnss_sigma_v_mps;

    estimation::TranslationalCovariance p0 =
        estimation::TranslationalCovariance::zero();
    for (std::size_t i = 0; i < 3; ++i) {
        p0(i, i) = (25.0e3) * (25.0e3);
        p0(i + 3, i + 3) = 5.0 * 5.0;
    }

    estimation::TranslationalEkf filter(
        filter_config,
        orbit::CartesianState{
            initial_truth.position + cfg.initial_position_error_m,
            initial_truth.velocity + cfg.initial_velocity_error_mps},
        p0);

    MultiRateResult result;
    result.initial_position_error_m = cfg.initial_position_error_m.norm();
    result.time_s.reserve(truth_samples.size());
    result.position_error_m.reserve(truth_samples.size());

    ReferenceAttitudePropagator reference_attitude{
        math::Quaternion::identity()
            * math::Quaternion::from_axis_angle(
                math::Vector3{1.0, 2.0, 3.0},
                cfg.assumed_attitude_error_rad),
        scenario_omega};

    double previous_timestamp_s = 0.0;

    for (std::size_t k = 0; k < truth_samples.size(); ++k) {
        const double t_s = truth_samples[k].time_s;
        result.time_s.push_back(t_s);

        // Advance the TRUTH attitude (kinematic) and the external
        // navigation-layer attitude estimate for this epoch.
        if (k > 0) {
            truth_sc.rotational =
                numerics::rk4_step(
                    0.0, truth_sc.rotational, cfg.imu_dt_s,
                    [](double, const attitude::RotationalState& s) {
                        return attitude::RotationalState{
                            attitude::quaternion_derivative(
                                s.orientation, s.angular_velocity_rad_per_s),
                            math::Vector3{},
                        };
                    })
                    .normalized();
            reference_attitude.advance(cfg.imu_dt_s);
            CHECK(truth_sc.rotational.orientation.w()
                  == reference_attitude.orientation().w());
        }

        // IMU measurement through the M12 sensor model.
        truth_sc.translational = truth_samples[k].state;
        const sensors::ImuMeasurement imu_sample =
            imu.measure(t_s, truth_sc);

        if (k > 0) {
            const double dt_s = t_s - previous_timestamp_s;  // actual interval
            estimation::ImuPredictionInput prediction_input;
            prediction_input.dt_s = dt_s;
            prediction_input.measured_specific_force_body_mps2 =
                imu_sample.specific_force_body_mps2;
            prediction_input.measured_angular_velocity_body_rad_s =
                imu_sample.angular_velocity_body_rad_s;
            prediction_input.attitude =
                estimation::NavigationAttitudeEstimate{
                    reference_attitude.orientation()};
            filter.predict_with_imu(prediction_input);
        }
        previous_timestamp_s = t_s;

        // GNSS update at 1 Hz when enabled and not in a dropout window.
        bool applied_measurement = false;
        if (cfg.gnss_enabled && (k % 100u == 0u)) {
            const sensors::GnssMeasurement gnss_sample =
                gnss.measure(t_s, truth_sc);
            if (gnss_sample.valid) {
                filter.update_gnss(
                    gnss_sample.position_eci_m, gnss_sample.velocity_eci_mps);
                applied_measurement = true;
                result.nis.push_back(
                    filter.last_update_diagnostics().normalized_innovation_squared);
                result.raw_position_error_m.push_back(
                    (gnss_sample.position_eci_m - truth_samples[k].state.position)
                        .norm());
                result.raw_velocity_error_m.push_back(
                    (gnss_sample.velocity_eci_mps
                     - truth_samples[k].state.velocity)
                        .norm());
            } else {
                result.nis.push_back(std::numeric_limits<double>::quiet_NaN());
                result.raw_position_error_m.push_back(
                    std::numeric_limits<double>::quiet_NaN());
                result.raw_velocity_error_m.push_back(
                    std::numeric_limits<double>::quiet_NaN());
            }
        }
        result.measured.push_back(applied_measurement);

        const orbit::CartesianState estimate = filter.estimated_state();
        result.final_estimate = estimate;
        result.position_error_m.push_back(
            (truth_samples[k].state.position - estimate.position).norm());
        result.velocity_error_m.push_back(
            (truth_samples[k].state.velocity - estimate.velocity).norm());
        result.covariance_trace.push_back([&filter]() {
            const auto diagonal = filter.covariance().diagonal();
            double trace = 0.0;
            for (std::size_t i = 0; i < 6; ++i) {
                trace += diagonal(i, 0);
            }
            return trace;
        }());
        result.nees.push_back(estimation::translational_nees(
            truth_samples[k].state, estimate, filter.covariance()));
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
    REQUIRE(values.size() > 0);
    const std::size_t n = values.size();
    if (n % 2 == 1) {
        return values[n / 2];
    }
    return 0.5 * (values[n / 2 - 1] + values[n / 2]);
}

}  // namespace

// ---------------------------------------------------------------------------
// Layer 1: specific-force frame transformation (DCM direction audit)
// ---------------------------------------------------------------------------

TEST_CASE("Specific force transforms BODY to ECI through the supplied attitude "
          "with correct rotation direction",
          "[imu_ekf][frames]") {
    // q rotates ECI<-BODY by +90 degrees about z: body x-hat maps to ECI
    // y-hat. A conjugated (inverted) DCM would map body x-hat to NEGATIVE
    // ECI y — this test fails on that classic direction error.
    const math::Quaternion q = math::Quaternion::from_axis_angle(
        math::Vector3{0.0, 0.0, 1.0}, constants::pi / 2.0);

    // Position along +x: g points along -x with magnitude mu/r^2.
    const double radius_m = 7.0e6;
    const math::Vector3 position{radius_m, 0.0, 0.0};
    const double gMagnitude = kMu / (radius_m * radius_m);

    const math::Vector3 acceleration =
        estimation::imu_inertial_acceleration_eci(
            math::Vector3{1.0, 0.0, 0.0},  // f_B along body x
            math::Vector3{},               // no bias
            q,
            position,
            kMu);

    CHECK_THAT(acceleration.x(), WithinAbs(-gMagnitude, 1.0e-12));
    INFO("a_eci = " << acceleration.x() << ", " << acceleration.y() << ", "
         << acceleration.z());
    CHECK_THAT(acceleration.y(), WithinAbs(1.0, 1.0e-15));
    CHECK_THAT(acceleration.z(), WithinAbs(0.0, 1.0e-15));
}

TEST_CASE("Known-attitude acceleration reconstruction matches an independent "
          "textbook DCM computation across nontrivial attitudes",
          "[imu_ekf][analytic]") {
    struct Case {
        math::Quaternion q;
        math::Vector3 f_body;
        math::Vector3 bias_body;
    };

    const math::Vector3 position{-2.3e6, 5.1e6, -4.7e6};  // LEO-scale, off-axis

    const std::vector<Case> cases{
        // Identity attitude: acceleration must equal f + g exactly.
        {math::Quaternion::identity(),
         math::Vector3{0.3, -1.2, 0.7}, math::Vector3{}},
        // Large rotations about non-axis-aligned unit axes.
        {math::Quaternion::from_axis_angle(
             math::Vector3{1.0, 1.0, 1.0}, 2.0 / 3.0 * constants::pi),
         math::Vector3{0.05, 0.02, -0.03}, math::Vector3{0.001, -0.002, 0.0005}},
        {math::Quaternion::from_axis_angle(
             math::Vector3{-2.0, 5.0, 3.0}, 2.4),
         math::Vector3{10.0, -4.0, 2.5}, math::Vector3{}},
        {math::Quaternion::from_axis_angle(
             math::Vector3{0.0, 1.0, 0.0}, 0.05),  // small tilt
         math::Vector3{0.1, 0.0, 0.0}, math::Vector3{0.0, 0.0, 0.001}},
    };

    for (const auto& item : cases) {
        const math::Matrix<3, 3> c = independent_dcm_eci_from_body(item.q);
        const math::Vector3 f_corrected = item.f_body - item.bias_body;
        const math::Vector3 g = dynamics::two_body_acceleration(position, kMu);
        const math::Vector3 expected =
            math::Vector3{
                g.x() + c(0, 0) * f_corrected.x() + c(0, 1) * f_corrected.y()
                    + c(0, 2) * f_corrected.z(),
                g.y() + c(1, 0) * f_corrected.x() + c(1, 1) * f_corrected.y()
                    + c(1, 2) * f_corrected.z(),
                g.z() + c(2, 0) * f_corrected.x() + c(2, 1) * f_corrected.y()
                    + c(2, 2) * f_corrected.z()};

        const math::Vector3 produced =
            estimation::imu_inertial_acceleration_eci(
                item.f_body, item.bias_body, item.q, position, kMu);

        const math::Vector3 diff = produced - expected;
        INFO("diff = " << diff.x() << ", " << diff.y() << ", " << diff.z());
        CHECK(diff.norm() < 1.0e-12 * expected.norm());
    }
}

// ---------------------------------------------------------------------------
// Layer 1: prediction model analytics
// ---------------------------------------------------------------------------

TEST_CASE("IMU prediction reduces bitwise to the M13A dynamics-only prediction "
          "when measured specific force and bias cancel exactly",
          "[imu_ekf][regression]") {
    using constants::earth_reference_radius_m;
    const double radius_m = earth_reference_radius_m + 500.0e3;
    const orbit::CartesianState initial{
        {radius_m, 0.0, 0.0},
        {0.0, orbit::circular_orbit_speed_m_per_s(kMu, radius_m), 0.0},
    };

    estimation::TranslationalEkfConfig config;
    config.gravitational_parameter_m3_per_s2 = kMu;
    config.imu.accelerometer_noise_std_mps2 = 0.0;
    estimation::TranslationalEkf imu_filter(config, initial,
                                            estimation::TranslationalCovariance::identity());
    estimation::TranslationalEkf plain_filter(config, initial,
                                              estimation::TranslationalCovariance::identity());

    const estimation::NavigationAttitudeEstimate attitude{
        math::Quaternion::from_axis_angle(math::Vector3{3.0, -1.0, 2.0}, 1.1)};
    const estimation::ImuPredictionInput input{
        1.0,
        math::Vector3{},  // free fall: measured specific force is exactly zero
        kTumbleRadS,
        attitude};

    for (int step = 0; step < 300; ++step) {
        imu_filter.predict_with_imu(input);
        plain_filter.predict(1.0);
        const orbit::CartesianState a = imu_filter.estimated_state();
        const orbit::CartesianState b = plain_filter.estimated_state();
        CHECK(a.position.x() == b.position.x());
        CHECK(a.position.y() == b.position.y());
        CHECK(a.position.z() == b.position.z());
        CHECK(a.velocity.x() == b.velocity.x());
        CHECK(a.velocity.y() == b.velocity.y());
        CHECK(a.velocity.z() == b.velocity.z());
    }
}

TEST_CASE("IMU prediction mean satisfies the first-order small-dt limit",
          "[imu_ekf][analytic]") {
    using constants::earth_reference_radius_m;
    const double radius_m = earth_reference_radius_m + 500.0e3;
    const orbit::CartesianState initial{
        {radius_m, 0.0, 0.0},
        {0.0, 7540.0, 1200.0},
    };

    estimation::TranslationalEkfConfig config;
    config.gravitational_parameter_m3_per_s2 = kMu;
    estimation::TranslationalEkf filter(config, initial,
                                        estimation::TranslationalCovariance::identity());

    const math::Quaternion q = math::Quaternion::from_axis_angle(
        math::Vector3{1.0, -2.0, 0.5}, 0.9);
    const math::Vector3 f_body{0.02, -0.01, 0.005};
    // dt large enough that first-order deltas sit far above double-precision
    // cancellation floors (~ulp of the state components), small enough that
    // O(dt^2) truncation stays below 1e-6 relative.
    const double dt = 1.0e-4;

    const math::Vector3 expected_acceleration =
        estimation::imu_inertial_acceleration_eci(f_body, {}, q, initial.position, kMu);

    const estimation::ImuPredictionDiagnostics diagnostics =
        filter.predict_with_imu({dt, f_body, kTumbleRadS, {q}});

    // Diagnostic acceleration must equal the analytic reconstruction exactly
    // (same function, same inputs, evaluated at the pre-step state).
    const math::Vector3 accel_diff =
        diagnostics.inertial_acceleration_eci_mps2 - expected_acceleration;
    CHECK(accel_diff.norm() == 0.0);

    // ODE solution to first order: dr = dt*v, dv = dt*a0 (+ O(dt^2)).
    const math::Vector3 delta_position =
        diagnostics.predicted_position_eci_m - initial.position;
    const math::Vector3 delta_velocity =
        diagnostics.predicted_velocity_eci_mps - initial.velocity;
    const math::Vector3 expected_delta_position = initial.velocity * dt;
    const math::Vector3 expected_delta_velocity = expected_acceleration * dt;

    CHECK((delta_position - expected_delta_position).norm()
          < 1.0e-6 * expected_delta_position.norm());
    CHECK((delta_velocity - expected_delta_velocity).norm()
          < 1.0e-6 * expected_delta_velocity.norm());
}

TEST_CASE("IMU-aided covariance prediction uses the accelerometer-driven Q, "
          "not the M13A process noise",
          "[imu_ekf][covariance]") {
    using constants::earth_reference_radius_m;
    const double radius_m = earth_reference_radius_m + 500.0e3;
    const orbit::CartesianState initial{
        {radius_m, 0.0, 0.0},
        {0.0, 7500.0, 500.0},
    };

    // Two DIFFERENT sigma values: the filter must wire the IMU one into Q.
    constexpr double kM13aSigma = 7.0e-3;
    constexpr double kImuSigma = 2.0e-3;

    estimation::TranslationalEkfConfig config;
    config.gravitational_parameter_m3_per_s2 = kMu;
    config.acceleration_noise_std_mps2 = kM13aSigma;
    config.imu.accelerometer_noise_std_mps2 = kImuSigma;

    estimation::TranslationalCovariance p0 =
        estimation::TranslationalCovariance::zero();
    for (std::size_t i = 0; i < 6; ++i) {
        p0(i, i) = static_cast<double>(i + 1);
    }

    estimation::TranslationalEkf reference(config, initial, p0);
    const double dt = 0.5;
    static_cast<void>(reference.predict_with_imu({dt, {}, {}, {}}));

    // Independent expectation: Phi P0 Phi^T + Q(sigma_imu).
    const estimation::TranslationalCovariance phi =
        estimation::discrete_state_transition(
            estimation::continuous_dynamics_matrix(initial.position, kMu), dt);
    const estimation::TranslationalCovariance q_imu =
        estimation::imu_prediction_process_noise(kImuSigma, dt);
    const estimation::TranslationalCovariance q_wrong =
        estimation::discrete_process_noise(kM13aSigma, dt);
    const estimation::TranslationalCovariance expected =
        (phi * p0 * phi.transpose() + q_imu).symmetrized();

    const math::Matrix<6, 6> difference =
        reference.covariance() - expected;
    CHECK(math::max_abs_difference(difference,
                                   estimation::TranslationalCovariance::zero())
          < 1.0e-18);

    // And the wrong-Q covariance would have been measurably different.
    const estimation::TranslationalCovariance expected_wrong =
        (phi * p0 * phi.transpose() + q_wrong).symmetrized();
    CHECK(math::max_abs_difference(expected - expected_wrong,
                                   estimation::TranslationalCovariance::zero())
          > 1.0e-6);

    // Q scaling contract: position block dt^3, velocity block dt, cross dt^2.
    const double variance = kImuSigma * kImuSigma;
    const auto q_half = estimation::imu_prediction_process_noise(kImuSigma, dt / 2.0);
    CHECK_THAT(q_half(0, 0),
               WithinAbs(variance * (dt / 2.0 * dt / 2.0 * dt / 2.0) / 3.0, 1.0e-30));
    CHECK_THAT(q_half(3, 3), WithinAbs(variance * (dt / 2.0), 1.0e-30));
    CHECK_THAT(q_half(0, 3),
               WithinAbs(variance * (dt / 2.0 * dt / 2.0) / 2.0, 1.0e-30));
    // Isotropic per-axis construction: axes decoupled inside each block.
    CHECK(q_half(0, 1) == 0.0);
    CHECK(q_half(3, 4) == 0.0);
}

TEST_CASE("GNSS update following IMU prediction matches direct ekf_update algebra",
          "[imu_ekf][regression]") {
    using constants::earth_reference_radius_m;
    const double radius_m = earth_reference_radius_m + 500.0e3;
    const orbit::CartesianState initial_estimate{
        {radius_m + 100.0, 0.0, 0.0},
        {0.0, 7550.0, 0.0},
    };

    estimation::TranslationalEkfConfig config;
    config.gravitational_parameter_m3_per_s2 = kMu;
    config.acceleration_noise_std_mps2 = 1.0e-3;
    config.imu.accelerometer_noise_std_mps2 = 1.0e-3;
    config.gnss_noise.position_variance_m2 = 100.0;
    config.gnss_noise.velocity_variance_m2_per_s2 = 2.5e-3;

    estimation::TranslationalCovariance p0 =
        estimation::TranslationalCovariance::identity() * 400.0;

    estimation::TranslationalEkf filter(config, initial_estimate, p0);
    filter.predict_with_imu({1.0, {}, {}, {}});

    // Capture prior from public accessors, compute expected posterior with
    // the generic kalman.hpp primitives, then compare with update_gnss().
    const orbit::CartesianState prior = filter.estimated_state();
    const estimation::TranslationalCovariance prior_p = filter.covariance();

    const math::Vector3 measured_position = prior.position + math::Vector3{3.0, -4.0, 5.0};
    const math::Vector3 measured_velocity = prior.velocity + math::Vector3{0.01, 0.02, -0.03};

    const auto expected = estimation::ekf_update(
        estimation::translational_state_from_cartesian(prior),
        prior_p,
        estimation::translational_state_from_cartesian(
            orbit::CartesianState{measured_position, measured_velocity}),
        estimation::translational_state_from_cartesian(prior),
        estimation::gnss_measurement_matrix(),
        config.gnss_noise.covariance());

    static_cast<void>(filter.update_gnss(measured_position, measured_velocity));

    const orbit::CartesianState posterior = filter.estimated_state();
    CHECK(posterior.position.x() == expected.posterior_state(0, 0));
    CHECK(posterior.velocity.z() == expected.posterior_state(5, 0));
    CHECK(math::max_abs_difference(filter.covariance() - expected.posterior_covariance,
                                   estimation::TranslationalCovariance::zero())
          < 1.0e-12);
}

// ---------------------------------------------------------------------------
// Layer 1: free-fall orbital propagation
// ---------------------------------------------------------------------------

TEST_CASE("Free-fall IMU propagation reproduces the M04/M10 two-body reference "
          "trajectory bitwise under a tumbling known attitude",
          "[imu_ekf][freefall]") {
    using constants::earth_reference_radius_m;
    const double radius_m = earth_reference_radius_m + 500.0e3;
    const orbit::CartesianState initial_truth{
        {radius_m, 0.0, 0.0},
        {0.0, orbit::circular_orbit_speed_m_per_s(kMu, radius_m), 0.0},
    };

    const double imu_dt = 0.01;
    const double duration_s = 60.0;  // 6000 steps @ 100 Hz

    // Reference: open-loop RK4 propagation at the same rate (truth model).
    const auto reference = numerics::propagate_fixed_step(
        0.0, duration_s, imu_dt, initial_truth,
        numerics::IntegrationMethod::classical_rk4,
        [](double, const orbit::CartesianState& s) {
            return orbit::two_body_state_derivative(0.0, s, kMu);
        });

    // Ideal IMU: zero noise, zero bias => measured specific force is EXACTLY
    // zero in free fall regardless of attitude.
    sensors::ImuConfig imu_cfg;
    imu_cfg.sample_period_s = imu_dt;
    sensors::ImuSensor imu(imu_cfg);

    estimation::TranslationalEkfConfig config;
    config.gravitational_parameter_m3_per_s2 = kMu;
    config.imu.accelerometer_noise_std_mps2 = 0.0;
    estimation::TranslationalEkf filter(config, initial_truth,
                                        estimation::TranslationalCovariance::identity());

    ReferenceAttitudePropagator reference_attitude{
        math::Quaternion::identity(), kTumbleRadS};

    spacecraft::SpacecraftState truth_sc{};
    truth_sc.translational = initial_truth;
    truth_sc.rotational = attitude::RotationalState{
        math::Quaternion::identity(), kTumbleRadS};

    for (std::size_t k = 1; k < reference.size(); ++k) {
        const double t_s = reference[k].time_s;
        truth_sc.translational = reference[k].state;
        reference_attitude.advance(imu_dt);

        const sensors::ImuMeasurement sample = imu.measure(t_s, truth_sc);
        filter.predict_with_imu(
            {t_s - reference[k - 1].time_s,
             sample.specific_force_body_mps2,
             sample.angular_velocity_body_rad_s,
             {reference_attitude.orientation()}});

        // Same integrator, same model, ZOH of exact-zero force => bitwise
        // equality with the open-loop truth propagation at checkpoints.
        if (k % 500u == 0u || k == reference.size() - 1u) {
            // Same integrator, same model; ZOH of exact-zero force makes the
            // paths algebraically identical. Residual nm-scale differences
            // come ONLY from timestamp subtraction t_k - t_{k-1} rounding to
            // values microscopically off nominal_dt in binary floating point
            // -- which the filter CORRECTLY honors per the timestamp-
            // discipline contract, so exact-bitwise equality across two
            // different dt sequences is neither expected nor desired.
            const orbit::CartesianState estimated = filter.estimated_state();
            const orbit::CartesianState expected = reference[k].state;
            const double position_gap =
                (estimated.position - expected.position).norm();
            const double velocity_gap =
                (estimated.velocity - expected.velocity).norm();
            INFO("k=" << k << " pos gap " << position_gap << " m, vel gap "
                 << velocity_gap);
            CHECK(position_gap < 1.0e-6);
            CHECK(velocity_gap < 1.0e-8);
        }
    }
}

// ---------------------------------------------------------------------------
// Layer 2: non-gravitational force (drag) propagation
// ---------------------------------------------------------------------------

TEST_CASE("IMU propagation captures atmospheric drag through measured specific "
          "force where gravity-only prediction diverges",
          "[imu_ekf][drag]") {
    using constants::earth_reference_radius_m;
    const double radius_m = earth_reference_radius_m + 500.0e3;
    const orbit::CartesianState initial{
        {radius_m, 0.0, 0.0},
        {0.0, orbit::circular_orbit_speed_m_per_s(kMu, radius_m), 0.0},
    };
    const attitude::PrincipalInertia inertia{100.0, 75.0, 50.0};
    const spacecraft::SpacecraftParameters sc_params{kMu, inertia};
    const double mass_kg = 500.0;

    environment::EnvironmentalParameters env_params;
    env_params.mass_kg = mass_kg;
    env_params.drag_coefficient_cd = 2.2;
    env_params.drag_reference_area_m2 = 2.0;
    env_params.atmosphere_ref_alt_m = 500.0e3;
    env_params.atmosphere_ref_density_kg_per_m3 = 4.0e-12;
    env_params.atmosphere_scale_height_m = 60.0e3;
    const environment::EnvironmentConfiguration env_config{false, true, false, false};

    const double imu_dt = 0.01;
    const double duration_s = 400.0;

    // Truth: environmental 6-DOF propagation WITH drag.
    const auto truth = environment::propagate_spacecraft_environmental(
        spacecraft::SpacecraftState{
            initial,
            attitude::RotationalState{math::Quaternion::identity(), kTumbleRadS}},
        duration_s, imu_dt, sc_params, env_params, env_config);

    // Drag magnitude sanity: ~1 micrometer/s^2 class at 500 km.
    const math::Vector3 drag_ref = environment::drag_acceleration_eci(
        initial.velocity, initial.position, mass_kg, 2.2, 2.0,
        constants::earth_reference_radius_m,
        constants::earth_rotation_rate_rad_per_s,
        500.0e3, 4.0e-12, 60.0e3);
    const double drag_scale = drag_ref.norm();
    REQUIRE(drag_scale > 5.0e-7);
    REQUIRE(drag_scale < 5.0e-6);

    // High-grade accelerometer: resolving nano-g-class drag demands an
    // instrument far quieter than the nominal 1e-3 m/s^2 unit (documented in
    // Lesson 014); this isolates the physics being tested.
    sensors::ImuConfig imu_cfg;
    imu_cfg.sample_period_s = imu_dt;
    imu_cfg.accel_noise_std_mps2 = 1.0e-8;
    imu_cfg.random_seed = 77u;
    sensors::ImuSensor imu(imu_cfg);

    estimation::TranslationalEkfConfig config;
    config.gravitational_parameter_m3_per_s2 = kMu;
    config.imu.accelerometer_noise_std_mps2 = 1.0e-8;
    estimation::TranslationalEkf imu_filter(config, initial,
                                            estimation::TranslationalCovariance::identity());
    estimation::TranslationalEkf coast_filter(config, initial,
                                              estimation::TranslationalCovariance::identity());

    ReferenceAttitudePropagator reference_attitude{
        math::Quaternion::identity(), kTumbleRadS};

    // Truth attitude: MUST advance in lock-step with the reference
    // propagator so the known-attitude idealization holds exactly; otherwise
    // the sensor measures f in one frame and the filter rotates it with
    // another, sweeping the applied force direction and self-cancelling.
    math::Quaternion truth_attitude = math::Quaternion::identity();

    double first_sample_specific_force_norm = -1.0;
    for (std::size_t k = 1; k < truth.size(); ++k) {
        const double t_s = truth[k].time_s;
        const spacecraft::SpacecraftState truth_k = truth[k].state;
        reference_attitude.advance(imu_dt);
        truth_attitude = advance_truth_attitude(truth_attitude, imu_dt);

        // Non-gravitational acceleration driving the sensor (ECI -> BODY
        // happens INSIDE the sensor model using the TRUTH attitude).
        const environment::EnvironmentalEffects effects =
            environment::evaluate_environmental_effects(
                truth_k, sc_params, env_params, env_config);
        spacecraft::SpacecraftState sensor_state = truth_k;
        sensor_state.rotational = attitude::RotationalState{
            truth_attitude, kTumbleRadS};
        const sensors::ImuMeasurement sample =
            imu.measure(t_s, sensor_state, effects.acceleration_eci_mps2);
        if (first_sample_specific_force_norm < 0.0) {
            first_sample_specific_force_norm =
                sample.specific_force_body_mps2.norm();
        }

        const estimation::ImuPredictionInput input{
            t_s - truth[k - 1].time_s,  // actual interval (timestamp discipline)
            sample.specific_force_body_mps2,
            sample.angular_velocity_body_rad_s,
            {reference_attitude.orientation()}};
        imu_filter.predict_with_imu(input);
        coast_filter.predict(imu_dt);
    }

    const math::Vector3 final_truth_position = truth.back().state.position();
    const double imu_aided_error =
        (final_truth_position - imu_filter.estimated_state().position).norm();
    const double coast_error =
        (final_truth_position - coast_filter.estimated_state().position).norm();

    INFO("drag scale " << drag_scale << " m/s^2; |f_meas(t0)| "
         << first_sample_specific_force_norm);
    INFO("imu-aided final error " << imu_aided_error << " m; coast error "
         << coast_error << " m");

    // The accelerometer saw the drag (rotated into BODY by the truth
    // attitude; rotation preserves norm). Tolerance covers the injected
    // 1e-8 m/s^2 axis noise plus a generous margin.
    CHECK_THAT(first_sample_specific_force_norm,
               WithinAbs(drag_scale, 0.05 * drag_scale));

    // Dead reckoning THROUGH the measurement stays accurate; ignoring the
    // measured force diverges quadratically (~80 m over 400 s).
    CHECK(imu_aided_error < 1.0);
    CHECK(coast_error > 20.0);
    CHECK(coast_error > 20.0 * imu_aided_error);
}

// ---------------------------------------------------------------------------
// Layer 2: multi-rate IMU+GNSS convergence and dropout behavior
// ---------------------------------------------------------------------------

TEST_CASE("Multi-rate IMU+GNSS filter converges from a wrong initialization and "
          "beats raw GNSS",
          "[imu_ekf][convergence]") {
    MultiRateConfig cfg;
    cfg.seed = 42u;
    cfg.duration_s = 1200.0;
    cfg.accel_noise_std_mps2 = 1.0e-3;
    cfg.configured_accel_noise_std_mps2 = 1.0e-3;

    const MultiRateResult run = run_multi_rate_scenario(cfg);
    const std::size_t n = run.position_error_m.size();

    // RMSE over the converged half of EACH stream: raw arrays hold one entry
    // per GNSS update (~n/100), not per IMU step.
    const double ekf_pos_rmse =
        root_mean_square(run.position_error_m, n / 2);
    const auto& raw_pos = run.raw_position_error_m;
    const double raw_pos_rmse = root_mean_square(raw_pos, raw_pos.size() / 2);
    const double ekf_vel_rmse =
        root_mean_square(run.velocity_error_m, n / 2);
    const auto& raw_vel = run.raw_velocity_error_m;
    const double raw_vel_rmse = root_mean_square(raw_vel, raw_vel.size() / 2);

    INFO("EKF pos RMSE " << ekf_pos_rmse << " vs raw " << raw_pos_rmse);
    INFO("EKF vel RMSE " << ekf_vel_rmse << " vs raw " << raw_vel_rmse);
    CHECK(run.initial_position_error_m > 60.0e3);
    CHECK(ekf_pos_rmse < raw_pos_rmse);
    CHECK(ekf_vel_rmse < raw_vel_rmse);
    CHECK(run.position_error_m[n - 1] < 100.0);

    double mean_nis = 0.0;
    std::size_t count = 0;
    for (std::size_t i = run.nis.size() / 2; i < run.nis.size(); ++i) {
        if (!std::isnan(run.nis[i])) {
            mean_nis += run.nis[i];
            ++count;
        }
    }
    mean_nis /= static_cast<double>(count);
    INFO("mean NIS converged: " << mean_nis);
    CHECK(mean_nis > 4.0);
    CHECK(mean_nis < 8.0);

    // Final NEES within a generous single-run bound.
    CHECK(run.nees[n - 1] < 60.0);
}

TEST_CASE("GNSS dropout: dead-reckoning error grows sub-linearly, covariance "
          "grows monotonically, recovery restores information",
          "[imu_ekf][dropout]") {
    MultiRateConfig cfg;
    cfg.seed = 1234u;
    cfg.duration_s = 1200.0;
    cfg.accel_noise_std_mps2 = 1.0e-3;
    cfg.configured_accel_noise_std_mps2 = 1.0e-3;
    cfg.apply_dropout = true;
    cfg.dropout_start_s = 600.0;
    cfg.dropout_end_s = 700.0;

    const MultiRateResult run = run_multi_rate_scenario(cfg);

    const auto index_at = [&](double t_s) {
        return static_cast<std::size_t>(t_s / cfg.imu_dt_s);
    };

    // Covariance grows monotonically during pure prediction...
    const double trace_before = run.covariance_trace[index_at(599.0)];
    const double trace_late = run.covariance_trace[index_at(699.9)];
    CHECK(trace_late > trace_before);
    for (std::size_t i = index_at(601.0); i <= index_at(699.0); ++i) {
        if (run.covariance_trace[i] <= run.covariance_trace[i - 1]) {
            FAIL("covariance trace decreased during dropout at index " << i);
        }
    }

    // ...and contracts after measurements return.
    const double trace_recovered = run.covariance_trace[index_at(900.0)];
    INFO("trace before/late/recovered: " << trace_before << " / " << trace_late
         << " / " << trace_recovered);
    CHECK(trace_recovered < trace_late);

    // Dead reckoning bounds the outage error well below GNSS-fix scale.
    double max_dropout_position_error = 0.0;
    double max_dropout_velocity_error = 0.0;
    for (std::size_t i = index_at(600.0); i <= index_at(700.0); ++i) {
        max_dropout_position_error =
            std::max(max_dropout_position_error, run.position_error_m[i]);
        max_dropout_velocity_error =
            std::max(max_dropout_velocity_error, run.velocity_error_m[i]);
    }
    INFO("max dropout pos err " << max_dropout_position_error << " m, vel err "
         << max_dropout_velocity_error << " m/s");
    CHECK(max_dropout_position_error < 200.0);
    CHECK(max_dropout_velocity_error < 1.0);

    // Post-recovery accuracy returns to nominal levels.
    CHECK(root_mean_square(run.position_error_m, index_at(900.0)) < 25.0);
}

TEST_CASE("Accelerometer bias creates growing drift with linear-in-bias, "
          "quadratic-in-time position error",
          "[imu_ekf][bias]") {
    // Prediction-only runs isolate the dead-reckoning drift law:
    //   delta_v = b*t,  delta_r ~= b*t^2/2   (linear in bias, quadratic in t).
    // The identity-attitude configuration keeps the bias vector fixed in ECI
    // (a tumbling frame would rotate the force direction and partially
    // cancel itself); horizon 300 s keeps gravity-gradient curvature well
    // below the flat-space law.
    const double duration_s = 300.0;
    const std::vector<double> biases{0.0, 1.0e-4, 5.0e-4};

    std::vector<double> final_position_drift;
    std::vector<double> final_velocity_drift;

    for (const double bias : biases) {
        MultiRateConfig cfg;
        cfg.seed = 2024u;
        cfg.duration_s = duration_s;
        cfg.gnss_enabled = false;  // pure dead reckoning
        cfg.true_accel_bias_body_mps2 = math::Vector3{bias, 0.0, 0.0};
        cfg.configured_accel_bias_body_mps2 = math::Vector3{};  // UNCALIBRATED
        cfg.initial_position_error_m = math::Vector3{};         // start exact
        cfg.initial_velocity_error_mps = math::Vector3{};
        cfg.accel_noise_std_mps2 = 0.0;
        cfg.configured_accel_noise_std_mps2 = 0.0;
        cfg.zero_tumble = true;

        const MultiRateResult run = run_multi_rate_scenario(cfg);
        final_position_drift.push_back(run.position_error_m.back());
        final_velocity_drift.push_back(run.velocity_error_m.back());
    }

    INFO("drifts: " << final_position_drift[0] << ", "
         << final_position_drift[1] << ", " << final_position_drift[2]);

    // Zero bias: dead reckoning tracks truth to the timestamp-roundoff floor
    // (see the free-fall test discussion).
    CHECK(final_position_drift[0] < 1.0e-5);

    // Position drift scales LINEARLY with bias at fixed horizon (the t^2 law
    // is in TIME, not in bias magnitude).
    const double ratio_large_over_small =
        final_position_drift[2] / final_position_drift[1];
    CHECK_THAT(ratio_large_over_small, WithinAbs(5.0, 0.05));

    const double velocity_ratio =
        final_velocity_drift[2] / final_velocity_drift[1];
    CHECK_THAT(velocity_ratio, WithinAbs(5.0, 0.05));

    // Time scaling at fixed bias: halving time quarters the drift.
    MultiRateConfig short_cfg;
    short_cfg.seed = 2024u;
    short_cfg.duration_s = duration_s / 2.0;
    short_cfg.gnss_enabled = false;
    short_cfg.true_accel_bias_body_mps2 = math::Vector3{5.0e-4, 0.0, 0.0};
    short_cfg.accel_noise_std_mps2 = 0.0;
    short_cfg.configured_accel_noise_std_mps2 = 0.0;
    short_cfg.initial_position_error_m = math::Vector3{};
    short_cfg.initial_velocity_error_mps = math::Vector3{};
    short_cfg.zero_tumble = true;
    const MultiRateResult short_run = run_multi_rate_scenario(short_cfg);
    const double time_ratio =
        final_position_drift[2] / std::max(short_run.position_error_m.back(), 1.0e-12);
    CHECK_THAT(time_ratio, WithinAbs(4.0, 0.3));  // (300/150)^2 = 4

    // Absolute scale sanity: flat-space 0.5*b*t^2 predicts 22.5 m at
    // b=5e-4, t=300 s; gravity-gradient coupling bends this modestly.
    CHECK(final_position_drift[2] > 15.0);
    CHECK(final_position_drift[2] < 40.0);
    CHECK(final_velocity_drift[2] > 0.10);
    CHECK(final_velocity_drift[2] < 0.25);
}

TEST_CASE("Attitude dependency: identical specific force under different "
          "assumed attitudes produces different inertial acceleration",
          "[imu_ekf][attitude]") {
    const math::Vector3 position{7.0e6, 1.0e6, -0.5e6};
    const math::Vector3 f_body{0.1, 0.0, 0.0};  // thruster-class maneuver

    const math::Quaternion q_nominal = math::Quaternion::identity();
    const double delta_theta = 0.01;  // 10 mrad assumed-attitude error
    const math::Quaternion q_biased = math::Quaternion::from_axis_angle(
        math::Vector3{0.0, 1.0, 0.0}, delta_theta);

    const math::Vector3 a_nominal = estimation::imu_inertial_acceleration_eci(
        f_body, {}, q_nominal, position, kMu);
    const math::Vector3 a_biased = estimation::imu_inertial_acceleration_eci(
        f_body, {}, q_biased, position, kMu);

    // Gravity cancels in the difference (same position): the disagreement is
    // purely the rotated specific force.
    const double delta_a = (a_biased - a_nominal).norm();
    const double expected = 2.0 * std::sin(delta_theta / 2.0) * f_body.norm();
    INFO("delta|a| = " << delta_a << ", expected " << expected);
    CHECK_THAT(delta_a, WithinAbs(expected, 1.0e-12));

    // Small-angle scaling: doubling the attitude error scales the chord by
    // sin(delta)/sin(delta/2) -- 2x up to an O(delta^2) correction.
    const math::Quaternion q_biased_2x = math::Quaternion::from_axis_angle(
        math::Vector3{0.0, 1.0, 0.0}, 2.0 * delta_theta);
    const math::Vector3 a_biased_2x = estimation::imu_inertial_acceleration_eci(
        f_body, {}, q_biased_2x, position, kMu);
    const double delta_a_2x = (a_biased_2x - a_nominal).norm();
    const double expected_ratio =
        std::sin(2.0 * delta_theta / 2.0) / std::sin(delta_theta / 2.0);
    CHECK_THAT(delta_a_2x / delta_a, WithinAbs(expected_ratio, 1.0e-9));

    // A 10 mrad attitude error on a 0.1 m/s^2 force corrupts acceleration by
    // ~1 mm/s^2 -- larger than the nominal IMU noise floor. Translational
    // inertial navigation ultimately requires ESTIMATED attitude (M13C).
    CHECK(delta_a > 5.0e-4);
}

// ---------------------------------------------------------------------------
// Layer 2: robustness, timing discipline, failure policy
// ---------------------------------------------------------------------------

TEST_CASE("Timestamp jitter is honored through actual measurement intervals",
          "[imu_ekf][timing]") {
    using constants::earth_reference_radius_m;
    const double radius_m = earth_reference_radius_m + 500.0e3;
    const orbit::CartesianState initial{
        {radius_m, 0.0, 0.0},
        {0.0, orbit::circular_orbit_speed_m_per_s(kMu, radius_m), 0.0},
    };

    estimation::TranslationalEkfConfig config;
    config.gravitational_parameter_m3_per_s2 = kMu;
    config.imu.accelerometer_noise_std_mps2 = 0.0;
    estimation::TranslationalEkf jittered(config, initial,
                                          estimation::TranslationalCovariance::identity());
    estimation::TranslationalEkf uniform(config, initial,
                                         estimation::TranslationalCovariance::identity());

    const estimation::NavigationAttitudeEstimate attitude{
        math::Quaternion::from_axis_angle(math::Vector3{1.0, 1.0, 0.0}, 0.7)};
    const math::Vector3 f_body{1.0e-3, -2.0e-3, 5.0e-4};

    const int steps = 500;
    const double nominal_dt = 0.01;
    double elapsed_jittered_time = 0.0;
    double total_abs_jitter = 0.0;

    for (int k = 1; k <= steps; ++k) {
        // Deterministic +-1 ms wobble around the nominal schedule; intervals
        // stay strictly positive because |jitter| < nominal_dt.
        const double jitter = 1.0e-4 * std::sin(static_cast<double>(k) * 0.7);
        const double dt_jittered = nominal_dt + jitter;
        REQUIRE(dt_jittered > 0.0);
        static_cast<void>(jittered.predict_with_imu(
            {dt_jittered, f_body, {}, attitude}));
        static_cast<void>(uniform.predict_with_imu(
            {nominal_dt, f_body, {}, attitude}));
        elapsed_jittered_time += dt_jittered;
        total_abs_jitter += std::abs(jitter);
    }

    // The jittered run advanced a slightly different total time (bounded by
    // the accumulated wobble), so its final state must differ from the
    // uniform run by roughly |v| * delta_t -- measurable but small.
    const double elapsed_gap =
        elapsed_jittered_time - static_cast<double>(steps) * nominal_dt;
    INFO("elapsed gap " << elapsed_gap << " s, total |jitter| " << total_abs_jitter);
    CHECK(std::abs(elapsed_gap) < 1.0e-3);

    const math::Vector3 position_gap =
        jittered.estimated_state().position - uniform.estimated_state().position;
    const math::Vector3 velocity_gap =
        jittered.estimated_state().velocity - uniform.estimated_state().velocity;
    INFO("position gap " << position_gap.norm() << " m, velocity gap "
         << velocity_gap.norm() << " m/s");
    CHECK(position_gap.norm() < 10.0);   // ~|v| * 1e-3 s bound, generous
    CHECK(velocity_gap.norm() < 0.05);   // ~|a| * 1e-3 s bound, generous
    CHECK(position_gap.norm() > 0.0);    // the timing actually mattered
}

TEST_CASE("Invalid IMU prediction inputs fail loudly instead of propagating "
          "corruption",
          "[imu_ekf][validation]") {
    using constants::earth_reference_radius_m;
    const double radius_m = earth_reference_radius_m + 500.0e3;
    const orbit::CartesianState initial{
        {radius_m, 0.0, 0.0},
        {0.0, 7500.0, 0.0},
    };
    estimation::TranslationalEkfConfig config;
    config.gravitational_parameter_m3_per_s2 = kMu;
    estimation::TranslationalEkf filter(config, initial,
                                        estimation::TranslationalCovariance::identity());

    const estimation::NavigationAttitudeEstimate good_attitude{};
    const double nan_value = std::nan("");

    REQUIRE_THROWS_AS(filter.predict_with_imu({0.0, {}, {}, good_attitude}),
                      std::domain_error);                       // dt == 0
    REQUIRE_THROWS_AS(filter.predict_with_imu({-0.01, {}, {}, good_attitude}),
                      std::domain_error);                       // negative dt
    REQUIRE_THROWS_AS(
        filter.predict_with_imu({nan_value, {}, {}, good_attitude}),
        std::domain_error);                                     // NaN dt
    REQUIRE_THROWS_AS(
        filter.predict_with_imu({0.01, math::Vector3{nan_value, 0.0, 0.0}, {},
                                 good_attitude}),
        std::domain_error);                                     // NaN force
    REQUIRE_THROWS_AS(
        filter.predict_with_imu(
            {0.01, {}, math::Vector3{0.0, nan_value, 0.0}, good_attitude}),
        std::domain_error);                                     // NaN gyro
    REQUIRE_THROWS_AS(
        filter.predict_with_imu(
            {0.01, {}, {},
             estimation::NavigationAttitudeEstimate{math::Quaternion{
                 2.0, 0.0, 0.0, 0.0}}}),
        std::domain_error);                                     // non-unit q

    // NaN accelerometer bias in configuration also rejects.
    config.imu.accelerometer_bias_body_mps2 = math::Vector3{nan_value, 0.0, 0.0};
    estimation::TranslationalEkf bad_bias_filter(config, initial,
                                                 estimation::TranslationalCovariance::identity());
    REQUIRE_THROWS_AS(bad_bias_filter.predict_with_imu({0.01, {}, {}, good_attitude}),
                      std::domain_error);
}

TEST_CASE("IMU unavailability falls back to documented dynamics-only "
          "prediction without divergence",
          "[imu_ekf][failure]") {
    using constants::earth_reference_radius_m;
    const double radius_m = earth_reference_radius_m + 500.0e3;
    const orbit::CartesianState initial_truth{
        {radius_m, 0.0, 0.0},
        {0.0, orbit::circular_orbit_speed_m_per_s(kMu, radius_m), 0.0},
    };
    const auto reference = numerics::propagate_fixed_step(
        0.0, 30.0, 0.01, initial_truth,
        numerics::IntegrationMethod::classical_rk4,
        [](double, const orbit::CartesianState& s) {
            return orbit::two_body_state_derivative(0.0, s, kMu);
        });

    estimation::TranslationalEkfConfig config;
    config.gravitational_parameter_m3_per_s2 = kMu;
    config.imu.accelerometer_noise_std_mps2 = 0.0;
    estimation::TranslationalEkf filter(config, initial_truth,
                                        estimation::TranslationalCovariance::identity());

    // IMU outage for the whole window: the harness fallback policy is the
    // documented dynamics-only predict(dt). No invented measurements. The
    // fallback reproduces the open-loop reference to within the same
    // timestamp-roundoff tolerance documented in the free-fall test.
    for (std::size_t k = 1; k < reference.size(); ++k) {
        filter.predict(reference[k].time_s - reference[k - 1].time_s);
    }

    const orbit::CartesianState estimate = filter.estimated_state();
    CHECK(math::is_finite(estimate.position));
    CHECK(math::is_finite(estimate.velocity));
    CHECK((estimate.position - reference.back().state.position).norm()
          < 1.0e-6);
}

TEST_CASE("Estimator runs are bitwise deterministic for identical seeds",
          "[imu_ekf][determinism]") {
    MultiRateConfig cfg;
    cfg.seed = 777u;
    cfg.duration_s = 200.0;
    cfg.accel_noise_std_mps2 = 1.0e-3;
    cfg.configured_accel_noise_std_mps2 = 1.0e-3;

    const MultiRateResult first = run_multi_rate_scenario(cfg);
    const MultiRateResult second = run_multi_rate_scenario(cfg);

    REQUIRE(first.position_error_m.size() == second.position_error_m.size());
    for (std::size_t i = 0; i < first.position_error_m.size(); ++i) {
        CHECK(first.time_s[i] == second.time_s[i]);
        CHECK(first.position_error_m[i] == second.position_error_m[i]);
        CHECK(first.velocity_error_m[i] == second.velocity_error_m[i]);
        CHECK(first.nees[i] == second.nees[i]);
    }
    REQUIRE(first.nis.size() == second.nis.size());
    for (std::size_t i = 0; i < first.nis.size(); ++i) {
        CHECK(first.nis[i] == second.nis[i]);
    }
}

// ---------------------------------------------------------------------------
// Layer 1: Jacobian audit through the full IMU acceleration map
// ---------------------------------------------------------------------------

TEST_CASE("Finite-difference audit of the IMU-aided dynamics gradient across "
          "nontrivial attitudes",
          "[imu_ekf][jacobian][audit]") {
    // The IMU-aided acceleration a(r) = C(q)(f - b) + g(r) has gravity as its
    // ONLY position dependence: da/dr = G(r) regardless of attitude or force.
    // The audit verifies the full production map (including the C(q) term)
    // against central differences, at multiple attitudes, with a documented
    // step sweep (truncation O(h^2) vs roundoff ~ eps/h).
    const std::vector<math::Quaternion> attitudes{
        math::Quaternion::identity(),
        math::Quaternion::from_axis_angle(math::Vector3{1.0, 2.0, -1.0}, 1.9),
        math::Quaternion::from_axis_angle(math::Vector3{0.0, 0.0, 1.0},
                                          2.0 * constants::pi / 3.0),
        math::Quaternion::from_axis_angle(math::Vector3{-1.0, 0.5, 3.0}, 0.35),
    };
    const std::vector<math::Vector3> positions{
        math::Vector3{7000.0e3, 0.0, 0.0},
        math::Vector3{-3000.0e3, 5000.0e3, 8000.0e3},
    };
    const math::Vector3 f_body{0.04, -0.11, 0.07};
    const math::Vector3 bias{0.002, -0.001, 0.0015};
    const std::vector<double> steps_m{1.0e6, 1.0e5, 1.0e4, 1.0e3, 1.0e2, 10.0, 1.0, 0.1};

    double worst_best_relative_error = 0.0;

    for (const auto& attitude : attitudes) {
        for (const auto& position : positions) {
            const math::Matrix<3, 3> analytical =
                estimation::gravity_position_jacobian(position, kMu);
            const double analytical_norm = frobenius_norm(analytical);

            double best_relative_error = std::numeric_limits<double>::infinity();
            double best_step_m = 0.0;
            for (const double h : steps_m) {
                math::Matrix<3, 3> numerical = math::Matrix<3, 3>::zero();
                for (std::size_t j = 0; j < 3; ++j) {
                    const auto accel_at = [&](const math::Vector3& r) {
                        return estimation::imu_inertial_acceleration_eci(
                            f_body, bias, attitude, r, kMu);
                    };
                    const math::Vector3 plus = accel_at(offset_along(position, j, h));
                    const math::Vector3 minus = accel_at(offset_along(position, j, -h));
                    const math::Vector3 column = (plus - minus) / (2.0 * h);
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
            INFO("attitude case, best step " << best_step_m << " rel err "
                 << best_relative_error);
            worst_best_relative_error =
                std::max(worst_best_relative_error, best_relative_error);
        }
    }

    CHECK(worst_best_relative_error < 1.0e-10);
}

// ---------------------------------------------------------------------------
// Layer 2: statistical consistency (lightweight Monte Carlo) & truth safety
// ---------------------------------------------------------------------------

TEST_CASE("Monte Carlo IMU+GNSS study keeps NEES consistent across seeds",
          "[imu_ekf][montecarlo]") {
    constexpr int kRuns = 30;
    std::vector<double> nees_finals;
    std::vector<double> final_errors;
    int beats_raw = 0;

    for (int run = 0; run < kRuns; ++run) {
        MultiRateConfig cfg;
        cfg.seed = static_cast<std::uint64_t>(50'000u + run);
        cfg.duration_s = 400.0;
        cfg.accel_noise_std_mps2 = 1.0e-3;
        cfg.configured_accel_noise_std_mps2 = 1.0e-3;

        const MultiRateResult result = run_multi_rate_scenario(cfg);
        const std::size_t n = result.position_error_m.size();
        final_errors.push_back(result.position_error_m[n - 1]);
        nees_finals.push_back(result.nees[n - 1]);

        // RMSE over the converged half of EACH stream: raw arrays hold one
        // entry per GNSS update (~n/100), not per IMU step.
        const double ekf_rmse =
            root_mean_square(result.position_error_m, n / 2);
        const auto& raw = result.raw_position_error_m;
        const double raw_rmse = root_mean_square(raw, raw.size() / 2);
        if (ekf_rmse < raw_rmse) {
            ++beats_raw;
        }
    }

    const double median_nees = median_of(nees_finals);
    const double median_error = median_of(final_errors);
    const double beat_fraction =
        static_cast<double>(beats_raw) / static_cast<double>(kRuns);

    INFO("median NEES " << median_nees << "; median final error "
         << median_error << "; beat fraction " << beat_fraction);
    CHECK(median_nees > estimation::k_chi2_6dof_95_lower);
    CHECK(median_nees < estimation::k_chi2_6dof_95_upper);
    CHECK(median_error < 17.32);
    CHECK(beat_fraction >= 0.9);
}

TEST_CASE("Truth non-interference: multi-rate estimation leaves all truth and "
          "measurement records bitwise untouched",
          "[imu_ekf][noninterference]") {
    using constants::earth_reference_radius_m;
    const double radius_m = earth_reference_radius_m + 500.0e3;
    const orbit::CartesianState initial_truth{
        {radius_m, 0.0, 0.0},
        {0.0, orbit::circular_orbit_speed_m_per_s(kMu, radius_m), 0.0},
    };
    const auto truth_samples = numerics::propagate_fixed_step(
        0.0, 100.0, 0.01, initial_truth,
        numerics::IntegrationMethod::classical_rk4,
        [](double, const orbit::CartesianState& s) {
            return orbit::two_body_state_derivative(0.0, s, kMu);
        });

    // Deep copies taken BEFORE estimation.
    const std::vector<numerics::StateSample<orbit::CartesianState>>
        truth_snapshot = truth_samples;

    MultiRateConfig cfg;
    cfg.seed = 999u;
    cfg.duration_s = 100.0;
    cfg.accel_noise_std_mps2 = 1.0e-3;
    cfg.configured_accel_noise_std_mps2 = 1.0e-3;
    static_cast<void>(run_multi_rate_scenario(cfg));

    // Truth telemetry bitwise unchanged after estimation ran to completion.
    for (std::size_t k = 0; k < truth_samples.size(); ++k) {
        CHECK(truth_samples[k].state.position.x()
              == truth_snapshot[k].state.position.x());
        CHECK(truth_samples[k].state.velocity.z()
              == truth_snapshot[k].state.velocity.z());
    }
}
