#pragma once

#include "dynamics/two_body.hpp"
#include "estimation/imu_prediction.hpp"
#include "estimation/kalman.hpp"
#include "estimation/range_update.hpp"
#include "math/linalg.hpp"
#include "numerics/integrators.hpp"
#include "orbit/cartesian_state.hpp"
#include "orbit/two_body_orbit.hpp"

#include <cmath>
#include <stdexcept>

namespace astradock::estimation {

// M13A — Translational Extended Kalman Filter.
//
// State (Earth-Centered Inertial frame, SI units, fixed ordering):
//
//     x = [ r_x, r_y, r_z, v_x, v_y, v_z ]^T   (m, m/s)
//
// Covariance P is 6x6 with the same row/column ordering:
//
//     P = E[(x_true - x_estimate)(x_true - x_estimate)^T]
//
//     index | quantity        | unit
//     ------+-----------------+----------
//      0-2  | position (ECI)  | m^2
//      3-5  | velocity (ECI)  | m^2/s^2
//     cross | pos-vel corr.   | m^2/s
//
// Dynamics model (identical to the M04/M12 truth model for baseline work):
//
//     dr/dt = v
//     dv/dt = -mu * r / ||r||^3
//
// The filter mean is advanced with the same RK4 integrator used by the truth
// propagator. The error-state covariance uses the first-order discretization
// Phi = I + F dt of the linearized dynamics; see discrete_state_transition().
// Process noise models unmodeled accelerations as continuous white noise.

inline constexpr std::size_t k_translational_state_dim = 6;
inline constexpr std::size_t k_translational_position_index = 0;
inline constexpr std::size_t k_translational_velocity_index = 3;

using TranslationalCovariance = math::Matrix<k_translational_state_dim, k_translational_state_dim>;
using TranslationalMeanState = math::ColVector<k_translational_state_dim>;

// Converts an orbital Cartesian state into the 6x1 filter state vector.
[[nodiscard]] inline TranslationalMeanState translational_state_from_cartesian(
    const orbit::CartesianState& state) noexcept {
    return TranslationalMeanState(std::array<double, 6>{
        state.position.x(),
        state.position.y(),
        state.position.z(),
        state.velocity.x(),
        state.velocity.y(),
        state.velocity.z(),
    });
}

// Extracts the orbital Cartesian state from a 6x1 filter state vector.
[[nodiscard]] inline orbit::CartesianState cartesian_from_translational_state(
    const TranslationalMeanState& x) {
    if (!x.all_finite()) {
        throw std::domain_error("Translational state must contain only finite values");
    }
    return orbit::CartesianState{
        math::Vector3{x(0, 0), x(1, 0), x(2, 0)},
        math::Vector3{x(3, 0), x(4, 0), x(5, 0)},
    };
}

// Continuous-time gravity gradient G = da/dr at position r:
//
//     a(r) = -mu * r / ||r||^3
//     d a_i / d r_j = -mu * ( delta_ij / r^3 - 3 r_i r_j / r^5 )
//
// so G = -mu/r^3 * ( I - 3 rhat rhat^T ).
//
// Properties verified by tests: G is symmetric, traceless (the Newtonian
// potential is harmonic away from the origin), and has eigenvalues
// (+2 mu/r^3) along the radial direction and (-mu/r^3) along both transverse
// directions: radial separation grows, tangential separation oscillates.
// This analytical Jacobian is independently audited against central finite
// differences in the test suite; finite-differencing the production function
// alone would prove nothing about its correctness.
[[nodiscard]] inline math::Matrix<3, 3> gravity_position_jacobian(
    const math::Vector3& position_eci_m,
    double gravitational_parameter_m3_per_s2) {
    if (!std::isfinite(gravitational_parameter_m3_per_s2)
        || gravitational_parameter_m3_per_s2 <= 0.0) {
        throw std::domain_error(
            "Gravity Jacobian gravitational parameter must be finite and positive");
    }
    if (!math::is_finite(position_eci_m)) {
        throw std::domain_error("Gravity Jacobian position must contain only finite values");
    }

    const double radius_m = position_eci_m.norm();
    if (!(radius_m > 0.0)) {
        throw std::domain_error("Gravity Jacobian is undefined at the origin");
    }
    const double inverse_radius3 =
        1.0 / (radius_m * radius_m * radius_m);
    const double scale = -gravitational_parameter_m3_per_s2 * inverse_radius3;

    // Outer product rhat rhat^T.
    const math::Vector3 unit = position_eci_m / radius_m;
    math::Matrix<3, 3> outer;
    outer(0, 0) = unit.x() * unit.x();
    outer(0, 1) = unit.x() * unit.y();
    outer(0, 2) = unit.x() * unit.z();
    outer(1, 0) = unit.y() * unit.x();
    outer(1, 1) = unit.y() * unit.y();
    outer(1, 2) = unit.y() * unit.z();
    outer(2, 0) = unit.z() * unit.x();
    outer(2, 1) = unit.z() * unit.y();
    outer(2, 2) = unit.z() * unit.z();

    return scale * (math::Matrix<3, 3>::identity() - 3.0 * outer);
}

// Continuous-time dynamics Jacobian F = df/dx of the two-body system:
//
//         [  0_{3x3}   I_{3x3} ]
//     F = [                    ]
//         [    G         0     ]
//
// evaluated at the current estimate's position.
[[nodiscard]] inline TranslationalCovariance continuous_dynamics_matrix(
    const math::Vector3& position_eci_m,
    double gravitational_parameter_m3_per_s2) {
    const math::Matrix<3, 3> gravity_gradient = gravity_position_jacobian(
        position_eci_m, gravitational_parameter_m3_per_s2);

    TranslationalCovariance f = TranslationalCovariance::zero();
    for (std::size_t i = 0; i < 3; ++i) {
        // Top-right identity block: dr/dt = v.
        f(i, i + 3) = 1.0;
        // Bottom-left gravity-gradient block.
        for (std::size_t j = 0; j < 3; ++j) {
            f(i + 3, j) = gravity_gradient(i, j);
        }
    }
    return f;
}

// Discrete state transition matrix for the linearized error dynamics over dt.
//
//     Phi ~ I + F dt        (first-order discretization)
//
// Valid when the linearized dynamics change little over one step, i.e. when
// omega*dt << 1 where omega = sqrt(mu/r^3) is the local mean motion. For LEO
// with dt = 1 s this gives omega*dt ~ 1e-3, so neglected terms are O(1e-6).
// The second-order term (1/2) F^2 dt^2 could be added cheaply later if a
// validation case demands it; it is intentionally omitted to keep the first
// estimator minimal and its approximations explicit.
[[nodiscard]] inline TranslationalCovariance discrete_state_transition(
    const TranslationalCovariance& continuous_jacobian,
    double dt_s) {
    if (!std::isfinite(dt_s)) {
        throw std::domain_error("State transition timestep must be finite");
    }
    return math::Matrix<6, 6>::identity() + continuous_jacobian * dt_s;
}

// Discrete process noise covariance from a continuous white-noise acceleration
// model with standard deviation sigma_a (m/s^2 per axis).
//
// Model: acceleration disturbance w(t) enters the velocity equation directly
// (G_q = [0; I]) with spectral density sigma_a^2. Exact integration through
// the double-integrator chain dr/dt = v gives the standard piecewise-constant
// result:
//
//     Q = sigma_a^2 * [ dt^3/3 I   dt^2/2 I ]
//                     [ dt^2/2 I   dt     I ]
//
// Q expresses *model uncertainty* (what the dynamics got wrong between
// measurements); R expresses *sensor uncertainty*. They must not be tuned as a
// single knob. sigma_a has units m/s^2; a value of 1e-3 m/s^2 comfortably
// covers residual unmodeled LEO accelerations for this idealized baseline.
[[nodiscard]] inline TranslationalCovariance discrete_process_noise(
    double acceleration_noise_std_mps2,
    double dt_s) {
    if (!std::isfinite(acceleration_noise_std_mps2) || acceleration_noise_std_mps2 < 0.0) {
        throw std::domain_error("Process noise standard deviation must be finite and >= 0");
    }
    if (!std::isfinite(dt_s)) {
        throw std::domain_error("Process noise timestep must be finite");
    }

    const double variance = acceleration_noise_std_mps2 * acceleration_noise_std_mps2;
    const double c_position = variance * dt_s * dt_s * dt_s / 3.0;
    const double c_cross = variance * dt_s * dt_s / 2.0;
    const double c_velocity = variance * dt_s;

    TranslationalCovariance q = TranslationalCovariance::zero();
    for (std::size_t i = 0; i < 3; ++i) {
        q(k_translational_position_index + i, k_translational_position_index + i) = c_position;
        q(k_translational_velocity_index + i, k_translational_velocity_index + i) = c_velocity;
        q(k_translational_position_index + i, k_translational_velocity_index + i) = c_cross;
        q(k_translational_velocity_index + i, k_translational_position_index + i) = c_cross;
    }
    return q;
}

// Discrete translational process noise for IMU-AIDED propagation (M13B).
//
// The residual acceleration error after bias removal is dominated by the
// accelerometer white noise (isotropic sigma per axis, m/s^2); the two-body
// gravity model itself is exact for the M13B scenarios, and unmodeled-
// perturbation margin is intentionally zero because non-gravitational
// accelerations now enter THROUGH the measurement instead of being absent
// from both truth and model. Should a mission profile demand extra margin,
// an RSS combination sigma_total^2 = sigma_accel^2 + sigma_unmodeled^2 is
// the documented extension point; arbitrary Q inflation to "make the filter
// look good" is forbidden.
//
// Algebraically identical to discrete_process_noise():
//
//     Q = sigma^2 [ dt^3/3 I   dt^2/2 I ]
//                 [ dt^2/2 I   dt     I ]
//
// State ordering [r(0:2), v(3:5)]; acceleration noise isotropic in ECI under
// the M13B known-attitude idealization. Kept as a named function (rather
// than a raw alias) so tests can audit the Q-vs-IMU-noise contract directly.
[[nodiscard]] inline TranslationalCovariance imu_prediction_process_noise(
    double accelerometer_noise_std_mps2,
    double dt_s) {
    return discrete_process_noise(accelerometer_noise_std_mps2, dt_s);
}

// GNSS measurement model for M13A: z = [r_meas; v_meas] with h(x) = H x and
//
//     H = [ I_3  0_3 ]
//         [ 0_3  I_3 ]  = I_6
//
// This linear measurement keeps the first filter analytically checkable;
// nonlinear measurement models (range, attitude) arrive after this baseline
// is fully validated.
[[nodiscard]] inline TranslationalCovariance gnss_measurement_matrix() noexcept {
    return TranslationalCovariance::identity();
}

// Builds the GNSS measurement noise covariance R (isotropic per block).
struct GnssNoiseModel {
    double position_variance_m2{100.0};       // (sigma_r)^2, default sigma_r = 10 m
    double velocity_variance_m2_per_s2{2.5e-3}; // (sigma_v)^2, default sigma_v = 0.05 m/s

    [[nodiscard]] TranslationalCovariance covariance() const {
        if (!std::isfinite(position_variance_m2) || position_variance_m2 <= 0.0
            || !std::isfinite(velocity_variance_m2_per_s2)
            || velocity_variance_m2_per_s2 <= 0.0) {
            throw std::domain_error(
                "GNSS noise variances must be finite and strictly positive "
                "(zero R makes the innovation covariance singular)");
        }
        TranslationalCovariance r = TranslationalCovariance::zero();
        for (std::size_t i = 0; i < 3; ++i) {
            r(k_translational_position_index + i, k_translational_position_index + i) =
                position_variance_m2;
            r(k_translational_velocity_index + i, k_translational_velocity_index + i) =
                velocity_variance_m2_per_s2;
        }
        return r;
    }
};

// Configuration for the translational EKF.
//
// Two prediction paths share one filter instance:
//   - predict(dt)               dynamics-only propagation (M13A baseline),
//   - predict_with_imu(input)   IMU dead-reckoning propagation (M13B).
// The M13A field `acceleration_noise_std_mps2` parameterizes the FIRST path;
// the IMU-specific process noise of the second path is driven by
// `imu.accelerometer_noise_std_mps2` (see ImuPredictionConfig). Keeping them
// separate makes each Q derivation explicit and prevents the two baselines
// from silently sharing a tuning knob.
struct TranslationalEkfConfig {
    double gravitational_parameter_m3_per_s2{0.0};
    double acceleration_noise_std_mps2{1.0e-3};
    GnssNoiseModel gnss_noise{};
    ImuPredictionConfig imu{};
};

// One complete IMU prediction step record for telemetry/diagnostics.
struct ImuPredictionRecord {
    double dt_s{0.0};
    // Acceleration actually integrated across this interval, evaluated at the
    // pre-step estimate: a = C(q_hat)(f_m - b_a) + g(r_hat).
    math::Vector3 inertial_acceleration_eci_mps2{};
    math::Vector3 measured_specific_force_body_mps2{};
    math::Vector3 measured_angular_velocity_body_rad_s{};
};

// One complete GNSS measurement update record for telemetry/diagnostics.
struct GnssUpdateDiagnostics {
    math::Vector3 innovation_position_eci_m{};
    math::Vector3 innovation_velocity_eci_mps{};
    double normalized_innovation_squared{0.0};
    double kalman_gain_max_abs{0.0};
};

// Translational Extended Kalman Filter (M13A).
//
// Architectural contract: the filter consumes ONLY measurements. No API below
// accepts a truth state. Truth-dependent statistics (NEES, estimation error)
// are computed by evaluation harnesses outside the estimation module.
class TranslationalEkf {
public:
    TranslationalEkf(
        const TranslationalEkfConfig& config,
        const orbit::CartesianState& initial_estimate,
        const TranslationalCovariance& initial_covariance)
        : config_(config),
          state_(translational_state_from_cartesian(initial_estimate)),
          covariance_(initial_covariance) {
        validate_config(config);
        validate_covariance(covariance_, "TranslationalEkf: initial covariance");
    }

    [[nodiscard]] const TranslationalEkfConfig& config() const noexcept {
        return config_;
    }

    [[nodiscard]] orbit::CartesianState estimated_state() const {
        return cartesian_from_translational_state(state_);
    }

    [[nodiscard]] const TranslationalCovariance& covariance() const noexcept {
        return covariance_;
    }

    [[nodiscard]] const GnssUpdateDiagnostics& last_update_diagnostics() const noexcept {
        return last_diagnostics_;
    }

    [[nodiscard]] const ImuPredictionRecord& last_imu_prediction_diagnostics()
        const noexcept {
        return last_imu_record_;
    }

    [[nodiscard]] const RangeUpdateDiagnostics& last_range_diagnostics() const noexcept {
        return last_range_diagnostics_;
    }

    void reset(
        const orbit::CartesianState& initial_estimate,
        const TranslationalCovariance& initial_covariance) {
        state_ = translational_state_from_cartesian(initial_estimate);
        covariance_ = initial_covariance;
        validate_covariance(covariance_, "TranslationalEkf: reset covariance");
        last_diagnostics_ = GnssUpdateDiagnostics{};
        last_imu_record_ = ImuPredictionRecord{};
    }

    // Time update across dt seconds.
    //
    // Mean: RK4 step of the nonlinear two-body dynamics (same integrator as
    // the truth propagator, so zero process-model mismatch tests isolate pure
    // estimation behavior). Covariance: Phi P Phi^T + Q with the first-order
    // transition documented above.
    void predict(double dt_s) {
        if (!std::isfinite(dt_s)) {
            throw std::domain_error("EKF prediction timestep must be finite");
        }

        const auto process_model =
            [mu = config_.gravitational_parameter_m3_per_s2, dt_s](
                const TranslationalMeanState& x) -> TranslationalMeanState {
            const orbit::CartesianState propagated = numerics::rk4_step(
                0.0,
                cartesian_from_translational_state(x),
                dt_s,
                [mu](double, const orbit::CartesianState& s) {
                    return orbit::two_body_state_derivative(0.0, s, mu);
                });
            return translational_state_from_cartesian(propagated);
        };

        const TranslationalCovariance f = continuous_dynamics_matrix(
            math::Vector3{state_(0, 0), state_(1, 0), state_(2, 0)},
            config_.gravitational_parameter_m3_per_s2);
        const TranslationalCovariance phi =
            discrete_state_transition(f, dt_s);
        const TranslationalCovariance q = discrete_process_noise(
            config_.acceleration_noise_std_mps2, dt_s);

        const EkfPredictionResult<k_translational_state_dim> prediction =
            ekf_predict(state_, covariance_, process_model, phi, q);
        state_ = prediction.prior_state;
        covariance_ = prediction.prior_covariance;
    }

    // IMU-aided time update across dt seconds (M13B).
    //
    // Mean: RK4 integration of the inertial-navigation ODE over [t, t + dt]
    // with ZERO-ORDER HOLD on both the measured specific force and the
    // supplied attitude estimate (both constant across the interval; gravity
    // varies continuously with the integrated position):
    //
    //     dr/dt = v
    //     dv/dt = C_I_B(q_hat)(f_m - b_a) + g(r)
    //
    // Covariance: Phi = I + F dt with F = [[0, I], [G, 0]] evaluated at the
    // pre-step estimate position, plus Q built from the CONFIGURED
    // accelerometer noise via imu_prediction_process_noise().
    //
    // Documented linearization limitation (deliberate M13B staging): the
    // acceleration's dependence on attitude error, dC(q)/d(delta theta) * f,
    // contributes NO columns to F because attitude is treated as an exogenous
    // known input. The Jacobian is therefore exact only for the known-
    // attitude idealization; in free fall (f ~= 0) the neglected coupling is
    // identically zero anyway, which is why this approximation is benign for
    // orbital coast phases and must be revisited by the coupled M13C filter
    // whenever sustained non-gravitational forces act on the spacecraft.
    ImuPredictionDiagnostics predict_with_imu(const ImuPredictionInput& input) {
        config_.imu.validate();
        if (!std::isfinite(input.dt_s) || input.dt_s <= 0.0) {
            throw std::domain_error(
                "IMU prediction timestep must be finite and strictly positive "
                "(measurements must arrive in nondecreasing time order)");
        }
        if (!math::is_finite(input.measured_specific_force_body_mps2)
            || !math::is_finite(input.measured_angular_velocity_body_rad_s)) {
            throw std::domain_error("IMU prediction measurements must contain only finite values");
        }
        const math::Quaternion q_attitude = input.attitude.validated_orientation();

        const math::Vector3 pre_step_position{
            state_(0, 0), state_(1, 0), state_(2, 0)};

        const auto process_model =
            [&input, &q_attitude, mu = config_.gravitational_parameter_m3_per_s2,
             bias = config_.imu.accelerometer_bias_body_mps2,
             dt_s = input.dt_s](const TranslationalMeanState& x) -> TranslationalMeanState {
            const orbit::CartesianState propagated = numerics::rk4_step(
                0.0,
                cartesian_from_translational_state(x),
                dt_s,
                [mu, &q_attitude, bias, &input](double, const orbit::CartesianState& s) {
                    return orbit::CartesianState{
                        s.velocity,
                        imu_inertial_acceleration_eci(
                            input.measured_specific_force_body_mps2,
                            bias,
                            q_attitude,
                            s.position,
                            mu),
                    };
                });
            return translational_state_from_cartesian(propagated);
        };

        const TranslationalCovariance f = continuous_dynamics_matrix(
            math::Vector3{state_(0, 0), state_(1, 0), state_(2, 0)},
            config_.gravitational_parameter_m3_per_s2);
        const TranslationalCovariance phi = discrete_state_transition(f, input.dt_s);
        const TranslationalCovariance q = imu_prediction_process_noise(
            config_.imu.accelerometer_noise_std_mps2, input.dt_s);

        const EkfPredictionResult<k_translational_state_dim> prediction =
            ekf_predict(state_, covariance_, process_model, phi, q);
        state_ = prediction.prior_state;
        covariance_ = prediction.prior_covariance;

        last_imu_record_ = ImuPredictionRecord{
            input.dt_s,
            imu_inertial_acceleration_eci(
                input.measured_specific_force_body_mps2,
                config_.imu.accelerometer_bias_body_mps2,
                q_attitude,
                pre_step_position,
                config_.gravitational_parameter_m3_per_s2),
            input.measured_specific_force_body_mps2,
            input.measured_angular_velocity_body_rad_s,
        };

        return ImuPredictionDiagnostics{
            input.dt_s,
            last_imu_record_.inertial_acceleration_eci_mps2,
            math::Vector3{state_(0, 0), state_(1, 0), state_(2, 0)},
            math::Vector3{state_(3, 0), state_(4, 0), state_(5, 0)},
        };
    }

    // Measurement update with one GNSS solution. Returns false when no update
    // was applied (invalid flag), true otherwise.
    bool update_gnss(
        const math::Vector3& measured_position_eci_m,
        const math::Vector3& measured_velocity_eci_mps) {
        if (!math::is_finite(measured_position_eci_m)
            || !math::is_finite(measured_velocity_eci_mps)) {
            throw std::domain_error("GNSS update measurements must contain only finite values");
        }

        const EkfUpdateResult<6, 6> result = ekf_update(
            state_,
            covariance_,
            translational_state_from_cartesian(
                orbit::CartesianState{measured_position_eci_m, measured_velocity_eci_mps}),
            state_,  // h(x) = H x = x for the identity GNSS model
            gnss_measurement_matrix(),
            config_.gnss_noise.covariance());

        state_ = result.posterior_state;
        covariance_ = result.posterior_covariance;
        last_diagnostics_ = GnssUpdateDiagnostics{
            math::Vector3{
                result.innovation(0, 0),
                result.innovation(1, 0),
                result.innovation(2, 0),
            },
            math::Vector3{
                result.innovation(3, 0),
                result.innovation(4, 0),
                result.innovation(5, 0),
            },
            result.normalized_innovation_squared,
            max_abs(result.kalman_gain),
        };
        return true;
    }

    // Measurement update with one scalar range observation (M13C).
    // Updates position estimates along the line of sight.
    bool update_range(
        double measured_range_m,
        const math::Vector3& target_position_eci_m,
        double range_noise_std_m) {
        if (!std::isfinite(measured_range_m) || !math::is_finite(target_position_eci_m)) {
            throw std::domain_error("Range update measurements must contain only finite values");
        }

        const auto result = execute_range_update(
            state_,
            covariance_,
            measured_range_m,
            target_position_eci_m,
            range_noise_std_m,
            &last_range_diagnostics_
        );

        state_ = result.posterior_state;
        covariance_ = result.posterior_covariance;
        return true;
    }

private:
    [[nodiscard]] static double max_abs(const TranslationalCovariance& matrix) noexcept {
        double worst = 0.0;
        for (const double value : matrix.data()) {
            const double magnitude = value < 0.0 ? -value : value;
            if (magnitude > worst) {
                worst = magnitude;
            }
        }
        return worst;
    }

    static void validate_config(const TranslationalEkfConfig& config) {
        if (!std::isfinite(config.gravitational_parameter_m3_per_s2)
            || config.gravitational_parameter_m3_per_s2 <= 0.0) {
            throw std::domain_error(
                "Translational EKF gravitational parameter must be finite and positive");
        }
        if (!std::isfinite(config.acceleration_noise_std_mps2)
            || config.acceleration_noise_std_mps2 < 0.0) {
            throw std::domain_error(
                "Translational EKF process noise standard deviation must be finite and >= 0");
        }
        // Validates positivity rules inside GnssNoiseModel.
        static_cast<void>(config.gnss_noise.covariance());
    }

    TranslationalEkfConfig config_;
    TranslationalMeanState state_{};
    TranslationalCovariance covariance_{};
    GnssUpdateDiagnostics last_diagnostics_{};
    ImuPredictionRecord last_imu_record_{};
    RangeUpdateDiagnostics last_range_diagnostics_{};
};

}  // namespace astradock::estimation
