#pragma once

#include "math/linalg.hpp"

#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>

namespace astradock::estimation {

// Generic Extended Kalman Filter (EKF) prediction and measurement update.
//
// Physical problem
// ----------------
// The navigation computer must infer a spacecraft state from imperfect
// measurements. The Kalman filter is the classical probabilistic answer for
// linear-Gaussian systems. Spacecraft dynamics (two-body gravity) and several
// measurement models are nonlinear, so the *Extended* KF linearizes about the
// current estimate:
//
//     x_{k+1} = f(x_k) + w_k          (nonlinear process)
//     z_k     = h(x_k) + v_k          (nonlinear measurement)
//
// with zero-mean, white, mutually independent process and measurement noise:
//
//     w ~ N(0, Q),   v ~ N(0, R)
//
// The filter maintains a Gaussian belief N(x, P):
//
//     x = E[x_true],  P = E[(x_true - x)(x_true - x)^T]
//
// Prediction (time update), with discrete transition Phi approximating the
// linearized dynamics over dt:
//
//     x_prior = f(x_posterior)
//     P_prior = Phi P Phi^T + Q
//
// Measurement update, with H the Jacobian of h at x_prior:
//
//     y = z - h(x_prior)                 (innovation: what the sensor saw minus
//                                         what the filter expected)
//     S = H P H^T + R                    (innovation covariance)
//     K = P H^T S^{-1}                   (Kalman gain)
//     x_post = x_prior + K y
//     P_post = (I - K H) P (I - K H)^T + K R K^T      (Joseph form)
//
// Numerical policy
// ----------------
// - All solves use Cholesky factorization; no explicit matrix inverse is ever
//   formed.
// - The Joseph form is preferred over the shorter (I - KH)P because it better
//   preserves symmetry and positive semidefiniteness when K is computed in
//   finite precision. Both agree in exact arithmetic; see Lesson 013.
// - Covariances are symmetrized after every justified numerical operation;
//   arbitrary entry clamping is forbidden.
// - Non-finite states, non-positive covariance diagonals, or non-SPD
//   innovation covariances throw immediately: a silently wrong filter is worse
//   than a stopped one.

// Maximum tolerated asymmetry when validating an incoming covariance matrix.
// The check is RELATIVE to the largest matrix element: covariance blocks mix
// unit scales (position variances reach ~(7.5e4 m)^2 ~ 5.6e9 m^2 while
// velocity variances are O(10)), so an absolute tolerance would either be
// meaningless or reject physically valid matrices.
inline constexpr double k_covariance_symmetry_relative_tolerance = 1.0e-10;
inline constexpr double k_covariance_symmetry_absolute_floor_m2 = 1.0e-12;

// Relative tolerance applied to covariance diagonals when rejecting negative
// variances. Tiny negative values (~1e-18 relative) can appear from roundoff
// before symmetrization and are repaired by symmetrization; anything beyond
// this tolerance indicates a broken filter.
inline constexpr double k_covariance_negative_diagonal_relative_tolerance = 1.0e-12;

// Validates that a candidate covariance matrix is finite, symmetric within
// tolerance, and has no negative diagonal entries beyond roundoff tolerance.
template <std::size_t N>
void validate_covariance(const math::Matrix<N, N>& covariance, const char* context) {
    if (!covariance.all_finite()) {
        throw std::domain_error(
            std::string(context) + ": covariance must contain only finite values");
    }
    double max_abs_element = 0.0;
    for (const double value : covariance.data()) {
        const double magnitude = value < 0.0 ? -value : value;
        if (magnitude > max_abs_element) {
            max_abs_element = magnitude;
        }
    }
    const double symmetry_tolerance =
        math::approximately_equal(max_abs_element, 0.0)
            ? k_covariance_symmetry_absolute_floor_m2
            : k_covariance_symmetry_relative_tolerance * max_abs_element;
    if (covariance.symmetry_error() > symmetry_tolerance) {
        throw std::domain_error(
            std::string(context) + ": covariance is not symmetric within tolerance");
    }
    const auto diagonal = covariance.diagonal();
    double scale = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        const double value = diagonal(i, 0);
        const double magnitude = value < 0.0 ? -value : value;
        if (magnitude > scale) {
            scale = magnitude;
        }
    }
    const double negative_tolerance =
        k_covariance_negative_diagonal_relative_tolerance * scale;
    for (std::size_t i = 0; i < N; ++i) {
        if (diagonal(i, 0) < -negative_tolerance) {
            throw std::domain_error(
                std::string(context) + ": covariance diagonal contains negative variance");
        }
    }
}

// Result of an EKF prediction step.
template <std::size_t N>
struct EkfPredictionResult {
    math::Matrix<N, 1> prior_state{};
    math::Matrix<N, N> prior_covariance{};
};

// Result of an EKF measurement update step, including full diagnostics.
template <std::size_t N, std::size_t M>
struct EkfUpdateResult {
    math::Matrix<N, 1> posterior_state{};
    math::Matrix<N, N> posterior_covariance{};

    // Diagnostics exposed so the learner can inspect filter behavior:
    math::Matrix<M, 1> innovation{};              // y = z - h(x_prior)
    math::Matrix<M, M> innovation_covariance{};   // S = H P H^T + R
    math::Matrix<N, M> kalman_gain{};             // K = P H^T S^{-1}
    double normalized_innovation_squared{0.0};    // NIS = y^T S^{-1} y
};

// EKF time update.
//
// process_model is any callable f: Matrix<N,1> -> Matrix<N,1> advancing the
// nominal state mean across dt (for M13A this wraps an RK4 two-body step).
// Phi is the discrete state transition matrix of the linearized error
// dynamics; Q the discrete process noise covariance. The caller owns both
// discretizations; this function only performs the covariance algebra and its
// numerical validation.
template <std::size_t N, typename ProcessModel>
[[nodiscard]] EkfPredictionResult<N> ekf_predict(
    const math::Matrix<N, 1>& posterior_state,
    const math::Matrix<N, N>& posterior_covariance,
    ProcessModel&& process_model,
    const math::Matrix<N, N>& state_transition,
    const math::Matrix<N, N>& process_noise) {
    validate_covariance(
        posterior_covariance, "EKF predict: prior (posterior-input) covariance");
    if (!posterior_state.all_finite()) {
        throw std::domain_error("EKF predict: state must contain only finite values");
    }
    if (!state_transition.all_finite() || !process_noise.all_finite()) {
        throw std::domain_error("EKF predict: Phi and Q must contain only finite values");
    }

    EkfPredictionResult<N> result;
    result.prior_state = std::invoke(
        std::forward<ProcessModel>(process_model), posterior_state);
    if (!result.prior_state.all_finite()) {
        throw std::overflow_error("EKF predict: process model produced a non-finite state");
    }

    result.prior_covariance =
        (state_transition * posterior_covariance * state_transition.transpose()
         + process_noise)
            .symmetrized();
    validate_covariance(result.prior_covariance, "EKF predict: predicted covariance");
    return result;
}

// EKF measurement update (§16 of the M13 specification).
//
// z               measurement vector
// z_predicted     h(x_prior), evaluated by the caller at the prior estimate
// H               measurement Jacobian dh/dx at x_prior (h linear => H constant)
// R               measurement noise covariance (must be positive definite so S
//                 stays invertible even for a perfectly certain prior)
template <std::size_t N, std::size_t M>
[[nodiscard]] EkfUpdateResult<N, M> ekf_update(
    const math::Matrix<N, 1>& prior_state,
    const math::Matrix<N, N>& prior_covariance,
    const math::Matrix<M, 1>& measurement,
    const math::Matrix<M, 1>& predicted_measurement,
    const math::Matrix<M, N>& measurement_jacobian,
    const math::Matrix<M, M>& measurement_noise) {
    if (!prior_state.all_finite() || !measurement.all_finite()
        || !predicted_measurement.all_finite() || !measurement_jacobian.all_finite()
        || !measurement_noise.all_finite()) {
        throw std::domain_error("EKF update: inputs must contain only finite values");
    }
    validate_covariance(prior_covariance, "EKF update: prior covariance");

    EkfUpdateResult<N, M> result;

    // Innovation: disagreement between sensor and prediction.
    result.innovation = measurement - predicted_measurement;

    // Innovation covariance S = H P H^T + R, symmetrized against roundoff.
    result.innovation_covariance =
        (measurement_jacobian * prior_covariance * measurement_jacobian.transpose()
         + measurement_noise)
            .symmetrized();

    // Kalman gain via SPD solve: K = P H^T S^{-1} = (S^{-1} H P)^T using the
    // symmetry of S^{-1} and P. No explicit inverse is formed.
    result.kalman_gain =
        math::solve_spd_multi(result.innovation_covariance,
                              measurement_jacobian * prior_covariance)
            .transpose();

    // Normalized innovation squared: scalar consistency statistic; for a
    // consistent filter NIS follows approximately chi-square(M).
    const math::Matrix<M, 1> whitened =
        math::solve_spd(result.innovation_covariance, result.innovation);
    double nis = 0.0;
    for (std::size_t i = 0; i < M; ++i) {
        nis += result.innovation(i, 0) * whitened(i, 0);
    }
    if (!std::isfinite(nis) || nis < 0.0) {
        throw std::overflow_error("EKF update: NIS computation produced an invalid value");
    }
    result.normalized_innovation_squared = nis;

    // Posterior mean.
    result.posterior_state = prior_state + result.kalman_gain * result.innovation;

    // Joseph-form posterior covariance: numerically safer than (I-KH)P.
    const math::Matrix<N, N> identity_minus_kh =
        math::Matrix<N, N>::identity() - result.kalman_gain * measurement_jacobian;
    result.posterior_covariance =
        (identity_minus_kh * prior_covariance * identity_minus_kh.transpose()
         + result.kalman_gain * measurement_noise * result.kalman_gain.transpose())
            .symmetrized();

    validate_covariance(result.posterior_covariance, "EKF update: posterior covariance");
    return result;
}

}  // namespace astradock::estimation
