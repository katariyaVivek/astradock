#pragma once

#include "math/constants.hpp"
#include "numerics/fixed_step_propagation.hpp"
#include "orbit/cartesian_state.hpp"
#include "orbit/two_body_orbit.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace astradock::orbit::diagnostics {

// M05 validation diagnostics for the ideal two-body propagator.
//
// These functions operate on states and trajectories produced by the C++ propagator
// and expose measurable error, invariant, and convergence metrics. They provide
// reusable, testable QA diagnostics without embedding analysis logic inside tests.

namespace detail {

inline void require_positive_finite(double value, const char* message) {
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::domain_error(message);
    }
}

}  // namespace detail

// ---------------------------------------------------------------------------
// Position and velocity closure errors
// ---------------------------------------------------------------------------

[[nodiscard]] inline double position_closure_error_m(
    const CartesianState& initial_state,
    const CartesianState& final_state) noexcept {
    return (final_state.position - initial_state.position).norm();
}

[[nodiscard]] inline double velocity_closure_error_m_per_s(
    const CartesianState& initial_state,
    const CartesianState& final_state) noexcept {
    return (final_state.velocity - initial_state.velocity).norm();
}

// ---------------------------------------------------------------------------
// Invariant errors and drift
// ---------------------------------------------------------------------------

// Absolute error between instantaneous specific energy and reference energy.
[[nodiscard]] inline double specific_energy_error_m2_per_s2(
    const CartesianState& state,
    double reference_energy_m2_per_s2,
    double gravitational_parameter_m3_per_s2) {
    const double energy = specific_orbital_energy_m2_per_s2(
        state, gravitational_parameter_m3_per_s2);
    return std::abs(energy - reference_energy_m2_per_s2);
}

// Relative error in specific mechanical energy |epsilon - epsilon_0| / |epsilon_0|.
[[nodiscard]] inline double relative_energy_error(
    const CartesianState& state,
    double reference_energy_m2_per_s2,
    double gravitational_parameter_m3_per_s2) {
    detail::require_positive_finite(
        std::abs(reference_energy_m2_per_s2),
        "Relative energy error: reference energy magnitude must be nonzero and finite");
    const double error = specific_energy_error_m2_per_s2(
        state, reference_energy_m2_per_s2, gravitational_parameter_m3_per_s2);
    return error / std::abs(reference_energy_m2_per_s2);
}

// Maximum relative energy drift over a full trajectory relative to the initial state.
[[nodiscard]] inline double max_relative_energy_drift(
    const std::vector<numerics::StateSample<CartesianState>>& samples,
    double gravitational_parameter_m3_per_s2) {
    if (samples.empty()) {
        return 0.0;
    }
    const double initial_energy = specific_orbital_energy_m2_per_s2(
        samples.front().state, gravitational_parameter_m3_per_s2);
    double max_drift = 0.0;
    for (const auto& sample : samples) {
        const double rel_err = relative_energy_error(
            sample.state, initial_energy, gravitational_parameter_m3_per_s2);
        max_drift = std::max(max_drift, rel_err);
    }
    return max_drift;
}

// Absolute error between instantaneous angular momentum magnitude and reference value.
[[nodiscard]] inline double specific_angular_momentum_magnitude_error_m2_per_s(
    const CartesianState& state,
    double reference_h_magnitude_m2_per_s) {
    const double h_mag = specific_angular_momentum_m2_per_s(state).norm();
    return std::abs(h_mag - reference_h_magnitude_m2_per_s);
}

// Relative error in specific angular momentum magnitude |h - h_0| / h_0.
[[nodiscard]] inline double relative_angular_momentum_magnitude_error(
    const CartesianState& state,
    double reference_h_magnitude_m2_per_s) {
    detail::require_positive_finite(
        reference_h_magnitude_m2_per_s,
        "Relative angular momentum error: reference magnitude must be positive and finite");
    const double error = specific_angular_momentum_magnitude_error_m2_per_s(
        state, reference_h_magnitude_m2_per_s);
    return error / reference_h_magnitude_m2_per_s;
}

// Maximum relative angular momentum magnitude drift over a trajectory.
[[nodiscard]] inline double max_relative_angular_momentum_drift(
    const std::vector<numerics::StateSample<CartesianState>>& samples) {
    if (samples.empty()) {
        return 0.0;
    }
    const double initial_h = specific_angular_momentum_m2_per_s(samples.front().state).norm();
    double max_drift = 0.0;
    for (const auto& sample : samples) {
        const double rel_err = relative_angular_momentum_magnitude_error(
            sample.state, initial_h);
        max_drift = std::max(max_drift, rel_err);
    }
    return max_drift;
}

// Maximum radial deviation |r(t) - r_ref| over a trajectory.
[[nodiscard]] inline double max_radius_deviation_m(
    const std::vector<numerics::StateSample<CartesianState>>& samples,
    double reference_radius_m) {
    detail::require_positive_finite(
        reference_radius_m,
        "Max radius deviation: reference radius must be positive and finite");
    double max_deviation = 0.0;
    for (const auto& sample : samples) {
        const double r = sample.state.position.norm();
        max_deviation = std::max(max_deviation, std::abs(r - reference_radius_m));
    }
    return max_deviation;
}

// ---------------------------------------------------------------------------
// Phase angle in the orbital (X-Y) plane
// ---------------------------------------------------------------------------
// For the simplified equatorial orbit used in M04/M05, the spacecraft stays in
// the X-Y plane. The phase angle is atan2(y, x), measured from the +X axis.
// Returns phase in radians in [-pi, pi].
[[nodiscard]] inline double phase_angle_rad(const CartesianState& state) {
    return std::atan2(state.position.y(), state.position.x());
}

// ---------------------------------------------------------------------------
// Analytical phase for a circular orbit
// ---------------------------------------------------------------------------
// For a circular orbit with mean motion n = sqrt(mu / r^3), the analytical
// phase at time t is theta(t) = n * t.
// The result is not wrapped: it grows linearly with time.
[[nodiscard]] inline double analytical_phase_rad(
    double time_s,
    double gravitational_parameter_m3_per_s2,
    double circular_radius_m) {
    detail::require_positive_finite(
        gravitational_parameter_m3_per_s2,
        "Analytical phase: gravitational parameter must be finite and positive");
    detail::require_positive_finite(
        circular_radius_m,
        "Analytical phase: orbital radius must be finite and positive");

    const double mean_motion_rad_per_s =
        std::sqrt(gravitational_parameter_m3_per_s2 / (circular_radius_m * circular_radius_m * circular_radius_m));
    if (!std::isfinite(mean_motion_rad_per_s)) {
        throw std::overflow_error("Analytical phase: mean motion is not representable");
    }
    return mean_motion_rad_per_s * time_s;
}

// Instantaneous phase error for a circular orbit at time t.
[[nodiscard]] inline double phase_error_rad(
    const CartesianState& state,
    double time_s,
    double gravitational_parameter_m3_per_s2,
    double circular_radius_m) {
    const double expected = analytical_phase_rad(
        time_s, gravitational_parameter_m3_per_s2, circular_radius_m);
    // Wrap expected to [-pi, pi] for comparison with wrapped instantaneous angle
    double wrapped_expected = std::fmod(expected + constants::pi, 2.0 * constants::pi);
    if (wrapped_expected < 0.0) {
        wrapped_expected += 2.0 * constants::pi;
    }
    wrapped_expected -= constants::pi;

    const double actual = phase_angle_rad(state);
    double diff = actual - wrapped_expected;
    while (diff > constants::pi) diff -= 2.0 * constants::pi;
    while (diff < -constants::pi) diff += 2.0 * constants::pi;
    return std::abs(diff);
}

// ---------------------------------------------------------------------------
// Unwrapped numerical phase from trajectory
// ---------------------------------------------------------------------------
// Starting from the initial sample phase, unwraps 2*pi phase discontinuities
// using standard step-to-step delta wrapping.
[[nodiscard]] inline std::vector<double> unwrapped_phase_rad(
    const std::vector<numerics::StateSample<CartesianState>>& samples) {
    if (samples.empty()) {
        return {};
    }

    std::vector<double> result;
    result.reserve(samples.size());
    double cumulative = phase_angle_rad(samples.front().state);
    result.push_back(cumulative);

    double prev_wrapped = cumulative;
    for (std::size_t i = 1; i < samples.size(); ++i) {
        const double curr_wrapped = phase_angle_rad(samples[i].state);
        double delta = curr_wrapped - prev_wrapped;
        while (delta > constants::pi) delta -= 2.0 * constants::pi;
        while (delta < -constants::pi) delta += 2.0 * constants::pi;
        cumulative += delta;
        prev_wrapped = curr_wrapped;
        result.push_back(cumulative);
    }

    return result;
}

// ---------------------------------------------------------------------------
// Phase error: max |numerical_unwrapped - analytical| over the trajectory
// ---------------------------------------------------------------------------
[[nodiscard]] inline double max_phase_error_rad(
    const std::vector<numerics::StateSample<CartesianState>>& samples,
    double gravitational_parameter_m3_per_s2,
    double circular_radius_m) {
    if (samples.size() < 2) {
        return 0.0;
    }

    const std::vector<double> phases = unwrapped_phase_rad(samples);
    double max_error = 0.0;

    for (std::size_t i = 0; i < samples.size(); ++i) {
        const double analytical =
            analytical_phase_rad(samples[i].time_s, gravitational_parameter_m3_per_s2, circular_radius_m);
        max_error = std::max(max_error, std::abs(phases[i] - analytical));
    }

    return max_error;
}

// ---------------------------------------------------------------------------
// Orbital period estimation from trajectory
// ---------------------------------------------------------------------------
// Detects the first return to the positive X-axis (crossing from y < 0 to y >= 0
// with x > 0) after at least 0.5 * analytical_period_hint_s.
// Uses linear interpolation between the two bracketing trajectory samples.
[[nodiscard]] inline double estimate_period_from_trajectory(
    const std::vector<numerics::StateSample<CartesianState>>& samples,
    double analytical_period_hint_s) {
    detail::require_positive_finite(
        analytical_period_hint_s,
        "Period estimation: analytical period hint must be finite and positive");

    if (samples.size() < 3) {
        throw std::domain_error(
            "Period estimation requires at least three trajectory samples");
    }

    const double min_duration = 0.5 * analytical_period_hint_s;
    if (samples.back().time_s < min_duration) {
        throw std::domain_error(
            "Trajectory is too short for period estimation");
    }

    for (std::size_t i = 1; i < samples.size(); ++i) {
        if (samples[i].time_s < min_duration) {
            continue;
        }

        const double y_prev = samples[i - 1].state.position.y();
        const double y_curr = samples[i].state.position.y();
        const double x_prev = samples[i - 1].state.position.x();
        const double x_curr = samples[i].state.position.x();

        // Detect crossing from y < 0 to y >= 0 on the positive-X half-plane (x > 0)
        if (y_prev < 0.0 && y_curr >= 0.0 && (x_prev > 0.0 || x_curr > 0.0)) {
            const double denom = y_curr - y_prev;
            if (denom <= 0.0) {
                continue;
            }
            const double fraction = -y_prev / denom;
            const double t0 = samples[i - 1].time_s;
            const double t1 = samples[i].time_s;
            const double estimated_period = t0 + fraction * (t1 - t0);
            if (std::isfinite(estimated_period) && estimated_period > 0.0) {
                return estimated_period;
            }
        }
    }

    throw std::domain_error("No orbital period crossing detected in trajectory");
}

// ---------------------------------------------------------------------------
// Empirical convergence order
// ---------------------------------------------------------------------------
// Given two timesteps h1 > h2 and their corresponding error values E1, E2:
// p = log(E1 / E2) / log(h1 / h2).
[[nodiscard]] inline double empirical_convergence_order(
    double h1,
    double error_at_h1,
    double h2,
    double error_at_h2) {
    if (!std::isfinite(h1) || !std::isfinite(h2) || h1 <= h2 || h2 <= 0.0) {
        throw std::domain_error(
            "Empirical convergence order: require h1 > h2 > 0 and both finite");
    }
    if (!std::isfinite(error_at_h1) || !std::isfinite(error_at_h2) ||
        error_at_h1 <= 0.0 || error_at_h2 <= 0.0) {
        throw std::domain_error(
            "Empirical convergence order: require positive, finite error values");
    }
    const double ratio_h = h1 / h2;
    const double ratio_e = error_at_h1 / error_at_h2;
    return std::log(ratio_e) / std::log(ratio_h);
}

// ---------------------------------------------------------------------------
// Maximum relative angular-momentum vector-direction drift
// ---------------------------------------------------------------------------
// Measures the maximum angle (in radians) between the initial angular momentum
// vector and the current angular momentum vector. This captures orbital plane
// tilting/precession caused by numerical error.
[[nodiscard]] inline double max_angular_momentum_direction_drift_rad(
    const std::vector<numerics::StateSample<CartesianState>>& samples) {
    if (samples.empty()) {
        return 0.0;
    }

    const math::Vector3 h0 = specific_angular_momentum_m2_per_s(samples.front().state);
    const double h0_norm = h0.norm();
    detail::require_positive_finite(
        h0_norm,
        "Angular momentum direction drift: initial angular momentum is zero or invalid");

    double max_angle = 0.0;
    for (const auto& sample : samples) {
        const math::Vector3 h = specific_angular_momentum_m2_per_s(sample.state);
        const double h_norm = h.norm();
        if (h_norm < 1.0e-200 || !std::isfinite(h_norm)) {
            continue;  // Skip degenerate/underflow vectors
        }
        const double dot_product = h0.dot(h) / (h0_norm * h_norm);
        const double clamped = std::clamp(dot_product, -1.0, 1.0);
        const double angle = std::acos(clamped);
        max_angle = std::max(max_angle, angle);
    }

    return max_angle;
}

}  // namespace astradock::orbit::diagnostics