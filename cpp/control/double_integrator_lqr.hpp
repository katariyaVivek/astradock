#pragma once

// AstraDock M16D — One small LQR example on a double-integrator axis.
//
// Physical problem (teaching goal):
//   PD gains say "how stiff / how damped" but never answer "stiff at what
//   cost". LQR makes the tradeoff explicit: minimize J = ∫(x'Qx + u'Ru) with
//   state cost Q and control cost R, and the optimal linear law u = -Kx drops
//   out of one Riccati equation. AstraDock teaches this once, on the smallest
//   problem where every number can be hand-checked: one translational axis
//   (double integrator), then uses the resulting K as a cross-check on PD
//   intuition — not as a general optimal-control framework (explicit non-goal).
//
// Plant: dx/dt = A x + B u with A = [[0,1],[0,0]], B = [0,1]' (1 kg mass, so
//   u is acceleration in m/s^2). Cost: Q = diag(q_pos, q_vel), R = r scalar.
// Closed-form solution (solve by hand in the lesson, verified in tests):
//   P = [[p11, p12],[p12, p22]], p12 = sqrt(q_pos * r),
//       p22 = sqrt(r (2 p12 + q_vel)), p11 = p12 p22 / r,
//   K = [p12, p22] / r.
// Stability guarantee: A - BK is Hurwitz for all q_pos, q_vel, r > 0, with
//   closed-loop damping ratio ζ = K_vel / (2 sqrt(K_pos)).
//
// Units: SI (m, m/s, m/s^2); Q in (1/m^2, 1/(m/s)^2)-weighted units, R in
//   1/(m/s^2)^2. These are tuning weights, not physical constants — the lesson
//   documents the engineering basis for the chosen values.

#include "math/vector3.hpp"

#include <array>
#include <cmath>
#include <stdexcept>

namespace astradock::control {

namespace detail {

inline void require_positive_lqr_weights(double q_pos, double q_vel, double r) {
    if (!std::isfinite(q_pos) || !std::isfinite(q_vel) || !std::isfinite(r)) {
        throw std::domain_error("LQR weights must contain only finite values");
    }
    if (q_pos <= 0.0 || q_vel <= 0.0 || r <= 0.0) {
        throw std::domain_error("LQR weights must be strictly positive");
    }
}

}  // namespace detail

// LQR tuning weights for one double-integrator axis.
struct DoubleIntegratorLqrWeights {
    double q_position{1.0};
    double q_velocity{1.0};
    double r_control{1.0};
};

// Optimal gain row K = [k_pos, k_vel] minimizing J for the weights.
struct DoubleIntegratorLqrGain {
    double k_position{1.0};
    double k_velocity{1.0};
};

[[nodiscard]] inline DoubleIntegratorLqrGain solve_double_integrator_lqr(
    const DoubleIntegratorLqrWeights& weights) {
    detail::require_positive_lqr_weights(weights.q_position, weights.q_velocity, weights.r_control);
    const double p12 = std::sqrt(weights.q_position * weights.r_control);
    const double p22 = std::sqrt(weights.r_control * (2.0 * p12 + weights.q_velocity));
    return {p12 / weights.r_control, p22 / weights.r_control};
}

// Applies u = -(K x) = -(k_pos * position + k_vel * velocity).
[[nodiscard]] inline double apply_double_integrator_lqr(
    const DoubleIntegratorLqrGain& gain,
    double position_m,
    double velocity_mps) {
    if (!std::isfinite(gain.k_position) || !std::isfinite(gain.k_velocity)) {
        throw std::domain_error("LQR gains must contain only finite values");
    }
    if (!std::isfinite(position_m) || !std::isfinite(velocity_mps)) {
        throw std::domain_error("LQR state must contain only finite values");
    }
    return -(gain.k_position * position_m + gain.k_velocity * velocity_mps);
}

// Closed-loop eigenvalues of A - BK (real parts must be negative). Returns the
// pair (lambda1, lambda2); complex conjugates share one real part.
[[nodiscard]] inline std::array<double, 2> double_integrator_closed_loop_eigenvalues(
    const DoubleIntegratorLqrGain& gain) {
    if (!std::isfinite(gain.k_position) || !std::isfinite(gain.k_velocity)
        || gain.k_position <= 0.0 || gain.k_velocity <= 0.0) {
        throw std::domain_error("LQR gains must be finite and strictly positive");
    }
    // s^2 + k_vel s + k_pos = 0.
    const double discriminant =
        gain.k_velocity * gain.k_velocity - 4.0 * gain.k_position;
    if (discriminant >= 0.0) {
        const double root = std::sqrt(discriminant);
        return {(-gain.k_velocity + root) / 2.0, (-gain.k_velocity - root) / 2.0};
    }
    return {-gain.k_velocity / 2.0, -gain.k_velocity / 2.0};
}

// Closed-loop damping ratio ζ = k_vel / (2 sqrt(k_pos)).
[[nodiscard]] inline double double_integrator_damping_ratio(const DoubleIntegratorLqrGain& gain) {
    if (!std::isfinite(gain.k_position) || !std::isfinite(gain.k_velocity)
        || gain.k_position <= 0.0 || gain.k_velocity <= 0.0) {
        throw std::domain_error("LQR gains must be finite and strictly positive");
    }
    return gain.k_velocity / (2.0 * std::sqrt(gain.k_position));
}

}  // namespace astradock::control
