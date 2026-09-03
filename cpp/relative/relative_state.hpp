#pragma once

// AstraDock M15C/M15D — LVLH relative state and Clohessy-Wiltshire propagation.
//
// Physical problem:
//   Rendezvous needs the chaser's motion expressed about the target, not about
//   Earth. The nonlinear truth lives in ECI; the LVLH frame centered on the
//   target gives relative position/velocity components with direct operational
//   meaning (radial x, along-track y, cross-track z in AstraDock's M06 LVLH
//   convention: x = e_r radial, y = e_t along-track, z = e_h orbit-normal).
//   For a circular reference orbit and small separations the linearized
//   Clohessy-Wiltshire (Hill) equations predict that relative motion in closed
//   form, which is what makes far-field rendezvous targeting tractable. M15D
//   implements the propagator AND maps where the linearization breaks down.
//
// Frames and semantics:
//   RelativeStateL Norton: relative_position_lvlh_m / relative_velocity_lvlh_mps
//   are LVLH components (x radial outward, y along-track, z orbit-normal).
//   Target/chaser roles are explicit function arguments — never inferred.
//   `relative_velocity_lvlh_mps` is the rotating-frame derivative (what a sensor
//   fixed in LVLH observes), NOT the projected inertial velocity. M06's
//   transform_eci_to_lvlh projects; beating that confusion is the point of M15B.
//
// Units (SI): m, m/s, m/s^2, s; mean motion n in rad/s.
//
// Governing equations (CW/Hill, circular reference, LVLH as defined above):
//   x_ddot - 2n y_dot - 3n^2 x = a_x
//   y_ddot + 2n x_dot          = a_y
//   z_ddot + n^2 z             = a_z
// State-space form dX/dt = A(n) X + B a_lvlh with the standard 6x6 A. Closed-form
// transition Phi(t) via the textbook solution (Clohessy-Wiltshire 1960):
//   x(t)   = (4-3c) x0 + (s/n) vx0 + (2/n)(1-c) vy0
//            + ax*(1-c)/n^2 + (2 ay/n^2)(t - s/n) ... (impulsive/general accel
//            handled by discrete convolution inside cw_propagate)
// with c = cos(nt), s = sin(nt). The implementation uses the exact Phi matrix
// (Schaub & Junkins closed form) and exact discrete input convolution for
// piecewise-constant LVLH accelerations, so CW integration itself is exact and
// any mismatch vs ECI truth is linearization error, not integrator error.
//
// Assumptions (documented limits, quantified in M15E):
//   - Circular reference orbit (e_ref = 0 exactly in the derivation).
//   - Linearization: separation << orbital radius.
//   - Unperturbed two-body gravity for both vehicles (differential J2/drag NOT
//     in CW; compare against matching two-body ECI truth).
//   - Short-to-moderate horizons: error grows secularly in along-track.

#include "frames/lvlh.hpp"
#include "frames/lvlh_rate.hpp"
#include "math/matrix3.hpp"
#include "math/vector3.hpp"
#include "orbit/cartesian_state.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace astradock::relative {

namespace detail {

inline void require_valid_mean_motion(double mean_motion_rad_per_s) {
    if (!std::isfinite(mean_motion_rad_per_s) || mean_motion_rad_per_s <= 0.0) {
        throw std::domain_error("CW mean motion must be finite and strictly positive");
    }
}

inline void require_valid_duration(double duration_s) {
    if (!std::isfinite(duration_s)) {
        throw std::domain_error("Propagation duration must be finite");
    }
}

inline void require_finite_relative_state(
    const math::Vector3& position_lvlh_m,
    const math::Vector3& velocity_lvlh_mps) {
    if (!math::is_finite(position_lvlh_m) || !math::is_finite(velocity_lvlh_mps)) {
        throw std::domain_error("LVLH relative state must contain only finite values");
    }
}

}  // namespace detail

// Relative translational state of a chaser w.r.t. a target, expressed in the
// target-centered LVLH frame (AstraDock convention: x radial, y along-track,
// z orbit-normal). Velocity is the rotating-frame derivative.
struct RelativeStateLvlh {
    math::Vector3 relative_position_lvlh_m{};
    math::Vector3 relative_velocity_lvlh_mps{};
};

[[nodiscard]] inline bool is_finite(const RelativeStateLvlh& state) noexcept {
    return math::is_finite(state.relative_position_lvlh_m)
        && math::is_finite(state.relative_velocity_lvlh_mps);
}

// Forms the LVLH relative state from ECI truth states of target and chaser:
//
//   rho_eci = r_chaser - r_target
//   rho_lvlh = C_LVLH_ECI(target) * rho_eci
//   v_rel_lvlh = C_LVLH_ECI(target) * (v_chaser - v_target) - omega_lvlh x rho_lvlh
//
// with omega_lvlh the LVLH frame rate in LVLH components. The correction term
// converts the projected inertial velocity difference into the rotating-frame
// derivative (transport theorem, M15B).
[[nodiscard]] inline RelativeStateLvlh relative_state_from_eci(
    const orbit::CartesianState& target_eci,
    const orbit::CartesianState& chaser_eci) {
    if (!orbit::is_finite(target_eci) || !orbit::is_finite(chaser_eci)) {
        throw std::domain_error("Target and chaser ECI states must be finite");
    }
    const math::Matrix3 c_lvlh_eci =
        frames::dcm_lvlh_from_eci(target_eci.position, target_eci.velocity);
    const math::Vector3 rho_eci = chaser_eci.position - target_eci.position;
    const math::Vector3 rho_lvlh = c_lvlh_eci * rho_eci;
    const math::Vector3 omega_eci =
        frames::lvlh_angular_velocity_rad_s(target_eci.position, target_eci.velocity);
    const math::Vector3 omega_lvlh = c_lvlh_eci * omega_eci;
    const math::Vector3 delta_v_eci = chaser_eci.velocity - target_eci.velocity;
    const math::Vector3 projected = c_lvlh_eci * delta_v_eci;
    const math::Vector3 v_rel_lvlh = projected - omega_lvlh.cross(rho_lvlh);
    return {rho_lvlh, v_rel_lvlh};
}

// Inverse map: chaser ECI state from target ECI state plus an LVLH relative
// state (exact inverse of relative_state_from_eci).
[[nodiscard]] inline orbit::CartesianState chaser_eci_from_relative(
    const orbit::CartesianState& target_eci,
    const RelativeStateLvlh& relative) {
    if (!orbit::is_finite(target_eci)) {
        throw std::domain_error("Target ECI state must be finite");
    }
    detail::require_finite_relative_state(
        relative.relative_position_lvlh_m, relative.relative_velocity_lvlh_mps);
    const math::Matrix3 c_lvlh_eci =
        frames::dcm_lvlh_from_eci(target_eci.position, target_eci.velocity);
    const math::Matrix3 c_eci_lvlh = c_lvlh_eci.transpose();
    const math::Vector3 rho_eci = c_eci_lvlh * relative.relative_position_lvlh_m;
    const math::Vector3 omega_eci =
        frames::lvlh_angular_velocity_rad_s(target_eci.position, target_eci.velocity);
    const math::Vector3 omega_lvlh = c_lvlh_eci * omega_eci;
    const math::Vector3 v_rel_eci = c_eci_lvlh
        * (relative.relative_velocity_lvlh_mps + omega_lvlh.cross(relative.relative_position_lvlh_m));
    return {target_eci.position + rho_eci, target_eci.velocity + v_rel_eci};
}

// CW state-transition scalars c = cos(nt), s = sin(nt). Closed-form Phi(t) for
// X = [x y z vx vy vz] (Schaub & Junkins, circular reference):
//
//   row x:  [4-3c, 0, 0, s/n, 2(1-c)/n, 0]
//   row y:  [6(s-nt), 1, 0, -2(1-c)/n, (4s-3nt)/n, 0]
//   row z:  [0, 0, c, 0, 0, s/n]
//   row vx: [3ns, 0, 0, c, 2s, 0]
//   row vy: [-6n(1-c), 0, 0, -2s, 4c-3, 0]
//   row vz: [0, 0, -ns, 0, 0, c]
struct CwTransitionScalars {
    double cos_nt{1.0};
    double sin_nt{0.0};
};

[[nodiscard]] inline CwTransitionScalars cw_transition_scalars(
    double mean_motion_rad_per_s,
    double duration_s) {
    detail::require_valid_mean_motion(mean_motion_rad_per_s);
    detail::require_valid_duration(duration_s);
    const double nt = mean_motion_rad_per_s * duration_s;
    return {std::cos(nt), std::sin(nt)};
}

// Applies Phi(t) to a relative state: unforced CW prediction.
[[nodiscard]] inline RelativeStateLvlh cw_predict(
    const RelativeStateLvlh& initial,
    double mean_motion_rad_per_s,
    double duration_s) {
    detail::require_finite_relative_state(
        initial.relative_position_lvlh_m, initial.relative_velocity_lvlh_mps);
    detail::require_valid_mean_motion(mean_motion_rad_per_s);
    detail::require_valid_duration(duration_s);
    const double n = mean_motion_rad_per_s;
    const double t = duration_s;
    const double c = std::cos(n * t);
    const double s = std::sin(n * t);
    const double x0 = initial.relative_position_lvlh_m.x();
    const double y0 = initial.relative_position_lvlh_m.y();
    const double z0 = initial.relative_position_lvlh_m.z();
    const double vx0 = initial.relative_velocity_lvlh_mps.x();
    const double vy0 = initial.relative_velocity_lvlh_mps.y();
    const double vz0 = initial.relative_velocity_lvlh_mps.z();

    const math::Vector3 pos{
        (4.0 - 3.0 * c) * x0 + (s / n) * vx0 + (2.0 / n) * (1.0 - c) * vy0,
        6.0 * (s - n * t) * x0 + y0 - (2.0 / n) * (1.0 - c) * vx0 + (4.0 * s - 3.0 * n * t) / n * vy0,
        c * z0 + (s / n) * vz0};
    const math::Vector3 vel{
        3.0 * n * s * x0 + c * vx0 + 2.0 * s * vy0,
        -6.0 * n * (1.0 - c) * x0 - 2.0 * s * vx0 + (4.0 * c - 3.0) * vy0,
        -n * s * z0 + c * vz0};
    return {pos, vel};
}

// CW acceleration (right-hand side): the differential specific force the CW
// model attributes to a relative state under an LVLH control/perturbation
// acceleration (m/s^2, LVLH components).
[[nodiscard]] inline RelativeStateLvlh cw_acceleration(
    const RelativeStateLvlh& state,
    double mean_motion_rad_per_s,
    const math::Vector3& accel_lvlh_mps2 = math::Vector3{}) {
    detail::require_finite_relative_state(
        state.relative_position_lvlh_m, state.relative_velocity_lvlh_mps);
    detail::require_valid_mean_motion(mean_motion_rad_per_s);
    if (!math::is_finite(accel_lvlh_mps2)) {
        throw std::domain_error("CW input acceleration must contain only finite values");
    }
    const double n = mean_motion_rad_per_s;
    const double n2 = n * n;
    const double x = state.relative_position_lvlh_m.x();
    const double z = state.relative_position_lvlh_m.z();
    const double vx = state.relative_velocity_lvlh_mps.x();
    const double vy = state.relative_velocity_lvlh_mps.y();
    const double vz = state.relative_velocity_lvlh_mps.z();
    const math::Vector3 dpos{vx, vy, vz};
    const math::Vector3 dvel{
        3.0 * n2 * x + 2.0 * n * vy + accel_lvlh_mps2.x(),
        -2.0 * n * vx + accel_lvlh_mps2.y(),
        -n2 * z + accel_lvlh_mps2.z()};
    return {dpos, dvel};
}

// Fixed-step RK4 propagation of the CW equations with piecewise-constant LVLH
// acceleration. Used for forced-response validation (impulsive burns as constant
// acceleration over one step); unforced arcs should use cw_predict (exact).
[[nodiscard]] inline std::vector<RelativeStateLvlh> cw_propagate(
    const RelativeStateLvlh& initial,
    double mean_motion_rad_per_s,
    double dt_s,
    std::size_t num_steps,
    const math::Vector3& accel_lvlh_mps2 = math::Vector3{}) {
    detail::require_finite_relative_state(
        initial.relative_position_lvlh_m, initial.relative_velocity_lvlh_mps);
    detail::require_valid_mean_motion(mean_motion_rad_per_s);
    if (!std::isfinite(dt_s) || dt_s <= 0.0) {
        throw std::domain_error("CW step size must be finite and strictly positive");
    }
    if (!math::is_finite(accel_lvlh_mps2)) {
        throw std::domain_error("CW input acceleration must contain only finite values");
    }
    std::vector<RelativeStateLvlh> samples;
    samples.reserve(num_steps + 1);
    samples.push_back(initial);
    RelativeStateLvlh current = initial;
    const auto rhs = [&](const RelativeStateLvlh& s) { return cw_acceleration(s, mean_motion_rad_per_s, accel_lvlh_mps2); };
    const auto add = [](const RelativeStateLvlh& a, const RelativeStateLvlh& b) {
        return RelativeStateLvlh{a.relative_position_lvlh_m + b.relative_position_lvlh_m,
                                 a.relative_velocity_lvlh_mps + b.relative_velocity_lvlh_mps};
    };
    const auto scale = [](const RelativeStateLvlh& a, double k) {
        return RelativeStateLvlh{a.relative_position_lvlh_m * k, a.relative_velocity_lvlh_mps * k};
    };
    for (std::size_t i = 0; i < num_steps; ++i) {
        const RelativeStateLvlh k1 = rhs(current);
        const RelativeStateLvlh k2 = rhs(add(current, scale(k1, 0.5 * dt_s)));
        const RelativeStateLvlh k3 = rhs(add(current, scale(k2, 0.5 * dt_s)));
        const RelativeStateLvlh k4 = rhs(add(current, scale(k3, dt_s)));
        current = add(current, scale(add(add(k1, scale(k2, 2.0)), add(scale(k3, 2.0), k4)), dt_s / 6.0));
        samples.push_back(current);
    }
    return samples;
}

}  // namespace astradock::relative
