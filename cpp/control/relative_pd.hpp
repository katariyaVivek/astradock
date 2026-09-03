#pragma once

// AstraDock M16C — LVLH relative-state PD translation control.
//
// Physical problem:
//   Hold a commanded LVLH offset (station keeping) or track a slow reference
//   trajectory using body thrust. The CW linearization (M15D) says the local
//   dynamics look like coupled oscillators plus drift, so a PD law on relative
//   position/velocity errors with the CW acceleration terms as feedforward gives
//   a physically grounded starting controller — strictly better than pretending
//   each axis is a double integrator, while still classical and transparent.
//
// Frames:
//   Errors are LVLH components (x radial, y along-track, z orbit-normal, M06/M15
//   convention). The commanded acceleration is LVLH; the caller rotates it into
//   the body frame for thruster allocation (M17), because allocation needs the
//   live attitude. Keeping this law in LVLH preserves the frame discipline that
//   a body-frame formulation would blur (LVLH error vs body actuation).
//
// Units (SI): m, m/s, m/s^2; gains Kp (1/s^2), Kd (1/s).
//
// Governing equations:
//   e_r = rho - rho_des,  e_v = v_rel - v_des
//   a_cmd = a_des + a_cw_ff(rho, v_rel, n) - Kp .* e_r - Kd .* e_v
// with feedforward a_cw_ff = -([3n^2 x + 2n vy], [-2n vx], [-n^2 z]) evaluated
// at the CURRENT state, which cancels the linearized natural dynamics so the
// PD terms only fight residual error (exact cancellation when CW is exact).
// With a_des = 0 and rho_des constant this is station keeping; nonzero v_des /
// a_des tracks a reference trajectory.
//
// Saturation discipline (M16E hook, enforced in the demo/tests via M14):
//   this header NEVER clips. It reports the vector a_cmd; desired-vs-achieved
//   accounting happens where thrust limits live (actuator assembly), because
//   only the allocator knows thruster geometry and mass.
//
// Assumptions: circular reference (CW feedforward), small separations, thrust
// achievable exactly along the commanded LVLH direction (allocation is M17).

#include "math/vector3.hpp"

#include <cmath>
#include <stdexcept>

namespace astradock::control {

namespace detail {

inline void require_finite_translation_inputs(
    const math::Vector3& rho,
    const math::Vector3& v_rel,
    const math::Vector3& rho_des,
    const math::Vector3& v_des,
    const math::Vector3& a_des,
    const math::Vector3& kp,
    const math::Vector3& kd) {
    if (!math::is_finite(rho) || !math::is_finite(v_rel) || !math::is_finite(rho_des)
        || !math::is_finite(v_des) || !math::is_finite(a_des)) {
        throw std::domain_error("Translation control states must contain only finite values");
    }
    if (!math::is_finite(kp) || !math::is_finite(kd)) {
        throw std::domain_error("Translation control gains must contain only finite values");
    }
    if (kp.x() <= 0.0 || kp.y() <= 0.0 || kp.z() <= 0.0 || kd.x() <= 0.0 || kd.y() <= 0.0
        || kd.z() <= 0.0) {
        throw std::domain_error("Translation control gains must be strictly positive per axis");
    }
}

}  // namespace detail

struct RelativePdGains {
    math::Vector3 kp_per_s2{1.0e-4, 1.0e-4, 1.0e-4};
    math::Vector3 kd_per_s{1.0e-2, 1.0e-2, 1.0e-2};
};

// CW natural-acceleration feedforward evaluated at the current relative state:
// the acceleration the linearized dynamics would produce without control. The
// controller subtracts it (see below) so PD acts on residual error only.
[[nodiscard]] inline math::Vector3 cw_feedforward_accel_lvlh_mps2(
    const math::Vector3& rho_lvlh_m,
    const math::Vector3& v_rel_lvlh_mps,
    double mean_motion_rad_per_s) {
    if (!math::is_finite(rho_lvlh_m) || !math::is_finite(v_rel_lvlh_mps)) {
        throw std::domain_error("Relative state must contain only finite values");
    }
    if (!std::isfinite(mean_motion_rad_per_s) || mean_motion_rad_per_s <= 0.0) {
        throw std::domain_error("Mean motion must be finite and strictly positive");
    }
    const double n = mean_motion_rad_per_s;
    return {
        3.0 * n * n * rho_lvlh_m.x() + 2.0 * n * v_rel_lvlh_mps.y(),
        -2.0 * n * v_rel_lvlh_mps.x(),
        -n * n * rho_lvlh_m.z(),
    };
}

// Full relative-state PD + CW feedforward law (LVLH acceleration command).
[[nodiscard]] inline math::Vector3 relative_pd_accel_lvlh_mps2(
    const math::Vector3& rho_lvlh_m,
    const math::Vector3& v_rel_lvlh_mps,
    const math::Vector3& rho_des_lvlh_m,
    const math::Vector3& v_des_lvlh_mps,
    const math::Vector3& a_des_lvlh_mps2,
    double mean_motion_rad_per_s,
    const RelativePdGains& gains) {
    detail::require_finite_translation_inputs(
        rho_lvlh_m, v_rel_lvlh_mps, rho_des_lvlh_m, v_des_lvlh_mps, a_des_lvlh_mps2,
        gains.kp_per_s2, gains.kd_per_s);
    if (!std::isfinite(mean_motion_rad_per_s) || mean_motion_rad_per_s <= 0.0) {
        throw std::domain_error("Mean motion must be finite and strictly positive");
    }
    const math::Vector3 e_r = rho_lvlh_m - rho_des_lvlh_m;
    const math::Vector3 e_v = v_rel_lvlh_mps - v_des_lvlh_mps;
    const math::Vector3 ff =
        cw_feedforward_accel_lvlh_mps2(rho_lvlh_m, v_rel_lvlh_mps, mean_motion_rad_per_s);
    return {
        a_des_lvlh_mps2.x() - ff.x() - gains.kp_per_s2.x() * e_r.x() - gains.kd_per_s.x() * e_v.x(),
        a_des_lvlh_mps2.y() - ff.y() - gains.kp_per_s2.y() * e_r.y() - gains.kd_per_s.y() * e_v.y(),
        a_des_lvlh_mps2.z() - ff.z() - gains.kp_per_s2.z() * e_r.z() - gains.kd_per_s.z() * e_v.z(),
    };
}

}  // namespace astradock::control
