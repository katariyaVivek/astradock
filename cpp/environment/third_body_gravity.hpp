#pragma once

#include "math/vector3.hpp"

#include <cmath>
#include <stdexcept>

namespace astradock::environment {

// Computes the tidal gravitational perturbation acceleration exerted by an external
// third body (e.g. Moon, Sun) on a spacecraft in an Earth-Centered Inertial (ECI) frame.
//
// Physical Model:
//   In an Earth-centered reference frame, Earth itself is accelerating toward the third
//   body under gravitational attraction. The net apparent acceleration of the spacecraft
//   relative to Earth consists of:
//     1. Direct gravitational acceleration of the spacecraft toward the third body:
//          a_direct = mu_3 * (r_3 - r) / |r_3 - r|^3
//     2. Indirect (inertial reaction) acceleration of Earth toward the third body:
//          a_indirect = - mu_3 * r_3 / |r_3|^3
//
//   Total third-body tidal acceleration:
//     a_3B = a_direct + a_indirect = mu_3 * [ (r_3 - r) / |r_3 - r|^3 - r_3 / |r_3|^3 ]
//
// Frame and Units:
//   - spacecraft_position_eci_m: Spacecraft position relative to Earth in ECI (m)
//   - third_body_position_eci_m: Third body position relative to Earth in ECI (m)
//   - third_body_mu_m3_per_s2: Gravitational parameter of the third body mu_3 (m^3/s^2)
//   - Returns: Net perturbing acceleration in ECI frame (m/s^2)
//
// Invariant & Analytical Properties:
//   - When r = 0 (spacecraft at Earth center), a_3B = 0 exactly.
//   - When mu_3 = 0, a_3B = 0 exactly.
//   - For |r| << |r_3|, expands to classical tidal dipole: a_3B ~ (mu_3 / |r_3|^3) * [ 3 * (r . r_hat_3) * r_hat_3 - r ]
[[nodiscard]] inline math::Vector3 third_body_acceleration_eci(
    const math::Vector3& spacecraft_position_eci_m,
    const math::Vector3& third_body_position_eci_m,
    double third_body_mu_m3_per_s2) {
    if (!math::is_finite(spacecraft_position_eci_m) || !math::is_finite(third_body_position_eci_m)) {
        throw std::domain_error("Spacecraft and third-body positions must be finite");
    }
    if (!std::isfinite(third_body_mu_m3_per_s2) || third_body_mu_m3_per_s2 < 0.0) {
        throw std::domain_error("Third-body gravitational parameter must be finite and non-negative");
    }

    if (third_body_mu_m3_per_s2 == 0.0) {
        return math::Vector3{0.0, 0.0, 0.0};
    }

    const double r3_sq = third_body_position_eci_m.squared_norm();
    if (r3_sq <= 0.0) {
        throw std::domain_error("Third-body position relative to Earth cannot be zero");
    }

    const math::Vector3 r_rel = third_body_position_eci_m - spacecraft_position_eci_m;
    const double r_rel_sq = r_rel.squared_norm();
    if (r_rel_sq <= 0.0) {
        throw std::domain_error("Spacecraft cannot be coincident with third body point mass");
    }

    const double r3 = std::sqrt(r3_sq);
    const double r3_cubed = r3_sq * r3;

    const double r_rel_norm = std::sqrt(r_rel_sq);
    const double r_rel_cubed = r_rel_sq * r_rel_norm;

    const math::Vector3 direct_term = r_rel * (third_body_mu_m3_per_s2 / r_rel_cubed);
    const math::Vector3 indirect_term = third_body_position_eci_m * (third_body_mu_m3_per_s2 / r3_cubed);

    return direct_term - indirect_term;
}

}  // namespace astradock::environment
