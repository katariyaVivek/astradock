#pragma once

#include "math/constants.hpp"
#include "math/vector3.hpp"

#include <cmath>
#include <stdexcept>

namespace astradock::environment {

// Computes the gravitational perturbation acceleration due to Earth's second zonal
// harmonic (J2 oblateness) in Earth-Centered Inertial (ECI) coordinates.
//
// Physical Model:
//   Earth is an oblate spheroid with an equatorial bulge due to planetary rotation.
//   The gravitational potential expanded in Legendre polynomials is:
//     U(r, phi) = (mu / r) * [ 1 - J2 * (R_E / r)^2 * (3 * sin^2(phi) - 1) / 2 ]
//   where sin(phi) = z / r is the geocentric latitude.
//
//   Taking the gradient of the J2 perturbing potential produces the acceleration:
//     a_J2 = - (3 * mu * J2 * R_E^2) / (2 * r^5) * [
//              x * (1 - 5 * z^2 / r^2),
//              y * (1 - 5 * z^2 / r^2),
//              z * (3 - 5 * z^2 / r^2)
//            ]^T
//
// Frame and Units:
//   - position_eci_m: Position vector in ECI frame (m)
//   - mu_m3_per_s2: Gravitational parameter (m^3/s^2)
//   - earth_reference_radius_m: Equatorial reference radius R_E (m)
//   - j2: Dimensionless second zonal harmonic coefficient (e.g. WGS 84 J2 ~ 1.08262668e-3)
//   - Returns: Perturbation acceleration in ECI frame (m/s^2)
//
// Special Cases & Physical Checks:
//   - Equatorial plane (z = 0): a_J2 = - (3 * mu * J2 * R_E^2) / (2 * r^4) * r_hat (inward radial)
//   - North/South Pole (x = y = 0, z = r): a_J2 = + (3 * mu * J2 * R_E^2) / r^4 * z_hat (outward/equatorward)
//   - Critical latitude (z^2 / r^2 = 1/5 => phi ~ 26.565°): equatorial perturbation term vanishes.
[[nodiscard]] inline math::Vector3 j2_acceleration_eci(
    const math::Vector3& position_eci_m,
    double mu_m3_per_s2 = constants::earth_gravitational_parameter_m3_per_s2,
    double earth_reference_radius_m = constants::earth_reference_radius_m,
    double j2 = constants::earth_j2) {
    if (!math::is_finite(position_eci_m)) {
        throw std::domain_error("Spacecraft position must contain only finite values for J2 acceleration");
    }
    if (!std::isfinite(mu_m3_per_s2) || mu_m3_per_s2 <= 0.0) {
        throw std::domain_error("Gravitational parameter must be finite and strictly positive for J2 acceleration");
    }
    if (!std::isfinite(earth_reference_radius_m) || earth_reference_radius_m <= 0.0) {
        throw std::domain_error("Earth reference radius must be finite and strictly positive for J2 acceleration");
    }
    if (!std::isfinite(j2)) {
        throw std::domain_error("J2 coefficient must be finite for J2 acceleration");
    }

    const double r2 = position_eci_m.squared_norm();
    if (r2 <= 0.0) {
        throw std::domain_error("Spacecraft position cannot be zero for J2 acceleration (singularity at center)");
    }

    const double r = std::sqrt(r2);
    const double r5 = r2 * r2 * r;
    const double z2_over_r2 = (position_eci_m.z() * position_eci_m.z()) / r2;

    const double factor = -1.5 * mu_m3_per_s2 * j2 * (earth_reference_radius_m * earth_reference_radius_m) / r5;

    const double ax = factor * position_eci_m.x() * (1.0 - 5.0 * z2_over_r2);
    const double ay = factor * position_eci_m.y() * (1.0 - 5.0 * z2_over_r2);
    const double az = factor * position_eci_m.z() * (3.0 - 5.0 * z2_over_r2);

    return math::Vector3{ax, ay, az};
}

}  // namespace astradock::environment
