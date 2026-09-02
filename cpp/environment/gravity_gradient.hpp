#pragma once

#include "attitude/attitude_state.hpp"
#include "attitude/principal_inertia.hpp"
#include "math/constants.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"

#include <cmath>
#include <stdexcept>

namespace astradock::environment {

// Computes the gravity-gradient environmental torque acting on a rigid spacecraft
// expressed in the spacecraft BODY frame.
//
// Physical Model:
//   Earth's gravitational field varies as 1/r^2 across the finite spatial dimensions
//   of a non-spherical spacecraft. The differential gravitational attraction across
//   mass elements dm produces a net restorative/destabilizing torque about the center of mass:
//
//     tau_gg = (3 * mu / r^3) * r_hat_B x (I * r_hat_B)
//
// Frame Conventions:
//   - position_eci_m: Spacecraft center-of-mass position in ECI frame (m)
//   - orientation_eci_from_body: Unit quaternion q_(I_B) mapping Body coordinates to ECI
//   - inertia: Principal moments of inertia I = diag(I_xx, I_yy, I_zz) in Body frame (kg*m^2)
//   - mu_m3_per_s2: Central body gravitational parameter (m^3/s^2)
//   - Returns: Environmental torque in spacecraft BODY frame (N*m)
//
// Derivation of Body-Frame Radial Vector:
//   1. ECI radial unit vector: r_hat_I = r_ECI / |r_ECI|
//   2. Body radial unit vector: r_hat_B = C_(B_I) * r_hat_I = C_(I_B)^T * r_hat_I
//      Using quaternion algebra: r_hat_B = rotate_vector_inertial_to_body(q, r_hat_I)
//                                       = q.conjugate().rotate(r_hat_I)
//
// Component-Wise Form in Principal Axes (r_hat_B = [u_x, u_y, u_z]^T):
//   tau_gg_x = (3 * mu / r^3) * (I_zz - I_yy) * u_y * u_z
//   tau_gg_y = (3 * mu / r^3) * (I_xx - I_zz) * u_x * u_z
//   tau_gg_z = (3 * mu / r^3) * (I_yy - I_xx) * u_x * u_y
[[nodiscard]] inline math::Vector3 gravity_gradient_torque_body(
    const math::Vector3& position_eci_m,
    const math::Quaternion& orientation_eci_from_body,
    const attitude::PrincipalInertia& inertia,
    double mu_m3_per_s2 = constants::earth_gravitational_parameter_m3_per_s2) {
    if (!math::is_finite(position_eci_m)) {
        throw std::domain_error("Spacecraft position must contain only finite values for gravity-gradient torque");
    }
    if (!math::is_finite(orientation_eci_from_body)) {
        throw std::domain_error("Spacecraft attitude quaternion must contain only finite values");
    }
    if (!orientation_eci_from_body.is_unit(1.0e-3)) {
        throw std::domain_error("Attitude quaternion must have unit norm for gravity-gradient calculation");
    }
    if (!attitude::is_valid(inertia)) {
        throw std::domain_error("Spacecraft principal inertia moments must be finite and strictly positive");
    }
    if (!std::isfinite(mu_m3_per_s2) || mu_m3_per_s2 <= 0.0) {
        throw std::domain_error("Gravitational parameter must be finite and strictly positive");
    }

    const double r2 = position_eci_m.squared_norm();
    if (r2 <= 0.0) {
        throw std::domain_error("Spacecraft position cannot be zero for gravity-gradient torque");
    }

    const double r = std::sqrt(r2);
    const double r3 = r2 * r;

    // Unit radial vector in ECI
    const math::Vector3 r_hat_eci = position_eci_m / r;

    // Transform ECI radial unit vector to Spacecraft Body frame
    // Active rotation q maps Body to ECI, so q.conjugate() maps ECI to Body
    const math::Vector3 r_hat_body = orientation_eci_from_body.conjugate().rotate_vector(r_hat_eci);

    const double ux = r_hat_body.x();
    const double uy = r_hat_body.y();
    const double uz = r_hat_body.z();

    const double factor = 3.0 * mu_m3_per_s2 / r3;

    const double tau_x = factor * (inertia.Izz() - inertia.Iyy()) * uy * uz;
    const double tau_y = factor * (inertia.Ixx() - inertia.Izz()) * ux * uz;
    const double tau_z = factor * (inertia.Iyy() - inertia.Ixx()) * ux * uy;

    return math::Vector3{tau_x, tau_y, tau_z};
}

}  // namespace astradock::environment
