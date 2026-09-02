#pragma once

#include "math/constants.hpp"
#include "math/vector3.hpp"

#include <cmath>
#include <stdexcept>

namespace astradock::environment {

// Computes the velocity of a spacecraft relative to Earth's co-rotating atmosphere
// in Earth-Centered Inertial (ECI) coordinates.
//
// Physical Model:
//   Earth rotates at angular rate omega_E about its inertial +Z axis:
//     omega_vector = [0, 0, omega_E]^T
//   The atmospheric wind velocity in ECI is:
//     v_atm = omega_vector x r_ECI = [-omega_E * y, omega_E * x, 0]^T
//   The spacecraft velocity relative to the local rotating atmosphere is:
//     v_rel = v_ECI - v_atm = [v_x + omega_E * y, v_y - omega_E * x, v_z]^T
[[nodiscard]] inline math::Vector3 relative_atmospheric_velocity_eci(
    const math::Vector3& velocity_eci_mps,
    const math::Vector3& position_eci_m,
    double earth_rotation_rate_rad_per_s = constants::earth_rotation_rate_rad_per_s) {
    if (!math::is_finite(velocity_eci_mps) || !math::is_finite(position_eci_m)) {
        throw std::domain_error("Velocity and position must be finite for relative atmospheric velocity");
    }
    if (!std::isfinite(earth_rotation_rate_rad_per_s)) {
        throw std::domain_error("Earth rotation rate must be finite");
    }

    const double v_rel_x = velocity_eci_mps.x() + earth_rotation_rate_rad_per_s * position_eci_m.y();
    const double v_rel_y = velocity_eci_mps.y() - earth_rotation_rate_rad_per_s * position_eci_m.x();
    const double v_rel_z = velocity_eci_mps.z();

    return math::Vector3{v_rel_x, v_rel_y, v_rel_z};
}

// Computes atmospheric density using an educational exponential scale-height model.
//
// Physical Model:
//   rho(h) = rho_0 * exp( - (h - h_0) / H )
//
// Default LEO Reference Parameters:
//   - h_0: 500 km (500,000 m) reference altitude
//   - rho_0: 6.967e-13 kg/m^3 (approximate mean solar-activity density at 500 km)
//   - H: 63.8 km (63,800 m) atmospheric scale height
[[nodiscard]] inline double exponential_atmospheric_density(
    double altitude_m,
    double reference_altitude_m = 500.0e3,
    double reference_density_kg_per_m3 = 6.967e-13,
    double scale_height_m = 63.8e3) {
    if (!std::isfinite(altitude_m) || !std::isfinite(reference_altitude_m)
        || !std::isfinite(reference_density_kg_per_m3) || !std::isfinite(scale_height_m)) {
        throw std::domain_error("All atmospheric density parameters must be finite");
    }
    if (reference_density_kg_per_m3 < 0.0) {
        throw std::domain_error("Reference atmospheric density must be non-negative");
    }
    if (scale_height_m <= 0.0) {
        throw std::domain_error("Atmospheric scale height must be strictly positive");
    }

    // Clamp effective altitude to zero if below Earth reference surface
    const double effective_alt = std::max(altitude_m, 0.0);
    const double exponent = -(effective_alt - reference_altitude_m) / scale_height_m;

    // Prevent floating-point underflow / overflow
    if (exponent < -100.0) {
        return 0.0;
    }
    if (exponent > 100.0) {
        return reference_density_kg_per_m3 * std::exp(100.0);
    }

    return reference_density_kg_per_m3 * std::exp(exponent);
}

// Computes aerodynamic drag acceleration in ECI coordinates.
//
// Physical Model:
//   a_drag = - 0.5 * C_D * (A / m) * rho(h) * |v_rel| * v_rel
//
// Units:
//   - velocity_eci_mps: Spacecraft inertial velocity (m/s)
//   - position_eci_m: Spacecraft inertial position (m)
//   - mass_kg: Spacecraft total mass (kg)
//   - drag_coefficient_cd: Dimensionless drag coefficient C_D (typically ~2.0 - 2.2 in free-molecular flow)
//   - drag_reference_area_m2: Cross-sectional drag area A (m^2)
//   - earth_reference_radius_m: Earth reference radius R_E (m)
//   - earth_rotation_rate_rad_per_s: Earth rotation rate omega_E (rad/s)
//   - Returns: Drag acceleration in ECI frame (m/s^2)
[[nodiscard]] inline math::Vector3 drag_acceleration_eci(
    const math::Vector3& velocity_eci_mps,
    const math::Vector3& position_eci_m,
    double mass_kg,
    double drag_coefficient_cd,
    double drag_reference_area_m2,
    double earth_reference_radius_m = constants::earth_reference_radius_m,
    double earth_rotation_rate_rad_per_s = constants::earth_rotation_rate_rad_per_s,
    double reference_altitude_m = 500.0e3,
    double reference_density_kg_per_m3 = 6.967e-13,
    double scale_height_m = 63.8e3) {
    if (!math::is_finite(velocity_eci_mps) || !math::is_finite(position_eci_m)) {
        throw std::domain_error("Velocity and position must be finite for drag acceleration");
    }
    if (!std::isfinite(mass_kg) || mass_kg <= 0.0) {
        throw std::domain_error("Spacecraft mass must be finite and strictly positive for drag acceleration");
    }
    if (!std::isfinite(drag_coefficient_cd) || drag_coefficient_cd < 0.0) {
        throw std::domain_error("Drag coefficient must be finite and non-negative");
    }
    if (!std::isfinite(drag_reference_area_m2) || drag_reference_area_m2 < 0.0) {
        throw std::domain_error("Drag reference area must be finite and non-negative");
    }

    const double r_norm = position_eci_m.norm();
    if (r_norm <= 0.0) {
        throw std::domain_error("Spacecraft position cannot be zero for drag calculation");
    }

    const double altitude_m = r_norm - earth_reference_radius_m;
    const double density = exponential_atmospheric_density(
        altitude_m,
        reference_altitude_m,
        reference_density_kg_per_m3,
        scale_height_m
    );

    if (density <= 0.0 || drag_coefficient_cd == 0.0 || drag_reference_area_m2 == 0.0) {
        return math::Vector3{0.0, 0.0, 0.0};
    }

    const math::Vector3 v_rel = relative_atmospheric_velocity_eci(
        velocity_eci_mps,
        position_eci_m,
        earth_rotation_rate_rad_per_s
    );

    const double v_rel_norm = v_rel.norm();
    if (v_rel_norm == 0.0) {
        return math::Vector3{0.0, 0.0, 0.0};
    }

    const double ballistic_factor = 0.5 * drag_coefficient_cd * (drag_reference_area_m2 / mass_kg);
    const double drag_mag_factor = -ballistic_factor * density * v_rel_norm;

    return v_rel * drag_mag_factor;
}

}  // namespace astradock::environment
