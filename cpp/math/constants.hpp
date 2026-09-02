#pragma once

#include <numbers>

namespace astradock::constants {

inline constexpr double pi = std::numbers::pi_v<double>;

// WGS 84 conventional Earth gravitational parameter, including the mass of
// the atmosphere. Units: m^3 / s^2.
inline constexpr double earth_gravitational_parameter_m3_per_s2 = 3.986004418e14;

// WGS 84 semi-major axis. AstraDock uses this as the spherical Earth reference
// radius for the introductory analytical examples. Units: m.
inline constexpr double earth_reference_radius_m = 6'378'137.0;

// WGS 84 second zonal harmonic (Earth oblateness coefficient, dimensionless).
inline constexpr double earth_j2 = 1.08262668e-3;

// Nominal Earth rotation rate relative to inertial space (WGS 84). Units: rad/s.
// (7.2921150e-5 rad/s corresponds to 1 sidereal day of 86164.0905 s).
inline constexpr double earth_rotation_rate_rad_per_s = 7.2921150e-5;

// Canonical lunar gravitational parameter (Moon). Units: m^3 / s^2.
inline constexpr double moon_gravitational_parameter_m3_per_s2 = 4.9048695e12;

// Canonical solar gravitational parameter (Sun). Units: m^3 / s^2.
inline constexpr double sun_gravitational_parameter_m3_per_s2 = 1.32712440018e20;

}  // namespace astradock::constants
