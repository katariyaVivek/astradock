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

}  // namespace astradock::constants
