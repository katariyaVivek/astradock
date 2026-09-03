#pragma once

// AstraDock M24 — High-fidelity environment: Earth rotation, geodetic frame,
// zonal gravity, and analytical third-body ephemeris.
//
// Physical problem:
//   Through M23 the Earth is a spherical point mass in an idealized inertial
//   frame with a fixed Moon direction and an exponential atmosphere. That was
//   the right staging — dynamics, estimation, and GNC first. M24 upgrades the
//   environment models that later missions (and the M25 capstone) build on,
//   each verified independently BEFORE composition:
//
//   M24A ECI/ECEF: explicit Earth rotation about the inertial +Z axis by
//     theta(t) = omega_E * t + theta_0 (documented time assumption below).
//   M24B Time: simulation elapsed (s) vs Julian-date civil reference vs
//     dynamical rate. No leap-second claims; the mapping is linear + offset.
//   M24C WGS-84 geodetic: ECEF <-> (lat, lon, alt) on the oblate ellipsoid
//     (Bowring iteration, polar/equatorial edge cases tested).
//   M24D Zonal gravity J3/J4: next Legendre terms beyond J2, same gradient
//     machinery, verified against finite differences of the potential.
//   M24E Ephemeris: analytical circular coplanar Moon/Sun directions with
//     documented mean elements (phase + rate), replacing the fixed vector.
//   M24F Atmosphere: retained exponential model with documented validity —
//     explicitly NOT claimed as operational (no NRLMSISE upgrade without a
//     justifiable source; the lesson records why).
//
// Time assumption (M24B, stated once, used everywhere): the simulator clock
// t_sim (s) maps to Earth rotation angle theta = omega_E * t_sim + theta_0
// with theta_0 = 0 (Greenwich meridian aligned with inertial +X at t = 0).
// Julian date JD = JD_REF + t_sim / 86400 with JD_REF = 2451545.0 (J2000).
// This is a SIMULATION time system, not UTC/TT: no leap seconds, no EOP.

#include "math/constants.hpp"
#include "math/matrix3.hpp"
#include "math/vector3.hpp"

#include <cmath>
#include <stdexcept>

namespace astradock::environment {

namespace detail {

inline void require_finite_m24(const math::Vector3& v, const char* message) {
    if (!math::is_finite(v)) {
        throw std::domain_error(message);
    }
}

}  // namespace detail

// WGS-84 ellipsoid shape constants (documented here, not scattered).
inline constexpr double k_wgs84_flattening = 1.0 / 298.257223563;
inline constexpr double k_wgs84_semimajor_m = 6378137.0;
inline constexpr double k_j2_ref = 1.08262668e-3;
inline constexpr double k_j3_ref = -2.532656485e-6;   // J3 (pear-shape, WGS-84/EGM)
inline constexpr double k_j4_ref = -1.619621591e-6;   // J4 (EGM96-class value)
inline constexpr double k_j2000_jd = 2451545.0;
inline constexpr double k_seconds_per_day = 86400.0;
// Mean lunar/solar ephemeris rates (circular coplanar analytical model).
inline constexpr double k_moon_distance_m = 384400.0e3;
inline constexpr double k_moon_period_s = 27.321661 * k_seconds_per_day;
inline constexpr double k_sun_distance_m = 149597870700.0;
inline constexpr double k_sun_period_s = 365.256363004 * k_seconds_per_day;

// M24B: simulation time -> Earth rotation angle (rad), theta_0 = 0.
[[nodiscard]] inline double earth_rotation_angle_rad(double sim_time_s) {
    if (!std::isfinite(sim_time_s)) {
        throw std::domain_error("Simulation time must be finite");
    }
    return constants::earth_rotation_rate_rad_per_s * sim_time_s;
}

// M24B: simulation time -> Julian date (civil reference, linear mapping).
[[nodiscard]] inline double sim_time_to_julian_date(double sim_time_s) {
    if (!std::isfinite(sim_time_s)) {
        throw std::domain_error("Simulation time must be finite");
    }
    return k_j2000_jd + sim_time_s / k_seconds_per_day;
}

// M24A: ECI -> ECEF rotation (principal rotation about +Z by theta).
[[nodiscard]] inline math::Matrix3 dcm_ecef_from_eci(double sim_time_s) {
    const double theta = earth_rotation_angle_rad(sim_time_s);
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    return math::Matrix3(c, s, 0.0, -s, c, 0.0, 0.0, 0.0, 1.0);
}

[[nodiscard]] inline math::Vector3 eci_to_ecef(
    const math::Vector3& position_eci_m,
    double sim_time_s) {
    detail::require_finite_m24(position_eci_m, "ECI position must be finite");
    return dcm_ecef_from_eci(sim_time_s) * position_eci_m;
}

[[nodiscard]] inline math::Vector3 ecef_to_eci(
    const math::Vector3& position_ecef_m,
    double sim_time_s) {
    detail::require_finite_m24(position_ecef_m, "ECEF position must be finite");
    return dcm_ecef_from_eci(sim_time_s).transpose() * position_ecef_m;
}

// M24C: geodetic coordinates on the WGS-84 ellipsoid.
struct GeodeticCoord {
    double latitude_rad{0.0};   // [-pi/2, pi/2]
    double longitude_rad{0.0};  // [-pi, pi]
    double altitude_m{0.0};     // above ellipsoid
};

// ECEF -> geodetic via Bowring iteration (converges for all non-polar points;
// polar axis handled by closed form).
[[nodiscard]] inline GeodeticCoord ecef_to_geodetic(const math::Vector3& ecef_m) {
    detail::require_finite_m24(ecef_m, "ECEF position must be finite");
    const double a = k_wgs84_semimajor_m;
    const double f = k_wgs84_flattening;
    const double e2 = 2.0 * f - f * f;
    const double x = ecef_m.x();
    const double y = ecef_m.y();
    const double z = ecef_m.z();
    const double p = std::sqrt(x * x + y * y);
    GeodeticCoord geo;
    geo.longitude_rad = std::atan2(y, x);
    if (p < 1.0) {
        // Polar axis: latitude is +/-90 deg, altitude from polar radius.
        geo.latitude_rad = (z >= 0.0 ? 1.0 : -1.0) * constants::pi / 2.0;
        geo.altitude_m = std::abs(z) - a * (1.0 - f);
        return geo;
    }
    double lat = std::atan2(z, p * (1.0 - e2));  // first guess
    for (int i = 0; i < 8; ++i) {
        const double sin_lat = std::sin(lat);
        const double prime = a / std::sqrt(1.0 - e2 * sin_lat * sin_lat);
        lat = std::atan2(z + e2 * prime * sin_lat, p);
    }
    const double sin_lat = std::sin(lat);
    const double prime = a / std::sqrt(1.0 - e2 * sin_lat * sin_lat);
    geo.latitude_rad = lat;
    geo.altitude_m = p / std::cos(lat) - prime;
    return geo;
}

// Geodetic -> ECEF (closed form).
[[nodiscard]] inline math::Vector3 geodetic_to_ecef(const GeodeticCoord& geo) {
    if (!std::isfinite(geo.latitude_rad) || !std::isfinite(geo.longitude_rad)
        || !std::isfinite(geo.altitude_m)) {
        throw std::domain_error("Geodetic coordinates must be finite");
    }
    if (std::abs(geo.latitude_rad) > constants::pi / 2.0) {
        throw std::domain_error("Geodetic latitude must lie within [-pi/2, pi/2]");
    }
    const double a = k_wgs84_semimajor_m;
    const double f = k_wgs84_flattening;
    const double e2 = 2.0 * f - f * f;
    const double sin_lat = std::sin(geo.latitude_rad);
    const double cos_lat = std::cos(geo.latitude_rad);
    const double prime = a / std::sqrt(1.0 - e2 * sin_lat * sin_lat);
    const double r_axial = (prime + geo.altitude_m) * cos_lat;
    return math::Vector3{
        r_axial * std::cos(geo.longitude_rad),
        r_axial * std::sin(geo.longitude_rad),
        ((1.0 - e2) * prime + geo.altitude_m) * sin_lat,
    };
}

// M24D: zonal acceleration (TOTAL: central + zonal perturbation) through J4.
// Potential (unnormalized):
//   U = mu/r [1 - J2 (R/r)^2 P2 - J3 (R/r)^3 P3 - J4 (R/r)^4 P4], s = z/r.
// Acceleration = grad U via the chain rule in (r, s):
//   a = (dU/dr) r_hat + (1/r)(dU/ds) (z_hat - s r_hat).
// Returns TOTAL acceleration (central -mu/r^2 term INCLUDED): this is what
// propagators integrate. Contrast the legacy M11 j2_acceleration_eci, which
// returns the PERTURBATION only — subtract the central term to compare.
// J2-only input reproduces (central + legacy J2) exactly (tested).
[[nodiscard]] inline math::Vector3 zonal_acceleration_eci(
    const math::Vector3& position_eci_m,
    double mu_m3_per_s2 = constants::earth_gravitational_parameter_m3_per_s2,
    double radius_m = constants::earth_reference_radius_m,
    double j2 = k_j2_ref,
    double j3 = 0.0,
    double j4 = 0.0) {
    detail::require_finite_m24(position_eci_m, "Position must be finite for zonal gravity");
    if (!std::isfinite(mu_m3_per_s2) || mu_m3_per_s2 <= 0.0) {
        throw std::domain_error("Gravitational parameter must be finite and positive");
    }
    if (!std::isfinite(radius_m) || radius_m <= 0.0) {
        throw std::domain_error("Reference radius must be finite and positive");
    }
    if (!std::isfinite(j2) || !std::isfinite(j3) || !std::isfinite(j4)) {
        throw std::domain_error("Zonal coefficients must be finite");
    }
    const double r2 = position_eci_m.squared_norm();
    if (r2 <= 0.0) {
        throw std::domain_error("Zonal gravity is singular at the origin");
    }
    const double r = std::sqrt(r2);
    const double s = position_eci_m.z() / r;
    const double rho = radius_m / r;
    // Legendre polynomials and derivatives in s.
    const double p2 = 0.5 * (3.0 * s * s - 1.0);
    const double p3 = 0.5 * (5.0 * s * s * s - 3.0 * s);
    const double p4 = 0.125 * (35.0 * s * s * s * s - 30.0 * s * s + 3.0);
    const double dp2 = 3.0 * s;
    const double dp3 = 0.5 * (15.0 * s * s - 3.0);
    const double dp4 = 0.125 * (140.0 * s * s * s - 60.0 * s);
    const double rho2 = rho * rho;
    const double rho3 = rho2 * rho;
    const double rho4 = rho3 * rho;
    // U = mu/r * S, S = 1 - (J2 rho^2 P2 + J3 rho^3 P3 + J4 rho^4 P4).
    // dS/dr = (2 J2 rho^2 P2 + 3 J3 rho^3 P3 + 4 J4 rho^4 P4)/r.
    const double dS_dr =
        (2.0 * j2 * rho2 * p2 + 3.0 * j3 * rho3 * p3 + 4.0 * j4 * rho4 * p4) / r;
    const double S = 1.0 - (j2 * rho2 * p2 + j3 * rho3 * p3 + j4 * rho4 * p4);
    const double dU_dr = -mu_m3_per_s2 / (r * r) * S + mu_m3_per_s2 / r * dS_dr;
    // dU/ds = mu/r * dS/ds, dS/ds = -(J2 rho^2 dP2 + ...).
    const double dS_ds = -(j2 * rho2 * dp2 + j3 * rho3 * dp3 + j4 * rho4 * dp4);
    const double dU_ds = mu_m3_per_s2 / r * dS_ds;
    const math::Vector3 r_hat = position_eci_m / r;
    const math::Vector3 z_hat{0.0, 0.0, 1.0};
    const math::Vector3 grad_s_dir = z_hat - r_hat * s;
    return r_hat * dU_dr + grad_s_dir * (dU_ds / r);
}

// M24E: analytical circular coplanar ephemeris. Moon and Sun orbit in the
// equatorial plane (documented simplification: real lunar inclination 5.14 deg
// and eccentricity 0.055 are NOT modeled) with fixed phase at t = 0.
[[nodiscard]] inline math::Vector3 moon_position_eci_m(double sim_time_s, double phase_rad = 0.0) {
    if (!std::isfinite(sim_time_s) || !std::isfinite(phase_rad)) {
        throw std::domain_error("Ephemeris time and phase must be finite");
    }
    const double angle = phase_rad + 2.0 * constants::pi * sim_time_s / k_moon_period_s;
    return math::Vector3{
        k_moon_distance_m * std::cos(angle), k_moon_distance_m * std::sin(angle), 0.0};
}

[[nodiscard]] inline math::Vector3 sun_position_eci_m(double sim_time_s, double phase_rad = 0.0) {
    if (!std::isfinite(sim_time_s) || !std::isfinite(phase_rad)) {
        throw std::domain_error("Ephemeris time and phase must be finite");
    }
    const double angle = phase_rad + 2.0 * constants::pi * sim_time_s / k_sun_period_s;
    return math::Vector3{
        k_sun_distance_m * std::cos(angle), k_sun_distance_m * std::sin(angle), 0.0};
}

}  // namespace astradock::environment
