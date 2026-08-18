#pragma once

#include "dynamics/two_body.hpp"
#include "math/angle.hpp"
#include "math/constants.hpp"
#include "math/matrix3.hpp"
#include "math/vector3.hpp"
#include "orbit/cartesian_state.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace astradock::orbit {

// Classical Keplerian Orbital Elements (Keplerian element set).
//
// These six parameters describe the size, shape, orientation, and instantaneous
// position of a spacecraft in an idealized two-body Keplerian orbit relative
// to an Earth-Centered Inertial (ECI) coordinate frame.
//
// Units:
//   semi_major_axis_m:          metres (m)
//   eccentricity:               dimensionless
//   inclination_rad:            radians (rad) in [0, pi]
//   raan_rad:                   radians (rad) in [0, 2pi)
//   argument_of_periapsis_rad:  radians (rad) in [0, 2pi)
//   true_anomaly_rad:           radians (rad) in [0, 2pi)
struct ClassicalOrbitalElements {
    double semi_major_axis_m;
    double eccentricity;
    double inclination_rad;
    double raan_rad;
    double argument_of_periapsis_rad;
    double true_anomaly_rad;
};

// Numerical classification for orbit shape.
enum class OrbitShape {
    circular,         // e < eccentricity_tolerance
    elliptic,         // eccentricity_tolerance <= e < 1.0
    parabolic_limit,  // e ≈ 1.0 (or epsilon ≈ 0)
    hyperbolic,       // e > 1.0 (epsilon > 0)
};

// Numerical classification for orbit inclination.
enum class OrbitInclination {
    equatorial_prograde,    // i ≈ 0
    inclined,               // 0 < i < pi/2 or pi/2 < i < pi
    polar,                  // i ≈ pi/2
    equatorial_retrograde,  // i ≈ pi
};

struct OrbitClassification {
    OrbitShape shape;
    OrbitInclination inclination;
    bool is_circular;
    bool is_equatorial;
};

// Centralized default numerical singularity thresholds
inline constexpr double default_eccentricity_tolerance = 1.0e-11;
inline constexpr double default_node_norm_tolerance = 1.0e-11;
inline constexpr double default_inclination_tolerance_rad = 1.0e-11;

// ---------------------------------------------------------------------------
// Analytical Auxiliary Functions
// ---------------------------------------------------------------------------

// Computes the semi-latus rectum p = a * (1 - e^2) in metres.
[[nodiscard]] inline double semi_latus_rectum_m(
    double semi_major_axis_m,
    double eccentricity) {
    if (!std::isfinite(semi_major_axis_m) || semi_major_axis_m <= 0.0) {
        throw std::domain_error("Semi-major axis must be finite and positive for bound orbits");
    }
    if (!std::isfinite(eccentricity) || eccentricity < 0.0 || eccentricity >= 1.0) {
        throw std::domain_error("Eccentricity must be in [0, 1) for elliptic orbits");
    }
    return semi_major_axis_m * (1.0 - eccentricity * eccentricity);
}

// Computes the periapsis radius r_p = a * (1 - e) in metres.
[[nodiscard]] inline double periapsis_radius_m(
    double semi_major_axis_m,
    double eccentricity) {
    if (!std::isfinite(semi_major_axis_m) || semi_major_axis_m <= 0.0) {
        throw std::domain_error("Semi-major axis must be finite and positive");
    }
    if (!std::isfinite(eccentricity) || eccentricity < 0.0 || eccentricity >= 1.0) {
        throw std::domain_error("Eccentricity must be in [0, 1) for elliptic orbits");
    }
    return semi_major_axis_m * (1.0 - eccentricity);
}

// Computes the apoapsis radius r_a = a * (1 + e) in metres.
[[nodiscard]] inline double apoapsis_radius_m(
    double semi_major_axis_m,
    double eccentricity) {
    if (!std::isfinite(semi_major_axis_m) || semi_major_axis_m <= 0.0) {
        throw std::domain_error("Semi-major axis must be finite and positive");
    }
    if (!std::isfinite(eccentricity) || eccentricity < 0.0 || eccentricity >= 1.0) {
        throw std::domain_error("Eccentricity must be in [0, 1) for elliptic orbits");
    }
    return semi_major_axis_m * (1.0 + eccentricity);
}

// Computes mean motion n = sqrt(mu / a^3) in radians per second.
[[nodiscard]] inline double mean_motion_rad_per_s(
    double semi_major_axis_m,
    double gravitational_parameter_m3_per_s2) {
    if (!std::isfinite(gravitational_parameter_m3_per_s2) || gravitational_parameter_m3_per_s2 <= 0.0) {
        throw std::domain_error("Gravitational parameter must be finite and positive");
    }
    if (!std::isfinite(semi_major_axis_m) || semi_major_axis_m <= 0.0) {
        throw std::domain_error("Semi-major axis must be finite and positive for bound orbits");
    }
    const double a3 = semi_major_axis_m * semi_major_axis_m * semi_major_axis_m;
    return std::sqrt(gravitational_parameter_m3_per_s2 / a3);
}

// Computes the orbital period T = 2 * pi * sqrt(a^3 / mu) in seconds.
[[nodiscard]] inline double orbital_period_s(
    double semi_major_axis_m,
    double gravitational_parameter_m3_per_s2) {
    const double n = mean_motion_rad_per_s(semi_major_axis_m, gravitational_parameter_m3_per_s2);
    return 2.0 * constants::pi / n;
}

// Classifies an orbit given its classical elements.
[[nodiscard]] inline OrbitClassification classify_orbit(
    const ClassicalOrbitalElements& elements,
    double eccentricity_tolerance = default_eccentricity_tolerance,
    double inclination_tolerance_rad = default_inclination_tolerance_rad) noexcept {
    OrbitShape shape = OrbitShape::elliptic;
    if (elements.eccentricity < eccentricity_tolerance) {
        shape = OrbitShape::circular;
    } else if (std::abs(elements.eccentricity - 1.0) < eccentricity_tolerance) {
        shape = OrbitShape::parabolic_limit;
    } else if (elements.eccentricity > 1.0) {
        shape = OrbitShape::hyperbolic;
    }

    OrbitInclination inc = OrbitInclination::inclined;
    if (elements.inclination_rad < inclination_tolerance_rad) {
        inc = OrbitInclination::equatorial_prograde;
    } else if (std::abs(elements.inclination_rad - constants::pi * 0.5) < inclination_tolerance_rad) {
        inc = OrbitInclination::polar;
    } else if (std::abs(elements.inclination_rad - constants::pi) < inclination_tolerance_rad) {
        inc = OrbitInclination::equatorial_retrograde;
    }

    const bool is_circ = (shape == OrbitShape::circular);
    const bool is_eq = (inc == OrbitInclination::equatorial_prograde || inc == OrbitInclination::equatorial_retrograde);

    return {shape, inc, is_circ, is_eq};
}

// ---------------------------------------------------------------------------
// State -> Classical Orbital Elements Conversion
// ---------------------------------------------------------------------------

// Converts a Cartesian position and velocity state (in ECI) into classical
// Keplerian orbital elements.
//
// Inputs:
//   state:                             Cartesian state in ECI (position in m, velocity in m/s).
//   gravitational_parameter_m3_per_s2: Central body gravitational parameter mu (m^3/s^2).
//   eccentricity_tolerance:            Threshold below which an orbit is classified circular.
//   node_norm_tolerance:               Threshold below which an orbit is classified equatorial.
//
// Singularities and Conventions:
//   - Circular orbit (e < tol):
//       Argument of periapsis is undefined (convention: 0.0 rad).
//       True anomaly is measured from the ascending node (argument of latitude u), or from the
//       +X axis if also equatorial (true longitude lambda).
//   - Equatorial orbit (n < tol):
//       RAAN is undefined (convention: 0.0 rad).
//       Argument of periapsis is measured from the +X axis (longitude of periapsis varpi).
//   - Circular Equatorial orbit:
//       RAAN = 0.0 rad, Argument of Periapsis = 0.0 rad, True Anomaly = lambda (true longitude).
[[nodiscard]] inline ClassicalOrbitalElements state_to_classical_elements(
    const CartesianState& state,
    double gravitational_parameter_m3_per_s2,
    double eccentricity_tolerance = default_eccentricity_tolerance,
    double node_norm_tolerance = default_node_norm_tolerance) {
    if (!std::isfinite(gravitational_parameter_m3_per_s2) || gravitational_parameter_m3_per_s2 <= 0.0) {
        throw std::domain_error("State to elements: gravitational parameter must be finite and positive");
    }
    if (!math::is_finite(state.position) || !math::is_finite(state.velocity)) {
        throw std::domain_error("State to elements: position and velocity must contain only finite components");
    }

    const double r_mag = state.position.norm();
    if (r_mag == 0.0) {
        throw std::domain_error("State to elements: position magnitude cannot be zero");
    }

    const double v_mag = state.velocity.norm();

    // 1. Specific angular momentum vector: h = r x v
    const math::Vector3 h_vec = state.position.cross(state.velocity);
    const double h_mag = h_vec.norm();
    if (h_mag == 0.0) {
        throw std::domain_error("State to elements: angular momentum is zero (rectilinear/radial trajectory)");
    }
    const math::Vector3 h_hat = h_vec / h_mag;

    // 2. Specific orbital energy: epsilon = v^2 / 2 - mu / r
    const double specific_energy = 0.5 * v_mag * v_mag - gravitational_parameter_m3_per_s2 / r_mag;

    // For M07, require bound elliptical motion (epsilon < 0)
    if (specific_energy >= 0.0) {
        throw std::domain_error("State to elements: non-elliptic/unbound trajectory (parabolic/hyperbolic not supported in M07)");
    }

    // 3. Semi-major axis: a = -mu / (2 * epsilon)
    const double semi_major_axis = -gravitational_parameter_m3_per_s2 / (2.0 * specific_energy);

    // 4. Eccentricity vector: e = (v x h) / mu - r / |r|
    const math::Vector3 v_cross_h = state.velocity.cross(h_vec);
    const math::Vector3 e_vec = (v_cross_h / gravitational_parameter_m3_per_s2) - (state.position / r_mag);
    const double e_mag = e_vec.norm();

    // 5. Inclination: i = acos(h_z / h) in [0, pi]
    const double cos_i = std::clamp(h_vec.z() / h_mag, -1.0, 1.0);
    const double inclination = std::acos(cos_i);

    // 6. Ascending node vector: n = k x h = [-h_y, h_x, 0]
    const math::Vector3 n_vec(-h_vec.y(), h_vec.x(), 0.0);
    const double n_mag = n_vec.norm();

    const bool is_circular = (e_mag < eccentricity_tolerance);
    const bool is_equatorial = (n_mag < node_norm_tolerance);

    double raan = 0.0;
    double argument_of_periapsis = 0.0;
    double true_anomaly = 0.0;

    // 7. Right Ascension of Ascending Node (RAAN, Omega)
    if (!is_equatorial) {
        raan = math::normalize_angle_2pi_rad(std::atan2(n_vec.y(), n_vec.x()));
    } else {
        // Convention for equatorial orbit: RAAN = 0
        raan = 0.0;
    }

    // 8. Argument of Periapsis (omega)
    if (!is_circular && !is_equatorial) {
        // Standard non-singular case: angle between node vector n and eccentricity vector e
        const double sin_omega = n_vec.cross(e_vec).dot(h_hat);
        const double cos_omega = n_vec.dot(e_vec);
        argument_of_periapsis = math::normalize_angle_2pi_rad(std::atan2(sin_omega, cos_omega));
    } else if (!is_circular && is_equatorial) {
        // Equatorial eccentric orbit: longitude of periapsis varpi = atan2(e_y, e_x)
        const double sin_varpi = e_vec.y();
        const double cos_varpi = e_vec.x();
        argument_of_periapsis = math::normalize_angle_2pi_rad(std::atan2(sin_varpi, cos_varpi));
    } else {
        // Circular orbit (equatorial or inclined): argument of periapsis is undefined (convention: 0.0)
        argument_of_periapsis = 0.0;
    }

    // 9. True Anomaly (nu)
    if (!is_circular) {
        // Standard eccentric case: angle between eccentricity vector e and position vector r
        const double sin_nu = e_vec.cross(state.position).dot(h_hat);
        const double cos_nu = e_vec.dot(state.position);
        true_anomaly = math::normalize_angle_2pi_rad(std::atan2(sin_nu, cos_nu));
    } else if (!is_equatorial) {
        // Circular inclined orbit: argument of latitude u = angle between n and r
        const double sin_u = n_vec.cross(state.position).dot(h_hat);
        const double cos_u = n_vec.dot(state.position);
        true_anomaly = math::normalize_angle_2pi_rad(std::atan2(sin_u, cos_u));
    } else {
        // Circular equatorial orbit: true longitude lambda = atan2(r_y, r_x)
        true_anomaly = math::normalize_angle_2pi_rad(std::atan2(state.position.y(), state.position.x()));
    }

    return {
        semi_major_axis,
        e_mag,
        inclination,
        raan,
        argument_of_periapsis,
        true_anomaly,
    };
}

// ---------------------------------------------------------------------------
// Classical Orbital Elements -> State Conversion
// ---------------------------------------------------------------------------

// Converts classical Keplerian orbital elements into a Cartesian position and
// velocity state in ECI.
//
// The conversion constructs coordinates in the Perifocal Frame (PQW), where:
//   P points toward periapsis in the orbital plane,
//   W points along the orbital angular momentum vector (orbit normal),
//   Q completes the right-handed triad in the orbital plane (W x P).
//
// Perifocal coordinates are rotated into ECI coordinates via the standard sequence:
//   C_ECI_PQW = R3(-Omega) * R1(-i) * R3(-omega)
//
// The columns of C_ECI_PQW are the unit vectors P, Q, W expressed in ECI:
//   P = [cos(Omega)*cos(omega) - sin(Omega)*sin(omega)*cos(i),
//        sin(Omega)*cos(omega) + cos(Omega)*sin(omega)*cos(i),
//        sin(omega)*sin(i)]
//   Q = [-cos(Omega)*sin(omega) - sin(Omega)*cos(omega)*cos(i),
//        -sin(Omega)*sin(omega) + cos(Omega)*cos(omega)*cos(i),
//        cos(omega)*sin(i)]
//   W = [sin(Omega)*sin(i),
//        -cos(Omega)*sin(i),
//        cos(i)]
[[nodiscard]] inline CartesianState classical_elements_to_state(
    const ClassicalOrbitalElements& elements,
    double gravitational_parameter_m3_per_s2) {
    if (!std::isfinite(gravitational_parameter_m3_per_s2) || gravitational_parameter_m3_per_s2 <= 0.0) {
        throw std::domain_error("Elements to state: gravitational parameter must be finite and positive");
    }
    if (!std::isfinite(elements.semi_major_axis_m) || elements.semi_major_axis_m <= 0.0) {
        throw std::domain_error("Elements to state: semi-major axis must be finite and positive");
    }
    if (!std::isfinite(elements.eccentricity) || elements.eccentricity < 0.0 || elements.eccentricity >= 1.0) {
        throw std::domain_error("Elements to state: eccentricity must be in [0, 1) for elliptic orbits");
    }
    if (!std::isfinite(elements.inclination_rad)
        || !std::isfinite(elements.raan_rad)
        || !std::isfinite(elements.argument_of_periapsis_rad)
        || !std::isfinite(elements.true_anomaly_rad)) {
        throw std::domain_error("Elements to state: angular elements must contain only finite values");
    }

    const double p = semi_latus_rectum_m(elements.semi_major_axis_m, elements.eccentricity);
    const double cos_nu = std::cos(elements.true_anomaly_rad);
    const double sin_nu = std::sin(elements.true_anomaly_rad);

    const double denom = 1.0 + elements.eccentricity * cos_nu;
    if (denom <= 0.0) {
        throw std::domain_error("Elements to state: non-physical orbit denominator (1 + e*cos(nu) <= 0)");
    }
    const double r_orbit = p / denom;

    // Perifocal coordinates:
    // r_PQW = [r * cos(nu), r * sin(nu), 0]
    // v_PQW = [sqrt(mu / p) * (-sin(nu)), sqrt(mu / p) * (e + cos(nu)), 0]
    const double v_coeff = std::sqrt(gravitational_parameter_m3_per_s2 / p);

    const double r_p = r_orbit * cos_nu;
    const double r_q = r_orbit * sin_nu;

    const double v_p = v_coeff * (-sin_nu);
    const double v_q = v_coeff * (elements.eccentricity + cos_nu);

    // Orientation angles
    const double cos_raan = std::cos(elements.raan_rad);
    const double sin_raan = std::sin(elements.raan_rad);
    const double cos_inc = std::cos(elements.inclination_rad);
    const double sin_inc = std::sin(elements.inclination_rad);
    const double cos_argp = std::cos(elements.argument_of_periapsis_rad);
    const double sin_argp = std::sin(elements.argument_of_periapsis_rad);

    // Unit vectors P and Q expressed in ECI:
    const math::Vector3 P_vec(
        cos_raan * cos_argp - sin_raan * sin_argp * cos_inc,
        sin_raan * cos_argp + cos_raan * sin_argp * cos_inc,
        sin_argp * sin_inc
    );

    const math::Vector3 Q_vec(
        -cos_raan * sin_argp - sin_raan * cos_argp * cos_inc,
        -sin_raan * sin_argp + cos_raan * cos_argp * cos_inc,
        cos_argp * sin_inc
    );

    // Transform Perifocal vectors to ECI
    const math::Vector3 position_eci = P_vec * r_p + Q_vec * r_q;
    const math::Vector3 velocity_eci = P_vec * v_p + Q_vec * v_q;

    return {position_eci, velocity_eci};
}

}  // namespace astradock::orbit
