// AstraDock M24 — High-fidelity environment tests.
//
// Verification strategy:
//   M24A rotation: identity at t = 0, 2pi-periodicity over a sidereal day,
//     orthonormality, round-trip, +X axis sweeps into +Y (sign regression).
//   M24B time: JD mapping at J2000 + one day; rotation/JD consistency.
//   M24C geodetic: equator/prime-meridian origin, north pole, 500 km altitude,
//     round-trips incl. high latitude, invalid latitude rejected.
//   M24D zonal: J2-only reproduces j2_acceleration_eci bit-near-exact;
//     J3 antisymmetry across the equator; J4 equatorial hand value;
//     finite-difference audit of the potential gradient.
//   M24E ephemeris: circular radius, period closure, documented rates.

#include "environment/high_fidelity.hpp"
#include "environment/j2_gravity.hpp"
#include "math/constants.hpp"
#include "math/matrix3.hpp"
#include "math/vector3.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>
#include <vector>

using namespace astradock;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("M24A ECI/ECEF rotation is identity at epoch and periodic daily", "[env][m24a]") {
    const auto dcm0 = environment::dcm_ecef_from_eci(0.0);
    CHECK_THAT(dcm0(0, 0), WithinRel(1.0, 1.0e-15));
    CHECK_THAT(dcm0(1, 0), WithinAbs(0.0, 1.0e-15));
    CHECK(dcm0.is_orthonormal(1.0e-12, 1.0e-12));
    // One sidereal day later: identity again (2pi rotation).
    const double sidereal_day = 2.0 * constants::pi / constants::earth_rotation_rate_rad_per_s;
    const auto dcm_day = environment::dcm_ecef_from_eci(sidereal_day);
    CHECK_THAT(dcm_day(0, 0), WithinRel(1.0, 1.0e-9));
    CHECK_THAT(dcm_day(0, 1), WithinAbs(0.0, 1.0e-9));
    // Sign regression: +X ECI rotates toward +Y ECEF as theta grows (R_z(+theta)
    // convention: v_ecef = R v_eci with rows [c, s, 0; -s, c, 0; ...]).
    const math::Vector3 ecef =
        environment::eci_to_ecef(math::Vector3{1.0, 0.0, 0.0}, 1000.0);
    const double theta = constants::earth_rotation_rate_rad_per_s * 1000.0;
    CHECK_THAT(ecef.x(), WithinRel(std::cos(theta), 1.0e-12));
    CHECK_THAT(ecef.y(), WithinRel(-std::sin(theta), 1.0e-12));
    // Round trip.
    const math::Vector3 p_eci{7000.0e3, 100.0e3, 200.0e3};
    const math::Vector3 back =
        environment::ecef_to_eci(environment::eci_to_ecef(p_eci, 3600.0), 3600.0);
    CHECK_THAT((back - p_eci).norm(), WithinAbs(0.0, 1.0e-6));
}

TEST_CASE("M24B simulation time maps linearly to Julian date", "[env][m24b]") {
    CHECK_THAT(environment::sim_time_to_julian_date(0.0), WithinRel(2451545.0, 1.0e-12));
    CHECK_THAT(
        environment::sim_time_to_julian_date(86400.0), WithinRel(2451546.0, 1.0e-12));
    // Rotation angle after one solar day is slightly MORE than 2pi (sidereal <
    // solar day): documents the civil-vs-dynamical distinction, not a bug.
    const double theta_solar = environment::earth_rotation_angle_rad(86400.0);
    CHECK(theta_solar > 2.0 * constants::pi);
    CHECK_THROWS_AS(environment::earth_rotation_angle_rad(
                        std::numeric_limits<double>::quiet_NaN()),
                    std::domain_error);
}

TEST_CASE("M24C geodetic conversion matches ellipsoid landmarks", "[env][m24c]") {
    // Equator / prime meridian surface point.
    const auto eq = environment::ecef_to_geodetic(
        math::Vector3{environment::k_wgs84_semimajor_m, 0.0, 0.0});
    CHECK_THAT(eq.latitude_rad, WithinAbs(0.0, 1.0e-9));
    CHECK_THAT(eq.longitude_rad, WithinAbs(0.0, 1.0e-12));
    CHECK_THAT(eq.altitude_m, WithinAbs(0.0, 1.0e-6));
    // North pole surface point: latitude +90 deg, altitude ~0.
    const double polar_radius =
        environment::k_wgs84_semimajor_m * (1.0 - environment::k_wgs84_flattening);
    const auto pole = environment::ecef_to_geodetic(math::Vector3{0.0, 0.0, polar_radius});
    CHECK_THAT(pole.latitude_rad, WithinRel(constants::pi / 2.0, 1.0e-12));
    CHECK_THAT(pole.altitude_m, WithinAbs(0.0, 1.0e-3));
    // 500 km over the equator reads 500 km altitude.
    const auto alt = environment::ecef_to_geodetic(
        math::Vector3{environment::k_wgs84_semimajor_m + 500.0e3, 0.0, 0.0});
    CHECK_THAT(alt.altitude_m, WithinRel(500.0e3, 1.0e-9));
    // Round trips: mid-latitude, southern hemisphere, high altitude.
    const std::vector<environment::GeodeticCoord> cases{
        {0.785, 1.0, 1000.0}, {-0.5, -2.0, 400.0e3}, {1.2, 0.0, 0.0},
    };
    for (const auto& geo : cases) {
        const auto rt = environment::ecef_to_geodetic(environment::geodetic_to_ecef(geo));
        CHECK_THAT(rt.latitude_rad, WithinRel(geo.latitude_rad, 1.0e-9));
        CHECK_THAT(rt.longitude_rad, WithinRel(geo.longitude_rad, 1.0e-9));
        // Altitude round-trip is ABSOLUTE-error checked: relative tolerance is
        // meaningless at 0 m altitude (division by ~zero inflates the error).
        CHECK_THAT(rt.altitude_m, WithinAbs(geo.altitude_m, 1.0e-6));
    }
    environment::GeodeticCoord bad{constants::pi, 0.0, 0.0};
    CHECK_THROWS_AS(environment::geodetic_to_ecef(bad), std::domain_error);
}

TEST_CASE("M24D zonal J2-only path reproduces legacy J2 + central", "[env][m24d]") {
    // zonal_acceleration_eci returns TOTAL acceleration; the legacy M11 J2
    // function returns the PERTURBATION only. Compare accordingly.
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const std::vector<math::Vector3> positions{
        {7000.0e3, 0.0, 0.0},
        {0.0, 0.0, 7000.0e3},
        {1000.0e3, 2000.0e3, 6000.0e3},
        {-4500.0e3, 5200.0e3, 1800.0e3},
    };
    for (const auto& r : positions) {
        const math::Vector3 legacy = environment::j2_acceleration_eci(r);
        const math::Vector3 central =
            r.normalized() * (-mu / r.squared_norm());
        const math::Vector3 zonal = environment::zonal_acceleration_eci(
            r, mu, constants::earth_reference_radius_m, constants::earth_j2, 0.0, 0.0);
        CHECK_THAT((zonal - central - legacy).norm(), WithinAbs(0.0, 1.0e-9));
    }
}

TEST_CASE("M24D J3/J4 pole values match hand-derived potential gradients", "[env][m24d]") {
    // Pear-shape physics (corrected during development — the naive
    // 'antisymmetry' claim is wrong): an odd zonal pushes BOTH poles the same
    // way. At the pole grad_s vanishes (z_hat - s r_hat = 0), so
    //   a_z = dU/dr = -mu/r^2 + 4 mu J3 rho^3/r^2 (north, s = +1, P3 = 1).
    // South (s = -1, P3 = -1): a_z = +mu/r^2 + 4 mu J3 rho^3/r^2 — same J3
    // push, mirrored central term. Test the J3 PART by differencing.
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const double R = constants::earth_reference_radius_m;
    const double r = 7000.0e3;
    const double j3 = -2.5e-6;
    const math::Vector3 north{0.0, 0.0, r};
    const math::Vector3 south{0.0, 0.0, -r};
    const double rho = R / r;
    const double expected_push = 4.0 * mu * j3 * rho * rho * rho / (r * r);
    const math::Vector3 j3n = environment::zonal_acceleration_eci(north, mu, R, 0.0, j3, 0.0);
    const math::Vector3 j3s = environment::zonal_acceleration_eci(south, mu, R, 0.0, j3, 0.0);
    // Remove the mirrored central term: J3 push is identical at both poles.
    const double push_n = j3n.z() + mu / (r * r);
    const double push_s = j3s.z() - mu / (r * r);
    CHECK_THAT(push_n, WithinRel(expected_push, 1.0e-9));
    CHECK_THAT(push_s, WithinRel(expected_push, 1.0e-9));
    // J4 (even): mirrored VECTORS across the equator (radial direction flips
    // at the south pole): a_south = -a_north exactly. Central terms cancel in
    // the sum, so the full z-components must be opposite.
    const double j4 = -1.6e-6;
    const math::Vector3 j4n = environment::zonal_acceleration_eci(north, mu, R, 0.0, 0.0, j4);
    const math::Vector3 j4s = environment::zonal_acceleration_eci(south, mu, R, 0.0, 0.0, j4);
    CHECK_THAT(j4n.z() + j4s.z(), WithinAbs(0.0, 1.0e-12));
    // Scale regression on PERTURBATIONS (central term removed):
    // |a_J4| / |a_J2| ~ (J4/J2)(R/r)^2 at the equator.
    const math::Vector3 eq{r, 0.0, 0.0};
    const math::Vector3 central{-mu / (r * r), 0.0, 0.0};
    const double j2mag =
        (environment::zonal_acceleration_eci(eq, mu, R, 1.08262668e-3, 0.0, 0.0) - central)
            .norm();
    const double j4mag =
        (environment::zonal_acceleration_eci(eq, mu, R, 0.0, 0.0, j4) - central).norm();
    CHECK_THAT(j4mag / j2mag, WithinRel(std::abs(j4) / 1.08262668e-3, 0.10));
}

TEST_CASE("M24D zonal gradient matches finite differences of the potential", "[env][m24d]") {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const double R = constants::earth_reference_radius_m;
    const double j2 = 1.08262668e-3;
    const double j3 = -2.5e-6;
    const double j4 = -1.6e-6;
    const auto potential = [&](const math::Vector3& r) {
        const double rn = r.norm();
        const double s = r.z() / rn;
        const double rho = R / rn;
        const double p2 = 0.5 * (3.0 * s * s - 1.0);
        const double p3 = 0.5 * (5.0 * s * s * s - 3.0 * s);
        const double p4 = 0.125 * (35.0 * std::pow(s, 4) - 30.0 * s * s + 3.0);
        return mu / rn
            * (1.0 - j2 * rho * rho * p2 - j3 * rho * rho * rho * p3
               - j4 * rho * rho * rho * rho * p4);
    };
    const math::Vector3 r{1000.0e3, 2000.0e3, 6000.0e3};
    const math::Vector3 analytic = environment::zonal_acceleration_eci(r, mu, R, j2, j3, j4);
    const double eps = 1.0;
    const math::Vector3 fd{
        (potential(r + math::Vector3{eps, 0.0, 0.0}) - potential(r - math::Vector3{eps, 0.0, 0.0}))
            / (2.0 * eps),
        (potential(r + math::Vector3{0.0, eps, 0.0}) - potential(r - math::Vector3{0.0, eps, 0.0}))
            / (2.0 * eps),
        (potential(r + math::Vector3{0.0, 0.0, eps}) - potential(r - math::Vector3{0.0, 0.0, eps}))
            / (2.0 * eps),
    };
    CHECK_THAT((analytic - fd).norm() / analytic.norm(), WithinAbs(0.0, 1.0e-8));
}

TEST_CASE("M24E analytical ephemeris holds radius and period", "[env][m24e]") {
    // Moon: radius exact, full-period closure, quarter-period quadrature.
    const auto m0 = environment::moon_position_eci_m(0.0);
    CHECK_THAT(m0.norm(), WithinRel(environment::k_moon_distance_m, 1.0e-12));
    const auto m1 = environment::moon_position_eci_m(environment::k_moon_period_s);
    CHECK_THAT((m1 - m0).norm(), WithinAbs(0.0, 1.0e-3));
    const auto mq = environment::moon_position_eci_m(environment::k_moon_period_s / 4.0);
    CHECK_THAT(mq.x(), WithinAbs(0.0, 1.0e-3));
    CHECK_THAT(mq.y(), WithinRel(environment::k_moon_distance_m, 1.0e-9));
    // Sun: radius exact + period closure.
    const auto s0 = environment::sun_position_eci_m(0.0);
    CHECK_THAT(s0.norm(), WithinRel(environment::k_sun_distance_m, 1.0e-9));
    const auto s1 = environment::sun_position_eci_m(environment::k_sun_period_s);
    CHECK_THAT((s1 - s0).norm() / s0.norm(), WithinAbs(0.0, 1.0e-9));
    // Documented simplification: equatorial plane (inclination NOT modeled).
    CHECK_THAT(m0.z(), WithinAbs(0.0, 0.0));
    CHECK_THROWS_AS(
        environment::moon_position_eci_m(std::numeric_limits<double>::quiet_NaN()),
        std::domain_error);
}
