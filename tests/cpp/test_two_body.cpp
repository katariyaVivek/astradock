#include "dynamics/two_body.hpp"
#include "math/constants.hpp"
#include "math/vector3.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace {

using astradock::constants::earth_gravitational_parameter_m3_per_s2;
using astradock::constants::earth_reference_radius_m;
using astradock::dynamics::two_body_acceleration;
using astradock::math::Vector3;
using Catch::Approx;

TEST_CASE("Two-body gravity has the analytical radial magnitude and inward direction",
          "[dynamics][two-body]") {
    constexpr double radius_m = 7.0e6;
    constexpr double mu_m3_per_s2 = earth_gravitational_parameter_m3_per_s2;
    const Vector3 acceleration =
        two_body_acceleration(Vector3{radius_m, 0.0, 0.0}, mu_m3_per_s2);
    const double expected_magnitude_m_per_s2 = mu_m3_per_s2 / (radius_m * radius_m);

    REQUIRE(acceleration.x() < 0.0);
    REQUIRE(acceleration.y() == Approx(0.0));
    REQUIRE(acceleration.z() == Approx(0.0));
    REQUIRE(acceleration.norm() == Approx(expected_magnitude_m_per_s2).epsilon(1.0e-12));
}

TEST_CASE("Two-body gravity points toward the origin on every coordinate axis",
          "[dynamics][two-body]") {
    constexpr double radius_m = 7.0e6;
    constexpr double mu_m3_per_s2 = earth_gravitational_parameter_m3_per_s2;

    const Vector3 acceleration_x =
        two_body_acceleration(Vector3{radius_m, 0.0, 0.0}, mu_m3_per_s2);
    const Vector3 acceleration_y =
        two_body_acceleration(Vector3{0.0, radius_m, 0.0}, mu_m3_per_s2);
    const Vector3 acceleration_z =
        two_body_acceleration(Vector3{0.0, 0.0, radius_m}, mu_m3_per_s2);
    const double expected_magnitude_m_per_s2 = mu_m3_per_s2 / (radius_m * radius_m);

    REQUIRE(astradock::math::approximately_equal(
        acceleration_x, Vector3{-expected_magnitude_m_per_s2, 0.0, 0.0}));
    REQUIRE(astradock::math::approximately_equal(
        acceleration_y, Vector3{0.0, -expected_magnitude_m_per_s2, 0.0}));
    REQUIRE(astradock::math::approximately_equal(
        acceleration_z, Vector3{0.0, 0.0, -expected_magnitude_m_per_s2}));
    REQUIRE(acceleration_x.norm() == Approx(acceleration_y.norm()));
    REQUIRE(acceleration_y.norm() == Approx(acceleration_z.norm()));
}

TEST_CASE("Two-body gravity is anti-parallel to an arbitrary position",
          "[dynamics][two-body]") {
    const Vector3 position_central_body_inertial_m{7.0e6, -2.0e6, 1.0e6};
    const Vector3 acceleration_central_body_inertial_m_per_s2 = two_body_acceleration(
        position_central_body_inertial_m, earth_gravitational_parameter_m3_per_s2);
    const Vector3 radial_cross_acceleration =
        position_central_body_inertial_m.cross(acceleration_central_body_inertial_m_per_s2);
    const double normalized_cross_magnitude =
        radial_cross_acceleration.norm()
        / (position_central_body_inertial_m.norm()
           * acceleration_central_body_inertial_m_per_s2.norm());

    REQUIRE(position_central_body_inertial_m.dot(acceleration_central_body_inertial_m_per_s2)
            < 0.0);
    REQUIRE(normalized_cross_magnitude == Approx(0.0).margin(1.0e-14));
    REQUIRE(astradock::math::approximately_equal(
        acceleration_central_body_inertial_m_per_s2.normalized(),
        -1.0 * position_central_body_inertial_m.normalized(),
        1.0e-13,
        1.0e-13));
}

TEST_CASE("Two-body gravity follows the inverse-square relationship",
          "[dynamics][two-body]") {
    constexpr double radius_m = 7.0e6;
    const double acceleration_at_radius_m_per_s2 =
        two_body_acceleration(
            Vector3{radius_m, 0.0, 0.0}, earth_gravitational_parameter_m3_per_s2)
            .norm();
    const double acceleration_at_twice_radius_m_per_s2 =
        two_body_acceleration(
            Vector3{2.0 * radius_m, 0.0, 0.0}, earth_gravitational_parameter_m3_per_s2)
            .norm();

    REQUIRE(acceleration_at_twice_radius_m_per_s2
            == Approx(acceleration_at_radius_m_per_s2 / 4.0).epsilon(1.0e-12));
}

TEST_CASE("Two-body Earth surface gravity is approximately 9.8 metres per square second",
          "[dynamics][two-body][earth]") {
    const double surface_gravity_m_per_s2 =
        two_body_acceleration(
            Vector3{earth_reference_radius_m, 0.0, 0.0},
            earth_gravitational_parameter_m3_per_s2)
            .norm();

    // This broad physical sanity check deliberately does not model Earth's
    // rotation, ellipsoidal shape, altitude variation, or higher gravity terms.
    REQUIRE(surface_gravity_m_per_s2 == Approx(9.8).margin(0.02));
}

TEST_CASE("Analytical circular speed at 500 kilometres is physically reasonable",
          "[dynamics][two-body][earth]") {
    constexpr double altitude_m = 500'000.0;
    const double orbital_radius_m = earth_reference_radius_m + altitude_m;
    const double circular_speed_m_per_s =
        std::sqrt(earth_gravitational_parameter_m3_per_s2 / orbital_radius_m);

    REQUIRE(circular_speed_m_per_s == Approx(7'612.6).margin(0.1));
    REQUIRE(circular_speed_m_per_s > 7'500.0);
    REQUIRE(circular_speed_m_per_s < 7'800.0);
}

TEST_CASE("Two-body gravity rejects the central-body origin", "[dynamics][two-body]") {
    REQUIRE_THROWS_AS(
        two_body_acceleration(Vector3{}, earth_gravitational_parameter_m3_per_s2),
        std::domain_error);
}

TEST_CASE("Two-body gravity rejects invalid physical inputs", "[dynamics][two-body]") {
    const Vector3 valid_position_m{earth_reference_radius_m, 0.0, 0.0};
    const double infinity = std::numeric_limits<double>::infinity();
    const double not_a_number = std::numeric_limits<double>::quiet_NaN();

    REQUIRE_THROWS_AS((two_body_acceleration(valid_position_m, 0.0)), std::domain_error);
    REQUIRE_THROWS_AS((two_body_acceleration(valid_position_m, -1.0)), std::domain_error);
    REQUIRE_THROWS_AS((two_body_acceleration(valid_position_m, infinity)), std::domain_error);
    REQUIRE_THROWS_AS((two_body_acceleration(valid_position_m, not_a_number)), std::domain_error);
    REQUIRE_THROWS_AS(
        (two_body_acceleration(Vector3{infinity, 0.0, 0.0},
                               earth_gravitational_parameter_m3_per_s2)),
        std::domain_error);
}

}  // namespace
