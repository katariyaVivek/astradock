#include "math/constants.hpp"
#include "math/vector3.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <numbers>
#include <stdexcept>

namespace {

using astradock::math::Vector3;
using Catch::Approx;

TEST_CASE("Vector3 addition is component-wise", "[math][vector3]") {
    const Vector3 result = Vector3{1.0, -2.0, 3.5} + Vector3{4.0, 5.0, -0.5};

    REQUIRE(result.x() == Approx(5.0));
    REQUIRE(result.y() == Approx(3.0));
    REQUIRE(result.z() == Approx(3.0));
}

TEST_CASE("Vector3 subtraction is component-wise", "[math][vector3]") {
    const Vector3 result = Vector3{4.0, 5.0, 6.0} - Vector3{1.0, 2.0, 3.0};

    REQUIRE(result.x() == Approx(3.0));
    REQUIRE(result.y() == Approx(3.0));
    REQUIRE(result.z() == Approx(3.0));
}

TEST_CASE("Vector3 supports scalar multiplication on both sides", "[math][vector3]") {
    const Vector3 vector{1.0, -2.0, 0.5};

    REQUIRE(astradock::math::approximately_equal(vector * 3.0, Vector3{3.0, -6.0, 1.5}));
    REQUIRE(astradock::math::approximately_equal(3.0 * vector, Vector3{3.0, -6.0, 1.5}));
}

TEST_CASE("Vector3 supports nonzero scalar division", "[math][vector3]") {
    const Vector3 result = Vector3{6.0, -3.0, 1.5} / 3.0;

    REQUIRE(astradock::math::approximately_equal(result, Vector3{2.0, -1.0, 0.5}));
    REQUIRE_THROWS_AS((Vector3{1.0, 2.0, 3.0} / 0.0), std::domain_error);
}

TEST_CASE("Vector3 dot product matches an analytical result", "[math][vector3]") {
    const Vector3 lhs{1.0, 2.0, 3.0};
    const Vector3 rhs{4.0, -5.0, 6.0};

    REQUIRE(lhs.dot(rhs) == Approx(12.0));
}

TEST_CASE("Vector3 cross product follows the right-hand rule", "[math][vector3]") {
    const Vector3 x_axis{1.0, 0.0, 0.0};
    const Vector3 y_axis{0.0, 1.0, 0.0};

    REQUIRE(astradock::math::approximately_equal(x_axis.cross(y_axis), Vector3{0.0, 0.0, 1.0}));
    REQUIRE(astradock::math::approximately_equal(y_axis.cross(x_axis), Vector3{0.0, 0.0, -1.0}));
}

TEST_CASE("Vector3 norm and squared norm match a 3-4-5 triangle", "[math][vector3]") {
    const Vector3 vector{3.0, 4.0, 0.0};

    REQUIRE(vector.squared_norm() == Approx(25.0));
    REQUIRE(vector.norm() == Approx(5.0));
}

TEST_CASE("Vector3 normalization produces a unit vector", "[math][vector3]") {
    const Vector3 normalized = Vector3{3.0, 4.0, 0.0}.normalized();

    REQUIRE(normalized.x() == Approx(0.6));
    REQUIRE(normalized.y() == Approx(0.8));
    REQUIRE(normalized.z() == Approx(0.0));
    REQUIRE(normalized.norm() == Approx(1.0));
}

TEST_CASE("Vector3 cross product is orthogonal to both operands", "[math][vector3]") {
    const Vector3 lhs{2.0, -3.0, 4.0};
    const Vector3 rhs{-1.0, 5.0, 2.0};
    const Vector3 cross = lhs.cross(rhs);

    REQUIRE(cross.dot(lhs) == Approx(0.0).margin(1.0e-12));
    REQUIRE(cross.dot(rhs) == Approx(0.0).margin(1.0e-12));
}

TEST_CASE("Normalizing the zero vector fails explicitly", "[math][vector3]") {
    REQUIRE_THROWS_AS(Vector3{}.normalized(), std::domain_error);
}

TEST_CASE("The minimal constants header exposes pi", "[math][constants]") {
    REQUIRE(astradock::constants::pi == Approx(std::numbers::pi_v<double>));
}

}  // namespace
