#include "frames/frame_basis.hpp"
#include "frames/lvlh.hpp"
#include "math/constants.hpp"
#include "math/matrix3.hpp"
#include "math/vector3.hpp"
#include "numerics/fixed_step_propagation.hpp"
#include "orbit/cartesian_state.hpp"
#include "orbit/two_body_orbit.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

using astradock::frames::compute_lvlh_basis;
using astradock::frames::dcm_eci_from_lvlh;
using astradock::frames::dcm_lvlh_from_eci;
using astradock::frames::FrameBasis;
using astradock::frames::is_orthonormal;
using astradock::frames::transform_eci_to_lvlh;
using astradock::frames::transform_lvlh_to_eci;
using astradock::math::approximately_equal;
using astradock::math::Matrix3;
using astradock::math::Vector3;
using astradock::numerics::IntegrationMethod;
using astradock::numerics::propagate_fixed_step;
using astradock::orbit::CartesianState;
using astradock::orbit::compute_circular_orbit_reference;
using astradock::orbit::two_body_state_derivative;
using Catch::Matchers::WithinAbs;

TEST_CASE("FrameBasis validation detects orthonormal and non-orthonormal triads", "[frames][basis]") {
    // Valid standard right-handed Cartesian basis
    const FrameBasis standard_basis{
        Vector3(1.0, 0.0, 0.0),
        Vector3(0.0, 1.0, 0.0),
        Vector3(0.0, 0.0, 1.0)
    };
    CHECK(is_orthonormal(standard_basis));

    // Non-unit norm basis vector
    const FrameBasis non_unit_basis{
        Vector3(1.01, 0.0, 0.0),
        Vector3(0.0, 1.0, 0.0),
        Vector3(0.0, 0.0, 1.0)
    };
    CHECK_FALSE(is_orthonormal(non_unit_basis));

    // Non-orthogonal vectors
    const FrameBasis non_orthogonal_basis{
        Vector3(1.0, 0.0, 0.0),
        Vector3(0.1, 0.99498743710662, 0.0),
        Vector3(0.0, 0.0, 1.0)
    };
    CHECK_FALSE(is_orthonormal(non_orthogonal_basis));

    // Left-handed basis (reflection)
    const FrameBasis left_handed_basis{
        Vector3(1.0, 0.0, 0.0),
        Vector3(0.0, 1.0, 0.0),
        Vector3(0.0, 0.0, -1.0)
    };
    CHECK_FALSE(is_orthonormal(left_handed_basis));

    // Non-finite basis
    const FrameBasis nan_basis{
        Vector3(std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0),
        Vector3(0.0, 1.0, 0.0),
        Vector3(0.0, 0.0, 1.0)
    };
    CHECK_FALSE(is_orthonormal(nan_basis));
}

TEST_CASE("Simple circular orbit produces identity-like LVLH basis", "[frames][lvlh]") {
    const double radius = 6878137.0;
    const double speed = 7612.608173;

    // Initial state: r along +X, v along +Y
    const Vector3 r_eci(radius, 0.0, 0.0);
    const Vector3 v_eci(0.0, speed, 0.0);

    const FrameBasis basis = compute_lvlh_basis(r_eci, v_eci);
    CHECK(is_orthonormal(basis));

    // Expected axes:
    // e_r = +X
    CHECK_THAT(basis.x.x(), WithinAbs(1.0, 1.0e-14));
    CHECK_THAT(basis.x.y(), WithinAbs(0.0, 1.0e-14));
    CHECK_THAT(basis.x.z(), WithinAbs(0.0, 1.0e-14));

    // e_t = +Y
    CHECK_THAT(basis.y.x(), WithinAbs(0.0, 1.0e-14));
    CHECK_THAT(basis.y.y(), WithinAbs(1.0, 1.0e-14));
    CHECK_THAT(basis.y.z(), WithinAbs(0.0, 1.0e-14));

    // e_h = +Z
    CHECK_THAT(basis.z.x(), WithinAbs(0.0, 1.0e-14));
    CHECK_THAT(basis.z.y(), WithinAbs(0.0, 1.0e-14));
    CHECK_THAT(basis.z.z(), WithinAbs(1.0, 1.0e-14));

    // DCM should equal identity
    const Matrix3 c_lvlh_eci = dcm_lvlh_from_eci(basis);
    CHECK(approximately_equal(c_lvlh_eci, Matrix3::identity()));
}

TEST_CASE("Circular orbit analytical quadrant check (quarter, half, three-quarter)", "[frames][lvlh]") {
    const double r = 6878137.0;
    const double v = 7612.608173;

    // 1. Quarter orbit (theta = 90 deg): r along +Y, v along -X
    {
        const Vector3 r_90(0.0, r, 0.0);
        const Vector3 v_90(-v, 0.0, 0.0);
        const FrameBasis basis = compute_lvlh_basis(r_90, v_90);

        CHECK(is_orthonormal(basis));
        // e_r points along +Y
        CHECK_THAT(basis.x.x(), WithinAbs(0.0, 1.0e-12));
        CHECK_THAT(basis.x.y(), WithinAbs(1.0, 1.0e-12));
        CHECK_THAT(basis.x.z(), WithinAbs(0.0, 1.0e-12));

        // e_t points along -X
        CHECK_THAT(basis.y.x(), WithinAbs(-1.0, 1.0e-12));
        CHECK_THAT(basis.y.y(), WithinAbs(0.0, 1.0e-12));
        CHECK_THAT(basis.y.z(), WithinAbs(0.0, 1.0e-12));

        // e_h points along +Z
        CHECK_THAT(basis.z.x(), WithinAbs(0.0, 1.0e-12));
        CHECK_THAT(basis.z.y(), WithinAbs(0.0, 1.0e-12));
        CHECK_THAT(basis.z.z(), WithinAbs(1.0, 1.0e-12));
    }

    // 2. Half orbit (theta = 180 deg): r along -X, v along -Y
    {
        const Vector3 r_180(-r, 0.0, 0.0);
        const Vector3 v_180(0.0, -v, 0.0);
        const FrameBasis basis = compute_lvlh_basis(r_180, v_180);

        CHECK(is_orthonormal(basis));
        // e_r points along -X
        CHECK_THAT(basis.x.x(), WithinAbs(-1.0, 1.0e-12));
        CHECK_THAT(basis.x.y(), WithinAbs(0.0, 1.0e-12));
        CHECK_THAT(basis.x.z(), WithinAbs(0.0, 1.0e-12));

        // e_t points along -Y
        CHECK_THAT(basis.y.x(), WithinAbs(0.0, 1.0e-12));
        CHECK_THAT(basis.y.y(), WithinAbs(-1.0, 1.0e-12));
        CHECK_THAT(basis.y.z(), WithinAbs(0.0, 1.0e-12));

        // e_h points along +Z
        CHECK_THAT(basis.z.x(), WithinAbs(0.0, 1.0e-12));
        CHECK_THAT(basis.z.y(), WithinAbs(0.0, 1.0e-12));
        CHECK_THAT(basis.z.z(), WithinAbs(1.0, 1.0e-12));
    }

    // 3. Three-quarter orbit (theta = 270 deg): r along -Y, v along +X
    {
        const Vector3 r_270(0.0, -r, 0.0);
        const Vector3 v_270(v, 0.0, 0.0);
        const FrameBasis basis = compute_lvlh_basis(r_270, v_270);

        CHECK(is_orthonormal(basis));
        // e_r points along -Y
        CHECK_THAT(basis.x.x(), WithinAbs(0.0, 1.0e-12));
        CHECK_THAT(basis.x.y(), WithinAbs(-1.0, 1.0e-12));
        CHECK_THAT(basis.x.z(), WithinAbs(0.0, 1.0e-12));

        // e_t points along +X
        CHECK_THAT(basis.y.x(), WithinAbs(1.0, 1.0e-12));
        CHECK_THAT(basis.y.y(), WithinAbs(0.0, 1.0e-12));
        CHECK_THAT(basis.y.z(), WithinAbs(0.0, 1.0e-12));

        // e_h points along +Z
        CHECK_THAT(basis.z.x(), WithinAbs(0.0, 1.0e-12));
        CHECK_THAT(basis.z.y(), WithinAbs(0.0, 1.0e-12));
        CHECK_THAT(basis.z.z(), WithinAbs(1.0, 1.0e-12));
    }
}

TEST_CASE("Arbitrary non-axis-aligned orbit constructs valid orthonormal right-handed basis", "[frames][lvlh]") {
    const Vector3 r_eci(3500000.0, -4200000.0, 5100000.0);
    const Vector3 v_eci(-4000.0, 3500.0, 5500.0);

    const FrameBasis basis = compute_lvlh_basis(r_eci, v_eci);
    CHECK(is_orthonormal(basis));

    // Orthonormality checks
    CHECK_THAT(basis.x.norm(), WithinAbs(1.0, 1.0e-12));
    CHECK_THAT(basis.y.norm(), WithinAbs(1.0, 1.0e-12));
    CHECK_THAT(basis.z.norm(), WithinAbs(1.0, 1.0e-12));

    CHECK_THAT(basis.x.dot(basis.y), WithinAbs(0.0, 1.0e-12));
    CHECK_THAT(basis.y.dot(basis.z), WithinAbs(0.0, 1.0e-12));
    CHECK_THAT(basis.z.dot(basis.x), WithinAbs(0.0, 1.0e-12));

    // Right-handed orientation: e_r x e_t = e_h
    const Vector3 cross_rt = basis.x.cross(basis.y);
    CHECK(approximately_equal(cross_rt, basis.z, 1.0e-12, 1.0e-12));
    CHECK_THAT(cross_rt.dot(basis.z), WithinAbs(1.0, 1.0e-12));

    // Check DCM properties
    const Matrix3 c_lvlh_eci = dcm_lvlh_from_eci(basis);
    const Matrix3 c_eci_lvlh = dcm_eci_from_lvlh(basis);

    CHECK(c_lvlh_eci.is_orthonormal());
    CHECK(c_eci_lvlh.is_orthonormal());
    CHECK_THAT(c_lvlh_eci.determinant(), WithinAbs(1.0, 1.0e-12));
    CHECK_THAT(c_eci_lvlh.determinant(), WithinAbs(1.0, 1.0e-12));

    // Transpose relationship
    CHECK(approximately_equal(c_lvlh_eci.transpose(), c_eci_lvlh));
    CHECK(approximately_equal(c_lvlh_eci * c_eci_lvlh, Matrix3::identity()));
}

TEST_CASE("DCM transforms position and velocity vectors correctly into LVLH", "[frames][transform]") {
    const double radius = 6878137.0;
    const double speed = 7612.608173;

    const Vector3 r_eci(radius, 0.0, 0.0);
    const Vector3 v_eci(0.0, speed, 0.0);

    const Matrix3 c_lvlh_eci = dcm_lvlh_from_eci(r_eci, v_eci);

    // Transforming position into LVLH:
    // In LVLH, spacecraft position is entirely along radial axis: [ ||r||, 0, 0 ]
    const Vector3 r_lvlh = transform_eci_to_lvlh(r_eci, c_lvlh_eci);
    CHECK_THAT(r_lvlh.x(), WithinAbs(radius, 1.0e-8));
    CHECK_THAT(r_lvlh.y(), WithinAbs(0.0, 1.0e-8));
    CHECK_THAT(r_lvlh.z(), WithinAbs(0.0, 1.0e-8));

    // Transforming circular velocity into LVLH:
    // For a circular orbit, velocity is purely along-track: [ 0, ||v||, 0 ]
    const Vector3 v_lvlh = transform_eci_to_lvlh(v_eci, c_lvlh_eci);
    CHECK_THAT(v_lvlh.x(), WithinAbs(0.0, 1.0e-8));
    CHECK_THAT(v_lvlh.y(), WithinAbs(speed, 1.0e-8));
    CHECK_THAT(v_lvlh.z(), WithinAbs(0.0, 1.0e-8));

    // Norm preservation: ||v_LVLH|| == ||v_ECI||
    CHECK_THAT(r_lvlh.norm(), WithinAbs(r_eci.norm(), 1.0e-8));
    CHECK_THAT(v_lvlh.norm(), WithinAbs(v_eci.norm(), 1.0e-8));

    // Round-trip transformation
    const Vector3 r_eci_round_trip = transform_lvlh_to_eci(r_lvlh, c_lvlh_eci);
    const Vector3 v_eci_round_trip = transform_lvlh_to_eci(v_lvlh, c_lvlh_eci);
    CHECK(approximately_equal(r_eci_round_trip, r_eci, 1.0e-8, 1.0e-8));
    CHECK(approximately_equal(v_eci_round_trip, v_eci, 1.0e-8, 1.0e-8));
}

TEST_CASE("Arbitrary vector transformation preserves norm and satisfies round-trip", "[frames][transform]") {
    const Vector3 r_eci(3500000.0, -4200000.0, 5100000.0);
    const Vector3 v_eci(-4000.0, 3500.0, 5500.0);
    const Matrix3 c_lvlh_eci = dcm_lvlh_from_eci(r_eci, v_eci);

    const Vector3 test_vec_eci(1234.5, -6789.1, 4321.0);
    const Vector3 test_vec_lvlh = transform_eci_to_lvlh(test_vec_eci, c_lvlh_eci);

    // Magnitude invariance
    CHECK_THAT(test_vec_lvlh.norm(), WithinAbs(test_vec_eci.norm(), 1.0e-10));

    // Round trip
    const Vector3 test_vec_reconstructed = transform_lvlh_to_eci(test_vec_lvlh, c_lvlh_eci);
    CHECK(approximately_equal(test_vec_reconstructed, test_vec_eci, 1.0e-10, 1.0e-10));
}

TEST_CASE("Degenerate and invalid orbital states are rejected", "[frames][degeneracy]") {
    const Vector3 valid_r(6878137.0, 0.0, 0.0);
    const Vector3 valid_v(0.0, 7612.608173, 0.0);

    // 1. Zero position magnitude
    CHECK_THROWS_AS(
        compute_lvlh_basis(Vector3(0.0, 0.0, 0.0), valid_v),
        std::domain_error);

    // 2. Zero velocity (angular momentum = 0)
    CHECK_THROWS_AS(
        compute_lvlh_basis(valid_r, Vector3(0.0, 0.0, 0.0)),
        std::domain_error);

    // 3. Collinear position and velocity (radial trajectory, angular momentum = 0)
    CHECK_THROWS_AS(
        compute_lvlh_basis(valid_r, Vector3(1000.0, 0.0, 0.0)),
        std::domain_error);
    CHECK_THROWS_AS(
        compute_lvlh_basis(valid_r, Vector3(-500.0, 0.0, 0.0)),
        std::domain_error);

    // 4. Non-finite values in position
    CHECK_THROWS_AS(
        compute_lvlh_basis(
            Vector3(std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0),
            valid_v),
        std::domain_error);
    CHECK_THROWS_AS(
        compute_lvlh_basis(
            Vector3(std::numeric_limits<double>::infinity(), 0.0, 0.0),
            valid_v),
        std::domain_error);

    // 5. Non-finite values in velocity
    CHECK_THROWS_AS(
        compute_lvlh_basis(
            valid_r,
            Vector3(0.0, std::numeric_limits<double>::quiet_NaN(), 0.0)),
        std::domain_error);
    CHECK_THROWS_AS(
        compute_lvlh_basis(
            valid_r,
            Vector3(0.0, std::numeric_limits<double>::infinity(), 0.0)),
        std::domain_error);
}

TEST_CASE("Propagated circular orbit maintains orthonormal LVLH basis across full period", "[frames][propagation]") {
    const double altitude_m = 500000.0;
    const double radius_m = astradock::constants::earth_reference_radius_m + altitude_m;
    const auto ref = compute_circular_orbit_reference(
        astradock::constants::earth_gravitational_parameter_m3_per_s2,
        radius_m);

    const CartesianState initial_state{
        Vector3(radius_m, 0.0, 0.0),
        Vector3(0.0, ref.speed_m_per_s, 0.0)
    };

    const auto trajectory = propagate_fixed_step(
        0.0,
        ref.period_s,
        10.0,
        initial_state,
        IntegrationMethod::classical_rk4,
        [](double t, const CartesianState& s) {
            return two_body_state_derivative(
                t, s, astradock::constants::earth_gravitational_parameter_m3_per_s2);
        });

    REQUIRE(trajectory.size() > 500);

    // Sample along the orbit and check LVLH properties at every step
    for (const auto& sample : trajectory) {
        const FrameBasis basis = compute_lvlh_basis(
            sample.state.position, sample.state.velocity);

        CHECK(is_orthonormal(basis, 1.0e-10, 1.0e-10));

        const Matrix3 dcm = dcm_lvlh_from_eci(basis);
        CHECK_THAT(dcm.determinant(), WithinAbs(1.0, 1.0e-10));

        // Radial vector in LVLH should be [ ||r||, 0, 0 ]
        const Vector3 r_lvlh = transform_eci_to_lvlh(sample.state.position, dcm);
        CHECK_THAT(r_lvlh.x(), WithinAbs(sample.state.position.norm(), 1.0e-6));
        CHECK_THAT(r_lvlh.y(), WithinAbs(0.0, 1.0e-6));
        CHECK_THAT(r_lvlh.z(), WithinAbs(0.0, 1.0e-6));

        // For circular orbit, radial velocity component is ~0 and along-track is ~speed
        const Vector3 v_lvlh = transform_eci_to_lvlh(sample.state.velocity, dcm);
        CHECK_THAT(v_lvlh.x(), WithinAbs(0.0, 1.0e-2));  // small radial velocity for near-circular
        CHECK_THAT(v_lvlh.y(), WithinAbs(ref.speed_m_per_s, 1.0e-2));
        CHECK_THAT(v_lvlh.z(), WithinAbs(0.0, 1.0e-6));  // zero out-of-plane velocity
    }
}
