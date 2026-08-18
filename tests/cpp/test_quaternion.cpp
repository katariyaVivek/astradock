#include "attitude/attitude_state.hpp"
#include "math/constants.hpp"
#include "math/euler_angles.hpp"
#include "math/matrix3.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <random>

using astradock::attitude::AttitudeState;
using astradock::constants::pi;
using astradock::math::approximately_equal;
using astradock::math::EulerAngles;
using astradock::math::euler_to_quaternion;
using astradock::math::euler_to_rotation_matrix;
using astradock::math::is_finite;
using astradock::math::Matrix3;
using astradock::math::Quaternion;
using astradock::math::quaternion_to_euler;
using astradock::math::represents_same_rotation;
using astradock::math::rotation_matrix_to_euler;
using astradock::math::Vector3;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// ===========================================================================
// Test 1: Quaternion Construction, Norm, Normalization, and Conjugate
// ===========================================================================
TEST_CASE("Quaternion basic primitives satisfy algebraic properties", "[attitude][quaternion]") {
    const Quaternion q_id = Quaternion::identity();
    CHECK_THAT(q_id.w(), WithinAbs(1.0, 1.0e-15));
    CHECK_THAT(q_id.x(), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(q_id.y(), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(q_id.z(), WithinAbs(0.0, 1.0e-15));
    CHECK(q_id.is_unit());

    const Quaternion q(1.0, 2.0, 3.0, 4.0);
    // norm = sqrt(1 + 4 + 9 + 16) = sqrt(30)
    const double expected_norm = std::sqrt(30.0);
    CHECK_THAT(q.norm(), WithinRel(expected_norm, 1.0e-14));
    CHECK_THAT(q.squared_norm(), WithinRel(30.0, 1.0e-14));
    CHECK_FALSE(q.is_unit());

    const Quaternion q_norm = q.normalized();
    CHECK(q_norm.is_unit());
    CHECK_THAT(q_norm.norm(), WithinRel(1.0, 1.0e-14));
    CHECK_THAT(q_norm.w(), WithinRel(1.0 / expected_norm, 1.0e-14));
    CHECK_THAT(q_norm.x(), WithinRel(2.0 / expected_norm, 1.0e-14));
    CHECK_THAT(q_norm.y(), WithinRel(3.0 / expected_norm, 1.0e-14));
    CHECK_THAT(q_norm.z(), WithinRel(4.0 / expected_norm, 1.0e-14));

    const Quaternion q_conj = q.conjugate();
    CHECK_THAT(q_conj.w(), WithinAbs(1.0, 1.0e-15));
    CHECK_THAT(q_conj.x(), WithinAbs(-2.0, 1.0e-15));
    CHECK_THAT(q_conj.y(), WithinAbs(-3.0, 1.0e-15));
    CHECK_THAT(q_conj.z(), WithinAbs(-4.0, 1.0e-15));

    // Double conjugate is original
    CHECK(approximately_equal(q_conj.conjugate(), q));
}

// ===========================================================================
// Test 2: Hamilton Product, Identity, Non-Commutativity, and Inverses
// ===========================================================================
TEST_CASE("Hamilton product and inverse operations behave correctly", "[attitude][quaternion]") {
    const Quaternion q1 = Quaternion::from_axis_angle(Vector3(1.0, 0.0, 0.0), 0.5);
    const Quaternion q2 = Quaternion::from_axis_angle(Vector3(0.0, 1.0, 0.0), 0.8);
    const Quaternion q_id = Quaternion::identity();

    // Identity multiplication: q * I = q and I * q = q
    CHECK(approximately_equal(q1 * q_id, q1));
    CHECK(approximately_equal(q_id * q1, q1));

    // Non-commutativity: q1 * q2 != q2 * q1 for orthogonal axes
    const Quaternion p12 = q1 * q2;
    const Quaternion p21 = q2 * q1;
    CHECK_FALSE(approximately_equal(p12, p21, 1.0e-4, 1.0e-4));

    // Norm multiplication property: ||q1 * q2|| == ||q1|| * ||q2||
    CHECK_THAT(p12.norm(), WithinRel(1.0, 1.0e-14));

    // General inverse: q * q^(-1) == q_id
    const Quaternion q_arb(2.0, -1.0, 3.0, -0.5);
    const Quaternion q_inv = q_arb.inverse();
    const Quaternion prod_right = q_arb * q_inv;
    const Quaternion prod_left = q_inv * q_arb;
    CHECK(approximately_equal(prod_right, q_id, 1.0e-12, 1.0e-12));
    CHECK(approximately_equal(prod_left, q_id, 1.0e-12, 1.0e-12));

    // For unit quaternion: q^(-1) == q*
    CHECK(approximately_equal(q1.inverse(), q1.conjugate(), 1.0e-14, 1.0e-14));
}

// ===========================================================================
// Test 3: Known Canonical Rotations (90 deg, 180 deg around X, Y, Z)
// ===========================================================================
TEST_CASE("Known canonical rotations transform orthogonal basis vectors", "[attitude][quaternion]") {
    const double half_pi = 0.5 * pi;

    // 1. 90 deg about +X: rotates [0, 1, 0] to [0, 0, 1]
    {
        const Quaternion q = Quaternion::from_axis_angle(Vector3(1.0, 0.0, 0.0), half_pi);
        const Vector3 v(0.0, 1.0, 0.0);
        const Vector3 v_rot = q.rotate_vector(v);
        CHECK(approximately_equal(v_rot, Vector3(0.0, 0.0, 1.0), 1.0e-14, 1.0e-14));
    }

    // 2. 90 deg about +Y: rotates [0, 0, 1] to [1, 0, 0]
    {
        const Quaternion q = Quaternion::from_axis_angle(Vector3(0.0, 1.0, 0.0), half_pi);
        const Vector3 v(0.0, 0.0, 1.0);
        const Vector3 v_rot = q.rotate_vector(v);
        CHECK(approximately_equal(v_rot, Vector3(1.0, 0.0, 0.0), 1.0e-14, 1.0e-14));
    }

    // 3. 90 deg about +Z: rotates [1, 0, 0] to [0, 1, 0]
    {
        const Quaternion q = Quaternion::from_axis_angle(Vector3(0.0, 0.0, 1.0), half_pi);
        const Vector3 v(1.0, 0.0, 0.0);
        const Vector3 v_rot = q.rotate_vector(v);
        CHECK(approximately_equal(v_rot, Vector3(0.0, 1.0, 0.0), 1.0e-14, 1.0e-14));
    }

    // 4. 180 deg about +X: rotates [0, 1, 0] to [0, -1, 0]
    {
        const Quaternion q = Quaternion::from_axis_angle(Vector3(1.0, 0.0, 0.0), pi);
        const Vector3 v(0.0, 1.0, 0.0);
        const Vector3 v_rot = q.rotate_vector(v);
        CHECK(approximately_equal(v_rot, Vector3(0.0, -1.0, 0.0), 1.0e-14, 1.0e-14));
    }

    // 5. 180 deg about +Y: rotates [1, 0, 0] to [-1, 0, 0]
    {
        const Quaternion q = Quaternion::from_axis_angle(Vector3(0.0, 1.0, 0.0), pi);
        const Vector3 v(1.0, 0.0, 0.0);
        const Vector3 v_rot = q.rotate_vector(v);
        CHECK(approximately_equal(v_rot, Vector3(-1.0, 0.0, 0.0), 1.0e-14, 1.0e-14));
    }

    // 6. 180 deg about +Z: rotates [1, 0, 0] to [-1, 0, 0]
    {
        const Quaternion q = Quaternion::from_axis_angle(Vector3(0.0, 0.0, 1.0), pi);
        const Vector3 v(1.0, 0.0, 0.0);
        const Vector3 v_rot = q.rotate_vector(v);
        CHECK(approximately_equal(v_rot, Vector3(-1.0, 0.0, 0.0), 1.0e-14, 1.0e-14));
    }
}

// ===========================================================================
// Test 4: Double-Cover Property (q and -q physical equivalence)
// ===========================================================================
TEST_CASE("Double-cover property guarantees physical equivalence of q and -q", "[attitude][quaternion]") {
    const Quaternion q = Quaternion::from_axis_angle(Vector3(1.0, 2.0, 3.0), 1.25);
    const Quaternion q_neg = -q;

    // Element-wise they are opposites
    CHECK_THAT(q.w() + q_neg.w(), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(q.x() + q_neg.x(), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(q.y() + q_neg.y(), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(q.z() + q_neg.z(), WithinAbs(0.0, 1.0e-15));

    // Represents same rotation utility
    CHECK(represents_same_rotation(q, q_neg));

    // Both produce identical vector rotations
    const Vector3 v(123.456, -789.012, 345.678);
    const Vector3 v_rot1 = q.rotate_vector(v);
    const Vector3 v_rot2 = q_neg.rotate_vector(v);
    CHECK(approximately_equal(v_rot1, v_rot2, 1.0e-12, 1.0e-12));

    // Both produce identical Direction Cosine Matrices
    const Matrix3 mat1 = q.to_rotation_matrix();
    const Matrix3 mat2 = q_neg.to_rotation_matrix();
    CHECK(approximately_equal(mat1, mat2, 1.0e-14, 1.0e-14));
}

// ===========================================================================
// Test 5: Quaternion <-> DCM Conversion and Round-Trip Robustness
// ===========================================================================
TEST_CASE("Quaternion and DCM conversions satisfy round-trip accuracy across regimes", "[attitude][quaternion]") {
    const std::vector<std::pair<Vector3, double>> test_rotations = {
        {Vector3(1.0, 0.0, 0.0), 0.0},                  // Identity (0 deg)
        {Vector3(1.0, 0.0, 0.0), 0.5 * pi},             // 90 deg X
        {Vector3(0.0, 1.0, 0.0), 0.5 * pi},             // 90 deg Y
        {Vector3(0.0, 0.0, 1.0), 0.5 * pi},             // 90 deg Z
        {Vector3(1.0, 0.0, 0.0), pi},                    // 180 deg X (Branch case 1)
        {Vector3(0.0, 1.0, 0.0), pi},                    // 180 deg Y (Branch case 2)
        {Vector3(0.0, 0.0, 1.0), pi},                    // 180 deg Z (Branch case 3)
        {Vector3(1.0, 1.0, 1.0), 0.25 * pi},            // 45 deg Arbitrary
        {Vector3(2.0, -1.0, 3.0), 179.99 * pi / 180.0}, // Near-180 deg Arbitrary
    };

    for (const auto& [axis, angle] : test_rotations) {
        const Quaternion q_orig = Quaternion::from_axis_angle(axis, angle);
        const Matrix3 dcm = q_orig.to_rotation_matrix();

        // Check DCM Orthonormality and proper rotation
        CHECK(dcm.is_orthonormal(1.0e-12));
        CHECK_THAT(dcm.determinant(), WithinAbs(1.0, 1.0e-12));

        // Reconstruct Quaternion from DCM (Shepperd algorithm)
        const Quaternion q_recon = Quaternion::from_rotation_matrix(dcm);
        CHECK(represents_same_rotation(q_orig, q_recon, 1.0e-12, 1.0e-12));

        // Reconstruct DCM from reconstructed Quaternion
        const Matrix3 dcm_recon = q_recon.to_rotation_matrix();
        CHECK(approximately_equal(dcm, dcm_recon, 1.0e-12, 1.0e-12));

        // Test vector rotation consistency between Quaternion and DCM
        const Vector3 v_sample(42.0, -17.5, 99.1);
        const Vector3 v_rot_q = q_orig.rotate_vector(v_sample);
        const Vector3 v_rot_dcm = dcm * v_sample;
        CHECK(approximately_equal(v_rot_q, v_rot_dcm, 1.0e-12, 1.0e-12));
    }
}

// ===========================================================================
// Test 6: Composition Consistency (Quaternion vs Matrix Multiplication)
// ===========================================================================
TEST_CASE("Rotation composition order matches Direction Cosine Matrix multiplication", "[attitude][quaternion]") {
    const Quaternion q1 = Quaternion::from_axis_angle(Vector3(1.0, 2.0, 0.5), 0.73);
    const Quaternion q2 = Quaternion::from_axis_angle(Vector3(-1.0, 0.5, 3.0), 1.42);

    // Composite Quaternion: q_composite = q1 * q2
    const Quaternion q_comp = q1 * q2;

    // Independent DCMs
    const Matrix3 c1 = q1.to_rotation_matrix();
    const Matrix3 c2 = q2.to_rotation_matrix();
    const Matrix3 c_comp_mat = c1 * c2; // Independent matrix product

    const Matrix3 c_comp_from_q = q_comp.to_rotation_matrix();

    CHECK(approximately_equal(c_comp_from_q, c_comp_mat, 1.0e-14, 1.0e-14));
}

// ===========================================================================
// Test 7: Randomized Scale-Aware Vector Norm Preservation Property
// ===========================================================================
TEST_CASE("Rotations preserve vector norm across wide scale variations", "[attitude][quaternion]") {
    std::mt19937_64 rng(42); // Fixed deterministic seed
    std::uniform_real_distribution<double> dist_axis(-1.0, 1.0);
    std::uniform_real_distribution<double> dist_angle(-pi, pi);

    const std::vector<double> scales = {1.0e-9, 1.0e-6, 1.0e-3, 1.0, 1.0e3, 1.0e6};

    for (int iter = 0; iter < 50; ++iter) {
        const Vector3 axis(dist_axis(rng), dist_axis(rng), dist_axis(rng));
        if (axis.norm() < 1.0e-3) continue;
        const double angle = dist_angle(rng);

        const Quaternion q = Quaternion::from_axis_angle(axis, angle);

        for (const double s : scales) {
            const Vector3 v_orig(dist_axis(rng) * s, dist_axis(rng) * s, dist_axis(rng) * s);
            const double orig_norm = v_orig.norm();

            const Vector3 v_rot = q.rotate_vector(v_orig);
            const double rot_norm = v_rot.norm();

            CHECK_THAT(rot_norm, WithinRel(orig_norm, 1.0e-12));
        }
    }
}

// ===========================================================================
// Test 8: ZYX Euler Angles and Gimbal Lock Detection
// ===========================================================================
TEST_CASE("Euler angle conversions and gimbal lock handling behave as documented", "[attitude][euler]") {
    // 1. Nominal case: yaw = 30 deg, pitch = 45 deg, roll = 60 deg
    const EulerAngles euler_nominal{30.0 * pi / 180.0, 45.0 * pi / 180.0, 60.0 * pi / 180.0};
    const Matrix3 dcm = euler_to_rotation_matrix(euler_nominal);
    const EulerAngles euler_recon = rotation_matrix_to_euler(dcm);

    CHECK_THAT(euler_recon.yaw_rad, WithinAbs(euler_nominal.yaw_rad, 1.0e-12));
    CHECK_THAT(euler_recon.pitch_rad, WithinAbs(euler_nominal.pitch_rad, 1.0e-12));
    CHECK_THAT(euler_recon.roll_rad, WithinAbs(euler_nominal.roll_rad, 1.0e-12));

    const Quaternion q_euler = euler_to_quaternion(euler_nominal);
    const EulerAngles euler_from_q = quaternion_to_euler(q_euler);
    CHECK_THAT(euler_from_q.yaw_rad, WithinAbs(euler_nominal.yaw_rad, 1.0e-12));
    CHECK_THAT(euler_from_q.pitch_rad, WithinAbs(euler_nominal.pitch_rad, 1.0e-12));
    CHECK_THAT(euler_from_q.roll_rad, WithinAbs(euler_nominal.roll_rad, 1.0e-12));

    // 2. Gimbal Lock at pitch = +90 deg
    const EulerAngles euler_gimbal{0.0, 0.5 * pi, 45.0 * pi / 180.0};
    const Matrix3 dcm_gimbal = euler_to_rotation_matrix(euler_gimbal);
    const EulerAngles euler_gimbal_recon = rotation_matrix_to_euler(dcm_gimbal);

    CHECK_THAT(euler_gimbal_recon.pitch_rad, WithinAbs(0.5 * pi, 1.0e-12));
    CHECK_THAT(euler_gimbal_recon.yaw_rad, WithinAbs(0.0, 1.0e-12)); // Resolved by convention
    CHECK_THAT(euler_gimbal_recon.roll_rad, WithinAbs(45.0 * pi / 180.0, 1.0e-12));
}

// ===========================================================================
// Test 9: Defensive Validation and Degenerate Input Handling
// ===========================================================================
TEST_CASE("Defensive validation rejects degenerate and non-finite inputs", "[attitude][quaternion]") {
    // 1. Zero axis in axis-angle
    CHECK_THROWS_AS(Quaternion::from_axis_angle(Vector3(0.0, 0.0, 0.0), 1.0), std::domain_error);

    // 2. Non-finite angle in axis-angle
    CHECK_THROWS_AS(Quaternion::from_axis_angle(Vector3(1.0, 0.0, 0.0), std::numeric_limits<double>::infinity()), std::domain_error);

    // 3. Non-finite vector in axis-angle
    CHECK_THROWS_AS(Quaternion::from_axis_angle(Vector3(std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0), 1.0), std::domain_error);

    // 4. Zero quaternion normalization and inversion
    const Quaternion q_zero(0.0, 0.0, 0.0, 0.0);
    CHECK_THROWS_AS(q_zero.normalized(), std::domain_error);
    CHECK_THROWS_AS(q_zero.inverse(), std::domain_error);

    // 5. Non-finite matrix to quaternion
    Matrix3 non_finite_mat = Matrix3::identity();
    non_finite_mat(0, 0) = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROWS_AS(Quaternion::from_rotation_matrix(non_finite_mat), std::domain_error);

    // 6. AttitudeState instantiation
    const AttitudeState state;
    CHECK(state.orientation.is_unit());
}
