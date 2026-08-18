#include "math/matrix3.hpp"
#include "math/vector3.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>

using astradock::math::approximately_equal;
using astradock::math::is_finite;
using astradock::math::Matrix3;
using astradock::math::Vector3;
using Catch::Matchers::WithinAbs;

TEST_CASE("Matrix3 default construction and factories initialize expected patterns", "[math][matrix3]") {
    const Matrix3 default_zero;
    const Matrix3 explicit_zero = Matrix3::zero();
    const Matrix3 identity = Matrix3::identity();

    for (std::size_t r = 0; r < 3; ++r) {
        for (std::size_t c = 0; c < 3; ++c) {
            CHECK(default_zero(r, c) == 0.0);
            CHECK(explicit_zero(r, c) == 0.0);
            if (r == c) {
                CHECK(identity(r, c) == 1.0);
            } else {
                CHECK(identity(r, c) == 0.0);
            }
        }
    }
}

TEST_CASE("Matrix3 element access and bounds checking behave correctly", "[math][matrix3]") {
    Matrix3 mat(
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0
    );

    CHECK(mat(0, 0) == 1.0);
    CHECK(mat(0, 1) == 2.0);
    CHECK(mat(0, 2) == 3.0);
    CHECK(mat(1, 0) == 4.0);
    CHECK(mat(1, 1) == 5.0);
    CHECK(mat(1, 2) == 6.0);
    CHECK(mat(2, 0) == 7.0);
    CHECK(mat(2, 1) == 8.0);
    CHECK(mat(2, 2) == 9.0);

    mat(1, 1) = 42.0;
    CHECK(mat(1, 1) == 42.0);

    CHECK_THROWS_AS(mat(3, 0), std::out_of_range);
    CHECK_THROWS_AS(mat(0, 3), std::out_of_range);
    CHECK_THROWS_AS(mat.row(3), std::out_of_range);
    CHECK_THROWS_AS(mat.column(3), std::out_of_range);
}

TEST_CASE("Matrix3 row and column construction and extraction match", "[math][matrix3]") {
    const Vector3 v0(1.0, 2.0, 3.0);
    const Vector3 v1(4.0, 5.0, 6.0);
    const Vector3 v2(7.0, 8.0, 9.0);

    const Matrix3 from_r = Matrix3::from_rows(v0, v1, v2);
    CHECK(from_r.row(0).x() == 1.0);
    CHECK(from_r.row(0).y() == 2.0);
    CHECK(from_r.row(0).z() == 3.0);
    CHECK(from_r.row(1).x() == 4.0);
    CHECK(from_r.row(2).z() == 9.0);

    const Matrix3 from_c = Matrix3::from_columns(v0, v1, v2);
    CHECK(from_c.column(0).x() == 1.0);
    CHECK(from_c.column(0).y() == 2.0);
    CHECK(from_c.column(0).z() == 3.0);
    CHECK(from_c.column(1).x() == 4.0);
    CHECK(from_c.column(2).z() == 9.0);

    // The transpose of from_rows should equal from_columns
    CHECK(approximately_equal(from_r.transpose(), from_c));
}

TEST_CASE("Matrix3 arithmetic operations satisfy algebraic identities", "[math][matrix3]") {
    const Matrix3 a(
        1.0, 2.0, 3.0,
        0.0, 1.0, 4.0,
        5.0, 6.0, 0.0
    );
    const Matrix3 b(
        2.0, 0.0, -1.0,
        1.0, 3.0,  2.0,
        0.0, 1.0,  1.0
    );

    const Matrix3 sum = a + b;
    CHECK(sum(0, 0) == 3.0);
    CHECK(sum(0, 1) == 2.0);
    CHECK(sum(0, 2) == 2.0);
    CHECK(sum(1, 0) == 1.0);
    CHECK(sum(1, 1) == 4.0);

    const Matrix3 diff = a - b;
    CHECK(diff(0, 0) == -1.0);
    CHECK(diff(0, 1) == 2.0);
    CHECK(diff(0, 2) == 4.0);

    const Matrix3 scaled = a * 2.0;
    CHECK(scaled(0, 0) == 2.0);
    CHECK(scaled(0, 2) == 6.0);
    CHECK(scaled(1, 2) == 8.0);

    const Matrix3 scalar_mult_left = 2.0 * a;
    CHECK(approximately_equal(scaled, scalar_mult_left));

    const Matrix3 div = scaled / 2.0;
    CHECK(approximately_equal(div, a));

    CHECK_THROWS_AS(a / 0.0, std::domain_error);
}

TEST_CASE("Matrix3 multiplication and identity properties hold", "[math][matrix3]") {
    const Matrix3 id = Matrix3::identity();
    const Vector3 v(3.0, -4.0, 5.0);

    // Identity * vector = vector
    const Vector3 iv = id * v;
    CHECK(iv.x() == v.x());
    CHECK(iv.y() == v.y());
    CHECK(iv.z() == v.z());

    const Matrix3 m(
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0
    );

    // M * Identity = M and Identity * M = M
    CHECK(approximately_equal(m * id, m));
    CHECK(approximately_equal(id * m, m));

    // Matrix vector product check
    const Vector3 mv = m * v;
    CHECK(mv.x() == 1.0 * 3.0 + 2.0 * (-4.0) + 3.0 * 5.0);  // 3 - 8 + 15 = 10
    CHECK(mv.y() == 4.0 * 3.0 + 5.0 * (-4.0) + 6.0 * 5.0);  // 12 - 20 + 30 = 22
    CHECK(mv.z() == 7.0 * 3.0 + 8.0 * (-4.0) + 9.0 * 5.0);  // 21 - 32 + 45 = 34
}

TEST_CASE("Matrix3 transpose and product transpose identity", "[math][matrix3]") {
    const Matrix3 a(
        1.0, 2.0, 3.0,
        0.0, 1.0, 4.0,
        5.0, 6.0, 0.0
    );
    const Matrix3 b(
        2.0, 0.0, -1.0,
        1.0, 3.0,  2.0,
        0.0, 1.0,  1.0
    );

    // (A^T)^T = A
    CHECK(approximately_equal(a.transpose().transpose(), a));

    // (A * B)^T = B^T * A^T
    const Matrix3 ab_t = (a * b).transpose();
    const Matrix3 bt_at = b.transpose() * a.transpose();
    CHECK(approximately_equal(ab_t, bt_at));
}

TEST_CASE("Matrix3 determinant and trace compute analytical values", "[math][matrix3]") {
    const Matrix3 id = Matrix3::identity();
    CHECK(id.determinant() == 1.0);
    CHECK(id.trace() == 3.0);

    const Matrix3 diag(
        2.0, 0.0, 0.0,
        0.0, 3.0, 0.0,
        0.0, 0.0, 4.0
    );
    CHECK(diag.determinant() == 24.0);
    CHECK(diag.trace() == 9.0);

    // Known 3x3 determinant:
    // [ 1, 2, 3 ]
    // [ 0, 1, 4 ]
    // [ 5, 6, 0 ]
    // det = 1*(0 - 24) - 2*(0 - 20) + 3*(0 - 5) = -24 + 40 - 15 = 1
    const Matrix3 test_mat(
        1.0, 2.0, 3.0,
        0.0, 1.0, 4.0,
        5.0, 6.0, 0.0
    );
    CHECK(test_mat.determinant() == 1.0);
    CHECK(test_mat.trace() == 2.0);
}

TEST_CASE("Matrix3 orthonormality and rotation properties", "[math][matrix3]") {
    // 90 degree rotation about Z axis
    const Matrix3 rot_z_90(
        0.0, -1.0, 0.0,
        1.0,  0.0, 0.0,
        0.0,  0.0, 1.0
    );

    CHECK(rot_z_90.is_orthonormal());
    CHECK_THAT(rot_z_90.determinant(), WithinAbs(1.0, 1.0e-14));

    // C * C^T = I
    const Matrix3 prod = rot_z_90 * rot_z_90.transpose();
    CHECK(approximately_equal(prod, Matrix3::identity()));

    // Norm preservation: ||C * v|| == ||v||
    const Vector3 v(3.0, 4.0, 5.0);
    const Vector3 rot_v = rot_z_90 * v;
    CHECK_THAT(rot_v.norm(), WithinAbs(v.norm(), 1.0e-12));
    CHECK_THAT(rot_v.x(), WithinAbs(-4.0, 1.0e-12));
    CHECK_THAT(rot_v.y(), WithinAbs(3.0, 1.0e-12));
    CHECK_THAT(rot_v.z(), WithinAbs(5.0, 1.0e-12));

    // Round-trip transformation: C^T * (C * v) = v
    const Vector3 round_trip = rot_z_90.transpose() * rot_v;
    CHECK(approximately_equal(round_trip, v));

    // Reflection matrix (det = -1) is orthogonal but NOT proper rotation
    const Matrix3 reflection(
        -1.0, 0.0, 0.0,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0
    );
    CHECK_THAT(reflection.determinant(), WithinAbs(-1.0, 1.0e-14));
    CHECK_FALSE(reflection.is_orthonormal());  // is_orthonormal requires det = +1
}

TEST_CASE("Matrix3 finite checks identify non-finite values", "[math][matrix3]") {
    const Matrix3 finite_mat = Matrix3::identity();
    CHECK(is_finite(finite_mat));

    Matrix3 nan_mat = Matrix3::identity();
    nan_mat(0, 1) = std::numeric_limits<double>::quiet_NaN();
    CHECK_FALSE(is_finite(nan_mat));

    Matrix3 inf_mat = Matrix3::identity();
    inf_mat(2, 2) = std::numeric_limits<double>::infinity();
    CHECK_FALSE(is_finite(inf_mat));
}
