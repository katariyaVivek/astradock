#pragma once

#include "math/vector3.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace astradock::math {

// A 3x3 double-precision matrix stored in row-major order:
//   [ m00, m01, m02 ]   indices: [ 0, 1, 2 ]
//   [ m10, m11, m12 ]            [ 3, 4, 5 ]
//   [ m20, m21, m22 ]            [ 6, 7, 8 ]
//
// Like Vector3, Matrix3 carries no implicit coordinate frame or convention.
// When used as a Direction Cosine Matrix (DCM), the frame transformation
// direction (e.g., C_to_from) must be explicitly stated by the calling API.
class Matrix3 {
public:
    constexpr Matrix3() noexcept : data_{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0} {}

    constexpr Matrix3(
        double m00, double m01, double m02,
        double m10, double m11, double m12,
        double m20, double m21, double m22) noexcept
        : data_{m00, m01, m02, m10, m11, m12, m20, m21, m22} {}

    // Constructs a matrix with rows set to the three given vectors.
    constexpr Matrix3(const Vector3& r0, const Vector3& r1, const Vector3& r2) noexcept
        : data_{
            r0.x(), r0.y(), r0.z(),
            r1.x(), r1.y(), r1.z(),
            r2.x(), r2.y(), r2.z()} {}

    [[nodiscard]] static constexpr Matrix3 zero() noexcept {
        return Matrix3{};
    }

    [[nodiscard]] static constexpr Matrix3 identity() noexcept {
        return Matrix3{
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0
        };
    }

    // Constructs a Matrix3 whose rows are the three given vectors:
    //   row 0 = r0, row 1 = r1, row 2 = r2
    [[nodiscard]] static constexpr Matrix3 from_rows(
        const Vector3& r0,
        const Vector3& r1,
        const Vector3& r2) noexcept {
        return Matrix3{r0, r1, r2};
    }

    // Constructs a Matrix3 whose columns are the three given vectors:
    //   col 0 = c0, col 1 = c1, col 2 = c2
    [[nodiscard]] static constexpr Matrix3 from_columns(
        const Vector3& c0,
        const Vector3& c1,
        const Vector3& c2) noexcept {
        return Matrix3{
            c0.x(), c1.x(), c2.x(),
            c0.y(), c1.y(), c2.y(),
            c0.z(), c1.z(), c2.z()
        };
    }

    // Element access by 0-indexed row and column.
    [[nodiscard]] constexpr double operator()(std::size_t row, std::size_t col) const {
        if (row >= 3 || col >= 3) {
            throw std::out_of_range("Matrix3 index out of range");
        }
        return data_[row * 3 + col];
    }

    [[nodiscard]] constexpr double& operator()(std::size_t row, std::size_t col) {
        if (row >= 3 || col >= 3) {
            throw std::out_of_range("Matrix3 index out of range");
        }
        return data_[row * 3 + col];
    }

    [[nodiscard]] constexpr Vector3 row(std::size_t r) const {
        if (r >= 3) {
            throw std::out_of_range("Matrix3 row index out of range");
        }
        return {data_[r * 3], data_[r * 3 + 1], data_[r * 3 + 2]};
    }

    [[nodiscard]] constexpr Vector3 column(std::size_t c) const {
        if (c >= 3) {
            throw std::out_of_range("Matrix3 column index out of range");
        }
        return {data_[c], data_[3 + c], data_[6 + c]};
    }

    [[nodiscard]] constexpr Matrix3 operator+(const Matrix3& other) const noexcept {
        return Matrix3{
            data_[0] + other.data_[0], data_[1] + other.data_[1], data_[2] + other.data_[2],
            data_[3] + other.data_[3], data_[4] + other.data_[4], data_[5] + other.data_[5],
            data_[6] + other.data_[6], data_[7] + other.data_[7], data_[8] + other.data_[8]
        };
    }

    [[nodiscard]] constexpr Matrix3 operator-(const Matrix3& other) const noexcept {
        return Matrix3{
            data_[0] - other.data_[0], data_[1] - other.data_[1], data_[2] - other.data_[2],
            data_[3] - other.data_[3], data_[4] - other.data_[4], data_[5] - other.data_[5],
            data_[6] - other.data_[6], data_[7] - other.data_[7], data_[8] - other.data_[8]
        };
    }

    [[nodiscard]] constexpr Matrix3 operator*(double scalar) const noexcept {
        return Matrix3{
            data_[0] * scalar, data_[1] * scalar, data_[2] * scalar,
            data_[3] * scalar, data_[4] * scalar, data_[5] * scalar,
            data_[6] * scalar, data_[7] * scalar, data_[8] * scalar
        };
    }

    [[nodiscard]] constexpr Matrix3 operator/(double scalar) const {
        if (scalar == 0.0) {
            throw std::domain_error("Matrix3 division by zero");
        }
        return Matrix3{
            data_[0] / scalar, data_[1] / scalar, data_[2] / scalar,
            data_[3] / scalar, data_[4] / scalar, data_[5] / scalar,
            data_[6] / scalar, data_[7] / scalar, data_[8] / scalar
        };
    }

    // Matrix-vector multiplication: result = M * v
    // Equivalent to dot product of each row with the vector v.
    [[nodiscard]] constexpr Vector3 operator*(const Vector3& v) const noexcept {
        return {
            data_[0] * v.x() + data_[1] * v.y() + data_[2] * v.z(),
            data_[3] * v.x() + data_[4] * v.y() + data_[5] * v.z(),
            data_[6] * v.x() + data_[7] * v.y() + data_[8] * v.z()
        };
    }

    // Matrix-matrix multiplication: result = M * other
    [[nodiscard]] constexpr Matrix3 operator*(const Matrix3& other) const noexcept {
        return Matrix3{
            data_[0] * other.data_[0] + data_[1] * other.data_[3] + data_[2] * other.data_[6],
            data_[0] * other.data_[1] + data_[1] * other.data_[4] + data_[2] * other.data_[7],
            data_[0] * other.data_[2] + data_[1] * other.data_[5] + data_[2] * other.data_[8],

            data_[3] * other.data_[0] + data_[4] * other.data_[3] + data_[5] * other.data_[6],
            data_[3] * other.data_[1] + data_[4] * other.data_[4] + data_[5] * other.data_[7],
            data_[3] * other.data_[2] + data_[4] * other.data_[5] + data_[5] * other.data_[8],

            data_[6] * other.data_[0] + data_[7] * other.data_[3] + data_[8] * other.data_[6],
            data_[6] * other.data_[1] + data_[7] * other.data_[4] + data_[8] * other.data_[7],
            data_[6] * other.data_[2] + data_[7] * other.data_[5] + data_[8] * other.data_[8]
        };
    }

    // Matrix transpose: M^T
    [[nodiscard]] constexpr Matrix3 transpose() const noexcept {
        return Matrix3{
            data_[0], data_[3], data_[6],
            data_[1], data_[4], data_[7],
            data_[2], data_[5], data_[8]
        };
    }

    // Matrix determinant: det(M)
    [[nodiscard]] constexpr double determinant() const noexcept {
        return data_[0] * (data_[4] * data_[8] - data_[5] * data_[7])
             - data_[1] * (data_[3] * data_[8] - data_[5] * data_[6])
             + data_[2] * (data_[3] * data_[7] - data_[4] * data_[6]);
    }

    // Matrix trace: sum of diagonal elements
    [[nodiscard]] constexpr double trace() const noexcept {
        return data_[0] + data_[4] + data_[8];
    }

    // Checks if the matrix represents an orthonormal rotation matrix:
    //   1. M * M^T ≈ I
    //   2. det(M) ≈ +1.0 (proper rotation, right-handed)
    [[nodiscard]] bool is_orthonormal(
        double absolute_tolerance = 1.0e-12,
        double relative_tolerance = 1.0e-12) const noexcept;

private:
    std::array<double, 9> data_;
};

[[nodiscard]] constexpr Matrix3 operator*(double scalar, const Matrix3& m) noexcept {
    return m * scalar;
}

[[nodiscard]] inline bool approximately_equal(
    const Matrix3& lhs,
    const Matrix3& rhs,
    double absolute_tolerance = 1.0e-12,
    double relative_tolerance = 1.0e-12) noexcept {
    for (std::size_t r = 0; r < 3; ++r) {
        for (std::size_t c = 0; c < 3; ++c) {
            if (!approximately_equal(lhs(r, c), rhs(r, c), absolute_tolerance, relative_tolerance)) {
                return false;
            }
        }
    }
    return true;
}

inline bool Matrix3::is_orthonormal(
    double absolute_tolerance,
    double relative_tolerance) const noexcept {
    const Matrix3 prod = (*this) * transpose();
    const Matrix3 id = identity();
    if (!approximately_equal(prod, id, absolute_tolerance, relative_tolerance)) {
        return false;
    }
    return approximately_equal(determinant(), 1.0, absolute_tolerance, relative_tolerance);
}

[[nodiscard]] inline bool is_finite(const Matrix3& m) noexcept {
    for (std::size_t r = 0; r < 3; ++r) {
        for (std::size_t c = 0; c < 3; ++c) {
            if (!std::isfinite(m(r, c))) {
                return false;
            }
        }
    }
    return true;
}

}  // namespace astradock::math
