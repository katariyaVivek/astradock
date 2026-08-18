#pragma once

#include "math/matrix3.hpp"
#include "math/vector3.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace astradock::math {

// A four-component, double-precision unit quaternion representing a 3D rotation.
//
// Convention:
//   Scalar-first: q = [w, x, y, z] = [w, v]
//   where w is the scalar (real) part and [x, y, z] is the vector (imaginary) part.
//
// Transformation / Action Convention:
//   A unit quaternion q rotates a physical 3D vector v via:
//     v' = q * [0, v] * q^(-1) = q * [0, v] * q*
//   This produces an active rotation matching the project's Matrix3 / DCM convention:
//     v' = C(q) * v
//   For frame transformations:
//     If q_A_B represents the orientation of frame B relative to frame A, it transforms
//     vector coordinates from frame B into frame A: v_A = q_A_B * v_B * q_A_B*.
class Quaternion {
public:
    // Default constructor initializes to the identity rotation: [1, 0, 0, 0].
    constexpr Quaternion() noexcept : w_(1.0), x_(0.0), y_(0.0), z_(0.0) {}

    // Direct component constructor (scalar-first: w, x, y, z).
    constexpr Quaternion(double w, double x, double y, double z) noexcept
        : w_(w), x_(x), y_(y), z_(z) {}

    [[nodiscard]] constexpr double w() const noexcept { return w_; }
    [[nodiscard]] constexpr double x() const noexcept { return x_; }
    [[nodiscard]] constexpr double y() const noexcept { return y_; }
    [[nodiscard]] constexpr double z() const noexcept { return z_; }

    [[nodiscard]] constexpr Vector3 vector_part() const noexcept {
        return {x_, y_, z_};
    }

    // Factory for the identity quaternion [1, 0, 0, 0].
    [[nodiscard]] static constexpr Quaternion identity() noexcept {
        return {1.0, 0.0, 0.0, 0.0};
    }

    // Constructs a unit quaternion from an axis and rotation angle in radians:
    //   q = [cos(theta/2), sin(theta/2) * u_hat]
    [[nodiscard]] static Quaternion from_axis_angle(const Vector3& axis, double angle_rad) {
        if (!std::isfinite(angle_rad)) {
            throw std::domain_error("Quaternion axis-angle: angle must be finite");
        }
        if (!is_finite(axis)) {
            throw std::domain_error("Quaternion axis-angle: axis components must be finite");
        }
        const double axis_norm = axis.norm();
        if (axis_norm == 0.0) {
            throw std::domain_error("Quaternion axis-angle: rotation axis magnitude cannot be zero");
        }

        const double half_angle = 0.5 * angle_rad;
        const double sin_half = std::sin(half_angle);
        const double cos_half = std::cos(half_angle);
        const Vector3 u_hat = axis / axis_norm;

        return {
            cos_half,
            sin_half * u_hat.x(),
            sin_half * u_hat.y(),
            sin_half * u_hat.z()
        };
    }

    // Squared Euclidean norm: w^2 + x^2 + y^2 + z^2.
    [[nodiscard]] constexpr double squared_norm() const noexcept {
        return w_ * w_ + x_ * x_ + y_ * y_ + z_ * z_;
    }

    // Euclidean norm (magnitude).
    [[nodiscard]] double norm() const noexcept {
        return std::sqrt(squared_norm());
    }

    // Checks if the quaternion is approximately a unit quaternion (||q|| ≈ 1).
    [[nodiscard]] bool is_unit(double tolerance = 1.0e-12) const noexcept {
        return std::abs(squared_norm() - 1.0) <= tolerance;
    }

    // Returns a normalized unit quaternion. Throws std::domain_error if norm is zero or non-finite.
    [[nodiscard]] Quaternion normalized(double min_norm_tolerance = 1.0e-14) const {
        const double n = norm();
        if (!std::isfinite(n) || n <= min_norm_tolerance) {
            throw std::domain_error("Cannot normalize a zero or near-zero quaternion");
        }
        return {w_ / n, x_ / n, y_ / n, z_ / n};
    }

    // Quaternion conjugate: q* = [w, -x, -y, -z].
    [[nodiscard]] constexpr Quaternion conjugate() const noexcept {
        return {w_, -x_, -y_, -z_};
    }

    // Quaternion inverse: q^(-1) = q* / ||q||^2.
    // For unit quaternions, q^(-1) == q*.
    [[nodiscard]] Quaternion inverse(double min_norm_tolerance = 1.0e-14) const {
        const double sq_norm = squared_norm();
        if (!std::isfinite(sq_norm) || sq_norm <= min_norm_tolerance) {
            throw std::domain_error("Cannot invert a zero or near-zero quaternion");
        }
        return {w_ / sq_norm, -x_ / sq_norm, -y_ / sq_norm, -z_ / sq_norm};
    }

    // Hamilton product: q_result = (*this) ⊗ other.
    // Non-commutative composition of rotations.
    [[nodiscard]] constexpr Quaternion operator*(const Quaternion& other) const noexcept {
        return {
            w_ * other.w_ - x_ * other.x_ - y_ * other.y_ - z_ * other.z_,
            w_ * other.x_ + x_ * other.w_ + y_ * other.z_ - z_ * other.y_,
            w_ * other.y_ - x_ * other.z_ + y_ * other.w_ + z_ * other.x_,
            w_ * other.z_ + x_ * other.y_ - y_ * other.x_ + z_ * other.w_
        };
    }

    // Rotates a 3D vector v via active rotation:
    //   v' = v + 2 * w * (u x v) + 2 * (u x (u x v))
    // where u = [x, y, z]. This is mathematically equivalent to q * [0, v] * q*.
    [[nodiscard]] Vector3 rotate_vector(const Vector3& vector) const {
        if (!is_finite(vector)) {
            throw std::domain_error("Cannot rotate non-finite vector");
        }
        const Vector3 u(x_, y_, z_);
        const Vector3 u_cross_v = u.cross(vector);
        const Vector3 u_cross_u_cross_v = u.cross(u_cross_v);

        return vector + (u_cross_v * (2.0 * w_)) + (u_cross_u_cross_v * 2.0);
    }

    // Converts this unit quaternion to an equivalent 3x3 Direction Cosine Matrix (DCM).
    // The matrix satisfies v' = C * v, matching rotate_vector().
    [[nodiscard]] Matrix3 to_rotation_matrix() const {
        const double xx = x_ * x_;
        const double yy = y_ * y_;
        const double zz = z_ * z_;

        const double wx = w_ * x_;
        const double wy = w_ * y_;
        const double wz = w_ * z_;

        const double xy = x_ * y_;
        const double xz = x_ * z_;
        const double yz = y_ * z_;

        return Matrix3(
            1.0 - 2.0 * (yy + zz), 2.0 * (xy - wz),       2.0 * (xz + wy),
            2.0 * (xy + wz),       1.0 - 2.0 * (xx + zz), 2.0 * (yz - wx),
            2.0 * (xz - wy),       2.0 * (yz + wx),       1.0 - 2.0 * (xx + yy)
        );
    }

    // Converts a 3x3 Direction Cosine Matrix into a unit quaternion using
    // Shepperd's numerically stable algorithm (avoids division by near-zero values).
    [[nodiscard]] static Quaternion from_rotation_matrix(const Matrix3& matrix) {
        if (!is_finite(matrix)) {
            throw std::domain_error("Cannot create quaternion from non-finite matrix");
        }

        const double r00 = matrix(0, 0);
        const double r11 = matrix(1, 1);
        const double r22 = matrix(2, 2);
        const double trace = matrix.trace();

        double w = 0.0;
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;

        if (trace > r00 && trace > r11 && trace > r22) {
            // Case 0: Trace is largest -> compute w first
            w = 0.5 * std::sqrt(1.0 + trace);
            const double s = 0.25 / w;
            x = (matrix(2, 1) - matrix(1, 2)) * s;
            y = (matrix(0, 2) - matrix(2, 0)) * s;
            z = (matrix(1, 0) - matrix(0, 1)) * s;
        } else if (r00 > r11 && r00 > r22) {
            // Case 1: r00 is largest -> compute x first
            x = 0.5 * std::sqrt(1.0 + 2.0 * r00 - trace);
            const double s = 0.25 / x;
            w = (matrix(2, 1) - matrix(1, 2)) * s;
            y = (matrix(0, 1) + matrix(1, 0)) * s;
            z = (matrix(0, 2) + matrix(2, 0)) * s;
        } else if (r11 > r22) {
            // Case 2: r11 is largest -> compute y first
            y = 0.5 * std::sqrt(1.0 + 2.0 * r11 - trace);
            const double s = 0.25 / y;
            w = (matrix(0, 2) - matrix(2, 0)) * s;
            x = (matrix(0, 1) + matrix(1, 0)) * s;
            z = (matrix(1, 2) + matrix(2, 1)) * s;
        } else {
            // Case 3: r22 is largest -> compute z first
            z = 0.5 * std::sqrt(1.0 + 2.0 * r22 - trace);
            const double s = 0.25 / z;
            w = (matrix(1, 0) - matrix(0, 1)) * s;
            x = (matrix(0, 2) + matrix(2, 0)) * s;
            y = (matrix(1, 2) + matrix(2, 1)) * s;
        }

        const Quaternion q(w, x, y, z);
        return q.normalized();
    }

    // Negation operator: -q (double cover).
    [[nodiscard]] constexpr Quaternion operator-() const noexcept {
        return {-w_, -x_, -y_, -z_};
    }

    // Scalar arithmetic operators
    [[nodiscard]] constexpr Quaternion operator+(const Quaternion& other) const noexcept {
        return {w_ + other.w_, x_ + other.x_, y_ + other.y_, z_ + other.z_};
    }

    [[nodiscard]] constexpr Quaternion operator-(const Quaternion& other) const noexcept {
        return {w_ - other.w_, x_ - other.x_, y_ - other.y_, z_ - other.z_};
    }

    [[nodiscard]] constexpr Quaternion operator*(double scalar) const noexcept {
        return {w_ * scalar, x_ * scalar, y_ * scalar, z_ * scalar};
    }

    [[nodiscard]] Quaternion operator/(double scalar) const {
        if (scalar == 0.0) {
            throw std::domain_error("Quaternion division by zero");
        }
        return {w_ / scalar, x_ / scalar, y_ / scalar, z_ / scalar};
    }

private:
    double w_{1.0};
    double x_{0.0};
    double y_{0.0};
    double z_{0.0};
};

[[nodiscard]] inline bool is_finite(const Quaternion& q) noexcept {
    return std::isfinite(q.w()) && std::isfinite(q.x()) && std::isfinite(q.y()) && std::isfinite(q.z());
}

// Checks if two quaternions are element-wise approximately equal: q1 ≈ q2.
[[nodiscard]] inline bool approximately_equal(
    const Quaternion& lhs,
    const Quaternion& rhs,
    double absolute_tolerance = 1.0e-12,
    double relative_tolerance = 1.0e-12) noexcept {
    return approximately_equal(lhs.w(), rhs.w(), absolute_tolerance, relative_tolerance)
        && approximately_equal(lhs.x(), rhs.x(), absolute_tolerance, relative_tolerance)
        && approximately_equal(lhs.y(), rhs.y(), absolute_tolerance, relative_tolerance)
        && approximately_equal(lhs.z(), rhs.z(), absolute_tolerance, relative_tolerance);
}

// Checks if two quaternions represent the same physical rotation (Double Cover: q ≈ +q' or q ≈ -q').
[[nodiscard]] inline bool represents_same_rotation(
    const Quaternion& lhs,
    const Quaternion& rhs,
    double absolute_tolerance = 1.0e-12,
    double relative_tolerance = 1.0e-12) noexcept {
    return approximately_equal(lhs, rhs, absolute_tolerance, relative_tolerance)
        || approximately_equal(lhs, -rhs, absolute_tolerance, relative_tolerance);
}

}  // namespace astradock::math
