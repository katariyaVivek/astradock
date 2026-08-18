#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace astradock::math {

// A three-component, double-precision vector. The type deliberately carries no
// implicit coordinate frame or unit; callers must document both at each API
// boundary until the project has a justified strongly typed frame system.
class Vector3 {
public:
    constexpr Vector3() noexcept = default;
    constexpr Vector3(double x, double y, double z) noexcept : x_(x), y_(y), z_(z) {}

    [[nodiscard]] constexpr double x() const noexcept { return x_; }
    [[nodiscard]] constexpr double y() const noexcept { return y_; }
    [[nodiscard]] constexpr double z() const noexcept { return z_; }

    [[nodiscard]] constexpr Vector3 operator+(const Vector3& other) const noexcept {
        return {x_ + other.x_, y_ + other.y_, z_ + other.z_};
    }

    [[nodiscard]] constexpr Vector3 operator-(const Vector3& other) const noexcept {
        return {x_ - other.x_, y_ - other.y_, z_ - other.z_};
    }

    [[nodiscard]] constexpr Vector3 operator*(double scalar) const noexcept {
        return {x_ * scalar, y_ * scalar, z_ * scalar};
    }

    [[nodiscard]] constexpr Vector3 operator/(double scalar) const {
        if (scalar == 0.0) {
            throw std::domain_error("Vector3 division by zero");
        }
        return {x_ / scalar, y_ / scalar, z_ / scalar};
    }

    [[nodiscard]] constexpr double dot(const Vector3& other) const noexcept {
        return x_ * other.x_ + y_ * other.y_ + z_ * other.z_;
    }

    [[nodiscard]] constexpr Vector3 cross(const Vector3& other) const noexcept {
        return {
            y_ * other.z_ - z_ * other.y_,
            z_ * other.x_ - x_ * other.z_,
            x_ * other.y_ - y_ * other.x_,
        };
    }

    [[nodiscard]] constexpr double squared_norm() const noexcept { return dot(*this); }

    [[nodiscard]] double norm() const noexcept { return std::hypot(x_, y_, z_); }

    [[nodiscard]] Vector3 normalized() const {
        const double magnitude = norm();
        if (magnitude == 0.0) {
            throw std::domain_error("Cannot normalize a zero Vector3");
        }
        return *this / magnitude;
    }

private:
    double x_{0.0};
    double y_{0.0};
    double z_{0.0};
};

[[nodiscard]] constexpr Vector3 operator*(double scalar, const Vector3& vector) noexcept {
    return vector * scalar;
}

[[nodiscard]] inline bool approximately_equal(
    double lhs,
    double rhs,
    double absolute_tolerance = 1.0e-12,
    double relative_tolerance = 1.0e-12) noexcept {
    if (absolute_tolerance < 0.0 || relative_tolerance < 0.0) {
        return false;
    }

    const double scale = std::max(std::abs(lhs), std::abs(rhs));
    return std::abs(lhs - rhs) <= absolute_tolerance + relative_tolerance * scale;
}

[[nodiscard]] inline bool approximately_equal(
    const Vector3& lhs,
    const Vector3& rhs,
    double absolute_tolerance = 1.0e-12,
    double relative_tolerance = 1.0e-12) noexcept {
    return approximately_equal(lhs.x(), rhs.x(), absolute_tolerance, relative_tolerance)
        && approximately_equal(lhs.y(), rhs.y(), absolute_tolerance, relative_tolerance)
        && approximately_equal(lhs.z(), rhs.z(), absolute_tolerance, relative_tolerance);
}

[[nodiscard]] inline bool is_finite(const Vector3& vector) noexcept {
    return std::isfinite(vector.x()) && std::isfinite(vector.y()) && std::isfinite(vector.z());
}

}  // namespace astradock::math
