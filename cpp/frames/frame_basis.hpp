#pragma once

#include "math/vector3.hpp"

#include <cmath>

namespace astradock::frames {

// An orthonormal basis triad (x, y, z) defining a 3D coordinate frame.
// Each basis vector is a unit vector expressed in a parent/reference frame.
//
// For example, in an LVLH frame resolved in ECI coordinates:
//   x = radial unit vector (e_r)
//   y = along-track unit vector (e_t)
//   z = orbit-normal unit vector (e_h)
// where each vector has length 1.0 and they form a right-handed triad: x x y = z.
struct FrameBasis {
    math::Vector3 x;
    math::Vector3 y;
    math::Vector3 z;
};

// Checks whether all components in the basis triad are finite real numbers.
[[nodiscard]] inline bool is_finite(const FrameBasis& basis) noexcept {
    return std::isfinite(basis.x.x()) && std::isfinite(basis.x.y()) && std::isfinite(basis.x.z())
        && std::isfinite(basis.y.x()) && std::isfinite(basis.y.y()) && std::isfinite(basis.y.z())
        && std::isfinite(basis.z.x()) && std::isfinite(basis.z.y()) && std::isfinite(basis.z.z());
}

// Validates that the frame basis satisfies orthonormality and right-handedness:
//   1. ||x|| ≈ 1, ||y|| ≈ 1, ||z|| ≈ 1
//   2. x · y ≈ 0, y · z ≈ 0, z · x ≈ 0
//   3. (x × y) · z ≈ +1 (right-handed orientation)
[[nodiscard]] inline bool is_orthonormal(
    const FrameBasis& basis,
    double absolute_tolerance = 1.0e-12,
    double relative_tolerance = 1.0e-12) noexcept {
    if (!is_finite(basis)) {
        return false;
    }

    // 1. Unit norms
    if (!math::approximately_equal(basis.x.norm(), 1.0, absolute_tolerance, relative_tolerance)
        || !math::approximately_equal(basis.y.norm(), 1.0, absolute_tolerance, relative_tolerance)
        || !math::approximately_equal(basis.z.norm(), 1.0, absolute_tolerance, relative_tolerance)) {
        return false;
    }

    // 2. Pairwise orthogonality
    if (!math::approximately_equal(basis.x.dot(basis.y), 0.0, absolute_tolerance, relative_tolerance)
        || !math::approximately_equal(basis.y.dot(basis.z), 0.0, absolute_tolerance, relative_tolerance)
        || !math::approximately_equal(basis.z.dot(basis.x), 0.0, absolute_tolerance, relative_tolerance)) {
        return false;
    }

    // 3. Right-handed orientation (cross product x × y should equal z, scalar triple product = 1)
    const math::Vector3 cross_xy = basis.x.cross(basis.y);
    if (!math::approximately_equal(cross_xy, basis.z, absolute_tolerance, relative_tolerance)) {
        return false;
    }

    return math::approximately_equal(cross_xy.dot(basis.z), 1.0, absolute_tolerance, relative_tolerance);
}

}  // namespace astradock::frames
