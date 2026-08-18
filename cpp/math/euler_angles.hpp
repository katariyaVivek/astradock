#pragma once

#include "math/angle.hpp"
#include "math/matrix3.hpp"
#include "math/quaternion.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace astradock::math {

// Euler Angles representation in the standard aerospace 3-2-1 (ZYX) rotation sequence:
//   yaw_rad:   rotation about +Z axis (psi) in [-pi, pi] or [0, 2pi)
//   pitch_rad: rotation about +Y' axis (theta) in [-pi/2, pi/2]
//   roll_rad:  rotation about +X'' axis (phi) in [-pi, pi] or [0, 2pi)
//
// Active Rotation Matrix:
//   C = R_z(yaw) * R_y(pitch) * R_x(roll)
struct EulerAngles {
    double yaw_rad;    // psi (about Z)
    double pitch_rad;  // theta (about Y)
    double roll_rad;   // phi (about X)
};

// Converts ZYX Euler angles into a Direction Cosine Matrix (DCM).
[[nodiscard]] inline Matrix3 euler_to_rotation_matrix(const EulerAngles& euler) {
    if (!std::isfinite(euler.yaw_rad) || !std::isfinite(euler.pitch_rad) || !std::isfinite(euler.roll_rad)) {
        throw std::domain_error("Euler angles must contain only finite values");
    }

    const double cy = std::cos(euler.yaw_rad);
    const double sy = std::sin(euler.yaw_rad);
    const double cp = std::cos(euler.pitch_rad);
    const double sp = std::sin(euler.pitch_rad);
    const double cr = std::cos(euler.roll_rad);
    const double sr = std::sin(euler.roll_rad);

    return Matrix3(
        cy * cp, cy * sp * sr - sy * cr, cy * sp * cr + sy * sr,
        sy * cp, sy * sp * sr + cy * cr, sy * sp * cr - cy * sr,
        -sp,     cp * sr,                cp * cr
    );
}

// Converts a Direction Cosine Matrix (DCM) into ZYX Euler angles.
// Detects gimbal lock at pitch ≈ ±90 deg (|pitch| ≈ pi/2) where yaw and roll become coupled.
[[nodiscard]] inline EulerAngles rotation_matrix_to_euler(const Matrix3& matrix, double gimbal_lock_tol = 1.0e-8) {
    if (!is_finite(matrix)) {
        throw std::domain_error("Cannot convert non-finite rotation matrix to Euler angles");
    }

    const double sp = -matrix(2, 0);
    double pitch = 0.0;
    double yaw = 0.0;
    double roll = 0.0;

    if (sp >= 1.0 - gimbal_lock_tol) {
        // Gimbal lock: pitch ≈ +90 deg (theta = +pi/2)
        pitch = 0.5 * constants::pi;
        yaw = 0.0; // Convention
        roll = std::atan2(matrix(0, 1), matrix(0, 2));
    } else if (sp <= -1.0 + gimbal_lock_tol) {
        // Gimbal lock: pitch ≈ -90 deg (theta = -pi/2)
        pitch = -0.5 * constants::pi;
        yaw = 0.0; // Convention
        roll = std::atan2(-matrix(0, 1), matrix(0, 2));
    } else {
        pitch = std::asin(std::clamp(sp, -1.0, 1.0));
        yaw = std::atan2(matrix(1, 0), matrix(0, 0));
        roll = std::atan2(matrix(2, 1), matrix(2, 2));
    }

    return {yaw, pitch, roll};
}

// Converts ZYX Euler angles to a unit quaternion.
[[nodiscard]] inline Quaternion euler_to_quaternion(const EulerAngles& euler) {
    const Matrix3 mat = euler_to_rotation_matrix(euler);
    return Quaternion::from_rotation_matrix(mat);
}

// Converts a unit quaternion to ZYX Euler angles.
[[nodiscard]] inline EulerAngles quaternion_to_euler(const Quaternion& q) {
    const Matrix3 mat = q.to_rotation_matrix();
    return rotation_matrix_to_euler(mat);
}

}  // namespace astradock::math
