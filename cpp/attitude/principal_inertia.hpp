#pragma once

#include "math/matrix3.hpp"

#include <cmath>
#include <stdexcept>

namespace astradock::attitude {

// Represents the principal moments of inertia of a rigid spacecraft body
// expressed in its principal-axis body-fixed frame.
//
// In this frame, products of inertia (Ixy, Ixz, Iyz) are identically zero,
// simplifying the inertia tensor to a positive-definite diagonal matrix:
//   I = diag(Ixx, Iyy, Izz)
//
// Units: kg * m^2 (kilogram square-metres)
struct PrincipalInertia {
    double Ixx_kg_m2{1.0};
    double Iyy_kg_m2{1.0};
    double Izz_kg_m2{1.0};

    // Default constructor creates a normalized identity inertia (1, 1, 1) kg*m^2.
    constexpr PrincipalInertia() noexcept = default;

    // Parameterized constructor with strict physical validation.
    PrincipalInertia(double i_xx, double i_yy, double i_zz) {
        if (!std::isfinite(i_xx) || !std::isfinite(i_yy) || !std::isfinite(i_zz)) {
            throw std::domain_error("Principal moments of inertia must be finite");
        }
        if (i_xx <= 0.0 || i_yy <= 0.0 || i_zz <= 0.0) {
            throw std::domain_error("Principal moments of inertia must be strictly positive");
        }
        Ixx_kg_m2 = i_xx;
        Iyy_kg_m2 = i_yy;
        Izz_kg_m2 = i_zz;
    }

    // Component accessors
    [[nodiscard]] constexpr double Ixx() const noexcept { return Ixx_kg_m2; }
    [[nodiscard]] constexpr double Iyy() const noexcept { return Iyy_kg_m2; }
    [[nodiscard]] constexpr double Izz() const noexcept { return Izz_kg_m2; }

    // Computes inverse principal moments (1/Ixx, 1/Iyy, 1/Izz).
    [[nodiscard]] PrincipalInertia inverse() const {
        return PrincipalInertia(1.0 / Ixx_kg_m2, 1.0 / Iyy_kg_m2, 1.0 / Izz_kg_m2);
    }

    // Converts the principal moments into a 3x3 diagonal Direction Cosine / Inertia Matrix.
    [[nodiscard]] constexpr math::Matrix3 to_diagonal_matrix() const noexcept {
        return math::Matrix3(
            Ixx_kg_m2, 0.0,       0.0,
            0.0,       Iyy_kg_m2, 0.0,
            0.0,       0.0,       Izz_kg_m2
        );
    }
};

[[nodiscard]] inline bool is_finite(const PrincipalInertia& inertia) noexcept {
    return std::isfinite(inertia.Ixx_kg_m2)
        && std::isfinite(inertia.Iyy_kg_m2)
        && std::isfinite(inertia.Izz_kg_m2);
}

[[nodiscard]] inline bool is_valid(const PrincipalInertia& inertia) noexcept {
    return is_finite(inertia)
        && (inertia.Ixx_kg_m2 > 0.0)
        && (inertia.Iyy_kg_m2 > 0.0)
        && (inertia.Izz_kg_m2 > 0.0);
}

}  // namespace astradock::attitude
