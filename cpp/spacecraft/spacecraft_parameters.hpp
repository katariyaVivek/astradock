#pragma once

#include "attitude/principal_inertia.hpp"

#include <cmath>

namespace astradock::spacecraft {

// Physical parameters for 6-DOF rigid-spacecraft simulation.
//
// In accordance with the point-mass gravity model, translational gravitational
// acceleration a = -mu * r / |r|^3 is independent of spacecraft mass.
// Therefore, spacecraft mass is not required for idealized two-body central motion.
// Spacecraft rotational dynamics depends on the principal moments of inertia.
struct SpacecraftParameters {
    double gravitational_parameter_m3_per_s2{3.986004418e14}; // WGS 84 default (m^3/s^2)
    attitude::PrincipalInertia inertia{10.0, 10.0, 10.0};       // Principal inertia moments (kg*m^2)

    constexpr SpacecraftParameters() noexcept = default;

    SpacecraftParameters(
        double mu_m3_per_s2,
        const attitude::PrincipalInertia& principal_inertia)
        : gravitational_parameter_m3_per_s2(mu_m3_per_s2),
          inertia(principal_inertia) {
        if (!std::isfinite(mu_m3_per_s2) || mu_m3_per_s2 <= 0.0) {
            throw std::domain_error("Gravitational parameter must be finite and strictly positive");
        }
        if (!attitude::is_valid(principal_inertia)) {
            throw std::domain_error("Principal inertia moments must be finite and strictly positive");
        }
    }
};

[[nodiscard]] inline bool is_valid(const SpacecraftParameters& params) noexcept {
    return std::isfinite(params.gravitational_parameter_m3_per_s2)
        && params.gravitational_parameter_m3_per_s2 > 0.0
        && attitude::is_valid(params.inertia);
}

}  // namespace astradock::spacecraft
