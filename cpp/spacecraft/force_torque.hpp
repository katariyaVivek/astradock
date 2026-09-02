#pragma once

#include "math/vector3.hpp"

namespace astradock::spacecraft {

// Generalized external force and torque input applied to the spacecraft.
//
// Frame Conventions:
//   - force_eci_N: Resultant external force vector expressed in the Earth-Centered
//                  Inertial (ECI) coordinate frame (N).
//   - torque_body_Nm: Resultant external torque vector about the spacecraft center of
//                     mass, expressed in the spacecraft principal BODY coordinate frame (N*m).
//
// Why Different Frames:
//   Translational orbital equations are naturally formulated in the inertial ECI frame.
//   Rotational Euler rigid-body equations are formulated in the body-fixed principal axes
//   where the inertia tensor I is constant and diagonal.
//
// Default State:
//   Zero external force and zero external torque (unforced two-body orbit + torque-free spin).
struct ForceTorqueInput {
    math::Vector3 force_eci_N{};
    math::Vector3 torque_body_Nm{};

    constexpr ForceTorqueInput() noexcept = default;

    constexpr ForceTorqueInput(
        const math::Vector3& force_eci,
        const math::Vector3& torque_body) noexcept
        : force_eci_N(force_eci),
          torque_body_Nm(torque_body) {}
};

[[nodiscard]] inline bool is_finite(const ForceTorqueInput& input) noexcept {
    return math::is_finite(input.force_eci_N) && math::is_finite(input.torque_body_Nm);
}

}  // namespace astradock::spacecraft
