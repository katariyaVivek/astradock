#pragma once

#include "frames/frame_basis.hpp"
#include "math/matrix3.hpp"
#include "math/vector3.hpp"

#include <cmath>
#include <stdexcept>

namespace astradock::frames {

namespace detail {

inline void require_finite_vector(const math::Vector3& v, const char* message) {
    if (!std::isfinite(v.x()) || !std::isfinite(v.y()) || !std::isfinite(v.z())) {
        throw std::domain_error(message);
    }
}

}  // namespace detail

// Computes the orthonormal LVLH (Local Vertical Local Horizontal) basis triad
// expressed in Earth-Centered Inertial (ECI) coordinates.
//
// Convention:
//   r_hat = r_eci / ||r_eci||          (radial direction, pointing from Earth center to satellite)
//   h     = r_eci x v_eci              (specific angular momentum vector)
//   h_hat = h / ||h||                  (orbit-normal direction, perpendicular to orbital plane)
//   t_hat = h_hat x r_hat              (along-track direction, in direction of orbital motion)
//
// Basis triad:
//   x = e_r = r_hat  (radial axis)
//   y = e_t = t_hat  (along-track axis)
//   z = e_h = h_hat  (orbit-normal axis)
//
// Handedness:
//   e_r x e_t = r_hat x (h_hat x r_hat) = h_hat = e_h (strictly right-handed).
//
// Rejects:
//   - Non-finite components in position or velocity.
//   - Zero position norm (orbital position undefined at central-body origin).
//   - Zero angular momentum norm (position and velocity are collinear / zero, so orbital plane is undefined).
[[nodiscard]] inline FrameBasis compute_lvlh_basis(
    const math::Vector3& position_eci,
    const math::Vector3& velocity_eci) {
    detail::require_finite_vector(
        position_eci,
        "Cannot compute LVLH basis: orbital position contains non-finite values");
    detail::require_finite_vector(
        velocity_eci,
        "Cannot compute LVLH basis: orbital velocity contains non-finite values");

    const double r_norm = position_eci.norm();
    if (r_norm == 0.0) {
        throw std::domain_error(
            "Cannot compute LVLH basis: orbital position magnitude is zero");
    }

    const math::Vector3 angular_momentum = position_eci.cross(velocity_eci);
    const double h_norm = angular_momentum.norm();
    if (h_norm == 0.0) {
        throw std::domain_error(
            "Cannot compute LVLH basis: angular momentum is zero (collinear position and velocity)");
    }

    const math::Vector3 e_r = position_eci / r_norm;
    const math::Vector3 e_h = angular_momentum / h_norm;
    const math::Vector3 e_t = e_h.cross(e_r);

    return FrameBasis{e_r, e_t, e_h};
}

// Direction Cosine Matrix (DCM) transforming vector components from ECI to LVLH:
//   v_lvlh = C_LVLH_ECI * v_eci
//
// Because v_lvlh components are the dot products with the LVLH basis vectors
// expressed in ECI:
//   v_lvlh.x = e_r · v_eci
//   v_lvlh.y = e_t · v_eci
//   v_lvlh.z = e_h · v_eci
//
// The rows of C_LVLH_ECI are the basis vectors e_r, e_t, e_h expressed in ECI:
//   C_LVLH_ECI = [ e_r^T ]
//                [ e_t^T ]
//                [ e_h^T ]
[[nodiscard]] inline math::Matrix3 dcm_lvlh_from_eci(const FrameBasis& lvlh_basis) {
    return math::Matrix3::from_rows(lvlh_basis.x, lvlh_basis.y, lvlh_basis.z);
}

[[nodiscard]] inline math::Matrix3 dcm_lvlh_from_eci(
    const math::Vector3& position_eci,
    const math::Vector3& velocity_eci) {
    return dcm_lvlh_from_eci(compute_lvlh_basis(position_eci, velocity_eci));
}

// Direction Cosine Matrix (DCM) transforming vector components from LVLH to ECI:
//   v_eci = C_ECI_LVLH * v_lvlh
//
// Since C is orthonormal, C_ECI_LVLH = (C_LVLH_ECI)^T.
// The columns of C_ECI_LVLH are the basis vectors e_r, e_t, e_h expressed in ECI:
//   C_ECI_LVLH = [ e_r, e_t, e_h ]
[[nodiscard]] inline math::Matrix3 dcm_eci_from_lvlh(const FrameBasis& lvlh_basis) {
    return math::Matrix3::from_columns(lvlh_basis.x, lvlh_basis.y, lvlh_basis.z);
}

[[nodiscard]] inline math::Matrix3 dcm_eci_from_lvlh(
    const math::Vector3& position_eci,
    const math::Vector3& velocity_eci) {
    return dcm_eci_from_lvlh(compute_lvlh_basis(position_eci, velocity_eci));
}

// Transforms coordinates of a physical vector expressed in ECI into coordinates
// expressed in LVLH:
//   v_lvlh = dcm_lvlh_from_eci * v_eci
//
// NOTE: This performs a geometric vector coordinate transformation.
// For velocity, this projects the inertial velocity onto the instantaneous LVLH axes.
// It is NOT the apparent relative velocity observed from the rotating LVLH frame,
// which requires accounting for frame angular velocity: v_rel = dr/dt - omega x r.
[[nodiscard]] inline math::Vector3 transform_eci_to_lvlh(
    const math::Vector3& vector_eci,
    const math::Matrix3& dcm_lvlh_from_eci) noexcept {
    return dcm_lvlh_from_eci * vector_eci;
}

// Transforms coordinates of a physical vector expressed in LVLH into coordinates
// expressed in ECI:
//   v_eci = (dcm_lvlh_from_eci)^T * v_lvlh = dcm_eci_from_lvlh * v_lvlh
[[nodiscard]] inline math::Vector3 transform_lvlh_to_eci(
    const math::Vector3& vector_lvlh,
    const math::Matrix3& dcm_lvlh_from_eci) noexcept {
    return dcm_lvlh_from_eci.transpose() * vector_lvlh;
}

}  // namespace astradock::frames
