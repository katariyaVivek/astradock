#pragma once

// AstraDock M22 — Synthetic camera model & projective geometry.
//
// Physical problem:
//   Before any pose estimation (classical M22, learned M23), the simulator
//   needs a camera that turns 3-D target geometry into 2-D pixel observations
//   through documented optics — no OpenCV, no renderer, just the pinhole
//   equations with explicit frames. Docking fiducials (known 3-D points on the
//   target body) project into the chaser camera; the pixel set IS the
//   measurement that M22-PnP and M23 models consume. Failure modes (missed /
//   false / occluded features, noise, blur-as-noise) are modeled at the pixel
//   level, deterministically.
//
// Frames (the heart of this header — four of them):
//   WORLD (target body): fiducial points p_w, meters, target-body coordinates.
//   CAMERA: origin at optical center, +Z along the boresight OUT of the lens
//     toward the scene, +X right, +Y down in the image (standard computer
//     vision convention). Extrinsics: p_c = R_cw * (p_w - t_w) with R_cw the
//     world-to-camera rotation and t_w the camera center in world coords.
//   NORMALIZED IMAGE: (x_n, y_n) = (X_c / Z_c, Y_c / Z_c), defined only for
//     Z_c > 0 (behind-camera is a hard rejection, not a projection).
//   PIXEL: u = fx * x_n + cx, v = fy * y_n + cy (+ optional single-parameter
//     radial distortion about the principal point).
//
// Units: meters (3-D), pixels (2-D), radians (angles), focal length in PIXELS
//   (fx, fy) so depth-to-pixel scaling needs no sensor-size constants.
//
// Governing equations:
//   p_c = R_cw (p_w - t_w);  visible iff Z_c > near_clip (default 0.1 m).
//   x_n = X_c/Z_c, y_n = Y_c/Z_c; r^2 = x_n^2 + y_n^2.
//   Radial distortion: s = 1 + k1 r^2; (x_d, y_d) = s (x_n, y_n).
//   u = fx x_d + cx, v = fy y_d + cy. In-bounds iff 0 <= u < W, 0 <= v < H.
//   Pixel noise: additive i.i.d. Gaussian per axis (M12 DeterministicRng
//     discipline for seeded runs — the RNG lives in the harness, this header
//     takes noise as an input vector to stay pure).
//
// Assumptions: pinhole (no skew, square-or-not pixels via fx != fy allowed);
// single radial coefficient (barrel/pincushion to first order); no tangential
// distortion, no rolling shutter, no motion blur beyond noise inflation;
// fiducials are ideal points (no extent, no perspective-foreshortening of the
// marker itself — M23 may add appearance).

#include "math/matrix3.hpp"
#include "math/vector3.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace astradock::vision {

namespace detail {

inline void require_finite_camera(double value, const char* message) {
    if (!std::isfinite(value)) {
        throw std::domain_error(message);
    }
}

}  // namespace detail

// Camera intrinsics: focal lengths (pixels), principal point (pixels),
// image size, single radial distortion coefficient, near clip plane.
struct CameraIntrinsics {
    double fx_px{800.0};
    double fy_px{800.0};
    double cx_px{320.0};
    double cy_px{240.0};
    int width_px{640};
    int height_px{480};
    double k1_radial{0.0};
    double near_clip_m{0.1};

    void validate() const {
        if (!std::isfinite(fx_px) || fx_px <= 0.0 || !std::isfinite(fy_px) || fy_px <= 0.0) {
            throw std::domain_error("Camera focal lengths must be finite and positive");
        }
        if (!std::isfinite(cx_px) || !std::isfinite(cy_px)) {
            throw std::domain_error("Camera principal point must be finite");
        }
        if (width_px <= 0 || height_px <= 0) {
            throw std::domain_error("Camera image size must be positive");
        }
        if (!std::isfinite(k1_radial)) {
            throw std::domain_error("Camera distortion coefficient must be finite");
        }
        if (!std::isfinite(near_clip_m) || near_clip_m <= 0.0) {
            throw std::domain_error("Camera near clip must be finite and positive");
        }
    }
};

// Camera extrinsics: world-to-camera rotation + camera center in world coords.
struct CameraExtrinsics {
    math::Matrix3 rotation_world_to_camera{math::Matrix3::identity()};
    math::Vector3 center_world_m{};

    void validate() const {
        if (!math::is_finite(rotation_world_to_camera) || !math::is_finite(center_world_m)) {
            throw std::domain_error("Camera extrinsics must be finite");
        }
        if (!rotation_world_to_camera.is_orthonormal(1.0e-9, 1.0e-9)) {
            throw std::domain_error("World-to-camera rotation must be orthonormal");
        }
    }
};

// One projected observation: pixel + visibility flags. A point can be
// invisible for exactly one reason — the flags say which (never silent).
struct PixelObservation {
    double u_px{0.0};
    double v_px{0.0};
    bool visible{false};
    bool behind_camera{false};
    bool out_of_bounds{false};
};

// 3-D point in world -> camera coordinates.
[[nodiscard]] inline math::Vector3 world_to_camera(
    const math::Vector3& point_world_m,
    const CameraExtrinsics& extrinsics) {
    extrinsics.validate();
    if (!math::is_finite(point_world_m)) {
        throw std::domain_error("World point must be finite");
    }
    return extrinsics.rotation_world_to_camera * (point_world_m - extrinsics.center_world_m);
}

// Pinhole projection of a camera-frame point to pixels (no noise).
[[nodiscard]] inline PixelObservation project_camera_point(
    const math::Vector3& point_camera_m,
    const CameraIntrinsics& intrinsics) {
    intrinsics.validate();
    if (!math::is_finite(point_camera_m)) {
        throw std::domain_error("Camera-frame point must be finite");
    }
    PixelObservation obs;
    if (point_camera_m.z() <= intrinsics.near_clip_m) {
        obs.behind_camera = true;
        return obs;
    }
    const double x_n = point_camera_m.x() / point_camera_m.z();
    const double y_n = point_camera_m.y() / point_camera_m.z();
    const double r2 = x_n * x_n + y_n * y_n;
    const double s = 1.0 + intrinsics.k1_radial * r2;
    const double u = intrinsics.fx_px * (s * x_n) + intrinsics.cx_px;
    const double v = intrinsics.fy_px * (s * y_n) + intrinsics.cy_px;
    obs.u_px = u;
    obs.v_px = v;
    if (u < 0.0 || u >= static_cast<double>(intrinsics.width_px) || v < 0.0
        || v >= static_cast<double>(intrinsics.height_px)) {
        obs.out_of_bounds = true;
        return obs;
    }
    obs.visible = true;
    return obs;
}

// Full pipeline: world -> camera -> pixels.
[[nodiscard]] inline PixelObservation project_world_point(
    const math::Vector3& point_world_m,
    const CameraExtrinsics& extrinsics,
    const CameraIntrinsics& intrinsics) {
    return project_camera_point(world_to_camera(point_world_m, extrinsics), intrinsics);
}

// Back-projection: pixel -> unit bearing vector in camera coordinates
// (undistorted first via one Newton refinement when k1 != 0, then normalized
// [x_n, y_n, 1]). Bearing (not range): the third component is exactly 1
// before normalization... returns the NORMALIZED direction.
[[nodiscard]] inline math::Vector3 pixel_to_bearing(
    double u_px,
    double v_px,
    const CameraIntrinsics& intrinsics) {
    intrinsics.validate();
    detail::require_finite_camera(u_px, "Pixel u must be finite");
    detail::require_finite_camera(v_px, "Pixel v must be finite");
    double x_d = (u_px - intrinsics.cx_px) / intrinsics.fx_px;
    double y_d = (v_px - intrinsics.cy_px) / intrinsics.fy_px;
    double x_n = x_d;
    double y_n = y_d;
    if (intrinsics.k1_radial != 0.0) {
        // Invert s = 1 + k1 r^2 iteratively (3 fixed-point steps suffice for
        // |k1| r^2 << 1 educational lenses; exactness is tested at k1 = 0).
        for (int i = 0; i < 8; ++i) {
            const double r2 = x_n * x_n + y_n * y_n;
            const double s = 1.0 + intrinsics.k1_radial * r2;
            x_n = x_d / s;
            y_n = y_d / s;
        }
    }
    const math::Vector3 bearing{x_n, y_n, 1.0};
    return bearing.normalized();
}

// Reprojection error in pixels between an observed pixel and a re-projected
// world point (the PnP cost function, M22 pose header consumes this).
[[nodiscard]] inline double reprojection_error_px(
    double observed_u_px,
    double observed_v_px,
    const math::Vector3& point_world_m,
    const CameraExtrinsics& extrinsics,
    const CameraIntrinsics& intrinsics) {
    const PixelObservation predicted = project_world_point(point_world_m, extrinsics, intrinsics);
    if (!predicted.visible) {
        throw std::domain_error("Reprojection target is not visible; error is undefined");
    }
    const double du = observed_u_px - predicted.u_px;
    const double dv = observed_v_px - predicted.v_px;
    return std::sqrt(du * du + dv * dv);
}

// Docking fiducial set: known 3-D points on the target body (world frame).
struct FiducialSet {
    std::vector<math::Vector3> points_world_m{};

    void validate() const {
        if (points_world_m.size() < 4) {
            throw std::domain_error("Fiducial set needs at least 4 non-coplanar points for PnP");
        }
        for (const auto& point : points_world_m) {
            if (!math::is_finite(point)) {
                throw std::domain_error("Fiducial points must be finite");
            }
        }
    }
};

// Default docking fiducial plate: 4 corners of a 1 m square at z = 0 plus a
// 0.25 m standoff center (non-coplanar 5th point breaks the planar ambiguity).
[[nodiscard]] inline FiducialSet default_fiducial_plate() {
    FiducialSet set;
    set.points_world_m = {
        math::Vector3{-0.5, -0.5, 0.0},
        math::Vector3{0.5, -0.5, 0.0},
        math::Vector3{0.5, 0.5, 0.0},
        math::Vector3{-0.5, 0.5, 0.0},
        math::Vector3{0.0, 0.0, 0.25},
    };
    return set;
}

}  // namespace astradock::vision
