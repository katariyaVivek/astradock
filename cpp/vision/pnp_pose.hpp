#pragma once

// AstraDock M22 — Geometric relative-pose estimation from fiducial pixels.
//
// Physical problem:
//   Given N 2-D pixel observations of KNOWN 3-D fiducial points (M22 camera),
//   recover the camera pose (extrinsics): rotation R_cw + center t_w. This is
//   the Perspective-n-Point problem, solved here CLASSICALLY (Gauss-Newton on
//   the reprojection-error least-squares cost) — the baseline every learned
//   method (M23) must beat, and the honest fallback when learning is
//   unavailable. Deterministic, no RANSAC library: robustness comes from an
//   explicit inlier gate + single re-fit, with failure REPORTED (pose_valid =
//   false) instead of a confident wrong answer.
//
// Frames/units: inherits the camera header (world meters, pixels). Pose
//   output is CameraExtrinsics (world-to-camera). Translation reported both
//   as center-in-world and range (center norm) for docking relevance.
//
// Governing method (Gauss-Newton on SO(3) x R^3):
//   cost = sum_i ||obs_i - project(p_i; R, t)||^2 over inlier pixels.
//   Parameterization: axis-angle delta (3) + translation delta (3), updated as
//     R <- exp([dphi]x) R,  t <- t + R^T ... (translation in world coords).
//   Jacobians: analytic pinhole chain rule; rotation part via the cross-product
//     linearization d(Rp)/dphi = -[Rp]x (world-frame perturbation of the
//     rotated point). Damped normal equations (Levenberg-Marquardt-lite with a
//     fixed small lambda) for conditioning; 6x6 SPD solve via Cholesky.
//   Init: true-extrinsics-perturbed seed supplied by the CALLER (script-side),
//     because global initialization (DLT/EPnP) is out of scope — documented.
//     In operations the seed is the filter prediction / acquisition estimate.
//   Outlier gate: per-point reprojection error vs pixel gate; one re-fit on
//     inliers; pose rejected when inliers < 4 or RMS above tolerance.
//
// Assumptions: known correspondences (association is given); distortion known;
// at least 4 non-coplanar inliers; seed within the basin (tested: 5 deg / 1 m).
// Failure modes are FIRST-CLASS outputs: divergent cost, singular normal
// equations, inlier starvation all yield pose_valid = false with a reason.

#include "math/linalg.hpp"
#include "math/matrix3.hpp"
#include "math/vector3.hpp"
#include "vision/camera_model.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace astradock::vision {

namespace detail {

inline math::Matrix3 axis_angle_to_matrix(const math::Vector3& phi) {
    const double angle = phi.norm();
    if (angle < 1.0e-12) {
        return math::Matrix3::identity();
    }
    const math::Vector3 axis = phi / angle;
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    const double x = axis.x();
    const double y = axis.y();
    const double z = axis.z();
    // Rodrigues rotation formula.
    return math::Matrix3(
        c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s,
        y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s,
        z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c));
}

inline math::Matrix3 skew(const math::Vector3& v) {
    return math::Matrix3(
        0.0, -v.z(), v.y(),
        v.z(), 0.0, -v.x(),
        -v.y(), v.x(), 0.0);
}

}  // namespace detail

// Pose estimation result: extrinsics + diagnostics. pose_valid == false means
// "no trustworthy pose" — downstream MUST check it (docking never consumes an
// invalid pose; the demo asserts this path explicitly).
struct PoseEstimate {
    CameraExtrinsics extrinsics{};
    bool pose_valid{false};
    std::string failure_reason{"uninitialized"};
    double rms_reprojection_px{0.0};
    std::size_t inliers{0};
    std::size_t iterations{0};
    double range_m{0.0};
};

struct PnpOptions {
    int max_iterations{50};
    double inlier_gate_px{3.0};
    double rms_accept_px{2.0};
    double damping{1.0e-3};  // relative to diag(J^T J): scale-aware LM damping
    double step_tolerance{1.0e-9};
    double huber_delta_px{4.0};  // Huber knee for pass-1 robustification
    bool robust_pass{true};      // Huber-weight pass 1, hard-gate + refit after
};

// Residual vector + analytic 2Nx6 Jacobian for the current pose. Rows: pixel
// (u, v) residuals; columns: [dphi_world(3), dt_world(3)]. Derivation: with
// q = R(p - t), dq/dphi = -[q]x (world-frame rotation perturbation) and
// dq/dt = -R. Pixel chain rule through the pinhole + radial distortion.
struct PnpLinearization {
    std::vector<double> residuals;  // length 2N
    math::Matrix<6, 6> normal{};
    std::array<double, 6> gradient{};
};

[[nodiscard]] inline PnpLinearization pnp_linearize(
    const std::vector<math::Vector3>& points_world_m,
    const std::vector<math::Vector3>& observed_px,  // (u, v, 0) carriers
    const CameraExtrinsics& extrinsics,
    const CameraIntrinsics& intrinsics,
    bool use_huber = false,
    double huber_delta_px = 4.0) {
    extrinsics.validate();
    intrinsics.validate();
    if (points_world_m.size() != observed_px.size() || points_world_m.empty()) {
        throw std::domain_error("PnP needs matching non-empty point/observation sets");
    }
    PnpLinearization lin;
    lin.residuals.reserve(2 * points_world_m.size());
    lin.normal = math::Matrix<6, 6>::zero();
    for (int k = 0; k < 6; ++k) {
        lin.gradient[k] = 0.0;
    }
    const math::Matrix3& R = extrinsics.rotation_world_to_camera;
    for (std::size_t i = 0; i < points_world_m.size(); ++i) {
        const math::Vector3 q = R * (points_world_m[i] - extrinsics.center_world_m);
        const double z = q.z();
        if (z <= intrinsics.near_clip_m) {
            throw std::domain_error("PnP linearization point behind near clip");
        }
        const double x_n = q.x() / z;
        const double y_n = q.y() / z;
        const double r2 = x_n * x_n + y_n * y_n;
        const double s = 1.0 + intrinsics.k1_radial * r2;
        // Distorted normalized coords and their derivatives w.r.t. (x_n, y_n).
        const double xd = s * x_n;
        const double yd = s * y_n;
        const double ds_dx = 2.0 * intrinsics.k1_radial * x_n;
        const double ds_dy = 2.0 * intrinsics.k1_radial * y_n;
        // d(pixel)/d(x_n, y_n): 2x2 chain.
        const double du_dxn = intrinsics.fx_px * (s + x_n * ds_dx);
        const double du_dyn = intrinsics.fx_px * (x_n * ds_dy);
        const double dv_dxn = intrinsics.fy_px * (y_n * ds_dx);
        const double dv_dyn = intrinsics.fy_px * (s + y_n * ds_dy);
        // d(x_n, y_n)/dq: 2x3 perspective division.
        const double z2 = z * z;
        // Rows of the 2x6 Jacobian for rotation (world perturbation) then translation.
        const math::Matrix3 neg_skew_q = detail::skew(q) * (-1.0);  // dq/dphi
        const math::Matrix3 neg_R = R * (-1.0);                     // dq/dt
        double J_phi[2][3];
        double J_t[2][3];
        for (int c = 0; c < 3; ++c) {
            const math::Vector3 dphi_col{neg_skew_q(0, c), neg_skew_q(1, c), neg_skew_q(2, c)};
            const double dxn = (dphi_col.x() * z - q.x() * dphi_col.z()) / z2;
            const double dyn = (dphi_col.y() * z - q.y() * dphi_col.z()) / z2;
            J_phi[0][c] = du_dxn * dxn + du_dyn * dyn;
            J_phi[1][c] = dv_dxn * dxn + dv_dyn * dyn;
            const math::Vector3 dt_col{neg_R(0, c), neg_R(1, c), neg_R(2, c)};
            const double txn = (dt_col.x() * z - q.x() * dt_col.z()) / z2;
            const double tyn = (dt_col.y() * z - q.y() * dt_col.z()) / z2;
            J_t[0][c] = du_dxn * txn + du_dyn * tyn;
            J_t[1][c] = dv_dxn * txn + dv_dyn * tyn;
        }
        const double pred_u = intrinsics.fx_px * xd + intrinsics.cx_px;
        const double pred_v = intrinsics.fy_px * yd + intrinsics.cy_px;
        const double res_u = observed_px[i].x() - pred_u;
        const double res_v = observed_px[i].y() - pred_v;
        // Huber weight from the CURRENT reprojection error of this point:
        // w = 1 inside the knee, delta/err outside (IRLS). Both rows of the
        // point share sqrt(w) so the normal equations see w * J'J, w * J'r.
        const double point_err = std::sqrt(res_u * res_u + res_v * res_v);
        double point_weight = 1.0;
        if (use_huber && point_err > huber_delta_px) {
            point_weight = huber_delta_px / point_err;
        }
        const double sqrt_w = std::sqrt(point_weight);
        lin.residuals.push_back(sqrt_w * res_u);
        lin.residuals.push_back(sqrt_w * res_v);
        const double* rows[2] = {nullptr, nullptr};
        (void)rows;
        double J[2][6];
        for (int c = 0; c < 3; ++c) {
            J[0][c] = J_phi[0][c];
            J[1][c] = J_phi[1][c];
            J[0][3 + c] = J_t[0][c];
            J[1][3 + c] = J_t[1][c];
        }
        const double res[2] = {sqrt_w * res_u, sqrt_w * res_v};
        double Jw[2][6];
        for (int r = 0; r < 2; ++r) {
            for (int a = 0; a < 6; ++a) {
                Jw[r][a] = sqrt_w * J[r][a];
            }
        }
        for (int r = 0; r < 2; ++r) {
            for (int a = 0; a < 6; ++a) {
                lin.gradient[a] += Jw[r][a] * res[r];
                for (int b = 0; b < 6; ++b) {
                    lin.normal(a, b) += Jw[r][a] * Jw[r][b];
                }
            }
        }
    }
    return lin;
}

// Gauss-Newton PnP with inlier gating + single re-fit. Correspondences are
// given (associated). Seed = caller-supplied initial extrinsics.
[[nodiscard]] inline PoseEstimate estimate_pose_pnp(
    const std::vector<math::Vector3>& points_world_m,
    const std::vector<math::Vector3>& observed_px,
    const CameraExtrinsics& seed_extrinsics,
    const CameraIntrinsics& intrinsics,
    const PnpOptions& options = PnpOptions{}) {
    if (points_world_m.size() != observed_px.size() || points_world_m.size() < 4) {
        throw std::domain_error("PnP needs >= 4 matched point/observation pairs");
    }
    seed_extrinsics.validate();
    intrinsics.validate();
    if (options.max_iterations <= 0 || options.inlier_gate_px <= 0.0 || options.rms_accept_px <= 0.0) {
        throw std::domain_error("PnP options must be positive");
    }

    PoseEstimate result;
    CameraExtrinsics current = seed_extrinsics;
    const auto solve_once = [&](const std::vector<std::size_t>& indices,
                                CameraExtrinsics start,
                                bool robust) -> std::pair<CameraExtrinsics, double> {
        CameraExtrinsics pose = start;
        double prev_cost = 0.0;
        // Adaptive Levenberg-Marquardt damping: shrink on success (approach
        // Gauss-Newton for quadratic finish), grow on failure. Fixed lambda
        // crawls near the optimum and stalls at millimeter accuracy.
        double lambda = options.damping;
        for (int iter = 0; iter < options.max_iterations; ++iter) {
            std::vector<math::Vector3> sub_pts;
            std::vector<math::Vector3> sub_obs;
            sub_pts.reserve(indices.size());
            sub_obs.reserve(indices.size());
            for (std::size_t k : indices) {
                sub_pts.push_back(points_world_m[k]);
                sub_obs.push_back(observed_px[k]);
            }
            PnpLinearization lin;
            try {
                lin = pnp_linearize(
                    sub_pts, sub_obs, pose, intrinsics, robust, options.huber_delta_px);
            } catch (const std::domain_error&) {
                return {pose, -1.0};
            }
            double cost = 0.0;
            for (double r : lin.residuals) {
                cost += r * r;
            }
            // Convergence: cost stalled AND the parameter step is tiny.
            // Cost-difference alone stalls prematurely under damping (slow
            // linear crawl with sizable steps still ahead); the step norm is
            // the honest certificate. Both must hold to quit iterating.
            bool cost_stalled = iter > 0
                && std::abs(prev_cost - cost) < options.step_tolerance * (1.0 + prev_cost);
            // Levenberg-Marquardt damping scaled to the normal-equation
            // diagonal (a fixed absolute lambda is meaningless next to J^T J
            // entries ~1e6 from 800 px focal lengths). Backtracking line
            // search halves the step until the cost decreases — standard
            // hardening against outlier-yanked Gauss-Newton overshoot.
            math::Matrix<6, 6> damped = lin.normal;
            for (int d = 0; d < 6; ++d) {
                damped(d, d) += lambda * lin.normal(d, d);
            }
            math::Matrix<6, 1> pos_grad;
            for (int k = 0; k < 6; ++k) {
                pos_grad(k, 0) = lin.gradient[k];
            }
            math::Matrix<6, 1> delta;
            try {
                // Gauss-Newton step for residuals r = observed - predicted:
                //   delta = +(J^T J)^-1 J^T r (descent; the negation would climb).
                delta = math::solve_spd(damped, pos_grad);
            } catch (const std::runtime_error&) {
                return {pose, -1.0};
            }
            const math::Vector3 dphi_full{delta(0, 0), delta(1, 0), delta(2, 0)};
            const math::Vector3 dt_full{delta(3, 0), delta(4, 0), delta(5, 0)};
            const double step_norm = std::sqrt(
                dphi_full.squared_norm() + dt_full.squared_norm());
            if (cost_stalled && step_norm < 1.0e-9) {
                prev_cost = cost;
                break;  // certified: flat cost AND negligible step
            }
            prev_cost = cost;
            bool improved = false;
            for (int halve = 0; halve < 10; ++halve) {
                const double scale = 1.0 / static_cast<double>(1 << halve);
                CameraExtrinsics trial = pose;
                trial.rotation_world_to_camera =
                    detail::axis_angle_to_matrix(dphi_full * scale) * trial.rotation_world_to_camera;
                trial.center_world_m = trial.center_world_m + dt_full * scale;
                double trial_cost = 0.0;
                bool ok = true;
                try {
                    const PnpLinearization trial_lin = pnp_linearize(
                        sub_pts, sub_obs, trial, intrinsics, robust, options.huber_delta_px);
                    for (double r : trial_lin.residuals) {
                        trial_cost += r * r;
                    }
                } catch (const std::domain_error&) {
                    ok = false;
                }
                if (ok && trial_cost < cost) {
                    pose = trial;
                    improved = true;
                    lambda = std::max(lambda * 0.1, 1.0e-9);
                    break;
                }
            }
            if (!improved) {
                lambda = std::min(lambda * 10.0, 1.0e7);
                if (lambda >= 1.0e7) {
                    break;  // converged as far as this basin allows
                }
                continue;
            }
            result.iterations++;
        }
        // Final cost.
        std::vector<math::Vector3> sub_pts;
        std::vector<math::Vector3> sub_obs;
        for (std::size_t k : indices) {
            sub_pts.push_back(points_world_m[k]);
            sub_obs.push_back(observed_px[k]);
        }
        double cost = 0.0;
        try {
            const PnpLinearization lin =
                pnp_linearize(sub_pts, sub_obs, pose, intrinsics, robust, options.huber_delta_px);
            for (double r : lin.residuals) {
                cost += r * r;
            }
        } catch (const std::domain_error&) {
            return {pose, -1.0};
        }
        return {pose, cost};
    };

    // Pass 1: all points (Huber-robust when enabled: outliers downweighted,
    // not obeyed). Pass 2 re-fits the hard-gated inliers WITHOUT robust
    // weights so the final RMS is an honest pixel statistic.
    std::vector<std::size_t> all_idx(points_world_m.size());
    for (std::size_t i = 0; i < all_idx.size(); ++i) {
        all_idx[i] = i;
    }
    auto [pose1, cost1] = solve_once(all_idx, current, options.robust_pass);
    if (cost1 < 0.0) {
        result.failure_reason = "diverged_or_singular";
        return result;
    }
    current = pose1;
    // Gate inliers on reprojection error.
    std::vector<std::size_t> inliers;
    double inlier_ss = 0.0;
    for (std::size_t i = 0; i < points_world_m.size(); ++i) {
        const PixelObservation pred =
            project_world_point(points_world_m[i], current, intrinsics);
        if (!pred.visible) {
            continue;
        }
        const double du = observed_px[i].x() - pred.u_px;
        const double dv = observed_px[i].y() - pred.v_px;
        const double err = std::sqrt(du * du + dv * dv);
        if (err <= options.inlier_gate_px) {
            inliers.push_back(i);
            inlier_ss += err * err;
        }
    }
    // Pass 2: re-fit on inliers when some were rejected (unweighted: the
    // outliers are gone, so Huber would only bias the honest RMS).
    if (inliers.size() < all_idx.size() && inliers.size() >= 4) {
        auto [pose2, cost2] = solve_once(inliers, current, false);
        if (cost2 >= 0.0) {
            current = pose2;
            inliers.clear();
            inlier_ss = 0.0;
            for (std::size_t i = 0; i < points_world_m.size(); ++i) {
                const PixelObservation pred =
                    project_world_point(points_world_m[i], current, intrinsics);
                if (!pred.visible) {
                    continue;
                }
                const double du = observed_px[i].x() - pred.u_px;
                const double dv = observed_px[i].y() - pred.v_px;
                const double err = std::sqrt(du * du + dv * dv);
                if (err <= options.inlier_gate_px) {
                    inliers.push_back(i);
                    inlier_ss += err * err;
                }
            }
        }
    }
    result.inliers = inliers.size();
    if (inliers.size() < 4) {
        result.failure_reason = "inlier_starvation";
        return result;
    }
    result.rms_reprojection_px = std::sqrt(inlier_ss / static_cast<double>(2 * inliers.size()));
    if (result.rms_reprojection_px > options.rms_accept_px) {
        result.failure_reason = "rms_above_tolerance";
        return result;
    }
    result.extrinsics = current;
    result.pose_valid = true;
    result.failure_reason = "none";
    result.range_m = current.center_world_m.norm();
    return result;
}

}  // namespace astradock::vision
