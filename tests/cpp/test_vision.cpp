// AstraDock M22 — Computer vision & optical navigation tests.
//
// Verification strategy:
//   Camera (M22A-C): principal-point centering hand case, off-axis analytic
//     pixels, behind-camera + out-of-bounds flags, distortion direction,
//     bearing round-trip, reprojection error, invalid rejection.
//   Fiducials: default plate size/count/non-coplanarity.
//   PnP (M22D-G): exact recovery from noiseless pixels; noisy recovery within
//     bounds; outlier rejection + re-fit; inlier starvation reports invalid;
//     pixel-noise sensitivity sweep; RMS gate enforcement.

#include "math/matrix3.hpp"
#include "math/vector3.hpp"
#include "vision/camera_model.hpp"
#include "vision/pnp_pose.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

using namespace astradock;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

vision::CameraIntrinsics test_cam() {
    vision::CameraIntrinsics cam;
    cam.fx_px = 800.0;
    cam.fy_px = 800.0;
    cam.cx_px = 320.0;
    cam.cy_px = 240.0;
    cam.width_px = 640;
    cam.height_px = 480;
    cam.k1_radial = 0.0;
    cam.near_clip_m = 0.1;
    return cam;
}

vision::CameraExtrinsics identity_pose() {
    vision::CameraExtrinsics extr;
    extr.rotation_world_to_camera = math::Matrix3::identity();
    extr.center_world_m = math::Vector3{};
    return extr;
}

}  // namespace

TEST_CASE("M22 optical axis projects to the principal point", "[vision][camera]") {
    const auto cam = test_cam();
    // Point on the boresight at 10 m: (320, 240) exactly.
    const auto obs = vision::project_camera_point(math::Vector3{0.0, 0.0, 10.0}, cam);
    CHECK(obs.visible);
    CHECK_THAT(obs.u_px, WithinRel(320.0, 1.0e-12));
    CHECK_THAT(obs.v_px, WithinRel(240.0, 1.0e-12));
    // Off-axis analytic: (1, -0.5, 10) -> u = 800*0.1+320 = 400, v = 800*(-0.05)+240 = 200.
    const auto off = vision::project_camera_point(math::Vector3{1.0, -0.5, 10.0}, cam);
    CHECK(off.visible);
    CHECK_THAT(off.u_px, WithinRel(400.0, 1.0e-12));
    CHECK_THAT(off.v_px, WithinRel(200.0, 1.0e-12));
}

TEST_CASE("M22 behind-camera and out-of-bounds flag exactly one reason", "[vision][camera]") {
    const auto cam = test_cam();
    const auto behind = vision::project_camera_point(math::Vector3{0.0, 0.0, -5.0}, cam);
    CHECK_FALSE(behind.visible);
    CHECK(behind.behind_camera);
    CHECK_FALSE(behind.out_of_bounds);
    const auto at_clip = vision::project_camera_point(math::Vector3{0.0, 0.0, 0.05}, cam);
    CHECK_FALSE(at_clip.visible);
    CHECK(at_clip.behind_camera);
    // Far off-axis: projects outside 640x480.
    const auto oob = vision::project_camera_point(math::Vector3{10.0, 0.0, 5.0}, cam);
    CHECK_FALSE(oob.visible);
    CHECK_FALSE(oob.behind_camera);
    CHECK(oob.out_of_bounds);
}

TEST_CASE("M22 radial distortion displaces outward for barrel lenses", "[vision][camera]") {
    auto cam = test_cam();
    cam.k1_radial = 0.1;
    const auto plain = test_cam();
    const math::Vector3 pt{2.0, 0.0, 10.0};
    const auto distorted = vision::project_camera_point(pt, cam);
    const auto undistorted = vision::project_camera_point(pt, plain);
    CHECK(distorted.visible);
    // x_n = 0.2, r^2 = 0.04, s = 1.004: u grows past the undistorted 480.
    CHECK_THAT(undistorted.u_px, WithinRel(480.0, 1.0e-12));
    CHECK(distorted.u_px > undistorted.u_px);
    CHECK_THAT(distorted.u_px, WithinRel(320.0 + 800.0 * 1.004 * 0.2, 1.0e-9));
    // Bearing round-trips through projection at k1 = 0.
    const math::Vector3 bearing = vision::pixel_to_bearing(400.0, 200.0, plain);
    CHECK_THAT(bearing.x() / bearing.z(), WithinRel(0.1, 1.0e-12));
    CHECK_THAT(bearing.y() / bearing.z(), WithinRel(-0.05, 1.0e-12));
}

TEST_CASE("M22 world-to-camera translation is subtracted before rotation", "[vision][camera]") {
    const auto cam = test_cam();
    vision::CameraExtrinsics extr;
    extr.rotation_world_to_camera = math::Matrix3::identity();
    extr.center_world_m = math::Vector3{0.0, 0.0, -10.0};
    // World origin sits 10 m along the boresight: principal point.
    const auto obs = vision::project_world_point(math::Vector3{}, extr, cam);
    CHECK(obs.visible);
    CHECK_THAT(obs.u_px, WithinRel(320.0, 1.0e-12));
    CHECK_THAT(obs.v_px, WithinRel(240.0, 1.0e-12));
    // Reprojection error of the exact projection is zero.
    CHECK_THAT(
        vision::reprojection_error_px(obs.u_px, obs.v_px, math::Vector3{}, extr, cam),
        WithinAbs(0.0, 1.0e-9));
    // A world point BEHIND the camera (camera-frame z < 0) is unprojectable.
    vision::CameraExtrinsics flipped = extr;
    flipped.center_world_m = math::Vector3{0.0, 0.0, 10.0};  // camera past the point
    CHECK_THROWS_AS(
        vision::reprojection_error_px(320.0, 240.0, math::Vector3{}, flipped, cam),
        std::domain_error);
}

TEST_CASE("M22 default fiducial plate is sized and non-coplanar", "[vision][fiducials]") {
    const auto plate = vision::default_fiducial_plate();
    plate.validate();
    REQUIRE(plate.points_world_m.size() == 5);
    // 1 m square corners: edge length exactly 1.
    CHECK_THAT(
        (plate.points_world_m[0] - plate.points_world_m[1]).norm(), WithinRel(1.0, 1.0e-15));
    // Standoff point breaks planarity (z = 0.25).
    CHECK_THAT(plate.points_world_m[4].z(), WithinRel(0.25, 1.0e-15));
    // Few-point sets rejected for PnP.
    vision::FiducialSet small;
    small.points_world_m = {math::Vector3{}, math::Vector3{1.0, 0.0, 0.0}};
    CHECK_THROWS_AS(small.validate(), std::domain_error);
}

TEST_CASE("M22 PnP recovers the exact pose from noiseless pixels", "[vision][pnp]") {
    const auto cam = test_cam();
    const auto plate = vision::default_fiducial_plate();
    // True pose: camera 8 m back on -Z looking +Z at the plate (R = I... plate
    // at z=0 faces +Z; camera center (0,0,-8) with identity rotation sees +Z).
    vision::CameraExtrinsics truth;
    truth.rotation_world_to_camera = math::Matrix3::identity();
    truth.center_world_m = math::Vector3{0.0, 0.0, -8.0};
    std::vector<math::Vector3> pixels;
    for (const auto& point : plate.points_world_m) {
        const auto obs = vision::project_world_point(point, truth, cam);
        REQUIRE(obs.visible);
        pixels.push_back(math::Vector3{obs.u_px, obs.v_px, 0.0});
    }
    // Seed: 3-deg rotation + 0.5 m translation perturbation.
    vision::CameraExtrinsics seed = truth;
    seed.center_world_m = seed.center_world_m + math::Vector3{0.3, -0.2, 0.4};
    const auto estimate = vision::estimate_pose_pnp(plate.points_world_m, pixels, seed, cam);
    CHECK(estimate.pose_valid);
    CHECK(estimate.inliers == 5);
    CHECK_THAT(estimate.rms_reprojection_px, WithinAbs(0.0, 1.0e-6));
    CHECK_THAT(
        (estimate.extrinsics.center_world_m - truth.center_world_m).norm(), WithinAbs(0.0, 1.0e-6));
    CHECK_THAT(estimate.range_m, WithinRel(8.0, 1.0e-9));
}

TEST_CASE("M22 PnP recovers pose under pixel noise within bounds", "[vision][pnp]") {
    const auto cam = test_cam();
    const auto plate = vision::default_fiducial_plate();
    vision::CameraExtrinsics truth;
    truth.rotation_world_to_camera = math::Matrix3::identity();
    truth.center_world_m = math::Vector3{0.5, -0.3, -12.0};
    // Deterministic 0.5 px noise pattern (no RNG in unit tests).
    const std::vector<math::Vector3> noise{
        {0.5, -0.3, 0.0}, {-0.4, 0.2, 0.0}, {0.3, 0.4, 0.0}, {-0.2, -0.5, 0.0}, {0.1, 0.1, 0.0}};
    std::vector<math::Vector3> pixels;
    for (std::size_t i = 0; i < plate.points_world_m.size(); ++i) {
        const auto obs = vision::project_world_point(plate.points_world_m[i], truth, cam);
        REQUIRE(obs.visible);
        pixels.push_back(
            math::Vector3{obs.u_px + noise[i].x(), obs.v_px + noise[i].y(), 0.0});
    }
    const auto estimate = vision::estimate_pose_pnp(plate.points_world_m, pixels, truth, cam);
    CHECK(estimate.pose_valid);
    CHECK(estimate.rms_reprojection_px < 2.0);
    // Depth dilution is physics: dz ~ z^2/(f*b) * dpx = 144/800 * 0.5 ~= 0.09 m
    // at 12 m, plus lateral terms — bound at 0.25 m, not 0.10.
    CHECK((estimate.extrinsics.center_world_m - truth.center_world_m).norm() < 0.25);
}

TEST_CASE("M22 PnP rejects an outlier and still converges", "[vision][outlier]") {
    const auto cam = test_cam();
    // 6-point set (plate + extra) so one outlier still leaves 5 inliers.
    std::vector<math::Vector3> points = vision::default_fiducial_plate().points_world_m;
    points.push_back(math::Vector3{0.25, 0.25, 0.1});
    vision::CameraExtrinsics truth;
    truth.rotation_world_to_camera = math::Matrix3::identity();
    truth.center_world_m = math::Vector3{0.0, 0.0, -8.0};
    std::vector<math::Vector3> pixels;
    for (std::size_t i = 0; i < points.size(); ++i) {
        const auto obs = vision::project_world_point(points[i], truth, cam);
        REQUIRE(obs.visible);
        if (i == 2) {
            pixels.push_back(math::Vector3{obs.u_px + 50.0, obs.v_px - 40.0, 0.0});
        } else {
            pixels.push_back(math::Vector3{obs.u_px, obs.v_px, 0.0});
        }
    }
    const auto estimate = vision::estimate_pose_pnp(points, pixels, truth, cam);
    CHECK(estimate.pose_valid);
    CHECK(estimate.inliers == 5);
    CHECK_THAT(
        (estimate.extrinsics.center_world_m - truth.center_world_m).norm(), WithinAbs(0.0, 1.0e-6));
}

TEST_CASE("M22 PnP reports invalid instead of a confident wrong answer", "[vision][pnp]") {
    const auto cam = test_cam();
    const auto plate = vision::default_fiducial_plate();
    vision::CameraExtrinsics truth;
    truth.rotation_world_to_camera = math::Matrix3::identity();
    truth.center_world_m = math::Vector3{0.0, 0.0, -8.0};
    // All pixels scrambled: inlier starvation -> invalid with reason.
    std::vector<math::Vector3> pixels;
    for (std::size_t i = 0; i < plate.points_world_m.size(); ++i) {
        pixels.push_back(math::Vector3{10.0 + 20.0 * i, 400.0 - 15.0 * i, 0.0});
    }
    const auto estimate = vision::estimate_pose_pnp(plate.points_world_m, pixels, truth, cam);
    CHECK_FALSE(estimate.pose_valid);
    CHECK(estimate.inliers < 4);
}

TEST_CASE("M22 PnP degrades gracefully with pixel-noise level", "[vision][sensitivity]") {
    const auto cam = test_cam();
    const auto plate = vision::default_fiducial_plate();
    vision::CameraExtrinsics truth;
    truth.rotation_world_to_camera = math::Matrix3::identity();
    truth.center_world_m = math::Vector3{0.0, 0.0, -10.0};
    double prev_err = 0.0;
    for (double sigma : {0.0, 0.25, 0.5, 1.0}) {
        std::vector<math::Vector3> pixels;
        for (std::size_t i = 0; i < plate.points_world_m.size(); ++i) {
            const auto obs = vision::project_world_point(plate.points_world_m[i], truth, cam);
            REQUIRE(obs.visible);
            const double s = (i % 2 == 0) ? sigma : -sigma;
            pixels.push_back(math::Vector3{obs.u_px + s, obs.v_px - s, 0.0});
        }
        const auto estimate = vision::estimate_pose_pnp(plate.points_world_m, pixels, truth, cam);
        REQUIRE(estimate.pose_valid);
        const double err =
            (estimate.extrinsics.center_world_m - truth.center_world_m).norm();
        CHECK(err >= prev_err - 1.0e-9);  // monotone non-decreasing in noise
        prev_err = err;
        if (sigma == 1.0) {
            // NOTE: the alternating ± pattern is a sawtooth SHEAR (systematic),
            // not zero-mean noise — the least-squares answer correctly absorbs
            // ~0.3 m at 10 m range. Random Gaussian noise does better (see the
            // demo sensitivity sweep). Bound documents shear response.
            CHECK(err < 0.5);
        }
    }
}

TEST_CASE("M22 camera inputs reject invalid configurations", "[vision][invalid]") {
    auto cam = test_cam();
    cam.fx_px = 0.0;
    CHECK_THROWS_AS(
        vision::project_camera_point(math::Vector3{0.0, 0.0, 5.0}, cam), std::domain_error);
    vision::CameraExtrinsics extr;
    extr.rotation_world_to_camera = math::Matrix3(
        1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0);  // not a rotation
    CHECK_THROWS_AS(
        vision::project_world_point(math::Vector3{}, extr, test_cam()), std::domain_error);
    const auto plate = vision::default_fiducial_plate();
    const std::vector<math::Vector3> one{math::Vector3{320.0, 240.0, 0.0}};
    CHECK_THROWS_AS(
        vision::estimate_pose_pnp(
            {plate.points_world_m[0]}, one, identity_pose(), test_cam()),
        std::domain_error);
}
