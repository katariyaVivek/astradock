// AstraDock M15 — Rotating-frame kinematics & relative orbital dynamics tests.
//
// Verification strategy (mirrors the M15 spec):
//   M15A LVLH rate: circular n = sqrt(mu/R^3), direction along orbit normal,
//     general-eccentricity identity omega = h/r^2, radial-velocity independence.
//   M15B transport: first-derivative split both directions, Coriolis/centrifugal/
//     Euler closed forms, full second-order composition, angular-acceleration
//     quotient rule, circular alpha = 0, rate-consistency residual.
//   M15C relative state: ECI->LVLH->ECI round trip, target/chaser antisymmetry,
//     zero separation, projected-vs-rotating velocity distinction.
//   M15D CW: analytical closed-form cases (along-track drift, radial oscillator,
//     cross-track oscillator, stationary hold), Phi(t=0) = identity, unforced
//     RK4-vs-closed-form agreement, forced constant-acceleration response,
//     invalid rejection.

#include "dynamics/two_body.hpp"
#include "frames/lvlh.hpp"
#include "frames/lvlh_rate.hpp"
#include "math/constants.hpp"
#include "math/matrix3.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"
#include "orbit/cartesian_state.hpp"
#include "orbit/two_body_orbit.hpp"
#include "relative/relative_state.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>

using namespace astradock;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

orbit::CartesianState circular_target(double radius_m, double mu) {
    const double v = std::sqrt(mu / radius_m);
    return {{radius_m, 0.0, 0.0}, {0.0, v, 0.0}};
}

}  // namespace

TEST_CASE("M15A LVLH rate equals circular mean motion along the orbit normal", "[relative][m15a]") {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const double R = constants::earth_reference_radius_m + 500.0e3;
    const orbit::CartesianState target = circular_target(R, mu);
    const math::Vector3 omega =
        frames::lvlh_angular_velocity_rad_s(target.position, target.velocity);
    const double n = std::sqrt(mu / (R * R * R));
    // Circular equatorial-equivalent start: motion in +Y at +X => normal is +Z.
    CHECK_THAT(omega.x(), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(omega.y(), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(omega.z(), WithinRel(n, 1.0e-12));
    CHECK_THAT(frames::lvlh_rate_rad_per_s(target.position, target.velocity), WithinRel(n, 1.0e-12));
    // Matches the canonical circular-orbit reference mean motion exactly.
    const auto ref = orbit::compute_circular_orbit_reference(mu, R);
    CHECK_THAT(omega.norm(), WithinRel(ref.mean_motion_rad_per_s, 1.0e-12));
}

TEST_CASE("M15A LVLH rate identity holds for eccentric states: omega = h / r^2", "[relative][m15a]") {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    // Eccentric LEO: perigee 400 km with 800 m/s super-circular horizontal velocity.
    const double rp = constants::earth_reference_radius_m + 400.0e3;
    const math::Vector3 r{rp, 0.0, 0.0};
    const math::Vector3 v{0.0, std::sqrt(mu / rp) + 800.0, 0.0};
    const math::Vector3 omega = frames::lvlh_angular_velocity_rad_s(r, v);
    const math::Vector3 h = r.cross(v);
    const math::Vector3 expected = h / r.squared_norm();
    CHECK_THAT(omega.x(), WithinAbs(expected.x(), 1.0e-18));
    CHECK_THAT(omega.y(), WithinAbs(expected.y(), 1.0e-18));
    CHECK_THAT(omega.z(), WithinRel(expected.z(), 1.0e-12));
    // Radial velocity does not affect the rate: add 500 m/s radial component.
    const math::Vector3 v_radial{500.0, v.y(), 0.0};
    const math::Vector3 omega2 = frames::lvlh_angular_velocity_rad_s(r, v_radial);
    CHECK_THAT(omega2.z(), WithinRel(omega.z(), 1.0e-12));
    // Inclined orbit: rate direction is parallel to h, not blindly +Z.
    const math::Vector3 ri{7000.0e3, 100.0e3, 200.0e3};
    const math::Vector3 vi{-100.0, 7400.0, 300.0};
    const math::Vector3 oi = frames::lvlh_angular_velocity_rad_s(ri, vi);
    const math::Vector3 hi = ri.cross(vi).normalized();
    CHECK_THAT(oi.normalized().dot(hi), WithinRel(1.0, 1.0e-12));
    // Degenerate states rejected (mirrors M06 LVLH policy).
    CHECK_THROWS_AS(
        frames::lvlh_angular_velocity_rad_s(math::Vector3{}, math::Vector3{0.0, 1.0, 0.0}),
        std::domain_error);
    CHECK_THROWS_AS(
        frames::lvlh_angular_velocity_rad_s(
            math::Vector3{7000.0e3, 0.0, 0.0}, math::Vector3{8000.0, 0.0, 0.0}),
        std::domain_error);
}

TEST_CASE("M15B transport theorem splits inertial and rotating derivatives", "[relative][m15b]") {
    // Hand case: omega = [0,0,1], v = [1,0,0], rotating derivative = 0.
    // Inertial derivative = omega x v = [0,1,0].
    const math::Vector3 omega{0.0, 0.0, 1.0};
    const math::Vector3 v{1.0, 0.0, 0.0};
    const math::Vector3 zero{};
    const math::Vector3 inertial = frames::transport_first_derivative(v, omega, zero);
    CHECK_THAT(inertial.x(), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(inertial.y(), WithinRel(1.0, 1.0e-15));
    CHECK_THAT(inertial.z(), WithinAbs(0.0, 1.0e-15));
    // Inverse split recovers the rotating derivative.
    const math::Vector3 back = frames::rotating_frame_derivative(v, omega, inertial);
    CHECK_THAT(back.norm(), WithinAbs(0.0, 1.0e-15));
    // Nonzero rotating derivative adds linearly.
    const math::Vector3 drot{0.0, 0.0, 2.0};
    const math::Vector3 inertial2 = frames::transport_first_derivative(v, omega, drot);
    CHECK_THAT(inertial2.y(), WithinRel(1.0, 1.0e-15));
    CHECK_THAT(inertial2.z(), WithinRel(2.0, 1.0e-15));
    // Coriolis / centrifugal / Euler closed forms.
    const math::Vector3 cor = frames::coriolis_acceleration(omega, math::Vector3{0.0, 1.0, 0.0});
    CHECK_THAT(cor.x(), WithinRel(-2.0, 1.0e-15));  // 2 z-hat x y-hat = -2 x-hat
    const math::Vector3 cen = frames::centrifugal_acceleration(omega, math::Vector3{1.0, 0.0, 0.0});
    CHECK_THAT(cen.x(), WithinRel(-1.0, 1.0e-15));  // z x (z x x) = -x
    const math::Vector3 eul =
        frames::euler_acceleration(math::Vector3{0.0, 0.0, 3.0}, math::Vector3{1.0, 0.0, 0.0});
    CHECK_THAT(eul.y(), WithinRel(3.0, 1.0e-15));  // alpha z x x = +y
}

TEST_CASE("M15B full acceleration composition matches the textbook identity", "[relative][m15b]") {
    const math::Vector3 omega{0.0, 0.0, 0.1};
    const math::Vector3 alpha{0.0, 0.0, 0.01};
    const math::Vector3 r_rel{100.0, 0.0, 0.0};
    const math::Vector3 v_rel{0.0, 5.0, 0.0};
    const math::Vector3 a_rel{0.1, 0.0, 0.0};
    const math::Vector3 a_origin{0.0, -8.0, 0.0};
    const math::Vector3 composed = frames::compose_inertial_acceleration(
        r_rel, v_rel, a_rel, omega, alpha, a_origin);
    // Independent term-by-term reconstruction:
    //   a_rel = [0.1,0,0]; alpha x r = [0,1,0]; 2w x v = [-1,0,0];
    //   w x (w x r) = [-1,0,0]; origin = [0,-8,0].
    const math::Vector3 expected{-1.9, -7.0, 0.0};
    CHECK_THAT(composed.x(), WithinRel(expected.x(), 1.0e-12));
    CHECK_THAT(composed.y(), WithinRel(expected.y(), 1.0e-12));
    CHECK_THAT(composed.z(), WithinAbs(0.0, 1.0e-15));
}

TEST_CASE("M15B LVLH angular acceleration vanishes on circular Keplerian arcs", "[relative][m15b]") {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const double R = constants::earth_reference_radius_m + 500.0e3;
    const orbit::CartesianState target = circular_target(R, mu);
    const math::Vector3 a_ref = dynamics::two_body_acceleration(target.position, mu);
    const math::Vector3 alpha = frames::lvlh_angular_acceleration_rad_per_s2(
        target.position, target.velocity, a_ref);
    CHECK_THAT(alpha.norm(), WithinAbs(0.0, 1.0e-15));
    // Rate-consistency residual: v = radial-rate + omega x r holds exactly.
    const math::Vector3 residual = frames::lvlh_rate_consistency_residual(target.position, target.velocity);
    CHECK_THAT(residual.norm(), WithinAbs(0.0, 1.0e-9));
    // Eccentric case: residual still zero (identity), alpha nonzero in general.
    const math::Vector3 re{7000.0e3, 0.0, 0.0};
    const math::Vector3 ve{1000.0, 8000.0, 0.0};
    const math::Vector3 ae = dynamics::two_body_acceleration(re, mu);
    CHECK_THAT(
        frames::lvlh_rate_consistency_residual(re, ve).norm(), WithinAbs(0.0, 1.0e-9));
    CHECK(frames::lvlh_angular_acceleration_rad_per_s2(re, ve, ae).norm() > 0.0);
}

TEST_CASE("M15C LVLH relative state round-trips ECI exactly", "[relative][m15c]") {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const orbit::CartesianState target = circular_target(6878137.0, mu);
    const orbit::CartesianState chaser{
        target.position + math::Vector3{100.0, -2000.0, 50.0},
        target.velocity + math::Vector3{-1.0, 3.0, 0.5}};
    const relative::RelativeStateLvlh rel = relative::relative_state_from_eci(target, chaser);
    REQUIRE(math::is_finite(rel.relative_position_lvlh_m));
    const orbit::CartesianState recovered = relative::chaser_eci_from_relative(target, rel);
    CHECK_THAT(recovered.position.x(), WithinRel(chaser.position.x(), 1.0e-12));
    CHECK_THAT(recovered.position.y(), WithinRel(chaser.position.y(), 1.0e-12));
    CHECK_THAT(recovered.position.z(), WithinRel(chaser.position.z(), 1.0e-12));
    CHECK_THAT(recovered.velocity.x(), WithinRel(chaser.velocity.x(), 1.0e-9));
    CHECK_THAT(recovered.velocity.y(), WithinRel(chaser.velocity.y(), 1.0e-9));
    CHECK_THAT(recovered.velocity.z(), WithinRel(chaser.velocity.z(), 1.0e-9));
    // Zero separation: chaser co-located with target.
    const relative::RelativeStateLvlh zero = relative::relative_state_from_eci(target, target);
    CHECK_THAT(zero.relative_position_lvlh_m.norm(), WithinAbs(0.0, 1.0e-9));
    CHECK_THAT(zero.relative_velocity_lvlh_mps.norm(), WithinAbs(0.0, 1.0e-9));
}

TEST_CASE("M15C relative velocity is the rotating-frame derivative, not the projection", "[relative][m15c]") {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const orbit::CartesianState target = circular_target(6878137.0, mu);
    // Chaser displaced +2 km along-track with identical inertial velocity.
    // Projection C*dv = 0, but the rotating derivative is -omega x rho != 0.
    const math::Matrix3 c = frames::dcm_lvlh_from_eci(target.position, target.velocity);
    const orbit::CartesianState chaser{target.position + c.transpose() * math::Vector3{0.0, 2000.0, 0.0},
                                       target.velocity};
    const relative::RelativeStateLvlh rel = relative::relative_state_from_eci(target, chaser);
    CHECK_THAT(rel.relative_position_lvlh_m.y(), WithinRel(2000.0, 1.0e-9));
    const math::Vector3 omega_lvlh = c * frames::lvlh_angular_velocity_rad_s(target.position, target.velocity);
    const math::Vector3 expected_v = omega_lvlh.cross(rel.relative_position_lvlh_m) * (-1.0);
    CHECK_THAT(rel.relative_velocity_lvlh_mps.x(), WithinRel(expected_v.x(), 1.0e-9));
    CHECK_THAT(rel.relative_velocity_lvlh_mps.y(), WithinRel(expected_v.y(), 1.0e-9));
    // Naive projection would report ~0; the corrected value is +n*rho in x
    // (the chaser carries excess radial velocity vs local circular: same
    // inertial velocity at a larger orbital angle means outward drift in LVLH).
    const double n = std::sqrt(mu / (6878137.0 * 6878137.0 * 6878137.0));
    CHECK_THAT(rel.relative_velocity_lvlh_mps.x(), WithinRel(n * 2000.0, 1.0e-6));
}

TEST_CASE("M15D CW closed form reproduces textbook analytical cases", "[relative][m15d]") {
    const double n = 1.0e-3;
    // Case 1: pure along-track offset drifts at -3n*x? No: pure y0 holds y, x=0 static.
    const relative::RelativeStateLvlh drift{{0.0, 1000.0, 0.0}, {0.0, 0.0, 0.0}};
    const auto drifted = relative::cw_predict(drift, n, 1000.0);
    CHECK_THAT(drifted.relative_position_lvlh_m.x(), WithinAbs(0.0, 1.0e-9));
    CHECK_THAT(drifted.relative_position_lvlh_m.y(), WithinRel(1000.0, 1.0e-12));
    // Case 2: radial offset x0 induces secular along-track drift -6n x0 (t - sin(nt)/n)... verify numerically.
    const relative::RelativeStateLvlh radial{{100.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    const auto grown = relative::cw_predict(radial, n, 1000.0);
    const double c = std::cos(n * 1000.0);
    const double s = std::sin(n * 1000.0);
    CHECK_THAT(grown.relative_position_lvlh_m.x(), WithinRel((4.0 - 3.0 * c) * 100.0, 1.0e-12));
    CHECK_THAT(
        grown.relative_position_lvlh_m.y(), WithinRel(6.0 * (s - n * 1000.0) * 100.0, 1.0e-12));
    // Case 3: cross-track oscillator z(t) = z0 cos(nt), vz(t) = -n z0 sin(nt).
    const relative::RelativeStateLvlh cross{{0.0, 0.0, 50.0}, {0.0, 0.0, 0.0}};
    const auto osc = relative::cw_predict(cross, n, 500.0);
    CHECK_THAT(osc.relative_position_lvlh_m.z(), WithinRel(50.0 * std::cos(n * 500.0), 1.0e-12));
    CHECK_THAT(
        osc.relative_velocity_lvlh_mps.z(), WithinRel(-n * 50.0 * std::sin(n * 500.0), 1.0e-9));
    // Case 4: bounded relative ellipse condition vy0 = -2n x0 gives periodic x/y.
    const relative::RelativeStateLvlh ellipse{{100.0, 0.0, 0.0}, {0.0, -2.0 * n * 100.0, 0.0}};
    const double period = 2.0 * constants::pi / n;
    const auto closed = relative::cw_predict(ellipse, n, period);
    CHECK_THAT(closed.relative_position_lvlh_m.x(), WithinRel(100.0, 1.0e-9));
    CHECK_THAT(closed.relative_position_lvlh_m.y(), WithinAbs(0.0, 1.0e-6));
    // Phi(0) is the identity.
    const auto same = relative::cw_predict(ellipse, n, 0.0);
    CHECK_THAT(same.relative_position_lvlh_m.x(), WithinRel(100.0, 1.0e-15));
    CHECK_THAT(same.relative_velocity_lvlh_mps.y(), WithinRel(-2.0 * n * 100.0, 1.0e-15));
}

TEST_CASE("M15D CW acceleration matches the governing differential equations", "[relative][m15d]") {
    const double n = 1.1e-3;
    const relative::RelativeStateLvlh s{{10.0, -20.0, 5.0}, {1.0, 2.0, -0.5}};
    const auto d = relative::cw_acceleration(s, n);
    // Position slot = state kinematics d(pos)/dt = vel.
    CHECK_THAT(d.relative_position_lvlh_m.x(), WithinRel(1.0, 1.0e-15));
    CHECK_THAT(d.relative_position_lvlh_m.y(), WithinRel(2.0, 1.0e-15));
    CHECK_THAT(d.relative_position_lvlh_m.z(), WithinRel(-0.5, 1.0e-15));
    // Velocity slot = acceleration: x_ddot = 3n^2 x + 2n vy; y_ddot = -2n vx; z_ddot = -n^2 z.
    CHECK_THAT(d.relative_velocity_lvlh_mps.x(), WithinRel(3.0 * n * n * 10.0 + 2.0 * n * 2.0, 1.0e-12));
    CHECK_THAT(d.relative_velocity_lvlh_mps.y(), WithinRel(-2.0 * n * 1.0, 1.0e-12));
    CHECK_THAT(d.relative_velocity_lvlh_mps.z(), WithinRel(-n * n * 5.0, 1.0e-12));
    // Input acceleration adds linearly to the velocity slot only.
    const math::Vector3 push{0.01, -0.02, 0.03};
    const auto d_forced = relative::cw_acceleration(s, n, push);
    const auto d_free = relative::cw_acceleration(s, n);
    CHECK_THAT(
        (d_forced.relative_position_lvlh_m - d_free.relative_position_lvlh_m).norm(),
        WithinAbs(0.0, 1.0e-15));
    // Forced RK4 propagation under constant along-track accel: pushing -Y for
    // 600 s must move the chaser -Y relative to the unforced arc.
    const auto traj = relative::cw_propagate(s, n, 1.0, 600, push);
    const auto drift = relative::cw_propagate(s, n, 1.0, 600);
    REQUIRE(traj.size() == 601);
    CHECK(traj.back().relative_position_lvlh_m.y() < drift.back().relative_position_lvlh_m.y() - 1.0);
}

TEST_CASE("M15D unforced CW RK4 agrees with the exact closed form", "[relative][m15d]") {
    const double n = 1.0e-3;
    const relative::RelativeStateLvlh s{{100.0, -50.0, 25.0}, {0.5, -0.2, 0.1}};
    const auto exact = relative::cw_predict(s, n, 3600.0);
    const auto traj = relative::cw_propagate(s, n, 1.0, 3600);
    const auto& num = traj.back();
    CHECK_THAT(num.relative_position_lvlh_m.x(), WithinRel(exact.relative_position_lvlh_m.x(), 1.0e-6));
    CHECK_THAT(num.relative_position_lvlh_m.y(), WithinRel(exact.relative_position_lvlh_m.y(), 1.0e-6));
    CHECK_THAT(num.relative_position_lvlh_m.z(), WithinRel(exact.relative_position_lvlh_m.z(), 1.0e-6));
    CHECK_THAT(
        num.relative_velocity_lvlh_mps.x(), WithinRel(exact.relative_velocity_lvlh_mps.x(), 1.0e-6));
}

TEST_CASE("M15 inputs reject invalid relative states and parameters", "[relative][invalid]") {
    const relative::RelativeStateLvlh s{{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    CHECK_THROWS_AS(relative::cw_predict(s, 0.0, 10.0), std::domain_error);
    CHECK_THROWS_AS(relative::cw_predict(s, -1.0e-3, 10.0), std::domain_error);
    CHECK_THROWS_AS(
        relative::cw_propagate(s, 1.0e-3, 0.0, 10),
        std::domain_error);
    CHECK_THROWS_AS(
        relative::cw_propagate(s, 1.0e-3, -1.0, 10),
        std::domain_error);
    const relative::RelativeStateLvlh bad{{1.0, 0.0, 0.0},
                                          {std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0}};
    CHECK_THROWS_AS(relative::cw_predict(bad, 1.0e-3, 10.0), std::domain_error);
}
