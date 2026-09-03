// AstraDock M16 — Guidance & Control Foundations tests.
//
// Verification strategy (mirrors the M16 spec):
//   M16A attitude PD: error-vector hand cases per axis, double-cover sign
//     alignment, geodesic angle, torque law linearity, never raw subtraction.
//   M16B rate damping: principal-axis analytical cases, command tracking.
//   M16C translation: CW feedforward closed form, PD error response, station
//     keeping at origin, acceleration linearity.
//   M16D LQR: closed-form Riccati solution, ARE residual, Hurwitz check,
//     damping ratio, invalid rejection.
//   M16E saturation: desired/achieved/deficit accounting, per-axis clamp.
//   Metrics: PD-driven 6-DOF detumble settles; saturated torque stays bounded.

#include "attitude/principal_inertia.hpp"
#include "attitude/rigid_body.hpp"
#include "control/attitude_pd.hpp"
#include "control/control_saturation.hpp"
#include "control/double_integrator_lqr.hpp"
#include "control/relative_pd.hpp"
#include "math/constants.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"
#include "relative/relative_state.hpp"
#include "spacecraft/six_dof_dynamics.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace astradock;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("M16A attitude error vector matches small-angle hand cases per axis", "[control][m16a]") {
    // 10 deg rotation about +X: q_err = [cos5°, sin5°, 0, 0], vector part [sin5°,0,0].
    const math::Quaternion q_des = math::Quaternion::identity();
    const math::Quaternion q_cur = math::Quaternion::from_axis_angle(
        math::Vector3{1.0, 0.0, 0.0}, 10.0 * constants::pi / 180.0);
    const math::Vector3 e = control::attitude_error_vector(q_cur, q_des);
    CHECK_THAT(e.x(), WithinRel(std::sin(5.0 * constants::pi / 180.0), 1.0e-12));
    CHECK_THAT(e.y(), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(e.z(), WithinAbs(0.0, 1.0e-15));
    // Sign: rotation about -Y gives negative Y error.
    const math::Quaternion q_neg = math::Quaternion::from_axis_angle(
        math::Vector3{0.0, 1.0, 0.0}, -6.0 * constants::pi / 180.0);
    const math::Vector3 e_neg = control::attitude_error_vector(q_neg, q_des);
    CHECK(e_neg.y() < 0.0);
    CHECK_THAT(e_neg.y(), WithinRel(-std::sin(3.0 * constants::pi / 180.0), 1.0e-12));
    // Identity: zero error.
    const math::Vector3 e_zero = control::attitude_error_vector(q_des, q_des);
    CHECK_THAT(e_zero.norm(), WithinAbs(0.0, 1.0e-15));
    // Geodesic angle: 10 deg error reads 10 deg.
    CHECK_THAT(
        control::attitude_error_angle_rad(q_cur, q_des),
        WithinRel(10.0 * constants::pi / 180.0, 1.0e-12));
}

TEST_CASE("M16A error vector is invariant under quaternion double cover", "[control][m16a]") {
    const math::Quaternion q_des = math::Quaternion::identity();
    const math::Quaternion q_cur = math::Quaternion::from_axis_angle(
        math::Vector3{0.0, 0.0, 1.0}, 30.0 * constants::pi / 180.0);
    const math::Quaternion q_flip(-q_cur.w(), -q_cur.x(), -q_cur.y(), -q_cur.z());
    const math::Vector3 e1 = control::attitude_error_vector(q_cur, q_des);
    const math::Vector3 e2 = control::attitude_error_vector(q_flip, q_des);
    CHECK_THAT((e1 - e2).norm(), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(
        control::attitude_error_angle_rad(q_flip, q_des),
        WithinRel(30.0 * constants::pi / 180.0, 1.0e-12));
}

TEST_CASE("M16A PD torque law is linear in attitude and rate errors", "[control][m16a]") {
    const control::AttitudePdGains gains{{2.0, 3.0, 4.0}, {0.5, 0.6, 0.7}};
    const math::Quaternion q_des = math::Quaternion::identity();
    const math::Quaternion q_cur = math::Quaternion::from_axis_angle(
        math::Vector3{1.0, 0.0, 0.0}, 4.0 * constants::pi / 180.0);
    const math::Vector3 w{0.1, -0.2, 0.3};
    const math::Vector3 w_des{};
    const math::Vector3 tau =
        control::attitude_pd_torque_body_Nm(q_cur, q_des, w, w_des, gains);
    const math::Vector3 e = control::attitude_error_vector(q_cur, q_des);
    CHECK_THAT(tau.x(), WithinRel(-2.0 * e.x() - 0.5 * 0.1, 1.0e-12));
    CHECK_THAT(tau.y(), WithinRel(-3.0 * e.y() - 0.6 * (-0.2), 1.0e-12));
    CHECK_THAT(tau.z(), WithinRel(-4.0 * e.z() - 0.7 * 0.3, 1.0e-12));
    // Zero error + zero rate error: zero torque.
    const math::Vector3 tau_zero =
        control::attitude_pd_torque_body_Nm(q_des, q_des, w_des, w_des, gains);
    CHECK_THAT(tau_zero.norm(), WithinAbs(0.0, 1.0e-15));
    // Restoring direction: positive X error demands negative X torque.
    CHECK(tau.x() < 0.0);
    // Gain suggestion has the documented engineering basis (zeta, settle time).
    const control::PdTuningRequest req{{10.0, 20.0, 30.0}, 0.9, 30.0};
    const control::AttitudePdGains suggested = control::suggest_pd_gains(req);
    const double wn = 4.0 / (0.9 * 30.0);
    CHECK_THAT(suggested.kp_Nm_per_rad.x(), WithinRel(2.0 * 10.0 * wn * wn, 1.0e-12));
    CHECK_THAT(suggested.kd_Nm_s_per_rad.z(), WithinRel(2.0 * 0.9 * wn * 30.0, 1.0e-12));
}

TEST_CASE("M16B rate damping matches principal-axis analytical cases", "[control][m16b]") {
    const math::Vector3 kw{1.5, 2.5, 3.5};
    const math::Vector3 w{0.2, -0.1, 0.15};
    const math::Vector3 zero{};
    const math::Vector3 tau = control::rate_damping_torque_body_Nm(w, zero, kw);
    CHECK_THAT(tau.x(), WithinRel(-1.5 * 0.2, 1.0e-15));
    CHECK_THAT(tau.y(), WithinRel(-2.5 * (-0.1), 1.0e-15));
    CHECK_THAT(tau.z(), WithinRel(-3.5 * 0.15, 1.0e-15));
    // Nonzero command tracks the difference.
    const math::Vector3 cmd{0.05, 0.05, 0.05};
    const math::Vector3 tau_cmd = control::rate_damping_torque_body_Nm(w, cmd, kw);
    CHECK_THAT(tau_cmd.x(), WithinRel(-1.5 * 0.15, 1.0e-15));
    // Zero error: zero torque.
    CHECK_THAT(control::rate_damping_torque_body_Nm(cmd, cmd, kw).norm(), WithinAbs(0.0, 1.0e-15));
    // Decelerating torque opposes the spin on every axis (energy removal).
    CHECK(tau.x() * w.x() < 0.0);
    CHECK(tau.y() * w.y() < 0.0);
    CHECK(tau.z() * w.z() < 0.0);
}

TEST_CASE("M16B PD with zero attitude error reduces to rate damping", "[control][m16b]") {
    const math::Quaternion q = math::Quaternion::identity();
    const math::Vector3 w{0.1, 0.2, -0.15};
    const math::Vector3 w_des{0.01, -0.01, 0.02};
    const control::AttitudePdGains gains{{5.0, 5.0, 5.0}, {1.5, 2.5, 3.5}};
    const math::Vector3 tau_pd = control::attitude_pd_torque_body_Nm(q, q, w, w_des, gains);
    const math::Vector3 tau_rd =
        control::rate_damping_torque_body_Nm(w, w_des, gains.kd_Nm_s_per_rad);
    CHECK_THAT((tau_pd - tau_rd).norm(), WithinAbs(0.0, 1.0e-15));
}

TEST_CASE("M16C CW feedforward cancels the linearized natural dynamics", "[control][m16c]") {
    const double n = 1.1e-3;
    const math::Vector3 rho{10.0, -20.0, 5.0};
    const math::Vector3 vel{1.0, 2.0, -0.5};
    const math::Vector3 ff = control::cw_feedforward_accel_lvlh_mps2(rho, vel, n);
    CHECK_THAT(ff.x(), WithinRel(3.0 * n * n * 10.0 + 2.0 * n * 2.0, 1.0e-12));
    CHECK_THAT(ff.y(), WithinRel(-2.0 * n * 1.0, 1.0e-12));
    CHECK_THAT(ff.z(), WithinRel(-n * n * 5.0, 1.0e-12));
    // At the setpoint with matching velocity, PD terms vanish and the command
    // equals a_des minus feedforward: exact cancellation of natural dynamics.
    const control::RelativePdGains gains{{1.0e-3, 1.0e-3, 1.0e-3}, {0.05, 0.05, 0.05}};
    const math::Vector3 cmd =
        control::relative_pd_accel_lvlh_mps2(rho, vel, rho, vel, math::Vector3{}, n, gains);
    CHECK_THAT((cmd + ff).norm(), WithinAbs(0.0, 1.0e-15));
    // Pure position error drives a restoring acceleration.
    const math::Vector3 cmd_off =
        control::relative_pd_accel_lvlh_mps2(rho, vel, math::Vector3{}, vel, math::Vector3{}, n, gains);
    CHECK(cmd_off.x() < 0.0);  // +10 m radial error -> inward command
    CHECK(cmd_off.z() < 0.0);  // +5 m cross-track error -> downward command
}

TEST_CASE("M16D LQR closed form satisfies the Riccati equation", "[control][m16d]") {
    const control::DoubleIntegratorLqrWeights w{2.0, 3.0, 0.5};
    const control::DoubleIntegratorLqrGain k = control::solve_double_integrator_lqr(w);
    // Hand values: p12 = sqrt(2*0.5) = 1; p22 = sqrt(0.5*(2+3)) = sqrt(2.5).
    CHECK_THAT(k.k_position, WithinRel(1.0 / 0.5, 1.0e-12));
    CHECK_THAT(k.k_velocity, WithinRel(std::sqrt(2.5) / 0.5, 1.0e-12));
    // ARE residuals for P = [[p11,p12],[p12,p22]] with K = B'P/R:
    //   (1,1): -p12^2/R + q_pos = 0; (2,2): 2 p22... check via K directly:
    //   A'P + PA - PBR^-1B'P + Q = 0  <=>  k_pos = p12/R, k_vel = p22/R with
    //   p12^2 = q_pos R and (2 p12 + q_vel) R = p22^2.
    const double p12 = k.k_position * w.r_control;
    const double p22 = k.k_velocity * w.r_control;
    CHECK_THAT(p12 * p12, WithinRel(w.q_position * w.r_control, 1.0e-12));
    CHECK_THAT(p22 * p22, WithinRel((2.0 * p12 + w.q_velocity) * w.r_control, 1.0e-12));
    // Closed loop is Hurwitz: both eigenvalues have negative real parts.
    const auto eig = control::double_integrator_closed_loop_eigenvalues(k);
    CHECK(eig[0] < 0.0);
    CHECK(eig[1] < 0.0);
    // Unit weights give the textbook K = [1, sqrt(3)].
    const control::DoubleIntegratorLqrGain k1 =
        control::solve_double_integrator_lqr(control::DoubleIntegratorLqrWeights{});
    CHECK_THAT(k1.k_position, WithinRel(1.0, 1.0e-12));
    CHECK_THAT(k1.k_velocity, WithinRel(std::sqrt(3.0), 1.0e-12));
    CHECK_THAT(control::double_integrator_damping_ratio(k1), WithinRel(std::sqrt(3.0) / 2.0, 1.0e-12));
    // Control action opposes displacement.
    CHECK(control::apply_double_integrator_lqr(k1, 1.0, 0.0) < 0.0);
}

TEST_CASE("M16E saturation reports desired vs achieved without silent clipping", "[control][m16e]") {
    const math::Vector3 desired{0.5, -0.02, 0.001};
    const math::Vector3 limit{0.1, 0.1, 0.1};
    const control::SaturatedControl sat = control::saturate_control_vector(desired, limit);
    CHECK(sat.saturated);
    CHECK_THAT(sat.achieved.x(), WithinRel(0.1, 1.0e-15));
    CHECK_THAT(sat.achieved.y(), WithinRel(-0.02, 1.0e-15));
    CHECK_THAT(sat.deficit.x(), WithinRel(0.4, 1.0e-12));
    CHECK_THAT((sat.desired - sat.achieved - sat.deficit).norm(), WithinAbs(0.0, 1.0e-15));
    // Within limits: achieved == desired, no flag.
    const control::SaturatedControl free =
        control::saturate_control_vector(math::Vector3{0.01, -0.02, 0.03}, limit);
    CHECK_FALSE(free.saturated);
    CHECK_THAT((free.achieved - free.desired).norm(), WithinAbs(0.0, 1.0e-15));
    CHECK_THAT(free.deficit.norm(), WithinAbs(0.0, 1.0e-15));
    // Negative side clamps symmetrically.
    const control::SaturatedControl neg =
        control::saturate_control_vector(math::Vector3{0.0, -0.5, 0.0}, limit);
    CHECK_THAT(neg.achieved.y(), WithinRel(-0.1, 1.0e-15));
}

TEST_CASE("M16 attitude PD detumbles the 6-DOF bus within the settling budget", "[control][metrics]") {
    // Tumble from the master spec: omega = [0.2, -0.1, 0.15] rad/s.
    attitude::PrincipalInertia inertia{10.0, 20.0, 30.0};
    attitude::RotationalState rot{
        math::Quaternion::identity(), math::Vector3{0.2, -0.1, 0.15}};
    const control::AttitudePdGains gains =
        control::suggest_pd_gains(control::PdTuningRequest{{10.0, 20.0, 30.0}, 0.9, 30.0});
    const math::Quaternion q_hold = math::Quaternion::identity();
    const math::Vector3 w_hold{};
    const double dt = 0.05;
    double t = 0.0;
    double settle_time = -1.0;
    double peak_rate = rot.angular_velocity_rad_per_s.norm();
    double peak_torque = 0.0;
    math::Vector3 torque{};
    for (int i = 0; i < 2000; ++i) {
        torque = control::attitude_pd_torque_body_Nm(
            rot.orientation, q_hold, rot.angular_velocity_rad_per_s, w_hold, gains);
        peak_torque = std::max(peak_torque, torque.norm());
        rot = attitude::rk4_step_rotational(t, rot, dt, inertia, torque);
        t += dt;
        peak_rate = std::max(peak_rate, rot.angular_velocity_rad_per_s.norm());
        const double err = control::attitude_error_angle_rad(rot.orientation, q_hold);
        if (err < 0.01 && rot.angular_velocity_rad_per_s.norm() < 0.005 && settle_time < 0.0) {
            settle_time = t;
        }
    }
    CHECK(peak_rate <= 0.27);  // bounded transient, no violent excursion
    CHECK(settle_time > 0.0);  // settled inside the 100 s window
    CHECK(settle_time < 60.0);
    CHECK(control::attitude_error_angle_rad(rot.orientation, q_hold) < 0.005);
    CHECK(rot.angular_velocity_rad_per_s.norm() < 0.002);
    CHECK(peak_torque < 2.0);  // well inside a 0.05 Nm-class wheel only if wheels sized so;
                               // documents why M14 torque limits force M17 allocation care
}

TEST_CASE("M16 translation PD closes a 100 m radial offset under CW dynamics", "[control][metrics]") {
    const double n = 1.1e-3;
    const control::RelativePdGains gains{{1.0e-5, 1.0e-5, 1.0e-5}, {6.0e-3, 6.0e-3, 6.0e-3}};
    math::Vector3 rho{100.0, 0.0, 0.0};
    math::Vector3 vel{};
    const double dt = 1.0;
    double max_accel = 0.0;
    for (int i = 0; i < 3600; ++i) {
        const math::Vector3 a_cmd = control::relative_pd_accel_lvlh_mps2(
            rho, vel, math::Vector3{}, math::Vector3{}, math::Vector3{}, n, gains);
        max_accel = std::max(max_accel, a_cmd.norm());
        // Exact CW discrete update over dt with piecewise-constant command.
        const relative::RelativeStateLvlh s{rho, vel};
        const auto unforced = relative::cw_predict(s, n, dt);
        // Forced correction via double-integral of constant accel rotated by CW is
        // second order; first-order Euler on the command is sufficient for a
        // controller-response (not truth-model) check at dt = 1 s.
        vel = unforced.relative_velocity_lvlh_mps + a_cmd * dt;
        rho = unforced.relative_position_lvlh_m + a_cmd * (0.5 * dt * dt);
    }
    CHECK(rho.norm() < 5.0);
    CHECK(vel.norm() < 0.05);
    CHECK(max_accel < 0.01);
}

TEST_CASE("M16 controllers reject invalid inputs", "[control][invalid]") {
    const math::Quaternion q = math::Quaternion::identity();
    const math::Vector3 w{};
    const control::AttitudePdGains gains{{1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}};
    CHECK_THROWS_AS(
        control::attitude_pd_torque_body_Nm(
            q, q, w, w, control::AttitudePdGains{{0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}}),
        std::domain_error);
    CHECK_THROWS_AS(
        control::rate_damping_torque_body_Nm(w, w, math::Vector3{1.0, 0.0, -1.0}),
        std::domain_error);
    CHECK_THROWS_AS(control::solve_double_integrator_lqr({1.0, 1.0, 0.0}), std::domain_error);
    CHECK_THROWS_AS(
        control::relative_pd_accel_lvlh_mps2(w, w, w, w, w, -1.0e-3, control::RelativePdGains{}),
        std::domain_error);
    CHECK_THROWS_AS(
        control::saturate_control_vector(w, math::Vector3{0.0, 0.1, 0.1}), std::domain_error);
    CHECK_THROWS_AS(
        control::attitude_error_vector(
            math::Quaternion{0.0, 0.0, 0.0, 0.0}, math::Quaternion::identity()),
        std::domain_error);
}
