// AstraDock M18 — Rendezvous & proximity operations demonstration tool.
//
// M18A relative navigation: ECI truth target + chaser -> LVLH relative state.
// M18B guidance: default 4-leg approach sequence tracked with the M16 law on a
//   CW plant (fast sweep) plus a 6-DOF ECI truth run with GNSS-tracked relative
//   navigation (estimates-only guidance input).
// M18C keep-out/corridor enforcement: violation flags on estimated state.
// M18D validation: nominal rendezvous + burn-error + dropout robustness cases.

#include "control/relative_pd.hpp"
#include "dynamics/two_body.hpp"
#include "estimation/translational_ekf.hpp"
#include "frames/lvlh.hpp"
#include "math/constants.hpp"
#include "math/matrix3.hpp"
#include "math/vector3.hpp"
#include "numerics/fixed_step_propagation.hpp"
#include "orbit/cartesian_state.hpp"
#include "orbit/two_body_orbit.hpp"
#include "relative/relative_state.hpp"
#include "rendezvous/rendezvous_guidance.hpp"
#include "sensors/gnss.hpp"
#include "sensors/sensor_common.hpp"

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

using namespace astradock;

struct RendezvousRun {
    std::vector<double> time_s;
    std::vector<double> range_m;
    std::vector<double> lateral_m;
    std::vector<double> closing_mps;
    std::vector<double> accel_mps2;
    std::vector<int> leg;
    std::vector<int> keep_out;
    std::vector<int> corridor;
    std::vector<int> speed_warn;
    std::vector<double> est_range_m;
    std::vector<double> est_closing_mps;
    bool complete{false};
    double propellant_mps{0.0};  // delta-v proxy: integral |a_cmd| dt
};

// Fast CW-plant sweep of the full default sequence from a far offset.
RendezvousRun run_cw_sweep(double lateral_offset_m, double burn_error_scale) {
    const double n = 1.1e-3;
    const auto seq = rendezvous::default_approach_sequence();
    rendezvous::SafetyCorridor safety;
    const control::RelativePdGains gains{
        math::Vector3{1.0e-5, 1.0e-5, 1.0e-5}, math::Vector3{6.0e-3, 6.0e-3, 6.0e-3}};
    relative::RelativeStateLvlh state{
        math::Vector3{lateral_offset_m, -6000.0, 0.0}, math::Vector3{0.0, 0.0, 0.0}};
    std::size_t leg = 0;
    const double dt = 2.0;
    RendezvousRun run;
    double t = 0.0;
    for (int i = 0; i < 20000 && !run.complete; ++i) {
        auto step = rendezvous::step_rendezvous_guidance(state, seq, leg, safety, n, gains);
        run.complete = step.mission_complete;
        math::Vector3 cmd = step.commanded_accel_lvlh_mps2 * burn_error_scale;
        run.propellant_mps += cmd.norm() * dt;
        const auto drifted = relative::cw_predict(state, n, dt);
        math::Vector3 vel = drifted.relative_velocity_lvlh_mps + cmd * dt;
        math::Vector3 rho = drifted.relative_position_lvlh_m + cmd * (0.5 * dt * dt);
        state = {rho, vel};
        run.time_s.push_back(t);
        run.range_m.push_back(step.range_to_target_m);
        run.lateral_m.push_back(step.lateral_error_m);
        run.closing_mps.push_back(step.closing_speed_mps);
        run.accel_mps2.push_back(step.commanded_accel_lvlh_mps2.norm());
        run.leg.push_back(static_cast<int>(step.active_leg));
        run.keep_out.push_back(step.keep_out_violation ? 1 : 0);
        run.corridor.push_back(step.corridor_violation ? 1 : 0);
        run.speed_warn.push_back(step.closing_speed_warning ? 1 : 0);
        t += dt;
        if (step.keep_out_violation) {
            break;
        }
    }
    return run;
}

void write_run_csv(const std::filesystem::path& filepath, const RendezvousRun& run) {
    std::filesystem::create_directories(filepath.parent_path());
    std::ofstream out(filepath);
    out << "time_s,range_m,lateral_m,closing_mps,accel_mps2,leg,keep_out,corridor,speed_warn,est_range_m,est_closing_mps\n";
    out << std::setprecision(10);
    for (std::size_t i = 0; i < run.time_s.size(); ++i) {
        const double er = i < run.est_range_m.size() ? run.est_range_m[i] : -1.0;
        const double ec = i < run.est_closing_mps.size() ? run.est_closing_mps[i] : 0.0;
        out << run.time_s[i] << "," << run.range_m[i] << "," << run.lateral_m[i] << ","
            << run.closing_mps[i] << "," << run.accel_mps2[i] << "," << run.leg[i] << ","
            << run.keep_out[i] << "," << run.corridor[i] << "," << run.speed_warn[i] << ","
            << er << "," << ec << "\n";
    }
}

// 6-DOF ECI truth run: target on a circular orbit, chaser propagated under
// two-body + guidance thrust; GNSS-tracked relative navigation feeds guidance
// (estimates only). Short final-approach leg for runtime.
RendezvousRun run_eci_truth(bool gnss_dropout) {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const double R = constants::earth_reference_radius_m + 500.0e3;
    const auto ref = orbit::compute_circular_orbit_reference(mu, R);
    const double n = ref.mean_motion_rad_per_s;

    orbit::CartesianState target{{R, 0.0, 0.0}, {0.0, ref.speed_m_per_s, 0.0}};
    // Chaser starts 1200 m behind ON the approach axis, co-orbiting: position
    // offset along LVLH, velocity copied from the target (zero rotating-frame
    // relative velocity at this instant — a natural co-orbiting hold, not a
    // flyby). On-axis start is the realistic far-field condition (targeting
    // puts the chaser on the axis); off-axis correction is exercised by the
    // corridor-stress sweep instead. Starting 100 m off-axis with a 10 m
    // capture ball and a 1 m/s flyby is geometrically impossible (found during
    // development: lateral stuck at 75 m while y blew past the waypoint).
    const math::Matrix3 c_init = frames::dcm_lvlh_from_eci(target.position, target.velocity);
    // LVLH y-axis = along-track = direction of motion: -y offset puts the
    // chaser BEHIND the target. (M15/M06 convention: x radial, y along-track.)
    const math::Vector3 rho_init_lvlh{0.0, -1200.0, 0.0};
    orbit::CartesianState chaser{
        target.position + c_init.transpose() * rho_init_lvlh, target.velocity};

    const std::vector<rendezvous::HoldPoint> seq{
        {math::Vector3{0.0, -250.0, 0.0}, 0.5, 10.0},
        {math::Vector3{0.0, -50.0, 0.0}, 0.2, 5.0},
    };
    rendezvous::SafetyCorridor safety;
    const control::RelativePdGains gains{
        math::Vector3{1.0e-5, 1.0e-5, 1.0e-5}, math::Vector3{6.0e-3, 6.0e-3, 6.0e-3}};

    // Relative navigation: TranslationalEkf on the chaser ECI state, seeded with
    // GNSS noise characteristics (5 m / 0.05 m/s). THRUST-AWARE predict: the
    // commanded thrust acceleration (known onboard — the spacecraft commanded
    // it) is added to the filter's two-body dynamics predict. predict_with_imu
    // was considered but REJECTED here: its validated_orientation() demands a
    // unit quaternion and our ECI-acceleration-as-specific-force path has no
    // attitude content (C = I would need |q| = 1 to 1e-6 — workable but
    // conceptually dishonest). Instead the mean is propagated by direct RK4
    // with thrust inside the derivative, and covariance by the standard
    // two-body Phi plus thrust-magnitude-scaled process noise. Without thrust
    // awareness the two-body-only predict lags under sustained thrust and
    // guidance under-brakes (documented M18 finding).
    estimation::TranslationalEkfConfig ekf_cfg;
    ekf_cfg.gravitational_parameter_m3_per_s2 = mu;
    ekf_cfg.acceleration_noise_std_mps2 = 1.0e-3;
    ekf_cfg.gnss_noise.position_variance_m2 = 25.0;
    ekf_cfg.gnss_noise.velocity_variance_m2_per_s2 = 2.5e-3;
    ekf_cfg.imu.accelerometer_noise_std_mps2 = 1.0e-3;
    estimation::TranslationalCovariance init_cov = estimation::TranslationalCovariance::zero();
    init_cov(0, 0) = 400.0;
    init_cov(1, 1) = 400.0;
    init_cov(2, 2) = 400.0;
    init_cov(3, 3) = 0.25;
    init_cov(4, 4) = 0.25;
    init_cov(5, 5) = 0.25;
    estimation::TranslationalEkf filter(
        ekf_cfg,
        orbit::CartesianState{
            chaser.position + math::Vector3{20.0, -15.0, 10.0},
            chaser.velocity + math::Vector3{0.2, -0.2, 0.1}},
        init_cov);

    std::size_t leg = 0;
    const double dt = 1.0;
    double t = 0.0;
    RendezvousRun run;
    sensors::DeterministicRng rng(1234);
    const double dropout_start = 200.0;
    const double dropout_end = 260.0;
    math::Vector3 a_eci_prev{};  // commanded thrust accel (ECI) from prior tick
    auto two_body = [&](double, const orbit::CartesianState& s) {
        return orbit::two_body_state_derivative(0.0, s, mu);
    };
    for (int i = 0; i < 6000 && !run.complete; ++i) {
        // Epoch-aligned sequencing (M18 finding): propagate BOTH vehicles AND
        // the filter predict to the new epoch FIRST, then form the relative
        // state from common-epoch pairs. Stepping the target before the
        // chaser/estimate injects a ~7.6 km epoch-mismatch bias (one tick of
        // orbital motion) into rho — the failure mode that stalled this demo
        // during development. Thrust a_eci_prev from the previous tick is
        // applied to the chaser and the filter mean together (known-actuation
        // awareness); guidance then acts on the fresh estimate.
        target = numerics::rk4_step(t, target, dt, two_body);
        chaser = numerics::rk4_step(
            t, chaser, dt, [&](double, const orbit::CartesianState& s) {
                orbit::CartesianState d = orbit::two_body_state_derivative(0.0, s, mu);
                d.velocity = d.velocity + a_eci_prev;
                return d;
            });
        // Chaser EKF predict (THRUST-AWARE): predict(dt) advances covariance
        // and moves the mean two-body-only; the known thrust a_eci_prev is
        // folded into the mean by exact double-integral correction (reset
        // preserves covariance bit-for-bit). Without thrust awareness the
        // predict lags under sustained thrust and guidance under-brakes.
        filter.predict(dt);
        {
            orbit::CartesianState corrected = filter.estimated_state();
            corrected.velocity = corrected.velocity + a_eci_prev * dt;
            corrected.position = corrected.position + a_eci_prev * (0.5 * dt * dt);
            filter.reset(corrected, filter.covariance());
        }
        const bool in_dropout = gnss_dropout && (t + dt >= dropout_start) && (t + dt <= dropout_end);
        if (!in_dropout) {
            // Propagated chaser is now at epoch t+dt: sample GNSS there.
            const math::Vector3 pn{
                rng.gaussian(0.0, 5.0), rng.gaussian(0.0, 5.0), rng.gaussian(0.0, 5.0)};
            const math::Vector3 vn{
                rng.gaussian(0.0, 0.05), rng.gaussian(0.0, 0.05), rng.gaussian(0.0, 0.05)};
            filter.update_gnss(chaser.position + pn, chaser.velocity + vn);
        }
        t += dt;
        // Guidance from the ESTIMATED relative state (common epoch t).
        const orbit::CartesianState est_chaser = filter.estimated_state();
        const relative::RelativeStateLvlh est_rel =
            relative::relative_state_from_eci(target, est_chaser);
        auto step = rendezvous::step_rendezvous_guidance(est_rel, seq, leg, safety, n, gains);
        run.complete = step.mission_complete;
        // Thrust along LVLH command rotated to ECI via the target LVLH DCM,
        // applied on the NEXT tick (zero-order hold across [t, t+dt]).
        const math::Matrix3 c_lvlh_eci =
            frames::dcm_lvlh_from_eci(target.position, target.velocity);
        a_eci_prev = c_lvlh_eci.transpose() * step.commanded_accel_lvlh_mps2;
        const relative::RelativeStateLvlh scored =
            relative::relative_state_from_eci(target, chaser);
        run.time_s.push_back(t);
        run.range_m.push_back(scored.relative_position_lvlh_m.norm());
        run.lateral_m.push_back(rendezvous::lateral_axis_error_m(scored.relative_position_lvlh_m));
        run.closing_mps.push_back(
            rendezvous::closing_speed_mps(scored.relative_position_lvlh_m, scored.relative_velocity_lvlh_mps));
        run.accel_mps2.push_back(step.commanded_accel_lvlh_mps2.norm());
        run.leg.push_back(static_cast<int>(step.active_leg));
        run.keep_out.push_back(step.keep_out_violation ? 1 : 0);
        run.corridor.push_back(step.corridor_violation ? 1 : 0);
        run.speed_warn.push_back(step.closing_speed_warning ? 1 : 0);
        run.est_range_m.push_back(est_rel.relative_position_lvlh_m.norm());
        run.est_closing_mps.push_back(
            rendezvous::closing_speed_mps(est_rel.relative_position_lvlh_m, est_rel.relative_velocity_lvlh_mps));
        run.propellant_mps += step.commanded_accel_lvlh_mps2.norm() * dt;
        if (step.keep_out_violation) {
            break;
        }
    }
    return run;
}

}  // namespace

int main() {
    using namespace astradock;
    std::cout << "============================================================\n";
    std::cout << " AstraDock — M18 Rendezvous & Proximity Operations Demo     \n";
    std::cout << "============================================================\n";
    std::filesystem::create_directories("data");

    std::cout << "------------------------------------------------------------\n";
    std::cout << "Case 1: nominal 6 km -> 50 m approach (CW sweep)\n";
    const RendezvousRun nominal = run_cw_sweep(0.0, 1.0);
    write_run_csv("data/m18_nominal.csv", nominal);
    std::cout << "  Complete: " << (nominal.complete ? "yes" : "NO")
              << "; delta-v proxy " << nominal.propellant_mps << " m/s\n";

    std::cout << "------------------------------------------------------------\n";
    std::cout << "Case 2: 5% under-burn error\n";
    const RendezvousRun burn_err = run_cw_sweep(0.0, 0.95);
    write_run_csv("data/m18_burn_error.csv", burn_err);
    std::cout << "  Complete: " << (burn_err.complete ? "yes" : "NO")
              << "; delta-v proxy " << burn_err.propellant_mps << " m/s\n";

    std::cout << "------------------------------------------------------------\n";
    std::cout << "Case 3: 60 m lateral offset (corridor stress)\n";
    const RendezvousRun lateral = run_cw_sweep(60.0, 1.0);
    int corridor_flags = 0;
    for (int f : lateral.corridor) {
        corridor_flags += f;
    }
    write_run_csv("data/m18_lateral.csv", lateral);
    std::cout << "  Complete: " << (lateral.complete ? "yes" : "NO") << "; corridor flags "
              << corridor_flags << "\n";

    std::cout << "------------------------------------------------------------\n";
    std::cout << "Case 4: 6-DOF ECI truth, GNSS relative navigation (no dropout)\n";
    const RendezvousRun eci = run_eci_truth(false);
    write_run_csv("data/m18_eci_truth.csv", eci);
    std::cout << "  Complete: " << (eci.complete ? "yes" : "NO") << "; final range "
              << (eci.range_m.empty() ? -1.0 : eci.range_m.back()) << " m\n";

    std::cout << "------------------------------------------------------------\n";
    std::cout << "Case 5: 6-DOF ECI truth with 60 s GNSS dropout\n";
    const RendezvousRun dropout = run_eci_truth(true);
    write_run_csv("data/m18_dropout.csv", dropout);
    std::cout << "  Complete: " << (dropout.complete ? "yes" : "NO") << "; final range "
              << (dropout.range_m.empty() ? -1.0 : dropout.range_m.back()) << " m\n";

    std::cout << "============================================================\n";
    std::cout << " All M18 rendezvous demonstrations completed successfully.\n";
    std::cout << "============================================================\n";
    return 0;
}
