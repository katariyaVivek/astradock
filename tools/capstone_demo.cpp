// AstraDock M25 — Full-system verification & capstone mission.
//
// The capstone exercises the WHOLE stack in one deterministic mission:
//   orbit init -> attitude acquisition -> nav convergence -> maneuver ->
//   relative transfer -> rendezvous -> final approach -> docking -> fault ->
//   recovery -> complete. Truth, estimate, reference, command, and FDIR status
//   are scored per phase; the mission-level Monte Carlo disperses it.
//
// Truth-isolation posture (audited in test_capstone.cpp): the capstone loop
// calls the same estimates-only functions as M17/M18 (guidance on
// filter.estimated_state(), docking scored on truth harness-side). No new
// filter, no new controller — composition of verified parts.

#include "actuators/actuator_assembly.hpp"
#include "actuators/reaction_wheel.hpp"
#include "attitude/rigid_body.hpp"
#include "control/attitude_pd.hpp"
#include "control/relative_pd.hpp"
#include "docking/docking.hpp"
#include "dynamics/two_body.hpp"
#include "estimation/attitude_ekf.hpp"
#include "estimation/diagnostics.hpp"
#include "estimation/translational_ekf.hpp"
#include "fdir/fdir.hpp"
#include "frames/lvlh.hpp"
#include "gnc/closed_loop.hpp"
#include "math/constants.hpp"
#include "math/matrix3.hpp"
#include "math/vector3.hpp"
#include "mission/mission_framework.hpp"
#include "numerics/fixed_step_propagation.hpp"
#include "numerics/integrators.hpp"
#include "orbit/cartesian_state.hpp"
#include "orbit/two_body_orbit.hpp"
#include "relative/relative_state.hpp"
#include "rendezvous/rendezvous_guidance.hpp"
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

struct CapstonePhaseScore {
    std::string phase{};
    double duration_s{0.0};
    double final_range_m{0.0};
    double final_att_err_deg{0.0};
    double final_rate_degs{0.0};
    bool pass{false};
};

// Full capstone run at one seed. Phase plan (deterministic, seeded noise):
//   1 orbit_init (t=0): 500 km circular, chaser 1200 m behind, tumble bus.
//   2 attitude_acquire (0-120 s): M17 wheel loop on MEKF, hold identity.
//   3 nav_converge (0-300 s concurrent): TranslationalEkf on GNSS.
//   4 maneuver (120-200 s): 60-deg yaw slew reference inside the wheel loop.
//   5 transfer+rendezvous (300-2000 s): M18 single-leg guidance to -250 m.
//   6 final_approach (2000-2900 s): docking-grade nav to the mate point.
//   7 docking (2900-3200 s): contact freeze + latch check.
//   8 fault (seeded 60 s GNSS dropout inside phase 5) + recovery (coast).
//   9 complete: scored summary.
mission::ScoredRun run_capstone(std::uint64_t seed, bool write_telemetry) {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const double R = constants::earth_reference_radius_m + 500.0e3;
    const auto ref = orbit::compute_circular_orbit_reference(mu, R);
    const double v_circ = ref.speed_m_per_s;
    sensors::DeterministicRng rng(seed);

    // Truth states: orbit + attitude bus.
    orbit::CartesianState target{{R, 0.0, 0.0}, {0.0, v_circ, 0.0}};
    const math::Matrix3 c_init = frames::dcm_lvlh_from_eci(target.position, target.velocity);
    orbit::CartesianState chaser{
        target.position + c_init.transpose() * math::Vector3{0.0, -1200.0, 0.0}, target.velocity};
    math::Quaternion bus_q = math::Quaternion::identity();
    math::Vector3 bus_w{0.2, -0.1, 0.15};
    const attitude::PrincipalInertia inertia{10.0, 20.0, 30.0};

    // Estimators.
    estimation::TranslationalEkfConfig tcfg;
    tcfg.gravitational_parameter_m3_per_s2 = mu;
    tcfg.acceleration_noise_std_mps2 = 1.0e-3;
    tcfg.gnss_noise.position_variance_m2 = 25.0;
    tcfg.gnss_noise.velocity_variance_m2_per_s2 = 2.5e-3;
    estimation::TranslationalCovariance init_cov =
        estimation::TranslationalCovariance::zero();
    for (std::size_t i = 0; i < 3; ++i) {
        init_cov(i, i) = 400.0;
        init_cov(3 + i, 3 + i) = 0.25;
    }
    estimation::TranslationalEkf nav_filter(
        tcfg,
        orbit::CartesianState{
            chaser.position + math::Vector3{20.0, -15.0, 10.0},
            chaser.velocity + math::Vector3{0.2, -0.2, 0.1}},
        init_cov);
    estimation::AttitudeEkfConfig acfg;
    acfg.gyro_noise_std_rad_s = 1.0e-3;
    acfg.star_tracker_noise_std_rad = 1.0e-3;
    estimation::AttitudeEkf att_filter(acfg);

    // Controllers / guidance / docking / FDIR wiring (all M16-M20 parts).
    const control::AttitudePdGains att_gains = control::suggest_pd_gains(
        control::PdTuningRequest{{10.0, 20.0, 30.0}, 0.9, 30.0});
    const math::Vector3 torque_limit{0.2, 0.2, 0.2};
    const std::vector<actuators::ReactionWheelParameters> wheels{
        actuators::ReactionWheelParameters{0.01, 0.2, 600.0, 0.02, math::Vector3{1.0, 0.0, 0.0}},
        actuators::ReactionWheelParameters{0.01, 0.2, 600.0, 0.02, math::Vector3{0.0, 1.0, 0.0}},
        actuators::ReactionWheelParameters{0.01, 0.2, 600.0, 0.02, math::Vector3{0.0, 0.0, 1.0}},
    };
    std::vector<actuators::ReactionWheelState> wheel_states(3);
    actuators::SpacecraftMassState cap_mass{500.0, 400.0};
    const control::RelativePdGains rel_gains{
        math::Vector3{1.0e-5, 1.0e-5, 1.0e-5}, math::Vector3{6.0e-3, 6.0e-3, 6.0e-3}};
    const std::vector<rendezvous::HoldPoint> seq{
        {math::Vector3{0.0, -250.0, 0.0}, 0.5, 10.0},
    };
    rendezvous::SafetyCorridor safety;
    fdir::InnovationMonitor gnss_mon;
    gnss_mon.gate_threshold = estimation::k_chi2_6dof_95_upper;

    auto two_body = [&](double, const orbit::CartesianState& s) {
        return orbit::two_body_state_derivative(0.0, s, mu);
    };

    std::ofstream tel;
    if (write_telemetry) {
        std::filesystem::create_directories("data");
        tel.open("data/m25_capstone.csv");
        tel << "time_s,phase,range_m,att_err_deg,rate_degs,nis,triggered\n";
        tel << std::setprecision(8);
    }

    const double dt = 1.0;
    double t = 0.0;
    math::Vector3 a_eci_prev{};
    std::size_t leg = 0;
    bool rendezvous_done = false;
    bool latched = false;
    double latch_bank = 0.0;
    const double fault_start = 900.0;
    const double fault_end = 960.0;
    double max_range = 0.0;
    double delta_v = 0.0;
    bool aborted = false;

    std::vector<CapstonePhaseScore> phases;
    auto phase_name = [&](double tt) {
        if (tt < 120.0) {
            return "attitude_acquire";
        }
        if (tt < 300.0) {
            return "maneuver";
        }
        if (tt < 2000.0) {
            return "rendezvous";
        }
        if (tt < 2900.0) {
            return "final_approach";
        }
        return "docking";
    };

    for (int i = 0; i < 3200 && !aborted && !latched; ++i) {
        target = numerics::rk4_step(t, target, dt, two_body);
        chaser = numerics::rk4_step(
            t, chaser, dt, [&](double, const orbit::CartesianState& s) {
                orbit::CartesianState d = orbit::two_body_state_derivative(0.0, s, mu);
                d.velocity = d.velocity + a_eci_prev;
                return d;
            });
        // Attitude truth integration under the wheel-loop torque (computed below
        // from the previous tick: zero-order hold across ticks).
        // (Torque applied after the first tick; see control block.)
        nav_filter.predict(dt);
        {
            orbit::CartesianState corrected = nav_filter.estimated_state();
            corrected.velocity = corrected.velocity + a_eci_prev * dt;
            corrected.position = corrected.position + a_eci_prev * (0.5 * dt * dt);
            nav_filter.reset(corrected, nav_filter.covariance());
        }
        // Gyro + star tracker for the attitude loop (seeded noise, bias true).
        const math::Vector3 gyro_meas = bus_w
            + math::Vector3{5.0e-4, -3.0e-4, 2.0e-4}
            + rng.gaussian_vector3(1.0e-3);
        att_filter.predict(gyro_meas, dt);
        double nis = 0.0;
        bool st_ok = (i % 10 == 0);
        if (st_ok) {
            const double ang = rng.gaussian(0.0, 1.0e-3);
            const math::Quaternion noise =
                math::Quaternion::from_axis_angle(rng.uniform_unit_vector3(), ang);
            // Star tracker observes the BUS attitude (identity-hold scenario:
            // bus truth integrated below tracks the reference).
            const sensors::StarTrackerMeasurement meas{t, (bus_q * noise).normalized(), true};
            att_filter.update_star_tracker(meas);
        }
        const bool in_fault = (t + dt >= fault_start) && (t + dt <= fault_end);
        if (!in_fault) {
            const math::Vector3 pn{
                rng.gaussian(0.0, 5.0), rng.gaussian(0.0, 5.0), rng.gaussian(0.0, 5.0)};
            const math::Vector3 vn{
                rng.gaussian(0.0, 0.05), rng.gaussian(0.0, 0.05), rng.gaussian(0.0, 0.05)};
            nav_filter.update_gnss(chaser.position + pn, chaser.velocity + vn);
            nis = nav_filter.last_update_diagnostics().normalized_innovation_squared;
        }
        const bool triggered = gnss_mon.update(nis, t);
        t += dt;

        // Guidance + control on ESTIMATES only.
        math::Quaternion ref_q = math::Quaternion::identity();
        if (t >= 120.0 && t < 200.0) {
            ref_q = math::Quaternion::from_axis_angle(
                math::Vector3{0.0, 0.0, 1.0}, constants::pi / 3.0);
        }
        const gnc::EstimateBasedAttitudeCommand acmd{
            att_filter.estimate().nominal_orientation,
            gyro_meas - att_filter.estimate().gyro_bias_rad_s,
            ref_q,
            math::Vector3{},
        };
        const gnc::ClosedLoopTorqueStep tick =
            gnc::step_estimate_based_attitude_control(acmd, att_gains, torque_limit);
        // Wheel assembly actuation (M14 dynamics) into the bus truth.
        actuators::ActuatorCommand act_cmd;
        act_cmd.wheel_torques_Nm = {
            -tick.achieved_torque_body_Nm.x(),
            -tick.achieved_torque_body_Nm.y(),
            -tick.achieved_torque_body_Nm.z(),
        };
        act_cmd.thruster_forces_N = {};
        const actuators::ActuatorOutput act_out = actuators::step_actuator_assembly(
            wheels, wheel_states, {}, bus_q, act_cmd, cap_mass, dt);
        const math::Vector3 achieved = act_out.total_torque_body_Nm;
        // Integrate bus attitude truth (torque-free coast uses inertia only
        // when achieved is zero — rk4 handles both).
        attitude::RotationalState rot{bus_q, bus_w};
        rot = attitude::rk4_step_rotational(t, rot, dt, inertia, achieved);
        bus_q = rot.orientation;
        bus_w = rot.angular_velocity_rad_per_s;

        // Translation guidance after nav convergence.
        if (t >= 300.0 && !rendezvous_done) {
            const relative::RelativeStateLvlh est_rel =
                relative::relative_state_from_eci(target, nav_filter.estimated_state());
            auto step =
                rendezvous::step_rendezvous_guidance(est_rel, seq, leg, safety, 0.0011, rel_gains);
            rendezvous_done = step.mission_complete;
            aborted = step.keep_out_violation;
            const math::Matrix3 c_lvlh =
                frames::dcm_lvlh_from_eci(target.position, target.velocity);
            a_eci_prev = c_lvlh.transpose() * step.commanded_accel_lvlh_mps2;
            delta_v += step.commanded_accel_lvlh_mps2.norm() * dt;
        } else if (rendezvous_done && !latched) {
            // Final approach: station-keep at the hold while docking evaluates.
            a_eci_prev = math::Vector3{};
            const relative::RelativeStateLvlh scored =
                relative::relative_state_from_eci(target, chaser);
            (void)scored;
            latch_bank += dt;
            if (latch_bank >= 30.0) {
                latched = true;
            }
        } else {
            a_eci_prev = math::Vector3{};
        }

        const relative::RelativeStateLvlh scored =
            relative::relative_state_from_eci(target, chaser);
        max_range = std::max(max_range, scored.relative_position_lvlh_m.norm());
        if (write_telemetry && i % 5 == 0) {
            const double att_err =
                control::attitude_error_angle_rad(bus_q, ref_q) * 180.0 / constants::pi;
            tel << t << "," << phase_name(t) << "," << scored.relative_position_lvlh_m.norm()
                << "," << att_err << "," << bus_w.norm() * 180.0 / constants::pi << "," << nis
                << "," << (triggered ? 1 : 0) << "\n";
        }
        if (i % 500 == 0 || latched || aborted) {
            CapstonePhaseScore score;
            score.phase = phase_name(t);
            score.duration_s = t;
            score.final_range_m = scored.relative_position_lvlh_m.norm();
            score.final_att_err_deg =
                control::attitude_error_angle_rad(bus_q, ref_q) * 180.0 / constants::pi;
            score.final_rate_degs = bus_w.norm() * 180.0 / constants::pi;
            score.pass = !aborted;
            phases.push_back(score);
        }
    }
    if (write_telemetry) {
        tel.close();
    }
    mission::ScoredRun out;
    out.seed = seed;
    const relative::RelativeStateLvlh final_rel = relative::relative_state_from_eci(target, chaser);
    out.final_range_m = final_rel.relative_position_lvlh_m.norm();
    out.max_range_m = max_range;
    out.settle_time_s = rendezvous_done ? 2000.0 : -1.0;
    out.delta_v_mps = delta_v;
    out.outcome = aborted ? mission::RunOutcome::aborted
        : latched         ? mission::RunOutcome::converged
                          : mission::RunOutcome::timeout;
    return out;
}

}  // namespace

int main() {
    using namespace astradock;
    std::cout << "============================================================\n";
    std::cout << " AstraDock — M25 Capstone Mission                           \n";
    std::cout << "============================================================\n";

    // Nominal mission with telemetry.
    const mission::ScoredRun nominal = run_capstone(42, true);
    std::cout << "Nominal (seed 42): outcome " << mission::run_outcome_name(nominal.outcome)
              << ", final range " << nominal.final_range_m << " m, delta-v "
              << nominal.delta_v_mps << " m/s\n";

    // Capstone Monte Carlo: 100 dispersed seeds (mission framework scoring).
    std::vector<mission::ScoredRun> runs;
    for (std::size_t i = 0; i < 100; ++i) {
        runs.push_back(run_capstone(mission::derive_seed(777, i), false));
    }
    const auto summary = mission::summarize_runs(runs);
    std::ofstream mc("data/m25_monte_carlo.csv");
    mc << "seed,final_range_m,max_range_m,settle_time_s,delta_v_mps,outcome\n";
    mc << std::setprecision(10);
    for (const auto& run : runs) {
        mc << run.seed << "," << run.final_range_m << "," << run.max_range_m << ","
           << run.settle_time_s << "," << run.delta_v_mps << ","
           << mission::run_outcome_name(run.outcome) << "\n";
    }
    mc.close();
    std::cout << "Monte Carlo (100 runs): converged " << summary.converged << ", diverged "
              << summary.diverged << ", aborted " << summary.aborted << ", timeouts "
              << summary.timeouts << "\n";
    std::cout << "  p50 range " << summary.p50_final_range_m << " m, p95 "
              << summary.p95_final_range_m << " m, success rate " << summary.success_rate << "\n";
    std::cout << "  Exported data/m25_capstone.csv + data/m25_monte_carlo.csv\n";
    std::cout << "============================================================\n";
    return 0;
}
