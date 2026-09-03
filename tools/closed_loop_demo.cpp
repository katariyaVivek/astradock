// AstraDock M17 — Integrated Closed-Loop GNC demonstration tool.
//
// Scenario (deterministic, seeded):
//   M17A attitude hold: tumbling bus, gyro + star tracker -> AttitudeEkf ->
//     estimate-based PD -> 3-axis reaction-wheel assembly -> 6-DOF truth.
//   M17B maneuver: reference timeline switches identity -> 60-deg yaw at t=40 s.
//
// Truth isolation: the controller receives ONLY filter estimates + timeline
// reference. Truth is read for sensor generation and telemetry scoring. The
// compile-time shape of EstimateBasedAttitudeCommand (no truth member) plus
// the M17 static audit test enforce this structurally.

#include "actuators/actuator_assembly.hpp"
#include "actuators/reaction_wheel.hpp"
#include "attitude/rigid_body.hpp"
#include "control/attitude_pd.hpp"
#include "estimation/attitude_ekf.hpp"
#include "estimation/diagnostics.hpp"
#include "gnc/closed_loop.hpp"
#include "math/constants.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"
#include "orbit/cartesian_state.hpp"
#include "orbit/two_body_orbit.hpp"
#include "sensors/sensor_common.hpp"
#include "sensors/star_tracker.hpp"
#include "spacecraft/force_torque.hpp"
#include "spacecraft/six_dof_dynamics.hpp"
#include "spacecraft/spacecraft_parameters.hpp"
#include "spacecraft/spacecraft_state.hpp"

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

struct ClosedLoopRun {
    std::vector<double> time_s;
    std::vector<double> truth_error_rad;
    std::vector<double> estimate_error_rad;
    std::vector<double> reference_error_rad;
    std::vector<double> rate_norm_rad_s;
    std::vector<double> desired_torque_Nm;
    std::vector<double> achieved_torque_Nm;
    std::vector<int> saturated;
    std::vector<int> st_valid;
    std::vector<double> nis;
    std::vector<double> wheel_speed_0;
    std::vector<double> wheel_momentum_total;
};

ClosedLoopRun run_closed_loop(
    const std::string& csv_path,
    std::uint64_t seed,
    const gnc::AttitudeReferenceTimeline& timeline,
    double duration_s,
    bool log_header) {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const double r_orbit = constants::earth_reference_radius_m + 500.0e3;
    const double v_circ = std::sqrt(mu / r_orbit);

    // Truth: 500 km circular orbit + master-spec tumble.
    spacecraft::SpacecraftState truth{
        orbit::CartesianState{
            math::Vector3{r_orbit, 0.0, 0.0}, math::Vector3{0.0, v_circ, 0.0}},
        attitude::RotationalState{
            math::Quaternion::identity(), math::Vector3{0.2, -0.1, 0.15}}};
    const spacecraft::SpacecraftParameters params{
        mu, attitude::PrincipalInertia{10.0, 20.0, 30.0}};

    // Sensors observe truth (seeded determinism).
    sensors::StarTrackerConfig st_cfg;
    st_cfg.sample_period_s = 0.1;
    st_cfg.noise_std_rad = 1.0e-3;
    st_cfg.random_seed = seed + 11;
    sensors::StarTrackerSensor star_tracker(st_cfg);
    sensors::DeterministicRng gyro_rng(seed + 77);
    const double gyro_noise_std = 1.0e-3;
    const math::Vector3 true_gyro_bias{5.0e-4, -3.0e-4, 2.0e-4};

    // Navigation: MEKF on gyro + star tracker (estimates only from here on).
    estimation::AttitudeEkfConfig ekf_cfg;
    ekf_cfg.gyro_noise_std_rad_s = gyro_noise_std;
    ekf_cfg.star_tracker_noise_std_rad = 1.0e-3;
    estimation::AttitudeEkf filter(ekf_cfg);

    // Control: M16-tuned PD + 0.2 Nm per-axis demo cap.
    const control::AttitudePdGains gains = control::suggest_pd_gains(
        control::PdTuningRequest{{10.0, 20.0, 30.0}, 0.9, 30.0});
    const math::Vector3 torque_limit{0.2, 0.2, 0.2};

    // Actuation: 3 orthogonal wheels, 0.2 Nm / 6000 rpm-class envelope.
    const std::vector<actuators::ReactionWheelParameters> wheels{
        actuators::ReactionWheelParameters{0.01, 0.2, 600.0, 0.02, math::Vector3{1.0, 0.0, 0.0}},
        actuators::ReactionWheelParameters{0.01, 0.2, 600.0, 0.02, math::Vector3{0.0, 1.0, 0.0}},
        actuators::ReactionWheelParameters{0.01, 0.2, 600.0, 0.02, math::Vector3{0.0, 0.0, 1.0}},
    };
    std::vector<actuators::ReactionWheelState> wheel_states(3);

    const double dt = 0.05;  // 20 Hz GNC tick; gyro sampled each tick
    const double st_period = 0.1;
    double next_st_time = 0.0;

    std::ofstream csv(csv_path);
    if (log_header) {
        csv << "time_s,scored_error_rad,estimate_error_rad,reference_error_rad,"
               "rate_norm_rad_s,desired_torque_Nm,achieved_torque_Nm,saturated,"
               "st_valid,nis,wheel_speed_0_rad_s,wheel_momentum_total_Nms\n";
    }
    csv << std::setprecision(10);

    ClosedLoopRun run;
    double t = 0.0;
    // Filter gyro-bias seed: start at zero (must converge to true bias online).
    while (t <= duration_s + 1.0e-12) {
        const auto& ref = timeline.reference_at(t);

        // 1. Sensor sampling (truth -> measurement).
        const math::Vector3 gyro_meas = truth.angular_velocity_rad_per_s()
            + true_gyro_bias + gyro_rng.gaussian_vector3(gyro_noise_std);
        filter.predict(gyro_meas, dt);

        bool st_ok = false;
        double nis = 0.0;
        if (t >= next_st_time - 1.0e-9) {
            const auto st_meas = star_tracker.measure(next_st_time, truth);
            if (st_meas.valid) {
                filter.update_star_tracker(st_meas);
                nis = filter.last_star_tracker_diagnostics().normalized_innovation_squared;
                st_ok = true;
            }
            next_st_time += st_period;
        }

        // 2-3. Guidance reference + estimate-based control (NO truth input).
        const gnc::EstimateBasedAttitudeCommand cmd{
            filter.estimate().nominal_orientation,
            gyro_meas - filter.estimate().gyro_bias_rad_s,
            ref.attitude_body_to_eci,
            ref.rate_body_rad_s,
        };
        const gnc::ClosedLoopTorqueStep tick =
            gnc::step_estimate_based_attitude_control(cmd, gains, torque_limit);

        // 4. Actuation: per-axis desired torque -> wheel assembly commands.
        actuators::ActuatorCommand act_cmd;
        act_cmd.wheel_torques_Nm = {
            -tick.achieved_torque_body_Nm.x(),
            -tick.achieved_torque_body_Nm.y(),
            -tick.achieved_torque_body_Nm.z(),
        };
        act_cmd.thruster_forces_N = {};
        actuators::SpacecraftMassState mass{500.0, 400.0};
        const actuators::ActuatorOutput act_out = actuators::step_actuator_assembly(
            wheels, wheel_states, {}, filter.estimate().nominal_orientation, act_cmd, mass, dt);

        // 5. Truth integration under achieved torque.
        truth = spacecraft::rk4_step_spacecraft(
            t, truth, dt, params,
            spacecraft::ForceTorqueInput{math::Vector3{}, act_out.total_torque_body_Nm});

        // 6. Telemetry (truth readable here: scoring only).
        t += dt;
        const double truth_err =
            control::attitude_error_angle_rad(truth.orientation(), ref.attitude_body_to_eci);
        const double est_err = estimation::attitude_error_vector_rad(
                                   truth.orientation(), filter.estimate().nominal_orientation)
                                   .norm();
        const double ref_err = tick.attitude_error_angle_rad;
        run.time_s.push_back(t);
        run.truth_error_rad.push_back(truth_err);
        run.estimate_error_rad.push_back(est_err);
        run.reference_error_rad.push_back(ref_err);
        run.rate_norm_rad_s.push_back(truth.angular_velocity_rad_per_s().norm());
        run.desired_torque_Nm.push_back(tick.desired_torque_body_Nm.norm());
        run.achieved_torque_Nm.push_back(tick.achieved_torque_body_Nm.norm());
        run.saturated.push_back(tick.saturated ? 1 : 0);
        run.st_valid.push_back(st_ok ? 1 : 0);
        run.nis.push_back(nis);
        run.wheel_speed_0.push_back(wheel_states[0].wheel_speed_rad_s);
        run.wheel_momentum_total.push_back(
            actuators::total_wheel_momentum_body_Nms(wheels, wheel_states).norm());
        csv << t << "," << truth_err << "," << est_err << "," << ref_err << ","
            << run.rate_norm_rad_s.back() << "," << run.desired_torque_Nm.back() << ","
            << run.achieved_torque_Nm.back() << "," << run.saturated.back() << ","
            << run.st_valid.back() << "," << nis << "," << run.wheel_speed_0.back() << ","
            << run.wheel_momentum_total.back() << "\n";
    }
    return run;
}

void write_monte_carlo_csv(
    const std::filesystem::path& filepath,
    const std::vector<std::uint64_t>& seeds,
    const std::vector<double>& settle_times,
    const std::vector<double>& final_errors,
    const std::vector<double>& final_rates,
    const std::vector<double>& max_errors,
    const std::vector<double>& sat_fractions,
    const std::vector<double>& rms_est_errors) {
    std::filesystem::create_directories(filepath.parent_path());
    std::ofstream out(filepath);
    out << "seed,settle_time_s,final_error_rad,final_rate_rad_s,max_error_rad,"
           "saturated_fraction,rms_estimate_error_rad,converged\n";
    out << std::setprecision(10);
    for (std::size_t i = 0; i < seeds.size(); ++i) {
        out << seeds[i] << "," << settle_times[i] << "," << final_errors[i] << ","
            << final_rates[i] << "," << max_errors[i] << "," << sat_fractions[i] << ","
            << rms_est_errors[i] << "," << (settle_times[i] >= 0.0 ? 1 : 0) << "\n";
    }
}

}  // namespace

int main() {
    using namespace astradock;

    std::cout << "============================================================\n";
    std::cout << " AstraDock — M17 Integrated Closed-Loop GNC Demo            \n";
    std::cout << "============================================================\n";
    std::filesystem::create_directories("data");

    // M17A: attitude hold timeline (identity, zero rate).
    gnc::AttitudeReferenceTimeline hold_timeline;
    hold_timeline.segments = {{0.0, math::Quaternion::identity(), math::Vector3{}}};

    std::cout << "------------------------------------------------------------\n";
    std::cout << "Scenario M17A: closed-loop attitude hold (120 s, seed 42)\n";
    const ClosedLoopRun hold =
        run_closed_loop("data/m17_attitude_hold.csv", 42, hold_timeline, 120.0, true);
    const gnc::ClosedLoopMetrics hold_metrics = gnc::summarize_closed_loop_run(
        hold.time_s, hold.truth_error_rad, hold.desired_torque_Nm, hold.estimate_error_rad,
        hold.saturated, 0.01);
    std::cout << "  Settled: " << (hold_metrics.converged ? "yes" : "NO") << " at "
              << hold_metrics.settle_time_s << " s; final truth error "
              << hold.truth_error_rad.back() * 180.0 / constants::pi << " deg\n";

    // M17B: maneuver timeline (identity -> 60-deg yaw at t = 40 s).
    gnc::AttitudeReferenceTimeline slew_timeline;
    slew_timeline.segments = {
        {0.0, math::Quaternion::identity(), math::Vector3{}},
        {40.0,
         math::Quaternion::from_axis_angle(math::Vector3{0.0, 0.0, 1.0}, constants::pi / 3.0),
         math::Vector3{}},
    };
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Scenario M17B: 60-deg yaw maneuver at t = 40 s (seed 42)\n";
    const ClosedLoopRun slew =
        run_closed_loop("data/m17_maneuver.csv", 42, slew_timeline, 120.0, true);
    // Post-slew convergence: error relative to the NEW reference over t >= 40.
    std::vector<double> post_t, post_err, post_cmd, post_est;
    std::vector<int> post_sat;
    for (std::size_t i = 0; i < slew.time_s.size(); ++i) {
        if (slew.time_s[i] >= 40.0) {
            post_t.push_back(slew.time_s[i] - 40.0);
            post_err.push_back(slew.truth_error_rad[i]);
            post_cmd.push_back(slew.desired_torque_Nm[i]);
            post_est.push_back(slew.estimate_error_rad[i]);
            post_sat.push_back(slew.saturated[i]);
        }
    }
    const gnc::ClosedLoopMetrics slew_metrics = gnc::summarize_closed_loop_run(
        post_t, post_err, post_cmd, post_est, post_sat, 0.01);
    std::cout << "  Post-slew settled: " << (slew_metrics.converged ? "yes" : "NO") << " at +"
              << slew_metrics.settle_time_s << " s; final truth error "
              << slew.truth_error_rad.back() * 180.0 / constants::pi << " deg\n";

    // M17D: 100-run Monte Carlo on the hold scenario (seeds 1..100).
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Scenario M17D: 100-seed Monte Carlo (hold, varying noise seeds)\n";
    std::vector<std::uint64_t> seeds;
    std::vector<double> settles, finals, rates, maxes, sats, GMC_est;
    int converged = 0;
    for (std::uint64_t seed = 1; seed <= 100; ++seed) {
        const ClosedLoopRun run = run_closed_loop(
            "data/m17_mc_tmp.csv", seed * 1013 + 7, hold_timeline, 120.0, seed == 1);
        const gnc::ClosedLoopMetrics m = gnc::summarize_closed_loop_run(
            run.time_s, run.truth_error_rad, run.desired_torque_Nm, run.estimate_error_rad,
            run.saturated, 0.01);
        seeds.push_back(seed);
        settles.push_back(m.settle_time_s);
        finals.push_back(run.truth_error_rad.back());
        rates.push_back(run.rate_norm_rad_s.back());
        maxes.push_back(m.max_attitude_error_rad);
        sats.push_back(m.saturated_fraction);
        GMC_est.push_back(m.rms_estimate_error_rad);
        converged += m.converged ? 1 : 0;
    }
    // Overwrite the temp CSV with the summary table.
    write_monte_carlo_csv(
        "data/m17_monte_carlo.csv", seeds, settles, finals, rates, maxes, sats, GMC_est);
    std::filesystem::remove("data/m17_mc_tmp.csv");
    double worst_final = 0.0;
    double worst_rate = 0.0;
    for (std::size_t i = 0; i < finals.size(); ++i) {
        worst_final = std::max(worst_final, finals[i]);
        worst_rate = std::max(worst_rate, rates[i]);
    }
    std::cout << "  Converged: " << converged << " / 100; worst final error "
              << worst_final * 180.0 / constants::pi << " deg; worst final rate "
              << worst_rate * 180.0 / constants::pi << " deg/s\n";
    std::cout << "  Exported: data/m17_attitude_hold.csv, data/m17_maneuver.csv,"
                 " data/m17_monte_carlo.csv\n";
    std::cout << "============================================================\n";
    std::cout << " All M17 closed-loop demonstrations completed successfully.\n";
    std::cout << "============================================================\n";
    return 0;
}
