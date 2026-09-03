// AstraDock M21 — Mission framework demonstration tool.
//
// Runs the canonical regression missions through the real ECI rendezvous loop
// (M18 runner, now config-driven): 20 seeded runs per mission with dispersed
// initial conditions, scored outcomes, percentile summaries, pass/fail verdicts
// with run identity + config hash stamped on every CSV.

#include "control/relative_pd.hpp"
#include "dynamics/two_body.hpp"
#include "estimation/translational_ekf.hpp"
#include "frames/lvlh.hpp"
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

// One config-driven rendezvous run: returns the scored outcome. Epoch-aligned
// + thrust-aware sequencing per M18 findings; GNSS fault optional per config.
mission::ScoredRun run_mission_case(const mission::MissionConfig& config, std::size_t run_index) {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const double R = config.orbit_radius_m;
    const auto ref = orbit::compute_circular_orbit_reference(mu, R);
    const double n = ref.mean_motion_rad_per_s;
    const double seed = static_cast<double>(mission::derive_seed(config.base_seed, run_index));
    sensors::DeterministicRng rng(static_cast<std::uint64_t>(seed));

    orbit::CartesianState target{{R, 0.0, 0.0}, {0.0, ref.speed_m_per_s, 0.0}};
    const math::Matrix3 c_init = frames::dcm_lvlh_from_eci(target.position, target.velocity);
    // Dispersed initial offset (deterministic per-run stream).
    const math::Vector3 rho_init{
        config.chaser_offset_x_m + rng.gaussian(0.0, config.mc_position_dispersion_m),
        config.chaser_offset_y_m + rng.gaussian(0.0, config.mc_position_dispersion_m),
        rng.gaussian(0.0, config.mc_position_dispersion_m)};
    orbit::CartesianState chaser{target.position + c_init.transpose() * rho_init, target.velocity};

    const std::vector<rendezvous::HoldPoint> seq{
        {math::Vector3{0.0, config.waypoint_y_m, 0.0}, config.cruise_speed_mps,
         config.capture_radius_m},
    };
    rendezvous::SafetyCorridor safety;
    const control::RelativePdGains gains{
        math::Vector3{config.kp, config.kp, config.kp},
        math::Vector3{config.kd, config.kd, config.kd}};

    estimation::TranslationalEkfConfig ekf_cfg;
    ekf_cfg.gravitational_parameter_m3_per_s2 = mu;
    ekf_cfg.acceleration_noise_std_mps2 = 1.0e-3;
    ekf_cfg.gnss_noise.position_variance_m2 =
        config.gnss_position_noise_std_m * config.gnss_position_noise_std_m;
    ekf_cfg.gnss_noise.velocity_variance_m2_per_s2 =
        config.gnss_velocity_noise_std_mps * config.gnss_velocity_noise_std_mps;
    estimation::TranslationalCovariance init_cov =
        estimation::TranslationalCovariance::zero();
    init_cov(0, 0) = 400.0;
    init_cov(1, 1) = 400.0;
    init_cov(2, 2) = 400.0;
    init_cov(3, 3) = 0.25;
    init_cov(4, 4) = 0.25;
    init_cov(5, 5) = 0.25;
    estimation::TranslationalEkf filter(
        ekf_cfg,
        orbit::CartesianState{
            chaser.position + math::Vector3{config.initial_position_error_m, 0.0, 0.0},
            chaser.velocity + math::Vector3{0.0, config.initial_velocity_error_mps, 0.0}},
        init_cov);

    auto two_body = [&](double, const orbit::CartesianState& s) {
        return orbit::two_body_state_derivative(0.0, s, mu);
    };
    std::size_t leg = 0;
    const double dt = config.dt_s;
    double t = 0.0;
    math::Vector3 a_eci_prev{};
    double delta_v = 0.0;
    double max_range = 0.0;
    double settle = -1.0;
    bool complete = false;
    bool keepout = false;
    const int steps = static_cast<int>(config.duration_s / dt);
    for (int i = 0; i < steps && !complete && !keepout; ++i) {
        target = numerics::rk4_step(t, target, dt, two_body);
        chaser = numerics::rk4_step(
            t, chaser, dt, [&](double, const orbit::CartesianState& s) {
                orbit::CartesianState d = orbit::two_body_state_derivative(0.0, s, mu);
                d.velocity = d.velocity + a_eci_prev;
                return d;
            });
        filter.predict(dt);
        {
            orbit::CartesianState corrected = filter.estimated_state();
            corrected.velocity = corrected.velocity + a_eci_prev * dt;
            corrected.position = corrected.position + a_eci_prev * (0.5 * dt * dt);
            filter.reset(corrected, filter.covariance());
        }
        const bool fault_now = config.fault_enabled && (t + dt >= config.fault_start_s)
            && (t + dt <= config.fault_end_s);
        if (!fault_now) {
            const math::Vector3 pn{
                rng.gaussian(0.0, config.gnss_position_noise_std_m),
                rng.gaussian(0.0, config.gnss_position_noise_std_m),
                rng.gaussian(0.0, config.gnss_position_noise_std_m)};
            const math::Vector3 vn{
                rng.gaussian(0.0, config.gnss_velocity_noise_std_mps),
                rng.gaussian(0.0, config.gnss_velocity_noise_std_mps),
                rng.gaussian(0.0, config.gnss_velocity_noise_std_mps)};
            filter.update_gnss(chaser.position + pn, chaser.velocity + vn);
        }
        t += dt;
        const relative::RelativeStateLvlh est_rel =
            relative::relative_state_from_eci(target, filter.estimated_state());
        auto step = rendezvous::step_rendezvous_guidance(est_rel, seq, leg, safety, n, gains);
        complete = step.mission_complete;
        keepout = step.keep_out_violation;
        if (complete && settle < 0.0) {
            settle = t;
        }
        const math::Matrix3 c_lvlh = frames::dcm_lvlh_from_eci(target.position, target.velocity);
        a_eci_prev = c_lvlh.transpose() * step.commanded_accel_lvlh_mps2;
        delta_v += step.commanded_accel_lvlh_mps2.norm() * dt;
        const relative::RelativeStateLvlh scored = relative::relative_state_from_eci(target, chaser);
        max_range = std::max(max_range, scored.relative_position_lvlh_m.norm());
    }
    const relative::RelativeStateLvlh final_rel = relative::relative_state_from_eci(target, chaser);
    mission::ScoredRun scored_run;
    scored_run.seed = mission::derive_seed(config.base_seed, run_index);
    scored_run.final_range_m = final_rel.relative_position_lvlh_m.norm();
    scored_run.max_range_m = max_range;
    scored_run.settle_time_s = settle;
    scored_run.delta_v_mps = delta_v;
    if (keepout) {
        scored_run.outcome = mission::RunOutcome::aborted;
    } else if (complete) {
        scored_run.outcome = mission::RunOutcome::converged;
    } else if (scored_run.final_range_m > 5000.0) {
        scored_run.outcome = mission::RunOutcome::diverged;
    } else {
        scored_run.outcome = mission::RunOutcome::timeout;
    }
    return scored_run;
}

void write_runs_csv(
    const std::filesystem::path& filepath,
    const mission::RunIdentity& identity,
    const std::vector<mission::ScoredRun>& runs) {
    std::filesystem::create_directories(filepath.parent_path());
    std::ofstream out(filepath);
    out << "# scenario_id," << identity.scenario_id << ",seed_base," << identity.seed
        << ",config_hash," << identity.config_hash_value << ",version," << identity.software_version
        << "\n";
    out << "seed,final_range_m,max_range_m,settle_time_s,delta_v_mps,outcome\n";
    out << std::setprecision(10);
    for (const auto& run : runs) {
        out << run.seed << "," << run.final_range_m << "," << run.max_range_m << ","
            << run.settle_time_s << "," << run.delta_v_mps << ","
            << mission::run_outcome_name(run.outcome) << "\n";
    }
}

}  // namespace

int main() {
    using namespace astradock;
    std::cout << "============================================================\n";
    std::cout << " AstraDock — M21 Mission & Monte Carlo Framework Demo       \n";
    std::cout << "============================================================\n";
    std::filesystem::create_directories("data");

    const auto missions = mission::default_regression_missions();
    std::ofstream scorecard("data/m21_scorecard.csv");
    scorecard << "mission_id,scenario_id,config_hash,runs,converged,success_rate,"
                 "p50_range_m,p95_range_m,mean_dv_mps,verdict\n";
    for (const auto& mission : missions) {
        mission.config.validate();
        const mission::RunIdentity identity{
            mission.config.scenario_id, mission.config.base_seed,
            mission::config_hash(mission.config), mission::k_astradock_version};
        std::cout << "------------------------------------------------------------\n";
        std::cout << "Mission " << mission.mission_id << " (" << identity.scenario_id
                  << ", hash " << identity.config_hash_value << ")\n";
        std::vector<mission::ScoredRun> runs;
        const std::size_t n_runs = std::min<std::size_t>(mission.config.mc_runs, 20);
        for (std::size_t i = 0; i < n_runs; ++i) {
            runs.push_back(run_mission_case(mission.config, i));
        }
        const auto summary = mission::summarize_runs(runs);
        const double metric_value = mission.metric == "success_rate" ? summary.success_rate
            : summary.p95_final_range_m;
        const bool pass = mission.passes(metric_value);
        write_runs_csv("data/m21_" + mission.config.scenario_id + ".csv", identity, runs);
        scorecard << mission.mission_id << "," << identity.scenario_id << ","
                  << identity.config_hash_value << "," << summary.runs << "," << summary.converged
                  << "," << summary.success_rate << "," << summary.p50_final_range_m << ","
                  << summary.p95_final_range_m << "," << summary.mean_delta_v_mps << ","
                  << (pass ? "PASS" : "FAIL") << "\n";
        std::cout << "  converged " << summary.converged << "/" << summary.runs << ", p50 "
                  << summary.p50_final_range_m << " m, p95 " << summary.p95_final_range_m
                  << " m, verdict " << (pass ? "PASS" : "FAIL") << "\n";
    }
    std::cout << "============================================================\n";
    std::cout << " All M21 mission demonstrations completed successfully.\n";
    std::cout << "============================================================\n";
    return 0;
}
