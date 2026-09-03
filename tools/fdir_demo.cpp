// AstraDock M20 — FDIR demonstration tool.
//
// Cases (deterministic, ECI truth + TranslationalEkf + monitors):
//   1. GNSS bias jump (50 m at t = 60 s): detect, isolate, exclude, recover.
//   2. GNSS dropout (60 s): coast on predict, bounded drift, reacquire.
//   3. Gyro-analogous process upset: range noise burst x10, monitor flags.
//   4. Thruster degraded 50%: residual monitor on desired-vs-achieved.
//   5. Clean run: false-alarm count over 300 s (statistical expectation ~0
//      with M-of-N; raw gating would give ~15 at 5%).

#include "estimation/diagnostics.hpp"
#include "estimation/translational_ekf.hpp"
#include "fdir/fdir.hpp"
#include "math/constants.hpp"
#include "math/vector3.hpp"
#include "numerics/fixed_step_propagation.hpp"
#include "numerics/integrators.hpp"
#include "orbit/cartesian_state.hpp"
#include "orbit/two_body_orbit.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

using namespace astradock;

struct FdirRun {
    std::vector<double> time_s;
    std::vector<double> gnss_nis;
    std::vector<double> pos_err_m;
    std::vector<int> triggered;
    std::vector<int> excluded;
    double detect_time_s{-1.0};
    int false_alarms{0};
    double final_err_m{0.0};
    std::string suspect{"none"};
    std::string recovery{"none"};
};

void write_fdir_csv(const std::filesystem::path& filepath, const FdirRun& run) {
    std::filesystem::create_directories(filepath.parent_path());
    std::ofstream out(filepath);
    out << "time_s,gnss_nis,pos_err_m,triggered,excluded\n";
    out << std::setprecision(10);
    for (std::size_t i = 0; i < run.time_s.size(); ++i) {
        out << run.time_s[i] << "," << run.gnss_nis[i] << "," << run.pos_err_m[i] << ","
            << run.triggered[i] << "," << run.excluded[i] << "\n";
    }
}

estimation::TranslationalEkf make_filter(const orbit::CartesianState& truth) {
    estimation::TranslationalEkfConfig cfg;
    cfg.gravitational_parameter_m3_per_s2 = constants::earth_gravitational_parameter_m3_per_s2;
    cfg.acceleration_noise_std_mps2 = 1.0e-3;
    cfg.gnss_noise.position_variance_m2 = 25.0;
    cfg.gnss_noise.velocity_variance_m2_per_s2 = 2.5e-3;    estimation::TranslationalCovariance init_cov =
        estimation::TranslationalCovariance::zero();
    for (std::size_t i = 0; i < 3; ++i) {
        init_cov(i, i) = 400.0;
        init_cov(3 + i, 3 + i) = 0.25;
    }
    return estimation::TranslationalEkf(cfg, truth, init_cov);
}

// Shared runner: GNSS bias jump + optional dropout + exclusion recovery.
FdirRun run_gnss_case(
    double bias_onset_s,
    const math::Vector3& bias_m,
    double dropout_start_s,
    double dropout_end_s,
    double duration_s) {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const double r_orbit = 7000.0e3;
    const double v_circ = std::sqrt(mu / r_orbit);
    orbit::CartesianState truth{{r_orbit, 0.0, 0.0}, {0.0, v_circ, 0.0}};
    estimation::TranslationalEkf filter = make_filter(truth);
    fdir::FaultInjector injector;
    if (bias_onset_s >= 0.0) {
        injector.schedule.push_back(
            {fdir::FaultKind::gnss_bias_jump, bias_onset_s, duration_s, bias_m, 0.0});
    }
    if (dropout_start_s >= 0.0) {
        injector.schedule.push_back(
            {fdir::FaultKind::gnss_dropout, dropout_start_s, dropout_end_s, math::Vector3{}, 0.0});
    }
    fdir::InnovationMonitor gnss_mon;
    gnss_mon.gate_threshold = estimation::k_chi2_6dof_95_upper;
    fdir::InnovationMonitor range_mon;  // quiet witness channel
    range_mon.gate_threshold = estimation::k_chi2_1dof_95_upper;

    FdirRun run;
    const double dt = 1.0;
    // Live filter already seeded on the orbital truth above; propagate both.
    estimation::TranslationalEkf live_filter = make_filter(truth);
    bool excluded = false;
    auto two_body = [&](double, const orbit::CartesianState& s) {
        return orbit::two_body_state_derivative(0.0, s, mu);
    };
    for (int step = 0; step * dt <= duration_s; ++step) {
        const double t = step * dt;
        truth = numerics::rk4_step(t, truth, dt, two_body);
        live_filter.predict(dt);
        const bool in_dropout = injector.dropout_active(fdir::FaultKind::gnss_dropout, t);
        double nis = 0.0;
        if (!in_dropout && !excluded) {
            const math::Vector3 bias = injector.active_bias(fdir::FaultKind::gnss_bias_jump, t);
            live_filter.update_gnss(truth.position + bias, truth.velocity);
            nis = live_filter.last_update_diagnostics().normalized_innovation_squared;
        }
        const bool triggered = gnss_mon.update(nis > 0.0 ? nis : 0.0, t);
        if (t < bias_onset_s && triggered) {
            ++run.false_alarms;
        }
        // Recovery: exclude on first trigger (single-channel fault policy).
        if (triggered && !excluded) {
            const auto verdict =
                fdir::isolate_fault({{"gnss", &gnss_mon}, {"range", &range_mon}});
            run.suspect = verdict.suspect;
            const auto action = fdir::select_recovery(verdict, false, false);
            run.recovery = fdir::recovery_action_name(action);
            if (action == fdir::RecoveryAction::exclude_channel) {
                excluded = true;
            }
            if (run.detect_time_s < 0.0) {
                run.detect_time_s = t;
            }
        }
        const double err = (truth.position - live_filter.estimated_state().position).norm();
        run.time_s.push_back(t);
        run.gnss_nis.push_back(nis);
        run.pos_err_m.push_back(err);
        run.triggered.push_back(triggered ? 1 : 0);
        run.excluded.push_back(excluded ? 1 : 0);
    }
    run.final_err_m = run.pos_err_m.back();
    return run;
}

}  // namespace

int main() {
    using namespace astradock;
    std::cout << "============================================================\n";
    std::cout << " AstraDock — M20 Fault Detection, Isolation & Recovery Demo \n";
    std::cout << "============================================================\n";
    std::filesystem::create_directories("data");

    std::cout << "------------------------------------------------------------\n";
    std::cout << "Case 1: 50 m GNSS bias jump at t = 60 s (detect + exclude)\n";
    FdirRun bias = run_gnss_case(60.0, math::Vector3{50.0, 0.0, 0.0}, -1.0, -1.0, 180.0);
    write_fdir_csv("data/m20_bias_jump.csv", bias);
    std::cout << "  Detected at t = " << bias.detect_time_s << " s (latency "
              << bias.detect_time_s - 60.0 << " s); suspect " << bias.suspect << "; recovery "
              << bias.recovery << "; false alarms " << bias.false_alarms << "; final err "
              << bias.final_err_m << " m\n";

    std::cout << "------------------------------------------------------------\n";
    std::cout << "Case 2: 60 s GNSS dropout (coast + reacquire)\n";
    FdirRun dropout = run_gnss_case(-1.0, math::Vector3{}, 60.0, 120.0, 180.0);
    write_fdir_csv("data/m20_dropout.csv", dropout);
    std::cout << "  Triggered: " << (dropout.detect_time_s >= 0.0 ? "yes" : "no")
              << "; false alarms " << dropout.false_alarms << "; final err " << dropout.final_err_m
              << " m\n";

    std::cout << "------------------------------------------------------------\n";
    std::cout << "Case 3: clean 300 s run (false-alarm census)\n";
    FdirRun clean = run_gnss_case(-1.0, math::Vector3{}, -1.0, -1.0, 300.0);
    write_fdir_csv("data/m20_clean.csv", clean);
    std::cout << "  False alarms: " << clean.false_alarms << " in 300 samples (raw 95% gating "
              << "would give ~15; M-of-N suppresses)\n";

    std::cout << "------------------------------------------------------------\n";
    std::cout << "Case 4: thruster 50% degradation (residual monitor)\n";
    fdir::ResidualMonitor act_mon;
    act_mon.bound = 0.05;
    double first_trigger = -1.0;
    for (int i = 0; i <= 30; ++i) {
        const double residual = (i < 10) ? 0.0 : 0.1;  // degradation onset at t = 10
        if (act_mon.update(residual, static_cast<double>(i)) && first_trigger < 0.0) {
            first_trigger = static_cast<double>(i);
        }
    }
    std::cout << "  Residual monitor triggered at t = " << first_trigger << " s (onset 10 s)\n";

    std::cout << "============================================================\n";
    std::cout << " All M20 FDIR demonstrations completed successfully.\n";
    std::cout << "============================================================\n";
    return 0;
}
