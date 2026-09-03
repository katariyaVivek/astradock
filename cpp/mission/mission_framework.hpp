#pragma once

// AstraDock M21 — Mission scenario framework: configuration, identity,
// timeline, Monte Carlo runner, regression missions.
//
// Physical problem:
//   M13–M20 demos each hard-code their scenario: orbit, noise, faults,
//   waypoints, gains, seeds scattered across main() bodies. That is fine for
//   teaching one subsystem and fatal for a mission program: nothing is
//   comparable, repeatable, or auditable across runs. M21 centralizes the
//   scenario as DATA (MissionConfig — every knob in one struct with
//   documented defaults), stamps every run with identity (scenario ID, seed,
//   config hash, software version), drives execution through a phase timeline
//   (mode changes, fault windows, guidance legs as data, not code branches),
//   and runs seeded Monte Carlo with percentile/failure statistics.
//
// Design decisions (documented, not incidental):
//   - Plain C++ structs, no YAML/JSON parser dependency. A config FILE format
//     would add a third-party parser for zero physics benefit; the struct IS
//     the schema, printable to CSV for audit. (Revisit only if cross-language
//     scenario exchange is ever needed.)
//   - FNV-1a 64-bit config hash: deterministic across platforms for the same
//     field values (fixed field order, integer bit patterns of doubles).
//     Detects "same seed, different config" — the silent Monte Carlo killer.
//   - Phases are time windows with mode tags; the runner exposes the active
//     phase. Guidance legs / fault schedules stay in their subsystem types and
//     are REFERENCED by the config (composed, not duplicated).
//   - Statistics are online (Welford) + sorted percentiles: no storage of all
//     runs, but exact p5/median/p95. Failure classification is an enum the
//     SCORING harness sets (converged / diverged / aborted / timeout).
//   - Determinism: seed + config fully determine a run; the runner derives
//     per-run seeds as seed(frame) = base_seed * stride + run_index so runs
//     are independent yet reproducible from one base seed.
//
// Units/frames: inherit subsystems. Time in seconds mission-elapsed.

#include "math/vector3.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace astradock::mission {

namespace detail {

// FNV-1a 64-bit over raw bytes: deterministic, platform-independent for
// identical IEEE-754 double bit patterns and field order.
inline void hash_bytes(std::uint64_t& hash, const void* data, std::size_t size) noexcept {
    const auto* bytes = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < size; ++i) {
        hash ^= static_cast<std::uint64_t>(bytes[i]);
        hash *= 1099511628211ULL;
    }
}

inline void hash_double(std::uint64_t& hash, double value) noexcept {
    hash_bytes(hash, &value, sizeof(value));
}

inline void require_finite_config(double value, const char* message) {
    if (!std::isfinite(value)) {
        throw std::domain_error(message);
    }
}

}  // namespace detail

inline constexpr const char* k_astradock_version = "0.1.0";
inline constexpr std::uint64_t k_seed_stride = 1000003ULL;

// Mission phase: a named time window with a mode tag. The runner reports the
// active phase; subsystem behavior per phase stays in subsystem code.
struct MissionPhase {
    std::string name{"cruise"};
    double start_time_s{0.0};
    double end_time_s{0.0};
};

// Complete scenario configuration: every knob in one place. Defaults describe
// the M18-style final-approach rendezvous used by the regression missions.
struct MissionConfig {
    std::string scenario_id{"m21_nominal_approach"};
    std::uint64_t base_seed{42};
    // Orbit / truth.
    double orbit_radius_m{6878137.0};
    double chaser_offset_y_m{-1200.0};
    double chaser_offset_x_m{0.0};
    double duration_s{3000.0};
    double dt_s{1.0};
    // Navigation.
    double gnss_position_noise_std_m{5.0};
    double gnss_velocity_noise_std_mps{0.05};
    double initial_position_error_m{20.0};
    double initial_velocity_error_mps{0.2};
    // Guidance / control.
    double waypoint_y_m{-250.0};
    double cruise_speed_mps{0.5};
    double capture_radius_m{10.0};
    double kp{1.0e-5};
    double kd{6.0e-3};
    // Fault injection (harness-side schedule, M20 kinds as plain data).
    bool fault_enabled{false};
    double fault_start_s{200.0};
    double fault_end_s{260.0};
    double fault_bias_m{50.0};
    // Monte Carlo dispersion (harness-side sampling widths).
    double mc_position_dispersion_m{10.0};
    double mc_velocity_dispersion_mps{0.1};
    std::size_t mc_runs{100};

    void validate() const {
        if (scenario_id.empty()) {
            throw std::domain_error("Mission config scenario_id must be non-empty");
        }
        if (!(orbit_radius_m > 6378137.0) || !std::isfinite(orbit_radius_m)) {
            throw std::domain_error("Mission orbit radius must be finite and above the surface");
        }
        detail::require_finite_config(duration_s, "Mission duration must be finite");
        detail::require_finite_config(dt_s, "Mission dt must be finite");
        if (duration_s <= 0.0 || dt_s <= 0.0) {
            throw std::domain_error("Mission duration and dt must be positive");
        }
        if (!(gnss_position_noise_std_m >= 0.0) || !std::isfinite(gnss_position_noise_std_m)) {
            throw std::domain_error("GNSS noise std must be finite and non-negative");
        }
        if (mc_runs == 0) {
            throw std::domain_error("Monte Carlo run count must be positive");
        }
    }
};

// Deterministic config hash (FNV-1a over ordered fields): same config ->
// same hash on any platform; any field change flips it.
[[nodiscard]] inline std::uint64_t config_hash(const MissionConfig& config) {
    std::uint64_t hash = 1469598103934665603ULL;
    for (char c : config.scenario_id) {
        detail::hash_bytes(hash, &c, 1);
    }
    detail::hash_bytes(hash, &config.base_seed, sizeof(config.base_seed));
    detail::hash_double(hash, config.orbit_radius_m);
    detail::hash_double(hash, config.chaser_offset_y_m);
    detail::hash_double(hash, config.chaser_offset_x_m);
    detail::hash_double(hash, config.duration_s);
    detail::hash_double(hash, config.dt_s);
    detail::hash_double(hash, config.gnss_position_noise_std_m);
    detail::hash_double(hash, config.gnss_velocity_noise_std_mps);
    detail::hash_double(hash, config.initial_position_error_m);
    detail::hash_double(hash, config.initial_velocity_error_mps);
    detail::hash_double(hash, config.waypoint_y_m);
    detail::hash_double(hash, config.cruise_speed_mps);
    detail::hash_double(hash, config.capture_radius_m);
    detail::hash_double(hash, config.kp);
    detail::hash_double(hash, config.kd);
    detail::hash_bytes(hash, &config.fault_enabled, sizeof(config.fault_enabled));
    detail::hash_double(hash, config.fault_start_s);
    detail::hash_double(hash, config.fault_end_s);
    detail::hash_double(hash, config.fault_bias_m);
    detail::hash_double(hash, config.mc_position_dispersion_m);
    detail::hash_double(hash, config.mc_velocity_dispersion_mps);
    detail::hash_bytes(hash, &config.mc_runs, sizeof(config.mc_runs));
    return hash;
}

// Run identity: stamped on every run's telemetry header.
struct RunIdentity {
    std::string scenario_id{};
    std::uint64_t seed{0};
    std::uint64_t config_hash_value{0};
    std::string software_version{k_astradock_version};
};

// Per-run seed derivation: independent reproducible streams from one base.
[[nodiscard]] inline std::uint64_t derive_seed(std::uint64_t base_seed, std::size_t run_index) {
    return base_seed * k_seed_stride + static_cast<std::uint64_t>(run_index);
}

// Phase timeline query: latest phase with start <= t (mirrors M17/M18 policy).
[[nodiscard]] inline const MissionPhase& active_phase(
    const std::vector<MissionPhase>& timeline,
    double time_s) {
    if (timeline.empty()) {
        throw std::domain_error("Mission timeline must contain at least one phase");
    }
    if (!std::isfinite(time_s)) {
        throw std::domain_error("Mission time must be finite");
    }
    const MissionPhase* active = &timeline.front();
    for (const auto& phase : timeline) {
        if (phase.start_time_s <= time_s) {
            active = &phase;
        }
    }
    return *active;
}

// Default rendezvous timeline: acquire -> approach -> terminal -> capture.
[[nodiscard]] inline std::vector<MissionPhase> default_rendezvous_timeline(double duration_s) {
    return {
        {"acquire", 0.0, 300.0},
        {"approach", 300.0, duration_s - 600.0},
        {"terminal", duration_s - 600.0, duration_s - 60.0},
        {"capture", duration_s - 60.0, duration_s},
    };
}

// Outcome classification the scoring harness assigns per run.
enum class RunOutcome {
    converged,
    diverged,
    aborted,
    timeout,
};

[[nodiscard]] inline const char* run_outcome_name(RunOutcome outcome) noexcept {
    switch (outcome) {
    case RunOutcome::converged:
        return "converged";
    case RunOutcome::diverged:
        return "diverged";
    case RunOutcome::aborted:
        return "aborted";
    case RunOutcome::timeout:
        return "timeout";
    }
    return "unknown";
}

// One scored run: scalar metrics + outcome (harness fills these in).
struct ScoredRun {
    std::uint64_t seed{0};
    double final_range_m{0.0};
    double max_range_m{0.0};
    double settle_time_s{-1.0};
    double delta_v_mps{0.0};
    RunOutcome outcome{RunOutcome::timeout};
};

// Monte Carlo summary: online mean/variance (Welford) + exact percentiles
// over stored scalars (N <= thousands: storing doubles is honest, not lazy).
struct MonteCarloSummary {
    std::size_t runs{0};
    std::size_t converged{0};
    std::size_t diverged{0};
    std::size_t aborted{0};
    std::size_t timeouts{0};
    double mean_final_range_m{0.0};
    double std_final_range_m{0.0};
    double p5_final_range_m{0.0};
    double p50_final_range_m{0.0};
    double p95_final_range_m{0.0};
    double mean_delta_v_mps{0.0};
    double success_rate{0.0};
};

[[nodiscard]] inline MonteCarloSummary summarize_runs(const std::vector<ScoredRun>& runs) {
    if (runs.empty()) {
        throw std::domain_error("Monte Carlo summary needs at least one scored run");
    }
    MonteCarloSummary summary;
    summary.runs = runs.size();
    double mean = 0.0;
    double m2 = 0.0;
    double mean_dv = 0.0;
    std::vector<double> finals;
    finals.reserve(runs.size());
    for (const auto& run : runs) {
        switch (run.outcome) {
        case RunOutcome::converged:
            ++summary.converged;
            break;
        case RunOutcome::diverged:
            ++summary.diverged;
            break;
        case RunOutcome::aborted:
            ++summary.aborted;
            break;
        case RunOutcome::timeout:
            ++summary.timeouts;
            break;
        }
        finals.push_back(run.final_range_m);
        const double delta = run.final_range_m - mean;
        mean += delta / static_cast<double>(summary.converged + summary.diverged + summary.aborted + summary.timeouts);
        m2 += delta * (run.final_range_m - mean);
        mean_dv += (run.delta_v_mps - mean_dv)
            / static_cast<double>(summary.converged + summary.diverged + summary.aborted + summary.timeouts);
    }
    summary.mean_final_range_m = mean;
    summary.std_final_range_m =
        runs.size() > 1 ? std::sqrt(m2 / static_cast<double>(runs.size() - 1)) : 0.0;
    std::sort(finals.begin(), finals.end());
    const auto percentile = [&](double p) {
        const double rank = p * static_cast<double>(finals.size() - 1);
        const std::size_t lo = static_cast<std::size_t>(rank);
        const std::size_t hi = std::min(lo + 1, finals.size() - 1);
        const double frac = rank - static_cast<double>(lo);
        return finals[lo] * (1.0 - frac) + finals[hi] * frac;
    };
    summary.p5_final_range_m = percentile(0.05);
    summary.p50_final_range_m = percentile(0.50);
    summary.p95_final_range_m = percentile(0.95);
    summary.mean_delta_v_mps = mean_dv;
    summary.success_rate = static_cast<double>(summary.converged) / static_cast<double>(runs.size());
    return summary;
}

// Canonical regression missions: fixed configs with acceptance criteria.
// Each mission pins scenario_id + key parameters; the harness checks the
// criterion on the summary (kept as data: metric name + bound + direction).
struct RegressionMission {
    std::string mission_id{};
    MissionConfig config{};
    std::string metric{"success_rate"};
    double bound{0.95};
    bool higher_is_better{true};

    [[nodiscard]] bool passes(double metric_value) const {
        if (!std::isfinite(metric_value)) {
            throw std::domain_error("Regression metric value must be finite");
        }
        return higher_is_better ? metric_value >= bound : metric_value <= bound;
    }
};

[[nodiscard]] inline std::vector<RegressionMission> default_regression_missions() {
    MissionConfig nominal;
    nominal.scenario_id = "rendezvous_nominal";
    MissionConfig dropout = nominal;
    dropout.scenario_id = "rendezvous_gnss_dropout";
    dropout.fault_enabled = true;
    MissionConfig hot = nominal;
    hot.scenario_id = "rendezvous_hot_close";
    hot.cruise_speed_mps = 1.0;
    MissionConfig docking_grade = nominal;
    docking_grade.scenario_id = "docking_grade_nav";
    docking_grade.gnss_position_noise_std_m = 0.05;
    docking_grade.waypoint_y_m = -50.0;
    return {
        {"m21_nominal_converges", nominal, "success_rate", 0.95, true},
        {"m21_dropout_rides_through", dropout, "success_rate", 0.90, true},
        {"m21_hot_close_bounded", hot, "p95_final_range_m", 500.0, false},
        {"m21_docking_grade_capture", docking_grade, "success_rate", 0.90, true},
    };
}

}  // namespace astradock::mission
