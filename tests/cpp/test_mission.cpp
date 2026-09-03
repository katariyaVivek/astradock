// AstraDock M21 — Mission framework tests.
//
// Verification strategy:
//   Config: defaults validate; each invalid field throws; hash deterministic +
//     sensitive to every field group (flip one field -> hash changes).
//   Identity/seeds: derive_seed streams distinct + reproducible.
//   Timeline: phase selection incl. boundaries; empty/singleton rejection.
//   Statistics: Welford mean/variance vs closed form; percentiles on a known
//     uniform grid; outcome counts; empty rejection.
//   Regression: passes() directions; mission table has 4 entries with sane
//     bounds; non-finite metric rejected.

#include "mission/mission_framework.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

using namespace astradock;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("M21 mission config validates and hashes deterministically", "[mission][config]") {
    mission::MissionConfig config;
    CHECK_NOTHROW(config.validate());
    CHECK(mission::config_hash(config) == mission::config_hash(mission::MissionConfig{}));
    // Sensitivity: flipping any field group changes the hash.
    mission::MissionConfig altered = config;
    altered.scenario_id = "other";
    CHECK(mission::config_hash(altered) != mission::config_hash(config));
    altered = config;
    altered.base_seed = 999;
    CHECK(mission::config_hash(altered) != mission::config_hash(config));
    altered = config;
    altered.orbit_radius_m += 1.0;
    CHECK(mission::config_hash(altered) != mission::config_hash(config));
    altered = config;
    altered.kp *= 2.0;
    CHECK(mission::config_hash(altered) != mission::config_hash(config));
    altered = config;
    altered.fault_enabled = true;
    CHECK(mission::config_hash(altered) != mission::config_hash(config));
    altered = config;
    altered.mc_runs = 50;
    CHECK(mission::config_hash(altered) != mission::config_hash(config));
    // Invalid fields rejected.
    mission::MissionConfig bad = config;
    bad.scenario_id = "";
    CHECK_THROWS_AS(bad.validate(), std::domain_error);
    bad = config;
    bad.orbit_radius_m = 6378137.0;
    CHECK_THROWS_AS(bad.validate(), std::domain_error);
    bad = config;
    bad.dt_s = 0.0;
    CHECK_THROWS_AS(bad.validate(), std::domain_error);
    bad = config;
    bad.mc_runs = 0;
    CHECK_THROWS_AS(bad.validate(), std::domain_error);
    bad = config;
    bad.gnss_position_noise_std_m = -1.0;
    CHECK_THROWS_AS(bad.validate(), std::domain_error);
}

TEST_CASE("M21 run identity and seed derivation are reproducible", "[mission][identity]") {
    mission::MissionConfig config;
    config.scenario_id = "m21_probe";
    config.base_seed = 7;
    const mission::RunIdentity identity{
        config.scenario_id, mission::derive_seed(config.base_seed, 0),
        mission::config_hash(config), mission::k_astradock_version};
    CHECK(identity.scenario_id == "m21_probe");
    CHECK(identity.seed == 7 * mission::k_seed_stride);
    CHECK(std::string(identity.software_version) == "0.1.0");
    // Streams distinct across runs, stable across calls.
    CHECK(mission::derive_seed(7, 0) != mission::derive_seed(7, 1));
    CHECK(mission::derive_seed(7, 1) == mission::derive_seed(7, 1));
    CHECK(mission::derive_seed(7, 5) != mission::derive_seed(8, 5));
}

TEST_CASE("M21 phase timeline selects deterministically", "[mission][timeline]") {
    const auto timeline = mission::default_rendezvous_timeline(3000.0);
    REQUIRE(timeline.size() == 4);
    CHECK(mission::active_phase(timeline, 0.0).name == "acquire");
    CHECK(mission::active_phase(timeline, 300.0).name == "approach");
    CHECK(mission::active_phase(timeline, 2400.0).name == "terminal");
    CHECK(mission::active_phase(timeline, 2940.0).name == "capture");
    CHECK(mission::active_phase(timeline, 2999.0).name == "capture");
    CHECK_THROWS_AS(mission::active_phase({}, 0.0), std::domain_error);
    const std::vector<mission::MissionPhase> singleton{{"only", 0.0, 10.0}};
    CHECK(mission::active_phase(singleton, 99.0).name == "only");
}

TEST_CASE("M21 Monte Carlo summary matches closed-form statistics", "[mission][statistics]") {
    // Uniform grid 1..100: mean 50.5, sample std sqrt(841.666..), p5/p50/p95 known.
    std::vector<mission::ScoredRun> runs;
    for (std::uint64_t i = 1; i <= 100; ++i) {
        const double v = static_cast<double>(i);
        runs.push_back({i, v, v, v, 0.1 * v,
                        i <= 95 ? mission::RunOutcome::converged : mission::RunOutcome::diverged});
    }
    const auto summary = mission::summarize_runs(runs);
    CHECK(summary.runs == 100);
    CHECK(summary.converged == 95);
    CHECK(summary.diverged == 5);
    CHECK_THAT(summary.mean_final_range_m, WithinRel(50.5, 1.0e-12));
    CHECK_THAT(summary.std_final_range_m, WithinRel(std::sqrt(100.0 * 101.0 / 12.0), 1.0e-9));
    CHECK_THAT(summary.p5_final_range_m, WithinRel(5.95, 1.0e-12));
    CHECK_THAT(summary.p50_final_range_m, WithinRel(50.5, 1.0e-12));
    CHECK_THAT(summary.p95_final_range_m, WithinRel(95.05, 1.0e-12));
    CHECK_THAT(summary.success_rate, WithinRel(0.95, 1.0e-15));
    CHECK_THAT(summary.mean_delta_v_mps, WithinRel(5.05, 1.0e-12));
    CHECK_THROWS_AS(mission::summarize_runs({}), std::domain_error);
}

TEST_CASE("M21 outcome counts classify every run exactly once", "[mission][statistics]") {
    const std::vector<mission::ScoredRun> runs{
        {0, 1.0, 1.0, 1.0, 0.0, mission::RunOutcome::converged},
        {1, 2.0, 2.0, -1.0, 0.0, mission::RunOutcome::diverged},
        {2, 3.0, 3.0, -1.0, 0.0, mission::RunOutcome::aborted},
        {3, 4.0, 4.0, -1.0, 0.0, mission::RunOutcome::timeout},
    };
    const auto summary = mission::summarize_runs(runs);
    CHECK(summary.converged == 1);
    CHECK(summary.diverged == 1);
    CHECK(summary.aborted == 1);
    CHECK(summary.timeouts == 1);
    CHECK_THAT(summary.success_rate, WithinRel(0.25, 1.0e-15));
    CHECK(std::string(mission::run_outcome_name(mission::RunOutcome::aborted)) == "aborted");
}

TEST_CASE("M21 regression missions carry pass/fail criteria", "[mission][regression]") {
    const auto missions = mission::default_regression_missions();
    REQUIRE(missions.size() == 4);
    CHECK(missions[0].mission_id == "m21_nominal_converges");
    CHECK(missions[0].passes(0.99));
    CHECK_FALSE(missions[0].passes(0.90));
    // Lower-is-better direction (p95 bound).
    CHECK(missions[2].passes(100.0));
    CHECK_FALSE(missions[2].passes(600.0));
    CHECK_THROWS_AS(missions[0].passes(std::numeric_limits<double>::quiet_NaN()), std::domain_error);
    for (const auto& mission : missions) {
        CHECK_NOTHROW(mission.config.validate());
    }
}
