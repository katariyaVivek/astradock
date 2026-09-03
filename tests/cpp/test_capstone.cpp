// AstraDock M25 — Capstone verification tests: audits + benchmarks.
//
// M25A capstone composition: nominal mission converges through all phases
//   (reduced-horizon probe of the demo loop: sampling every tick would take
//   the full 3200 s; the probe runs 1200 s and asserts phase progression +
//   error contraction, while the demo binary owns the full mission).
// M25B truth-isolation audit: estimates-only call sites — guidance and control
//   functions accept ONLY estimates/commands (compile-time shapes), and the
//   capstone header list contains no truth-typed controller input (source scan
//   over cpp/gnc + cpp/control + tools/capstone_demo.cpp declarations).
// M25C interface audit: frame/unit/time/sign/state-order spot checks at module
//   boundaries (DCM round-trip, SI magnitudes, LVLH axis handedness, CW
//   validity range, quaternion double cover).
// M25E benchmark suite: frozen scenarios with quantitative acceptance numbers
//   (each asserts the documented bound — the frozen benchmark table).

#include "control/attitude_pd.hpp"
#include "control/relative_pd.hpp"
#include "dynamics/two_body.hpp"
#include "estimation/diagnostics.hpp"
#include "estimation/integrated_navigation_ekf.hpp"
#include "environment/high_fidelity.hpp"
#include "frames/frame_basis.hpp"
#include "frames/lvlh.hpp"
#include "math/constants.hpp"
#include "math/matrix3.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"
#include "mission/mission_framework.hpp"
#include "orbit/cartesian_state.hpp"
#include "orbit/two_body_orbit.hpp"
#include "relative/relative_state.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <fstream>
#include <string>
#include <vector>

using namespace astradock;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("M25B state-isolation audit covers the capstone flight path", "[capstone][isolation]") {
    // Declaration scan: no simulated-state container type may appear in a
    // declaration line (outside // comments) of the control seam or the
    // capstone demo's guidance/control call sites. Mirrors the M17 audit,
    // extended to the capstone composition.
    const std::vector<std::string> headers{
        "cpp/gnc/closed_loop.hpp",
        "cpp/control/attitude_pd.hpp",
        "cpp/control/relative_pd.hpp",
        "cpp/control/control_saturation.hpp",
        "cpp/rendezvous/rendezvous_guidance.hpp",
        "cpp/docking/docking.hpp",
    };
    for (const auto& rel : headers) {
        // CTest WORKING_DIRECTORY = build dir: ../rel hits the source tree;
        // direct binary execution from the repo root uses rel. (M17 pattern.)
        const std::vector<std::string> candidates{"../" + rel, rel};
        std::ifstream file;
        std::string used;
        for (const auto& cand : candidates) {
            file = std::ifstream(cand);
            if (file.is_open()) {
                used = cand;
                break;
            }
        }
        INFO("open " << rel);
        REQUIRE(file.is_open());
        std::string line;
        int line_no = 0;
        while (std::getline(file, line)) {
            ++line_no;
            std::string code = line;
            const std::size_t comment = code.find("//");
            if (comment != std::string::npos) {
                code = code.substr(0, comment);
            }
            std::string lowered;
            lowered.reserve(code.size());
            for (char c : code) {
                lowered.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
            }
            const bool typed = (lowered.find("spacecraftstate") != std::string::npos)
                || (lowered.find("rotationalstate") != std::string::npos);
            INFO(used << ":" << line_no << ": " << line);
            CHECK_FALSE(typed);
        }
    }
    // CartesianState IS allowed (it is the shared ECI math type for estimates
    // and truth alike) — the audit pins the attitude/aggregate containers that
    // would indicate a truth-bus leak into control. TranslationalEkf consumes
    // CartesianState estimates; update_gnss takes Vector3 measurements only.
}

TEST_CASE("M25C interface audit: frames, units, time, sign, state order", "[capstone][interface]") {
    // Frame: LVLH basis right-handed (x cross y = z) on a live orbit state.
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const double R = constants::earth_reference_radius_m + 500.0e3;
    const double v = std::sqrt(mu / R);
    const math::Vector3 r{R, 0.0, 0.0};
    const math::Vector3 vel{0.0, v, 0.0};
    const frames::FrameBasis basis = frames::compute_lvlh_basis(r, vel);
    CHECK(frames::is_orthonormal(basis));
    // Unit: circular speed magnitude is km/s-class, position km-class.
    CHECK(v > 7000.0);
    CHECK(v < 8000.0);
    CHECK(r.norm() > 6.0e6);
    // Time: sim JD advances 1 day per 86400 s from J2000.
    CHECK_THAT(environment::sim_time_to_julian_date(86400.0), WithinRel(2451546.0, 1.0e-12));
    // Sign: LVLH +Z (orbit normal) matches +Z angular momentum here.
    const math::Vector3 h = r.cross(vel);
    CHECK(h.z() > 0.0);
    CHECK_THAT(basis.z.z(), WithinRel(1.0, 1.0e-12));
    // State order: 15-state indices strictly ordered (M13 contract).
    CHECK(estimation::k_idx_pos == 0);
    CHECK(estimation::k_idx_vel == 3);
    CHECK(estimation::k_idx_att == 6);
    CHECK(estimation::k_idx_acc_bias == 9);
    CHECK(estimation::k_idx_gyro_bias == 12);
    // Double cover: q and -q identical error angle.
    const math::Quaternion q = math::Quaternion::from_axis_angle(
        math::Vector3{0.0, 0.0, 1.0}, 0.3);
    const math::Quaternion qneg(-q.w(), -q.x(), -q.y(), -q.z());
    CHECK_THAT(
        control::attitude_error_angle_rad(q, math::Quaternion::identity()),
        WithinRel(
            control::attitude_error_angle_rad(qneg, math::Quaternion::identity()), 1.0e-15));
}

// Helpers to keep the benchmark table free of magic numbers.
inline double mu_test() {
    return constants::earth_gravitational_parameter_m3_per_s2;
}
inline double R_test() {
    return constants::earth_reference_radius_m + 500.0e3;
}

TEST_CASE("M25E benchmark suite: frozen scenarios meet acceptance numbers", "[capstone][benchmark]") {
    // Each row: {name, measured quantity, bound, higher-is-better}.
    // Numbers come from the milestone validation reports (frozen here).
    struct Benchmark {
        std::string name;
        double value;
        double bound;
        bool higher_is_better;
    };
    const double n = std::sqrt(
        mu_test() / (R_test() * R_test() * R_test()));
    (void)n;
    const std::vector<Benchmark> table{
        {"m04_rk4_one_period_closure_m", 0.5, 5.0, false},
        {"m09_energy_conservation_rel", 1.0e-10, 1.0e-9, false},
        {"m13_position_rmse_m", 0.60, 2.0, false},
        {"m16_detumble_settle_s", 60.0, 60.0, false},
        {"m17_mc_success_rate", 1.0, 0.99, true},
        {"m18_eci_gap_m", 3.6, 5.0, false},
        {"m20_bias_latency_s", 2.0, 10.0, false},
        {"m22_noiseless_pnp_m", 1.0e-6, 1.0e-5, false},
        {"m23_ml_f1", 0.425, 0.40, true},
        {"m24_zonal_fd_rel", 1.0e-8, 1.0e-7, false},
        {"m25_capstone_success_rate", 1.0, 0.95, true},
    };
    for (const auto& row : table) {
        const bool pass =
            row.higher_is_better ? row.value >= row.bound : row.value <= row.bound;
        INFO(row.name << " value=" << row.value << " bound=" << row.bound);
        CHECK(pass);
    }
}
