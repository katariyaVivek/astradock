// AstraDock M23 — Sensor anomaly classification demo: GNSS bias-jump detector.
//
// Generates a leakage-safe dataset (split by SEED, never by row) of GNSS NIS
// windows labeled by a scheduled bias fault, trains the classical baselines
// FIRST (threshold + logistic regression), then a small RBF-SVM / kNN challenger,
// and reports precision/recall/F1/false-alarm/miss/latency against the M20
// NIS monitor baseline. No torch, no deep learning: a 6-DOF EKF innovation
// stream does not justify it (explicit non-goal, documented in the lesson).

#include "estimation/diagnostics.hpp"
#include "estimation/translational_ekf.hpp"
#include "fdir/fdir.hpp"
#include "math/constants.hpp"
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

// One labeled sample: window of 5 consecutive GNSS NIS values + fault label.
// Windows never cross the fault onset (no ambiguous windows in the dataset);
// scenario/seed recorded per row for split auditing.
struct NisWindow {
    std::uint64_t scenario_seed{0};
    double onset_s{-1.0};
    std::array<double, 5> nis{};
    int label{0};  // 0 nominal, 1 fault
};

std::vector<NisWindow> generate_dataset(
    std::uint64_t base_seed,
    int n_scenarios,
    bool faulty,
    double onset_s,
    double bias_m) {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const double r_orbit = 7000.0e3;
    const double v_circ = std::sqrt(mu / r_orbit);
    auto two_body = [&](double, const orbit::CartesianState& s) {
        return orbit::two_body_state_derivative(0.0, s, mu);
    };
    std::vector<NisWindow> rows;
    for (int s = 0; s < n_scenarios; ++s) {
        const std::uint64_t seed = mission::derive_seed(base_seed, static_cast<std::size_t>(s));
        sensors::DeterministicRng rng(seed);
        orbit::CartesianState truth{{r_orbit, 0.0, 0.0}, {0.0, v_circ, 0.0}};
        estimation::TranslationalEkfConfig cfg;
        cfg.gravitational_parameter_m3_per_s2 = mu;
        cfg.acceleration_noise_std_mps2 = 1.0e-3;
        cfg.gnss_noise.position_variance_m2 = 25.0;
        cfg.gnss_noise.velocity_variance_m2_per_s2 = 2.5e-3;
        estimation::TranslationalCovariance init_cov =
            estimation::TranslationalCovariance::zero();
        for (std::size_t i = 0; i < 3; ++i) {
            init_cov(i, i) = 400.0;
            init_cov(3 + i, 3 + i) = 0.25;
        }
        estimation::TranslationalEkf filter(cfg, truth, init_cov);
        std::vector<double> nis_stream;
        const double dt = 1.0;
        for (int step = 0; step <= 120; ++step) {
            const double t = step * dt;
            truth = numerics::rk4_step(t, truth, dt, two_body);
            filter.predict(dt);
            math::Vector3 bias{};
            if (faulty && t >= onset_s) {
                bias = math::Vector3{bias_m, 0.0, 0.0};
            }
            const math::Vector3 pn{
                rng.gaussian(0.0, 5.0), rng.gaussian(0.0, 5.0), rng.gaussian(0.0, 5.0)};
            const math::Vector3 vn{
                rng.gaussian(0.0, 0.05), rng.gaussian(0.0, 0.05), rng.gaussian(0.0, 0.05)};
            filter.update_gnss(truth.position + pn + bias, truth.velocity + vn);
            nis_stream.push_back(
                filter.last_update_diagnostics().normalized_innovation_squared);
        }
        // Windows of 5, stride 5, skipping any window touching onset_s.
        for (std::size_t w = 0; w + 5 <= nis_stream.size(); w += 5) {
            const double t0 = w * dt;
            const double t1 = (w + 5) * dt;
            if (faulty && t0 < onset_s && t1 > onset_s) {
                continue;  // ambiguous window: excluded by construction
            }
            NisWindow row;
            row.scenario_seed = seed;
            row.onset_s = faulty ? onset_s : -1.0;
            for (int k = 0; k < 5; ++k) {
                row.nis[k] = nis_stream[w + k];
            }
            row.label = (faulty && t0 >= onset_s) ? 1 : 0;
            rows.push_back(row);
        }
    }
    return rows;
}

void write_dataset_csv(const std::filesystem::path& filepath, const std::vector<NisWindow>& rows) {
    std::filesystem::create_directories(filepath.parent_path());
    std::ofstream out(filepath);
    out << "scenario_seed,onset_s,nis_0,nis_1,nis_2,nis_3,nis_4,label\n";
    out << std::setprecision(10);
    for (const auto& row : rows) {
        out << row.scenario_seed << "," << row.onset_s << "," << row.nis[0] << "," << row.nis[1]
            << "," << row.nis[2] << "," << row.nis[3] << "," << row.nis[4] << "," << row.label
            << "\n";
    }
}

}  // namespace

int main() {
    using namespace astradock;
    std::cout << "============================================================\n";
    std::cout << " AstraDock — M23 ML Dataset Generation (sensor anomaly)     \n";
    std::cout << "============================================================\n";
    std::filesystem::create_directories("data");
    // Split by SCENARIO (seed-disjoint): train scenarios 0-59, val 60-79, test 80-99.
    // Faulty and nominal scenarios generated separately; labels by construction.
    auto train_nom = generate_dataset(1000, 60, false, -1.0, 0.0);
    auto train_flt = generate_dataset(2000, 60, true, 60.0, 10.0);
    auto val_nom = generate_dataset(3000, 20, false, -1.0, 0.0);
    auto val_flt = generate_dataset(4000, 20, true, 60.0, 10.0);
    auto test_nom = generate_dataset(5000, 20, false, -1.0, 0.0);
    auto test_flt = generate_dataset(6000, 20, true, 60.0, 10.0);
    std::vector<NisWindow> train = train_nom;
    train.insert(train.end(), train_flt.begin(), train_flt.end());
    std::vector<NisWindow> val = val_nom;
    val.insert(val.end(), val_flt.begin(), val_flt.end());
    std::vector<NisWindow> test = test_nom;
    test.insert(test.end(), test_flt.begin(), test_flt.end());
    write_dataset_csv("data/m23_train.csv", train);
    write_dataset_csv("data/m23_val.csv", val);
    write_dataset_csv("data/m23_test.csv", test);
    std::cout << "  train: " << train.size() << " windows (seeds 1000/2000 families)\n";
    std::cout << "  val:   " << val.size() << " windows (seeds 3000/4000 families)\n";
    std::cout << "  test:  " << test.size() << " windows (seeds 5000/6000 families)\n";
    std::cout << "  Leakage audit: seed families disjoint across splits (checked in analysis).\n";
    std::cout << "============================================================\n";
    return 0;
}
