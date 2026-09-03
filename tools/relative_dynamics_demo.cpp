#include "dynamics/two_body.hpp"
#include "frames/lvlh.hpp"
#include "frames/lvlh_rate.hpp"
#include "math/constants.hpp"
#include "math/vector3.hpp"
#include "numerics/fixed_step_propagation.hpp"
#include "orbit/cartesian_state.hpp"
#include "orbit/two_body_orbit.hpp"
#include "relative/relative_state.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

using namespace astradock;

void write_breakdown_csv(
    const std::filesystem::path& filepath,
    const std::vector<double>& sep_m,
    const std::vector<double>& err_1orbit_m,
    const std::vector<double>& err_3orbit_m) {
    std::filesystem::create_directories(filepath.parent_path());
    std::ofstream out(filepath);
    if (!out.is_open()) {
        std::cerr << "Failed to open CSV output file: " << filepath.string() << "\n";
        return;
    }
    out << "separation_m,cw_error_1orbit_m,cw_error_3orbits_m\n";
    out << std::setprecision(14);
    for (std::size_t i = 0; i < sep_m.size(); ++i) {
        out << sep_m[i] << "," << err_1orbit_m[i] << "," << err_3orbit_m[i] << "\n";
    }
}

void write_trajectory_csv(
    const std::filesystem::path& filepath,
    const std::vector<double>& time_s,
    const std::vector<relative::RelativeStateLvlh>& truth,
    const std::vector<relative::RelativeStateLvlh>& predicted) {
    std::filesystem::create_directories(filepath.parent_path());
    std::ofstream out(filepath);
    if (!out.is_open()) {
        std::cerr << "Failed to open CSV output file: " << filepath.string() << "\n";
        return;
    }
    out << "time_s,truth_x_m,truth_y_m,truth_z_m,pred_x_m,pred_y_m,pred_z_m,"
           "err_x_m,err_y_m,err_z_m,err_norm_m\n";
    out << std::setprecision(14);
    for (std::size_t i = 0; i < time_s.size(); ++i) {
        const math::Vector3 err = truth[i].relative_position_lvlh_m - predicted[i].relative_position_lvlh_m;
        out << time_s[i] << "," << truth[i].relative_position_lvlh_m.x() << ","
            << truth[i].relative_position_lvlh_m.y() << "," << truth[i].relative_position_lvlh_m.z() << ","
            << predicted[i].relative_position_lvlh_m.x() << "," << predicted[i].relative_position_lvlh_m.y()
            << "," << predicted[i].relative_position_lvlh_m.z() << "," << err.x() << "," << err.y() << ","
            << err.z() << "," << err.norm() << "\n";
    }
}

// Nonlinear two-body ECI truth for target + chaser, sampled at dt.
std::vector<orbit::CartesianState> propagate_eci(
    const orbit::CartesianState& initial, double mu, double dt, std::size_t steps) {
    auto deriv = [&](double, const orbit::CartesianState& s) {
        return orbit::two_body_state_derivative(0.0, s, mu);
    };
    auto samples = numerics::propagate_fixed_step(
        0.0, dt * static_cast<double>(steps), dt, initial, numerics::IntegrationMethod::classical_rk4, deriv);
    std::vector<orbit::CartesianState> states;
    states.reserve(samples.size());
    for (const auto& s : samples) {
        states.push_back(s.state);
    }
    return states;
}

}  // namespace

int main() {
    using namespace astradock;
    using namespace astradock::math;

    std::cout << "============================================================\n";
    std::cout << " AstraDock — M15 Relative Dynamics & CW Validation Demo     \n";
    std::cout << "============================================================\n";

    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const double R = constants::earth_reference_radius_m + 500.0e3;
    const auto ref = orbit::compute_circular_orbit_reference(mu, R);
    const double n = ref.mean_motion_rad_per_s;

    // Scenario 1: linearization breakdown sweep (M15E). Chaser starts on a
    // bounded relative ellipse (vy0 = -2n x0) at increasing along-track scales;
    // CW error vs two-body ECI truth after 1 and 3 orbits.
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Scenario 1: CW breakdown vs separation scale\n";
    const std::vector<double> separations{100.0, 500.0, 1000.0, 5000.0, 10000.0, 25000.0};
    std::vector<double> err_1, err_3;
    const double dt = 10.0;
    for (double sep : separations) {
        const orbit::CartesianState target0{{R, 0.0, 0.0}, {0.0, ref.speed_m_per_s, 0.0}};
        // Along-track offset sep with the bounded-ellipse condition in LVLH x/y.
        const relative::RelativeStateLvlh rel0{{0.0, -sep, 0.0}, {0.0, 0.0, 0.0}};
        const orbit::CartesianState chaser0 = relative::chaser_eci_from_relative(target0, rel0);
        const std::size_t steps_3 = static_cast<std::size_t>(std::round(3.0 * ref.period_s / dt));
        const auto tgt = propagate_eci(target0, mu, dt, steps_3);
        const auto cha = propagate_eci(chaser0, mu, dt, steps_3);
        const std::size_t idx_1 = static_cast<std::size_t>(std::round(ref.period_s / dt));
        double e1 = 0.0;
        double e3 = 0.0;
        for (std::size_t k : {idx_1, steps_3}) {
            const relative::RelativeStateLvlh truth = relative::relative_state_from_eci(tgt[k], cha[k]);
            const relative::RelativeStateLvlh pred =
                relative::cw_predict(rel0, n, dt * static_cast<double>(k));
            const double err =
                (truth.relative_position_lvlh_m - pred.relative_position_lvlh_m).norm();
            if (k == idx_1) {
                e1 = err;
            } else {
                e3 = err;
            }
        }
        err_1.push_back(e1);
        err_3.push_back(e3);
        std::cout << "  sep " << sep << " m: 1-orbit err " << e1 << " m, 3-orbit err " << e3 << " m\n";
    }
    write_breakdown_csv("data/m15_cw_breakdown.csv", separations, err_1, err_3);

    // Scenario 2: trajectory overlay at 5 km (truth vs CW over 1 orbit).
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Scenario 2: 5 km bounded-ellipse overlay (1 orbit)\n";
    const orbit::CartesianState target0{{R, 0.0, 0.0}, {0.0, ref.speed_m_per_s, 0.0}};
    const relative::RelativeStateLvlh rel0{{500.0, -5000.0, 200.0}, {0.0, -2.0 * n * 500.0, 0.0}};
    const orbit::CartesianState chaser0 = relative::chaser_eci_from_relative(target0, rel0);
    const std::size_t steps = static_cast<std::size_t>(std::round(ref.period_s / dt));
    const auto tgt = propagate_eci(target0, mu, dt, steps);
    const auto cha = propagate_eci(chaser0, mu, dt, steps);
    std::vector<double> time_v;
    std::vector<relative::RelativeStateLvlh> truth_v, pred_v;
    for (std::size_t k = 0; k <= steps; ++k) {
        time_v.push_back(dt * static_cast<double>(k));
        truth_v.push_back(relative::relative_state_from_eci(tgt[k], cha[k]));
        pred_v.push_back(relative::cw_predict(rel0, n, dt * static_cast<double>(k)));
    }
    write_trajectory_csv("data/m15_cw_trajectory_5km.csv", time_v, truth_v, pred_v);
    const double final_err =
        (truth_v.back().relative_position_lvlh_m - pred_v.back().relative_position_lvlh_m).norm();
    std::cout << "  Final 1-orbit CW position error at 5 km scale: " << final_err << " m\n";
    std::cout << "  Exported to: data/m15_cw_breakdown.csv, data/m15_cw_trajectory_5km.csv\n";

    std::cout << "============================================================\n";
    std::cout << " All M15 relative-dynamics demonstrations completed.\n";
    std::cout << "============================================================\n";
    return 0;
}
