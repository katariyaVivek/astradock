#include "math/constants.hpp"
#include "numerics/fixed_step_propagation.hpp"
#include "orbit/cartesian_state.hpp"
#include "orbit/orbital_diagnostics.hpp"
#include "orbit/two_body_orbit.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using astradock::numerics::IntegrationMethod;
using astradock::numerics::StateSample;
using astradock::orbit::CartesianState;
using astradock::orbit::CircularOrbitReference;

// Timestep-sweep results for one integrator at one dt.
struct SweepResult {
    std::string integrator_name;
    double dt_s;
    double duration_s;
    double final_position_error_m;
    double final_velocity_error_m_per_s;
    double max_radius_error_m;
    double max_relative_energy_error;
    double max_relative_angular_momentum_error;
    double max_angular_momentum_direction_drift_rad;
    double phase_error_rad;
    double min_altitude_m;
    double max_altitude_m;
    double estimated_period_s;
    double analytical_period_s;
    double period_error_s;
};

struct CircularScenario {
    CircularOrbitReference reference;
    CartesianState initial_state;
    double earth_reference_radius_m;
};

CircularScenario make_scenario(
    double altitude_m,
    double gravitational_parameter_m3_per_s2,
    double earth_reference_radius_m) {
    const double orbital_radius_m = earth_reference_radius_m + altitude_m;
    const CircularOrbitReference ref =
        astradock::orbit::compute_circular_orbit_reference(
            gravitational_parameter_m3_per_s2, orbital_radius_m);
    const CartesianState initial_state{
        {orbital_radius_m, 0.0, 0.0},
        {0.0, ref.speed_m_per_s, 0.0},
    };

    return {
        ref,
        initial_state,
        earth_reference_radius_m,
    };
}

SweepResult run_one_orbit(
    const CircularScenario& scenario,
    double dt_s,
    IntegrationMethod method,
    const std::string& integrator_name) {
    const auto derivative = [](double t, const CartesianState& state) {
        return astradock::orbit::two_body_state_derivative(
            t,
            state,
            astradock::constants::earth_gravitational_parameter_m3_per_s2);
    };

    const auto samples = astradock::numerics::propagate_fixed_step(
        0.0,
        scenario.reference.period_s,
        dt_s,
        scenario.initial_state,
        method,
        derivative);

    if (samples.size() < 2) {
        throw std::runtime_error("Propagation returned too few samples");
    }

    const CartesianState& final_state = samples.back().state;
    const double final_position_error_m =
        astradock::orbit::diagnostics::position_closure_error_m(
            scenario.initial_state, final_state);
    const double final_velocity_error_m_per_s =
        astradock::orbit::diagnostics::velocity_closure_error_m_per_s(
            scenario.initial_state, final_state);

    // Invariant diagnostics
    const double max_radius_error_m =
        astradock::orbit::diagnostics::max_radius_deviation_m(
            samples, scenario.reference.orbital_radius_m);
    const double max_relative_energy_error =
        astradock::orbit::diagnostics::max_relative_energy_drift(
            samples, astradock::constants::earth_gravitational_parameter_m3_per_s2);
    const double max_relative_h_error =
        astradock::orbit::diagnostics::max_relative_angular_momentum_drift(samples);
    const double max_h_direction_drift_rad =
        astradock::orbit::diagnostics::max_angular_momentum_direction_drift_rad(samples);

    double min_radius = samples.front().state.position.norm();
    double max_radius = min_radius;
    for (const auto& sample : samples) {
        const double r = sample.state.position.norm();
        min_radius = std::min(min_radius, r);
        max_radius = std::max(max_radius, r);
    }

    // Phase error
    const double phase_error_rad =
        astradock::orbit::diagnostics::max_phase_error_rad(
            samples,
            astradock::constants::earth_gravitational_parameter_m3_per_s2,
            scenario.reference.orbital_radius_m);

    // Period estimation
    double estimated_period_s = 0.0;
    double period_error_s = 0.0;
    try {
        estimated_period_s =
            astradock::orbit::diagnostics::estimate_period_from_trajectory(
                samples,
                scenario.reference.period_s);
        period_error_s = estimated_period_s - scenario.reference.period_s;
    } catch (const std::domain_error&) {
        estimated_period_s = std::numeric_limits<double>::quiet_NaN();
        period_error_s = std::numeric_limits<double>::quiet_NaN();
    }

    return {
        integrator_name,
        dt_s,
        scenario.reference.period_s,
        final_position_error_m,
        final_velocity_error_m_per_s,
        max_radius_error_m,
        max_relative_energy_error,
        max_relative_h_error,
        max_h_direction_drift_rad,
        phase_error_rad,
        min_radius - scenario.earth_reference_radius_m,
        max_radius - scenario.earth_reference_radius_m,
        estimated_period_s,
        scenario.reference.period_s,
        period_error_s,
    };
}

void write_convergence_csv(
    const std::filesystem::path& path,
    const std::vector<SweepResult>& results) {
    std::ofstream out(path, std::ios::trunc);
    if (!out) {
        throw std::runtime_error("Cannot open CSV: " + path.string());
    }
    out << std::setprecision(17);
    out << "integrator,dt_s,duration_s,"
        << "final_position_error_m,final_velocity_error_m_per_s,"
        << "max_radius_error_m,"
        << "max_relative_energy_error,max_relative_angular_momentum_error,"
        << "max_angular_momentum_direction_drift_rad,"
        << "phase_error_rad,"
        << "min_altitude_m,max_altitude_m,"
        << "estimated_period_s,analytical_period_s,period_error_s\n";

    for (const auto& r : results) {
        out << r.integrator_name << ',' << r.dt_s << ',' << r.duration_s << ','
            << r.final_position_error_m << ',' << r.final_velocity_error_m_per_s << ','
            << r.max_radius_error_m << ','
            << r.max_relative_energy_error << ',' << r.max_relative_angular_momentum_error << ','
            << r.max_angular_momentum_direction_drift_rad << ','
            << r.phase_error_rad << ','
            << r.min_altitude_m << ',' << r.max_altitude_m << ','
            << r.estimated_period_s << ',' << r.analytical_period_s << ','
            << r.period_error_s << '\n';
    }
    if (!out) {
        throw std::runtime_error("CSV write failed: " + path.string());
    }
}

void print_result(const SweepResult& r) {
    std::cout << std::setprecision(10)
              << std::setw(6) << r.integrator_name << " dt=" << std::setw(4) << r.dt_s << "s"
              << "  pos_err=" << std::setw(14) << r.final_position_error_m << " m"
              << "  vel_err=" << std::setw(14) << r.final_velocity_error_m_per_s << " m/s"
              << "  rad_err=" << std::setw(14) << r.max_radius_error_m << " m"
              << "  dE/E=" << std::setw(14) << r.max_relative_energy_error
              << "  dh/h=" << std::setw(14) << r.max_relative_angular_momentum_error
              << "  phase_err=" << std::setw(14) << r.phase_error_rad << " rad"
              << "  period_err=" << std::setw(14) << r.period_error_s << " s"
              << '\n';
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        const std::filesystem::path output_directory =
            (argc >= 2) ? std::filesystem::path{argv[1]}
                        : std::filesystem::path{"artifacts/data"};
        std::filesystem::create_directories(output_directory);

        constexpr double mu =
            astradock::constants::earth_gravitational_parameter_m3_per_s2;
        constexpr double R_earth =
            astradock::constants::earth_reference_radius_m;
        constexpr double altitude_m = 500'000.0;

        const CircularScenario scenario = make_scenario(altitude_m, mu, R_earth);

        std::cout << "=== M05 Propagator Validation ===" << '\n'
                  << "altitude_m                   = " << altitude_m << " m\n"
                  << "orbital_radius_m             = " << scenario.reference.orbital_radius_m << " m\n"
                  << "circular_speed_m_per_s       = " << scenario.reference.speed_m_per_s << " m/s\n"
                  << "analytical_period_s          = " << scenario.reference.period_s << " s\n"
                  << "analytical_energy_m2_per_s2  = " << scenario.reference.specific_energy_m2_per_s2 << " m^2/s^2\n"
                  << "analytical_h_m2_per_s        = " << scenario.reference.specific_angular_momentum_m2_per_s << " m^2/s\n"
                  << "gravitational_accel_m_per_s2 = " << scenario.reference.gravitational_acceleration_m_per_s2 << " m/s^2\n"
                  << '\n';

        // Timestep sweep for RK4
        constexpr std::array rk4_dts = {40.0, 20.0, 10.0, 5.0, 2.5};
        std::cout << "--- RK4 Timestep Sweep ---" << '\n';
        std::vector<SweepResult> rk4_results;

        for (const double dt : rk4_dts) {
            try {
                const auto result = run_one_orbit(
                    scenario, dt, IntegrationMethod::classical_rk4, "rk4");
                print_result(result);
                rk4_results.push_back(result);
            } catch (const std::exception& e) {
                std::cerr << "RK4 dt=" << dt << " failed: " << e.what() << '\n';
            }
        }

        // High-accuracy numerical reference (RK4 dt = 0.25 s)
        std::cout << '\n' << "--- High-Accuracy Numerical Reference (RK4 dt=0.25s) ---" << '\n';
        const auto rk4_ref = run_one_orbit(
            scenario, 0.25, IntegrationMethod::classical_rk4, "rk4_ref");
        print_result(rk4_ref);

        // Timestep sweep for Euler
        constexpr std::array euler_dts = {40.0, 20.0, 10.0, 5.0, 2.5};
        std::cout << '\n' << "--- Euler Timestep Sweep ---" << '\n';
        std::vector<SweepResult> euler_results;

        for (const double dt : euler_dts) {
            try {
                const auto result = run_one_orbit(
                    scenario, dt, IntegrationMethod::forward_euler, "euler");
                print_result(result);
                euler_results.push_back(result);
            } catch (const std::exception& e) {
                std::cerr << "Euler dt=" << dt << " failed: " << e.what() << '\n';
            }
        }

        // Combine for CSV
        std::vector<SweepResult> all_results;
        all_results.insert(all_results.end(), euler_results.begin(), euler_results.end());
        all_results.insert(all_results.end(), rk4_results.begin(), rk4_results.end());

        const std::filesystem::path convergence_csv =
            output_directory / "m05_convergence.csv";
        write_convergence_csv(convergence_csv, all_results);
        std::cout << '\n' << "wrote " << convergence_csv.string() << '\n';

        // Write period validation CSV
        const std::filesystem::path period_csv =
            output_directory / "m05_period_validation.csv";
        {
            std::ofstream out(period_csv, std::ios::trunc);
            if (!out) {
                throw std::runtime_error("Cannot open: " + period_csv.string());
            }
            out << std::setprecision(17);
            out << "integrator,dt_s,estimated_period_s,analytical_period_s,period_error_s\n";
            for (const auto& r : all_results) {
                out << r.integrator_name << ',' << r.dt_s << ','
                    << r.estimated_period_s << ',' << r.analytical_period_s << ','
                    << r.period_error_s << '\n';
            }
            if (!out) {
                throw std::runtime_error("Write failed: " + period_csv.string());
            }
        }
        std::cout << "wrote " << period_csv.string() << '\n';

        // Write integrator comparison CSV (for Python plots)
        const std::filesystem::path comparison_csv =
            output_directory / "m05_integrator_comparison.csv";
        {
            std::ofstream out(comparison_csv, std::ios::trunc);
            if (!out) {
                throw std::runtime_error("Cannot open: " + comparison_csv.string());
            }
            out << std::setprecision(17);
            out << "integrator,dt_s,final_position_error_m,final_velocity_error_m_per_s,"
                << "max_radius_error_m,max_relative_energy_error,"
                << "max_relative_angular_momentum_error,max_angular_momentum_direction_drift_rad,phase_error_rad\n";
            for (const auto& r : all_results) {
                out << r.integrator_name << ',' << r.dt_s << ','
                    << r.final_position_error_m << ',' << r.final_velocity_error_m_per_s << ','
                    << r.max_radius_error_m << ','
                    << r.max_relative_energy_error << ','
                    << r.max_relative_angular_momentum_error << ','
                    << r.max_angular_momentum_direction_drift_rad << ','
                    << r.phase_error_rad << '\n';
            }
            if (!out) {
                throw std::runtime_error("Write failed: " + comparison_csv.string());
            }
        }
        std::cout << "wrote " << comparison_csv.string() << '\n';

        // Compute and print empirical convergence orders
        std::cout << '\n' << "--- Empirical Convergence Orders (RK4) ---" << '\n';
        std::cout << std::setprecision(6);
        for (std::size_t i = 0; i + 1 < rk4_results.size(); ++i) {
            const auto& coarser = rk4_results[i];
            const auto& finer = rk4_results[i + 1];
            double order_pos = std::numeric_limits<double>::quiet_NaN();
            if (coarser.final_position_error_m > 0 && finer.final_position_error_m > 0) {
                order_pos = astradock::orbit::diagnostics::empirical_convergence_order(
                    coarser.dt_s, coarser.final_position_error_m,
                    finer.dt_s, finer.final_position_error_m);
            }
            double order_energy = std::numeric_limits<double>::quiet_NaN();
            if (coarser.max_relative_energy_error > 0 && finer.max_relative_energy_error > 0) {
                order_energy = astradock::orbit::diagnostics::empirical_convergence_order(
                    coarser.dt_s, coarser.max_relative_energy_error,
                    finer.dt_s, finer.max_relative_energy_error);
            }
            std::cout << coarser.dt_s << "s -> " << finer.dt_s << "s: "
                      << "order_pos=" << order_pos
                      << "  order_energy=" << order_energy << '\n';
        }

        std::cout << '\n' << "--- Empirical Convergence Orders (Euler) ---" << '\n';
        for (std::size_t i = 0; i + 1 < euler_results.size(); ++i) {
            const auto& coarser = euler_results[i];
            const auto& finer = euler_results[i + 1];
            double order_pos = std::numeric_limits<double>::quiet_NaN();
            if (coarser.final_position_error_m > 0 && finer.final_position_error_m > 0) {
                order_pos = astradock::orbit::diagnostics::empirical_convergence_order(
                    coarser.dt_s, coarser.final_position_error_m,
                    finer.dt_s, finer.final_position_error_m);
            }
            double order_energy = std::numeric_limits<double>::quiet_NaN();
            if (coarser.max_relative_energy_error > 0 && finer.max_relative_energy_error > 0) {
                order_energy = astradock::orbit::diagnostics::empirical_convergence_order(
                    coarser.dt_s, coarser.max_relative_energy_error,
                    finer.dt_s, finer.max_relative_energy_error);
            }
            std::cout << coarser.dt_s << "s -> " << finer.dt_s << "s: "
                      << "order_pos=" << order_pos
                      << "  order_energy=" << order_energy << '\n';
        }

        // Repeatability check
        std::cout << '\n' << "--- Repeatability Check (RK4, dt=10s, 1 orbit) ---" << '\n';
        {
            const auto result1 = run_one_orbit(
                scenario, 10.0, IntegrationMethod::classical_rk4, "rk4");
            const auto result2 = run_one_orbit(
                scenario, 10.0, IntegrationMethod::classical_rk4, "rk4");
            const bool identical = (result1.final_position_error_m == result2.final_position_error_m)
                && (result1.final_velocity_error_m_per_s == result2.final_velocity_error_m_per_s)
                && (result1.max_radius_error_m == result2.max_radius_error_m)
                && (result1.max_relative_energy_error == result2.max_relative_energy_error)
                && (result1.phase_error_rad == result2.phase_error_rad);
            std::cout << "two runs identical: " << (identical ? "yes (deterministic)" : "NO") << '\n';
        }

    } catch (const std::exception& error) {
        std::cerr << "astradock_validation_demo: " << error.what() << '\n';
        return 1;
    }
    return 0;
}