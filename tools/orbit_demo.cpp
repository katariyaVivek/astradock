#include "math/constants.hpp"
#include "numerics/fixed_step_propagation.hpp"
#include "orbit/cartesian_state.hpp"
#include "orbit/two_body_orbit.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using astradock::numerics::IntegrationMethod;
using astradock::numerics::StateSample;
using astradock::orbit::CartesianState;

struct PropagationMetrics {
    double final_position_closure_error_m;
    double final_velocity_closure_error_m_per_s;
    double maximum_radius_deviation_m;
    double maximum_relative_energy_drift;
    double maximum_relative_angular_momentum_drift;
};

[[nodiscard]] PropagationMetrics calculate_metrics(
    const std::vector<StateSample<CartesianState>>& samples,
    const CartesianState& initial_state,
    double gravitational_parameter_m3_per_s2) {
    const double initial_radius_m = initial_state.position.norm();
    const double initial_energy_m2_per_s2 =
        astradock::orbit::specific_orbital_energy_m2_per_s2(
            initial_state,
            gravitational_parameter_m3_per_s2);
    const double initial_angular_momentum_m2_per_s =
        astradock::orbit::specific_angular_momentum_m2_per_s(initial_state).norm();

    double maximum_radius_deviation_m = 0.0;
    double maximum_relative_energy_drift = 0.0;
    double maximum_relative_angular_momentum_drift = 0.0;
    for (const auto& sample : samples) {
        const double energy_m2_per_s2 =
            astradock::orbit::specific_orbital_energy_m2_per_s2(
                sample.state,
                gravitational_parameter_m3_per_s2);
        const double angular_momentum_m2_per_s =
            astradock::orbit::specific_angular_momentum_m2_per_s(sample.state).norm();

        maximum_radius_deviation_m = std::max(
            maximum_radius_deviation_m,
            std::abs(sample.state.position.norm() - initial_radius_m));
        maximum_relative_energy_drift = std::max(
            maximum_relative_energy_drift,
            std::abs(
                (energy_m2_per_s2 - initial_energy_m2_per_s2)
                / std::abs(initial_energy_m2_per_s2)));
        maximum_relative_angular_momentum_drift = std::max(
            maximum_relative_angular_momentum_drift,
            std::abs(
                (angular_momentum_m2_per_s - initial_angular_momentum_m2_per_s)
                / initial_angular_momentum_m2_per_s));
    }

    const CartesianState& final_state = samples.back().state;
    return {
        (final_state.position - initial_state.position).norm(),
        (final_state.velocity - initial_state.velocity).norm(),
        maximum_radius_deviation_m,
        maximum_relative_energy_drift,
        maximum_relative_angular_momentum_drift,
    };
}

void write_rows(
    std::ofstream& output,
    IntegrationMethod method,
    const std::vector<StateSample<CartesianState>>& samples,
    double gravitational_parameter_m3_per_s2,
    double earth_reference_radius_m) {
    const double initial_energy_m2_per_s2 =
        astradock::orbit::specific_orbital_energy_m2_per_s2(
            samples.front().state,
            gravitational_parameter_m3_per_s2);
    const double initial_angular_momentum_m2_per_s =
        astradock::orbit::specific_angular_momentum_m2_per_s(samples.front().state).norm();

    for (const auto& sample : samples) {
        const double radius_m = sample.state.position.norm();
        const double energy_m2_per_s2 =
            astradock::orbit::specific_orbital_energy_m2_per_s2(
                sample.state,
                gravitational_parameter_m3_per_s2);
        const double angular_momentum_m2_per_s =
            astradock::orbit::specific_angular_momentum_m2_per_s(sample.state).norm();

        output << astradock::numerics::integration_method_name(method) << ','
               << sample.time_s << ',' << sample.state.position.x() << ','
               << sample.state.position.y() << ',' << sample.state.position.z() << ','
               << sample.state.velocity.x() << ',' << sample.state.velocity.y() << ','
               << sample.state.velocity.z() << ',' << radius_m << ','
               << radius_m - earth_reference_radius_m << ','
               << sample.state.velocity.norm() << ',' << energy_m2_per_s2 << ','
               << angular_momentum_m2_per_s << ','
               << (energy_m2_per_s2 - initial_energy_m2_per_s2)
                      / std::abs(initial_energy_m2_per_s2)
               << ','
               << (angular_momentum_m2_per_s - initial_angular_momentum_m2_per_s)
                      / initial_angular_momentum_m2_per_s
               << '\n';
    }
}

void write_comparison_csv(
    const std::filesystem::path& output_path,
    const std::vector<StateSample<CartesianState>>& euler_samples,
    const std::vector<StateSample<CartesianState>>& rk4_samples,
    double gravitational_parameter_m3_per_s2,
    double earth_reference_radius_m) {
    std::ofstream output(output_path, std::ios::trunc);
    if (!output) {
        throw std::runtime_error("Could not open orbit CSV output: " + output_path.string());
    }

    output << std::setprecision(17);
    output
        << "integrator,time_s,x_m,y_m,z_m,vx_mps,vy_mps,vz_mps,radius_m,altitude_m,"
           "speed_mps,specific_energy_m2_s2,specific_angular_momentum_m2_s,"
           "relative_energy_error,relative_angular_momentum_error\n";
    write_rows(
        output,
        IntegrationMethod::forward_euler,
        euler_samples,
        gravitational_parameter_m3_per_s2,
        earth_reference_radius_m);
    write_rows(
        output,
        IntegrationMethod::classical_rk4,
        rk4_samples,
        gravitational_parameter_m3_per_s2,
        earth_reference_radius_m);

    if (!output) {
        throw std::runtime_error("Failed while writing orbit CSV output: " + output_path.string());
    }
}

void print_metrics(
    const char* experiment,
    const char* method,
    const PropagationMetrics& metrics) {
    std::cout << experiment << '.' << method
              << ".final_position_closure_error_m="
              << metrics.final_position_closure_error_m << '\n'
              << experiment << '.' << method
              << ".final_velocity_closure_error_m_per_s="
              << metrics.final_velocity_closure_error_m_per_s << '\n'
              << experiment << '.' << method
              << ".maximum_radius_deviation_m=" << metrics.maximum_radius_deviation_m << '\n'
              << experiment << '.' << method
              << ".maximum_relative_energy_drift="
              << metrics.maximum_relative_energy_drift << '\n'
              << experiment << '.' << method
              << ".maximum_relative_angular_momentum_drift="
              << metrics.maximum_relative_angular_momentum_drift << '\n';
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        if (argc > 2) {
            throw std::invalid_argument("Usage: astradock_orbit_demo [output_directory]");
        }

        const std::filesystem::path output_directory =
            argc == 2 ? std::filesystem::path{argv[1]}
                      : std::filesystem::path{"artifacts/data"};
        std::filesystem::create_directories(output_directory);

        constexpr double altitude_m = 500'000.0;
        constexpr double nominal_dt_s = 10.0;
        constexpr double gravitational_parameter_m3_per_s2 =
            astradock::constants::earth_gravitational_parameter_m3_per_s2;
        constexpr double earth_reference_radius_m =
            astradock::constants::earth_reference_radius_m;
        constexpr double orbital_radius_m = earth_reference_radius_m + altitude_m;

        const double circular_speed_m_per_s = astradock::orbit::circular_orbit_speed_m_per_s(
            gravitational_parameter_m3_per_s2,
            orbital_radius_m);
        const double orbital_period_s = astradock::orbit::circular_orbit_period_s(
            gravitational_parameter_m3_per_s2,
            orbital_radius_m);
        const CartesianState initial_state{
            {orbital_radius_m, 0.0, 0.0},
            {0.0, circular_speed_m_per_s, 0.0},
        };
        const double gravitational_acceleration_m_per_s2 =
            astradock::dynamics::two_body_acceleration(
                initial_state.position,
                gravitational_parameter_m3_per_s2)
                .norm();
        const double specific_energy_m2_per_s2 =
            astradock::orbit::specific_orbital_energy_m2_per_s2(
                initial_state,
                gravitational_parameter_m3_per_s2);
        const double specific_angular_momentum_m2_per_s =
            astradock::orbit::specific_angular_momentum_m2_per_s(initial_state).norm();
        const auto derivative = [](double time_s, const CartesianState& state) {
            return astradock::orbit::two_body_state_derivative(
                time_s,
                state,
                gravitational_parameter_m3_per_s2);
        };

        const auto one_orbit_euler = astradock::numerics::propagate_fixed_step(
            0.0,
            orbital_period_s,
            nominal_dt_s,
            initial_state,
            IntegrationMethod::forward_euler,
            derivative);
        const auto one_orbit_rk4 = astradock::numerics::propagate_fixed_step(
            0.0,
            orbital_period_s,
            nominal_dt_s,
            initial_state,
            IntegrationMethod::classical_rk4,
            derivative);
        const auto five_orbit_euler = astradock::numerics::propagate_fixed_step(
            0.0,
            5.0 * orbital_period_s,
            nominal_dt_s,
            initial_state,
            IntegrationMethod::forward_euler,
            derivative);
        const auto five_orbit_rk4 = astradock::numerics::propagate_fixed_step(
            0.0,
            5.0 * orbital_period_s,
            nominal_dt_s,
            initial_state,
            IntegrationMethod::classical_rk4,
            derivative);

        const std::filesystem::path one_orbit_csv =
            output_directory / "orbit_500km_one_orbit.csv";
        const std::filesystem::path five_orbit_csv =
            output_directory / "orbit_500km_five_orbits.csv";
        write_comparison_csv(
            one_orbit_csv,
            one_orbit_euler,
            one_orbit_rk4,
            gravitational_parameter_m3_per_s2,
            earth_reference_radius_m);
        write_comparison_csv(
            five_orbit_csv,
            five_orbit_euler,
            five_orbit_rk4,
            gravitational_parameter_m3_per_s2,
            earth_reference_radius_m);

        std::cout << std::setprecision(17)
                  << "scenario.altitude_m=" << altitude_m << '\n'
                  << "scenario.orbital_radius_m=" << orbital_radius_m << '\n'
                  << "scenario.circular_speed_m_per_s=" << circular_speed_m_per_s << '\n'
                  << "scenario.analytical_period_s=" << orbital_period_s << '\n'
                  << "scenario.gravitational_acceleration_m_per_s2="
                  << gravitational_acceleration_m_per_s2 << '\n'
                  << "scenario.specific_energy_m2_per_s2=" << specific_energy_m2_per_s2
                  << '\n'
                  << "scenario.specific_angular_momentum_m2_per_s="
                  << specific_angular_momentum_m2_per_s << '\n'
                  << "scenario.nominal_dt_s=" << nominal_dt_s << '\n';
        print_metrics(
            "one_orbit",
            "euler",
            calculate_metrics(
                one_orbit_euler,
                initial_state,
                gravitational_parameter_m3_per_s2));
        print_metrics(
            "one_orbit",
            "rk4",
            calculate_metrics(
                one_orbit_rk4,
                initial_state,
                gravitational_parameter_m3_per_s2));
        print_metrics(
            "five_orbits",
            "euler",
            calculate_metrics(
                five_orbit_euler,
                initial_state,
                gravitational_parameter_m3_per_s2));
        print_metrics(
            "five_orbits",
            "rk4",
            calculate_metrics(
                five_orbit_rk4,
                initial_state,
                gravitational_parameter_m3_per_s2));
        std::cout << "output.one_orbit_csv=" << one_orbit_csv.string() << '\n'
                  << "output.five_orbit_csv=" << five_orbit_csv.string() << '\n';
    } catch (const std::exception& error) {
        std::cerr << "astradock_orbit_demo: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
