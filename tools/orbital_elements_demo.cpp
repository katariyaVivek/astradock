#include "dynamics/two_body.hpp"
#include "math/angle.hpp"
#include "math/constants.hpp"
#include "math/vector3.hpp"
#include "numerics/fixed_step_propagation.hpp"
#include "orbit/cartesian_state.hpp"
#include "orbit/classical_elements.hpp"
#include "orbit/two_body_orbit.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using astradock::math::angular_distance_rad;
using astradock::constants::pi;
using astradock::numerics::IntegrationMethod;
using astradock::numerics::propagate_fixed_step;
using astradock::orbit::CartesianState;
using astradock::orbit::classical_elements_to_state;
using astradock::orbit::ClassicalOrbitalElements;
using astradock::orbit::state_to_classical_elements;
using astradock::orbit::two_body_state_derivative;

int main() {
    std::cout << "=================================================================\n";
    std::cout << " AstraDock — M07 Classical Orbital Elements Demo & Propagation \n";
    std::cout << "=================================================================\n\n";

    const double mu = astradock::constants::earth_gravitational_parameter_m3_per_s2;

    // Nominal Case D: Inclined Elliptical Orbit
    const ClassicalOrbitalElements initial_elem{
        10'000'000.0,            // semi_major_axis_m: 10,000 km
        0.2,                     // eccentricity: 0.2
        45.0 * pi / 180.0,       // inclination_rad: 45 deg
        120.0 * pi / 180.0,      // raan_rad: 120 deg
        60.0 * pi / 180.0,       // argument_of_periapsis_rad: 60 deg
        30.0 * pi / 180.0,       // true_anomaly_rad: 30 deg
    };

    const double period_s = astradock::orbit::orbital_period_s(initial_elem.semi_major_axis_m, mu);
    const CartesianState state0 = classical_elements_to_state(initial_elem, mu);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Initial Orbit Definition (Case D — Inclined Elliptic):\n";
    std::cout << "  Semi-Major Axis (a):         " << initial_elem.semi_major_axis_m / 1000.0 << " km\n";
    std::cout << "  Eccentricity (e):            " << initial_elem.eccentricity << "\n";
    std::cout << "  Inclination (i):             " << initial_elem.inclination_rad * 180.0 / pi << " deg\n";
    std::cout << "  RAAN (Omega):                " << initial_elem.raan_rad * 180.0 / pi << " deg\n";
    std::cout << "  Arg of Periapsis (omega):    " << initial_elem.argument_of_periapsis_rad * 180.0 / pi << " deg\n";
    std::cout << "  True Anomaly (nu):           " << initial_elem.true_anomaly_rad * 180.0 / pi << " deg\n";
    std::cout << "  Analytical Orbital Period:   " << period_s << " s (" << period_s / 60.0 << " min)\n\n";

    std::cout << "Initial Cartesian State in ECI:\n";
    std::cout << "  Position (r):                [" << state0.position.x() << ", " << state0.position.y() << ", " << state0.position.z() << "] m\n";
    std::cout << "  Velocity (v):                [" << state0.velocity.x() << ", " << state0.velocity.y() << ", " << state0.velocity.z() << "] m/s\n";
    std::cout << "  Radius (|r|):                " << state0.position.norm() / 1000.0 << " km\n";
    std::cout << "  Speed (|v|):                 " << state0.velocity.norm() << " m/s\n\n";

    // Propagate for 3 full orbital periods with dt = 10.0 s
    const double dt = 10.0;
    const double duration = 3.0 * period_s;

    std::cout << "Propagating 3 orbits (" << duration << " s, dt = " << dt << " s) with classical RK4...\n";
    const auto trajectory = propagate_fixed_step(
        0.0, duration, dt, state0,
        IntegrationMethod::classical_rk4,
        [mu](double t, const CartesianState& s) {
            return two_body_state_derivative(t, s, mu);
        });

    std::cout << "Generated " << trajectory.size() << " trajectory samples.\n\n";

    // Export to CSV
    std::filesystem::create_directories("artifacts/data");
    const std::string csv_path = "artifacts/data/m07_elements_propagation.csv";
    std::ofstream csv(csv_path);
    if (!csv.is_open()) {
        std::cerr << "Failed to open CSV file for writing: " << csv_path << "\n";
        return 1;
    }

    csv << "time_s,x_m,y_m,z_m,vx_mps,vy_mps,vz_mps,a_m,eccentricity,inclination_rad,raan_rad,argument_of_periapsis_rad,true_anomaly_rad,specific_energy_m2_s2,specific_angular_momentum_m2_s,roundtrip_pos_err_m,roundtrip_vel_err_mps\n";

    double max_a_drift_m = 0.0;
    double max_e_drift = 0.0;
    double max_i_drift_rad = 0.0;
    double max_raan_drift_rad = 0.0;
    double max_argp_drift_rad = 0.0;
    double max_roundtrip_pos_err_m = 0.0;
    double max_roundtrip_vel_err_mps = 0.0;

    for (const auto& sample : trajectory) {
        const ClassicalOrbitalElements elem = state_to_classical_elements(sample.state, mu);
        const CartesianState reconstructed = classical_elements_to_state(elem, mu);

        const double energy = astradock::orbit::specific_orbital_energy_m2_per_s2(sample.state, mu);
        const double h_mag = astradock::orbit::specific_angular_momentum_m2_per_s(sample.state).norm();

        const double pos_err = (reconstructed.position - sample.state.position).norm();
        const double vel_err = (reconstructed.velocity - sample.state.velocity).norm();

        max_a_drift_m = std::max(max_a_drift_m, std::abs(elem.semi_major_axis_m - initial_elem.semi_major_axis_m));
        max_e_drift = std::max(max_e_drift, std::abs(elem.eccentricity - initial_elem.eccentricity));
        max_i_drift_rad = std::max(max_i_drift_rad, angular_distance_rad(elem.inclination_rad, initial_elem.inclination_rad));
        max_raan_drift_rad = std::max(max_raan_drift_rad, angular_distance_rad(elem.raan_rad, initial_elem.raan_rad));
        max_argp_drift_rad = std::max(max_argp_drift_rad, angular_distance_rad(elem.argument_of_periapsis_rad, initial_elem.argument_of_periapsis_rad));

        max_roundtrip_pos_err_m = std::max(max_roundtrip_pos_err_m, pos_err);
        max_roundtrip_vel_err_mps = std::max(max_roundtrip_vel_err_mps, vel_err);

        csv << std::setprecision(10)
            << sample.time_s << ","
            << sample.state.position.x() << ","
            << sample.state.position.y() << ","
            << sample.state.position.z() << ","
            << sample.state.velocity.x() << ","
            << sample.state.velocity.y() << ","
            << sample.state.velocity.z() << ","
            << elem.semi_major_axis_m << ","
            << elem.eccentricity << ","
            << elem.inclination_rad << ","
            << elem.raan_rad << ","
            << elem.argument_of_periapsis_rad << ","
            << elem.true_anomaly_rad << ","
            << energy << ","
            << h_mag << ","
            << pos_err << ","
            << vel_err << "\n";
    }

    csv.close();
    std::cout << "Saved trajectory to: " << csv_path << "\n\n";

    std::cout << "Multi-Orbit Element Invariance & Round-Trip Summary (3 Orbits):\n";
    std::cout << "  Max Semi-Major Axis Drift:   " << std::scientific << max_a_drift_m << " m\n";
    std::cout << "  Max Eccentricity Drift:       " << max_e_drift << "\n";
    std::cout << "  Max Inclination Drift:        " << max_i_drift_rad << " rad (" << max_i_drift_rad * 180.0 / pi << " deg)\n";
    std::cout << "  Max RAAN Drift:               " << max_raan_drift_rad << " rad (" << max_raan_drift_rad * 180.0 / pi << " deg)\n";
    std::cout << "  Max Arg of Periapsis Drift:   " << max_argp_drift_rad << " rad (" << max_argp_drift_rad * 180.0 / pi << " deg)\n";
    std::cout << "  Max Round-Trip Position Err:  " << max_roundtrip_pos_err_m << " m\n";
    std::cout << "  Max Round-Trip Velocity Err:  " << max_roundtrip_vel_err_mps << " m/s\n\n";

    std::cout << "=================================================================\n";
    std::cout << " Demo execution complete. \n";
    std::cout << "=================================================================\n";

    return 0;
}
