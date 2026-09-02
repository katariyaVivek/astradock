#include "environment/environment_models.hpp"
#include "math/constants.hpp"
#include "orbit/classical_elements.hpp"
#include "orbit/two_body_orbit.hpp"
#include "spacecraft/six_dof_dynamics.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace astradock;

namespace {

void ensure_data_directory_exists(const std::filesystem::path& dir_path) {
    if (!std::filesystem::exists(dir_path)) {
        std::filesystem::create_directories(dir_path);
    }
}

void export_j2_scenario(const std::filesystem::path& output_path) {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const attitude::PrincipalInertia inertia{10.0, 10.0, 10.0};
    const spacecraft::SpacecraftParameters params{mu, inertia};

    // Orbit: a = 10,000 km, e = 0.2, i = 45 deg, RAAN = 120 deg, omega = 60 deg, nu = 30 deg
    const orbit::ClassicalOrbitalElements el0{
        10000.0e3,
        0.2,
        45.0 * constants::pi / 180.0,
        120.0 * constants::pi / 180.0,
        60.0 * constants::pi / 180.0,
        30.0 * constants::pi / 180.0
    };
    const orbit::CartesianState cart0 = orbit::classical_elements_to_state(el0, mu);
    const spacecraft::SpacecraftState initial_state{cart0, attitude::RotationalState{}};

    const double period = 2.0 * constants::pi * std::sqrt(std::pow(el0.semi_major_axis_m, 3) / mu);
    const double sim_duration = 5.0 * period; // 5 orbital periods (~15.6 hours)
    const double dt = 5.0; // 5 second timestep

    environment::EnvironmentConfiguration config_j2;
    config_j2.enable_j2 = true;
    const environment::EnvironmentalParameters env_params;

    const auto traj = environment::propagate_spacecraft_environmental(
        initial_state, sim_duration, dt, params, env_params, config_j2
    );

    std::ofstream file(output_path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open " + output_path.string() + " for writing");
    }

    file << "time_s,pos_x_m,pos_y_m,pos_z_m,vel_x_mps,vel_y_mps,vel_z_mps,"
         << "semi_major_axis_m,eccentricity,inclination_deg,raan_deg,argument_of_periapsis_deg,true_anomaly_deg\n";
    file << std::fixed << std::setprecision(8);

    for (const auto& sample : traj) {
        const auto cart = sample.state.translational;
        const auto el = orbit::state_to_classical_elements(cart, mu);

        file << sample.time_s << ","
             << cart.position.x() << "," << cart.position.y() << "," << cart.position.z() << ","
             << cart.velocity.x() << "," << cart.velocity.y() << "," << cart.velocity.z() << ","
             << el.semi_major_axis_m << ","
             << el.eccentricity << ","
             << (el.inclination_rad * 180.0 / constants::pi) << ","
             << (el.raan_rad * 180.0 / constants::pi) << ","
             << (el.argument_of_periapsis_rad * 180.0 / constants::pi) << ","
             << (el.true_anomaly_rad * 180.0 / constants::pi) << "\n";
    }

    std::cout << "Exported J2 scenario to: " << output_path << " (" << traj.size() << " samples)\n";
}

void export_drag_scenario(const std::filesystem::path& output_path) {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const attitude::PrincipalInertia inertia{10.0, 10.0, 10.0};
    const spacecraft::SpacecraftParameters params{mu, inertia};

    // Circular orbit at 300 km altitude
    const double alt = 300.0e3;
    const double r_orb = constants::earth_reference_radius_m + alt;
    const double v_circ = std::sqrt(mu / r_orb);

    const spacecraft::SpacecraftState initial_state{
        orbit::CartesianState{math::Vector3{r_orb, 0.0, 0.0}, math::Vector3{0.0, v_circ, 0.0}},
        attitude::RotationalState{}
    };

    const double period = 2.0 * constants::pi * std::sqrt(std::pow(r_orb, 3) / mu);
    const double sim_duration = 10.0 * period; // 10 orbits (~15 hours)
    const double dt = 2.0;

    environment::EnvironmentConfiguration config_drag;
    config_drag.enable_drag = true;

    environment::EnvironmentalParameters env_params;
    env_params.mass_kg = 250.0;
    env_params.drag_reference_area_m2 = 3.0;
    env_params.drag_coefficient_cd = 2.2;
    env_params.atmosphere_ref_alt_m = 300.0e3;
    env_params.atmosphere_ref_density_kg_per_m3 = 2.4e-11;
    env_params.atmosphere_scale_height_m = 45.0e3;

    const auto traj = environment::propagate_spacecraft_environmental(
        initial_state, sim_duration, dt, params, env_params, config_drag
    );

    std::ofstream file(output_path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open " + output_path.string() + " for writing");
    }

    file << "time_s,altitude_km,semi_major_axis_km,specific_energy_m2_s2,drag_accel_mps2\n";
    file << std::fixed << std::setprecision(8);

    for (const auto& sample : traj) {
        const double r = sample.state.position().norm();
        const double current_alt_km = (r - constants::earth_reference_radius_m) / 1000.0;
        const double energy = orbit::specific_orbital_energy_m2_per_s2(sample.state.translational, mu);
        const double sma_km = -mu / (2.0 * energy) / 1000.0;

        const math::Vector3 a_drag = environment::drag_acceleration_eci(
            sample.state.velocity(),
            sample.state.position(),
            env_params.mass_kg,
            env_params.drag_coefficient_cd,
            env_params.drag_reference_area_m2,
            env_params.earth_reference_radius_m,
            env_params.earth_rotation_rate_rad_per_s,
            env_params.atmosphere_ref_alt_m,
            env_params.atmosphere_ref_density_kg_per_m3,
            env_params.atmosphere_scale_height_m
        );

        file << sample.time_s << ","
             << current_alt_km << ","
             << sma_km << ","
             << energy << ","
             << a_drag.norm() << "\n";
    }

    std::cout << "Exported Drag scenario to: " << output_path << " (" << traj.size() << " samples)\n";
}

void export_third_body_scenario(const std::filesystem::path& output_path) {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const attitude::PrincipalInertia inertia{10.0, 10.0, 10.0};
    const spacecraft::SpacecraftParameters params{mu, inertia};

    // Geostationary orbit (r = 42,164 km)
    const double r_geo = 42164.0e3;
    const double v_geo = std::sqrt(mu / r_geo);

    const spacecraft::SpacecraftState initial_state{
        orbit::CartesianState{math::Vector3{r_geo, 0.0, 0.0}, math::Vector3{0.0, v_geo, 0.0}},
        attitude::RotationalState{}
    };

    const double period = 86164.0905; // 1 sidereal day
    const double sim_duration = 3.0 * period; // 3 days
    const double dt = 30.0; // 30 second step

    environment::EnvironmentConfiguration config_tb;
    config_tb.enable_third_body = true;
    const environment::EnvironmentalParameters env_params;

    const auto traj = environment::propagate_spacecraft_environmental(
        initial_state, sim_duration, dt, params, env_params, config_tb
    );

    std::ofstream file(output_path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open " + output_path.string() + " for writing");
    }

    file << "time_s,pos_x_m,pos_y_m,pos_z_m,a_3b_x_mps2,a_3b_y_mps2,a_3b_z_mps2,a_3b_mag_mps2\n";
    file << std::scientific << std::setprecision(12);

    for (const auto& sample : traj) {
        const math::Vector3 a_3b = environment::third_body_acceleration_eci(
            sample.state.position(),
            env_params.third_body_position_eci_m,
            env_params.third_body_mu_m3_per_s2
        );

        file << sample.time_s << ","
             << sample.state.position().x() << ","
             << sample.state.position().y() << ","
             << sample.state.position().z() << ","
             << a_3b.x() << "," << a_3b.y() << "," << a_3b.z() << ","
             << a_3b.norm() << "\n";
    }

    std::cout << "Exported Third-Body scenario to: " << output_path << " (" << traj.size() << " samples)\n";
}

void export_gravity_gradient_scenario(const std::filesystem::path& output_path) {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    // Dumbbell satellite (Ixx = 10, Iyy = 50, Izz = 50 kg*m^2)
    const attitude::PrincipalInertia inertia{10.0, 50.0, 50.0};
    const spacecraft::SpacecraftParameters params{mu, inertia};

    const double r_orb = constants::earth_reference_radius_m + 500.0e3;
    const double v_circ = std::sqrt(mu / r_orb);

    // Initial pitch offset of 10 degrees about +Y axis
    const double pitch_deg0 = 10.0;
    const double pitch_rad0 = pitch_deg0 * constants::pi / 180.0;
    const math::Quaternion q0{std::cos(pitch_rad0 / 2.0), 0.0, std::sin(pitch_rad0 / 2.0), 0.0};

    const spacecraft::SpacecraftState initial_state{
        orbit::CartesianState{math::Vector3{r_orb, 0.0, 0.0}, math::Vector3{0.0, v_circ, 0.0}},
        attitude::RotationalState{q0, math::Vector3{0.0, 0.0, 0.0}}
    };

    const double sim_duration = 3600.0; // 1 hour
    const double dt = 0.5;

    environment::EnvironmentConfiguration config_gg;
    config_gg.enable_gravity_gradient = true;
    const environment::EnvironmentalParameters env_params;

    const auto traj = environment::propagate_spacecraft_environmental(
        initial_state, sim_duration, dt, params, env_params, config_gg
    );

    std::ofstream file(output_path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open " + output_path.string() + " for writing");
    }

    file << "time_s,pitch_deg,omega_y_rad_s,torque_y_Nm,torque_mag_Nm\n";
    file << std::fixed << std::setprecision(8);

    for (const auto& sample : traj) {
        // Compute current pitch angle from body radial vector
        const math::Vector3 r_hat_eci = sample.state.position().normalized();
        const math::Vector3 r_hat_body = sample.state.orientation().conjugate().rotate_vector(r_hat_eci);
        const double pitch_deg = std::atan2(-r_hat_body.z(), r_hat_body.x()) * 180.0 / constants::pi;

        const math::Vector3 tau_gg = environment::gravity_gradient_torque_body(
            sample.state.position(),
            sample.state.orientation(),
            params.inertia,
            mu
        );

        file << sample.time_s << ","
             << pitch_deg << ","
             << sample.state.angular_velocity_rad_per_s().y() << ","
             << tau_gg.y() << ","
             << tau_gg.norm() << "\n";
    }

    std::cout << "Exported Gravity-Gradient scenario to: " << output_path << " (" << traj.size() << " samples)\n";
}

} // namespace

int main() {
    try {
        const std::filesystem::path data_dir = std::filesystem::path("data");
        ensure_data_directory_exists(data_dir);

        std::cout << "======================================================\n";
        std::cout << "  AstraDock Milestone M11 — Spacecraft Environment Demo\n";
        std::cout << "======================================================\n";

        export_j2_scenario(data_dir / "environment_j2_orbit.csv");
        export_drag_scenario(data_dir / "environment_drag_decay.csv");
        export_third_body_scenario(data_dir / "environment_third_body.csv");
        export_gravity_gradient_scenario(data_dir / "environment_gravity_gradient.csv");

        std::cout << "======================================================\n";
        std::cout << "  All M11 environmental scenarios exported successfully.\n";
        std::cout << "======================================================\n";

        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "FATAL ERROR: " << ex.what() << "\n";
        return 1;
    }
}
