#include "attitude/rigid_body.hpp"
#include "frames/frame_basis.hpp"
#include "frames/lvlh.hpp"
#include "math/constants.hpp"
#include "math/matrix3.hpp"
#include "orbit/classical_elements.hpp"
#include "orbit/orbital_diagnostics.hpp"
#include "orbit/two_body_orbit.hpp"
#include "spacecraft/force_torque.hpp"
#include "spacecraft/six_dof_dynamics.hpp"
#include "spacecraft/spacecraft_parameters.hpp"
#include "spacecraft/spacecraft_state.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

void write_six_dof_csv(
    const std::filesystem::path& filepath,
    const std::vector<astradock::numerics::StateSample<astradock::spacecraft::SpacecraftState>>& samples,
    const astradock::spacecraft::SpacecraftParameters& params) {
    std::filesystem::create_directories(filepath.parent_path());
    std::ofstream out(filepath);
    if (!out.is_open()) {
        std::cerr << "Failed to open CSV output file: " << filepath.string() << "\n";
        return;
    }

    out << "time_s,"
        << "position_eci_x_m,position_eci_y_m,position_eci_z_m,"
        << "velocity_eci_x_mps,velocity_eci_y_mps,velocity_eci_z_mps,"
        << "quaternion_w,quaternion_x,quaternion_y,quaternion_z,"
        << "angular_velocity_body_x_rad_s,angular_velocity_body_y_rad_s,angular_velocity_body_z_rad_s,"
        << "specific_orbital_energy_m2_s2,specific_angular_momentum_m2_s,"
        << "rotational_kinetic_energy_J,quaternion_norm\n";

    out << std::setprecision(14);

    for (const auto& sample : samples) {
        const auto& t = sample.time_s;
        const auto& state = sample.state;

        const double specific_energy = astradock::orbit::specific_orbital_energy_m2_per_s2(
            state.translational, params.gravitational_parameter_m3_per_s2);
        const double specific_h = astradock::orbit::specific_angular_momentum_m2_per_s(
            state.translational).norm();
        const double rot_energy = astradock::attitude::rotational_kinetic_energy_J(
            state.angular_velocity_rad_per_s(), params.inertia);
        const double q_norm = state.orientation().norm();

        out << t << ","
            << state.position().x() << "," << state.position().y() << "," << state.position().z() << ","
            << state.velocity().x() << "," << state.velocity().y() << "," << state.velocity().z() << ","
            << state.orientation().w() << "," << state.orientation().x() << ","
            << state.orientation().y() << "," << state.orientation().z() << ","
            << state.angular_velocity_rad_per_s().x() << ","
            << state.angular_velocity_rad_per_s().y() << ","
            << state.angular_velocity_rad_per_s().z() << ","
            << specific_energy << ","
            << specific_h << ","
            << rot_energy << ","
            << q_norm << "\n";
    }
}

}  // namespace

int main() {
    using namespace astradock;
    using namespace astradock::spacecraft;
    using namespace astradock::math;

    std::cout << "============================================================\n";
    std::cout << " AstraDock — M10 Integrated 6-DOF Spacecraft State Demo     \n";
    std::cout << "============================================================\n";

    // 1. Canonical 500 km Orbit Reference Setup
    const double r_orbit_m = constants::earth_reference_radius_m + 500.0e3; // 6,878,137 m
    const auto ref = orbit::compute_circular_orbit_reference(
        constants::earth_gravitational_parameter_m3_per_s2, r_orbit_m);

    const SpacecraftParameters params(
        ref.gravitational_parameter_m3_per_s2,
        attitude::PrincipalInertia(10.0, 20.0, 30.0) // Asymmetric inertia
    );

    // Initial 6-DOF state
    const Quaternion q_init = Quaternion::from_axis_angle(Vector3(0.0, 0.0, 1.0), constants::pi / 4.0); // 45 deg Z
    const Vector3 omega_init(0.05, 0.08, 0.02); // Multi-axis rate (rad/s)

    const SpacecraftState canonical_state0{
        orbit::CartesianState{
            Vector3(ref.orbital_radius_m, 0.0, 0.0),
            Vector3(0.0, ref.speed_m_per_s, 0.0)
        },
        attitude::RotationalState{
            q_init,
            omega_init
        }
    };

    // Scenario 1: Canonical Decoupled 6-DOF Motion (1 Full Orbit)
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Scenario 1: Canonical 500 km Circular Orbit + Tumbling Attitude\n";
    std::cout << "  Orbit Radius: " << ref.orbital_radius_m << " m (Altitude: 500 km)\n";
    std::cout << "  Orbital Speed:  " << ref.speed_m_per_s << " m/s\n";
    std::cout << "  Orbital Period: " << ref.period_s << " s (~" << ref.period_s / 60.0 << " min)\n";
    std::cout << "  Inertia: [" << params.inertia.Ixx() << ", " << params.inertia.Iyy() << ", " << params.inertia.Izz() << "] kg*m^2\n";
    std::cout << "  Initial Attitude: 45 deg about +Z\n";
    std::cout << "  Initial Body Rates: [" << omega_init.x() << ", " << omega_init.y() << ", " << omega_init.z() << "] rad/s\n";

    const double dt_canon = 1.0; // 1 second fixed step
    const auto canon_samples = propagate_spacecraft_fixed_step(
        0.0, ref.period_s, dt_canon, canonical_state0, params, ForceTorqueInput{});

    write_six_dof_csv("data/six_dof_canonical_orbit.csv", canon_samples, params);
    std::cout << "  Propagated " << canon_samples.size() << " samples.\n";
    std::cout << "  Exported to: data/six_dof_canonical_orbit.csv\n";

    // Invariant checks on canonical run
    const double init_orb_energy = orbit::specific_orbital_energy_m2_per_s2(
        canonical_state0.translational, params.gravitational_parameter_m3_per_s2);
    const double final_orb_energy = orbit::specific_orbital_energy_m2_per_s2(
        canon_samples.back().state.translational, params.gravitational_parameter_m3_per_s2);
    const double init_rot_energy = attitude::rotational_kinetic_energy_J(
        canonical_state0.angular_velocity_rad_per_s(), params.inertia);
    const double final_rot_energy = attitude::rotational_kinetic_energy_J(
        canon_samples.back().state.angular_velocity_rad_per_s(), params.inertia);

    std::cout << "  Relative Orbital Energy Drift:   "
              << std::abs(final_orb_energy - init_orb_energy) / std::abs(init_orb_energy) << "\n";
    std::cout << "  Relative Rotational Energy Drift: "
              << std::abs(final_rot_energy - init_rot_energy) / init_rot_energy << "\n";

    // Scenario 2: Constant External Torque (Attitude Driven, Orbit Undisturbed)
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Scenario 2: Constant Body Torque (+Z Axis: 0.2 N*m)\n";
    const SpacecraftState torque_state0{
        canonical_state0.translational,
        attitude::RotationalState{
            Quaternion::identity(),
            Vector3(0.0, 0.0, 0.0)
        }
    };
    const ForceTorqueInput torque_input(Vector3{}, Vector3(0.0, 0.0, 0.2));
    const double torque_duration = 200.0;
    const double dt_torque = 0.1;

    const auto torque_samples = propagate_spacecraft_fixed_step(
        0.0, torque_duration, dt_torque, torque_state0, params, torque_input);

    write_six_dof_csv("data/six_dof_constant_torque.csv", torque_samples, params);
    std::cout << "  Propagated " << torque_samples.size() << " samples.\n";
    std::cout << "  Final omega_z: " << torque_samples.back().state.angular_velocity_rad_per_s().z()
              << " rad/s (Analytical: " << (0.2 / params.inertia.Izz()) * torque_duration << " rad/s)\n";
    std::cout << "  Exported to: data/six_dof_constant_torque.csv\n";

    // Report Key Timestamps Frame Triad (Section 18)
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Integrated Frame Telemetry at Key Timestamps (ECI & Body):\n";
    const std::vector<double> report_fractions = {0.0, 0.25, 0.50, 0.75, 1.0};
    for (double frac : report_fractions) {
        const double t_target = frac * ref.period_s;
        // Find closest sample
        size_t idx = static_cast<size_t>(std::round(t_target / dt_canon));
        if (idx >= canon_samples.size()) idx = canon_samples.size() - 1;
        const auto& s = canon_samples[idx];
        const Matrix3 dcm_eci_from_b = s.state.orientation().to_rotation_matrix();

        std::cout << "  [t = " << std::fixed << std::setprecision(1) << s.time_s << " s (" << frac * 100 << "% T)]\n"
                  << "    Pos ECI (km): [" << s.state.position().x() / 1000.0 << ", "
                  << s.state.position().y() / 1000.0 << ", " << s.state.position().z() / 1000.0 << "]\n"
                  << "    Vel ECI (m/s): [" << s.state.velocity().x() << ", "
                  << s.state.velocity().y() << ", " << s.state.velocity().z() << "]\n"
                  << "    Body Rate (rad/s): [" << s.state.angular_velocity_rad_per_s().x() << ", "
                  << s.state.angular_velocity_rad_per_s().y() << ", " << s.state.angular_velocity_rad_per_s().z() << "]\n"
                  << "    DCM C_ECI_B col0: [" << dcm_eci_from_b(0, 0) << ", " << dcm_eci_from_b(1, 0) << ", " << dcm_eci_from_b(2, 0) << "]\n";
    }

    std::cout << "============================================================\n";
    std::cout << " All 6-DOF demonstrations completed successfully.\n";
    std::cout << "============================================================\n";

    return 0;
}
