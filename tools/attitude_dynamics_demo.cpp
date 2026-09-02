#include "attitude/principal_inertia.hpp"
#include "attitude/rigid_body.hpp"
#include "attitude/rotational_state.hpp"
#include "math/matrix3.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

using astradock::attitude::body_angular_momentum_kg_m2_per_s;
using astradock::attitude::inertial_angular_momentum_kg_m2_per_s;
using astradock::attitude::PrincipalInertia;
using astradock::attitude::rk4_step_rotational;
using astradock::attitude::rotational_kinetic_energy_J;
using astradock::attitude::RotationalState;
using astradock::math::Quaternion;
using astradock::math::Vector3;

void run_scenario(
    const std::string& scenario_name,
    const std::string& csv_filename,
    const RotationalState& initial_state,
    const PrincipalInertia& inertia,
    const Vector3& torque_body_Nm,
    double duration_s,
    double dt_s) {
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Scenario: " << scenario_name << "\n";
    std::cout << "Inertia: [" << inertia.Ixx_kg_m2 << ", " << inertia.Iyy_kg_m2
              << ", " << inertia.Izz_kg_m2 << "] kg*m^2\n";
    std::cout << "Torque:  [" << torque_body_Nm.x() << ", " << torque_body_Nm.y()
              << ", " << torque_body_Nm.z() << "] N*m\n";
    std::cout << "Initial omega: [" << initial_state.angular_velocity_rad_per_s.x() << ", "
              << initial_state.angular_velocity_rad_per_s.y() << ", "
              << initial_state.angular_velocity_rad_per_s.z() << "] rad/s\n";

    std::filesystem::create_directories("data");
    const std::string csv_path = "data/" + csv_filename;
    std::ofstream csv(csv_path);
    if (!csv.is_open()) {
        std::cerr << "Failed to open output CSV: " << csv_path << "\n";
        return;
    }

    csv << "time_s,qw,qx,qy,qz,wx_rad_s,wy_rad_s,wz_rad_s,"
        << "Hx_body_Nms,Hy_body_Nms,Hz_body_Nms,"
        << "Hx_inertial_Nms,Hy_inertial_Nms,Hz_inertial_Nms,"
        << "rotational_energy_J,quaternion_norm\n";

    RotationalState state = initial_state;
    const int total_steps = static_cast<int>(std::round(duration_s / dt_s));

    const double initial_energy = rotational_kinetic_energy_J(state.angular_velocity_rad_per_s, inertia);
    const Vector3 initial_h_inertial = inertial_angular_momentum_kg_m2_per_s(state, inertia);
    const double initial_h_norm = initial_h_inertial.norm();

    double max_rel_energy_err = 0.0;
    double max_h_inertial_err = 0.0;
    double max_norm_err = 0.0;

    for (int step = 0; step <= total_steps; ++step) {
        const double t = step * dt_s;

        const Vector3 h_body = body_angular_momentum_kg_m2_per_s(state.angular_velocity_rad_per_s, inertia);
        const Vector3 h_inertial = inertial_angular_momentum_kg_m2_per_s(state, inertia);
        const double energy = rotational_kinetic_energy_J(state.angular_velocity_rad_per_s, inertia);
        const double q_norm = state.orientation.norm();

        const double norm_err = std::abs(q_norm - 1.0);
        if (norm_err > max_norm_err) max_norm_err = norm_err;

        if (initial_energy > 0.0 && torque_body_Nm.norm() == 0.0) {
            const double rel_e_err = std::abs(energy - initial_energy) / initial_energy;
            if (rel_e_err > max_rel_energy_err) max_rel_energy_err = rel_e_err;
        }

        if (initial_h_norm > 0.0 && torque_body_Nm.norm() == 0.0) {
            const double h_err = (h_inertial - initial_h_inertial).norm() / initial_h_norm;
            if (h_err > max_h_inertial_err) max_h_inertial_err = h_err;
        }

        csv << std::setprecision(12)
            << t << ","
            << state.orientation.w() << ","
            << state.orientation.x() << ","
            << state.orientation.y() << ","
            << state.orientation.z() << ","
            << state.angular_velocity_rad_per_s.x() << ","
            << state.angular_velocity_rad_per_s.y() << ","
            << state.angular_velocity_rad_per_s.z() << ","
            << h_body.x() << "," << h_body.y() << "," << h_body.z() << ","
            << h_inertial.x() << "," << h_inertial.y() << "," << h_inertial.z() << ","
            << energy << ","
            << q_norm << "\n";

        if (step < total_steps) {
            state = rk4_step_rotational(t, state, dt_s, inertia, torque_body_Nm, true);
        }
    }

    csv.close();
    std::cout << "Simulation completed (" << total_steps << " steps).\n";
    std::cout << "Output exported to: " << csv_path << "\n";
    std::cout << "Max quaternion norm residual: " << max_norm_err << "\n";
    if (torque_body_Nm.norm() == 0.0 && initial_energy > 0.0) {
        std::cout << "Max relative energy drift:   " << max_rel_energy_err << "\n";
        std::cout << "Max inertial H vector drift: " << max_h_inertial_err << "\n";
    }
    std::cout << "Final omega: [" << state.angular_velocity_rad_per_s.x() << ", "
              << state.angular_velocity_rad_per_s.y() << ", "
              << state.angular_velocity_rad_per_s.z() << "] rad/s\n";
}

}  // namespace

int main() {
    std::cout << "============================================================\n";
    std::cout << " AstraDock — M09 Attitude Dynamics & Quaternion Kinematics  \n";
    std::cout << "============================================================\n";

    // 1. Principal-Axis Spin
    run_scenario(
        "Principal-Axis Spin (+Z Axis, Torque-Free)",
        "attitude_dynamics_principal_spin.csv",
        RotationalState{Quaternion::identity(), Vector3(0.0, 0.0, 0.5)},
        PrincipalInertia(10.0, 15.0, 20.0),
        Vector3{},
        20.0,
        0.01
    );

    // 2. Constant Torque About Principal X Axis
    run_scenario(
        "Constant Torque About Principal X-Axis",
        "attitude_dynamics_constant_torque.csv",
        RotationalState{Quaternion::identity(), Vector3{}},
        PrincipalInertia(10.0, 25.0, 40.0),
        Vector3(0.5, 0.0, 0.0),
        10.0,
        0.01
    );

    // 3. Torque-Free Asymmetric Tumbling
    run_scenario(
        "Torque-Free Asymmetric Rigid-Body Tumbling",
        "attitude_dynamics_asymmetric_tumble.csv",
        RotationalState{Quaternion::identity(), Vector3(0.2, 0.3, 0.1)},
        PrincipalInertia(10.0, 20.0, 30.0),
        Vector3{},
        30.0,
        0.01
    );

    std::cout << "============================================================\n";
    std::cout << " All attitude dynamics demonstrations executed successfully.\n";
    std::cout << "============================================================\n";

    return 0;
}
