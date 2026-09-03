#include "actuators/actuator_assembly.hpp"
#include "actuators/reaction_wheel.hpp"
#include "actuators/thruster.hpp"
#include "attitude/rigid_body.hpp"
#include "math/constants.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"
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

using namespace astradock;

void write_actuator_csv(
    const std::filesystem::path& filepath,
    const std::vector<double>& time_s,
    const std::vector<double>& wheel_speed,
    const std::vector<double>& wheel_momentum,
    const std::vector<double>& cmd_torque,
    const std::vector<double>& achieved_torque,
    const std::vector<double>& reaction_torque,
    const std::vector<int>& torque_sat,
    const std::vector<int>& speed_sat,
    const std::vector<int>& authority_lost) {
    std::filesystem::create_directories(filepath.parent_path());
    std::ofstream out(filepath);
    if (!out.is_open()) {
        std::cerr << "Failed to open CSV output file: " << filepath.string() << "\n";
        return;
    }
    out << "time_s,wheel_speed_rad_s,stored_momentum_Nms,commanded_torque_Nm,"
           "achieved_torque_Nm,reaction_torque_body_Nm,torque_saturated,"
           "speed_saturated,authority_lost\n";
    out << std::setprecision(14);
    for (std::size_t i = 0; i < time_s.size(); ++i) {
        out << time_s[i] << "," << wheel_speed[i] << "," << wheel_momentum[i] << ","
            << cmd_torque[i] << "," << achieved_torque[i] << "," << reaction_torque[i] << ","
            << torque_sat[i] << "," << speed_sat[i] << "," << authority_lost[i] << "\n";
    }
}

void write_thruster_csv(
    const std::filesystem::path& filepath,
    const std::vector<double>& time_s,
    const std::vector<double>& thrust,
    const std::vector<double>& force_x,
    const std::vector<double>& torque_z,
    const std::vector<double>& mass,
    const std::vector<double>& used) {
    std::filesystem::create_directories(filepath.parent_path());
    std::ofstream out(filepath);
    if (!out.is_open()) {
        std::cerr << "Failed to open CSV output file: " << filepath.string() << "\n";
        return;
    }
    out << "time_s,thrust_achieved_N,force_body_x_N,torque_body_z_Nm,"
           "spacecraft_mass_kg,propellant_used_kg\n";
    out << std::setprecision(14);
    for (std::size_t i = 0; i < time_s.size(); ++i) {
        out << time_s[i] << "," << thrust[i] << "," << force_x[i] << "," << torque_z[i]
            << "," << mass[i] << "," << used[i] << "\n";
    }
}

}  // namespace

int main() {
    using namespace astradock;
    using namespace astradock::actuators;
    using namespace astradock::math;

    std::cout << "============================================================\n";
    std::cout << " AstraDock — M14 Spacecraft Actuator Dynamics Demo          \n";
    std::cout << "============================================================\n";

    // Scenario 1: deterministic momentum buildup to saturation (M14B).
    // I = 1e-3 kg*m^2, max speed 100 rad/s -> 0.1 Nms capacity; 0.05 Nm -> 2 s.
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Scenario 1: Momentum buildup to wheel saturation\n";
    ReactionWheelParameters wheel{1.0e-3, 0.05, 100.0, 0.0, Vector3{0.0, 0.0, 1.0}};
    ReactionWheelState wheel_state{};
    const double dt = 0.05;
    const double duration = 3.0;
    std::vector<double> t_v, w_v, h_v, cmd_v, ach_v, react_v;
    std::vector<int> tsat_v, ssat_v, alost_v;
    double t = 0.0;
    while (t <= duration + 1.0e-12) {
        const double cmd = (t < 2.5) ? 0.05 : -0.05;  // saturate, then recover
        const ReactionWheelStep step = step_reaction_wheel(wheel, wheel_state, cmd, dt);
        t_v.push_back(t);
        w_v.push_back(step.state.wheel_speed_rad_s);
        h_v.push_back(step.stored_momentum_Nms);
        cmd_v.push_back(cmd);
        ach_v.push_back(step.state.achieved_motor_torque_Nm);
        react_v.push_back(step.reaction_torque_body_Nm.z());
        tsat_v.push_back(step.torque_saturated ? 1 : 0);
        ssat_v.push_back(step.speed_saturated ? 1 : 0);
        alost_v.push_back(step.authority_lost ? 1 : 0);
        wheel_state = step.state;
        t += dt;
    }
    write_actuator_csv(
        "data/m14_wheel_saturation.csv", t_v, w_v, h_v, cmd_v, ach_v, react_v, tsat_v, ssat_v, alost_v);
    std::cout << "  Final wheel speed: " << wheel_state.wheel_speed_rad_s << " rad/s\n";
    std::cout << "  Exported to: data/m14_wheel_saturation.csv\n";

    // Scenario 2: offset thruster burn with propellant bookkeeping (M14C/M14D).
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Scenario 2: Offset thruster burn with mass depletion\n";
    ThrusterParameters thruster{2.0, 220.0, Vector3{0.0, 0.5, 0.0}, Vector3{1.0, 0.0, 0.0}};
    const Quaternion attitude = Quaternion::identity();
    double wet_mass = 500.0;
    const double dry_mass = 400.0;
    const double burn_dt = 0.5;
    const double burn_duration = 20.0;
    std::vector<double> bt_v, th_v, fx_v, tz_v, m_v, u_v;
    double bt = 0.0;
    double cumulative_used = 0.0;
    while (bt <= burn_duration + 1.0e-12) {
        const double cmd = (bt < 10.0) ? 2.0 : 0.0;  // 10 s burn, then coast
        const ThrusterStep step = step_thruster(thruster, attitude, cmd, wet_mass, dry_mass, burn_dt);
        bt_v.push_back(bt);
        th_v.push_back(step.thrust_achieved_N);
        fx_v.push_back(step.force_body_N.x());
        tz_v.push_back(step.torque_body_Nm.z());
        wet_mass = step.spacecraft_mass_after_kg;
        cumulative_used += step.propellant_used_kg;
        m_v.push_back(wet_mass);
        u_v.push_back(cumulative_used);
        bt += burn_dt;
    }
    write_thruster_csv("data/m14_thruster_burn.csv", bt_v, th_v, fx_v, tz_v, m_v, u_v);
    std::cout << "  Propellant used: " << cumulative_used << " kg\n";
    std::cout << "  Final wet mass: " << wet_mass << " kg\n";
    std::cout << "  Exported to: data/m14_thruster_burn.csv\n";

    // Scenario 3: wheel-only torque on the 6-DOF bus — angular momentum exchange.
    // A +Z wheel torque must spin the bus -Z; total (bus + wheel) momentum conserved.
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Scenario 3: Wheel/bus angular momentum exchange on 6-DOF truth\n";
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const double r_orbit = constants::earth_reference_radius_m + 500.0e3;
    const double v_circ = std::sqrt(mu / r_orbit);
    spacecraft::SpacecraftParameters params{mu, attitude::PrincipalInertia{10.0, 10.0, 10.0}};
    spacecraft::SpacecraftState bus_state{
        orbit::CartesianState{Vector3{r_orbit, 0.0, 0.0}, Vector3{0.0, v_circ, 0.0}},
        attitude::RotationalState{Quaternion::identity(), Vector3{0.0, 0.0, 0.0}}};
    ReactionWheelParameters bus_wheel{0.01, 1.0, 1.0e6, 0.0, Vector3{0.0, 0.0, 1.0}};
    ReactionWheelState bus_wheel_state{};
    const double bus_dt = 0.1;
    const double bus_duration = 5.0;
    const double wheel_cmd = 0.1;
    double bus_t = 0.0;
    // Total inertial Z momentum at t=0 is zero (bus at rest, wheel at rest).
    while (bus_t < bus_duration - 1.0e-12) {
        const ReactionWheelStep wstep = step_reaction_wheel(bus_wheel, bus_wheel_state, wheel_cmd, bus_dt);
        bus_wheel_state = wstep.state;
        const spacecraft::ForceTorqueInput input{Vector3{}, wstep.reaction_torque_body_Nm};
        bus_state = spacecraft::rk4_step_spacecraft(bus_t, bus_state, bus_dt, params, input);
        bus_t += bus_dt;
    }
    const Vector3 h_bus_body =
        attitude::body_angular_momentum_kg_m2_per_s(bus_state.angular_velocity_rad_per_s(), params.inertia);
    const Vector3 h_bus_eci = bus_state.orientation().rotate_vector(h_bus_body);
    const double h_wheel_scalar = wheel_momentum_scalar_Nms(bus_wheel, bus_wheel_state.wheel_speed_rad_s);
    // Wheel axis is body +Z; express wheel momentum in ECI through the bus attitude.
    const Vector3 h_wheel_eci =
        bus_state.orientation().rotate_vector(Vector3{0.0, 0.0, h_wheel_scalar});
    const double total_z = h_bus_eci.z() + h_wheel_eci.z();
    std::cout << "  Bus omega_z: " << bus_state.angular_velocity_rad_per_s().z() << " rad/s\n";
    std::cout << "  Wheel speed: " << bus_wheel_state.wheel_speed_rad_s << " rad/s\n";
    std::cout << "  Total inertial Z momentum (bus + wheel): " << total_z << " Nms\n";

    std::cout << "============================================================\n";
    std::cout << " All M14 actuator demonstrations completed successfully.\n";
    std::cout << "============================================================\n";
    return 0;
}
