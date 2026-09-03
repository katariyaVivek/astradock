#include "attitude/principal_inertia.hpp"
#include "attitude/rigid_body.hpp"
#include "control/attitude_pd.hpp"
#include "control/control_saturation.hpp"
#include "control/relative_pd.hpp"
#include "math/constants.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"
#include "relative/relative_state.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

using namespace astradock;

void write_attitude_csv(
    const std::filesystem::path& filepath,
    const std::vector<double>& time_s,
    const std::vector<double>& err_deg,
    const std::vector<double>& rate_dps,
    const std::vector<double>& torque_Nm,
    const std::vector<double>& saturated) {
    std::filesystem::create_directories(filepath.parent_path());
    std::ofstream out(filepath);
    if (!out.is_open()) {
        std::cerr << "Failed to open CSV output file: " << filepath.string() << "\n";
        return;
    }
    out << "time_s,attitude_error_deg,rate_norm_deg_s,torque_norm_Nm,saturated\n";
    out << std::setprecision(14);
    for (std::size_t i = 0; i < time_s.size(); ++i) {
        out << time_s[i] << "," << err_deg[i] << "," << rate_dps[i] << "," << torque_Nm[i] << ","
            << saturated[i] << "\n";
    }
}

void write_translation_csv(
    const std::filesystem::path& filepath,
    const std::vector<double>& time_s,
    const std::vector<double>& pos_m,
    const std::vector<double>& vel_mps,
    const std::vector<double>& accel_mps2) {
    std::filesystem::create_directories(filepath.parent_path());
    std::ofstream out(filepath);
    if (!out.is_open()) {
        std::cerr << "Failed to open CSV output file: " << filepath.string() << "\n";
        return;
    }
    out << "time_s,position_norm_m,velocity_norm_mps,accel_norm_mps2\n";
    out << std::setprecision(14);
    for (std::size_t i = 0; i < time_s.size(); ++i) {
        out << time_s[i] << "," << pos_m[i] << "," << vel_mps[i] << "," << accel_mps2[i] << "\n";
    }
}

}  // namespace

int main() {
    using namespace astradock;
    using namespace astradock::control;
    using namespace astradock::math;

    std::cout << "============================================================\n";
    std::cout << " AstraDock — M16 Guidance & Control Foundations Demo        \n";
    std::cout << "============================================================\n";

    // Scenario 1: detumble the master-spec tumble + hold identity attitude.
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Scenario 1: detumble omega = [0.2, -0.1, 0.15] rad/s\n";
    attitude::PrincipalInertia inertia{10.0, 20.0, 30.0};
    attitude::RotationalState rot{
        Quaternion::identity(), Vector3{0.2, -0.1, 0.15}};
    const AttitudePdGains gains =
        suggest_pd_gains(PdTuningRequest{{10.0, 20.0, 30.0}, 0.9, 30.0});
    const Quaternion q_hold = Quaternion::identity();
    const Vector3 w_hold{};
    const Vector3 torque_limit{0.5, 0.5, 0.5};  // M16E accounting demo bound
    const double dt = 0.05;
    const double duration = 100.0;
    std::vector<double> t_v, err_v, rate_v, tq_v, sat_v;
    double t = 0.0;
    while (t <= duration + 1.0e-12) {
        const Vector3 desired = attitude_pd_torque_body_Nm(
            rot.orientation, q_hold, rot.angular_velocity_rad_per_s, w_hold, gains);
        const SaturatedControl sat = saturate_control_vector(desired, torque_limit);
        rot = attitude::rk4_step_rotational(t, rot, dt, inertia, sat.achieved);
        t += dt;
        t_v.push_back(t);
        err_v.push_back(attitude_error_angle_rad(rot.orientation, q_hold) * 180.0 / constants::pi);
        rate_v.push_back(
            rot.angular_velocity_rad_per_s.norm() * 180.0 / constants::pi);
        tq_v.push_back(sat.achieved.norm());
        sat_v.push_back(sat.saturated ? 1.0 : 0.0);
    }
    write_attitude_csv("data/m16_attitude_hold.csv", t_v, err_v, rate_v, tq_v, sat_v);
    std::cout << "  Final error: " << err_v.back() << " deg, final rate: " << rate_v.back()
              << " deg/s\n";
    std::cout << "  Exported to: data/m16_attitude_hold.csv\n";

    // Scenario 2: 90-deg slew about +Z with the same tuned law.
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Scenario 2: 90-deg yaw slew\n";
    attitude::RotationalState slew{Quaternion::identity(), Vector3{}};
    const Quaternion q_slew =
        Quaternion::from_axis_angle(Vector3{0.0, 0.0, 1.0}, constants::pi / 2.0);
    t = 0.0;
    std::vector<double> st_v, serr_v, srate_v, stq_v, ssat_v;
    while (t <= duration + 1.0e-12) {
        const Vector3 desired = attitude_pd_torque_body_Nm(
            slew.orientation, q_slew, slew.angular_velocity_rad_per_s, w_hold, gains);
        const SaturatedControl sat = saturate_control_vector(desired, torque_limit);
        slew = attitude::rk4_step_rotational(t, slew, dt, inertia, sat.achieved);
        t += dt;
        st_v.push_back(t);
        serr_v.push_back(attitude_error_angle_rad(slew.orientation, q_slew) * 180.0 / constants::pi);
        srate_v.push_back(slew.angular_velocity_rad_per_s.norm() * 180.0 / constants::pi);
        stq_v.push_back(sat.achieved.norm());
        ssat_v.push_back(sat.saturated ? 1.0 : 0.0);
    }
    write_attitude_csv("data/m16_attitude_slew.csv", st_v, serr_v, srate_v, stq_v, ssat_v);
    std::cout << "  Final error: " << serr_v.back() << " deg\n";
    std::cout << "  Exported to: data/m16_attitude_slew.csv\n";

    // Scenario 3: translation station keeping from 100 m radial offset.
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Scenario 3: 100 m radial station keeping (CW plant)\n";
    const double n = 1.1e-3;
    const RelativePdGains tp_gains{{1.0e-5, 1.0e-5, 1.0e-5}, {6.0e-3, 6.0e-3, 6.0e-3}};
    Vector3 rho{100.0, 0.0, 0.0};
    Vector3 vel{};
    const double tdt = 1.0;
    std::vector<double> ct_v, cp_v, cv_v, ca_v;
    t = 0.0;
    for (int i = 0; i <= 3600; ++i) {
        const Vector3 a_cmd = relative_pd_accel_lvlh_mps2(
            rho, vel, Vector3{}, Vector3{}, Vector3{}, n, tp_gains);
        ct_v.push_back(t);
        cp_v.push_back(rho.norm());
        cv_v.push_back(vel.norm());
        ca_v.push_back(a_cmd.norm());
        const relative::RelativeStateLvlh s{rho, vel};
        const auto unforced = relative::cw_predict(s, n, tdt);
        vel = unforced.relative_velocity_lvlh_mps + a_cmd * tdt;
        rho = unforced.relative_position_lvlh_m + a_cmd * (0.5 * tdt * tdt);
        t += tdt;
    }
    write_translation_csv("data/m16_station_keeping.csv", ct_v, cp_v, cv_v, ca_v);
    std::cout << "  Final offset: " << cp_v.back() << " m\n";
    std::cout << "  Exported to: data/m16_station_keeping.csv\n";

    std::cout << "============================================================\n";
    std::cout << " All M16 control demonstrations completed successfully.\n";
    std::cout << "============================================================\n";
    return 0;
}
