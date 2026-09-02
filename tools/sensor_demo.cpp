#include "environment/atmospheric_drag.hpp"
#include "environment/environment_models.hpp"
#include "math/constants.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"
#include "sensors/gnss.hpp"
#include "sensors/imu.hpp"
#include "sensors/range_sensor.hpp"
#include "sensors/star_tracker.hpp"
#include "spacecraft/six_dof_dynamics.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using namespace astradock;

int main() {
    std::cout << "=================================================================\n";
    std::cout << " AstraDock - Milestone M12: Sensor Simulation & Telemetry Demo   \n";
    std::cout << "=================================================================\n\n";

    fs::create_directories("data");

    // Scenario Configuration:
    // 500 km LEO Orbit with atmospheric drag and tumbling spacecraft bus
    const double altitude_m = 500.0e3;
    const double r_orbit_m = constants::earth_reference_radius_m + altitude_m;
    const double v_circ_mps = std::sqrt(constants::earth_gravitational_parameter_m3_per_s2 / r_orbit_m);

    const spacecraft::SpacecraftState initial_state{
        orbit::CartesianState{
            math::Vector3{r_orbit_m, 0.0, 0.0},
            math::Vector3{0.0, v_circ_mps, 0.0}
        },
        attitude::RotationalState{
            math::Quaternion::identity(),
            math::Vector3{0.01, 0.02, 0.015} // Tumbling rates (rad/s)
        }
    };

    const attitude::PrincipalInertia inertia{100.0, 75.0, 50.0};
    const double mass_kg = 500.0;
    const spacecraft::SpacecraftParameters sc_params{
        constants::earth_gravitational_parameter_m3_per_s2,
        inertia
    };

    // Environmental Parameters (atmospheric drag)
    environment::EnvironmentalParameters env_params;
    env_params.mass_kg = mass_kg;
    env_params.drag_reference_area_m2 = 2.0;
    env_params.drag_coefficient_cd = 2.2;
    env_params.atmosphere_ref_density_kg_per_m3 = 4.0e-12; // 500 km
    env_params.atmosphere_ref_alt_m = 500.0e3;
    env_params.atmosphere_scale_height_m = 60.0e3;

    environment::EnvironmentConfiguration env_config;
    env_config.enable_drag = true;

    // Simulation Setup: 1 Orbit (~5677 s), integration step dt = 0.01 s (100 Hz)
    const double orbital_period_s = 2.0 * constants::pi * std::sqrt(r_orbit_m * r_orbit_m * r_orbit_m / constants::earth_gravitational_parameter_m3_per_s2);
    const double dt_s = 0.01;
    const double sim_duration_s = orbital_period_s;

    // Target Spacecraft for Range Sensor: circular orbit leading by ~10 km
    const double target_lead_angle_rad = 10.0e3 / r_orbit_m;

    // Configure Sensors
    const std::vector<sensors::DropoutWindow> dropouts = {{2000.0, 2100.0}};

    // 1. IMU: 100 Hz (dt = 0.01 s)
    sensors::ImuConfig imu_cfg;
    imu_cfg.sample_period_s = 0.01;
    imu_cfg.gyro_bias_rad_s = math::Vector3{0.0005, -0.0003, 0.0002};
    imu_cfg.gyro_noise_std_rad_s = 0.001; // 1 mrad/s
    imu_cfg.accel_bias_mps2 = math::Vector3{0.002, -0.001, 0.0015};
    imu_cfg.accel_noise_std_mps2 = 0.005; // 5 mm/s^2
    imu_cfg.dropouts = dropouts;
    imu_cfg.random_seed = 1001;
    sensors::ImuSensor imu(imu_cfg);

    // 2. GNSS: 1 Hz (dt = 1.0 s)
    sensors::GnssConfig gnss_cfg;
    gnss_cfg.sample_period_s = 1.0;
    gnss_cfg.position_bias_eci_m = math::Vector3{2.5, -1.5, 3.0};
    gnss_cfg.position_noise_std_m = 3.0; // 3 m
    gnss_cfg.velocity_bias_eci_mps = math::Vector3{0.02, -0.01, 0.015};
    gnss_cfg.velocity_noise_std_mps = 0.03; // 3 cm/s
    gnss_cfg.dropouts = dropouts;
    gnss_cfg.random_seed = 2002;
    sensors::GnssSensor gnss(gnss_cfg);

    // 3. Star Tracker: 10 Hz (dt = 0.1 s)
    sensors::StarTrackerConfig st_cfg;
    st_cfg.sample_period_s = 0.1;
    st_cfg.bias_axis = math::Vector3{0.0, 0.0, 1.0};
    st_cfg.bias_angle_rad = 0.0005; // 0.5 mrad boresight bias
    st_cfg.noise_std_rad = 0.0002;  // 0.2 mrad (41 arcsec)
    st_cfg.dropouts = dropouts;
    st_cfg.random_seed = 3003;
    sensors::StarTrackerSensor star_tracker(st_cfg);

    // 4. Range Sensor: 10 Hz (dt = 0.1 s)
    sensors::RangeSensorConfig range_cfg;
    range_cfg.sample_period_s = 0.1;
    range_cfg.bias_m = 0.5;
    range_cfg.noise_std_m = 0.1; // 10 cm
    range_cfg.dropouts = dropouts;
    range_cfg.random_seed = 4004;
    sensors::RangeSensor range_sensor(range_cfg);

    // Open Output CSV Files with full 16-digit double precision
    std::ofstream imu_file("data/m12_imu.csv");
    imu_file << std::setprecision(16);
    imu_file << "time_s,valid,truth_omega_x_rad_s,truth_omega_y_rad_s,truth_omega_z_rad_s,"
             << "meas_omega_x_rad_s,meas_omega_y_rad_s,meas_omega_z_rad_s,"
             << "truth_spec_force_x_mps2,truth_spec_force_y_mps2,truth_spec_force_z_mps2,"
             << "meas_spec_force_x_mps2,meas_spec_force_y_mps2,meas_spec_force_z_mps2\n";

    std::ofstream gnss_file("data/m12_gnss.csv");
    gnss_file << std::setprecision(16);
    gnss_file << "time_s,valid,truth_r_x_m,truth_r_y_m,truth_r_z_m,"
              << "truth_v_x_mps,truth_v_y_mps,truth_v_z_mps,"
              << "meas_r_x_m,meas_r_y_m,meas_r_z_m,"
              << "meas_v_x_mps,meas_v_y_mps,meas_v_z_mps\n";

    std::ofstream st_file("data/m12_star_tracker.csv");
    st_file << std::setprecision(16);
    st_file << "time_s,valid,truth_q_w,truth_q_x,truth_q_y,truth_q_z,"
            << "meas_q_w,meas_q_x,meas_q_y,meas_q_z,orientation_error_rad\n";

    std::ofstream range_file("data/m12_range.csv");
    range_file << std::setprecision(16);
    range_file << "time_s,valid,truth_range_m,meas_range_m,"
               << "target_r_x_m,target_r_y_m,target_r_z_m\n";

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Starting Simulation:\n";
    std::cout << "  Duration:          " << sim_duration_s << " s (1 orbit)\n";
    std::cout << "  Timestep:          " << dt_s << " s\n";
    std::cout << "  Sampling Rates:    IMU 100 Hz, Star Tracker 10 Hz, Range 10 Hz, GNSS 1 Hz\n";
    std::cout << "  Dropout Window:    [2000.0 s, 2100.0 s]\n\n";

    double current_t_s = 0.0;
    spacecraft::SpacecraftState current_state = initial_state;

    std::size_t count_imu = 0;
    std::size_t count_gnss = 0;
    std::size_t count_st = 0;
    std::size_t count_range = 0;

    while (current_t_s <= sim_duration_s + 1.0e-9) {
        // Evaluate non-gravitational acceleration in ECI (atmospheric drag)
        const math::Vector3 a_drag_eci = environment::drag_acceleration_eci(
            current_state.velocity(),
            current_state.position(),
            env_params.mass_kg,
            env_params.drag_coefficient_cd,
            env_params.drag_reference_area_m2,
            env_params.earth_reference_radius_m,
            env_params.earth_rotation_rate_rad_per_s,
            env_params.atmosphere_ref_alt_m,
            env_params.atmosphere_ref_density_kg_per_m3,
            env_params.atmosphere_scale_height_m
        );

        // Target Spacecraft Position at current time (in circular orbit)
        const double mean_motion = 2.0 * constants::pi / orbital_period_s;
        const double target_nu = mean_motion * current_t_s + target_lead_angle_rad;
        const math::Vector3 target_pos_eci{
            r_orbit_m * std::cos(target_nu),
            r_orbit_m * std::sin(target_nu),
            0.0
        };

        // 1. Sample IMU (100 Hz)
        if (auto imu_sample = imu.sample(current_t_s, current_state, a_drag_eci)) {
            const math::Vector3 truth_omega = current_state.angular_velocity_rad_per_s();
            const math::Vector3 truth_spec_force = current_state.orientation().normalized().conjugate().rotate_vector(a_drag_eci);

            imu_file << imu_sample->timestamp_s << ","
                     << (imu_sample->valid ? 1 : 0) << ","
                     << truth_omega.x() << "," << truth_omega.y() << "," << truth_omega.z() << ","
                     << imu_sample->angular_velocity_body_rad_s.x() << ","
                     << imu_sample->angular_velocity_body_rad_s.y() << ","
                     << imu_sample->angular_velocity_body_rad_s.z() << ","
                     << truth_spec_force.x() << "," << truth_spec_force.y() << "," << truth_spec_force.z() << ","
                     << imu_sample->specific_force_body_mps2.x() << ","
                     << imu_sample->specific_force_body_mps2.y() << ","
                     << imu_sample->specific_force_body_mps2.z() << "\n";
            ++count_imu;
        }

        // 2. Sample GNSS (1 Hz)
        if (auto gnss_sample = gnss.sample(current_t_s, current_state)) {
            gnss_file << gnss_sample->timestamp_s << ","
                      << (gnss_sample->valid ? 1 : 0) << ","
                      << current_state.position().x() << ","
                      << current_state.position().y() << ","
                      << current_state.position().z() << ","
                      << current_state.velocity().x() << ","
                      << current_state.velocity().y() << ","
                      << current_state.velocity().z() << ","
                      << gnss_sample->position_eci_m.x() << ","
                      << gnss_sample->position_eci_m.y() << ","
                      << gnss_sample->position_eci_m.z() << ","
                      << gnss_sample->velocity_eci_mps.x() << ","
                      << gnss_sample->velocity_eci_mps.y() << ","
                      << gnss_sample->velocity_eci_mps.z() << "\n";
            ++count_gnss;
        }

        // 3. Sample Star Tracker (10 Hz)
        if (auto st_sample = star_tracker.sample(current_t_s, current_state)) {
            const math::Quaternion truth_q = current_state.orientation().normalized();
            const double err_rad = st_sample->valid ?
                spacecraft::quaternion_orientation_error_rad(truth_q, st_sample->orientation_eci_from_body) : 0.0;

            st_file << st_sample->timestamp_s << ","
                    << (st_sample->valid ? 1 : 0) << ","
                    << truth_q.w() << "," << truth_q.x() << "," << truth_q.y() << "," << truth_q.z() << ","
                    << st_sample->orientation_eci_from_body.w() << ","
                    << st_sample->orientation_eci_from_body.x() << ","
                    << st_sample->orientation_eci_from_body.y() << ","
                    << st_sample->orientation_eci_from_body.z() << ","
                    << err_rad << "\n";
            ++count_st;
        }

        // 4. Sample Range Sensor (10 Hz)
        if (auto range_sample = range_sensor.sample(current_t_s, current_state, target_pos_eci)) {
            const double truth_range = (target_pos_eci - current_state.position()).norm();

            range_file << range_sample->timestamp_s << ","
                       << (range_sample->valid ? 1 : 0) << ","
                       << truth_range << ","
                       << range_sample->range_m << ","
                       << target_pos_eci.x() << ","
                       << target_pos_eci.y() << ","
                       << target_pos_eci.z() << "\n";
            ++count_range;
        }

        // Advance simulation truth state with environmental RK4
        current_state = environment::rk4_step_spacecraft_environmental(
            current_t_s,
            current_state,
            dt_s,
            sc_params,
            env_params,
            env_config
        );
        current_t_s += dt_s;
    }

    std::cout << "Simulation Complete. Export Summary:\n";
    std::cout << "  IMU Telemetry:          data/m12_imu.csv (" << count_imu << " samples)\n";
    std::cout << "  GNSS Telemetry:         data/m12_gnss.csv (" << count_gnss << " samples)\n";
    std::cout << "  Star Tracker Telemetry: data/m12_star_tracker.csv (" << count_st << " samples)\n";
    std::cout << "  Range Telemetry:        data/m12_range.csv (" << count_range << " samples)\n";
    std::cout << "=================================================================\n";

    return 0;
}
