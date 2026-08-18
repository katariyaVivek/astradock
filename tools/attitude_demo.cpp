#include "attitude/attitude_state.hpp"
#include "math/constants.hpp"
#include "math/euler_angles.hpp"
#include "math/matrix3.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using astradock::attitude::AttitudeState;
using astradock::constants::pi;
using astradock::math::EulerAngles;
using astradock::math::euler_to_quaternion;
using astradock::math::Matrix3;
using astradock::math::Quaternion;
using astradock::math::quaternion_to_euler;
using astradock::math::Vector3;

int main() {
    std::cout << "=================================================================\n";
    std::cout << " AstraDock — M08 Attitude Representation & Quaternion Demo \n";
    std::cout << "=================================================================\n\n";

    // 1. Initial Spacecraft Attitude State
    const Vector3 rot_axis(1.0, 1.0, 1.0);
    const double rot_angle_rad = 45.0 * pi / 180.0;
    const Quaternion q_spacecraft = Quaternion::from_axis_angle(rot_axis, rot_angle_rad);
    const AttitudeState attitude{q_spacecraft};

    std::cout << "Spacecraft Attitude State (45 deg about [1, 1, 1]):\n";
    std::cout << "  Quaternion (w, x, y, z):     [" 
              << std::fixed << std::setprecision(6)
              << attitude.orientation.w() << ", "
              << attitude.orientation.x() << ", "
              << attitude.orientation.y() << ", "
              << attitude.orientation.z() << "]\n";
    std::cout << "  Quaternion Norm:             " << attitude.orientation.norm() << "\n";

    const Matrix3 dcm = attitude.orientation.to_rotation_matrix();
    std::cout << "  Direction Cosine Matrix:\n";
    std::cout << "    [" << std::setw(9) << dcm(0, 0) << ", " << std::setw(9) << dcm(0, 1) << ", " << std::setw(9) << dcm(0, 2) << "]\n";
    std::cout << "    [" << std::setw(9) << dcm(1, 0) << ", " << std::setw(9) << dcm(1, 1) << ", " << std::setw(9) << dcm(1, 2) << "]\n";
    std::cout << "    [" << std::setw(9) << dcm(2, 0) << ", " << std::setw(9) << dcm(2, 1) << ", " << std::setw(9) << dcm(2, 2) << "]\n";
    std::cout << "  DCM Determinant:             " << dcm.determinant() << "\n";

    const EulerAngles euler = quaternion_to_euler(attitude.orientation);
    std::cout << "  Equivalent ZYX Euler Angles: Yaw=" 
              << euler.yaw_rad * 180.0 / pi << " deg, Pitch="
              << euler.pitch_rad * 180.0 / pi << " deg, Roll="
              << euler.roll_rad * 180.0 / pi << " deg\n\n";

    // 2. Export discrete rotation sequence telemetry for analysis & visualization
    std::ofstream csv_file("artifacts/data/m08_attitude_telemetry.csv");
    if (!csv_file.is_open()) {
        std::cerr << "Error: Could not open artifacts/data/m08_attitude_telemetry.csv for writing.\n";
        return 1;
    }

    csv_file << "step,angle_deg,q_w,q_x,q_y,q_z,dcm_00,dcm_01,dcm_02,dcm_10,dcm_11,dcm_12,dcm_20,dcm_21,dcm_22,body_x_x,body_x_y,body_x_z,body_y_x,body_y_y,body_y_z,body_z_x,body_z_y,body_z_z\n";

    const int num_steps = 360;
    const Vector3 spin_axis(0.0, 0.0, 1.0); // Z spin
    for (int step = 0; step <= num_steps; ++step) {
        const double angle_deg = static_cast<double>(step);
        const double angle_rad = angle_deg * pi / 180.0;
        const Quaternion q = Quaternion::from_axis_angle(spin_axis, angle_rad);
        const Matrix3 m = q.to_rotation_matrix();

        const Vector3 body_x = q.rotate_vector(Vector3(1.0, 0.0, 0.0));
        const Vector3 body_y = q.rotate_vector(Vector3(0.0, 1.0, 0.0));
        const Vector3 body_z = q.rotate_vector(Vector3(0.0, 0.0, 1.0));

        csv_file << step << ","
                 << angle_deg << ","
                 << q.w() << "," << q.x() << "," << q.y() << "," << q.z() << ","
                 << m(0, 0) << "," << m(0, 1) << "," << m(0, 2) << ","
                 << m(1, 0) << "," << m(1, 1) << "," << m(1, 2) << ","
                 << m(2, 0) << "," << m(2, 1) << "," << m(2, 2) << ","
                 << body_x.x() << "," << body_x.y() << "," << body_x.z() << ","
                 << body_y.x() << "," << body_y.y() << "," << body_y.z() << ","
                 << body_z.x() << "," << body_z.y() << "," << body_z.z() << "\n";
    }
    csv_file.close();

    std::cout << "Saved attitude telemetry to: artifacts/data/m08_attitude_telemetry.csv\n\n";
    std::cout << "=================================================================\n";
    std::cout << " Demo execution complete. \n";
    std::cout << "=================================================================\n";

    return 0;
}
