#include "frames/frame_basis.hpp"
#include "frames/lvlh.hpp"
#include "math/constants.hpp"
#include "math/matrix3.hpp"
#include "math/vector3.hpp"
#include "numerics/fixed_step_propagation.hpp"
#include "orbit/cartesian_state.hpp"
#include "orbit/two_body_orbit.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using astradock::frames::compute_lvlh_basis;
using astradock::frames::dcm_lvlh_from_eci;
using astradock::frames::FrameBasis;
using astradock::frames::is_orthonormal;
using astradock::frames::transform_eci_to_lvlh;
using astradock::frames::transform_lvlh_to_eci;
using astradock::math::Matrix3;
using astradock::math::Vector3;
using astradock::numerics::IntegrationMethod;
using astradock::numerics::propagate_fixed_step;
using astradock::orbit::CartesianState;
using astradock::orbit::compute_circular_orbit_reference;
using astradock::orbit::two_body_state_derivative;

int main() {
    std::cout << "=================================================================\n";
    std::cout << " AstraDock — M06: Coordinate Frames & LVLH Reference Frames Demo \n";
    std::cout << "=================================================================\n\n";

    const double altitude_m = 500000.0;
    const double radius_m = astradock::constants::earth_reference_radius_m + altitude_m;
    const auto ref = compute_circular_orbit_reference(
        astradock::constants::earth_gravitational_parameter_m3_per_s2,
        radius_m);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Scenario Configuration:\n";
    std::cout << "  Altitude:            " << altitude_m / 1000.0 << " km\n";
    std::cout << "  Orbital Radius:      " << radius_m / 1000.0 << " km\n";
    std::cout << "  Circular Speed:      " << ref.speed_m_per_s << " m/s\n";
    std::cout << "  Orbital Period:      " << ref.period_s << " s (" << ref.period_s / 60.0 << " min)\n\n";

    const CartesianState initial_state{
        Vector3(radius_m, 0.0, 0.0),
        Vector3(0.0, ref.speed_m_per_s, 0.0)
    };

    const double dt_s = 10.0;
    const auto trajectory = propagate_fixed_step(
        0.0,
        ref.period_s,
        dt_s,
        initial_state,
        IntegrationMethod::classical_rk4,
        [](double t, const CartesianState& s) {
            return two_body_state_derivative(
                t, s, astradock::constants::earth_gravitational_parameter_m3_per_s2);
        });

    std::cout << "Propagated " << trajectory.size() << " samples across one full orbit.\n\n";

    // Setup CSV export
    const std::filesystem::path output_dir = "artifacts/data";
    std::filesystem::create_directories(output_dir);
    const std::filesystem::path csv_path = output_dir / "m06_frames.csv";
    std::ofstream csv(csv_path);

    if (!csv.is_open()) {
        std::cerr << "Failed to open output CSV: " << csv_path << "\n";
        return 1;
    }

    csv << "time_s,"
        << "position_eci_x_m,position_eci_y_m,position_eci_z_m,"
        << "velocity_eci_x_m_per_s,velocity_eci_y_m_per_s,velocity_eci_z_m_per_s,"
        << "lvlh_x_eci_x,lvlh_x_eci_y,lvlh_x_eci_z,"
        << "lvlh_y_eci_x,lvlh_y_eci_y,lvlh_y_eci_z,"
        << "lvlh_z_eci_x,lvlh_z_eci_y,lvlh_z_eci_z,"
        << "r_lvlh_x_m,r_lvlh_y_m,r_lvlh_z_m,"
        << "v_lvlh_x_m_per_s,v_lvlh_y_m_per_s,v_lvlh_z_m_per_s,"
        << "det_c,orthonormality_max_err,norm_diff_m,round_trip_err_m\n";

    double max_det_err = 0.0;
    double max_ortho_err = 0.0;
    double max_norm_err = 0.0;
    double max_round_trip_err = 0.0;

    std::cout << "Orbital Phase Samples (ECI vs LVLH Components):\n";
    std::cout << "----------------------------------------------------------------------------------------------------\n";
    std::cout << " Phase  | Time (s) |   ECI Position (km)   |   LVLH Position (km)  |   ECI Velocity (km/s) | LVLH Vel (km/s)\n";
    std::cout << "----------------------------------------------------------------------------------------------------\n";

    for (std::size_t i = 0; i < trajectory.size(); ++i) {
        const auto& sample = trajectory[i];
        const double t = sample.time_s;
        const Vector3 r_eci = sample.state.position;
        const Vector3 v_eci = sample.state.velocity;

        const FrameBasis basis = compute_lvlh_basis(r_eci, v_eci);
        const Matrix3 c_lvlh_eci = dcm_lvlh_from_eci(basis);

        const Vector3 r_lvlh = transform_eci_to_lvlh(r_eci, c_lvlh_eci);
        const Vector3 v_lvlh = transform_eci_to_lvlh(v_eci, c_lvlh_eci);

        const Vector3 r_recon = transform_lvlh_to_eci(r_lvlh, c_lvlh_eci);

        // Verification diagnostics
        const double det = c_lvlh_eci.determinant();
        const double det_err = std::abs(det - 1.0);
        max_det_err = std::max(max_det_err, det_err);

        const double ortho_err = std::max({
            std::abs(basis.x.norm() - 1.0),
            std::abs(basis.y.norm() - 1.0),
            std::abs(basis.z.norm() - 1.0),
            std::abs(basis.x.dot(basis.y)),
            std::abs(basis.y.dot(basis.z)),
            std::abs(basis.z.dot(basis.x))
        });
        max_ortho_err = std::max(max_ortho_err, ortho_err);

        const double norm_diff = std::abs(r_lvlh.norm() - r_eci.norm());
        max_norm_err = std::max(max_norm_err, norm_diff);

        const double round_trip_err = (r_recon - r_eci).norm();
        max_round_trip_err = std::max(max_round_trip_err, round_trip_err);

        csv << std::setprecision(6) << t << ","
            << std::setprecision(8)
            << r_eci.x() << "," << r_eci.y() << "," << r_eci.z() << ","
            << v_eci.x() << "," << v_eci.y() << "," << v_eci.z() << ","
            << basis.x.x() << "," << basis.x.y() << "," << basis.x.z() << ","
            << basis.y.x() << "," << basis.y.y() << "," << basis.y.z() << ","
            << basis.z.x() << "," << basis.z.y() << "," << basis.z.z() << ","
            << r_lvlh.x() << "," << r_lvlh.y() << "," << r_lvlh.z() << ","
            << v_lvlh.x() << "," << v_lvlh.y() << "," << v_lvlh.z() << ","
            << std::setprecision(14)
            << det << "," << ortho_err << "," << norm_diff << "," << round_trip_err << "\n";

        // Print key phase points (0, 90, 180, 270, 360 deg)
        const double frac = t / ref.period_s;
        if (i == 0 || (i == trajectory.size() / 4) || (i == trajectory.size() / 2)
            || (i == 3 * trajectory.size() / 4) || (i == trajectory.size() - 1)) {
            const double deg = frac * 360.0;
            std::cout << std::setw(5) << std::setprecision(0) << deg << "° | "
                      << std::setw(8) << std::setprecision(1) << t << " | "
                      << "[" << std::setw(7) << std::setprecision(1) << r_eci.x() / 1e3 << ", "
                             << std::setw(7) << r_eci.y() / 1e3 << ", "
                             << std::setw(4) << r_eci.z() / 1e3 << "] | "
                      << "[" << std::setw(7) << std::setprecision(1) << r_lvlh.x() / 1e3 << ", "
                             << std::setw(4) << r_lvlh.y() / 1e3 << ", "
                             << std::setw(4) << r_lvlh.z() / 1e3 << "] | "
                      << "[" << std::setw(6) << std::setprecision(2) << v_eci.x() / 1e3 << ", "
                             << std::setw(6) << v_eci.y() / 1e3 << ", "
                             << std::setw(4) << v_eci.z() / 1e3 << "] | "
                      << "[" << std::setw(4) << std::setprecision(2) << v_lvlh.x() / 1e3 << ", "
                             << std::setw(6) << v_lvlh.y() / 1e3 << ", "
                             << std::setw(4) << v_lvlh.z() / 1e3 << "]\n";
        }
    }

    std::cout << "----------------------------------------------------------------------------------------------------\n\n";
    std::cout << "Verification Summary Across Orbit:\n";
    std::cout << "  Max Determinant Error (|det(C) - 1|): " << std::scientific << max_det_err << "\n";
    std::cout << "  Max Orthonormality Error:             " << max_ortho_err << "\n";
    std::cout << "  Max Position Norm Diff (|r_L| - |r_E|):" << max_norm_err << " m\n";
    std::cout << "  Max Round-Trip Error (||C^T C r - r||):" << max_round_trip_err << " m\n\n";

    std::cout << "Exported frame data to: " << csv_path.string() << "\n";
    std::cout << "M06 C++ frame verification complete.\n";

    return 0;
}
