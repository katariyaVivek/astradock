// AstraDock M24 — High-fidelity environment demonstration tool.
//
// Compares legacy vs upgraded environment telemetry on one LEO orbit:
//   fixed vs ephemeris Moon position, J2-only vs J2/J3/J4 acceleration,
//   ECI/ECEF ground-track + geodetic coordinates, Julian-date time tags.

#include "environment/high_fidelity.hpp"
#include "environment/j2_gravity.hpp"
#include "environment/third_body_gravity.hpp"
#include "math/constants.hpp"
#include "math/vector3.hpp"
#include "numerics/fixed_step_propagation.hpp"
#include "numerics/integrators.hpp"
#include "orbit/cartesian_state.hpp"
#include "orbit/two_body_orbit.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

int main() {
    using namespace astradock;
    std::cout << "============================================================\n";
    std::cout << " AstraDock — M24 High-Fidelity Environment Demo             \n";
    std::cout << "============================================================\n";
    std::filesystem::create_directories("data");

    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const double R = constants::earth_reference_radius_m;
    const double r_orbit = R + 500.0e3;
    const double v_circ = std::sqrt(mu / r_orbit);
    orbit::CartesianState state{{r_orbit, 0.0, 0.0}, {0.0, v_circ, 0.0}};

    std::ofstream out("data/m24_environment.csv");
    out << "time_s,julian_date,lat_deg,lon_deg,alt_m,j2_perturb_mag,j234_total_mag,moon_fixed_err_km,moon_range_km\n";
    out << std::setprecision(10);
    auto two_body = [&](double, const orbit::CartesianState& s) {
        return orbit::two_body_state_derivative(0.0, s, mu);
    };
    const math::Vector3 moon_fixed{384400.0e3, 0.0, 0.0};
    const double dt = 60.0;
    double t = 0.0;
    for (int i = 0; i <= 100; ++i) {
        const math::Vector3 ecef = environment::eci_to_ecef(state.position, t);
        const auto geo = environment::ecef_to_geodetic(ecef);
        const math::Vector3 a_j2 = environment::j2_acceleration_eci(state.position);
        const math::Vector3 a_j234 = environment::zonal_acceleration_eci(
            state.position, mu, R, environment::k_j2_ref, environment::k_j3_ref,
            environment::k_j4_ref);
        const math::Vector3 moon = environment::moon_position_eci_m(t);
        out << t << "," << environment::sim_time_to_julian_date(t) << ","
            << geo.latitude_rad * 180.0 / constants::pi << ","
            << geo.longitude_rad * 180.0 / constants::pi << "," << geo.altitude_m << ","
            << a_j2.norm() << "," << a_j234.norm() << ","
            << (moon - moon_fixed).norm() / 1000.0 << "," << moon.norm() / 1000.0 << "\n";
        state = numerics::rk4_step(t, state, dt, two_body);
        t += dt;
    }
    std::cout << "  Exported 100-minute LEO telemetry to data/m24_environment.csv\n";
    std::cout << "============================================================\n";
    return 0;
}
