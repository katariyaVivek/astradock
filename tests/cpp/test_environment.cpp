#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "environment/atmospheric_drag.hpp"
#include "environment/environment_models.hpp"
#include "environment/gravity_gradient.hpp"
#include "environment/j2_gravity.hpp"
#include "environment/third_body_gravity.hpp"
#include "math/constants.hpp"
#include "math/vector3.hpp"
#include "orbit/classical_elements.hpp"
#include "orbit/two_body_orbit.hpp"
#include "spacecraft/six_dof_dynamics.hpp"

#include <cmath>

using namespace astradock;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("M10 disabled-environment regression preserves M10 bitwise behavior", "[environment][m10_regression]") {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const attitude::PrincipalInertia inertia{10.0, 20.0, 30.0};
    const spacecraft::SpacecraftParameters params{mu, inertia};

    const math::Vector3 r0{6878137.0, 0.0, 0.0};
    const math::Vector3 v0{0.0, 7612.608119, 0.0};
    const math::Quaternion q0{std::cos(constants::pi / 8.0), 0.0, 0.0, std::sin(constants::pi / 8.0)};
    const math::Vector3 w0{0.05, 0.08, 0.02};

    const spacecraft::SpacecraftState initial_state{
        orbit::CartesianState{r0, v0},
        attitude::RotationalState{q0, w0}
    };

    const environment::EnvironmentConfiguration disabled_config; // all false
    const environment::EnvironmentalParameters env_params;

    // 1. Derivative check
    const spacecraft::SpacecraftState m10_deriv = spacecraft::spacecraft_state_derivative(
        0.0, initial_state, params
    );
    const spacecraft::SpacecraftState m11_deriv = environment::spacecraft_environmental_derivative(
        0.0, initial_state, params, env_params, disabled_config
    );

    CHECK(m10_deriv.position().x() == m11_deriv.position().x());
    CHECK(m10_deriv.position().y() == m11_deriv.position().y());
    CHECK(m10_deriv.position().z() == m11_deriv.position().z());

    CHECK(m10_deriv.velocity().x() == m11_deriv.velocity().x());
    CHECK(m10_deriv.velocity().y() == m11_deriv.velocity().y());
    CHECK(m10_deriv.velocity().z() == m11_deriv.velocity().z());

    CHECK(m10_deriv.orientation().w() == m11_deriv.orientation().w());
    CHECK(m10_deriv.orientation().x() == m11_deriv.orientation().x());
    CHECK(m10_deriv.orientation().y() == m11_deriv.orientation().y());
    CHECK(m10_deriv.orientation().z() == m11_deriv.orientation().z());

    CHECK(m10_deriv.angular_velocity_rad_per_s().x() == m11_deriv.angular_velocity_rad_per_s().x());
    CHECK(m10_deriv.angular_velocity_rad_per_s().y() == m11_deriv.angular_velocity_rad_per_s().y());
    CHECK(m10_deriv.angular_velocity_rad_per_s().z() == m11_deriv.angular_velocity_rad_per_s().z());

    // 2. Trajectory propagation check (100 steps)
    const auto m10_traj = spacecraft::propagate_spacecraft_fixed_step(
        0.0, 100.0, 1.0, initial_state, params
    );
    const auto m11_traj = environment::propagate_spacecraft_environmental(
        initial_state, 100.0, 1.0, params, env_params, disabled_config
    );

    REQUIRE(m10_traj.size() == m11_traj.size());
    for (std::size_t i = 0; i < m10_traj.size(); ++i) {
        CHECK(m10_traj[i].state.position().x() == m11_traj[i].state.position().x());
        CHECK(m10_traj[i].state.position().y() == m11_traj[i].state.position().y());
        CHECK(m10_traj[i].state.position().z() == m11_traj[i].state.position().z());

        CHECK(m10_traj[i].state.velocity().x() == m11_traj[i].state.velocity().x());
        CHECK(m10_traj[i].state.velocity().y() == m11_traj[i].state.velocity().y());
        CHECK(m10_traj[i].state.velocity().z() == m11_traj[i].state.velocity().z());

        CHECK(m10_traj[i].state.orientation().w() == m11_traj[i].state.orientation().w());
        CHECK(m10_traj[i].state.orientation().x() == m11_traj[i].state.orientation().x());
        CHECK(m10_traj[i].state.orientation().y() == m11_traj[i].state.orientation().y());
        CHECK(m10_traj[i].state.orientation().z() == m11_traj[i].state.orientation().z());

        CHECK(m10_traj[i].state.angular_velocity_rad_per_s().x() == m11_traj[i].state.angular_velocity_rad_per_s().x());
        CHECK(m10_traj[i].state.angular_velocity_rad_per_s().y() == m11_traj[i].state.angular_velocity_rad_per_s().y());
        CHECK(m10_traj[i].state.angular_velocity_rad_per_s().z() == m11_traj[i].state.angular_velocity_rad_per_s().z());
    }
}

TEST_CASE("J2 gravitational perturbation acceleration satisfies analytical special cases and symmetries", "[environment][j2]") {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const double R_E = constants::earth_reference_radius_m;
    const double J2 = constants::earth_j2;
    const double r = 7000.0e3;

    // 1. Equatorial plane (z = 0, x = r, y = 0)
    // a_J2_eq = - 1.5 * mu * J2 * R_E^2 / r^4 in -X direction
    const math::Vector3 r_eq{r, 0.0, 0.0};
    const math::Vector3 a_eq = environment::j2_acceleration_eci(r_eq, mu, R_E, J2);
    const double expected_a_eq_mag = -1.5 * mu * J2 * (R_E * R_E) / (r * r * r * r);

    CHECK_THAT(a_eq.x(), WithinRel(expected_a_eq_mag, 1.0e-14));
    CHECK_THAT(a_eq.y(), WithinAbs(0.0, 1.0e-20));
    CHECK_THAT(a_eq.z(), WithinAbs(0.0, 1.0e-20));

    // 2. North Pole (x = y = 0, z = r)
    // a_J2_pole = - 1.5 * mu * J2 * R_E^2 / r^5 * z * (3 - 5) = + 3.0 * mu * J2 * R_E^2 / r^4 in +Z direction
    const math::Vector3 r_pole{0.0, 0.0, r};
    const math::Vector3 a_pole = environment::j2_acceleration_eci(r_pole, mu, R_E, J2);
    const double expected_a_pole_z = 3.0 * mu * J2 * (R_E * R_E) / (r * r * r * r);

    CHECK_THAT(a_pole.x(), WithinAbs(0.0, 1.0e-20));
    CHECK_THAT(a_pole.y(), WithinAbs(0.0, 1.0e-20));
    CHECK_THAT(a_pole.z(), WithinRel(expected_a_pole_z, 1.0e-14));

    // 3. South Pole (x = y = 0, z = -r) => a_z must point in -Z direction (equatorward)
    const math::Vector3 r_spole{0.0, 0.0, -r};
    const math::Vector3 a_spole = environment::j2_acceleration_eci(r_spole, mu, R_E, J2);
    CHECK_THAT(a_spole.z(), WithinRel(-expected_a_pole_z, 1.0e-14));

    // 4. Critical latitude where equatorial acceleration component vanishes (z^2/r^2 = 1/5)
    // sin(phi) = 1/sqrt(5) => z = r / sqrt(5), x = 2*r / sqrt(5)
    const double z_crit = r / std::sqrt(5.0);
    const double x_crit = 2.0 * r / std::sqrt(5.0);
    const math::Vector3 r_crit{x_crit, 0.0, z_crit};
    const math::Vector3 a_crit = environment::j2_acceleration_eci(r_crit, mu, R_E, J2);
    CHECK_THAT(a_crit.x(), WithinAbs(0.0, 1.0e-14));

    // 5. Distance scaling: r -> 2*r => a_J2 scales by 1 / 16 (1 / 2^4)
    const math::Vector3 r_2r{2.0 * r, 0.0, 0.0};
    const math::Vector3 a_2r = environment::j2_acceleration_eci(r_2r, mu, R_E, J2);
    CHECK_THAT(a_2r.x(), WithinRel(a_eq.x() / 16.0, 1.0e-14));

    // 6. Defensive checks
    CHECK_THROWS_AS(environment::j2_acceleration_eci(math::Vector3{0.0, 0.0, 0.0}), std::domain_error);
    CHECK_THROWS_AS(environment::j2_acceleration_eci(r_eq, -1.0), std::domain_error);
    CHECK_THROWS_AS(environment::j2_acceleration_eci(r_eq, mu, -1.0), std::domain_error);
}

TEST_CASE("J2 gravity causes secular RAAN regression and apsidal precession over inclined orbit", "[environment][j2][precession]") {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const attitude::PrincipalInertia inertia{10.0, 10.0, 10.0};
    const spacecraft::SpacecraftParameters params{mu, inertia};

    // Orbit: a = 10,000 km, e = 0.1, i = 45 deg, RAAN = 60 deg, omega = 30 deg, nu = 0 deg
    const orbit::ClassicalOrbitalElements el0{
        10000.0e3,
        0.1,
        45.0 * constants::pi / 180.0,
        60.0 * constants::pi / 180.0,
        30.0 * constants::pi / 180.0,
        0.0
    };
    const orbit::CartesianState cart0 = orbit::classical_elements_to_state(el0, mu);

    const spacecraft::SpacecraftState initial_state{
        cart0,
        attitude::RotationalState{}
    };

    const double period = 2.0 * constants::pi * std::sqrt(std::pow(el0.semi_major_axis_m, 3) / mu);
    const double sim_duration = 3.0 * period; // 3 orbital periods (~8.3 hours)
    const double dt = 5.0; // 5 second step

    environment::EnvironmentConfiguration config_j2;
    config_j2.enable_j2 = true;
    const environment::EnvironmentalParameters env_params;

    const auto traj = environment::propagate_spacecraft_environmental(
        initial_state, sim_duration, dt, params, env_params, config_j2
    );

    const auto final_cart = traj.back().state.translational;
    const auto final_el = orbit::state_to_classical_elements(final_cart, mu);

    CHECK_THAT(final_el.semi_major_axis_m, WithinRel(el0.semi_major_axis_m, 1.0e-3));
    CHECK_THAT(final_el.inclination_rad, WithinRel(el0.inclination_rad, 1.0e-3));

    // RAAN must strictly regress (decrease over time for prograde orbit)
    CHECK(final_el.raan_rad < el0.raan_rad);

    // Argument of periapsis must advance (increase over time for i = 45 deg)
    CHECK(final_el.argument_of_periapsis_rad > el0.argument_of_periapsis_rad);
}

TEST_CASE("Atmospheric drag relative velocity, exponential density, and force scaling", "[environment][drag]") {
    const double omega_E = constants::earth_rotation_rate_rad_per_s;
    const double R_E = constants::earth_reference_radius_m;

    // 1. Prograde equatorial state (x = r, y = 0, z = 0, vx = 0, vy = 7500, vz = 0)
    const math::Vector3 r_eq{R_E + 300.0e3, 0.0, 0.0};
    const math::Vector3 v_prograde{0.0, 7500.0, 0.0};

    const math::Vector3 v_rel_pro = environment::relative_atmospheric_velocity_eci(v_prograde, r_eq, omega_E);
    // Atmosphere speed at this radius: v_atm = omega_E * r
    const double v_atm = omega_E * r_eq.x();
    CHECK_THAT(v_rel_pro.x(), WithinAbs(0.0, 1.0e-12));
    CHECK_THAT(v_rel_pro.y(), WithinAbs(7500.0 - v_atm, 1.0e-12));
    CHECK_THAT(v_rel_pro.z(), WithinAbs(0.0, 1.0e-12));

    // 2. Exponential density scaling: Altitude h = 500 km matches reference density exactly
    const double rho_500 = environment::exponential_atmospheric_density(500.0e3, 500.0e3, 6.967e-13, 63.8e3);
    CHECK_THAT(rho_500, WithinRel(6.967e-13, 1.0e-14));

    // Higher altitude has lower density
    const double rho_600 = environment::exponential_atmospheric_density(600.0e3, 500.0e3, 6.967e-13, 63.8e3);
    CHECK(rho_600 < rho_500);

    // 3. Drag acceleration direction opposes v_rel
    const double mass = 500.0;
    const double Cd = 2.2;
    const double area = 2.0;

    const math::Vector3 a_drag = environment::drag_acceleration_eci(
        v_prograde, r_eq, mass, Cd, area, R_E, omega_E
    );

    // Drag must oppose v_rel
    CHECK(a_drag.y() < 0.0);
    CHECK_THAT(a_drag.x(), WithinAbs(0.0, 1.0e-18));
    CHECK_THAT(a_drag.z(), WithinAbs(0.0, 1.0e-18));

    // 4. Area and mass scaling: Double mass => Halves drag
    const math::Vector3 a_drag_2m = environment::drag_acceleration_eci(
        v_prograde, r_eq, 2.0 * mass, Cd, area, R_E, omega_E
    );
    CHECK_THAT(a_drag_2m.y(), WithinRel(a_drag.y() / 2.0, 1.0e-14));
}

TEST_CASE("Atmospheric drag causes orbital energy and semi-major axis secular decay", "[environment][drag][decay]") {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const attitude::PrincipalInertia inertia{10.0, 10.0, 10.0};
    const spacecraft::SpacecraftParameters params{mu, inertia};

    // Very low circular orbit (250 km) where drag is prominent
    const double alt = 250.0e3;
    const double r_orb = constants::earth_reference_radius_m + alt;
    const double v_circ = std::sqrt(mu / r_orb);

    const spacecraft::SpacecraftState initial_state{
        orbit::CartesianState{math::Vector3{r_orb, 0.0, 0.0}, math::Vector3{0.0, v_circ, 0.0}},
        attitude::RotationalState{}
    };

    const double period = 2.0 * constants::pi * std::sqrt(std::pow(r_orb, 3) / mu);
    const double sim_duration = 2.0 * period; // 2 orbits (~3 hours)

    environment::EnvironmentConfiguration config_drag;
    config_drag.enable_drag = true;

    environment::EnvironmentalParameters env_params;
    env_params.mass_kg = 200.0;
    env_params.drag_reference_area_m2 = 5.0; // Large area-to-mass to produce distinct decay
    env_params.atmosphere_ref_alt_m = 250.0e3;
    env_params.atmosphere_ref_density_kg_per_m3 = 6.0e-11; // 250 km typical density
    env_params.atmosphere_scale_height_m = 40.0e3;

    const auto traj = environment::propagate_spacecraft_environmental(
        initial_state, sim_duration, 1.0, params, env_params, config_drag
    );

    const double e0 = orbit::specific_orbital_energy_m2_per_s2(traj.front().state.translational, mu);
    const double e_final = orbit::specific_orbital_energy_m2_per_s2(traj.back().state.translational, mu);

    // Specific orbital energy MUST strictly decrease (become more negative) under drag
    CHECK(e_final < e0);

    // Final radius must be smaller than initial radius
    const double r_final = traj.back().state.position().norm();
    CHECK(r_final < r_orb);
}

TEST_CASE("Third-body tidal gravitational acceleration satisfies analytical limits and tidal dipole properties", "[environment][third_body]") {
    const double mu_moon = constants::moon_gravitational_parameter_m3_per_s2;
    const math::Vector3 r_moon{384400.0e3, 0.0, 0.0}; // Moon on +X axis
    const double r_sc = 7000.0e3; // LEO spacecraft

    // 1. Zero mass third body produces zero acceleration
    const math::Vector3 a_zero_mu = environment::third_body_acceleration_eci(
        math::Vector3{r_sc, 0.0, 0.0}, r_moon, 0.0
    );
    CHECK_THAT(a_zero_mu.norm(), WithinAbs(0.0, 1.0e-20));

    // 2. Spacecraft at Earth center (r = 0) experiences exactly zero tidal acceleration
    const math::Vector3 a_center = environment::third_body_acceleration_eci(
        math::Vector3{0.0, 0.0, 0.0}, r_moon, mu_moon
    );
    CHECK_THAT(a_center.norm(), WithinAbs(0.0, 1.0e-20));

    // 3. Sub-lunar point (spacecraft between Earth and Moon: +X)
    // Direct term > Indirect term => Net acceleration toward Moon (+X)
    const math::Vector3 a_sublunar = environment::third_body_acceleration_eci(
        math::Vector3{r_sc, 0.0, 0.0}, r_moon, mu_moon
    );
    CHECK(a_sublunar.x() > 0.0);
    CHECK_THAT(a_sublunar.y(), WithinAbs(0.0, 1.0e-20));
    CHECK_THAT(a_sublunar.z(), WithinAbs(0.0, 1.0e-20));

    // 4. Anti-lunar point (spacecraft on opposite side from Moon: -X)
    // Indirect acceleration of Earth toward Moon > Direct acceleration of spacecraft
    // Net acceleration points AWAY from Earth center (-X)
    const math::Vector3 a_antilunar = environment::third_body_acceleration_eci(
        math::Vector3{-r_sc, 0.0, 0.0}, r_moon, mu_moon
    );
    CHECK(a_antilunar.x() < 0.0);

    // 5. Exact analytical 1D tidal acceleration along line of centers
    const double expected_exact_x = mu_moon * (1.0 / ((r_moon.x() - r_sc) * (r_moon.x() - r_sc)) - 1.0 / (r_moon.x() * r_moon.x()));
    CHECK_THAT(a_sublunar.x(), WithinRel(expected_exact_x, 1.0e-14));

    // Tidal dipole analytical approximation: a_tidal ~ 2 * mu_3 * r / r_3^3 (accurate to ~2.7% due to quadrupole)
    const double expected_tidal_x = 2.0 * mu_moon * r_sc / std::pow(r_moon.x(), 3);
    CHECK_THAT(a_sublunar.x(), WithinRel(expected_tidal_x, 0.05));
}

TEST_CASE("Gravity-gradient torque satisfies spherical symmetry, principal axis alignment, and analytical misalignment", "[environment][gravity_gradient]") {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const double r = 7000.0e3;
    const math::Vector3 r_eci{r, 0.0, 0.0};
    const math::Quaternion q_identity{1.0, 0.0, 0.0, 0.0};

    // 1. Spherical inertia body (Ixx = Iyy = Izz) produces exactly zero torque
    const attitude::PrincipalInertia spherical_inertia{20.0, 20.0, 20.0};
    const math::Vector3 tau_spherical = environment::gravity_gradient_torque_body(
        r_eci, q_identity, spherical_inertia, mu
    );
    CHECK_THAT(tau_spherical.norm(), WithinAbs(0.0, 1.0e-20));

    // 2. Asymmetric inertia aligned with principal axis produces zero torque
    const attitude::PrincipalInertia asymmetric_inertia{10.0, 20.0, 30.0};
    const math::Vector3 tau_aligned_x = environment::gravity_gradient_torque_body(
        r_eci, q_identity, asymmetric_inertia, mu
    );
    CHECK_THAT(tau_aligned_x.norm(), WithinAbs(0.0, 1.0e-20));

    // 3. 45 deg pitch offset about +Y (q = [cos(22.5°), 0, sin(22.5°), 0])
    // Rotating r_hat_I = [1, 0, 0]^T by q* (i.e. -45° about +Y):
    // r_hat_B = [cos(45°), 0, sin(45°)] = [1/sqrt(2), 0, 1/sqrt(2)]
    // tau_y = 3 * mu / r^3 * (Ixx - Izz) * ux * uz
    //       = 3 * mu / r^3 * (10 - 30) * (1/sqrt(2)) * (1/sqrt(2))
    //       = 3 * mu / r^3 * (-20) * (0.5) = - 30 * mu / r^3 (restoring pitch torque)
    const math::Quaternion q_pitch45{std::cos(constants::pi / 8.0), 0.0, std::sin(constants::pi / 8.0), 0.0};
    const math::Vector3 tau_pitch = environment::gravity_gradient_torque_body(
        r_eci, q_pitch45, asymmetric_inertia, mu
    );

    const double expected_tau_y = 3.0 * mu / (r * r * r) * (10.0 - 30.0) * (1.0 / std::sqrt(2.0)) * (1.0 / std::sqrt(2.0));

    CHECK_THAT(tau_pitch.x(), WithinAbs(0.0, 1.0e-20));
    CHECK_THAT(tau_pitch.y(), WithinRel(expected_tau_y, 1.0e-14));
    CHECK_THAT(tau_pitch.z(), WithinAbs(0.0, 1.0e-20));

    // Defensive checks
    CHECK_THROWS_AS(environment::gravity_gradient_torque_body(math::Vector3{0.0, 0.0, 0.0}, q_identity, asymmetric_inertia), std::domain_error);
    CHECK_THROWS_AS(environment::gravity_gradient_torque_body(r_eci, math::Quaternion{0.0, 0.0, 0.0, 0.0}, asymmetric_inertia), std::domain_error);
}

TEST_CASE("Gravity-gradient torque induces attitude libration in circular orbit", "[environment][gravity_gradient][libration]") {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const double r = 7000.0e3;
    const double v_circ = std::sqrt(mu / r);

    // Gravity-gradient stable dumbbell satellite: Ixx < Iyy = Izz
    const attitude::PrincipalInertia inertia{10.0, 50.0, 50.0};
    const spacecraft::SpacecraftParameters params{mu, inertia};

    // Small 5 degree initial pitch offset
    const double pitch_angle = 5.0 * constants::pi / 180.0;
    const math::Quaternion q0{std::cos(pitch_angle / 2.0), 0.0, std::sin(pitch_angle / 2.0), 0.0};
    const math::Vector3 w0{0.0, 0.0, 0.0}; // Initially at rest

    const spacecraft::SpacecraftState initial_state{
        orbit::CartesianState{math::Vector3{r, 0.0, 0.0}, math::Vector3{0.0, v_circ, 0.0}},
        attitude::RotationalState{q0, w0}
    };

    environment::EnvironmentConfiguration config_gg;
    config_gg.enable_gravity_gradient = true;
    const environment::EnvironmentalParameters env_params;

    const double sim_duration = 500.0; // 500 seconds
    const auto traj = environment::propagate_spacecraft_environmental(
        initial_state, sim_duration, 0.5, params, env_params, config_gg
    );

    // Spacecraft should develop periodic pitch angular velocity oscillations (librations)
    double max_wy = 0.0;
    for (const auto& sample : traj) {
        max_wy = std::max(max_wy, std::abs(sample.state.angular_velocity_rad_per_s().y()));
    }

    CHECK(max_wy > 1.0e-5); // Clear restorative oscillation induced by gravity-gradient torque
}

TEST_CASE("Environmental configuration and deterministic repeatability", "[environment][determinism]") {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const attitude::PrincipalInertia inertia{10.0, 20.0, 30.0};
    const spacecraft::SpacecraftParameters params{mu, inertia};

    const spacecraft::SpacecraftState state0{
        orbit::CartesianState{math::Vector3{6878137.0, 0.0, 0.0}, math::Vector3{0.0, 7612.608, 0.0}},
        attitude::RotationalState{math::Quaternion{1.0, 0.0, 0.0, 0.0}, math::Vector3{0.01, 0.02, 0.03}}
    };

    environment::EnvironmentConfiguration full_config{true, true, true, true};
    environment::EnvironmentalParameters env_params;

    const auto run1 = environment::propagate_spacecraft_environmental(
        state0, 200.0, 1.0, params, env_params, full_config
    );
    const auto run2 = environment::propagate_spacecraft_environmental(
        state0, 200.0, 1.0, params, env_params, full_config
    );

    REQUIRE(run1.size() == run2.size());
    for (std::size_t i = 0; i < run1.size(); ++i) {
        CHECK(run1[i].time_s == run2[i].time_s);
        CHECK(run1[i].state.position().x() == run2[i].state.position().x());
        CHECK(run1[i].state.velocity().y() == run2[i].state.velocity().y());
        CHECK(run1[i].state.orientation().w() == run2[i].state.orientation().w());
        CHECK(run1[i].state.angular_velocity_rad_per_s().z() == run2[i].state.angular_velocity_rad_per_s().z());
    }
}
