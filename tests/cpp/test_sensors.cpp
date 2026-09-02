#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "math/constants.hpp"
#include "math/quaternion.hpp"
#include "math/vector3.hpp"
#include "sensors/gnss.hpp"
#include "sensors/imu.hpp"
#include "sensors/range_sensor.hpp"
#include "sensors/sensor_common.hpp"
#include "sensors/star_tracker.hpp"
#include "spacecraft/six_dof_dynamics.hpp"

#include <cmath>
#include <cstdint>
#include <numeric>
#include <vector>

using namespace astradock;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("DeterministicRng produces reproducible sequences and valid Gaussian statistics", "[sensors][rng]") {
    // 1. Identical seeds produce bitwise identical sequences
    sensors::DeterministicRng rng1(12345);
    sensors::DeterministicRng rng2(12345);

    for (int i = 0; i < 50; ++i) {
        CHECK(rng1.gaussian(0.0, 1.0) == rng2.gaussian(0.0, 1.0));
    }

    // 2. Different seeds produce distinct sequences
    sensors::DeterministicRng rng3(54321);
    CHECK(rng1.gaussian(0.0, 1.0) != rng3.gaussian(0.0, 1.0));

    // 3. Statistical mean and standard deviation over large sample size
    sensors::DeterministicRng rng_stat(42);
    const std::size_t N = 10000;
    const double target_mean = 5.0;
    const double target_std = 2.0;

    double sum = 0.0;
    double sum_sq = 0.0;

    for (std::size_t i = 0; i < N; ++i) {
        const double x = rng_stat.gaussian(target_mean, target_std);
        sum += x;
        sum_sq += x * x;
    }

    const double sample_mean = sum / static_cast<double>(N);
    const double sample_var = (sum_sq / static_cast<double>(N)) - (sample_mean * sample_mean);
    const double sample_std = std::sqrt(std::max(0.0, sample_var));

    // For N = 10,000, 3*sigma_mean = 3 * (2.0 / sqrt(10000)) = 0.06
    CHECK_THAT(sample_mean, WithinAbs(target_mean, 0.08));
    CHECK_THAT(sample_std, WithinAbs(target_std, 0.08));

    // 4. Uniform unit vector on S^2 always has norm 1.0
    for (int i = 0; i < 20; ++i) {
        const math::Vector3 u = rng_stat.uniform_unit_vector3();
        CHECK_THAT(u.norm(), WithinAbs(1.0, 1.0e-14));
    }

    // 5. Zero standard deviation returns exact mean
    CHECK(rng_stat.gaussian(3.14159, 0.0) == 3.14159);
    CHECK(rng_stat.gaussian_vector3(0.0).norm() == 0.0);
}

TEST_CASE("SensorSchedule and DropoutWindow manage timing and availability deterministically", "[sensors][timing]") {
    sensors::SensorSchedule sched(0.1, 0.0); // 10 Hz schedule

    CHECK(sched.sample_period_s() == 0.1);
    CHECK(sched.next_sample_time_s() == 0.0);
    CHECK(sched.should_sample(0.0));
    CHECK_FALSE(sched.should_sample(-0.01));

    // Advance to t = 0.0
    sched.advance(0.0);
    CHECK_THAT(sched.next_sample_time_s(), WithinAbs(0.1, 1.0e-12));
    CHECK_FALSE(sched.should_sample(0.05));
    CHECK(sched.should_sample(0.10));

    // Dropout windows
    std::vector<sensors::DropoutWindow> dropouts = {
        {10.0, 20.0},
        {50.0, 60.0}
    };

    CHECK_FALSE(sensors::is_in_dropout_window(9.99, dropouts));
    CHECK(sensors::is_in_dropout_window(10.0, dropouts));
    CHECK(sensors::is_in_dropout_window(15.0, dropouts));
    CHECK(sensors::is_in_dropout_window(20.0, dropouts));
    CHECK_FALSE(sensors::is_in_dropout_window(20.01, dropouts));
    CHECK(sensors::is_in_dropout_window(55.0, dropouts));
}

TEST_CASE("IMU Gyroscope measures body angular velocity with bias, noise, and frame fidelity", "[sensors][imu][gyro]") {
    const math::Vector3 true_omega{0.05, -0.02, 0.10}; // rad/s in Body frame
    const spacecraft::SpacecraftState state{
        orbit::CartesianState{math::Vector3{7000.0e3, 0.0, 0.0}, math::Vector3{0.0, 7500.0, 0.0}},
        attitude::RotationalState{math::Quaternion::identity(), true_omega}
    };

    // 1. Zero noise, zero bias -> exact body rate passthrough
    sensors::ImuConfig cfg_ideal;
    cfg_ideal.gyro_bias_rad_s = math::Vector3{0.0, 0.0, 0.0};
    cfg_ideal.gyro_noise_std_rad_s = 0.0;
    sensors::ImuSensor imu_ideal(cfg_ideal);

    const auto meas_ideal = imu_ideal.measure(0.0, state);
    CHECK(meas_ideal.valid);
    CHECK_THAT(meas_ideal.angular_velocity_body_rad_s.x(), WithinAbs(true_omega.x(), 1.0e-15));
    CHECK_THAT(meas_ideal.angular_velocity_body_rad_s.y(), WithinAbs(true_omega.y(), 1.0e-15));
    CHECK_THAT(meas_ideal.angular_velocity_body_rad_s.z(), WithinAbs(true_omega.z(), 1.0e-15));

    // 2. Constant bias addition
    sensors::ImuConfig cfg_bias;
    const math::Vector3 bias{0.001, -0.002, 0.003};
    cfg_bias.gyro_bias_rad_s = bias;
    cfg_bias.gyro_noise_std_rad_s = 0.0;
    sensors::ImuSensor imu_bias(cfg_bias);

    const auto meas_bias = imu_bias.measure(0.0, state);
    CHECK_THAT(meas_bias.angular_velocity_body_rad_s.x(), WithinAbs(true_omega.x() + bias.x(), 1.0e-15));
    CHECK_THAT(meas_bias.angular_velocity_body_rad_s.y(), WithinAbs(true_omega.y() + bias.y(), 1.0e-15));
    CHECK_THAT(meas_bias.angular_velocity_body_rad_s.z(), WithinAbs(true_omega.z() + bias.z(), 1.0e-15));

    // 3. Noise statistics
    sensors::ImuConfig cfg_noise;
    cfg_noise.gyro_bias_rad_s = bias;
    cfg_noise.gyro_noise_std_rad_s = 0.005; // 5 mrad/s noise
    cfg_noise.random_seed = 101;
    sensors::ImuSensor imu_noise(cfg_noise);

    const std::size_t N = 5000;
    double sum_err_x = 0.0;
    double sum_sq_err_x = 0.0;

    for (std::size_t i = 0; i < N; ++i) {
        const auto m = imu_noise.measure(static_cast<double>(i) * 0.01, state);
        const double err_x = m.angular_velocity_body_rad_s.x() - true_omega.x();
        sum_err_x += err_x;
        sum_sq_err_x += err_x * err_x;
    }

    const double mean_err_x = sum_err_x / static_cast<double>(N);
    const double var_err_x = (sum_sq_err_x / static_cast<double>(N)) - (mean_err_x * mean_err_x);
    const double std_err_x = std::sqrt(std::max(0.0, var_err_x));

    CHECK_THAT(mean_err_x, WithinAbs(bias.x(), 0.0005));
    CHECK_THAT(std_err_x, WithinAbs(cfg_noise.gyro_noise_std_rad_s, 0.0005));
}

TEST_CASE("IMU Accelerometer measures zero specific force in free fall and detects non-gravitational acceleration", "[sensors][imu][accel]") {
    // 1. Spacecraft in pure free-fall orbital motion (no non-gravitational forces)
    // Orientation: 45 degree roll about +X
    const math::Quaternion q_att = math::Quaternion::from_axis_angle(math::Vector3{1.0, 0.0, 0.0}, constants::pi / 4.0);
    const spacecraft::SpacecraftState state{
        orbit::CartesianState{math::Vector3{6878137.0, 0.0, 0.0}, math::Vector3{0.0, 7612.6, 0.0}},
        attitude::RotationalState{q_att, math::Vector3{}}
    };

    sensors::ImuConfig cfg_freefall;
    cfg_freefall.accel_bias_mps2 = math::Vector3{0.0, 0.0, 0.0};
    cfg_freefall.accel_noise_std_mps2 = 0.0;
    sensors::ImuSensor imu(cfg_freefall);

    // Free fall: specific force MUST be identically zero
    const auto meas_freefall = imu.measure(0.0, state, math::Vector3{0.0, 0.0, 0.0});
    CHECK(meas_freefall.valid);
    CHECK_THAT(meas_freefall.specific_force_body_mps2.norm(), WithinAbs(0.0, 1.0e-15));

    // 2. Non-gravitational acceleration in ECI (e.g. atmospheric drag opposing velocity in -Y)
    const math::Vector3 a_drag_eci{0.0, -0.05, 0.0}; // -0.05 m/s^2 along ECI Y
    const auto meas_drag = imu.measure(0.0, state, a_drag_eci);

    // Specific force in body frame: f_B = q* ⊗ f_I ⊗ q
    // With 45 deg roll about X, ECI Y rotates into Body (Y cos 45° - Z sin 45°, etc.)
    const math::Vector3 expected_f_body = q_att.conjugate().rotate_vector(a_drag_eci);
    CHECK_THAT(meas_drag.specific_force_body_mps2.x(), WithinAbs(expected_f_body.x(), 1.0e-15));
    CHECK_THAT(meas_drag.specific_force_body_mps2.y(), WithinAbs(expected_f_body.y(), 1.0e-15));
    CHECK_THAT(meas_drag.specific_force_body_mps2.z(), WithinAbs(expected_f_body.z(), 1.0e-15));
    CHECK_THAT(meas_drag.specific_force_body_mps2.norm(), WithinAbs(0.05, 1.0e-15));

    // 3. Dropout behavior
    sensors::ImuConfig cfg_dropout;
    cfg_dropout.dropouts = {{5.0, 10.0}};
    sensors::ImuSensor imu_drop(cfg_dropout);

    CHECK(imu_drop.measure(4.9, state).valid);
    const auto dropout_sample = imu_drop.measure(7.5, state);
    CHECK_FALSE(dropout_sample.valid);
    CHECK(dropout_sample.specific_force_body_mps2.norm() == 0.0);
    CHECK(dropout_sample.angular_velocity_body_rad_s.norm() == 0.0);
    CHECK(imu_drop.measure(10.1, state).valid);
}

TEST_CASE("GNSS Sensor measures absolute ECI position and velocity with bias, noise, and sampling schedule", "[sensors][gnss]") {
    const math::Vector3 true_r{6878137.0, 1000.0, -500.0};
    const math::Vector3 true_v{-10.0, 7612.6, 15.0};
    const spacecraft::SpacecraftState state{
        orbit::CartesianState{true_r, true_v},
        attitude::RotationalState{}
    };

    // 1. Ideal zero noise/bias passthrough
    sensors::GnssConfig cfg_ideal;
    cfg_ideal.sample_period_s = 1.0;
    sensors::GnssSensor gnss_ideal(cfg_ideal);

    const auto m_ideal = gnss_ideal.measure(0.0, state);
    CHECK(m_ideal.valid);
    CHECK_THAT(m_ideal.position_eci_m.x(), WithinAbs(true_r.x(), 1.0e-15));
    CHECK_THAT(m_ideal.position_eci_m.y(), WithinAbs(true_r.y(), 1.0e-15));
    CHECK_THAT(m_ideal.position_eci_m.z(), WithinAbs(true_r.z(), 1.0e-15));
    CHECK_THAT(m_ideal.velocity_eci_mps.x(), WithinAbs(true_v.x(), 1.0e-15));
    CHECK_THAT(m_ideal.velocity_eci_mps.y(), WithinAbs(true_v.y(), 1.0e-15));
    CHECK_THAT(m_ideal.velocity_eci_mps.z(), WithinAbs(true_v.z(), 1.0e-15));

    // 2. Bias and noise statistics
    sensors::GnssConfig cfg_noisy;
    cfg_noisy.position_bias_eci_m = math::Vector3{2.0, -1.0, 0.5};
    cfg_noisy.position_noise_std_m = 5.0; // 5 m 1-sigma position
    cfg_noisy.velocity_bias_eci_mps = math::Vector3{0.01, -0.02, 0.0};
    cfg_noisy.velocity_noise_std_mps = 0.05; // 5 cm/s 1-sigma velocity
    cfg_noisy.random_seed = 202;
    sensors::GnssSensor gnss_noisy(cfg_noisy);

    const std::size_t N = 5000;
    double sum_pos_err = 0.0;
    double sum_sq_pos_err = 0.0;

    for (std::size_t i = 0; i < N; ++i) {
        const auto m = gnss_noisy.measure(static_cast<double>(i), state);
        const double err = m.position_eci_m.x() - true_r.x();
        sum_pos_err += err;
        sum_sq_pos_err += err * err;
    }

    const double mean_pos_err = sum_pos_err / static_cast<double>(N);
    const double var_pos_err = (sum_sq_pos_err / static_cast<double>(N)) - (mean_pos_err * mean_pos_err);
    const double std_pos_err = std::sqrt(std::max(0.0, var_pos_err));

    CHECK_THAT(mean_pos_err, WithinAbs(cfg_noisy.position_bias_eci_m.x(), 0.3));
    CHECK_THAT(std_pos_err, WithinAbs(cfg_noisy.position_noise_std_m, 0.3));

    // 3. 1 Hz sampling schedule: samples only once per integer second
    gnss_ideal.reset();
    CHECK(gnss_ideal.sample(0.0, state).has_value());
    CHECK_FALSE(gnss_ideal.sample(0.2, state).has_value());
    CHECK_FALSE(gnss_ideal.sample(0.8, state).has_value());
    CHECK(gnss_ideal.sample(1.0, state).has_value());
}

TEST_CASE("Star Tracker measures attitude quaternion via valid physical SO(3) rotations", "[sensors][star_tracker]") {
    // 45 deg pitch about +Y
    const math::Quaternion true_q = math::Quaternion::from_axis_angle(math::Vector3{0.0, 1.0, 0.0}, constants::pi / 4.0);
    const spacecraft::SpacecraftState state{
        orbit::CartesianState{},
        attitude::RotationalState{true_q, math::Vector3{}}
    };

    // 1. Zero noise, zero bias -> exact truth quaternion
    sensors::StarTrackerConfig cfg_ideal;
    cfg_ideal.sample_period_s = 0.1;
    sensors::StarTrackerSensor st_ideal(cfg_ideal);

    const auto m_ideal = st_ideal.measure(0.0, state);
    CHECK(m_ideal.valid);
    CHECK_THAT(m_ideal.orientation_eci_from_body.norm(), WithinAbs(1.0, 1.0e-14));

    const double orientation_err = spacecraft::quaternion_orientation_error_rad(
        true_q, m_ideal.orientation_eci_from_body
    );
    CHECK_THAT(orientation_err, WithinAbs(0.0, 1.0e-14));

    // 2. Known small boresight bias rotation
    sensors::StarTrackerConfig cfg_bias;
    const double bias_angle = 0.01; // 10 mrad
    cfg_bias.bias_axis = math::Vector3{0.0, 0.0, 1.0};
    cfg_bias.bias_angle_rad = bias_angle;
    sensors::StarTrackerSensor st_bias(cfg_bias);

    const auto m_bias = st_bias.measure(0.0, state);
    CHECK_THAT(m_bias.orientation_eci_from_body.norm(), WithinAbs(1.0, 1.0e-14));
    const double bias_err = spacecraft::quaternion_orientation_error_rad(
        true_q, m_bias.orientation_eci_from_body
    );
    CHECK_THAT(bias_err, WithinAbs(bias_angle, 1.0e-12));

    // 3. Noise perturbation generates valid unit quaternions
    sensors::StarTrackerConfig cfg_noise;
    cfg_noise.noise_std_rad = 0.001; // 1 mrad 1-sigma attitude error
    cfg_noise.random_seed = 303;
    sensors::StarTrackerSensor st_noise(cfg_noise);

    for (int i = 0; i < 100; ++i) {
        const auto m = st_noise.measure(static_cast<double>(i) * 0.1, state);
        CHECK(m.valid);
        CHECK_THAT(m.orientation_eci_from_body.norm(), WithinAbs(1.0, 1.0e-14));
    }
}

TEST_CASE("Range Sensor measures scalar distance with translation invariance, bias, and noise", "[sensors][range]") {
    const math::Vector3 pos_sc{7000.0e3, 0.0, 0.0};
    const math::Vector3 pos_target{7000.0e3, 500.0, 0.0}; // 500 m separation along +Y
    const double true_dist = 500.0;

    const spacecraft::SpacecraftState state{
        orbit::CartesianState{pos_sc, math::Vector3{}},
        attitude::RotationalState{}
    };

    // 1. Ideal measurement
    sensors::RangeSensorConfig cfg_ideal;
    sensors::RangeSensor range_ideal(cfg_ideal);

    const auto m_ideal = range_ideal.measure(0.0, state, pos_target);
    CHECK(m_ideal.valid);
    CHECK_THAT(m_ideal.range_m, WithinAbs(true_dist, 1.0e-12));

    // 2. Translation Invariance: Shift both spacecraft and target by arbitrary vector c
    const math::Vector3 shift{1.234e7, -4.567e6, 8.910e5};
    const spacecraft::SpacecraftState state_shifted{
        orbit::CartesianState{pos_sc + shift, math::Vector3{}},
        attitude::RotationalState{}
    };
    const auto m_shifted = range_ideal.measure(0.0, state_shifted, pos_target + shift);
    CHECK_THAT(m_shifted.range_m, WithinAbs(true_dist, 1.0e-12));

    // 3. Bias and noise
    sensors::RangeSensorConfig cfg_noisy;
    cfg_noisy.bias_m = 1.5;
    cfg_noisy.noise_std_m = 0.2;
    cfg_noisy.random_seed = 404;
    sensors::RangeSensor range_noisy(cfg_noisy);

    const auto m_noisy = range_noisy.measure(0.0, state, pos_target);
    CHECK(m_noisy.valid);
    // Measured range should be close to 500 + 1.5 = 501.5 m
    CHECK_THAT(m_noisy.range_m, WithinAbs(501.5, 1.0));
}

TEST_CASE("Sensor simulation obeys truth non-interference contract", "[sensors][truth_regression]") {
    const double mu = constants::earth_gravitational_parameter_m3_per_s2;
    const attitude::PrincipalInertia inertia{10.0, 20.0, 30.0};
    const spacecraft::SpacecraftParameters params{mu, inertia};

    const spacecraft::SpacecraftState initial_state{
        orbit::CartesianState{math::Vector3{6878137.0, 0.0, 0.0}, math::Vector3{0.0, 7612.608, 0.0}},
        attitude::RotationalState{math::Quaternion{1.0, 0.0, 0.0, 0.0}, math::Vector3{0.01, 0.02, 0.015}}
    };

    const double dt = 0.1;
    const int total_steps = 1000;

    // 1. Run truth simulation standalone
    std::vector<spacecraft::SpacecraftState> truth_traj_standalone;
    {
        double current_t = 0.0;
        spacecraft::SpacecraftState state = initial_state;
        truth_traj_standalone.push_back(state);
        for (int step = 0; step < total_steps; ++step) {
            state = spacecraft::rk4_step_spacecraft(current_t, state, dt, params);
            current_t += dt;
            truth_traj_standalone.push_back(state);
        }
    }

    // 2. Initialize full sensor suite
    sensors::ImuSensor imu(sensors::ImuConfig{});
    sensors::GnssSensor gnss(sensors::GnssConfig{});
    sensors::StarTrackerSensor star_tracker(sensors::StarTrackerConfig{});
    sensors::RangeSensor range_sensor(sensors::RangeSensorConfig{});
    const math::Vector3 target_pos{7000.0e3, 0.0, 0.0};

    // 3. Run truth simulation step-by-step while sampling sensors
    std::vector<spacecraft::SpacecraftState> truth_traj_with_sensors;
    {
        double current_t = 0.0;
        spacecraft::SpacecraftState state = initial_state;
        truth_traj_with_sensors.push_back(state);
        for (int step = 0; step < total_steps; ++step) {
            // Take sensor measurements
            static_cast<void>(imu.sample(current_t, state));
            static_cast<void>(gnss.sample(current_t, state));
            static_cast<void>(star_tracker.sample(current_t, state));
            static_cast<void>(range_sensor.sample(current_t, state, target_pos));

            state = spacecraft::rk4_step_spacecraft(current_t, state, dt, params);
            current_t += dt;
            truth_traj_with_sensors.push_back(state);
        }
    }

    // 4. Verify truth states are 100% BITWISE IDENTICAL
    REQUIRE(truth_traj_standalone.size() == truth_traj_with_sensors.size());
    for (std::size_t i = 0; i < truth_traj_standalone.size(); ++i) {
        const auto& s_alone = truth_traj_standalone[i];
        const auto& s_sens = truth_traj_with_sensors[i];

        CHECK(s_alone.position().x() == s_sens.position().x());
        CHECK(s_alone.position().y() == s_sens.position().y());
        CHECK(s_alone.position().z() == s_sens.position().z());

        CHECK(s_alone.velocity().x() == s_sens.velocity().x());
        CHECK(s_alone.velocity().y() == s_sens.velocity().y());
        CHECK(s_alone.velocity().z() == s_sens.velocity().z());

        CHECK(s_alone.orientation().w() == s_sens.orientation().w());
        CHECK(s_alone.orientation().x() == s_sens.orientation().x());
        CHECK(s_alone.orientation().y() == s_sens.orientation().y());
        CHECK(s_alone.orientation().z() == s_sens.orientation().z());

        CHECK(s_alone.angular_velocity_rad_per_s().x() == s_sens.angular_velocity_rad_per_s().x());
        CHECK(s_alone.angular_velocity_rad_per_s().y() == s_sens.angular_velocity_rad_per_s().y());
        CHECK(s_alone.angular_velocity_rad_per_s().z() == s_sens.angular_velocity_rad_per_s().z());
    }
}

TEST_CASE("Defensive validation rejects invalid sensor configurations and non-finite truth states", "[sensors][defense]") {
    // 1. Negative sample periods
    sensors::ImuConfig invalid_imu;
    invalid_imu.sample_period_s = -0.01;
    CHECK_THROWS_AS(sensors::ImuSensor(invalid_imu), std::domain_error);

    sensors::GnssConfig invalid_gnss;
    invalid_gnss.sample_period_s = 0.0;
    CHECK_THROWS_AS(sensors::GnssSensor(invalid_gnss), std::domain_error);

    sensors::StarTrackerConfig invalid_st;
    invalid_st.sample_period_s = -1.0;
    CHECK_THROWS_AS(sensors::StarTrackerSensor(invalid_st), std::domain_error);

    sensors::RangeSensorConfig invalid_range;
    invalid_range.sample_period_s = -0.1;
    CHECK_THROWS_AS(sensors::RangeSensor(invalid_range), std::domain_error);

    // 2. Negative noise standard deviations
    sensors::ImuConfig invalid_noise_imu;
    invalid_noise_imu.gyro_noise_std_rad_s = -0.1;
    CHECK_THROWS_AS(sensors::ImuSensor(invalid_noise_imu), std::domain_error);

    // 3. Non-finite truth states
    sensors::ImuSensor valid_imu(sensors::ImuConfig{});
    spacecraft::SpacecraftState nan_state;
    nan_state.translational.position = math::Vector3{std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0};
    CHECK_THROWS_AS(valid_imu.measure(0.0, nan_state), std::domain_error);
}
