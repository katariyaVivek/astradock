cmake_minimum_required(VERSION 3.24)

project(
    AstraDock
    VERSION 0.1.0
    DESCRIPTION "An educational autonomous orbital rendezvous and docking simulator"
    LANGUAGES CXX
)

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_library(astradock_math INTERFACE)
add_library(AstraDock::math ALIAS astradock_math)

target_compile_features(astradock_math INTERFACE cxx_std_20)
target_include_directories(
    astradock_math
    INTERFACE
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/cpp>"
)

add_library(astradock_dynamics INTERFACE)
add_library(AstraDock::dynamics ALIAS astradock_dynamics)
target_link_libraries(astradock_dynamics INTERFACE AstraDock::math)

add_library(astradock_numerics INTERFACE)
add_library(AstraDock::numerics ALIAS astradock_numerics)
target_link_libraries(astradock_numerics INTERFACE AstraDock::math)

add_library(astradock_orbit INTERFACE)
add_library(AstraDock::orbit ALIAS astradock_orbit)
target_link_libraries(
    astradock_orbit
    INTERFACE
        AstraDock::dynamics
        AstraDock::numerics
)

add_library(astradock_frames INTERFACE)
add_library(AstraDock::frames ALIAS astradock_frames)
target_link_libraries(astradock_frames INTERFACE AstraDock::math)

add_library(astradock_attitude INTERFACE)
add_library(AstraDock::attitude ALIAS astradock_attitude)
target_link_libraries(
    astradock_attitude
    INTERFACE
        AstraDock::math
        AstraDock::numerics
)

add_library(astradock_spacecraft INTERFACE)
add_library(AstraDock::spacecraft ALIAS astradock_spacecraft)
target_link_libraries(
    astradock_spacecraft
    INTERFACE
        AstraDock::orbit
        AstraDock::attitude
        AstraDock::frames
)

add_library(astradock_environment INTERFACE)
add_library(AstraDock::environment ALIAS astradock_environment)
target_link_libraries(
    astradock_environment
    INTERFACE
        AstraDock::spacecraft
        AstraDock::orbit
        AstraDock::attitude
        AstraDock::frames
        AstraDock::dynamics
        AstraDock::numerics
        AstraDock::math
)

add_library(astradock_sensors INTERFACE)
add_library(AstraDock::sensors ALIAS astradock_sensors)
target_link_libraries(
    astradock_sensors
    INTERFACE
        AstraDock::spacecraft
        AstraDock::environment
        AstraDock::frames
        AstraDock::attitude
        AstraDock::orbit
        AstraDock::dynamics
        AstraDock::numerics
        AstraDock::math
)

add_library(astradock_estimation INTERFACE)
add_library(AstraDock::estimation ALIAS astradock_estimation)
target_link_libraries(
    astradock_estimation
    INTERFACE
        AstraDock::orbit
        AstraDock::math
)

add_executable(astradock_orbit_demo tools/orbit_demo.cpp)
target_link_libraries(astradock_orbit_demo PRIVATE AstraDock::orbit)

add_executable(astradock_validation_demo tools/validation_demo.cpp)
target_link_libraries(astradock_validation_demo PRIVATE AstraDock::orbit)

add_executable(astradock_frame_demo tools/frame_demo.cpp)
target_link_libraries(
    astradock_frame_demo
    PRIVATE
        AstraDock::frames
        AstraDock::orbit
)

add_executable(astradock_elements_demo tools/orbital_elements_demo.cpp)
target_link_libraries(
    astradock_elements_demo
    PRIVATE
        AstraDock::orbit
)

add_executable(astradock_attitude_demo tools/attitude_demo.cpp)
target_link_libraries(
    astradock_attitude_demo
    PRIVATE
        AstraDock::attitude
        AstraDock::frames
)

add_executable(astradock_attitude_dynamics_demo tools/attitude_dynamics_demo.cpp)
target_link_libraries(
    astradock_attitude_dynamics_demo
    PRIVATE
        AstraDock::attitude
        AstraDock::frames
)

add_executable(astradock_6dof_demo tools/six_dof_demo.cpp)
target_link_libraries(
    astradock_6dof_demo
    PRIVATE
        AstraDock::spacecraft
)

add_executable(astradock_environment_demo tools/environment_demo.cpp)
target_link_libraries(
    astradock_environment_demo
    PRIVATE
        AstraDock::environment
)

add_executable(astradock_sensor_demo tools/sensor_demo.cpp)
target_link_libraries(
    astradock_sensor_demo
    PRIVATE
        AstraDock::sensors
)

add_executable(astradock_ekf_demo tools/ekf_demo.cpp)
target_link_libraries(
    astradock_ekf_demo
    PRIVATE
        AstraDock::estimation
        AstraDock::sensors
)

add_executable(astradock_imu_ekf_demo tools/imu_ekf_demo.cpp)
target_link_libraries(
    astradock_imu_ekf_demo
    PRIVATE
        AstraDock::estimation
        AstraDock::sensors
        AstraDock::environment
)

add_executable(astradock_attitude_ekf_demo tools/attitude_ekf_demo.cpp)
target_link_libraries(
    astradock_attitude_ekf_demo
    PRIVATE
        AstraDock::estimation
        AstraDock::sensors
        AstraDock::attitude
        AstraDock::orbit
        AstraDock::spacecraft
)

add_executable(astradock_integrated_navigation_demo tools/integrated_navigation_demo.cpp)
target_link_libraries(
    astradock_integrated_navigation_demo
    PRIVATE
        AstraDock::estimation
        AstraDock::sensors
        AstraDock::attitude
        AstraDock::orbit
        AstraDock::spacecraft
)

if(MSVC)
    target_compile_options(astradock_orbit_demo PRIVATE /W4 /permissive-)
    target_compile_options(astradock_validation_demo PRIVATE /W4 /permissive-)
    target_compile_options(astradock_frame_demo PRIVATE /W4 /permissive-)
    target_compile_options(astradock_elements_demo PRIVATE /W4 /permissive-)
    target_compile_options(astradock_attitude_demo PRIVATE /W4 /permissive-)
    target_compile_options(astradock_attitude_dynamics_demo PRIVATE /W4 /permissive-)
    target_compile_options(astradock_6dof_demo PRIVATE /W4 /permissive-)
    target_compile_options(astradock_environment_demo PRIVATE /W4 /permissive-)
    target_compile_options(astradock_sensor_demo PRIVATE /W4 /permissive-)
    target_compile_options(astradock_ekf_demo PRIVATE /W4 /permissive-)
    target_compile_options(astradock_imu_ekf_demo PRIVATE /W4 /permissive-)
    target_compile_options(astradock_attitude_ekf_demo PRIVATE /W4 /permissive-)
    target_compile_options(astradock_integrated_navigation_demo PRIVATE /W4 /permissive-)
else()
    target_compile_options(
        astradock_orbit_demo
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
    )
    target_compile_options(
        astradock_validation_demo
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
    )
    target_compile_options(
        astradock_frame_demo
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
    )
    target_compile_options(
        astradock_elements_demo
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
    )
    target_compile_options(
        astradock_attitude_demo
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
    )
    target_compile_options(
        astradock_attitude_dynamics_demo
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
    )
    target_compile_options(
        astradock_6dof_demo
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
    )
    target_compile_options(
        astradock_attitude_ekf_demo
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
    )
endif()

include(CTest)

if(BUILD_TESTING)
    find_package(Catch2 3 CONFIG QUIET)

    if(NOT Catch2_FOUND)
        include(FetchContent)
        FetchContent_Declare(
            Catch2
            GIT_REPOSITORY https://github.com/catchorg/Catch2.git
            GIT_TAG v3.7.1
            GIT_SHALLOW TRUE
        )
        FetchContent_MakeAvailable(Catch2)
        list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/extras")
    endif()

    add_executable(astradock_vector3_tests tests/cpp/test_vector3.cpp)
    target_link_libraries(
        astradock_vector3_tests
        PRIVATE
            AstraDock::math
            Catch2::Catch2WithMain
    )

    add_executable(astradock_two_body_tests tests/cpp/test_two_body.cpp)
    target_link_libraries(
        astradock_two_body_tests
        PRIVATE
            AstraDock::dynamics
            Catch2::Catch2WithMain
    )

    add_executable(astradock_integrator_tests tests/cpp/test_integrators.cpp)
    target_link_libraries(
        astradock_integrator_tests
        PRIVATE
            AstraDock::numerics
            Catch2::Catch2WithMain
    )

    add_executable(astradock_orbit_tests tests/cpp/test_orbit_propagation.cpp)
    target_link_libraries(
        astradock_orbit_tests
        PRIVATE
            AstraDock::orbit
            Catch2::Catch2WithMain
    )

    add_executable(astradock_validation_tests tests/cpp/test_validation.cpp)
    target_link_libraries(
        astradock_validation_tests
        PRIVATE
            AstraDock::orbit
            Catch2::Catch2WithMain
    )

    add_executable(astradock_matrix3_tests tests/cpp/test_matrix3.cpp)
    target_link_libraries(
        astradock_matrix3_tests
        PRIVATE
            AstraDock::math
            Catch2::Catch2WithMain
    )

    add_executable(astradock_frame_tests tests/cpp/test_coordinate_frames.cpp)
    target_link_libraries(
        astradock_frame_tests
        PRIVATE
            AstraDock::frames
            AstraDock::orbit
            Catch2::Catch2WithMain
    )

    add_executable(astradock_audit_tests tests/cpp/test_m01_m06_audit.cpp)
    target_link_libraries(
        astradock_audit_tests
        PRIVATE
            AstraDock::frames
            AstraDock::orbit
            Catch2::Catch2WithMain
    )

    add_executable(astradock_elements_tests tests/cpp/test_classical_elements.cpp)
    target_link_libraries(
        astradock_elements_tests
        PRIVATE
            AstraDock::orbit
            Catch2::Catch2WithMain
    )

    add_executable(astradock_quaternion_tests tests/cpp/test_quaternion.cpp)
    target_link_libraries(
        astradock_quaternion_tests
        PRIVATE
            AstraDock::attitude
            AstraDock::frames
            Catch2::Catch2WithMain
    )

    add_executable(astradock_attitude_dynamics_tests tests/cpp/test_attitude_dynamics.cpp)
    target_link_libraries(
        astradock_attitude_dynamics_tests
        PRIVATE
            AstraDock::attitude
            AstraDock::frames
            Catch2::Catch2WithMain
    )

    add_executable(astradock_6dof_tests tests/cpp/test_6dof.cpp)
    target_link_libraries(
        astradock_6dof_tests
        PRIVATE
            AstraDock::spacecraft
            Catch2::Catch2WithMain
    )

    add_executable(astradock_environment_tests tests/cpp/test_environment.cpp)
    target_link_libraries(
        astradock_environment_tests
        PRIVATE
            AstraDock::environment
            Catch2::Catch2WithMain
    )

    add_executable(astradock_sensor_tests tests/cpp/test_sensors.cpp)
    target_link_libraries(
        astradock_sensor_tests
        PRIVATE
            AstraDock::sensors
            Catch2::Catch2WithMain
    )

    add_executable(astradock_estimation_tests tests/cpp/test_estimation.cpp)
    target_link_libraries(
        astradock_estimation_tests
        PRIVATE
            AstraDock::estimation
            AstraDock::sensors
            Catch2::Catch2WithMain
    )

    add_executable(astradock_imu_ekf_tests tests/cpp/test_imu_ekf.cpp)
    target_link_libraries(
        astradock_imu_ekf_tests
        PRIVATE
            AstraDock::estimation
            AstraDock::sensors
            AstraDock::environment
            Catch2::Catch2WithMain
    )

    add_executable(astradock_attitude_ekf_tests tests/cpp/test_attitude_ekf.cpp)
    target_link_libraries(
        astradock_attitude_ekf_tests
        PRIVATE
            AstraDock::estimation
            AstraDock::sensors
            AstraDock::attitude
            Catch2::Catch2WithMain
    )

    add_executable(astradock_range_tests tests/cpp/test_range_update.cpp)
    target_link_libraries(
        astradock_range_tests
        PRIVATE
            AstraDock::estimation
            AstraDock::sensors
            AstraDock::orbit
            Catch2::Catch2WithMain
    )

    add_executable(astradock_integrated_navigation_tests tests/cpp/test_integrated_navigation.cpp)
    target_link_libraries(
        astradock_integrated_navigation_tests
        PRIVATE
            AstraDock::estimation
            AstraDock::sensors
            AstraDock::attitude
            AstraDock::orbit
            Catch2::Catch2WithMain
    )

    if(MSVC)
        target_compile_options(astradock_vector3_tests PRIVATE /W4 /permissive-)
        target_compile_options(astradock_two_body_tests PRIVATE /W4 /permissive-)
        target_compile_options(astradock_integrator_tests PRIVATE /W4 /permissive-)
        target_compile_options(astradock_orbit_tests PRIVATE /W4 /permissive-)
        target_compile_options(astradock_validation_tests PRIVATE /W4 /permissive-)
        target_compile_options(astradock_matrix3_tests PRIVATE /W4 /permissive-)
        target_compile_options(astradock_frame_tests PRIVATE /W4 /permissive-)
        target_compile_options(astradock_audit_tests PRIVATE /W4 /permissive-)
        target_compile_options(astradock_elements_tests PRIVATE /W4 /permissive-)
        target_compile_options(astradock_quaternion_tests PRIVATE /W4 /permissive-)
        target_compile_options(astradock_attitude_dynamics_tests PRIVATE /W4 /permissive-)
        target_compile_options(astradock_6dof_tests PRIVATE /W4 /permissive-)
        target_compile_options(astradock_environment_tests PRIVATE /W4 /permissive-)
        target_compile_options(astradock_sensor_tests PRIVATE /W4 /permissive-)
        target_compile_options(astradock_estimation_tests PRIVATE /W4 /permissive-)
        target_compile_options(astradock_imu_ekf_tests PRIVATE /W4 /permissive-)
        target_compile_options(astradock_attitude_ekf_tests PRIVATE /W4 /permissive-)
        target_compile_options(astradock_range_tests PRIVATE /W4 /permissive-)
        target_compile_options(astradock_integrated_navigation_tests PRIVATE /W4 /permissive-)
    else()
        target_compile_options(
            astradock_vector3_tests
            PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
        )
        target_compile_options(
            astradock_two_body_tests
            PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
        )
        target_compile_options(
            astradock_integrator_tests
            PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
        )
        target_compile_options(
            astradock_orbit_tests
            PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
        )
        target_compile_options(
            astradock_validation_tests
            PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
        )
        target_compile_options(
            astradock_matrix3_tests
            PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
        )
        target_compile_options(
            astradock_frame_tests
            PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
        )
        target_compile_options(
            astradock_audit_tests
            PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
        )
        target_compile_options(
            astradock_elements_tests
            PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
        )
        target_compile_options(
            astradock_quaternion_tests
            PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
        )
        target_compile_options(
            astradock_attitude_dynamics_tests
            PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
        )
        target_compile_options(
            astradock_6dof_tests
            PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
        )
        target_compile_options(
            astradock_environment_tests
            PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
        )
        target_compile_options(
            astradock_sensor_tests
            PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
        )
        target_compile_options(
            astradock_estimation_tests
            PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
        )
        target_compile_options(
            astradock_imu_ekf_tests
            PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
        )
        target_compile_options(
            astradock_attitude_ekf_tests
            PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
        )
        target_compile_options(
            astradock_range_tests
            PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wsign-conversion
        )
    endif()

    include(Catch)
    catch_discover_tests(astradock_vector3_tests)
    catch_discover_tests(astradock_two_body_tests)
    catch_discover_tests(astradock_integrator_tests)
    catch_discover_tests(astradock_orbit_tests)
    catch_discover_tests(astradock_validation_tests)
    catch_discover_tests(astradock_matrix3_tests)
    catch_discover_tests(astradock_frame_tests)
    catch_discover_tests(astradock_audit_tests)
    catch_discover_tests(astradock_elements_tests)
    catch_discover_tests(astradock_quaternion_tests)
    catch_discover_tests(astradock_attitude_dynamics_tests)
    catch_discover_tests(astradock_6dof_tests)
    catch_discover_tests(astradock_environment_tests)
    catch_discover_tests(astradock_sensor_tests)
    catch_discover_tests(astradock_estimation_tests)
    catch_discover_tests(astradock_imu_ekf_tests)
    catch_discover_tests(astradock_attitude_ekf_tests)
    catch_discover_tests(astradock_range_tests)
    catch_discover_tests(astradock_integrated_navigation_tests)
endif()
