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

add_executable(astradock_orbit_demo tools/orbit_demo.cpp)
target_link_libraries(astradock_orbit_demo PRIVATE AstraDock::orbit)

if(MSVC)
    target_compile_options(astradock_orbit_demo PRIVATE /W4 /permissive-)
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

    if(MSVC)
        target_compile_options(astradock_vector3_tests PRIVATE /W4 /permissive-)
        target_compile_options(astradock_two_body_tests PRIVATE /W4 /permissive-)
        target_compile_options(astradock_integrator_tests PRIVATE /W4 /permissive-)
        target_compile_options(astradock_orbit_tests PRIVATE /W4 /permissive-)
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
    endif()

    include(Catch)
    catch_discover_tests(astradock_vector3_tests)
    catch_discover_tests(astradock_two_body_tests)
    catch_discover_tests(astradock_integrator_tests)
    catch_discover_tests(astradock_orbit_tests)
endif()
