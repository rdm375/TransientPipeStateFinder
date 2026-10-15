cmake_minimum_required(VERSION 3.20)
project(TransientPipeStateFinder VERSION 0.1.0 LANGUAGES CXX)
include(CTest)

set(PIPE_SIM_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../isothermal-pipe-sim"
    CACHE PATH "Path to the validated isothermal-pipe-sim checkout")
if(NOT EXISTS "${PIPE_SIM_SOURCE_DIR}/CMakeLists.txt")
  message(FATAL_ERROR "isothermal-pipe-sim missing at ${PIPE_SIM_SOURCE_DIR}; provide -DPIPE_SIM_SOURCE_DIR=/path/to/checkout")
endif()

# No physics equations are copied or reimplemented here.
add_subdirectory("${PIPE_SIM_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/pipe_sim")
if(NOT TARGET pipe_sim_cpp_core)
  message(FATAL_ERROR "Expected target pipe_sim_cpp_core missing from forward solver")
endif()

add_library(tpsf_forward STATIC
    src/forward.cpp
    src/scheduled_forward.cpp
    src/objective.cpp
    src/full_space_observations.cpp
    src/physical_boundary_channels.cpp
)
target_include_directories(tpsf_forward PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(tpsf_forward PUBLIC pipe_sim_cpp_core)
target_compile_features(tpsf_forward PUBLIC cxx_std_20)

add_executable(tpsf_smoke tests/cpp/forward_smoke.cpp)
target_link_libraries(tpsf_smoke PRIVATE tpsf_forward)
if(BUILD_TESTING)
  add_test(NAME tpsf_forward_smoke COMMAND tpsf_smoke)
endif()

# E2.1 — Scaled initial-state representation
add_executable(tpsf_scaled_state_test
    tests/cpp/scaled_state_test.cpp)

target_link_libraries(tpsf_scaled_state_test
    PRIVATE tpsf_forward)

if(BUILD_TESTING)
    add_test(
        NAME tpsf_scaled_state
        COMMAND tpsf_scaled_state_test)
endif()

# E2.2 — Boundary observation model






# E2.2 — Four-channel boundary observations

# E2.2 — Four-channel boundary observations
add_executable(tpsf_observations_test
    tests/cpp/observations_test.cpp)

target_link_libraries(tpsf_observations_test
    PRIVATE tpsf_forward)

if(BUILD_TESTING)
    add_test(
        NAME tpsf_observations
        COMMAND tpsf_observations_test)
endif()


# E2.3 — Scheduled boundary integration
add_executable(tpsf_scheduled_forward_test
    tests/cpp/scheduled_forward_test.cpp)

target_link_libraries(tpsf_scheduled_forward_test
    PRIVATE tpsf_forward)

if(BUILD_TESTING)
    add_test(
        NAME tpsf_scheduled_forward
        COMMAND tpsf_scheduled_forward_test)
endif()


# E2.4a — Four-channel objective evaluation
add_executable(tpsf_objective_test
    tests/cpp/objective_test.cpp)

target_link_libraries(tpsf_objective_test
    PRIVATE tpsf_forward)

if(BUILD_TESTING)
    add_test(
        NAME tpsf_objective
        COMMAND tpsf_objective_test)
endif()


# E2.4b.1 — Full-space observation variables
add_executable(tpsf_full_space_observations_test
    tests/cpp/full_space_observations_test.cpp)

target_link_libraries(tpsf_full_space_observations_test
    PRIVATE tpsf_forward)

if(BUILD_TESTING)
    add_test(
        NAME tpsf_full_space_observations
        COMMAND tpsf_full_space_observations_test)
endif()


# E2.4b.2 — Physical boundary-channel extraction
add_executable(tpsf_physical_boundary_channels_test
    tests/cpp/physical_boundary_channels_test.cpp)

target_link_libraries(tpsf_physical_boundary_channels_test
    PRIVATE tpsf_forward)

if(BUILD_TESTING)
    add_test(
        NAME tpsf_physical_boundary_channels
        COMMAND tpsf_physical_boundary_channels_test)
endif()
