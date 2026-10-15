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
    src/forecast_continuation.cpp
    src/objective.cpp
    src/full_space_observations.cpp
    src/physical_boundary_channels.cpp
    src/full_space_constraints.cpp
    src/full_space_tangent.cpp
    src/full_space_adjoint.cpp
    src/full_space_trajectory.cpp
    src/full_space_trajectory_derivatives.cpp
    src/full_space_measurement_objective.cpp
    src/full_space_problem.cpp
    src/matrix_free_kkt.cpp
    src/minres.cpp
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


# E2.4b.3 — Coupled full-space constraints
add_executable(tpsf_full_space_constraints_test
    tests/cpp/full_space_constraints_test.cpp)

target_link_libraries(tpsf_full_space_constraints_test
    PRIVATE tpsf_forward)

if(BUILD_TESTING)
    add_test(
        NAME tpsf_full_space_constraints
        COMMAND tpsf_full_space_constraints_test)
endif()

# E2.4b.4 — Exact full-space discrete tangent-linear verification
add_executable(tpsf_full_space_tangent_test
    tests/cpp/full_space_tangent_test.cpp)

target_link_libraries(tpsf_full_space_tangent_test
    PRIVATE tpsf_forward)

if(BUILD_TESTING)
    add_test(
        NAME tpsf_full_space_tangent
        COMMAND tpsf_full_space_tangent_test)
endif()

# E2.4b.5 — Exact discrete full-space adjoint verification
add_executable(tpsf_full_space_adjoint_test
    tests/cpp/full_space_adjoint_test.cpp)

target_link_libraries(tpsf_full_space_adjoint_test
    PRIVATE tpsf_forward)

if(BUILD_TESTING)
    add_test(
        NAME tpsf_full_space_adjoint
        COMMAND tpsf_full_space_adjoint_test)
endif()

# E2.4b.6 — Trajectory-level full-space constraints
add_executable(tpsf_full_space_trajectory_test
    tests/cpp/full_space_trajectory_test.cpp)

target_link_libraries(tpsf_full_space_trajectory_test
    PRIVATE tpsf_forward)

if(BUILD_TESTING)
    add_test(
        NAME tpsf_full_space_trajectory
        COMMAND tpsf_full_space_trajectory_test)
endif()


# E2.4b.7 — Scored forward trajectory recording
add_executable(tpsf_scored_trajectory_test
    tests/cpp/scored_trajectory_test.cpp)

target_link_libraries(tpsf_scored_trajectory_test
    PRIVATE tpsf_forward)

if(BUILD_TESTING)
    add_test(
        NAME tpsf_scored_trajectory
        COMMAND tpsf_scored_trajectory_test)
endif()


# E2.4b.8 — Exact trajectory-level tangent and adjoint
add_executable(tpsf_full_space_trajectory_derivatives_test
    tests/cpp/full_space_trajectory_derivatives_test.cpp)

target_link_libraries(
    tpsf_full_space_trajectory_derivatives_test
    PRIVATE tpsf_forward)

if(BUILD_TESTING)
    add_test(
        NAME tpsf_full_space_trajectory_derivatives
        COMMAND tpsf_full_space_trajectory_derivatives_test)
endif()


# E2.4b.9 — Masked weighted boundary measurement objective
add_executable(tpsf_full_space_measurement_objective_test
    tests/cpp/full_space_measurement_objective_test.cpp)

target_link_libraries(
    tpsf_full_space_measurement_objective_test
    PRIVATE tpsf_forward)

if(BUILD_TESTING)
    add_test(
        NAME tpsf_full_space_measurement_objective
        COMMAND tpsf_full_space_measurement_objective_test)
endif()


# E2.4b.10 — Unified scaled full-space optimization problem
add_executable(tpsf_full_space_problem_test
    tests/cpp/full_space_problem_test.cpp)

target_link_libraries(
    tpsf_full_space_problem_test
    PRIVATE tpsf_forward)

if(BUILD_TESTING)
    add_test(
        NAME tpsf_full_space_problem
        COMMAND tpsf_full_space_problem_test)
endif()


# E2.4b.11a — Dense reference constrained least-squares solver.
# Eigen is intentionally not a dependency of tpsf_forward.
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_library(tpsf_dense_solver STATIC
    src/dense_constrained_least_squares.cpp
    src/dense_full_space_subproblem.cpp
    src/dense_full_space_sqp.cpp)

target_include_directories(tpsf_dense_solver
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(tpsf_dense_solver
    PUBLIC
        Eigen3::Eigen
        tpsf_forward)

target_compile_features(tpsf_dense_solver
    PUBLIC cxx_std_20)


add_executable(tpsf_dense_constrained_least_squares_test
    tests/cpp/dense_constrained_least_squares_test.cpp)

target_link_libraries(tpsf_dense_constrained_least_squares_test
    PRIVATE tpsf_dense_solver)

if(BUILD_TESTING)
    add_test(
        NAME tpsf_dense_constrained_least_squares
        COMMAND tpsf_dense_constrained_least_squares_test)
endif()


# E2.4b.11b.1 — Verified dense full-space SQP assembly.
add_executable(tpsf_dense_full_space_subproblem_test
    tests/cpp/dense_full_space_subproblem_test.cpp)

target_link_libraries(tpsf_dense_full_space_subproblem_test
    PRIVATE tpsf_dense_solver tpsf_forward)

if(BUILD_TESTING)
    add_test(
        NAME tpsf_dense_full_space_subproblem
        COMMAND tpsf_dense_full_space_subproblem_test)
endif()

# E2.4b.11b.2 — Nonlinear dense SQP reference solver.
add_executable(tpsf_dense_full_space_sqp_test
    tests/cpp/dense_full_space_sqp_test.cpp)

target_link_libraries(tpsf_dense_full_space_sqp_test
    PRIVATE tpsf_dense_solver)

if(BUILD_TESTING)
    add_test(
        NAME tpsf_dense_full_space_sqp
        COMMAND tpsf_dense_full_space_sqp_test)
endif()

# E2.4b.11b.3 — Forward-generated synthetic reconstruction.
add_executable(tpsf_synthetic_pipeline_reconstruction_test
    tests/cpp/synthetic_pipeline_reconstruction_test.cpp)

target_link_libraries(
    tpsf_synthetic_pipeline_reconstruction_test
    PRIVATE
        tpsf_dense_solver
        tpsf_forward)

if(BUILD_TESTING)
    add_test(
        NAME tpsf_synthetic_pipeline_reconstruction
        COMMAND tpsf_synthetic_pipeline_reconstruction_test)
endif()

# E2.4b.11b.4a — Unknown initial outlet-pressure reconstruction.
add_executable(tpsf_unknown_initial_outlet_pressure_test
    tests/cpp/unknown_initial_outlet_pressure_test.cpp)

target_link_libraries(
    tpsf_unknown_initial_outlet_pressure_test
    PRIVATE
        tpsf_dense_solver
        tpsf_forward)

if(BUILD_TESTING)
    add_test(
        NAME tpsf_unknown_initial_outlet_pressure
        COMMAND tpsf_unknown_initial_outlet_pressure_test)
endif()

# E2.4b.11b.4b — Independent initial-pressure reconstruction.
add_executable(tpsf_independent_initial_pressure_test
    tests/cpp/independent_initial_pressure_test.cpp)

target_link_libraries(
    tpsf_independent_initial_pressure_test
    PRIVATE
        tpsf_dense_solver
        tpsf_forward)

if(BUILD_TESTING)
    add_test(
        NAME tpsf_independent_initial_pressure_direct
        COMMAND tpsf_independent_initial_pressure_test 0)

    add_test(
        NAME tpsf_independent_initial_pressure_delayed
        COMMAND tpsf_independent_initial_pressure_test 2)
endif()


# E3 — Forecast continuation from reconstructed states
add_executable(tpsf_forecast_continuation_test
    tests/cpp/forecast_continuation_test.cpp)

target_link_libraries(tpsf_forecast_continuation_test
    PRIVATE tpsf_forward)

if(BUILD_TESTING)
    add_test(
        NAME tpsf_forecast_continuation
        COMMAND tpsf_forecast_continuation_test)
endif()

# E4.2b — Dense diagnostic oracle for Jacobian structure,
# numerical rank, and local identifiability.
add_executable(tpsf_kkt_structure_diagnostic
    tests/cpp/kkt_structure_diagnostic.cpp)

target_link_libraries(tpsf_kkt_structure_diagnostic
    PRIVATE tpsf_dense_solver)

if(BUILD_TESTING)
    add_test(
        NAME tpsf_kkt_structure_diagnostic
        COMMAND tpsf_kkt_structure_diagnostic)
endif()

# E4.3b.1 — Matrix-free symmetric Krylov solver.
add_executable(tpsf_minres_test
    tests/cpp/minres_test.cpp)

target_link_libraries(tpsf_minres_test
    PRIVATE tpsf_forward Eigen3::Eigen)

if(BUILD_TESTING)
    add_test(
        NAME tpsf_minres
        COMMAND tpsf_minres_test)
endif()
