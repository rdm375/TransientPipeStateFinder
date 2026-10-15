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

add_library(tpsf_forward STATIC src/forward.cpp)
target_include_directories(tpsf_forward PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(tpsf_forward PUBLIC pipe_sim_cpp_core)
target_compile_features(tpsf_forward PUBLIC cxx_std_20)

add_executable(tpsf_smoke tests/cpp/forward_smoke.cpp)
target_link_libraries(tpsf_smoke PRIVATE tpsf_forward)
if(BUILD_TESTING)
  add_test(NAME tpsf_forward_smoke COMMAND tpsf_smoke)
endif()
