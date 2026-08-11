cmake_minimum_required(VERSION 3.10)
project(FEM_2d VERSION 1.0.1 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Generate version header from template
configure_file(
  ${CMAKE_SOURCE_DIR}/version.h.in
  ${CMAKE_BINARY_DIR}/generated/version.h
  @ONLY
)

# Build standalone executable
add_executable(fem_run barsystem.cpp)
target_include_directories(fem_run PRIVATE ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR}/generated)
if (MSVC)
  target_compile_options(fem_run PRIVATE /utf-8)
endif()

# Unit tests build with UNIT_TEST defined (excludes main)
add_executable(unit_tests tests/test_basic.cpp barsystem.cpp)
target_compile_definitions(unit_tests PRIVATE UNIT_TEST)
target_include_directories(unit_tests PRIVATE ${CMAKE_SOURCE_DIR})
if (MSVC)
  target_compile_options(unit_tests PRIVATE /utf-8)
endif()

enable_testing()
add_test(NAME basic_tests COMMAND unit_tests)
