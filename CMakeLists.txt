cmake_minimum_required(VERSION 3.10)
project(FEM_2d CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(fem_run barsystem.cpp)

# Unit tests build with UNIT_TEST defined (excludes main)
add_executable(unit_tests tests/test_basic.cpp barsystem.cpp)
target_compile_definitions(unit_tests PRIVATE UNIT_TEST)
target_include_directories(fem_run PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(unit_tests PRIVATE ${CMAKE_SOURCE_DIR})
if (MSVC)
  target_compile_options(fem_run PRIVATE /utf-8)
  target_compile_options(unit_tests PRIVATE /utf-8)
endif()

enable_testing()
add_test(NAME basic_tests COMMAND unit_tests)
