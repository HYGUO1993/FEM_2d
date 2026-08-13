cmake_minimum_required(VERSION 3.10)
project(FEM_2d CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# femcore 静态库:求解内核,不含 main(fem_run/femcli/unit_tests 共用)
add_library(femcore STATIC barsystem.cpp)
target_include_directories(femcore PUBLIC ${CMAKE_SOURCE_DIR})
# 排除 barsystem.cpp 中的 main, 避免与 femcli/fem_run 的 main 冲突
target_compile_definitions(femcore PRIVATE FEM_NO_MAIN)
if (MSVC)
  target_compile_options(femcore PRIVATE /utf-8)
endif()

# Build standalone executable (txt 模式, 向后兼容)
add_executable(fem_run barsystem.cpp)
target_include_directories(fem_run PRIVATE ${CMAKE_SOURCE_DIR})
if (MSVC)
  target_compile_options(fem_run PRIVATE /utf-8)
endif()

# femcli: JSON 模式命令行工具
add_executable(femcli femcli.cpp)
target_link_libraries(femcli PRIVATE femcore)
target_include_directories(femcli PRIVATE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/third_party)
if (MSVC)
  target_compile_options(femcli PRIVATE /utf-8)
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
