cmake_minimum_required(VERSION 3.10)
project(FEM_2d CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Option to build shared library
option(BUILD_SHARED_LIB "Build FEM as a shared library" ON)
option(BUILD_PYTHON_BINDINGS "Build Python bindings with pybind11" OFF)

# Core FEM library sources
set(FEM_CORE_SOURCES
    barsystem.cpp
    src/fem_core.cpp
    src/fem_api.cpp
)

set(FEM_CORE_HEADERS
    barsystem.h
    src/fem_types.h
    src/fem_core.h
    src/fem_api.h
)

# Build shared library if enabled
if(BUILD_SHARED_LIB)
    add_library(fem_core SHARED ${FEM_CORE_SOURCES} ${FEM_CORE_HEADERS})
    target_include_directories(fem_core PUBLIC
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/src
    )
    if (MSVC)
        target_compile_options(fem_core PRIVATE /utf-8)
    endif()

    # Set library properties
    set_target_properties(fem_core PROPERTIES
        VERSION 1.0.0
        SOVERSION 1
        PUBLIC_HEADER "${FEM_CORE_HEADERS}"
    )

    # Install library
    install(TARGETS fem_core
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        PUBLIC_HEADER DESTINATION include/fem
    )
endif()

# Build standalone executable (original version)
add_executable(fem_run barsystem.cpp)
target_compile_definitions(fem_run PRIVATE)
target_include_directories(fem_run PRIVATE ${CMAKE_SOURCE_DIR})
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

# Python bindings with pybind11
if(BUILD_PYTHON_BINDINGS)
    find_package(pybind11 CONFIG)
    if(pybind11_FOUND)
        pybind11_add_module(pyfem python/pyfem_bindings.cpp ${FEM_CORE_SOURCES})
        target_include_directories(pyfem PRIVATE
            ${CMAKE_SOURCE_DIR}
            ${CMAKE_SOURCE_DIR}/src
        )
        if (MSVC)
            target_compile_options(pyfem PRIVATE /utf-8)
        endif()

        install(TARGETS pyfem DESTINATION python)
    else()
        message(WARNING "pybind11 not found. Python bindings will not be built.")
    endif()
endif()

enable_testing()
add_test(NAME basic_tests COMMAND unit_tests)
