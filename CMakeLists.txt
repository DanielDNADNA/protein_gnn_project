cmake_minimum_required(VERSION 3.20)

project(protein_spatial_gnn CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(
    $<$<CONFIG:Release>:-O3>
    -march=native
)

# -----------------------------
# Dependencies
# -----------------------------

include(FetchContent)

set(GEMMI_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GEMMI_EXECUTABLES OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)

FetchContent_Declare(
    gemmi
    GIT_REPOSITORY https://github.com/project-gemmi/gemmi.git
    GIT_TAG v0.6.5
)

FetchContent_MakeAvailable(gemmi)

find_package(ZLIB REQUIRED)

# Python only needed for pybind module
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG v2.13.6
)

FetchContent_MakeAvailable(pybind11)


# -----------------------------
# Core C++ library
# -----------------------------

add_library(spatial_core
    src/kd_tree.cpp
    src/cell_list.cpp
)

target_include_directories(spatial_core PUBLIC
    src
)


# -----------------------------
# Pipeline test
# -----------------------------

add_executable(pipeline_test
    tests/pipeline_test_main.cpp
    src/protein_gnn_project.cpp
)

target_link_libraries(pipeline_test PRIVATE
    spatial_core
    gemmi_cpp
    ZLIB::ZLIB
)

target_include_directories(pipeline_test PRIVATE
    src
)


add_custom_command(TARGET pipeline_test POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/data/benchmark_raw
        $<TARGET_FILE_DIR:pipeline_test>/benchmark_raw
)


# -----------------------------
# Evaluation benchmark
# -----------------------------

add_executable(evaluate
    tests/evaluation.cpp
    src/protein_gnn_project.cpp
)

target_link_libraries(evaluate PRIVATE
    spatial_core
    gemmi_cpp
    ZLIB::ZLIB
)

target_include_directories(evaluate PRIVATE
    src
)


# Copy preprocessed binary benchmark files
add_custom_command(TARGET evaluate POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/data/benchmark_bin
        $<TARGET_FILE_DIR:evaluate>/benchmark_bin
)


# -----------------------------
# Python extension module
# -----------------------------

pybind11_add_module(spatial_engine
    src/bindings.cpp
    src/protein_gnn_project.cpp
)

target_link_libraries(spatial_engine PRIVATE
    spatial_core
    pybind11::module
    Python3::Module
    gemmi_cpp
    ZLIB::ZLIB
)

target_include_directories(spatial_engine PRIVATE
    src
)

set_target_properties(spatial_engine PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}
)


# -----------------------------
# Debug info
# -----------------------------

message(STATUS "Project: ${CMAKE_PROJECT_NAME}")

get_property(allTargets DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)

message(STATUS "Targets: ${allTargets}")