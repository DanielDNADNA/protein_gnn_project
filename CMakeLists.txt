cmake_minimum_required(VERSION 3.20)

project(protein_spatial_gnn CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(
    $<$<CONFIG:Release>:-O3>
    -march=native
)

find_package(OpenMP REQUIRED)

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

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG v2.13.6
)

FetchContent_MakeAvailable(pybind11)

# -----------------------------
# Core library
# -----------------------------

add_library(spatial_core
    src/kd_tree.cpp
    src/cell_list.cpp
)

target_include_directories(spatial_core PUBLIC
    src
)

target_link_libraries(spatial_core
    PUBLIC
    OpenMP::OpenMP_CXX
)

# -----------------------------
# Helper function for executables
# -----------------------------

function(add_spatial_executable target)

    add_executable(${target}
        ${ARGN}
        src/protein_gnn_project.cpp
    )

    target_link_libraries(${target}
        PRIVATE
        spatial_core
        gemmi_cpp
        ZLIB::ZLIB
        pybind11::embed
        Python3::Python
    )

    target_include_directories(${target}
        PRIVATE
        src
        ${Python3_INCLUDE_DIRS}
    )

endfunction()

# -----------------------------
# Pipeline test
# -----------------------------

add_spatial_executable(
    pipeline_test
    tests/pipeline_test_main.cpp
)

add_custom_command(TARGET pipeline_test POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/data/benchmark_raw
        $<TARGET_FILE_DIR:pipeline_test>/benchmark_raw
)

# -----------------------------
# Evaluation benchmark
# -----------------------------

add_spatial_executable(
    evaluate
    tests/evaluation.cpp
    tests/loader.cpp
)

add_custom_command(TARGET evaluate POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/data/benchmark_bin
        $<TARGET_FILE_DIR:evaluate>/benchmark_bin
)

# -----------------------------
# Save benchmark parser
# -----------------------------

add_spatial_executable(
    save_pdb_files
    tests/save_pdb_files.cpp
    tests/evaluation.cpp
)

add_custom_command(TARGET save_pdb_files POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/data/benchmark_raw
        $<TARGET_FILE_DIR:save_pdb_files>/benchmark_raw
)

target_compile_definitions(save_pdb_files PRIVATE
    PROJECT_DATA_DIR="${CMAKE_SOURCE_DIR}/data"
    BENCHMARK_RAW_DIR="${CMAKE_SOURCE_DIR}/data/benchmark_raw"
    BENCHMARK_BIN_DIR="${CMAKE_SOURCE_DIR}/data/benchmark_bin"
)

# -----------------------------
# Python extension
# -----------------------------

pybind11_add_module(
    spatial_engine
    src/bindings.cpp
    src/protein_gnn_project.cpp
)

target_link_libraries(spatial_engine
    PRIVATE
    spatial_core
    gemmi_cpp
    ZLIB::ZLIB
    pybind11::module
    Python3::Module
)

target_include_directories(spatial_engine
    PRIVATE
    src
    ${Python3_INCLUDE_DIRS}
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