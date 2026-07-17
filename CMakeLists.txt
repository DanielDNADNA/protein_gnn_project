cmake_minimum_required(VERSION 3.20)

project(protein_spatial_gnn CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# -----------------------------
# Core C++ library
# -----------------------------

add_library(spatial_core
    src/kd_tree.cpp
)

target_include_directories(spatial_core PUBLIC
    src
)


# -----------------------------
# C++ KD-tree test executable
# -----------------------------

add_executable(kd_tree_test
    tests/main.cpp
)

add_custom_command(TARGET kd_tree_test POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/tests/sample.txt
    $<TARGET_FILE_DIR:kd_tree_test>
)

target_link_libraries(kd_tree_test PRIVATE
    spatial_core
)


# -----------------------------
# Dependencies for Python module
# -----------------------------

include(FetchContent)

set(GEMMI_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GEMMI_EXECUTABLES OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(USE_PYTHON ON CACHE BOOL "" FORCE)


FetchContent_Declare(
    gemmi
    GIT_REPOSITORY https://github.com/project-gemmi/gemmi.git
    GIT_TAG v0.6.5
)

FetchContent_MakeAvailable(gemmi)


find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(ZLIB REQUIRED)


FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG v2.13.6
)

FetchContent_MakeAvailable(pybind11)


# -----------------------------
# Python extension module
# -----------------------------

set(PYTHON_SOURCES
    src/protein_gnn_project.cpp
    src/protein_gnn_project.h
    src/bindings.cpp
    src/atom.h
    src/edges.h
)

pybind11_add_module(spatial_engine
    ${PYTHON_SOURCES}
)


target_link_libraries(spatial_engine PRIVATE
    spatial_core
    pybind11::module
    Python3::Module
    gemmi_cpp
    ZLIB::ZLIB
)


# Output .so into project root
set_target_properties(spatial_engine PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}
)


# Debug info
message(STATUS "Targets: ${CMAKE_PROJECT_NAME}")
get_property(allTargets DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)
message(STATUS "All targets: ${allTargets}")