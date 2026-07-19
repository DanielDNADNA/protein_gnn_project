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
# C++ test executables
# -----------------------------

add_executable(kd_tree_test
    tests/kd_tree_test_main.cpp
)

add_executable(pipeline_test
    tests/pipeline_test_main.cpp
    src/protein_gnn_project.cpp    # Needed if pipeline_test uses Parse_mmCIF
)

add_custom_command(TARGET kd_tree_test POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_SOURCE_DIR}/tests/sample.txt
        $<TARGET_FILE_DIR:kd_tree_test>
)

set(TEST_DATA
    1CRN.cif.gz
    1TIM.cif.gz
    6QNR.cif
    6QNR.cif.gz
    8H0R.cif.gz
    3J3Q.cif.gz
    3J3Y.cif.gz
    4V9D.cif.gz
    3E7Y.cif.gz
    1DSB.cif.gz
)

foreach(FILE ${TEST_DATA})
    add_custom_command(TARGET pipeline_test POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_SOURCE_DIR}/tests/data/${FILE}
            $<TARGET_FILE_DIR:pipeline_test>
    )
endforeach()

target_link_libraries(kd_tree_test PRIVATE
    spatial_core
)

# -----------------------------
# Dependencies
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

target_link_libraries(pipeline_test PRIVATE
    spatial_core
    gemmi_cpp
    ZLIB::ZLIB
    pybind11::embed
    Python3::Python
)

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