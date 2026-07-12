cmake_minimum_required(VERSION 3.20)

# 1. Point to the vcpkg toolchain

# 2. Define the single project name
project(protein_spatial_gnn CXX)

# 3. Set global C++20 standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(GEMMI_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

# 4. Fetch the Gemmi library
include(FetchContent)

set(GEMMI_EXECUTABLES OFF CACHE BOOL "Don't build gemmi programs" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "Don't build tests" FORCE)

FetchContent_Declare(
    gemmi
    GIT_REPOSITORY https://github.com/project-gemmi/gemmi.git
    GIT_TAG        v0.6.5 
) 

set(USE_PYTHON ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(gemmi)

# 5. Find our external packages
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(ZLIB REQUIRED)

FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG v2.13.6
)
FetchContent_MakeAvailable(pybind11)


# 6. Define our files
set(SOURCE_FILES
    src/protein_gnn_project.cpp
    src/protein_gnn_project.h
    src/kd_tree.cpp
    src/kd_tree.h
    src/bindings.cpp
    src/atom.h
    src/edges.h
)

# 7. Build the Pybind11 Python extension module
pybind11_add_module(spatial_engine ${SOURCE_FILES})

message(STATUS "Targets: ${CMAKE_PROJECT_NAME}")
get_property(allTargets DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)
message(STATUS "All targets: ${allTargets}")

# 8. Link the external libraries to module
target_link_libraries(spatial_engine PRIVATE 
    pybind11::module
    Python3::Module
    gemmi_cpp
    ZLIB::ZLIB
)

# Tells CMake to output the final .so file to your project root
set_target_properties(spatial_engine PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})