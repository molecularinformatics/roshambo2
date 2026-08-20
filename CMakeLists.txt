cmake_minimum_required(VERSION 3.24)  # 3.24+ for CMAKE_CUDA_ARCHITECTURES=native
option(BUILD_WITH_CUDA "Build the CUDA functionaility" ON)

# NOTE: CUDA is deliberately not listed as a project language here. Enabling a
# language commits CMake to its own default for CMAKE_CUDA_ARCHITECTURES at
# that moment, so setting the variable afterwards is silently ignored -- the
# guard `if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)` is already false by then.
# CUDA is enabled below, after our own default has been set.
project(_roshambo2 LANGUAGES CXX)

# Default to Release. Without this the manual cmake build in DEV_GUIDE.md
# leaves CMAKE_BUILD_TYPE empty, which compiles the host code at -O0 with
# asserts enabled.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

# Opt-in for fast local iteration: target only the GPUs present on this
# machine. Off by default because the result is not portable.
option(ROSHAMBO_ARCH_NATIVE "Target only the CUDA architecture(s) detected on the build machine (fast dev builds, not portable)" OFF)

if (BUILD_WITH_CUDA)
    # Must be set before CUDA is enabled, see the note above project().
    # SASS for each listed architecture plus compute_120 PTX, so the driver can
    # JIT for GPUs newer than any listed here. Override with
    # -DCMAKE_CUDA_ARCHITECTURES=... to target something else.
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        if(ROSHAMBO_ARCH_NATIVE)
            set(CMAKE_CUDA_ARCHITECTURES native)
        else()
            set(CMAKE_CUDA_ARCHITECTURES "80-real;86-real;89-real;90-real;120-real;120-virtual")
        endif()
    endif()
    enable_language(CUDA)
endif (BUILD_WITH_CUDA)

# Find OpenMP package
find_package(OpenMP REQUIRED)

# Find pybind
find_package(pybind11 REQUIRED)

if (BUILD_WITH_CUDA)
    # CUDAToolkit replaces the deprecated FindCUDA module. The CUDA::cudart
    # target carries the toolkit include directories, so cuda_interface.cpp
    # finds the CUDA headers without a manual include_directories() call --
    # which was being handed an empty string on newer CMake.
    find_package(CUDAToolkit REQUIRED)
endif (BUILD_WITH_CUDA)

# TODO: optional GPU build
pybind11_add_module(_roshambo2_cpp roshambo2/backends/cpp_src/cpp_functions.cpp roshambo2/backends/cpp_src/cpp_helper_functions.cpp)
if (BUILD_WITH_CUDA)
    pybind11_add_module(_roshambo2_cuda roshambo2/backends/cuda_src/cuda_functions.cu roshambo2/backends/cuda_src/cuda_interface.cpp)
    target_link_libraries(_roshambo2_cuda PRIVATE OpenMP::OpenMP_CXX CUDA::cudart)
endif (BUILD_WITH_CUDA)

target_link_libraries(_roshambo2_cpp PRIVATE OpenMP::OpenMP_CXX)

install(TARGETS _roshambo2_cpp
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX})

if (BUILD_WITH_CUDA)
install(TARGETS _roshambo2_cuda
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX})
endif (BUILD_WITH_CUDA)
