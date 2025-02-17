cmake_minimum_required(VERSION 3.18)
project(_roshambo2 LANGUAGES CXX CUDA)

# Find OpenMP package
find_package(OpenMP REQUIRED)

# Find CUDA
find_package(CUDA REQUIRED)

# Find pybind
find_package(pybind11)

include_directories("${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")

# change here for newer or older cards
set(CMAKE_CUDA_ARCHITECTURES 60 70 80)


# TODO: optional GPU build
pybind11_add_module(_roshambo2_cpp roshambo2/backends/cpp_src/cpp_functions.cpp)
pybind11_add_module(_roshambo2_cuda roshambo2/backends/cuda_src/cuda_functions.cu roshambo2/backends/cuda_src/cuda_interface.cpp)

# Link the OpenMP library for both
target_link_libraries(_roshambo2_cpp PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(_roshambo2_cuda PRIVATE OpenMP::OpenMP_CXX)


install(TARGETS _roshambo2_cpp
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX})


install(TARGETS _roshambo2_cuda
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX})