cmake_minimum_required(VERSION 3.15)
project(mm_strip_reconstruction)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ROOT detection
#set(ROOT_DIR "$ENV{ROOTSYS}/cmake/ROOTConfig.cmake")
message(STATUS "Looking for ROOT... ROOT_DIR='${ROOT_DIR}'")
find_package(ROOT REQUIRED COMPONENTS Tree RIO Core)
message(STATUS "Found ROOT version: ${ROOT_VERSION}")
include(${ROOT_USE_FILE})

# Add subdirectories
add_subdirectory(common)
add_subdirectory(decoder)
add_subdirectory(waveform_analysis)
add_subdirectory(clusterizer)
