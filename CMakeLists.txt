cmake_minimum_required(VERSION 3.19)

# Project name
project(scaleICP)

# LIBRARIES
include_directories( $ENV{EIGEN3_INCLUDE_DIR} )
add_library(lib STATIC "lib/mesh.h" "lib/mesh.cpp" "lib/icp.h" "lib/icp.cpp")

# Define an executable with a name and the files
add_executable(scaleICP "main.cpp")

target_include_directories(scaleICP PUBLIC "lib")  # Make the lib folder an include directory
target_link_libraries(scaleICP PUBLIC lib)  # Link this library to our executable