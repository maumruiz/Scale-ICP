cmake_minimum_required(VERSION 3.19)

# Arbitrary name
project(scaleICP)

# LIBRARIES
SET( EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}" )
LIST(APPEND LIB_FILES "lib/mesh.h" "lib/mesh.cpp")
add_library(lib STATIC ${LIB_FILES})

# EXECUTABLES
# LIST(APPEND HEADER_FILES "mesh.h")
LIST(APPEND SOURCE_FILES "main.cpp")

# Define an executable with a name and the files
add_executable(scaleICP ${SOURCE_FILES} ${HEADER_FILES})

target_include_directories(scaleICP PUBLIC "lib")  # Make the lib folder an include directory
target_include_directories(scaleICP PUBLIC ${EIGEN3_INCLUDE_DIR})
target_link_libraries(scaleICP PUBLIC lib)  # Link this library to our executable

# Redirect standard output
execute_process(
    COMMAND ./scaleICP ${SOURCE_FILES} ${HEADER_FILES}
)