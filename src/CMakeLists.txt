cmake_minimum_required(VERSION 2.8.9)

# Maybe stop from CMAKEing in the wrong place
if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "Source and build directories cannot be the same. Go use the /build directory.")
endif()

# Collect all .cpp files in the subdirectories
file(GLOB_RECURSE SRCS *.cpp)
# file(GLOB_RECURSE HEADS *.h *.hpp *.ipp)

# Create a single library for the project
add_library(
        geometry-central 
        ${SRCS}
        )

# If we (optionally) resolved SPQR support, we want to link against that as well
if(HAVE_SUITESPARSE)
    target_link_libraries(geometry-central ${SUITESPARSE_LIBRARIES})
endif(HAVE_SUITESPARSE)
