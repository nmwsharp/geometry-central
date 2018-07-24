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

target_include_directories(geometry-central PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../include")
target_include_directories(geometry-central PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../deps/eigen")
target_include_directories(geometry-central PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../deps/RectangleBinPack/include")
target_include_directories(geometry-central PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../deps/nanort/include")
target_include_directories(geometry-central PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../deps/nanoflann/include")
target_include_directories(geometry-central PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../deps/glm")
target_include_directories(geometry-central PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../deps/tinyply/include")
target_include_directories(geometry-central PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../deps/happly")
if(SUITESPARSE_FOUND)
  target_include_directories(geometry-central PUBLIC "${SUITESPARSE_INCLUDE_DIRS}")
endif()

# If we (optionally) resolved suitesparse support, we want to link against that as well
if(HAVE_SUITESPARSE)
    target_link_libraries(geometry-central ${SUITESPARSE_LIBRARIES})
endif(HAVE_SUITESPARSE)
