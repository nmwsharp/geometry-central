cmake_minimum_required(VERSION 3.8.0)

# Set up testing
enable_testing()

project(geometry-central-tests)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

### Configure the compiler

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

  # using Clang (linux or apple) or GCC
  message("Using clang/gcc compiler flags")
  SET(BASE_CXX_FLAGS "-std=c++11 -Wall -Wextra -g3")
  SET(CMAKE_CXX_FLAGS_DEBUG          "")
  SET(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces -Wno-unused-private-field")
  SET(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")

  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message("Setting clang-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=3 -fcolor-diagnostics")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-fsanitize=address -fno-limit-debug-info")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message("Setting gcc-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wno-maybe-uninitialized -Wno-format-zero-length -Wno-unused-but-set-parameter -Wno-unused-but-set-variable")
  endif()

  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  SET(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -g3")
  SET(CMAKE_CXX_FLAGS_RELEASE        "-O3 -march=native -DNDEBUG")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
  message("Using Visual Studio compiler flags")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /W4")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /MP") # parallel build
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4267\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4244\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4305\"")  # ignore truncation on initialization
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")

  add_definitions(/D "_CRT_SECURE_NO_WARNINGS")
  add_definitions (-DNOMINMAX)
else()
  # unrecognized
  message( FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]" )
endif()


### Download and unpack googletest at configure time
# (from the googletest docs)
if(${CMAKE_VERSION} VERSION_LESS 3.14)
  include(FetchContentLocal)
else()
  include(FetchContent)
endif()

FetchContent_Declare(
    googletest
    GIT_REPOSITORY    https://github.com/google/googletest.git
    GIT_TAG           release-1.10.0
    SOURCE_DIR        "${CMAKE_CURRENT_BINARY_DIR}/googletest-src"
    BINARY_DIR        "${CMAKE_CURRENT_BINARY_DIR}/googletest-build"
)

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

set_property(TARGET gtest PROPERTY CXX_STANDARD 11)

# silence warning about deprecated copies in gtest
# https://github.com/abseil/abseil-cpp/issues/948#issuecomment-846592426
if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  target_compile_options(gtest INTERFACE -Wno-deprecated-copy)
endif()


### Set (hardcode) a variable for where to find the assets
set(REL_ASSETS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/assets/")
get_filename_component(ABS_ASSETS_DIR "${REL_ASSETS_DIR}" ABSOLUTE)
add_definitions(-DGC_TEST_ASSETS_ABS_PATH=\"${ABS_ASSETS_DIR}\")


# Build the tests
set(TEST_SRCS
  src/main_test.cpp
  src/load_test_meshes.cpp
  src/eigen_interop_helpers_test.cpp
  src/halfedge_mesh_test.cpp
  src/halfedge_mutation_test.cpp
  src/halfedge_geometry_test.cpp
  src/surface_misc_test.cpp
  src/point_cloud_test.cpp
  src/linear_algebra_test.cpp
  src/stl_reader_test.cpp
  src/intrinsic_triangulation_test.cpp
)

add_executable(geometry-central-test "${TEST_SRCS}")
target_include_directories(geometry-central-test PUBLIC "${gtest_SOURCE_DIR}/include")
target_include_directories(geometry-central-test PRIVATE "include/")
target_link_libraries(geometry-central-test gtest_main geometry-central)

# Add geometry central as a subproject
add_subdirectory(../ "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
if(GC_HAVE_SUITESPARSE)
    add_definitions(-DGC_HAVE_SUITESPARSE)
endif()
