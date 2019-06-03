cmake_minimum_required(VERSION 2.8.9)

project(geometry-central)

### Policy settings
cmake_policy(SET CMP0054 NEW)   # don't implicitly dereference inside if()


### Process settings
option(BUILD_SHARED_LIBS "Build the shared library" FALSE)
if(BUILD_SHARED_LIBS)
    message("-- Building SHARED libraries")
else()
    message("-- Building STATIC libraries")
endif()


### Resolve optional dependencies

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake") # look for stuff in the /cmake directory
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/deps/eigen/cmake") # look for stuff in the Eigen/cmake directory
include(UpdateCacheVariable)

# Work with non-standard homebrew installations
# (from ceres build system)
if (CMAKE_SYSTEM_NAME MATCHES "Darwin")
  find_program(HOMEBREW_EXECUTABLE brew)
  mark_as_advanced(FORCE HOMEBREW_EXECUTABLE)
  if (HOMEBREW_EXECUTABLE)
    # Detected a Homebrew install, query for its install prefix.
    execute_process(COMMAND ${HOMEBREW_EXECUTABLE} --prefix
      OUTPUT_VARIABLE HOMEBREW_INSTALL_PREFIX
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    message(STATUS "Detected Homebrew with install prefix: "
      "${HOMEBREW_INSTALL_PREFIX}, adding to CMake search paths.")
    list(APPEND CMAKE_PREFIX_PATH "${HOMEBREW_INSTALL_PREFIX}")
  endif()
endif()


# Look for suitesparse (adapted from ceres build system)
option(SUITESPARSE "Enable SuiteSparse." ON)
if (SUITESPARSE)

  # Check for SuiteSparse and dependencies.
  find_package(SuiteSparse)
  if (SUITESPARSE_FOUND)
    # By default, if all of SuiteSparse's dependencies are found, GC is
    # built with SuiteSparse support.
    message("-- Found SuiteSparse ${SUITESPARSE_VERSION}, "
            "building with SuiteSparse.")
  else (SUITESPARSE_FOUND)
    # Disable use of SuiteSparse if it cannot be found and continue.
    message("-- Did not find all SuiteSparse dependencies, disabling "
      "SuiteSparse support.")
    update_cache_variable(SUITESPARSE OFF)
  endif (SUITESPARSE_FOUND)
else (SUITESPARSE)
  message("-- Building without SuiteSparse.")
endif (SUITESPARSE)

# HACK BY NICK
# The SuiteSparse logic above doesn't look for UMFpack, but we need it.
# This code attempts to find if by assuming that it will be in the same place
# as cholmod
if(SUITESPARSE_FOUND)
    string(REGEX REPLACE "cholmod" "umfpack" UMFPACK_LIBRARY ${CHOLMOD_LIBRARY})
    message("-- Guesstimated umfpack location as: ${UMFPACK_LIBRARY}")
    if(EXISTS ${UMFPACK_LIBRARY})
        list(APPEND SUITESPARSE_LIBRARIES ${UMFPACK_LIBRARY}) 
    else()
        message(WARNING "UMFPack guess failed, so we don't actually have SUITESPARSE support.")
        set(SUITESPARSE_FOUND FALSE)
    endif()
endif()

if(SUITESPARSE_FOUND)
    SET(HAVE_SUITESPARSE TRUE)
    SET(HAVE_SUITESPARSE TRUE PARENT_SCOPE)
    add_definitions(-DHAVE_SUITESPARSE)
    include_directories(${SUITESPARSE_INCLUDE_DIRS})
endif()

### Handle windows-specific fixes
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  add_definitions (-DNOMINMAX)  # don't use weird windows built-in min/max
  add_definitions (-D_USE_MATH_DEFINES)  # match unix behavior of constants in cmath
endif()

### Do anything needed for dependencies and bring their stuff in to scope
add_subdirectory(deps)

### Includes from dependencies
list(APPEND GEOMETRY_CENTRAL_INCLUDE_DIRS "deps/eigen")
list(APPEND GEOMETRY_CENTRAL_INCLUDE_DIRS "deps/RectangleBinPack/include")
list(APPEND GEOMETRY_CENTRAL_INCLUDE_DIRS "deps/nanort/include")
list(APPEND GEOMETRY_CENTRAL_INCLUDE_DIRS "deps/nanoflann/include")
list(APPEND GEOMETRY_CENTRAL_INCLUDE_DIRS "deps/glm")
list(APPEND GEOMETRY_CENTRAL_INCLUDE_DIRS "deps/tinyply/include")
list(APPEND GEOMETRY_CENTRAL_INCLUDE_DIRS "deps/happly")
if(SUITESPARSE_FOUND)
  list(APPEND GEOMETRY_CENTRAL_INCLUDE_DIRS ${SUITESPARSE_INCLUDE_DIRS})
endif()

### Includes from this project
list(APPEND GEOMETRY_CENTRAL_INCLUDE_DIRS "include")

include_directories(${GEOMETRY_CENTRAL_INCLUDE_DIRS})

### Compiler options
set( CMAKE_EXPORT_COMPILE_COMMANDS 1 ) # Emit a compile flags file to support completion engines (YouCompleteMe in my case)

### Compiler options
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using Clang (linux or apple) or GCC
  message("Using clang/gcc compiler flags")
  SET(BASE_CXX_FLAGS "-std=c++11 -Wall -Wextra -Werror -g3")
  SET(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-int-in-bool-context -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations")
  SET(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")

  # clang-specific things
  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=5 -fcolor-diagnostics")
  endif()

  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TRACE_INCLUDES}") # uncomment if you need to track down where something is getting included from
  #SET(CMAKE_CXX_FLAGS_DEBUG          "-g3")
  SET(CMAKE_CXX_FLAGS_DEBUG          "-g3 -fsanitize=address")
  SET(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
  SET(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
  SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
  message("Using Visual Studio compiler flags")
  set(BASE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  set(BASE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP") # parallel build
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


### Recurse to the source core
add_subdirectory(src)


### Export variables up the stack

# Convert paths to absolute
foreach(rel_path ${GEOMETRY_CENTRAL_INCLUDE_DIRS})
  get_filename_component(abs_path "${rel_path}" REALPATH BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
  list(APPEND GEOMETRY_CENTRAL_ABS_INCLUDE_DIRS ${abs_path})
endforeach(rel_path)

set(GEOMETRY_CENTRAL_INCLUDE_DIRS ${GEOMETRY_CENTRAL_ABS_INCLUDE_DIRS} PARENT_SCOPE)
