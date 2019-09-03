cmake_minimum_required(VERSION 3.8.0)

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
if(SUITESPARSE AND SUITESPARSE_FOUND)
    string(REGEX REPLACE "cholmod" "umfpack" UMFPACK_LIBRARY ${CHOLMOD_LIBRARY})
    message("-- Guesstimated umfpack location as: ${UMFPACK_LIBRARY}")
    if(EXISTS ${UMFPACK_LIBRARY})
        list(APPEND SUITESPARSE_LIBRARIES ${UMFPACK_LIBRARY}) 
    else()
        message(WARNING "UMFPack guess failed, so we don't actually have SUITESPARSE support.")
        set(SUITESPARSE_FOUND FALSE)
    endif()
endif()

if(SUITESPARSE AND SUITESPARSE_FOUND)
    SET(GC_HAVE_SUITESPARSE TRUE)
    SET(GC_HAVE_SUITESPARSE TRUE PARENT_SCOPE)
    add_definitions(-DGC_HAVE_SUITESPARSE)
    include_directories(${SUITESPARSE_INCLUDE_DIRS})
else()
    SET(GC_HAVE_SUITESPARSE FALSE)
    SET(GC_HAVE_SUITESPARSE FALSE PARENT_SCOPE)
endif()

### Handle windows-specific fixes
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  add_definitions (-DNOMINMAX)  # don't use weird windows built-in min/max
  add_definitions (-D_USE_MATH_DEFINES)  # match unix behavior of constants in cmath
endif()

### Do anything needed for dependencies and bring their stuff in to scope
add_subdirectory(deps)

### Recurse to the source code
add_subdirectory(src)
