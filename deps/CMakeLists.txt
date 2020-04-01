cmake_minimum_required(VERSION 3.8)

# We will accumulate all depenency includes and link libraries here
# (can be relative paths when building list, we will convert all to absolute at end)
SET(GC_DEP_LIBS "")

# Make sure we don't try to compile anything in the dependencies with Werror on
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error")
endif()

# === Resolve Eigen dependency

# This script uses the following policies in order to find Eigen:
#  (1) If GC_EIGEN_LOCATION is set, use the copy of Eigen there (an error will be thrown if Eigen is not found there)
#  (2) Invoke find_package() (if no suitable Eigen is found, will proceed to next step)
#  (3) Download Eigen in to the source tree

message("--  resolving Eigen dependency:")

set(GC_EIGEN_LOCATION "" CACHE STRING "Path to Eigen installation. Only needed if Eigen is not already available in CMAKE, will not be found by find_package(), and should not be auto-downloaded.")

# Try (1) loading from specified path
if(GC_EIGEN_LOCATION)
  message("--    variable GC_EIGEN_LOCATION was set with value ${GC_EIGEN_LOCATION}")

  include(EigenChecker)
  eigen3checker(${GC_EIGEN_LOCATION} 3.3)

  if(EIGEN3_FOUND AND TARGET Eigen3::Eigen)
    message("--    success. Found suitable Eigen from GC_EIGEN_LOCATION at ${EIGEN3_INCLUDE_DIR}")
  else()
    message(FATAL_ERROR "Eigen could not be found at GC_EIGEN_LOCATION = ${GC_EIGEN_LOCATION}.")
  endif()

# Try (2 and 3) find_package() or download
else()

  # Try (2) find_package() Eigen
  find_package (Eigen3 3.3 QUIET)

  if(EIGEN3_FOUND)
    message("--    success. find_package() succeeded with Eigen at ${EIGEN3_INCLUDE_DIR}")

  # The find_package() failed, manually download the package (3)
  else()

    message("--    Could not find Eigen anywhere. Downloading it (or using cached download).")

    if(${CMAKE_VERSION} VERSION_LESS 3.11)
      include(FetchContentLocal)
    else()
      include(FetchContent)
    endif()
    FetchContent_Declare(
        eigen
        GIT_REPOSITORY  https://github.com/eigenteam/eigen-git-mirror.git
        GIT_TAG         3.3.7
        SOURCE_DIR      "${CMAKE_CURRENT_BINARY_DIR}/eigen-src"
        BINARY_DIR      "${CMAKE_CURRENT_BINARY_DIR}/eigen-build"
    )
    FetchContent_GetProperties(eigen)
    if(NOT eigen_POPULATED)
        FetchContent_Populate(eigen)
        # Create an interface library target
        add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
        target_include_directories(Eigen3::Eigen INTERFACE $<BUILD_INTERFACE:${eigen_SOURCE_DIR}>)
    endif()
    message("--    success. Found/Downloaded Eigen in ${eigen_SOURCE_DIR}")
  endif()
endif()

if(TARGET Eigen3::Eigen)
  set_target_properties(Eigen3::Eigen PROPERTIES IMPORTED_GLOBAL TRUE)
  list(APPEND GC_DEP_LIBS Eigen3::Eigen)
endif()

# ===  done finding Eigen


# ===  Try to find suitesparse

# Look for suitesparse (adapted from ceres build system)
option(SUITESPARSE "Enable SuiteSparse." ON)
if (SUITESPARSE)
  message("-- Searching for SuiteSparse.")

  # Check for SuiteSparse and dependencies.
  find_package(SuiteSparse)
  if (SUITESPARSE_FOUND)
    # By default, if all of SuiteSparse's dependencies are found, GC is
    # built with SuiteSparse support.
    message("-- Found SuiteSparse ${SUITESPARSE_VERSION}, building with SuiteSparse.")
  else (SUITESPARSE_FOUND)
    # Disable use of SuiteSparse if it cannot be found and continue.
    message("-- Did not find all SuiteSparse dependencies, disabling SuiteSparse support.")
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
    set(SUITESPARSE_INCLUDE_DIRS "${SUITESPARSE_INCLUDE_DIRS}" PARENT_SCOPE)
    set(SUITESPARSE_LIBRARIES    "${SUITESPARSE_LIBRARIES}"    PARENT_SCOPE)
else()
    SET(GC_HAVE_SUITESPARSE FALSE)
    SET(GC_HAVE_SUITESPARSE FALSE PARENT_SCOPE)
endif()

# ===  done finding Suitesparse

add_library(nanort INTERFACE)
target_include_directories(nanort INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/nanort/include>)

add_library(nanoflann INTERFACE)
target_include_directories(nanoflann INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/nanoflann/include>)

add_library(happly INTERFACE)
target_include_directories(happly INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/happly>)

# Find other simpler dependencies
list(APPEND GC_DEP_LIBS nanort)
list(APPEND GC_DEP_LIBS nanoflann)
list(APPEND GC_DEP_LIBS happly)

# == Ensure all paths are absolute

# # Helper function
# macro(convert_paths_to_absolute files)
#   foreach(path ${${files}})
#     get_filename_component(abs_path "${path}" ABSOLUTE)
#     list(REMOVE_ITEM ${files} "${path}")
#     list(APPEND ${files} "${abs_path}")
#   endforeach()
# endmacro()

# convert_paths_to_absolute(GC_DEP_LIBS)
# convert_paths_to_absolute(GC_DEP_LIBRARIES)

# Propagate dep variables up to parent
SET(GC_DEP_LIBS ${GC_DEP_LIBS} PARENT_SCOPE)
# SET(GC_DEP_LIBRARIES ${GC_DEP_LIBRARIES} PARENT_SCOPE)
