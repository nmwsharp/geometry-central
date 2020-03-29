cmake_minimum_required(VERSION 3.8)

# We will accumulate all depenency includes and link libraries here
# (can be relative paths when building list, we will convert all to absolute at end)
SET(GC_DEP_INCLUDES "")
SET(GC_DEP_LIBRARIES "")

# Make sure we don't try to compile anything in the dependencies with Werror on
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error")
endif()


# Directory in which to place downloaded libraries
set(GC_ARCHIVE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/downloads/)

# === Resolve Eigen dependency

# This script uses the following policies in order to find Eigen:
#  (1) If GC_EIGEN_LOCATION is set, use the copy of Eigen there (an error will be thrown if Eigen is not found there)
#  (2) If the usual Eigen variables EIGEN3_INCLUDE_DIR and EIGEN3_VERSION are set, CMAKE already has Eigen resolved, so use that copy
#  (3) Invoke find_package() (if no suitable Eigen is found, will proceed to next step)
#  (4) Download Eigen in to the source tree

message("--  resolving Eigen dependency:")

set(EIGEN_LOCATION "" CACHE STRING "Path to Eigen installation. Only needed if Eigen is not already available in CMAKE, will not be found by find_package(), and should not be auto-downloaded.")

# Try (1) loading from specified path
if(GC_EIGEN_LOCATION)
  message("--    variable EIGEN_LOCATION was set with value ${GC_EIGEN_LOCATION}")

  # verify that the directory seems to acrtually include Eigen
  # NOTE should also hold Eigen with suitable version but we don't check that
  set(GC_EIGEN_FILE_CHECK "${GC_EIGEN_LOCATION}/Eigen/Core")
  if(NOT EXISTS ${GC_EIGEN_FILE_CHECK})
    message(FATAL_ERROR "Path from GC_EIGEN_LOCATION does not seem to point to Eigen. Could not find file ${GC_EIGEN_FILE_CHECK}")
  endif()

  set(EIGEN3_INCLUDE_DIR ${GC_EIGEN_LOCATION})
  set(EIGEN3_VERSION "unknown") # would be great to read this. feels like we should be able to put find_package() to work

  message("--    success. Found suitable Eigen from EIGEN_LOCATION at ${EIGEN3_INCLUDE_DIR}")

# Try (2) hoping Eigen is already resolved
elseif((${EIGEN3_INCLUDE_DIR}) AND (NOT ${EIGEN3_VERSION} VERSION_LESS 3.3))
  message("--    success. CMake environment already has suitable Eigen defined at ${EIGEN3_INCLUDE_DIR}")

# Try (3 and 4) find_package() or download
else()

  # Try (3) find_package() Eigen
  find_package(Eigen3 3.3)

  # The find_package() succeeded!
  if(EIGEN3_FOUND)
    message("--    success. find_package() succeeded with Eigen at ${EIGEN3_INCLUDE_DIR}")

  # The find_package() failed, manually download the package (4)
  else()

    message("--    Could not find Eigen anywhere. Downloading it (or using cached download).")

    include(DownloadProject)
    download_project(
      PROJ EIGEN3
      PREFIX ${GC_ARCHIVE_PATH}
      GIT_REPOSITORY https://github.com/eigenteam/eigen-git-mirror.git
      GIT_TAG 3.3.7
		)

    set(EIGEN3_INCLUDE_DIR ${EIGEN3_SOURCE_DIR})

    # verify that the directory seems to actually include Eigen
    set(GC_EIGEN_FILE_CHECK "${EIGEN3_INCLUDE_DIR}/Eigen/Core")
    if(NOT EXISTS ${GC_EIGEN_FILE_CHECK})
      message(FATAL_ERROR "Path from GC_EIGEN_LOCATION does not seem to point to Eigen. Could not find file ${GC_EIGEN_FILE_CHECK}")
    endif()
    message("--    success. Found/Downloaded Eigen in ${EIGEN3_INCLUDE_DIR}")

  endif()

endif()

# Propagate Eigen variables to parent scope in case we every want to use them for something else
SET(EIGEN3_INCLUDE_DIR ${EIGEN3_INCLUDE_DIR} PARENT_SCOPE)
SET(EIGEN3_VERSION ${EIGEN3_VERSION} PARENT_SCOPE)

# Add Eigen to the list of deps to include
list(APPEND GC_DEP_INCLUDES "${EIGEN3_INCLUDE_DIR}")

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
    add_definitions(-DGC_HAVE_SUITESPARSE)
    include_directories(${SUITESPARSE_INCLUDE_DIRS})
else()
    SET(GC_HAVE_SUITESPARSE FALSE)
    SET(GC_HAVE_SUITESPARSE FALSE PARENT_SCOPE)
endif()

if(GC_HAVE_SUITESPARSE)
  list(APPEND GC_DEP_INCLUDES "${SUITESPARSE_INCLUDE_DIRS}")
  list(APPEND GC_DEP_LIBRARIES "${SUITESPARSE_LIBRARIES}")
endif()

# ===  done finding Suitesparse


# Find other simpler dependencies
list(APPEND GC_DEP_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/nanort/include")
list(APPEND GC_DEP_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/nanoflann/include")
list(APPEND GC_DEP_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/happly")


# == Ensure all paths are absolute

# Helper function
macro(convert_paths_to_absolute files)
  foreach(path ${${files}})
    get_filename_component(abs_path "${path}" ABSOLUTE)
    list(REMOVE_ITEM ${files} "${path}")
    list(APPEND ${files} "${abs_path}")
  endforeach()
endmacro()

convert_paths_to_absolute(GC_DEP_INCLUDES)
convert_paths_to_absolute(GC_DEP_LIBRARIES)

# Propagate dep variables up to parten
SET(GC_DEP_INCLUDES ${GC_DEP_INCLUDES} PARENT_SCOPE)
SET(GC_DEP_LIBRARIES ${GC_DEP_LIBRARIES} PARENT_SCOPE)
