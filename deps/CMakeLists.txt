cmake_minimum_required(VERSION 2.8.9)

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

message(INFO "  resolving Eigen dependency:")

set(EIGEN_LOCATION "" CACHE STRING "Path to Eigen installation. Only needed if Eigen is not already available in CMAKE, will not be found by find_package(), and should not be auto-downloaded.")

# Try (1) loading from specified path
if(GC_EIGEN_LOCATION)
  message(INFO "    variable EIGEN_LOCATION was set with value ${GC_EIGEN_LOCATION}")
 
  # verify that the directory seems to acrtually include Eigen
  # NOTE should also hold Eigen with suitable version but we don't check that
  set(GC_EIGEN_FILE_CHECK "${GC_EIGEN_LOCATION}/Eigen/Core")
  if(NOT EXISTS ${GC_EIGEN_FILE_CHECK})
    message(FATAL_ERROR "Path from GC_EIGEN_LOCATION does not seem to point to Eigen. Could not find file ${GC_EIGEN_FILE_CHECK}")
  endif()

  set(EIGEN3_INCLUDE_DIR ${GC_EIGEN_LOCATION})
  set(EIGEN3_VERSION "unknown") # would be great to read this. feels like we should be able to put find_package() to work
  
  message(INFO "    success. Found suitable Eigen from EIGEN_LOCATION at ${EIGEN3_INCLUDE_DIR}")

# Try (2) hoping Eigen is already resolved
elseif((${EIGEN3_INCLUDE_DIR}) AND (NOT ${EIGEN3_VERSION} VERSION_LESS 3.3))
  message(INFO "    success. CMake environment already has suitable Eigen defined at ${EIGEN3_INCLUDE_DIR}")

# Try (3 and 4) find_package() or download
else()

  # Try (3) find_package() Eigen
  find_package(Eigen3 3.3)

  # The find_package() succeeded!
  if(EIGEN3_FOUND)
    message(INFO "    success. find_package() succeeded with Eigen at ${EIGEN3_INCLUDE_DIR}")

  # The find_package() failed, manually download the package (4)
  else()
    
    message(INFO "    Could not find Eigen anywhere. Downloading it (or using cached download).")

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
    message(INFO "    success. Found/Downloaded Eigen in ${EIGEN3_INCLUDE_DIR}")

  endif()

endif()

# Propagate Eigen variables to parent scope
SET(EIGEN3_INCLUDE_DIR ${EIGEN3_INCLUDE_DIR} PARENT_SCOPE)
SET(EIGEN3_VERSION ${EIGEN3_VERSION} PARENT_SCOPE)

# ===  done finding Eigen


# TODO fill out a single GC_DEP_INCLUDES string 
