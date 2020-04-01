

function(eigen3checker GC_EIGEN_LOCATION EIGEN3_FIND_VERSION)
  # Function looks for Eigen3 at GC_EIGEN_LOCATION
  # Eigen version much be >= EIGEN3_FIND_VERSION
  #
  # If Eigen is found it sets the following variables on PARENT_SCOPE
  #
  #  EIGEN3_FOUND       - system has eigen lib with correct version
  #  EIGEN3_INCLUDE_DIR - the eigen include directory
  #  EIGEN3_VERSION     - eigen version
  #
  # It defines the following target:
  #
  #  Eigen3::Eigen
  #

  if(NOT EIGEN3_INCLUDE_DIR)
    set(EIGEN3_FOUND false PARENT_SCOPE)

    # Search for the signature_of_eigen3_matrix_library
    find_path(EIGEN3_INCLUDE_DIR NAMES signature_of_eigen3_matrix_library
        PATHS
          ${GC_EIGEN_LOCATION}
        PATH_SUFFIXES
          eigen3
          eigen
      )

    # message(STATUS "Eigen Include Dir: ${EIGEN3_INCLUDE_DIR}")
    if(EXISTS "${EIGEN3_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h")

      # Parse version from Macros.h
      file(READ "${EIGEN3_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h" _eigen3_version_header)

      string(REGEX MATCH "define[ \t]+EIGEN_WORLD_VERSION[ \t]+([0-9]+)" _eigen3_world_version_match "${_eigen3_version_header}")
      set(EIGEN3_WORLD_VERSION "${CMAKE_MATCH_1}")
      string(REGEX MATCH "define[ \t]+EIGEN_MAJOR_VERSION[ \t]+([0-9]+)" _eigen3_major_version_match "${_eigen3_version_header}")
      set(EIGEN3_MAJOR_VERSION "${CMAKE_MATCH_1}")
      string(REGEX MATCH "define[ \t]+EIGEN_MINOR_VERSION[ \t]+([0-9]+)" _eigen3_minor_version_match "${_eigen3_version_header}")
      set(EIGEN3_MINOR_VERSION "${CMAKE_MATCH_1}")

      set(EIGEN3_VERSION ${EIGEN3_WORLD_VERSION}.${EIGEN3_MAJOR_VERSION}.${EIGEN3_MINOR_VERSION})

      # message(STATUS "Eigen Version: ${EIGEN3_VERSION}")
      if(${EIGEN3_VERSION} VERSION_LESS ${EIGEN3_FIND_VERSION})
        # Complain about the inadequate version
        message(STATUS "Eigen3 version ${EIGEN3_VERSION} found in ${EIGEN3_INCLUDE_DIR}, " "but at least version ${EIGEN3_FIND_VERSION} is required")
      else(${EIGEN3_VERSION} VERSION_LESS ${EIGEN3_FIND_VERSION})
        # Suitable version of Eigen3
        set(EIGEN3_VERSION ${EIGEN3_VERSION} PARENT_SCOPE)
        set(EIGEN3_FOUND true PARENT_SCOPE)

        # In theory Eigen3::Eigen should not have been defined... but
        # it may have been from some other project...
        if(NOT TARGET Eigen3::Eigen)
          add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
          set_target_properties(Eigen3::Eigen PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
        endif()
      endif(${EIGEN3_VERSION} VERSION_LESS ${EIGEN3_FIND_VERSION})
    endif()
  else(NOT EIGEN3_INCLUDE_DIR)
    # EIGEN3_INCLUDE_DIR is already set
    message(WARNING "EIGEN3_INCLUDE_DIR is already set to ${EIGEN3_INCLUDE_DIR} ignoring GC_EIGEN_LOCATION.")
  endif(NOT EIGEN3_INCLUDE_DIR)
endfunction()