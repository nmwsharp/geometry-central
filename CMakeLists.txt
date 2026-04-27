cmake_minimum_required(VERSION 3.14.0)

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


list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake") # look for stuff in the /cmake directory
include(UpdateCacheVariable)
# Install a standard CMake config package so downstream users do not need a custom FindGeometryCentral.cmake module.
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(GC_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/GeometryCentral")

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


### Handle windows-specific fixes
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  add_definitions (-DNOMINMAX)  # don't use weird windows built-in min/max
  add_definitions (-D_USE_MATH_DEFINES)  # match unix behavior of constants in cmath
endif()

### Do anything needed for dependencies and bring their stuff in to scope
add_subdirectory(deps)

# copy variables set by deps upward
SET(GC_HAVE_SUITESPARSE ${GC_HAVE_SUITESPARSE} PARENT_SCOPE)

if(GC_EIGEN_WAS_DOWNLOADED)
  message(WARNING "Eigen was auto-downloaded into the build tree. Consumers of the installed package will "
                  "need Eigen 3.3+ available as a system package (e.g. installed via their package manager).")
endif()

### Recurse to the source code
add_subdirectory(src)

# Install the library and export it as an imported target for find_package().
install(
  TARGETS geometry-central
  EXPORT GeometryCentralTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING
  PATTERN "*.h"
  PATTERN "*.ipp")

# These vendored headers are referenced by installed public headers, so ship them
# in the install tree rather than relying on downstream postInstall copy steps.
install(
  FILES "${CMAKE_CURRENT_SOURCE_DIR}/deps/happly/happly.h"
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/GeometryCentralConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/GeometryCentralConfig.cmake"
  INSTALL_DESTINATION "${GC_INSTALL_CMAKEDIR}")

# Install the exported target definitions and the package entry point side-by-side.
install(
  EXPORT GeometryCentralTargets
  NAMESPACE geometry-central::
  DESTINATION "${GC_INSTALL_CMAKEDIR}")

install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/GeometryCentralConfig.cmake"
  DESTINATION "${GC_INSTALL_CMAKEDIR}")
