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

### Get software version
find_package(Git)
include(GetProjectVersionFromGit)
get_project_version_from_git("${PROJECT_NAME}" "${CMAKE_CURRENT_SOURCE_DIR}" "${GIT_EXECUTABLE}")

### Handle windows-specific fixes
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  add_definitions (-DNOMINMAX)  # don't use weird windows built-in min/max
  add_definitions (-D_USE_MATH_DEFINES)  # match unix behavior of constants in cmath
endif()

# Set up our directory structure for output libraries and binaries
include(GNUInstallDirs)

### Do anything needed for dependencies and bring their stuff in to scope
add_subdirectory(deps)

### Recurse to the source code
add_subdirectory(src)

# install
install(
  TARGETS geometry-central
  EXPORT GeometryCentral
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(
  DIRECTORY ${CMAKE_SOURCE_DIR}/include/
  DESTINATION include
  FILES_MATCHING
  PATTERN "*.h"
  PATTERN "*.ipp")

install(EXPORT GeometryCentral
  FILE
    geometrycentral.cmake
  NAMESPACE
    geometrycentral::
  DESTINATION
  ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/geometrycentral-config-install.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/geometrycentral-config.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/geometrycentral-config-version.cmake"
    VERSION ${${PROJECT_NAME}_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/geometrycentral-config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/geometrycentral-config-version.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)
