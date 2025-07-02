cmake_minimum_required(VERSION 3.14.0)
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

project(geometry-central)
set(PROJECT_VERSION 0.1.0)

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


### Handle windows-specific fixes
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  add_definitions (-DNOMINMAX)  # don't use weird windows built-in min/max
  add_definitions (-D_USE_MATH_DEFINES)  # match unix behavior of constants in cmath
endif()

### Do anything needed for dependencies and bring their stuff in to scope
add_subdirectory(deps)

# copy variables set by deps upward
SET(GC_HAVE_SUITESPARSE ${GC_HAVE_SUITESPARSE} PARENT_SCOPE)

### Recurse to the source code
add_subdirectory(src)

# install

configure_package_config_file(geometry-central-config.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/geometry-central-config.cmake INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/geometry-central/cmake)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/geometry-centralConfigVersion.cmake COMPATIBILITY SameMajorVersion)
install(TARGETS geometry-central eigen nanort nanoflann happly EXPORT geometry-central-targets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(EXPORT geometry-central-targets FILE geometry-centralTargets.cmake DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/geometry-central/cmake)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/geometry-central-config.cmake DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/geometry-central/cmake)
install(
  DIRECTORY ${CMAKE_SOURCE_DIR}/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING
  PATTERN "*.h"
  PATTERN "*.ipp")
