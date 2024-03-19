# This function attempts to determine the project version from git tags.
# It sets variables with a given prefix to reflect the version information.
# - version_prefix: The prefix for version-related variables.
# - project_source_dir: The directory containing the project's source code.
# - git_exe: The path to the git executable.
#
# On success, it sets variables indicating the version (major, minor, patch, etc.).
# On failure, a status message is printed.

function(get_project_version_from_git version_prefix project_source_dir git_exe)
  if(EXISTS "${git_exe}" AND EXISTS "${project_source_dir}/.git")
    execute_process(
      COMMAND "${git_exe}" describe --tags --always
      WORKING_DIRECTORY "${project_source_dir}"
      RESULT_VARIABLE git_result
      OUTPUT_VARIABLE version_output
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)

    if(git_result EQUAL 0)
      string(REGEX MATCH "([0-9]+)\\.([0-9]+)\\.([0-9]+)" version_parts "${version_output}")
      if(version_parts)
        # Extract version numbers
        set(version_major "${CMAKE_MATCH_1}")
        set(version_minor "${CMAKE_MATCH_2}")
        set(version_patch "${CMAKE_MATCH_3}")
        set(full_version "${version_major}.${version_minor}.${version_patch}")

        # Set variables in the parent scope
        set("${version_prefix}_VERSION" "${full_version}" PARENT_SCOPE)
        set("${version_prefix}_MAJOR" "${version_major}" PARENT_SCOPE)
        set("${version_prefix}_MINOR" "${version_minor}" PARENT_SCOPE)
        set("${version_prefix}_PATCH" "${version_patch}" PARENT_SCOPE)
      endif()
    else()
      message(STATUS "Unable to determine project version from git tags.")
    endif()
  else()
    message(STATUS "Git not available or .git directory not found.")
  endif()
endfunction()
