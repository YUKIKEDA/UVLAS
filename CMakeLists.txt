cmake_minimum_required(VERSION 3.25)

project(UVLAS VERSION 0.1.0 LANGUAGES CXX)

set(UVLAS_VERSION "${PROJECT_VERSION}-dev")
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include(cmake/CompilerPolicy.cmake)
include(cmake/ProjectOptions.cmake)
include(cmake/Dependencies.cmake)

# The platform key under which numerical results are qualified (QA-07).
string(TOLOWER "${CMAKE_SYSTEM_NAME}" UVLAS_OS_KEY)
string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" UVLAS_ARCH_KEY)
if(UVLAS_ARCH_KEY MATCHES "^(amd64|x86_64)$")
    set(UVLAS_ARCH_KEY "x64")
endif()
string(TOLOWER "${CMAKE_CXX_COMPILER_ID}" _uvlas_compiler)
string(REGEX MATCH "^[0-9]+\\.[0-9]+" _uvlas_compiler_version "${CMAKE_CXX_COMPILER_VERSION}")
set(UVLAS_COMPILER_KEY "${_uvlas_compiler}${_uvlas_compiler_version}")
configure_file(cmake/build_config.hpp.in generated/uvlas/build_config.hpp @ONLY)

# The commit of the checkout, refreshed at every build (QA-06).
add_custom_target(uvlas_git_revision ALL
    COMMAND "${CMAKE_COMMAND}" "-DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
        "-DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/generated/uvlas/git_revision.hpp"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/GitRevision.cmake"
    BYPRODUCTS "${CMAKE_CURRENT_BINARY_DIR}/generated/uvlas/git_revision.hpp"
    VERBATIM)

add_library(uvlas_core STATIC src/core/build_info.cpp)
uvlas_configure_target(uvlas_core)
target_include_directories(uvlas_core
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include"
    PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/generated")
target_compile_definitions(uvlas_core PRIVATE UVLAS_BUILD_CONFIGURATION="$<CONFIG>")
# build_info reports the versions of the numerical dependencies.
target_link_libraries(uvlas_core PRIVATE
    uvlas::eigen Spectra::Spectra Boost::headers nlohmann_json::nlohmann_json)
add_dependencies(uvlas_core uvlas_git_revision)

# QA-03: the first baseline case. Built with the product flags, compared by CTest.
add_executable(uvlas_numeric_canary qa/numeric_canary/main.cpp)
uvlas_configure_target(uvlas_numeric_canary)
target_link_libraries(uvlas_numeric_canary PRIVATE
    uvlas_core uvlas::eigen Boost::multiprecision)

add_executable(uvlas_cli src/app/main.cpp)
uvlas_configure_target(uvlas_cli)
set_target_properties(uvlas_cli PROPERTIES OUTPUT_NAME uvlas)
target_link_libraries(uvlas_cli PRIVATE uvlas_core)

include(CTest)
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()

