cmake_minimum_required(VERSION 3.25)

project(UVLAS VERSION 0.1.0 LANGUAGES CXX)

set(UVLAS_VERSION "${PROJECT_VERSION}-dev")
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include(cmake/CompilerPolicy.cmake)
include(cmake/ProjectOptions.cmake)
# Resolve CTest's default BUILD_TESTING=ON before selecting required notices.
include(CTest)
option(UVLAS_BUILD_BENCHMARKS "Build the separately linked kernel benchmarks" OFF)
option(UVLAS_BUILD_PYTHON "Build the low-level CPython product bindings" OFF)
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

set(UVLAS_CORE_SOURCES src/core/build_info.cpp src/core/error_certificate.cpp
    src/core/checked_dimensions.cpp src/core/state.cpp src/core/instrumentation.cpp)
add_library(uvlas_core STATIC ${UVLAS_CORE_SOURCES})
uvlas_configure_target(uvlas_core)
target_include_directories(uvlas_core
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include"
    PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/generated")
target_compile_definitions(uvlas_core PRIVATE UVLAS_BUILD_CONFIGURATION="$<CONFIG>")
# build_info reports the versions of the numerical dependencies.
target_link_libraries(uvlas_core PRIVATE
    uvlas::eigen Spectra::Spectra Boost::headers Boost::hash2 nlohmann_json::nlohmann_json)
add_dependencies(uvlas_core uvlas_git_revision)
if(WIN32)
    target_link_libraries(uvlas_core PRIVATE psapi)
endif()

set(UVLAS_NUMERICS_SOURCES src/numerics/storage.cpp src/numerics/precision.cpp
    src/numerics/eigen_common.cpp src/numerics/eigen_metrics.cpp src/numerics/eigen_action.cpp
    src/numerics/backend_memory.cpp src/numerics/linear_metrics.cpp src/numerics/linear_result.cpp src/numerics/linear_validation.cpp
    src/numerics/gmres_solve.cpp src/numerics/gmres_residual_detail.cpp
    src/numerics/eigen_arnoldi.cpp
    src/numerics/eigen_full.cpp src/numerics/dense_solve.cpp src/numerics/sparse_spd.cpp src/numerics/sparse_triangular.cpp src/numerics/sparse_solve.cpp src/numerics/eigen_generalized.cpp)
add_library(uvlas_numerics STATIC ${UVLAS_NUMERICS_SOURCES})
uvlas_configure_target(uvlas_numerics)
target_include_directories(uvlas_numerics PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(uvlas_numerics PUBLIC uvlas_core PRIVATE uvlas::eigen Spectra::Spectra Boost::headers Boost::hash2)

if(UVLAS_BUILD_PYTHON)
    if(NOT "python" IN_LIST VCPKG_MANIFEST_FEATURES)
        message(FATAL_ERROR "UVLAS_BUILD_PYTHON requires the vcpkg python feature")
    endif()
    set_target_properties(uvlas_core uvlas_numerics PROPERTIES POSITION_INDEPENDENT_CODE ON)
    add_subdirectory(src/python)
endif()

# QA-03: the first baseline case. Built with the product flags, compared by CTest.
add_executable(uvlas_numeric_canary qa/numeric_canary/main.cpp)
uvlas_configure_target(uvlas_numeric_canary)
target_link_libraries(uvlas_numeric_canary PRIVATE
    uvlas_core uvlas::eigen Boost::multiprecision)

add_executable(uvlas_eigen_arnoldi_canary qa/eigen_arnoldi/main.cpp)
uvlas_configure_target(uvlas_eigen_arnoldi_canary)
target_link_libraries(uvlas_eigen_arnoldi_canary PRIVATE uvlas_numerics)

add_executable(uvlas_eigen_common_canary qa/eigen_common/main.cpp)
uvlas_configure_target(uvlas_eigen_common_canary)
target_link_libraries(uvlas_eigen_common_canary PRIVATE uvlas_numerics)

add_executable(uvlas_eigen_full_canary qa/eigen_full/main.cpp)
uvlas_configure_target(uvlas_eigen_full_canary)
target_link_libraries(uvlas_eigen_full_canary PRIVATE uvlas_numerics)

add_executable(uvlas_cli src/app/main.cpp)
uvlas_configure_target(uvlas_cli)
set_target_properties(uvlas_cli PROPERTIES OUTPUT_NAME uvlas)
target_link_libraries(uvlas_cli PRIVATE uvlas_core)

add_executable(uvlas_eigen_generalized_canary qa/eigen_generalized/main.cpp)
uvlas_configure_target(uvlas_eigen_generalized_canary)
target_link_libraries(uvlas_eigen_generalized_canary PRIVATE uvlas_numerics)

if(BUILD_TESTING)
    include(cmake/AllocationAudit.cmake)
    # Audits recompile product sources; mixing normal Eigen inline definitions is forbidden.
    add_library(uvlas_core_audit STATIC ${UVLAS_CORE_SOURCES})
    uvlas_configure_target(uvlas_core_audit)
    uvlas_enable_allocation_audit(uvlas_core_audit)
    target_include_directories(uvlas_core_audit PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include"
        PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/generated")
    target_compile_definitions(uvlas_core_audit PRIVATE UVLAS_BUILD_CONFIGURATION="$<CONFIG>")
    target_link_libraries(uvlas_core_audit PRIVATE uvlas::eigen Spectra::Spectra
        Boost::headers Boost::hash2 nlohmann_json::nlohmann_json)
    if(WIN32)
        target_link_libraries(uvlas_core_audit PRIVATE psapi)
    endif()
    add_dependencies(uvlas_core_audit uvlas_git_revision)
    uvlas_require_audited_dependencies(uvlas_core_audit)
    add_library(uvlas_numerics_audit STATIC ${UVLAS_NUMERICS_SOURCES})
    uvlas_configure_target(uvlas_numerics_audit)
    uvlas_enable_allocation_audit(uvlas_numerics_audit)
    target_include_directories(uvlas_numerics_audit PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
    target_link_libraries(uvlas_numerics_audit PUBLIC uvlas_core_audit PRIVATE uvlas::eigen Spectra::Spectra Boost::headers Boost::hash2)
    uvlas_require_audited_dependencies(uvlas_numerics_audit)
    add_subdirectory(tests)
endif()

if(UVLAS_BUILD_BENCHMARKS)
    if(NOT "benchmarks" IN_LIST VCPKG_MANIFEST_FEATURES)
        message(FATAL_ERROR "UVLAS_BUILD_BENCHMARKS requires the vcpkg benchmarks feature")
    endif()
    find_package(benchmark CONFIG REQUIRED)
    add_subdirectory(benchmarks)
endif()

