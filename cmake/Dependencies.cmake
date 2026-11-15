# SPDX-License-Identifier: Apache-2.0
#
# Third-party packages come from the vcpkg manifest (vcpkg.json): the baseline
# commit and the overrides pin every version (N4-01, ADR 0001). The presets
# pass the vcpkg toolchain; nothing is fetched outside it.

find_package(Eigen3 CONFIG REQUIRED)
find_package(spectra CONFIG REQUIRED)
find_package(boost_multiprecision CONFIG REQUIRED)
find_package(boost_hash2 CONFIG REQUIRED)
find_package(boost_headers CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)

# One target carries the adopted numerical build policy for Eigen:
# MPL-2.0 parts only, and the single-threaded kernel as the reference path.
add_library(uvlas_eigen INTERFACE)
target_link_libraries(uvlas_eigen INTERFACE Eigen3::Eigen)
target_compile_definitions(uvlas_eigen INTERFACE EIGEN_MPL2_ONLY EIGEN_DONT_PARALLELIZE)
add_library(uvlas::eigen ALIAS uvlas_eigen)

# Carry the upstream notices verbatim, including transitive package notices.
# Distribution of these files with actual artifacts is implemented in IM-06.
set(UVLAS_THIRD_PARTY_SHARE "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/share")
set(UVLAS_THIRD_PARTY_NOTICES "${PROJECT_BINARY_DIR}/third-party-notices")
set(_uvlas_required_notices eigen3 spectra boost-multiprecision boost-hash2 nlohmann-json)
if(BUILD_TESTING)
    list(APPEND _uvlas_required_notices catch2)
endif()
if("python" IN_LIST VCPKG_MANIFEST_FEATURES)
    list(APPEND _uvlas_required_notices pybind11)
endif()
if(UVLAS_BUILD_BENCHMARKS)
    list(APPEND _uvlas_required_notices benchmark)
endif()
foreach(_uvlas_package IN LISTS _uvlas_required_notices)
    if(NOT EXISTS "${UVLAS_THIRD_PARTY_SHARE}/${_uvlas_package}/copyright")
        message(FATAL_ERROR "Missing third-party notice for ${_uvlas_package}")
    endif()
endforeach()
file(GLOB _uvlas_copyrights CONFIGURE_DEPENDS "${UVLAS_THIRD_PARTY_SHARE}/*/copyright")
if(NOT UVLAS_BUILD_BENCHMARKS)
    file(REMOVE "${UVLAS_THIRD_PARTY_NOTICES}/benchmark.txt")
endif()
foreach(_uvlas_copyright IN LISTS _uvlas_copyrights)
    get_filename_component(_uvlas_package_dir "${_uvlas_copyright}" DIRECTORY)
    get_filename_component(_uvlas_package "${_uvlas_package_dir}" NAME)
    if(_uvlas_package STREQUAL "benchmark" AND NOT UVLAS_BUILD_BENCHMARKS)
        continue()
    endif()
    configure_file("${_uvlas_copyright}" "${UVLAS_THIRD_PARTY_NOTICES}/${_uvlas_package}.txt" COPYONLY)
endforeach()
