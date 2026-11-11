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
