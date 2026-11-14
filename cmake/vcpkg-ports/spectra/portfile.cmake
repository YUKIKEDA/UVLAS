# Spectra 1.2.0 is the adopted version (N4-07). The vcpkg registry stops at
# 1.0.1, so this overlay port pins the release archive. Its SHA-256 is
# fdfccacba1c77d9b4ffefae7258c760c99e3c8a2823ca87ea5b11a50d297a73b, the archive
# used by the numerical spike of design work 4.
vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO yixuan/spectra
    REF v1.2.0
    SHA512 dfd8824d33853dc23c476af4ca8458393284e1d50b8524c000b14109df0095257c3975881d01286bdd057950366d0ff498afb63fbde02f9f313122212ddf62a6
    HEAD_REF master
)

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
)

vcpkg_cmake_install()
# 1.2.0 installs its package files to share/cmake/Spectra (1.0.1 used share/spectra/cmake).
vcpkg_cmake_config_fixup(PACKAGE_NAME spectra CONFIG_PATH share/cmake/Spectra)

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug")

file(INSTALL "${SOURCE_PATH}/LICENSE" DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}" RENAME copyright)
