# SPDX-License-Identifier: Apache-2.0
# This component constructs a private verification image from native Release targets.
if(NOT UVLAS_BUILD_PYTHON OR CMAKE_CROSSCOMPILING OR NOT UVLAS_ARCH_KEY STREQUAL "x64" OR
   NOT (WIN32 OR CMAKE_SYSTEM_NAME STREQUAL "Linux"))
    message(FATAL_ERROR "Private runtime image requires a native x64 Windows or Linux Python build")
endif()

set(_uvlas_private_component uvlas-private-runtime)
install(TARGETS uvlas_cli
    RUNTIME DESTINATION bin CONFIGURATIONS Release COMPONENT ${_uvlas_private_component} EXCLUDE_FROM_ALL)
install(TARGETS uvlas_python
    LIBRARY DESTINATION python CONFIGURATIONS Release COMPONENT ${_uvlas_private_component} EXCLUDE_FROM_ALL
    RUNTIME DESTINATION python CONFIGURATIONS Release COMPONENT ${_uvlas_private_component} EXCLUDE_FROM_ALL)
install(FILES LICENSE THIRD_PARTY_NOTICES.md
    DESTINATION share/uvlas CONFIGURATIONS Release
    COMPONENT ${_uvlas_private_component} EXCLUDE_FROM_ALL)
install(DIRECTORY "${UVLAS_THIRD_PARTY_NOTICES}/"
    DESTINATION share/uvlas/third-party-notices CONFIGURATIONS Release
    COMPONENT ${_uvlas_private_component} EXCLUDE_FROM_ALL
    FILES_MATCHING PATTERN "*.txt")
