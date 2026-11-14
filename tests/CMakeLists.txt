find_package(Catch2 3 CONFIG REQUIRED)

# uvlas_add_test(<name> LABELS <label>... SOURCES <file>... [LIBS <target>...] [TIMEOUT <s>])
# One executable is one CTest test. Labels: unit | property | integration |
# verification | negative, plus the acceptance case it serves (for example V-PG).
function(uvlas_add_test name)
    cmake_parse_arguments(arg "" "TIMEOUT" "LABELS;SOURCES;LIBS" ${ARGN})
    if(NOT arg_LABELS OR NOT arg_SOURCES)
        message(FATAL_ERROR "uvlas_add_test(${name}) needs LABELS and SOURCES")
    endif()
    string(REPLACE "." "_" target "uvlas_test_${name}")
    add_executable(${target} ${arg_SOURCES})
    uvlas_configure_target(${target})
    target_link_libraries(${target} PRIVATE Catch2::Catch2WithMain ${arg_LIBS})
    add_test(NAME ${name} COMMAND ${target})
    if(NOT arg_TIMEOUT)
        set(arg_TIMEOUT 60)
    endif()
    set_tests_properties(${name} PROPERTIES LABELS "${arg_LABELS}" TIMEOUT ${arg_TIMEOUT})
endfunction()

foreach(case IN ITEMS version help missing-option unknown-option extra-option)
    add_test(NAME cli.${case}
        COMMAND "${CMAKE_COMMAND}"
            "-DUVLAS_EXECUTABLE=$<TARGET_FILE:uvlas_cli>"
            "-DCASE=${case}"
            "-DEXPECTED_VERSION=${UVLAS_VERSION}"
            "-DEXPECTED_CONFIGURATION=$<CONFIG>"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/check_cli.cmake")
    set_tests_properties(cli.${case} PROPERTIES TIMEOUT 10 LABELS "integration")
endforeach()

# QA-04: results must equal the baseline of this platform bit for bit.
find_package(Python3 3.10 REQUIRED COMPONENTS Interpreter)
add_test(NAME qa.numeric_canary
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/qa/results.py" check
        --case numeric-canary --run "$<TARGET_FILE:uvlas_numeric_canary>")
set_tests_properties(qa.numeric_canary PROPERTIES TIMEOUT 120 LABELS "verification")

add_test(NAME cli.build_info COMMAND uvlas_cli --build-info)
set_tests_properties(cli.build_info PROPERTIES TIMEOUT 10 LABELS "integration"
    PASS_REGULAR_EXPRESSION "\"schema\": \"uvlas.build-info/1\".*\"commit\": \"[0-9a-f]+\".*\"platform_key\": \"[a-z0-9.-]+\".*\"eigen\": \"5\\.0\\.1\"")

uvlas_add_test(toolchain.pinned_dependencies
    LABELS unit
    SOURCES toolchain/pinned_dependencies_test.cpp
    LIBS uvlas::eigen Spectra::Spectra Boost::multiprecision Boost::hash2
         nlohmann_json::nlohmann_json)
target_compile_definitions(uvlas_test_toolchain_pinned_dependencies PRIVATE
    UVLAS_MANIFEST_PATH="${PROJECT_SOURCE_DIR}/vcpkg.json")

uvlas_add_test(toolchain.linear_algebra
    LABELS unit verification
    SOURCES toolchain/linear_algebra_test.cpp
    LIBS uvlas::eigen)
uvlas_add_test(toolchain.arnoldi
    LABELS unit verification
    SOURCES toolchain/arnoldi_test.cpp
    LIBS uvlas::eigen Spectra::Spectra)
uvlas_add_test(toolchain.precise_arithmetic
    LABELS unit verification
    SOURCES toolchain/precise_arithmetic_test.cpp
    LIBS Boost::multiprecision Boost::hash2)

add_test(NAME toolchain.third_party_notices
    COMMAND "${CMAKE_COMMAND}"
        "-DSHARE_DIR=${UVLAS_THIRD_PARTY_SHARE}"
        "-DNOTICES_DIR=${UVLAS_THIRD_PARTY_NOTICES}"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/check_notices.cmake")
set_tests_properties(toolchain.third_party_notices PROPERTIES TIMEOUT 10 LABELS "integration")

add_test(NAME toolchain.notice_configuration
    COMMAND "${CMAKE_COMMAND}"
        "-DSOURCE_ROOT=${PROJECT_SOURCE_DIR}"
        "-DBINARY_ROOT=${CMAKE_CURRENT_BINARY_DIR}"
        "-DDEPENDENCY_PREFIX=${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}"
        "-DTRIPLET=${VCPKG_TARGET_TRIPLET}"
        "-DCXX_COMPILER=${CMAKE_CXX_COMPILER}"
        "-DPYTHON_EXECUTABLE=${Python3_EXECUTABLE}"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/check_notice_configuration.cmake")
set_tests_properties(toolchain.notice_configuration PROPERTIES TIMEOUT 180 LABELS "integration;negative")
