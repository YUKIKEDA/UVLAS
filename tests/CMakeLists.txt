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
