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
        "-DBUILD_BENCHMARKS=${UVLAS_BUILD_BENCHMARKS}"
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
        "-DCTEST_COMMAND=${CMAKE_CTEST_COMMAND}"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/check_notice_configuration.cmake")
set_tests_properties(toolchain.notice_configuration PROPERTIES TIMEOUT 180 LABELS "integration;negative")

uvlas_add_test(core.failure
    LABELS unit negative
    SOURCES core/failure_test.cpp
    LIBS uvlas_core)
uvlas_add_test(core.error_certificate
    LABELS unit negative
    SOURCES core/error_certificate_test.cpp
    LIBS uvlas_core)

uvlas_add_test(core.checked_dimensions
    LABELS unit negative
    SOURCES core/checked_dimensions_test.cpp
    LIBS uvlas_core)
uvlas_add_test(core.resource_plan
    LABELS unit negative
    SOURCES core/resource_plan_test.cpp
    LIBS uvlas_core)

add_test(NAME core.identifier_compilation
    COMMAND "${CMAKE_COMMAND}"
        -S "${CMAKE_CURRENT_SOURCE_DIR}/core/identifier_compile"
        -B "${CMAKE_CURRENT_BINARY_DIR}/identifier-compile-$<CONFIG>"
        "-DPROJECT_INCLUDE=${PROJECT_SOURCE_DIR}/include"
        "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
set_tests_properties(core.identifier_compilation PROPERTIES TIMEOUT 120 LABELS "negative")

uvlas_add_test(core.instrumentation LABELS unit integration negative
    SOURCES core/instrumentation_test.cpp LIBS uvlas_core nlohmann_json::nlohmann_json)
uvlas_add_test(core.allocation_audit LABELS property negative
    SOURCES core/allocation_audit_test.cpp support/allocation_audit.cpp
    LIBS uvlas_core_audit uvlas::eigen)
uvlas_enable_allocation_audit(uvlas_test_core_allocation_audit)
uvlas_require_audited_dependencies(uvlas_test_core_allocation_audit)
add_executable(uvlas_eigen_allocation_probe core/eigen_allocation_probe.cpp support/allocation_audit.cpp)
uvlas_configure_target(uvlas_eigen_allocation_probe)
uvlas_enable_allocation_audit(uvlas_eigen_allocation_probe)
target_link_libraries(uvlas_eigen_allocation_probe PRIVATE uvlas::eigen)
add_test(NAME core.eigen_allocation_guard COMMAND "${CMAKE_COMMAND}"
    "-DPROBE=$<TARGET_FILE:uvlas_eigen_allocation_probe>"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/core/check_eigen_allocation.cmake")
set_tests_properties(core.eigen_allocation_guard PROPERTIES TIMEOUT 40 LABELS "negative")
add_test(NAME core.audit_configuration COMMAND "${CMAKE_COMMAND}"
    "-DSOURCE=${CMAKE_CURRENT_SOURCE_DIR}/core/audit_configuration"
    "-DBINARY=${CMAKE_CURRENT_BINARY_DIR}/audit-configuration-$<CONFIG>"
    "-DAUDIT_MODULE=${PROJECT_SOURCE_DIR}/cmake/AllocationAudit.cmake"
    "-DCOMPILER=${CMAKE_CXX_COMPILER}"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/core/check_audit_configuration.cmake")
set_tests_properties(core.audit_configuration PROPERTIES TIMEOUT 100 LABELS "negative")

find_package(Threads REQUIRED)
uvlas_add_test(core.state
    LABELS unit negative
    SOURCES core/state_test.cpp
    LIBS uvlas_core Threads::Threads)

add_test(NAME tools.precision_witness_export
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/qa/export_precision_witness_test.py")
set_tests_properties(tools.precision_witness_export PROPERTIES TIMEOUT 10 LABELS "unit;negative")

uvlas_add_test(numerics.linear_validation
    LABELS unit negative verification
    SOURCES numerics/linear_validation_test.cpp
    LIBS uvlas_numerics)
target_include_directories(uvlas_test_numerics_linear_validation PRIVATE "${PROJECT_SOURCE_DIR}/src/numerics")

uvlas_add_test(numerics.linear_result
    LABELS unit negative verification
    SOURCES numerics/linear_result_test.cpp
    LIBS uvlas_numerics)
target_include_directories(uvlas_test_numerics_linear_result PRIVATE "${PROJECT_SOURCE_DIR}/src/numerics")

uvlas_add_test(numerics.linear_metrics
    LABELS unit negative verification
    SOURCES numerics/linear_metrics_test.cpp
    LIBS uvlas_numerics Boost::multiprecision)
target_include_directories(uvlas_test_numerics_linear_metrics PRIVATE "${PROJECT_SOURCE_DIR}/src/numerics")

uvlas_add_test(numerics.backend_memory
    LABELS unit negative property
    SOURCES numerics/backend_memory_test.cpp support/allocation_audit.cpp
    LIBS uvlas_numerics_audit uvlas::eigen)
target_include_directories(uvlas_test_numerics_backend_memory PRIVATE "${PROJECT_SOURCE_DIR}/src/numerics")
uvlas_enable_allocation_audit(uvlas_test_numerics_backend_memory)
uvlas_require_audited_dependencies(uvlas_test_numerics_backend_memory)

uvlas_add_test(numerics.precision
    LABELS unit negative verification
    SOURCES numerics/precision_test.cpp
    LIBS uvlas_numerics Boost::multiprecision nlohmann_json::nlohmann_json)
target_compile_definitions(uvlas_test_numerics_precision PRIVATE
    UVLAS_PRECISION_WITNESS_PATH="${CMAKE_CURRENT_SOURCE_DIR}/numerics/fixtures/precision-4096-row-v1.json")

uvlas_add_test(numerics.storage LABELS unit negative
    SOURCES numerics/storage_test.cpp LIBS uvlas_numerics)
uvlas_add_test(numerics.eigen_views LABELS integration negative
    SOURCES numerics/eigen_views_test.cpp LIBS uvlas_numerics uvlas::eigen)
target_include_directories(uvlas_test_numerics_eigen_views PRIVATE "${PROJECT_SOURCE_DIR}/src")

add_test(NAME numerics.storage_compilation
    COMMAND "${CMAKE_COMMAND}"
        -S "${CMAKE_CURRENT_SOURCE_DIR}/numerics/borrow_compile"
        -B "${CMAKE_CURRENT_BINARY_DIR}/storage-compile-$<CONFIG>"
        "-DPROJECT_ROOT=${PROJECT_SOURCE_DIR}"
        "-DCMAKE_PREFIX_PATH=${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}"
        "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
set_tests_properties(numerics.storage_compilation PROPERTIES TIMEOUT 120 LABELS "negative")

uvlas_add_test(numerics.sparse_spd
    LABELS unit verification negative
    SOURCES numerics/sparse_spd_test.cpp
    LIBS uvlas_numerics Boost::multiprecision)
target_include_directories(uvlas_test_numerics_sparse_spd PRIVATE "${PROJECT_SOURCE_DIR}/src/numerics")

uvlas_add_test(numerics.sparse_triangular
    LABELS unit verification negative
    SOURCES numerics/sparse_triangular_test.cpp
    LIBS uvlas_numerics uvlas::eigen)
target_include_directories(uvlas_test_numerics_sparse_triangular PRIVATE "${PROJECT_SOURCE_DIR}/src/numerics")

uvlas_add_test(numerics.sparse_solve
    LABELS unit integration verification negative
    SOURCES numerics/sparse_solve_test.cpp
    LIBS uvlas_numerics Boost::multiprecision)

uvlas_add_test(numerics.sparse_capacity
    LABELS capacity verification
    SOURCES numerics/sparse_capacity_test.cpp
    LIBS uvlas_numerics)

uvlas_add_test(numerics.sparse_allocation
    LABELS property verification
    SOURCES numerics/sparse_allocation_test.cpp support/allocation_audit.cpp
    LIBS uvlas_numerics_audit uvlas::eigen)
target_include_directories(uvlas_test_numerics_sparse_allocation PRIVATE "${PROJECT_SOURCE_DIR}/src/numerics")
uvlas_enable_allocation_audit(uvlas_test_numerics_sparse_allocation)
uvlas_require_audited_dependencies(uvlas_test_numerics_sparse_allocation)

uvlas_add_test(numerics.dense_solve LABELS unit verification negative
 SOURCES numerics/dense_solve_test.cpp LIBS uvlas_numerics Boost::multiprecision)

uvlas_add_test(numerics.dense_contract LABELS unit negative verification
 SOURCES numerics/dense_contract_test.cpp LIBS uvlas_numerics)

uvlas_add_test(numerics.dense_allocations LABELS property verification
 SOURCES numerics/dense_allocations_test.cpp support/allocation_audit.cpp
 LIBS uvlas_numerics_audit uvlas::eigen)
uvlas_enable_allocation_audit(uvlas_test_numerics_dense_allocations)
uvlas_require_audited_dependencies(uvlas_test_numerics_dense_allocations)

uvlas_add_test(numerics.gmres_contract LABELS unit negative verification
    SOURCES numerics/gmres_contract_test.cpp LIBS uvlas_numerics Boost::multiprecision)
uvlas_add_test(numerics.gmres_solve LABELS unit verification
    SOURCES numerics/gmres_solve_test.cpp LIBS uvlas_numerics)
target_include_directories(uvlas_test_numerics_gmres_solve PRIVATE "${PROJECT_SOURCE_DIR}/src/numerics")
uvlas_add_test(numerics.gmres_residual LABELS unit negative verification
    SOURCES numerics/gmres_residual_test.cpp LIBS uvlas_numerics Boost::multiprecision)
target_include_directories(uvlas_test_numerics_gmres_residual PRIVATE "${PROJECT_SOURCE_DIR}/src/numerics")
uvlas_add_test(numerics.gmres_allocations LABELS unit negative verification
    SOURCES numerics/gmres_allocations_test.cpp support/allocation_audit.cpp
    LIBS uvlas_numerics_audit uvlas::eigen)
uvlas_enable_allocation_audit(uvlas_test_numerics_gmres_allocations)
uvlas_require_audited_dependencies(uvlas_test_numerics_gmres_allocations)
uvlas_add_test(numerics.foundation_reference
    LABELS integration verification negative
    SOURCES numerics/foundation_reference_test.cpp
    LIBS uvlas_numerics Boost::multiprecision Boost::hash2 nlohmann_json::nlohmann_json)
target_compile_definitions(uvlas_test_numerics_foundation_reference PRIVATE
    UVLAS_FOUNDATION_REFERENCE_ROOT="${CMAKE_CURRENT_SOURCE_DIR}/numerics/fixtures/j1-pg-free-history-24-v1")

add_test(NAME tools.foundation_reference
    COMMAND "${Python3_EXECUTABLE}" -m unittest tools.qa.foundation_reference_test)
add_test(NAME tools.foundation_reference_verify
    COMMAND "${Python3_EXECUTABLE}" -m tools.qa.foundation_reference verify
        "${CMAKE_CURRENT_SOURCE_DIR}/numerics/fixtures/j1-pg-free-history-24-v1")
set_tests_properties(tools.foundation_reference tools.foundation_reference_verify
    PROPERTIES WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}" TIMEOUT 30 LABELS "unit;verification;negative")

if(WIN32)
    find_program(UVLAS_PWSH_EXECUTABLE NAMES pwsh REQUIRED)
    add_test(NAME tools.msvc_arguments
        COMMAND "${UVLAS_PWSH_EXECUTABLE}" -NoProfile -File
            "${PROJECT_SOURCE_DIR}/tools/dev/msvc_test.ps1")
    set_tests_properties(tools.msvc_arguments PROPERTIES TIMEOUT 75 LABELS "unit;negative")
endif()

uvlas_add_test(numerics.linear_candidate LABELS unit negative verification
    SOURCES numerics/linear_candidate_test.cpp LIBS uvlas_numerics)
target_include_directories(uvlas_test_numerics_linear_candidate PRIVATE "${PROJECT_SOURCE_DIR}/src/numerics")
uvlas_add_test(numerics.direct_progress LABELS unit negative verification
    SOURCES numerics/direct_progress_test.cpp LIBS uvlas_numerics Boost::multiprecision)
