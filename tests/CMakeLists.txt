find_package(Catch2 3 CONFIG REQUIRED)

# uvlas_add_test(<name> LABELS <label>... SOURCES <file>... [LIBS <target>...]
#                [TIMEOUT <s>] [CONFIGURATIONS <configuration>...])
# One executable is one CTest test. Labels: unit | property | integration |
# verification | negative, plus the acceptance case it serves (for example V-PG).
function(uvlas_add_test name)
    cmake_parse_arguments(arg "" "TIMEOUT" "LABELS;SOURCES;LIBS;CONFIGURATIONS" ${ARGN})
    if(NOT arg_LABELS OR NOT arg_SOURCES)
        message(FATAL_ERROR "uvlas_add_test(${name}) needs LABELS and SOURCES")
    endif()
    string(REPLACE "." "_" target "uvlas_test_${name}")
    add_executable(${target} ${arg_SOURCES})
    uvlas_configure_target(${target})
    target_link_libraries(${target} PRIVATE Catch2::Catch2WithMain ${arg_LIBS})
    set(configuration_args)
    if(arg_CONFIGURATIONS)
        if(CMAKE_CONFIGURATION_TYPES)
            list(APPEND configuration_args CONFIGURATIONS ${arg_CONFIGURATIONS})
        elseif(NOT CMAKE_BUILD_TYPE IN_LIST arg_CONFIGURATIONS)
            return()
        endif()
    endif()
    # Single-config CTest normally has no -C argument. CONFIGURATIONS there
    # would silently omit even the matching Release test; select at configure time.
    add_test(NAME ${name} COMMAND ${target} ${configuration_args})
    if(NOT arg_TIMEOUT)
        set(arg_TIMEOUT 60)
    endif()
    set_tests_properties(${name} PROPERTIES LABELS "${arg_LABELS}" TIMEOUT ${arg_TIMEOUT})
endfunction()

uvlas_add_test(io.input_scalar
    LABELS unit property negative verification
    SOURCES io/input_scalar_test.cpp
    LIBS uvlas_io Threads::Threads)
uvlas_add_test(io.input_scalar_allocations
    LABELS property negative
    SOURCES io/input_scalar_allocations_test.cpp support/allocation_audit.cpp
    LIBS uvlas_io_audit uvlas::eigen)
uvlas_enable_allocation_audit(uvlas_test_io_input_scalar_allocations)
uvlas_require_audited_dependencies(uvlas_test_io_input_scalar_allocations)

uvlas_add_test(io.input_document
    LABELS unit property negative verification
    SOURCES io/input_document_test.cpp
    LIBS uvlas_io Threads::Threads)
target_include_directories(uvlas_test_io_input_document PRIVATE "${PROJECT_SOURCE_DIR}/src/io")
uvlas_add_test(io.input_document_resources
    LABELS property negative verification
    SOURCES io/input_document_resource_test.cpp
    LIBS uvlas_io_audit uvlas::eigen)
target_include_directories(uvlas_test_io_input_document_resources PRIVATE "${PROJECT_SOURCE_DIR}/src/io")
uvlas_enable_allocation_audit(uvlas_test_io_input_document_resources)
uvlas_require_audited_dependencies(uvlas_test_io_input_document_resources)

uvlas_add_test(io.reference_path
    LABELS unit property negative verification
    SOURCES io/reference_path_test.cpp support/allocation_audit.cpp
    LIBS uvlas_io_audit uvlas::eigen nlohmann_json::nlohmann_json)
target_include_directories(uvlas_test_io_reference_path PRIVATE "${PROJECT_SOURCE_DIR}/src/io")
target_compile_definitions(uvlas_test_io_reference_path PRIVATE
    UVLAS_REFERENCE_PATH_FIXTURE="${PROJECT_SOURCE_DIR}/docs/verification/IM-06/reference-snapshot-design-cases.json")
uvlas_enable_allocation_audit(uvlas_test_io_reference_path)
uvlas_require_audited_dependencies(uvlas_test_io_reference_path)

uvlas_add_test(io.reference_file_snapshot
    LABELS unit property negative verification
    SOURCES io/reference_file_snapshot_test.cpp support/allocation_audit.cpp
    LIBS uvlas_io_audit uvlas::eigen)
target_include_directories(uvlas_test_io_reference_file_snapshot PRIVATE "${PROJECT_SOURCE_DIR}/src/io")
uvlas_enable_allocation_audit(uvlas_test_io_reference_file_snapshot)
uvlas_require_audited_dependencies(uvlas_test_io_reference_file_snapshot)

if(WIN32 OR CMAKE_SYSTEM_NAME STREQUAL "Linux")
    uvlas_add_test(io.reference_input
        LABELS integration property negative verification
        SOURCES io/reference_input_test.cpp
        LIBS uvlas_io nlohmann_json::nlohmann_json)
    target_include_directories(uvlas_test_io_reference_input PRIVATE "${PROJECT_SOURCE_DIR}/src/io")
endif()

if(WIN32)
    uvlas_add_test(io.windows_reference_root
        LABELS unit property negative verification
        SOURCES io/windows_reference_root_test.cpp
        LIBS uvlas_io)
    target_include_directories(uvlas_test_io_windows_reference_root PRIVATE "${PROJECT_SOURCE_DIR}/src/io")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    uvlas_add_test(io.linux_reference_root
        LABELS unit property negative verification
        SOURCES io/linux_reference_root_test.cpp
        LIBS uvlas_io)
    target_include_directories(uvlas_test_io_linux_reference_root PRIVATE "${PROJECT_SOURCE_DIR}/src/io")
endif()

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
add_test(NAME tools.trace
    COMMAND "${Python3_EXECUTABLE}" -m unittest tools.qa.trace_test)
set_tests_properties(tools.trace PROPERTIES TIMEOUT 30 LABELS "unit;verification"
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}" ENVIRONMENT "PYTHONUTF8=1")
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

# Isolate injected allocation failures from the test runner and exercise
# optional-NRVO-disabled MSVC Debug paths through the installed JSON package.
add_executable(uvlas_json_allocation_probe toolchain/json_allocation_probe.cpp)
uvlas_configure_target(uvlas_json_allocation_probe)
target_link_libraries(uvlas_json_allocation_probe PRIVATE nlohmann_json::nlohmann_json)
if(MSVC)
    target_compile_options(uvlas_json_allocation_probe PRIVATE /Zc:nrvo-)
endif()

add_test(NAME toolchain.json_allocation_exceptions
    COMMAND "${Python3_EXECUTABLE}" -m tools.qa.json_allocation_probe_test
        --probe "$<TARGET_FILE:uvlas_json_allocation_probe>")
set_tests_properties(toolchain.json_allocation_exceptions PROPERTIES
    TIMEOUT 180 LABELS "negative;verification"
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}" ENVIRONMENT "PYTHONUTF8=1")

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

uvlas_add_test(aero.exact_triangle_predicates
    LABELS unit property negative verification V-AERO
    SOURCES aero/exact_triangle_predicates_test.cpp
    LIBS uvlas_aero Boost::multiprecision Threads::Threads)
target_include_directories(uvlas_test_aero_exact_triangle_predicates
    PRIVATE "${PROJECT_SOURCE_DIR}/src")

uvlas_add_test(aero.exact_triangle_allocations
    LABELS property negative verification
    SOURCES aero/exact_triangle_allocations_test.cpp support/allocation_audit.cpp
    LIBS uvlas_aero_audit uvlas::eigen Boost::headers)
target_include_directories(uvlas_test_aero_exact_triangle_allocations
    PRIVATE "${PROJECT_SOURCE_DIR}/src")
uvlas_enable_allocation_audit(uvlas_test_aero_exact_triangle_allocations)
uvlas_require_audited_dependencies(uvlas_test_aero_exact_triangle_allocations)

uvlas_add_test(aero.exact_triangle_resources
    LABELS property negative verification
    SOURCES aero/exact_triangle_resource_test.cpp
    LIBS uvlas_aero_audit Boost::headers)
target_include_directories(uvlas_test_aero_exact_triangle_resources
    PRIVATE "${PROJECT_SOURCE_DIR}/src")
uvlas_enable_allocation_audit(uvlas_test_aero_exact_triangle_resources)
uvlas_require_audited_dependencies(uvlas_test_aero_exact_triangle_resources)

uvlas_add_test(aero.geometry_snapshot
    LABELS unit property negative verification V-AERO
    SOURCES aero/geometry_snapshot_test.cpp
    LIBS uvlas_aero Boost::multiprecision Threads::Threads)
target_include_directories(uvlas_test_aero_geometry_snapshot
    PRIVATE "${PROJECT_SOURCE_DIR}/src")
uvlas_add_test(aero.geometry_snapshot_resources
    LABELS property negative verification V-AERO
    SOURCES aero/geometry_snapshot_resources_test.cpp
    LIBS uvlas_aero_audit Boost::headers)
target_include_directories(uvlas_test_aero_geometry_snapshot_resources
    PRIVATE "${PROJECT_SOURCE_DIR}/src")
uvlas_enable_allocation_audit(uvlas_test_aero_geometry_snapshot_resources)
uvlas_require_audited_dependencies(uvlas_test_aero_geometry_snapshot_resources)

uvlas_add_test(aero.local_geometry_packet
    LABELS unit property negative verification V-AERO
    SOURCES aero/local_geometry_packet_test.cpp
    LIBS uvlas_aero Boost::multiprecision Threads::Threads)
target_include_directories(uvlas_test_aero_local_geometry_packet
    PRIVATE "${PROJECT_SOURCE_DIR}/src")
uvlas_add_test(aero.local_geometry_packet_resources
    LABELS property negative verification V-AERO
    SOURCES aero/local_geometry_packet_resources_test.cpp
    LIBS uvlas_aero_audit Boost::headers)
target_include_directories(uvlas_test_aero_local_geometry_packet_resources
    PRIVATE "${PROJECT_SOURCE_DIR}/src")
uvlas_enable_allocation_audit(uvlas_test_aero_local_geometry_packet_resources)
uvlas_require_audited_dependencies(uvlas_test_aero_local_geometry_packet_resources)

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

add_test(NAME tools.precision_system_export
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/qa/export_precision_system_test.py")
set_tests_properties(tools.precision_system_export PROPERTIES TIMEOUT 30 LABELS "unit;negative"
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}" ENVIRONMENT "PYTHONUTF8=1")

configure_file(numerics/fixtures/precision-4096-system-v1.json
    "${CMAKE_CURRENT_BINARY_DIR}/precision-4096-system-v1.json" COPYONLY)
uvlas_add_test(numerics.full_precision_control
    LABELS unit verification negative
    SOURCES numerics/full_precision_control_test.cpp
    LIBS Boost::multiprecision Boost::hash2 nlohmann_json::nlohmann_json)
target_compile_definitions(uvlas_test_numerics_full_precision_control PRIVATE
    UVLAS_FULL_PRECISION_FIXTURE="${CMAKE_CURRENT_BINARY_DIR}/precision-4096-system-v1.json")
uvlas_add_test(numerics.foundation_precision_system
    LABELS capacity integration verification negative
    SOURCES numerics/foundation_precision_system_test.cpp
    LIBS uvlas_numerics Boost::multiprecision Boost::hash2 nlohmann_json::nlohmann_json
    CONFIGURATIONS Release
    TIMEOUT 3600)
target_compile_definitions(uvlas_test_numerics_foundation_precision_system PRIVATE
    UVLAS_FULL_PRECISION_FIXTURE="${CMAKE_CURRENT_BINARY_DIR}/precision-4096-system-v1.json")

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

uvlas_add_test(numerics.foundation_precision_rank
    LABELS integration verification negative
    SOURCES numerics/foundation_precision_rank_test.cpp
    LIBS uvlas_numerics Boost::multiprecision)

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

uvlas_add_test(structure.exact_polynomial
    LABELS unit property negative verification V-PG
    SOURCES structure/exact_polynomial_test.cpp
    LIBS uvlas_structure Boost::multiprecision Threads::Threads nlohmann_json::nlohmann_json)
target_include_directories(uvlas_test_structure_exact_polynomial
    PRIVATE "${PROJECT_SOURCE_DIR}/src/structure")
target_compile_definitions(uvlas_test_structure_exact_polynomial PRIVATE
    UVLAS_EXACT_POLYNOMIAL_FIXTURE=\"${PROJECT_SOURCE_DIR}/tests/structure/fixtures/exact-polynomial-v1.json\")

uvlas_add_test(structure.exact_polynomial_resources
    LABELS property negative verification
    SOURCES structure/exact_polynomial_resource_test.cpp
    LIBS uvlas_structure_audit)
uvlas_enable_allocation_audit(uvlas_test_structure_exact_polynomial_resources)
uvlas_require_audited_dependencies(uvlas_test_structure_exact_polynomial_resources)

uvlas_add_test(structure.pg_graph
    LABELS unit property negative verification V-PG
    SOURCES structure/pg_graph_test.cpp
    LIBS uvlas_structure Boost::multiprecision Threads::Threads nlohmann_json::nlohmann_json)
target_include_directories(uvlas_test_structure_pg_graph
    PRIVATE "${PROJECT_SOURCE_DIR}/src/structure")
target_compile_definitions(uvlas_test_structure_pg_graph PRIVATE
    UVLAS_PG_GRAPH_FIXTURE=\"${PROJECT_SOURCE_DIR}/tests/structure/fixtures/pg-graph-v1.json\")
uvlas_add_test(structure.pg_graph_resources
    LABELS property negative verification V-PG
    SOURCES structure/pg_graph_resource_test.cpp
    LIBS uvlas_structure_audit)
uvlas_enable_allocation_audit(uvlas_test_structure_pg_graph_resources)
uvlas_require_audited_dependencies(uvlas_test_structure_pg_graph_resources)

uvlas_add_test(structure.pg_coefficients
    TIMEOUT 300
    LABELS unit property negative verification V-PG
    SOURCES structure/pg_coefficients_test.cpp structure/pg_coefficients_reference_test.cpp
    LIBS uvlas_structure Boost::multiprecision Threads::Threads nlohmann_json::nlohmann_json)
# The instrumented exact-coefficient test exceeded the ordinary 300-second deadline.
# Keep the normal limit and a separate finite budget for instrumented execution.
if(UVLAS_SANITIZERS AND NOT MSVC)
    set_tests_properties(structure.pg_coefficients PROPERTIES TIMEOUT 900)
endif()
target_include_directories(uvlas_test_structure_pg_coefficients
    PRIVATE "${PROJECT_SOURCE_DIR}/src/structure")
target_compile_definitions(uvlas_test_structure_pg_coefficients PRIVATE
    UVLAS_PG_COEFFICIENTS_FIXTURE=\"${PROJECT_SOURCE_DIR}/tests/structure/fixtures/pg-coefficients-v1.json\")
uvlas_add_test(structure.pg_coefficients_resources
    LABELS property negative verification V-PG
    SOURCES structure/pg_coefficients_resource_test.cpp
    LIBS uvlas_structure_audit)
target_include_directories(uvlas_test_structure_pg_coefficients_resources
    PRIVATE "${PROJECT_SOURCE_DIR}/src/structure")
uvlas_enable_allocation_audit(uvlas_test_structure_pg_coefficients_resources)
uvlas_require_audited_dependencies(uvlas_test_structure_pg_coefficients_resources)

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

uvlas_add_test(numerics.eigen_generalized LABELS unit negative verification
    SOURCES numerics/eigen_generalized_test.cpp LIBS uvlas_numerics Boost::multiprecision)

add_test(NAME qa.eigen_generalized
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/qa/results.py" check
        --case eigen-generalized --run "$<TARGET_FILE:uvlas_eigen_generalized_canary>")
set_tests_properties(qa.eigen_generalized PROPERTIES TIMEOUT 30 LABELS "verification")

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

uvlas_add_test(numerics.eigen_full_phase_audit LABELS property verification
    SOURCES numerics/eigen_full_phase_audit_test.cpp support/allocation_audit.cpp
    LIBS uvlas_numerics_audit uvlas::eigen Boost::headers)
target_include_directories(uvlas_test_numerics_eigen_full_phase_audit PRIVATE "${PROJECT_SOURCE_DIR}/src/numerics")
uvlas_enable_allocation_audit(uvlas_test_numerics_eigen_full_phase_audit)
uvlas_require_audited_dependencies(uvlas_test_numerics_eigen_full_phase_audit)

uvlas_add_test(numerics.eigen_full_contract LABELS unit negative verification
    SOURCES numerics/eigen_full_contract_test.cpp LIBS uvlas_numerics)

uvlas_add_test(numerics.eigen_full LABELS unit negative verification
    SOURCES numerics/eigen_full_test.cpp LIBS uvlas_numerics Boost::multiprecision nlohmann_json::nlohmann_json)
target_include_directories(uvlas_test_numerics_eigen_full PRIVATE "${PROJECT_SOURCE_DIR}/src/numerics")
target_compile_definitions(uvlas_test_numerics_eigen_full PRIVATE
    UVLAS_FULL_FIXTURE_PATH="${PROJECT_SOURCE_DIR}/tests/numerics/fixtures/j1-pg-free-history-24-v1/fixture.json")

uvlas_add_test(numerics.dense_contract LABELS unit negative verification
 SOURCES numerics/dense_contract_test.cpp LIBS uvlas_numerics)

uvlas_add_test(numerics.dense_allocations LABELS property verification
 SOURCES numerics/dense_allocations_test.cpp support/allocation_audit.cpp
 LIBS uvlas_numerics_audit uvlas::eigen)
uvlas_enable_allocation_audit(uvlas_test_numerics_dense_allocations)
uvlas_require_audited_dependencies(uvlas_test_numerics_dense_allocations)

uvlas_add_test(numerics.eigen_arnoldi_allocations LABELS unit negative verification
    SOURCES numerics/eigen_arnoldi_allocations_test.cpp support/allocation_audit.cpp
    LIBS uvlas_numerics_audit uvlas::eigen Spectra::Spectra)
target_include_directories(uvlas_test_numerics_eigen_arnoldi_allocations PRIVATE "${PROJECT_SOURCE_DIR}/src/numerics")
uvlas_enable_allocation_audit(uvlas_test_numerics_eigen_arnoldi_allocations)
uvlas_require_audited_dependencies(uvlas_test_numerics_eigen_arnoldi_allocations)
uvlas_add_test(numerics.eigen_arnoldi LABELS unit negative verification
    SOURCES numerics/eigen_arnoldi_test.cpp LIBS uvlas_numerics uvlas::eigen Spectra::Spectra Boost::multiprecision nlohmann_json::nlohmann_json)
target_include_directories(uvlas_test_numerics_eigen_arnoldi PRIVATE "${PROJECT_SOURCE_DIR}/src/numerics")
target_compile_definitions(uvlas_test_numerics_eigen_arnoldi PRIVATE
    UVLAS_ARNOLDI_REFERENCE_ROOT="${CMAKE_CURRENT_SOURCE_DIR}/numerics/fixtures/j1-pg-free-history-24-v1")
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

uvlas_add_test(numerics.eigen_common LABELS unit negative verification
    SOURCES numerics/eigen_common_test.cpp LIBS uvlas_numerics)
uvlas_add_test(numerics.eigen_metrics LABELS unit negative verification
    SOURCES numerics/eigen_metrics_test.cpp LIBS uvlas_numerics Boost::multiprecision)
uvlas_add_test(numerics.eigen_action LABELS unit negative verification
    SOURCES numerics/eigen_action_test.cpp LIBS uvlas_numerics)

add_test(NAME qa.eigen_common
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/qa/results.py" check
        --case eigen-common --run "$<TARGET_FILE:uvlas_eigen_common_canary>")
set_tests_properties(qa.eigen_common PROPERTIES TIMEOUT 30 LABELS "verification")

add_test(NAME qa.eigen_full
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/qa/results.py" check
        --case eigen-full --run "$<TARGET_FILE:uvlas_eigen_full_canary>")
set_tests_properties(qa.eigen_full PROPERTIES TIMEOUT 30 LABELS "verification")

add_test(NAME qa.eigen_arnoldi
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/qa/results.py" check
        --case eigen-arnoldi --run "$<TARGET_FILE:uvlas_eigen_arnoldi_canary>")
set_tests_properties(qa.eigen_arnoldi PROPERTIES TIMEOUT 30 LABELS "verification")


uvlas_add_test(structure.pg_integration_partition
    LABELS unit property negative verification V-PG
    SOURCES structure/pg_integration_partition_test.cpp
    LIBS uvlas_structure Boost::multiprecision Threads::Threads nlohmann_json::nlohmann_json)
target_include_directories(uvlas_test_structure_pg_integration_partition
    PRIVATE "${PROJECT_SOURCE_DIR}/src/structure")
target_compile_definitions(uvlas_test_structure_pg_integration_partition PRIVATE
    UVLAS_PG_PARTITION_FIXTURE="${PROJECT_SOURCE_DIR}/tests/structure/fixtures/pg-partition-v1.json")
uvlas_add_test(structure.pg_integration_partition_resources
    LABELS property negative verification V-PG
    SOURCES structure/pg_integration_partition_resource_test.cpp
    LIBS uvlas_structure_audit)
target_include_directories(uvlas_test_structure_pg_integration_partition_resources
    PRIVATE "${PROJECT_SOURCE_DIR}/src/structure")
uvlas_enable_allocation_audit(uvlas_test_structure_pg_integration_partition_resources)
uvlas_require_audited_dependencies(uvlas_test_structure_pg_integration_partition_resources)

uvlas_add_test(structure.pg_rank_preparation
    LABELS property negative verification V-PG
    TIMEOUT 120
    SOURCES structure/pg_rank_preparation_test.cpp
    LIBS uvlas_structure Threads::Threads)
uvlas_add_test(structure.pg_rank_preparation_resources
    LABELS property negative verification V-PG
    SOURCES structure/pg_rank_preparation_resource_test.cpp
    LIBS uvlas_structure_audit)
uvlas_enable_allocation_audit(uvlas_test_structure_pg_rank_preparation_resources)
uvlas_require_audited_dependencies(uvlas_test_structure_pg_rank_preparation_resources)

uvlas_add_test(structure.pg_reference_geometry
    LABELS unit property negative verification V-PG
    SOURCES structure/pg_reference_geometry_test.cpp
    LIBS uvlas_structure Boost::multiprecision)
target_include_directories(uvlas_test_structure_pg_reference_geometry
    PRIVATE "${PROJECT_SOURCE_DIR}/src/structure")

uvlas_add_test(aero.p1_induction
    LABELS unit property negative verification V-AERO
    SOURCES aero/p1_induction_test.cpp aero/p1_area_reference_test.cpp
    LIBS uvlas_aero Boost::multiprecision Threads::Threads nlohmann_json::nlohmann_json)
target_include_directories(uvlas_test_aero_p1_induction PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_compile_definitions(uvlas_test_aero_p1_induction PRIVATE
    UVLAS_P1_GAUSS_FIXTURE="${PROJECT_SOURCE_DIR}/tests/aero/fixtures/p1-gauss-legendre-v1.json"
    UVLAS_P1_AREA_FIXTURE="${PROJECT_SOURCE_DIR}/tests/aero/fixtures/p1-area-reference-v1.json"
    UVLAS_P1_AREA_PRECISE_FIXTURE="${PROJECT_SOURCE_DIR}/tests/aero/fixtures/p1-area-precision-v1.json")
if(UVLAS_SANITIZERS AND NOT MSVC)
    # The full 13-column interval/independent-reference suite exceeded 60 s under
    # ASan/UBSan in PR129. This runner limit does not alter numerical or performance acceptance.
    set_tests_properties(aero.p1_induction PROPERTIES TIMEOUT 300)
endif()
uvlas_add_test(aero.p1_induction_resources
    LABELS property negative verification V-AERO
    SOURCES aero/p1_induction_resources_test.cpp
    LIBS uvlas_aero_audit Boost::headers)
target_include_directories(uvlas_test_aero_p1_induction_resources PRIVATE "${PROJECT_SOURCE_DIR}/src")
uvlas_enable_allocation_audit(uvlas_test_aero_p1_induction_resources)
uvlas_require_audited_dependencies(uvlas_test_aero_p1_induction_resources)

uvlas_add_test(aero.p1_shared_action
    LABELS unit property negative verification V-AERO
    SOURCES aero/p1_shared_action_test.cpp
    LIBS uvlas_aero Boost::multiprecision nlohmann_json::nlohmann_json)
target_include_directories(uvlas_test_aero_p1_shared_action PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_compile_definitions(uvlas_test_aero_p1_shared_action PRIVATE
    UVLAS_P1_SHARED_FIXTURE="${PROJECT_SOURCE_DIR}/tests/aero/fixtures/p1-shared-action-v1.json")

uvlas_add_test(aero.global_geometry
    LABELS unit property negative verification V-AERO
    SOURCES aero/global_geometry_test.cpp
    LIBS uvlas_aero Boost::multiprecision Threads::Threads nlohmann_json::nlohmann_json)
target_include_directories(uvlas_test_aero_global_geometry PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_compile_definitions(uvlas_test_aero_global_geometry PRIVATE
    UVLAS_GLOBAL_GEOMETRY_FIXTURE="${PROJECT_SOURCE_DIR}/tests/aero/fixtures/global-geometry-reference-v1.json")
uvlas_add_test(aero.global_geometry_resources
    LABELS property negative verification V-AERO
    SOURCES aero/global_geometry_resources_test.cpp
    LIBS uvlas_aero_audit Boost::headers)
target_include_directories(uvlas_test_aero_global_geometry_resources PRIVATE "${PROJECT_SOURCE_DIR}/src")
uvlas_enable_allocation_audit(uvlas_test_aero_global_geometry_resources)
uvlas_require_audited_dependencies(uvlas_test_aero_global_geometry_resources)
uvlas_add_test(aero.pchip_rational_slots
    LABELS "V-AERO;verification;negative"
    SOURCES aero/pchip_rational_slots_test.cpp
    LIBS uvlas_aero Boost::multiprecision Threads::Threads)
target_include_directories(uvlas_test_aero_pchip_rational_slots PRIVATE "${PROJECT_SOURCE_DIR}/src")
uvlas_add_test(aero.pchip_rational_slots_resources
    LABELS property negative verification V-AERO
    SOURCES aero/pchip_rational_slots_resources_test.cpp
    LIBS uvlas_aero_audit)
target_include_directories(uvlas_test_aero_pchip_rational_slots_resources PRIVATE "${PROJECT_SOURCE_DIR}/src")
uvlas_enable_allocation_audit(uvlas_test_aero_pchip_rational_slots_resources)
uvlas_require_audited_dependencies(uvlas_test_aero_pchip_rational_slots_resources)
uvlas_add_test(aero.pchip_coefficients
    LABELS V-AERO property negative verification
    SOURCES aero/pchip_coefficients_test.cpp
    LIBS uvlas_aero Boost::multiprecision Threads::Threads nlohmann_json::nlohmann_json)
target_include_directories(uvlas_test_aero_pchip_coefficients PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_compile_definitions(uvlas_test_aero_pchip_coefficients PRIVATE
    UVLAS_PCHIP_COEFFICIENT_REFERENCE="${PROJECT_SOURCE_DIR}/tests/aero/fixtures/pchip-coefficients-v1.json")
uvlas_add_test(aero.pchip_coefficients_resources
    LABELS V-AERO property negative verification
    SOURCES aero/pchip_coefficients_resources_test.cpp
    LIBS uvlas_aero_audit)
target_include_directories(uvlas_test_aero_pchip_coefficients_resources PRIVATE "${PROJECT_SOURCE_DIR}/src")
uvlas_enable_allocation_audit(uvlas_test_aero_pchip_coefficients_resources)
uvlas_require_audited_dependencies(uvlas_test_aero_pchip_coefficients_resources)
uvlas_add_test(aero.pchip_interval
    LABELS V-AERO property negative verification
    SOURCES aero/pchip_interval_jet_test.cpp
    LIBS uvlas_aero Boost::multiprecision Threads::Threads nlohmann_json::nlohmann_json)
target_include_directories(uvlas_test_aero_pchip_interval PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_compile_definitions(uvlas_test_aero_pchip_interval PRIVATE
    UVLAS_PCHIP_COEFFICIENT_REFERENCE="${PROJECT_SOURCE_DIR}/tests/aero/fixtures/pchip-coefficients-v1.json"
    UVLAS_PCHIP_INTERVAL_REFERENCE="${PROJECT_SOURCE_DIR}/tests/aero/fixtures/pchip-interval-values-v1.json")
uvlas_add_test(aero.pchip_locator
    LABELS V-AERO property negative verification
    SOURCES aero/pchip_knot_locator_test.cpp
    LIBS uvlas_aero Boost::multiprecision Threads::Threads nlohmann_json::nlohmann_json)
target_include_directories(uvlas_test_aero_pchip_locator PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_compile_definitions(uvlas_test_aero_pchip_locator PRIVATE
    UVLAS_PCHIP_COEFFICIENT_REFERENCE="${PROJECT_SOURCE_DIR}/tests/aero/fixtures/pchip-coefficients-v1.json"
    UVLAS_PCHIP_INTERVAL_REFERENCE="${PROJECT_SOURCE_DIR}/tests/aero/fixtures/pchip-interval-values-v1.json")

uvlas_add_test(aero.pchip_domain
    LABELS V-AERO property negative verification
    SOURCES aero/pchip_domain_route_test.cpp
    LIBS uvlas_aero Boost::multiprecision Threads::Threads nlohmann_json::nlohmann_json)
target_include_directories(uvlas_test_aero_pchip_domain PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_compile_definitions(uvlas_test_aero_pchip_domain PRIVATE
    UVLAS_PCHIP_COEFFICIENT_REFERENCE="${PROJECT_SOURCE_DIR}/tests/aero/fixtures/pchip-coefficients-v1.json")

uvlas_add_test(aero.pchip_boundary
    LABELS V-AERO property negative verification
    SOURCES aero/pchip_boundary_test.cpp
    LIBS uvlas_aero Boost::multiprecision Threads::Threads)
target_include_directories(uvlas_test_aero_pchip_boundary PRIVATE "${PROJECT_SOURCE_DIR}/src")

uvlas_add_test(aero.pchip_side
    LABELS V-AERO property negative verification
    SOURCES aero/pchip_side_test.cpp support/allocation_audit.cpp
    LIBS uvlas_aero_audit uvlas::eigen Boost::headers Threads::Threads)
target_include_directories(uvlas_test_aero_pchip_side PRIVATE "${PROJECT_SOURCE_DIR}/src")
uvlas_enable_allocation_audit(uvlas_test_aero_pchip_side)
uvlas_require_audited_dependencies(uvlas_test_aero_pchip_side)


uvlas_add_test(aero.pchip_preparation_input
    LABELS V-AERO property negative verification
    SOURCES aero/pchip_preparation_input_test.cpp
    LIBS uvlas_aero Boost::multiprecision Threads::Threads)
target_include_directories(uvlas_test_aero_pchip_preparation_input PRIVATE "${PROJECT_SOURCE_DIR}/src")

uvlas_add_test(aero.pchip_curve
    LABELS V-AERO property negative verification
    SOURCES aero/pchip_curve_test.cpp
    LIBS uvlas_aero Boost::multiprecision nlohmann_json::nlohmann_json Threads::Threads)
target_include_directories(uvlas_test_aero_pchip_curve PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_compile_definitions(uvlas_test_aero_pchip_curve PRIVATE
    UVLAS_PCHIP_COEFFICIENT_REFERENCE="${PROJECT_SOURCE_DIR}/tests/aero/fixtures/pchip-coefficients-v1.json"
    UVLAS_PCHIP_INTERVAL_REFERENCE="${PROJECT_SOURCE_DIR}/tests/aero/fixtures/pchip-interval-values-v1.json")

uvlas_add_test(io.case_schema
    LABELS unit property negative verification
    SOURCES io/case_schema_test.cpp io/case_schema_intrinsic_test.cpp io/case_schema_cross_test.cpp io/case_schema_fragments_test.cpp io/case_schema_inventory_test.cpp
    LIBS uvlas_io nlohmann_json::nlohmann_json)
target_include_directories(uvlas_test_io_case_schema PRIVATE "${PROJECT_SOURCE_DIR}/src/io")
uvlas_add_test(io.case_schema_resources
    LABELS property negative verification
    SOURCES io/case_schema_resource_test.cpp
    LIBS uvlas_io_audit Threads::Threads)
target_include_directories(uvlas_test_io_case_schema_resources PRIVATE "${PROJECT_SOURCE_DIR}/src/io")
uvlas_enable_allocation_audit(uvlas_test_io_case_schema_resources)
uvlas_require_audited_dependencies(uvlas_test_io_case_schema_resources)
if(UVLAS_SANITIZERS AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
    # Both full schema suites reached the default 60 s on GCC 14 with
    # ASan/UBSan. Give these runners a separate, bounded test time.
    set_tests_properties(io.case_schema io.case_schema_resources PROPERTIES TIMEOUT 300)
endif()

uvlas_add_test(io.reference_root_plan
    LABELS integration property negative verification
    SOURCES io/reference_root_plan_test.cpp
    LIBS uvlas_io nlohmann_json::nlohmann_json)
target_include_directories(uvlas_test_io_reference_root_plan PRIVATE "${PROJECT_SOURCE_DIR}/src/io")
uvlas_add_test(io.reference_dependency_walk
    LABELS unit property negative verification
    SOURCES io/reference_dependency_walk_test.cpp
    LIBS uvlas_io)
target_include_directories(uvlas_test_io_reference_dependency_walk
    PRIVATE "${PROJECT_SOURCE_DIR}/src/io")
uvlas_add_test(io.reference_declarations
    LABELS unit property negative verification
    SOURCES io/reference_declarations_test.cpp
    LIBS uvlas_io nlohmann_json::nlohmann_json)
target_include_directories(uvlas_test_io_reference_declarations PRIVATE "${PROJECT_SOURCE_DIR}/src/io")
uvlas_add_test(io.reference_declarations_resources
    LABELS property negative verification
    SOURCES io/reference_declarations_resource_test.cpp
    LIBS uvlas_io_audit nlohmann_json::nlohmann_json Threads::Threads)
target_include_directories(uvlas_test_io_reference_declarations_resources PRIVATE "${PROJECT_SOURCE_DIR}/src/io")
uvlas_enable_allocation_audit(uvlas_test_io_reference_declarations_resources)
uvlas_require_audited_dependencies(uvlas_test_io_reference_declarations_resources)

uvlas_add_test(io.si_interval
    LABELS unit property negative verification
    SOURCES io/si_case_interval_test.cpp
    LIBS uvlas_core Boost::headers)
target_include_directories(uvlas_test_io_si_interval PRIVATE "${PROJECT_SOURCE_DIR}/src/io")

uvlas_add_test(io.si_case
    LABELS unit property negative verification
    SOURCES io/si_case_units_test.cpp io/si_case_ownership_test.cpp io/si_case_fragments_test.cpp io/si_case_fields_test.cpp
    LIBS uvlas_io Boost::headers)
target_include_directories(uvlas_test_io_si_case PRIVATE "${PROJECT_SOURCE_DIR}/src/io")

uvlas_add_test(io.si_case_resources
    LABELS property negative verification
    SOURCES io/si_case_resource_test.cpp
    LIBS uvlas_io_audit Threads::Threads Boost::headers)
target_include_directories(uvlas_test_io_si_case_resources PRIVATE "${PROJECT_SOURCE_DIR}/src/io")
uvlas_enable_allocation_audit(uvlas_test_io_si_case_resources)
uvlas_require_audited_dependencies(uvlas_test_io_si_case_resources)
