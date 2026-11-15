# N4-01: exercise fresh configurations of the actual product sources, including
# CTest's default ON, without changing installed dependencies or building C++.
file(REAL_PATH "${BINARY_ROOT}" binary_root)
string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef token)
set(work "${binary_root}/notice-configuration-${token}")
if(EXISTS "${work}")
    message(FATAL_ERROR "Refusing to reuse a notice test directory")
endif()
file(MAKE_DIRECTORY "${work}/source")
file(REAL_PATH "${work}" work)
cmake_path(IS_PREFIX binary_root "${work}" NORMALIZE inside)
if(NOT inside)
    message(FATAL_ERROR "Notice test directory escaped its build directory")
endif()
function(fail_case text)
    file(REMOVE_RECURSE "${work}")
    message(FATAL_ERROR "${text}")
endfunction()

file(COPY "${SOURCE_ROOT}/CMakeLists.txt" "${SOURCE_ROOT}/cmake"
    "${SOURCE_ROOT}/src" "${SOURCE_ROOT}/include" "${SOURCE_ROOT}/tests"
    "${SOURCE_ROOT}/qa" "${SOURCE_ROOT}/benchmarks" DESTINATION "${work}/source")
set(installed "${work}/installed")
foreach(package IN ITEMS eigen3 spectra boost-multiprecision boost-hash2 nlohmann-json catch2)
    configure_file("${DEPENDENCY_PREFIX}/share/${package}/copyright"
        "${installed}/${TRIPLET}/share/${package}/copyright" COPYONLY)
endforeach()

set(count 0)
foreach(mode IN ITEMS default ON OFF)
    foreach(present IN ITEMS OFF ON)
        set(notice "${installed}/${TRIPLET}/share/catch2/copyright")
        if(present)
            configure_file("${DEPENDENCY_PREFIX}/share/catch2/copyright" "${notice}" COPYONLY)
        else()
            file(REMOVE "${notice}")
        endif()
        set(options)
        if(NOT mode STREQUAL "default")
            list(APPEND options "-DBUILD_TESTING:BOOL=${mode}")
        endif()
        set(build "${work}/${mode}-${present}")
        execute_process(COMMAND "${CMAKE_COMMAND}" -S "${work}/source" -B "${build}" -G Ninja
            "-DCMAKE_CXX_COMPILER=${CXX_COMPILER}"
            "-DCMAKE_PREFIX_PATH=${DEPENDENCY_PREFIX}"
            "-DVCPKG_INSTALLED_DIR=${installed}" "-DVCPKG_TARGET_TRIPLET=${TRIPLET}"
            "-DPython3_EXECUTABLE=${PYTHON_EXECUTABLE}" ${options}
            RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error TIMEOUT 60)
        if(NOT mode STREQUAL "OFF" AND NOT present)
            if(result EQUAL 0 OR NOT "${output}\n${error}" MATCHES "Missing third-party notice for catch2")
                fail_case("${mode}/missing Catch2 did not reject the missing notice: ${result}\n${output}\n${error}")
            endif()
        else()
            if(NOT result EQUAL 0)
                fail_case("${mode}/notice=${present} failed: ${result}\n${output}\n${error}")
            endif()
            if(mode STREQUAL "OFF")
                set(expected OFF)
            else()
                set(expected ON)
            endif()
            file(STRINGS "${build}/CMakeCache.txt" actual REGEX "^BUILD_TESTING:BOOL=")
            if(NOT actual STREQUAL "BUILD_TESTING:BOOL=${expected}")
                fail_case("${mode}/notice=${present}: unexpected testing mode ${actual}")
            endif()
        endif()
        math(EXPR count "${count}+1")
    endforeach()
endforeach()

# Only package discovery and the actual registered notice CTest run here. This
# interface fixture never supplies headers or libraries to a compiled target.
set(benchmark_package "${work}/benchmark-package")
file(MAKE_DIRECTORY "${benchmark_package}")
file(WRITE "${benchmark_package}/benchmarkConfig.cmake"
    "add_library(benchmark::benchmark INTERFACE IMPORTED)\n")
set(benchmark_notice "${installed}/${TRIPLET}/share/benchmark/copyright")
file(WRITE "${benchmark_notice}" "Benchmark upstream notice fixture\n")

function(configure_benchmark name mode)
    set(options)
    if(NOT mode STREQUAL "default")
        list(APPEND options "-DUVLAS_BUILD_BENCHMARKS:BOOL=${mode}")
    endif()
    execute_process(COMMAND "${CMAKE_COMMAND}" -S "${work}/source" -B "${work}/${name}" -G Ninja
        "-DCMAKE_CXX_COMPILER=${CXX_COMPILER}"
        "-DCMAKE_PREFIX_PATH=${DEPENDENCY_PREFIX};${benchmark_package}"
        "-DVCPKG_INSTALLED_DIR=${installed}" "-DVCPKG_TARGET_TRIPLET=${TRIPLET}"
        "-DVCPKG_MANIFEST_FEATURES=benchmarks" "-DBUILD_TESTING:BOOL=ON"
        "-DPython3_EXECUTABLE=${PYTHON_EXECUTABLE}" ${options}
        RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error TIMEOUT 60)
    if(NOT result EQUAL 0)
        fail_case("Benchmark ${name}/${mode} configure failed: ${result}\n${output}\n${error}")
    endif()
endfunction()

function(check_registered_notices name expect_success)
    execute_process(COMMAND "${CTEST_COMMAND}" --test-dir "${work}/${name}"
        --output-on-failure -R "^toolchain.third_party_notices$"
        RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error TIMEOUT 30)
    if(expect_success)
        if(NOT result EQUAL 0 OR NOT output MATCHES "1/1 Test")
            fail_case("Benchmark ${name} registered CTest failed: ${result}\n${output}\n${error}")
        endif()
    elseif(result EQUAL 0 OR NOT "${output}\n${error}" MATCHES "Missing copied notice: benchmark")
        fail_case("Benchmark ${name} missing copy was not rejected: ${result}\n${output}\n${error}")
    endif()
endfunction()

foreach(mode IN ITEMS default OFF ON)
    configure_benchmark("benchmark-${mode}" "${mode}")
    set(copy "${work}/benchmark-${mode}/third-party-notices/benchmark.txt")
    if(mode STREQUAL "ON")
        if(NOT EXISTS "${copy}")
            fail_case("Enabled Benchmark notice was not copied")
        endif()
    elseif(EXISTS "${copy}")
        fail_case("Disabled Benchmark notice was copied for ${mode}")
    endif()
    check_registered_notices("benchmark-${mode}" ON)
endforeach()
# The ON path must reject a lost copy, rather than simply skipping Benchmark.
file(REMOVE "${work}/benchmark-ON/third-party-notices/benchmark.txt")
check_registered_notices("benchmark-ON" OFF)
configure_benchmark("benchmark-ON" ON)
configure_benchmark("benchmark-ON" OFF)
if(EXISTS "${work}/benchmark-ON/third-party-notices/benchmark.txt")
    fail_case("ON to OFF retained the Benchmark notice")
endif()
check_registered_notices("benchmark-ON" ON)

# Required upstream notices must still be enforced when the feature is enabled.
file(REMOVE "${benchmark_notice}")
execute_process(COMMAND "${CMAKE_COMMAND}" -S "${work}/source"
    -B "${work}/benchmark-missing" -G Ninja "-DCMAKE_CXX_COMPILER=${CXX_COMPILER}"
    "-DCMAKE_PREFIX_PATH=${DEPENDENCY_PREFIX};${benchmark_package}"
    "-DVCPKG_INSTALLED_DIR=${installed}" "-DVCPKG_TARGET_TRIPLET=${TRIPLET}"
    "-DVCPKG_MANIFEST_FEATURES=benchmarks" "-DBUILD_TESTING:BOOL=ON"
    "-DUVLAS_BUILD_BENCHMARKS:BOOL=ON" "-DPython3_EXECUTABLE=${PYTHON_EXECUTABLE}"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error TIMEOUT 60)
if(result EQUAL 0 OR NOT "${output}\n${error}" MATCHES "Missing third-party notice for benchmark")
    fail_case("Enabled Benchmark accepted a missing upstream notice: ${result}\n${output}\n${error}")
endif()
file(REMOVE_RECURSE "${work}")
message(STATUS "${count} Catch2 configurations and Benchmark default/OFF/ON/toggle/negative registered CTest checks passed (N4-01)")
