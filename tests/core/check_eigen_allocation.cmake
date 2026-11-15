execute_process(COMMAND "${PROBE}" clean RESULT_VARIABLE clean OUTPUT_VARIABLE out ERROR_VARIABLE err TIMEOUT 10)
if(NOT "${clean}" STREQUAL "0")
    message(FATAL_ERROR "Eigen allocation positive control failed: ${out} ${err}")
endif()
foreach(mode IN ITEMS allocate reallocate)
    execute_process(COMMAND "${PROBE}" ${mode} RESULT_VARIABLE result
        OUTPUT_VARIABLE out ERROR_VARIABLE err TIMEOUT 10)
    if("${result}" STREQUAL "0" OR "${result}" STREQUAL "4" OR
       NOT "${out}${err}" MATCHES "heap allocation is forbidden")
        message(FATAL_ERROR "Eigen ${mode} was not rejected by the expected guard: ${result} ${out} ${err}")
    endif()
endforeach()
