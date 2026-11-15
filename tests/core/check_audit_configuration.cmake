foreach(case IN ITEMS good bad hidden object-good object-bad raw-source raw-library source-expression
                      imported-good imported-bad imported-source-good imported-source-bad directory-source)
    execute_process(COMMAND "${CMAKE_COMMAND}" -S "${SOURCE}" -B "${BINARY}/${case}"
        "-DCASE=${case}" "-DAUDIT_MODULE=${AUDIT_MODULE}" "-DCMAKE_CXX_COMPILER=${COMPILER}"
        RESULT_VARIABLE code OUTPUT_VARIABLE out ERROR_VARIABLE err TIMEOUT 30)
    if(case MATCHES "(^good$|-good$)")
        if(NOT "${code}" STREQUAL "0")
            message(FATAL_ERROR "Positive audit configuration failed: ${out} ${err}")
        endif()
    elseif("${code}" STREQUAL "0" OR NOT "${out}${err}" MATCHES
           "unaudited product target product|opaque compiled artifact|explicit source dependencies")
        message(FATAL_ERROR "Unaudited dependency ${case} was not rejected: ${out} ${err}")
    endif()
endforeach()
