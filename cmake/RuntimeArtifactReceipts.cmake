# SPDX-License-Identifier: Apache-2.0

set(uvlas_runtime_receipt_dir "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles")
set(uvlas_runtime_revision_header "${CMAKE_CURRENT_BINARY_DIR}/generated/uvlas/git_revision.hpp")

function(uvlas_record_runtime_link target)
    set(receipt "${uvlas_runtime_receipt_dir}/${target}-$<CONFIG>.sha256")
    add_custom_command(TARGET ${target} POST_BUILD
        COMMAND "${CMAKE_COMMAND}"
            "-DARTIFACT=$<TARGET_FILE:${target}>"
            "-DREVISION_HEADER=${uvlas_runtime_revision_header}"
            "-DOUTPUT=${receipt}"
            -P "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/RecordRuntimeArtifactLink.cmake"
        BYPRODUCTS "${receipt}"
        VERBATIM)
endfunction()
