# SPDX-License-Identifier: Apache-2.0

option(UVLAS_WARNINGS_AS_ERRORS "Treat compiler warnings in UVLAS targets as errors" OFF)
set(UVLAS_SANITIZERS "" CACHE STRING
    "Sanitizers for UVLAS targets and tests on GCC/Clang, for example address,undefined")

# The floating-point model, reported by build_info (QA-06).
if(MSVC)
    set(UVLAS_FLOATING_POINT "/fp:strict")
else()
    set(UVLAS_FLOATING_POINT "-ffp-contract=off -fno-fast-math")
endif()

# uvlas_configure_target(<target>)
# Warning policy and floating-point model shared by every UVLAS target.
function(uvlas_configure_target target)
    target_compile_features(${target} PUBLIC cxx_std_23)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /permissive- /utf-8 /Zc:__cplusplus)
        # The error-free transformations and compensated residuals of N4-06
        # need IEEE evaluation without contraction or value-changing
        # optimization. Same flags as the design-work-4 numerical spike.
        target_compile_options(${target} PRIVATE /fp:strict)
        target_compile_definitions(${target} PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
        if(UVLAS_WARNINGS_AS_ERRORS)
            target_compile_options(${target} PRIVATE /WX)
        endif()
    else()
        target_compile_options(${target} PRIVATE
            -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Wshadow)
        target_compile_options(${target} PRIVATE -ffp-contract=off -fno-fast-math)
        if(UVLAS_WARNINGS_AS_ERRORS)
            target_compile_options(${target} PRIVATE -Werror)
        endif()
        if(UVLAS_SANITIZERS)
            target_compile_options(${target} PRIVATE
                -fsanitize=${UVLAS_SANITIZERS} -fno-omit-frame-pointer -fno-sanitize-recover=all)
            target_link_options(${target} PRIVATE -fsanitize=${UVLAS_SANITIZERS})
        endif()
    endif()
endfunction()
