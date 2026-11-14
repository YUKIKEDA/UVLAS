# SPDX-License-Identifier: Apache-2.0
#
# Writes OUTPUT with the commit of SOURCE_DIR and whether tracked files differ
# from it (QA-06). Run at build time, so the header follows the checkout; the
# file is rewritten only when its content changes.
execute_process(COMMAND git -C "${SOURCE_DIR}" rev-parse HEAD
    OUTPUT_VARIABLE commit OUTPUT_STRIP_TRAILING_WHITESPACE RESULT_VARIABLE failed ERROR_QUIET)
if(failed)
    set(commit "unknown")
    set(dirty 1)
else()
    execute_process(COMMAND git -C "${SOURCE_DIR}" status --porcelain --untracked-files=no
        OUTPUT_VARIABLE changes OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(changes STREQUAL "")
        set(dirty 0)
    else()
        set(dirty 1)
    endif()
endif()
set(content "#pragma once\n\n#define UVLAS_GIT_COMMIT \"${commit}\"\n#define UVLAS_GIT_DIRTY ${dirty}\n")
if(EXISTS "${OUTPUT}")
    file(READ "${OUTPUT}" previous)
else()
    set(previous "")
endif()
if(NOT previous STREQUAL content)
    file(WRITE "${OUTPUT}" "${content}")
endif()
