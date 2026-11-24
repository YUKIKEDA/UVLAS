find_package(Python3 3.12 EXACT REQUIRED COMPONENTS Interpreter Development.Module NumPy)
if(Python3_NumPy_VERSION VERSION_LESS 2 OR NOT Python3_NumPy_VERSION VERSION_LESS 3)
    message(FATAL_ERROR "The validated Python binding requires NumPy 2.x")
endif()
set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 3.1.0 EXACT CONFIG REQUIRED)

pybind11_add_module(uvlas_python MODULE NO_EXTRAS module.cpp)
uvlas_configure_target(uvlas_python)
target_link_libraries(uvlas_python PRIVATE uvlas_core uvlas_numerics)
set_target_properties(uvlas_python PROPERTIES OUTPUT_NAME _uvlas)
if(NOT CMAKE_CROSSCOMPILING)
    uvlas_record_runtime_link(uvlas_python)
endif()

if(BUILD_TESTING)
    add_test(NAME python.ownership
        COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/python/test_ownership.py"
            "$<TARGET_FILE_DIR:uvlas_python>")
    set_tests_properties(python.ownership PROPERTIES LABELS "integration;negative" TIMEOUT 90)
    add_test(NAME python.case_file_operations
        COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/python/test_case_file_operations.py"
            "$<TARGET_FILE_DIR:uvlas_python>" "$<TARGET_FILE:uvlas_cli>"
            "${PROJECT_SOURCE_DIR}/docs/specs/foundation/case-schema-examples.json")
    set_tests_properties(python.case_file_operations PROPERTIES
        LABELS "integration;negative;verification" TIMEOUT 30)
endif()
