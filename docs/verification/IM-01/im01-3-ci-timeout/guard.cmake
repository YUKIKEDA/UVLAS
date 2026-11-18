if(UVLAS_SANITIZERS AND NOT MSVC)
    set_tests_properties(structure.pg_coefficients PROPERTIES TIMEOUT 900)
endif()
