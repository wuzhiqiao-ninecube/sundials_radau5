# Helper function to add a radau5 test
function(radau5_add_test test_name)
    add_executable(${test_name} ${test_name}.c)
    target_link_libraries(${test_name} PRIVATE radau5 m)
    target_include_directories(${test_name} PRIVATE
        ${CMAKE_SOURCE_DIR}/src
    )
    add_test(NAME ${test_name} COMMAND ${test_name})
endfunction()

# Unit tests
radau5_add_test(test_radau5_constants)
radau5_add_test(test_radau5_api)
radau5_add_test(test_radau5_build_e1)
radau5_add_test(test_radau5_build_e2)
radau5_add_test(test_radau5_dq_jac)
radau5_add_test(test_radau5_dahlquist)
radau5_add_test(test_radau5_colgroup)
radau5_add_test(test_radau5_dq_jac_sparse)
