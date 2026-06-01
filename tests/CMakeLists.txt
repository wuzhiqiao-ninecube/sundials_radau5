# Helper function to add a radau5 test
function(radau5_add_test test_name)
    add_executable(${test_name} ${test_name}.c)
    target_link_libraries(${test_name} PRIVATE radau5 m)
    target_include_directories(${test_name} PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/nvector/complex_serial
        ${CMAKE_SOURCE_DIR}/src/sunmatrix/complex_dense
        ${CMAKE_SOURCE_DIR}/src/sunmatrix/complex_band
        ${CMAKE_SOURCE_DIR}/src/sunmatrix/complex_sparse
        ${CMAKE_SOURCE_DIR}/src/sunlinsol/complex_dense
        ${CMAKE_SOURCE_DIR}/src/sunlinsol/complex_band
        ${CMAKE_SOURCE_DIR}/src/sunlinsol/complex_sparse
    )
    add_test(NAME ${test_name} COMMAND ${test_name})
endfunction()

# Unit tests
radau5_add_test(test_radau5_constants)
radau5_add_test(test_radau5_constants_high)
radau5_add_test(test_radau5_api)
radau5_add_test(test_radau5_build_e1)
radau5_add_test(test_radau5_build_e2)
radau5_add_test(test_radau5_dq_jac)
radau5_add_test(test_radau5_colgroup)
radau5_add_test(test_radau5_dq_jac_sparse)

# Complex module tests
radau5_add_test(test_nvector_complex_serial)
radau5_add_test(test_sunmatrix_complex_dense)
radau5_add_test(test_sunmatrix_complex_band)
radau5_add_test(test_sunmatrix_complex_sparse)
