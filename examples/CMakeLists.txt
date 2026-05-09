# Helper function to add a radau5 example
function(radau5_add_example example_name)
    add_executable(${example_name} ${example_name}.c)
    target_link_libraries(${example_name} PRIVATE radau5 m)
endfunction()

# Helper for examples that need the sparse matrix module
function(radau5_add_sparse_example example_name)
    add_executable(${example_name} ${example_name}.c)
    target_link_libraries(${example_name} PRIVATE
        radau5
        SUNDIALS::sunmatrixsparse
        m
    )
endfunction()

# All examples
radau5_add_example(radau5_vdpol)
radau5_add_example(radau5_rober)
radau5_add_example(radau5_orego)
radau5_add_example(radau5_e5)
radau5_add_example(radau5_hires)
radau5_add_example(radau5_plei)
radau5_add_example(radau5_ringmod)
radau5_add_example(radau5_heat1d)
radau5_add_example(radau5_medakzo)
radau5_add_sparse_example(radau5_pollu)
radau5_add_example(radau5_chemakzo)
radau5_add_example(radau5_transamp)
radau5_add_sparse_example(radau5_transamp_sparse)
radau5_add_example(radau5_andrews)
radau5_add_example(radau5_ark_analytic)
radau5_add_example(radau5_vdpolm)
radau5_add_example(radau5_caraxis)
radau5_add_example(radau5_beam)
radau5_add_example(radau5_crank)
radau5_add_example(radau5_water)
radau5_add_example(radau5_fekete)
radau5_add_example(radau5_tba)
radau5_add_example(radau5_bounce)
radau5_add_example(radau5_roberts_root)
radau5_add_example(radau5_rocket)
radau5_add_example(radau5_orbit)
radau5_add_example(radau5_kepler)
radau5_add_example(radau5_kneeode)

# ---------------------------------------------------------------------------
# CTest: each problem with default args, eigenvalue mode and Schur mode
# ---------------------------------------------------------------------------
function(radau5_add_tests prob rt at h)
    add_test(NAME radau5_${prob}_eigen
             COMMAND radau5_${prob} ${rt} ${at} ${h} 0)
    add_test(NAME radau5_${prob}_schur
             COMMAND radau5_${prob} ${rt} ${at} ${h} 1)
endfunction()

radau5_add_tests(vdpol    1e-6   1e-6   1e-6)
radau5_add_tests(rober    1e-10  1e-14  1e-12)
radau5_add_tests(orego    1e-6   1e-6   1e-6)
radau5_add_tests(e5       1e-4   1e-24  1e-6)
radau5_add_tests(hires    1e-6   1e-6   1e-6)
radau5_add_tests(plei     1e-6   1e-6   1e-6)
radau5_add_tests(ringmod  1e-7   1e-7   1e-6)
radau5_add_tests(heat1d   1e-10  1e-12  1e-6)
radau5_add_tests(medakzo  1e-7   1e-7   1e-6)
radau5_add_tests(pollu    1e-6   1e-6   1e-6)
radau5_add_tests(chemakzo 1e-7   1e-7   1e-10)
radau5_add_tests(transamp 1e-6   1e-6   1e-6)
radau5_add_tests(transamp_sparse 1e-6 1e-6 1e-6)
radau5_add_tests(andrews  1e-6   1e-6   1e-10)
radau5_add_tests(vdpolm   1e-6   1e-6   1e-6)
radau5_add_tests(caraxis  1e-6   1e-6   1e-10)
radau5_add_tests(beam     1e-6   1e-6   1e-6)
radau5_add_tests(crank    1e-6   1e-6   1e-10)
radau5_add_tests(water    1e-7   1e-7   1e-6)
radau5_add_tests(fekete   1e-6   1e-6   1e-6)
radau5_add_tests(tba      1e-5   1e-5   4e-5)

# Rootfinding examples (no Schur variant needed — just run once)
add_test(NAME radau5_bounce COMMAND radau5_bounce)
add_test(NAME radau5_roberts_root COMMAND radau5_roberts_root)
add_test(NAME radau5_rocket COMMAND radau5_rocket)
add_test(NAME radau5_orbit COMMAND radau5_orbit)
add_test(NAME radau5_kepler COMMAND radau5_kepler)
add_test(NAME radau5_kneeode COMMAND radau5_kneeode)
