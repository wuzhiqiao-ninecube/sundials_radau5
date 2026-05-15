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
radau5_add_example(radau5_tba2)
radau5_add_example(radau5_bounce)
radau5_add_example(radau5_roberts_root)
radau5_add_example(radau5_rocket)
radau5_add_example(radau5_orbit)
radau5_add_example(radau5_kepler)
radau5_add_example(radau5_kneeode)

# ---------------------------------------------------------------------------
# CTest: each problem × 2 transform modes × 4 order configurations
#
# Order configurations:
#   ns3   nsmin=3 nsmax=3  fixed order 5  (classic RADAU5)
#   ns5   nsmin=5 nsmax=5  fixed order 9
#   ns7   nsmin=7 nsmax=7  fixed order 13
#   nsvar nsmin=3 nsmax=7  variable order
# ---------------------------------------------------------------------------

# Standard problems: args are rtol atol h0 use_schur nsmin nsmax
function(radau5_add_tests prob rt at h)
    set(schur_tags eigen schur)
    set(schur_vals 0 1)
    set(ns_tags ns3 ns5 ns7 nsvar)
    set(ns_mins 3 5 7 3)
    set(ns_maxs 3 5 7 7)
    foreach(schur_tag schur_val IN ZIP_LISTS schur_tags schur_vals)
        foreach(ns_tag ns_min ns_max IN ZIP_LISTS ns_tags ns_mins ns_maxs)
            add_test(NAME radau5_${prob}_${schur_tag}_${ns_tag}
                     COMMAND radau5_${prob} ${rt} ${at} ${h} ${schur_val}
                             ${ns_min} ${ns_max})
            set_tests_properties(radau5_${prob}_${schur_tag}_${ns_tag}
                                 PROPERTIES TIMEOUT 600)
        endforeach()
    endforeach()
endfunction()

# DAE-3 problems: ns=7 is known to fail, only test ns3/ns5/nsvar
function(radau5_add_dae3_tests prob rt at h)
    set(schur_tags eigen schur)
    set(schur_vals 0 1)
    set(ns_tags ns3 ns5 nsvar)
    set(ns_mins 3 5 3)
    set(ns_maxs 3 5 7)
    foreach(schur_tag schur_val IN ZIP_LISTS schur_tags schur_vals)
        foreach(ns_tag ns_min ns_max IN ZIP_LISTS ns_tags ns_mins ns_maxs)
            add_test(NAME radau5_${prob}_${schur_tag}_${ns_tag}
                     COMMAND radau5_${prob} ${rt} ${at} ${h} ${schur_val}
                             ${ns_min} ${ns_max})
            set_tests_properties(radau5_${prob}_${schur_tag}_${ns_tag}
                                 PROPERTIES TIMEOUT 600)
        endforeach()
    endforeach()
endfunction()

radau5_add_tests(vdpol    1e-6   1e-6   1e-6)
radau5_add_tests(rober    1e-10  1e-14  1e-12)
radau5_add_tests(orego    1e-6   1e-6   1e-6)
radau5_add_tests(e5       1e-4   1e-24  1e-6)
radau5_add_tests(hires    1e-9   1e-9   1e-9)
radau5_add_tests(plei     1e-6   1e-6   1e-6)
radau5_add_tests(ringmod  1e-7   1e-7   1e-6)
radau5_add_tests(heat1d   1e-10  1e-12  1e-6)
radau5_add_tests(medakzo  1e-7   1e-7   1e-6)
radau5_add_tests(chemakzo 1e-7   1e-7   1e-10)
radau5_add_tests(transamp 1e-6   1e-6   1e-6)
radau5_add_tests(transamp_sparse 1e-6 1e-6 1e-6)
radau5_add_dae3_tests(andrews  1e-6   1e-6   1e-10)
radau5_add_tests(vdpolm   1e-6   1e-6   1e-6)
radau5_add_dae3_tests(caraxis  1e-6   1e-6   1e-10)
radau5_add_tests(beam     1e-6   1e-6   1e-6)
radau5_add_tests(crank    1e-6   1e-6   1e-10)
radau5_add_tests(water    1e-7   1e-7   1e-6)
radau5_add_tests(fekete   1e-6   1e-6   1e-6)
radau5_add_dae3_tests(tba      1e-5   1e-5   4e-5)
radau5_add_tests(ark_analytic 1e-5 1e-10 1e-4)

# Pollu has extra use_dq arg: rtol atol h0 use_schur use_dq nsmin nsmax
function(radau5_add_pollu_tests)
    set(schur_tags eigen schur)
    set(schur_vals 0 1)
    set(ns_tags ns3 ns5 ns7 nsvar)
    set(ns_mins 3 5 7 3)
    set(ns_maxs 3 5 7 7)
    foreach(schur_tag schur_val IN ZIP_LISTS schur_tags schur_vals)
        foreach(ns_tag ns_min ns_max IN ZIP_LISTS ns_tags ns_mins ns_maxs)
            add_test(NAME radau5_pollu_${schur_tag}_${ns_tag}
                     COMMAND radau5_pollu 1e-6 1e-6 1e-6 ${schur_val} 0
                             ${ns_min} ${ns_max})
            set_tests_properties(radau5_pollu_${schur_tag}_${ns_tag}
                                 PROPERTIES TIMEOUT 600)
        endforeach()
    endforeach()
endfunction()
radau5_add_pollu_tests()

# Rootfinding examples: args are nsmin nsmax (no rtol/atol/h0/schur)
function(radau5_add_root_tests prob)
    set(ns_tags ns3 ns5 ns7 nsvar)
    set(ns_mins 3 5 7 3)
    set(ns_maxs 3 5 7 7)
    foreach(ns_tag ns_min ns_max IN ZIP_LISTS ns_tags ns_mins ns_maxs)
        add_test(NAME radau5_${prob}_${ns_tag}
                 COMMAND radau5_${prob} ${ns_min} ${ns_max})
        set_tests_properties(radau5_${prob}_${ns_tag}
                             PROPERTIES TIMEOUT 600)
    endforeach()
endfunction()

radau5_add_root_tests(bounce)
radau5_add_root_tests(roberts_root)
radau5_add_root_tests(rocket)
radau5_add_root_tests(orbit)
radau5_add_root_tests(kepler)
radau5_add_root_tests(kneeode)
