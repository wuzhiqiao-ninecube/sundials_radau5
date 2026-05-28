cmake_minimum_required(VERSION 3.18)
project(radau5 C)

# Find SUNDIALS (user must set SUNDIALS_DIR or CMAKE_PREFIX_PATH)
find_package(SUNDIALS REQUIRED)

# Find LAPACK (required for complex E2 solver: zgetrf/zgetrs)
find_package(LAPACK REQUIRED)

# Find KLU directly (for klu_z_* complex functions)
find_package(PkgConfig QUIET)
if(TARGET SUNDIALS::sunlinsolklu OR TARGET SUNDIALS::sundials_sunlinsolklu)
  # KLU is available via SUNDIALS — find the headers and library directly
  find_path(KLU_INCLUDE_DIR klu.h
    HINTS /usr/include/suitesparse /usr/local/include/suitesparse
          ${KLU_DIR}/include ${SUITESPARSE_DIR}/include)
  find_library(KLU_LIBRARY NAMES klu
    HINTS /usr/lib /usr/local/lib /usr/lib/x86_64-linux-gnu
          ${KLU_DIR}/lib ${SUITESPARSE_DIR}/lib)
  set(RADAU5_HAVE_KLU ON)
endif()

# Source files
set(RADAU5_SOURCES
    src/radau5.c
    src/radau5_linsys.c
    src/radau5_newt.c
    src/radau5_estrad.c
    src/radau5_step.c
    src/radau5_contr.c
    src/radau5_ic.c
    src/radau5_colgroup.c
    src/radau5_root.c
)

# Library
add_library(radau5 ${RADAU5_SOURCES})
target_include_directories(radau5
    PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(radau5
    PUBLIC
        SUNDIALS::nvecserial
        SUNDIALS::sunmatrixdense
        SUNDIALS::sunlinsoldense
        SUNDIALS::core
    PRIVATE
        ${LAPACK_LIBRARIES}
)

# Optional KLU support
if(RADAU5_HAVE_KLU)
    if(TARGET SUNDIALS::sunlinsolklu)
        target_link_libraries(radau5 PUBLIC SUNDIALS::sunlinsolklu)
    elseif(TARGET SUNDIALS::sundials_sunlinsolklu)
        target_link_libraries(radau5 PUBLIC SUNDIALS::sundials_sunlinsolklu)
    endif()
    # Direct KLU link for klu_z_* complex functions
    if(KLU_INCLUDE_DIR)
        target_include_directories(radau5 PRIVATE ${KLU_INCLUDE_DIR})
    endif()
    if(KLU_LIBRARY)
        target_link_libraries(radau5 PRIVATE ${KLU_LIBRARY})
    endif()
    target_compile_definitions(radau5 PRIVATE RADAU5_HAVE_KLU)
endif()

# Optional band support
target_link_libraries(radau5 PUBLIC
    SUNDIALS::sunmatrixband
    SUNDIALS::sunlinsolband
)

# Math library
target_link_libraries(radau5 PRIVATE m)

# Tests
option(RADAU5_BUILD_TESTS "Build radau5 tests" ON)
if(RADAU5_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Examples
option(RADAU5_BUILD_EXAMPLES "Build radau5 examples" ON)
if(RADAU5_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
