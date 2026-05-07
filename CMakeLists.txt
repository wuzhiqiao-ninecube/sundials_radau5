cmake_minimum_required(VERSION 3.18)
project(radau5 C)

# Find SUNDIALS (user must set SUNDIALS_DIR or CMAKE_PREFIX_PATH)
find_package(SUNDIALS REQUIRED)

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
)

# Optional KLU support
if(TARGET SUNDIALS::sunlinsolklu)
    target_link_libraries(radau5 PUBLIC SUNDIALS::sunlinsolklu)
    target_compile_definitions(radau5 PRIVATE RADAU5_HAVE_KLU)
elseif(TARGET SUNDIALS::sundials_sunlinsolklu)
    target_link_libraries(radau5 PUBLIC SUNDIALS::sunlinsolklu)
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
