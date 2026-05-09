#!/bin/bash
set -e
SUNDIALS_SRC=/mnt/d/workon/sundials
SUNDIALS_BUILD=/mnt/d/workon/sundials/build
RADAU5=/mnt/d/workon/sundials/thirdpart/radau5
CC=clang

INC="-I${SUNDIALS_BUILD}/include -I${SUNDIALS_SRC}/include -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/include/suitesparse -I${RADAU5}/include -I${RADAU5}/src"

SLIBS="${SUNDIALS_BUILD}/src/nvector/serial/libsundials_nvecserial.a \
${SUNDIALS_BUILD}/src/sunmatrix/dense/libsundials_sunmatrixdense.a \
${SUNDIALS_BUILD}/src/sunmatrix/band/libsundials_sunmatrixband.a \
${SUNDIALS_BUILD}/src/sunmatrix/sparse/libsundials_sunmatrixsparse.a \
${SUNDIALS_BUILD}/src/sunlinsol/dense/libsundials_sunlinsoldense.a \
${SUNDIALS_BUILD}/src/sunlinsol/band/libsundials_sunlinsolband.a \
${SUNDIALS_BUILD}/src/sunlinsol/klu/libsundials_sunlinsolklu.a \
${SUNDIALS_BUILD}/src/sundials/libsundials_core.a"

MPI_LIB="-L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi"
KLU_LIB="-lklu -lamd -lcolamd -lbtf -lsuitesparseconfig"

RADAU5_SRCS="${RADAU5}/src/radau5.c ${RADAU5}/src/radau5_linsys.c ${RADAU5}/src/radau5_newt.c ${RADAU5}/src/radau5_estrad.c ${RADAU5}/src/radau5_step.c ${RADAU5}/src/radau5_contr.c ${RADAU5}/src/radau5_ic.c ${RADAU5}/src/radau5_colgroup.c ${RADAU5}/src/radau5_root.c"

mkdir -p bin

build_example() {
    local name=$1
    echo "Building ${name}..."
    ${CC} -O2 -o bin/${name} ${INC} examples/${name}.c ${RADAU5_SRCS} ${SLIBS} ${MPI_LIB} ${KLU_LIB} -lm 2>&1
}

for ex in "$@"; do
    build_example "$ex"
done
