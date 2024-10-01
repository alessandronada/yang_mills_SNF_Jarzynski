#!/bin/bash

NCOLORS=3
NTHREADS=1
compile_targets='yang_mills_local yang_mills_local_pt yang_mills_local_pt_multicanonic conf_check yang_mills_t0 yang_mills_gradient_flow yang_mills_local_jarzynski_bc yang_mills_local_jarzynski_beta'

# 1=yes, 0=no
ENABLE_THETA=0
ENABLE_OPENMP=0

flag_theta=''
flag_openmp=''
if [ "${ENABLE_THETA}" -eq '1' ]; then flag_theta='--enable-use-theta'; fi
if [ "${ENABLE_OPENMP}" -eq '1' ]; then flag_openmp='--enable-use-openmp'; fi

chmod +x configure
if [ -d lib/.deps ]; then make clean; fi

# standard configuration with gcc + standard -O3 optimizations
#./configure N_c=${NCOLORS} Num_threads=${NTHREADS} CC=gcc CFLAGS='-O3' ${flag_openmp} ${flag_theta}

# optimized configuration for Marconi (optimized compilation with Intel compiler on Intel Skylake processors)
./configure N_c=${NCOLORS} Num_threads=${NTHREADS} CC=icc CFLAGS='-O3 -axCORE-AVX512 -mtune=skylake -ip -ipo' LIBS="-ldl -lz -lc" ${flag_openmp} ${flag_theta}

# optimized configuration for Galileo100 (optimized compilation with Intel compiler on Intel Cascadelake processors)
#./configure N_c=${NCOLORS} Num_threads=${NTHREADS} CC=icc CFLAGS='-O3 -axCORE-AVX512 -mtune=cascadelake -ip -ipo' LIBS="-ldl -lz -lc" ${flag_openmp} ${flag_theta}

# compile
make ${compile_targets} # -j 18
