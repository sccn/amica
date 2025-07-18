# BUILDING FOR UNIVERSITY OF FLORIDA'S HIPERGATOR (REHL9)

module purge
module load gcc/12.2.0 openmpi/5.0.7 openblas lapack mkl

cd /location/of/amica/repo/
mpifort -cpp -fopenmp -O3 funmod2.f90 amica17.f90 \
  -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 \
  -Wl,--start-group -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -Wl,--end-group \
  -fopenmp -lpthread -lm -ldl -o amica17ub

# VARIOUS COMPILATION COMMANDS
## For GCC + MKL (sequential, no threading):
mpifort -cpp -fopenmp -O3 funmod2.f90 amica17.f90 \
  -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 \
  -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group \
  -lpthread -lm -ldl -o amica17ub

## For GCC + MKL (OpenMP threading):
mpifort -cpp -fopenmp -O3 funmod2.f90 amica17.f90 \
  -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 \
  -Wl,--start-group -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -Wl,--end-group \
  -fopenmp -lpthread -lm -ldl -o amica17ub

## For GCC + MKL (runtime, recommended for simplicity):
mpifort -cpp -fopenmp -O3 funmod2.f90 amica17.f90 \
  -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 \
  -lmkl_rt -lpthread -lm -ldl -o amica17ub

# NOTES & EDITS

## (07/17/2025) JS EDITS TO AMICA15,
removed dependencies of mkl libraries "vrda" "vml" and used native fortran functions. The code can be compiled without mkl dependencies, hence the addition of "_nothread".

compile command:
mpifort -cpp -fopenmp -O3 -DMKL funmod2.f90 amica15.f90 -L/apps/gcc/12.2.0/openblas/0.3.28/lib64 \
 -L/apps/gcc/12.2.0/lapack/3.11.0/lib -lopenblas -llapack -o amica15ub



## (07/17/2025) JS EDITS TO AMICA17, 
updating to use newer compiler. Needed to directly link mkl libraries inside intel directroy. HPG support team
highly discouraged just using the mkl module directly. some optinos for parallel computing compiler mpifort:

-liomp5: Intel OpenMP library for parallel processing
    -lmkl_sequential: Sequential Intel MKL library for basic linear algebra operations (option)
    -lmkl_intel_thread: Intel MKL library for multi-threaded operations (option)
-lmkl_intel_lp64: Intel MKL library for 64-bit linear algebra
-lmkl_core: Core Intel MKL library for basic linear algebra operations
-Wl,--start-group and -Wl,--end-group: These options are used to group the libraries together for linking, which is necessary when there are circular dependencies between libraries.

EDITS:
Needed to edit amica17.f90 for line character limits (132 ch per line) added '&' signs for line coninuation
Needed to edit amica17_header.f90 to include seed_size, seed_array(:) definitions. In coordination, had to
    add lines to amica17.f90 to update call to random_seed()

