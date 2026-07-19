![Screenshot 2024-07-25 at 14 12 51](https://github.com/user-attachments/assets/e5880f21-51ae-4ebf-a3de-bc040d6f4aab)

# The AMICA algorithm and EEGLAB plugin

Adaptive Mixture Independent Component Analysis (AMICA) is a program (for Linux, Mac, and Windows) that performs an independent component analysis (ICA) decomposition on input data, potentially with multiple ICA models. It can be run standalone, or from MATLAB.

Code for AMICA: Adaptive Mixture ICA with shared component

Refer to the [Amica wiki](https://github.com/japalmer29/amica/wiki) for documentation or the other menus if you are on the EEGLAB website.

Refer also to Jason Palmer's [AMICA page](https://sccn.ucsd.edu/~jason/amica_web.html).

# Compilation

Prebuilt binaries are included in this repository: `amica15mac` (Intel macOS), `amica15ub`
(Ubuntu/Linux), `amica15mkl.exe` (Windows), and `amica15_macos_arm64` (Apple Silicon, native,
built with the gfortran path below; it links Homebrew's Open MPI and LAPACK, so running it needs
`brew install gcc open-mpi lapack`). If you already have that toolchain, building from source is
just as easy:

## Build from source with gfortran + Open MPI (macOS and Linux, no Intel/MKL)

`amica15.f90` is an MPI + OpenMP + LAPACK program, so it needs an MPI Fortran wrapper
(`mpif90`) plus LAPACK/BLAS, not plain `gfortran`.

Install the toolchain:

- **macOS / Apple Silicon (Homebrew, no sudo):** `brew install gcc open-mpi lapack`
- **Debian/Ubuntu:** `sudo apt-get install -y gfortran libopenmpi-dev openmpi-bin liblapack-dev libblas-dev`

Then build and run:

```bash
bash build_gfortran.sh                 # -> ./amica15
OMP_NUM_THREADS=4 ./amica15 input.param
```

The binary runs as a single MPI rank (Open MPI singleton, no `mpirun` needed); OpenMP thread
count is set with `OMP_NUM_THREADS`. `build_gfortran.sh` selects LAPACK per platform (Homebrew
`lapack`, falling back to the Accelerate framework, on macOS; `-llapack -lblas` on Linux) and
passes the flags gfortran needs for this Intel-oriented source: `-cpp` (resolve the `#ifdef MKL`
guards so the MKL include is skipped), `-std=legacy -fallow-argument-mismatch` (relax obsolescent
constructs), and `-ffree-line-length-none`. `vmath_shim.c` supplies `vrda_exp`/`vrda_log` (the
non-MKL branch's vectorized exp/log, otherwise from AMD LibM) as libm loops, so no vendor math
library is required.

For Intel OneAPI / Windows builds, see the [Amica wiki](https://github.com/japalmer29/amica/wiki).

# Version history

1.7 - Update documentation for pop_runamica and add test file

1.6.1 - Modify Windows compilation instructions. Intel OneAPI should be tested for Mac and Ubuntu.

1.6 - Deprecate Comet and replace with Expanse supercomputer executable

1.5.2 - Comet supercomputer executable

