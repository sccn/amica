#!/usr/bin/env bash
# Build amica15 with gfortran + Open MPI (no Intel/MKL), on macOS or Linux.
#
# amica15.f90 is an MPI + OpenMP + LAPACK program, so this needs an MPI Fortran
# wrapper (mpif90) plus LAPACK/BLAS. The resulting binary runs as a single MPI
# rank (Open MPI singleton -- no mpirun needed); set threads with OMP_NUM_THREADS.
#
# Dependencies:
#   macOS (Homebrew):   brew install gcc open-mpi lapack
#   Debian/Ubuntu:      sudo apt-get install -y gfortran libopenmpi-dev openmpi-bin liblapack-dev libblas-dev
#
# Usage:
#   bash build_gfortran.sh
#   FC=mpiifort bash build_gfortran.sh                     # override the compiler
#   LAPACK_LIBS="-framework Accelerate" bash build_gfortran.sh   # override linkage
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fc="${FC:-mpif90}"
cc="${CC:-cc}"
out="${OUT:-$here/amica15}"
os="$(uname -s)"

# LAPACK/BLAS linkage: apt libs on Linux; Homebrew lapack (falling back to the
# Accelerate framework) on macOS. Override with LAPACK_LIBS=... .
if [[ -n "${LAPACK_LIBS:-}" ]]; then
  lapack_libs="$LAPACK_LIBS"
elif [[ "$os" == "Darwin" ]]; then
  if brew_lapack="$(brew --prefix lapack 2>/dev/null)" && [[ -d "$brew_lapack/lib" ]]; then
    lapack_libs="-L$brew_lapack/lib -llapack -lblas"
  else
    lapack_libs="-framework Accelerate"
  fi
else
  lapack_libs="-llapack -lblas"
fi

echo "amica gfortran build: $os $(uname -m), FC=$fc, LAPACK=$lapack_libs -> $out"

if ! command -v "$fc" >/dev/null 2>&1; then
  echo "ERROR: '$fc' not found -- install an MPI Fortran compiler (see the header)." >&2
  exit 1
fi

# funmod2 first (amica15 does `use funmod2`). -cpp resolves the #ifdef MKL guards
# (MKL undefined -> the plain-Fortran / vrda branch is used, resolved by
# vmath_shim.c). -std=legacy + -fallow-argument-mismatch relax obsolescent
# constructs; -ffree-line-length-none lifts gfortran's 132-column cap.
"$cc" -O3 -c "$here/vmath_shim.c" -o "$here/vmath_shim.o"
# shellcheck disable=SC2086  # $lapack_libs must word-split into separate flags
"$fc" -O3 -fopenmp -cpp -ffree-line-length-none -std=legacy -fallow-argument-mismatch \
  "$here/funmod2.f90" "$here/amica15.f90" "$here/vmath_shim.o" \
  -o "$out" $lapack_libs
rm -f "$here/vmath_shim.o" "$here"/*.mod

echo "built: $out"
echo "run:   OMP_NUM_THREADS=4 $out input.param"
