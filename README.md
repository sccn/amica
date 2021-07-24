# Amica
Code for AMICA: Adaptive Mixture ICA with shared component

Refer to the [Amica wiki](https://github.com/japalmer29/amica/wiki) for documentation.

## TO COMPILE WITH INTEL FORTRAN ON MAC

1. Install Intel Fortran compiler for Mac/Linux (free demo).
   See https://software.intel.com/en-us/intel-parallel-studio-xe

2. Compile MPICH2 setting environmental vars CC, CXX, FC, and F77 to icc and ifort. Set $FBIN to Intel Fortran bin directory.

   i) Download the mpich-3.2 code from: http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz
   
   ii) Compile mpich-3.2:
   
         $ cp /Users/$USER/downloads/mpich-3.2.tar.gz .   
         $ setenv CC $FBIN/icc
         $ setenv CXX $FBIN/icc
         $ setenv F77 $FBIN/ifort
         $ setenv FC $FBIN/ifort   
         $ tar xvf mpich-3.2.tar.gz   
         $ cd mpich-3.2   
         $ ./configure --prefix=/Users/$USER/mpich-3.2-install   
         $ make   
         $ make install

3. Compile Amica with the command:

         $ ~/mpich-3.2-install/bin/mpif90 -L/Users/$USER/mpich-3.2-install/lib/ -I/Users/$USER/mpich-3.2-install/include/ -qopenmp -mkl -static-intel -O3 -fpp -DMKL amica15.f90 funmod2.f90 -o amica15mac

4. Test:

   i) Download Sample EEG Data (Memorize.fdt and amicadefs.param) from: https://sccn.ucsd.edu/~jason/amica_web.html
   
   ii) Test binary:
   
         $ ./amica15mac ./amicadefs.param

## TO COMPILE WITH INTEL FORTRAN ON WINDOWS

1. Install Intel Fortran compiler for Windows.
2. Install MPICH2 library (fmpich2.lib) for Windows.
3. In a cmd.exe windows: Run compilervars.bat with argument (e.g. intel64): 

         > "c:\Program Files (x86)\Intel\Composer XE 2011 SP1\bin\compilervars.bat" intel64

4. Compile Amica with the command (/F sets the stack size):

         > ifort /Qopenmp /Qmkl /F2147483648 /DMKL /fpp  /O3 /exe:amica15mkl.exe funmod2.f90 amica15.f90 fmpich2.lib

5. Test:

         > .\amica15mkl.exe .\amicadefs.param

## TO COMPILE WITH INTEL FORTRAN ON UBUNTU

1. Install Intel Fortran compiler for Linux.
2. Compile MPICH2 setting environmental vars CC, CXX, FC, and F77 to icc and ifort.
3. Compile Amica with the command:

         $ /home/jason/mpich2-3.2-install/bin/mpif90 -I/opt/intel/mkl/include/ -fpp -qopenmp -O3 -mkl -static -static-intel -DMKL funmod2.f90 amica15.f90 -o amica15ub

4. Test:

         $ ./amica15ub ./amicadefs.param
         
## TO COMPILE WITH INTEL FORTRAN ON EXPANSE SUPERCOMPUTER

1. load appropriate modules:

```
   module purge
   module load cpu
   module load intel
   module load intel-mkl
   module load mvapich2

   mpif90 -static-intel -fpp -O3 -march=core-avx2 -heap-arrays \
       -qopenmp -mkl -DMKL -o amica15ex funmod2.f90 amica15.f90
```

2. Compile Amica with the command:

```
   mpif90 -static-intel -fpp -O3 -march=core-avx2 -heap-arrays \
       -qopenmp -mkl -DMKL -o amica15ex funmod2.f90 amica15.f90
```

3. Test:
```
   $ ./amica15ex ./amicadefs.param
```

## VERSION HISTORY

1.6 - Deprecate Comet and replace with Expanse supercomputer executable

1.5.2 - Comet supercomputer executable

