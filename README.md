# amica
Code for AMICA: Adaptive Mixture ICA with shared component

TO COMPILE WITH INTEL FORTRAN ON WINDOWS

1. Install Intel Fortran compiler for Windows.
2. Install MPICH2 library (fmpich2.lib) for Windows.
3. In a cmd.exe windows: Run compilervars.bat with argument (e.g. intel64): 

   \> "c:\Program Files (x86)\Intel\Composer XE 2011 SP1\bin\compilervars.bat" intel64

4. Compile Amica with command like (/F sets the stack size):

   \> ifort   /Qopenmp /Qmkl  /F2147483648 /DMKL /fpp  /O3  /exe:amica15mkl.exe  funmod2.f90 amica15.f90 fmpich2.lib

5. Test:
   \> .\amica15mkl.exe .\amicadefs.param



TO COMPILE WITH INTEL FORTRAN ON UBUNTU

1. Install Intel Fortran compiler for Linux.
2. Compile MPICH2 setting environmental vars CC, CXX, FC, and F77 to icc and ifort.
3. Compile Amica with command like:

   $ /home/jason/mpich2-3.2-install/bin/mpif90 -I/opt/intel/mkl/include/ -fpp -qopenmp -O3 -mkl -static -static-intel -DMKL funmod2.f90 amica15.f90 -o amica15ub

4. Test:

   $ ./amica15ub ./amicadefs.param


TO COMPILE WITH INTEL FORTRAN ON MAC

1. Install Intel Fortran compiler for Mac/Linux.
2. Compile MPICH2 setting environmental vars CC, CXX, FC, and F77 to icc and ifort.
3. Compile Amica with command like:

   $ ~/mpich-3.2-install/bin/mpif90 -L/Users/jason/mpich-3.2-install/lib/ -I/Users/jason/mpich-3.2-install/include/ -qopenmp -mkl -static-intel -O3 -fpp -DMKL amica15.f90 funmod2.f90 -o amica15mac
   

4. Test:

   $ ./amica15mac ./amicadefs.param
