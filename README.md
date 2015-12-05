# amica
Code for AMICA: Adaptive Mixture ICA with shared component

TO COMPILE WITH INTEL FORTRAN ON WINDOWS

1. In a cmd.exe windows: Run compilervars.bat with argument (e.g. intel64): 

   "c:\Program Files (x86)\Intel\Composer XE 2011 SP1\bin\compilervars.bat" intel64

2. Compile with command like (/F sets the stack size):

   ifort   /Qopenmp /Qmkl  /F2147483648 /DMKL /fpp  /O3  /exe:amica15mkl.exe  funmod2.f90 amica15.f90 fmpich2.lib

3. Test:

   .\amica15mkl.exe .\amicadefs.param
