heat-equation
=============

Using parallel algorithms (MPI) and cellular automatas to the problem of heat equation

Compilation
===========

Start "make" in the root directory creates several executables, which can be used with

mpirun -nb nb_procs ./executable file

Where :
- "nb_procs" is the number of processors used
- "executable" is the name of the executable (ex : average)
- "file" is the file where the program has to read the data
