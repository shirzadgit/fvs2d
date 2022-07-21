ifort ios_unstrc.f90 tecplot_unstrc.f90  grid_procs.f90 ios2tecplot.f90 -r8 -O3 -convert big_endian -traceback -o ios2tecplot.x -L ./ -ltecio -lstdc++
rm -f *.o *.mod
