cd src
rm -f *.o *.mod aggregation
gfortran -O3 -c potential.f90 integrators.f90
gfortran -O3 -o aggregation aggregation.f90 *.o
