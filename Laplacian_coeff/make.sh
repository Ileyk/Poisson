mpif90 -c ../src/mod_linal.f90 laplacian_coeff.f90
mpif90 mod_linal.o laplacian_coeff.o -o laplacian_coeff
rm *.o *.mod
