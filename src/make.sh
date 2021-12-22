mpif90 -c mod_basic_types.f90 mod_csts.f90 mod_params.f90 mod_mpi.f90 mod_io.f90 mod_init.f90 mod_grid.f90 mod_poisson.f90 main.f90
mpif90 mod_basic_types.o mod_csts.f90 mod_params.o mod_mpi.o mod_io.o mod_init.o mod_grid.o mod_poisson.o main.o -o exe
mv exe ../scratch/.
rm *.o *.mod
