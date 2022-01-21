program poisson

use mod_grid
use mod_init
use mod_mpi
use mod_poisson
use mod_io
use mod_bc

use mpi_f08

implicit none

double precision, allocatable :: x(:), w(:)

call init_MPI
if (id==0) print*, 'Let it snow on', NPROC, 'processors'

if (NPROC/=1) then
  call mpistop("Not MPI parallelized yet")
endif

call check_par

call init_par

call make_grid(N1_I,NGC,iImin1,iImax1,min1,max1,x)

call init(N1_I,min1,max1,x,w)

! call get_bc(N1_I,NGC,iOmin1,iOmax1,bc_type,w)

call save_vec(N1_I,0,'w_',w)

call solver(N1_I,NGC,iOmin1,iOmax1,min1,max1,pencil,solver_type,bc_type,x,w)

call finalize_MPI

if (id==0) print*, 'Done'

end program poisson
