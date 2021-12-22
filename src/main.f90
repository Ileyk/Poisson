program poisson

use mod_grid
use mod_init
use mod_mpi
use mod_poisson
use mod_io

! use mod_share
! use mod_rdm
! use class_pcl
! use mod_push
! use mod_io

use mpi_f08

implicit none

double precision, parameter :: min1=-0.5d0, max1=0.5d0
integer, parameter :: N1=100
double precision, allocatable :: x(:), w(:), f(:)

! Initialize the MPI environment
call init_MPI
if (id==0) print*, 'Let it snow on', NPROC, 'processors'

if (NPROC/=1) then
  call mpistop("Not MPI parallelized yet")
endif

call make_grid(N1,min1,max1,x)
call init(N1,min1,max1,x,w)

call save_vec(N1,0,'w_',w)

call frst_gss(N1,min1,max1,x,w,f)

call save_vec(N1,0,'f_',f)

call solver(N1,min1,max1,x,w,f)

call finalize_MPI

if (id==0) print*, 'Done'

end program poisson
