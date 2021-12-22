module mod_mpi
use mpi_f08
implicit none

private

public :: init_MPI
public :: finalize_MPI
public :: mpistop

!> The number of MPI tasks
integer, public :: NPROC
!> The rank of the current MPI task
integer, public :: id
!> The MPI communicator
! integer :: icomm ! old syntax, for "use mpi" and "include 'mpif.h'"
type(MPI_Comm), public :: icomm ! for "use mpi_f08"
! type(MPI_Info) :: iinfo
!> A global MPI error return code
integer, public :: ierr
! !> coordinates of the process id in the cartesian topology,
! !> i.e. coordinates of the block, from (0,0) to (NPR-1,NTH-1)
! integer, dimension(2), public :: coords
! !> neighbor array (in 2D, 8 neighbors)
! integer, dimension(8), public :: ngh

contains


!=======================================================================
!> Initialize the MPI environment
!=======================================================================
subroutine init_MPI
use mod_basic_types, only : size_real, size_double, size_int, size_logical
! ----------------------------------------------------------------------
integer(kind=MPI_ADDRESS_KIND) :: lb
integer(kind=MPI_ADDRESS_KIND) :: sizes
! ----------------------------------------------------------------------

! Initialize MPI
call MPI_INIT(ierr)

! Test whether MPI is working properly or not
IF (ierr.NE.MPI_SUCCESS) THEN
   print*,'Error starting MPI program. Terminating.'
   call MPI_ABORT(MPI_COMM_WORLD,ierr)
END IF

! Get the ID number of each process.
! Each process stores its rank, which ranges from 0 to NPROC-1,
! where NPROC is the number of processes.
call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierr)

! Get the number of processes NPROC
call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,ierr)

! Use the default communicator, which contains all the processes
icomm = MPI_COMM_WORLD

! Get size of double/integer
call MPI_TYPE_GET_EXTENT(MPI_REAL,lb,sizes,ierr)
if (sizes /= size_real) call mpistop("Incompatible real size")
call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,sizes,ierr)
if (sizes /= size_double) call mpistop("Incompatible double size")
call MPI_TYPE_GET_EXTENT(MPI_INTEGER,lb,sizes,ierr)
if (sizes /= size_int) call mpistop("Incompatible integer size")
call MPI_TYPE_GET_EXTENT(MPI_LOGICAL,lb,sizes,ierr)
if (sizes /= size_logical) call mpistop("Incompatible logical size")

end subroutine init_MPI
!=======================================================================


!=======================================================================
!> Finalize (or shutdown) the MPI environment
!=======================================================================
subroutine finalize_MPI
! ----------------------------------------------------------------------

call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call MPI_FINALIZE(ierr)

end subroutine finalize_MPI
!=======================================================================


!=======================================================================
!> Exit MPI-AMRVAC with an error message
!=======================================================================
subroutine mpistop(message)
! ----------------------------------------------------------------------
character(len=*), intent(in) :: message !< The error message to be printed
integer                      :: ierrcode
! ----------------------------------------------------------------------

write(*, *) "ERROR for processor", id, ":"
write(*, *) trim(message)

call MPI_ABORT(icomm,ierrcode,ierr)

end subroutine mpistop
!=======================================================================

end module mod_mpi
