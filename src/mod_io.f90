!> Input/Output
module mod_io

implicit none

public :: save_vec
public :: save_grid

contains

! -----------------------------------------------------------------------------
!> Serial writing of outputs
! -----------------------------------------------------------------------------
subroutine save_grid(N1,N2,it,base_filename,x)
integer, intent(in) :: N1, N2, it
double precision, intent(in) :: x(N1,N2,2)
character(len=*), intent(in) :: base_filename
integer :: i, j
character(len=256) :: filename
character(len=*), parameter :: data_fmt='.dat', output_folder='output/'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

write(filename,"(A100,I0.6,A4)") adjustl(trim(output_folder))//adjustl(trim(base_filename)), it, adjustl(data_fmt)
open(1,file=adjustl(filename))
do i=1,N1
  write(1,*) x(i,1,1)
enddo
do j=1,N2
  write(1,*) x(1,j,2)
enddo
close(1)

end subroutine save_grid
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Serial writing of outputs
! -----------------------------------------------------------------------------
subroutine save_vec(N1,N2,it,base_filename,vec)
integer, intent(in) :: N1, N2, it
double precision, intent(in) :: vec(N1,N2)
character(len=*), intent(in) :: base_filename
integer :: i, j
character(len=256) :: filename
character(len=*), parameter :: data_fmt='.dat', output_folder='output/'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

write(filename,"(A100,I0.6,A4)") adjustl(trim(output_folder))//adjustl(trim(base_filename)), it, adjustl(data_fmt)
open(1,file=adjustl(filename))
do i=1,N1
do j=1,N2
  write(1,*) vec(i,j)
enddo
enddo
close(1)

end subroutine save_vec
! -----------------------------------------------------------------------------

end module mod_io
