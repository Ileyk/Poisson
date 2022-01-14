!> Input/Output
module mod_io

implicit none

public :: save_vec

contains

! -----------------------------------------------------------------------------
!> Serial writing of outputs
! -----------------------------------------------------------------------------
subroutine save_vec(N1,it,base_filename,vec)
integer, intent(in) :: N1, it
double precision, intent(in) :: vec(N1)
character(len=*), intent(in) :: base_filename
integer :: i
character(len=256) :: filename
character(len=*), parameter :: data_fmt='.dat', output_folder='output/'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

write(filename,"(A100,I0.6,A4)") adjustl(trim(output_folder))//adjustl(trim(base_filename)), it, adjustl(data_fmt)
open(1,file=adjustl(filename))
do i=1,N1
  write(1,*) vec(i)
enddo
close(1)

end subroutine save_vec
! -----------------------------------------------------------------------------

end module mod_io
