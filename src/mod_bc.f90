module mod_bc

implicit none

public :: get_bc

contains

! -----------------------------------------------------------------------------
!> Set physical boundary conditions for scalar field y
!>
!> Fix the value of y in ghost cells
!> @todo add new types of boundary conditions
! -----------------------------------------------------------------------------
subroutine get_bc(N1,NGC,iOmin1,iOmax1,bc_type,y)
integer, intent(in) :: N1, NGC
integer, intent(in) :: iOmin1, iOmax1
character(len=*), intent(in) :: bc_type
double precision, intent(inout) :: y(N1)
integer :: i
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if (bc_type=='periodic') then
  do i=1,NGC
    y(i)=y(iOmax1-(NGC-i))
  enddo
  do i=iOmax1+1,N1
    y(i)=y(iOmin1+(i-(iOmax1+1)))
  enddo
endif

end subroutine get_bc
! -----------------------------------------------------------------------------

end module mod_bc
