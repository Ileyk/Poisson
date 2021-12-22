module mod_grid

implicit none

public :: make_grid

contains

! -----------------------------------------------------------------------------
!> Make 1D uniform Cartesian mesh
!> Inputs
!>  - N1: number of cells in direction # 1
!>  - min1: minimum edge of direction # 1
!>  - max1: maximum edge of direction # 1
! -----------------------------------------------------------------------------
subroutine make_grid(N1,min1,max1,x)
integer, intent(in) :: N1
double precision, intent(in) :: min1, max1
double precision, allocatable, intent(out) :: x(:)
integer :: i
double precision :: dx
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

allocate(x(N1))
dx=(max1-min1)/dble(N1)
do i=1,N1
  x(i)=min1+dx*(i-0.5)
enddo

end subroutine make_grid
! -----------------------------------------------------------------------------

end module mod_grid
