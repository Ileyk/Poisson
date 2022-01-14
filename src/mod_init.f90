!> Initialization of quantities on the grid
module mod_init
use mod_params

implicit none

public :: init

contains

! -----------------------------------------------------------------------------
!> Initialize w scalar field on grid
!> The spatial averaged (i) does not play any role in the Poisson
!> equation (since its Laplacian is zero) and (ii) needs to be subtracted
!> for the Jacobi method to work
! -----------------------------------------------------------------------------
subroutine init(N1,iOmin1,iOmax1,min1,max1,x,w)
integer, intent(in) :: N1
integer, intent(in) :: iOmin1, iOmax1
double precision, intent(in) :: min1, max1
double precision, intent(in) :: x(N1)
double precision, allocatable, intent(out) :: w(:)
integer :: i
double precision :: w_avg
double precision :: dpi=dacos(-1.d0)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

allocate(w(N1))
do i=iOmin1,iOmax1
  w(i)=dsin(2.d0*dpi*kwv*((x(i)-min1)/(max1-min1)))
enddo

! w_avg=sum(w)/N1
! w=w-w_avg

! w(N1/2+1)=1.d0

end subroutine init
! -----------------------------------------------------------------------------

end module mod_init
