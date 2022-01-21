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
subroutine init(N1,min1,max1,x,w)
integer, intent(in) :: N1
double precision, intent(in) :: min1, max1
double precision, intent(in) :: x(N1)
double precision, allocatable, intent(out) :: w(:)
integer :: i
double precision :: w_avg
double precision :: dpi=dacos(-1.d0)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

allocate(w(N1))
w=0.d0
do i=1,N1
  ! w(i)=dsin(2.d0*dpi*kwv*((x(i)-min1)/(max1-min1)))
  ! w(i)= dsin(2.d0*dpi*kwv*((x(i)-min1)/(max1-min1)))       *dexp(dsin(2.d0*dpi*kwv*((x(i)-min1)/(max1-min1))))-&
  !      (dcos(2.d0*dpi*kwv*((x(i)-min1)/(max1-min1))))**2.d0*dexp(dsin(2.d0*dpi*kwv*((x(i)-min1)/(max1-min1))))
  ! w(i)= dsin(2.d0*dpi*kwv*((x(i)-min1)/(max1-min1)))/(2.d0+dsin(2.d0*dpi*kwv*((x(i)-min1)/(max1-min1))))**2.d0+&
  !       2.d0*(dcos(2.d0*dpi*kwv*((x(i)-min1)/(max1-min1))))**2.d0/(2.d0+dsin(2.d0*dpi*kwv*((x(i)-min1)/(max1-min1))))**3.d0
  w(i)=x(i)**2.d0 ! -1.5d0
enddo

! w_avg=sum(w)/N1
! print*, w_avg
! w=w-w_avg

! w(N1/2+1)=1.d0

end subroutine init
! -----------------------------------------------------------------------------

end module mod_init
