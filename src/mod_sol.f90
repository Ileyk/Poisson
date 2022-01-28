module mod_sol

implicit none

contains

! -----------------------------------------------------------------------------
! det(gm)
! -----------------------------------------------------------------------------
function analytic_ (x,y)
use mod_params
use mod_csts, only : dpi
double precision :: analytic_
double precision, intent(in) :: x, y
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! analytic_=-(1.d0/(2.d0*(2.d0*dpi*kwv/(max1-min1))**2.d0))*&
!          dsin(2.d0*dpi*kwv*((x-min1)/(max1-min1)))*&
!          dsin(2.d0*dpi*kwv*((y-min2)/(max2-min2)))
! analytic_=-(1.d0/(2.d0*dpi*kwv/(max1-min1))**2.d0)*&
!          dexp(dsin(2.d0*dpi*kwv*((x-min1)/(max1-min1))))*&
!          dexp(dcos(2.d0*dpi*kwv*((y-min2)/(max2-min2))))
! analytic_=x**4.d0/12.d0+y**4.d0/12.d0
analytic_=1.d0/x+1.d0

return
end function analytic_
! -----------------------------------------------------------------------------

end module mod_sol
