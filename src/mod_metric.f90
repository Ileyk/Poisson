! -----------------------------------------------------------------------------
!> Module where you specify the functions corresponding (i) to each element
!> of the invert of the space metric gm (i.e. gm^{rr}, gm^{rt}, etc)
!> and (ii) to the determinant of the space metric.
!> We also specify the required derivatives.
! -----------------------------------------------------------------------------
module mod_metric

implicit none

contains

! -----------------------------------------------------------------------------
! det(gm)
! -----------------------------------------------------------------------------
function det_root_ (r,t)
double precision :: det_root_
double precision, intent(in) :: r, t
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
det_root_=r**2.d0*(dsin(t))**2.d0
return
end function det_root_
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! d(det(gm))/dr
! -----------------------------------------------------------------------------
function det_root_dr_ (r,t)
double precision :: det_root_dr_
double precision, intent(in) :: r, t
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
det_root_dr_=2.d0*r*(dsin(t))**2.d0
return
end function det_root_dr_
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! gm^{rr}
! -----------------------------------------------------------------------------
function rr_ (r)
double precision :: rr_
double precision, intent(in) :: r
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rr_=1.d0
return
end function rr_
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! d(gm^{rr})/dr
! -----------------------------------------------------------------------------
function rr_dr_ (r)
double precision :: rr_dr_
double precision, intent(in) :: r
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rr_dr_=0.d0
return
end function rr_dr_
! -----------------------------------------------------------------------------

end module mod_metric
