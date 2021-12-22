module class_pcl
implicit none

! -----------------------------------------------------------------------------
type, public :: type_pcl
  ! private
  double precision, dimension(3) :: pos
  double precision, dimension(3) :: vel
end type type_pcl
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
type, extends(type_pcl) :: type_electrons
end type type_electrons
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
type, extends(type_pcl) :: type_positrons
end type type_positrons
! -----------------------------------------------------------------------------

! type, public :: type_pcls
!   type(type_pcl), allocatable :: pcl
! end type type_pcls

contains

! Necessary subroutines if attributes of the class are private.
! But why would we do such a thing?...

! subroutine set_pos(x,pcl)
!   double precision, dimension(3), intent(in) :: x
!   class(type_pcl) :: pcl
!   pcl%pos=x
! end subroutine set_pos
!
! subroutine set_vel(v,pcl)
!   double precision, dimension(3), intent(in) :: v
!   class(type_pcl) :: pcl
!   pcl%vel=v
! end subroutine set_vel
!
! function get_pos(pcl) result(x)
!   class(type_pcl), intent(in) :: pcl
!   double precision, dimension(3) :: x
!   x=pcl%pos
! end function get_pos
!
! function get_vel(pcl) result(v)
!   class(type_pcl), intent(in) :: pcl
!   double precision, dimension(3) :: v
!   v=pcl%vel
! end function get_vel

end module class_pcl
