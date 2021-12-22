module mod_push
use class_pcl
implicit none

public :: push

contains

! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
subroutine push(dt,g,pcls)
double precision, intent(in) :: dt, g
! "class" instead of "type" to say that it could be any
! extended type from the type_pcl type
class(type_pcl), dimension(:), intent(inout) :: pcls
! double precision, dimension(3) :: x, v
integer :: N_pcl, i_pcl
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
N_pcl=size(pcls)
! I need to be able to do this
! pcls%pos=pcls%pos+dt*pcls%vel
! pcls%vel=pcls%vel+dt*g
do i_pcl=1,N_pcl
  ! x = get_pos(pcls(i_pcl))
  ! v = get_vel(pcls(i_pcl))
  ! x   =max(0.d0,x+dt*v)
  ! v(3)=v(3)+dt*g
  ! call set_pos(x,pcls(i_pcl))
  ! call set_vel(v,pcls(i_pcl))
  ! Which dynamic equation for particles? Depends on their specific type,
  ! extended from type_pcl.
  select type (pcls)
  class is (type_electrons)
    pcls(i_pcl)%pos   =max(0.d0,pcls(i_pcl)%pos +dt*pcls(i_pcl)%vel)
    pcls(i_pcl)%vel(3)=       pcls(i_pcl)%vel(3)+dt*g
  class is (type_positrons)
    pcls(i_pcl)%pos   =min(1.d0,pcls(i_pcl)%pos +dt*pcls(i_pcl)%vel)
    pcls(i_pcl)%vel(3)=       pcls(i_pcl)%vel(3)-dt*g
  class default
    stop 'Unknown type'
  end select
enddo
end subroutine push
! -----------------------------------------------------------------------------

! ! -----------------------------------------------------------------------------
! ! -----------------------------------------------------------------------------
! subroutine remove(pcls)
! ! type(type_pcl), dimension(:), intent(inout) :: pcls
! class(type_pcl), dimension(:), allocatable, intent(inout) :: pcls
! class(type_pcl), dimension(:), allocatable :: pcls_tmp
! integer :: N_pcl, i_pcl, N_rm
! logical, dimension(:), allocatable :: mask
! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! N_pcl=size(pcls)
! allocate(mask(N_pcl))
!
! ! mask true if pcls to be kept
! select type (pcls)
! class is (type_electrons)
!   mask=pcls(:)%pos(3)>0.d0
! class is (type_positrons)
!   mask=pcls(:)%pos(3)<1.d0
! class default
!   stop 'Unknown type'
! end select
!
! N_rm=N_pcl-count(mask)
! print*, N_rm
! ! where(mask)
!
! pcls_tmp = pack(pcls,mask)
! deallocate(pcls)
! call move_alloc(pcls_tmp,pcls)
! deallocate(pcls_tmp)
!
! ! print*, N_rm
! ! N_pcl=size(pcls)
!
! ! print*, pack(pcls(:)%pos(3), pcls(:)%pos(3)>1.d0)
! ! print*, count(pcls(:)%pos(3)>1.d0-1.d-12)
!
! ! do i_pcl=1,N_pcl
! !   if (pcls(i_pcl)%x(3)>1.d0 .or. pcls(i_pcl)%x(3))
! ! enddo
!
! deallocate(mask)
!
! end subroutine remove
! ! -----------------------------------------------------------------------------

end module mod_push
