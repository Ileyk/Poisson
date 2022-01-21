module mod_bc

implicit none

public :: get_bc
public :: get_bc_BiCGSTAB_p
public :: get_bc_BiCGSTAB_s

contains

! -----------------------------------------------------------------------------
!> Set physical boundary conditions for scalar field y
!>
!> @todo add new types of boundary conditions
! -----------------------------------------------------------------------------
subroutine get_bc(N1,NGC,iOmin1,iOmax1,bc_type,x,w,y)
integer, intent(in) :: N1, NGC
integer, intent(in) :: iOmin1, iOmax1
character(len=*), intent(in) :: bc_type
double precision, intent(in) :: x(N1)
double precision, intent(in) :: w(N1)
double precision, intent(inout) :: y(N1)
integer :: i
double precision :: dpi=dacos(-1.d0)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if     (bc_type=='periodic') then
  do i=1,NGC
    y(i)=y(iOmax1-(NGC-i))
  enddo
  do i=iOmax1+1,N1
    y(i)=y(iOmin1+(i-(iOmax1+1)))
  enddo

else if (bc_type=='usr') then
  ! no rm_avg
  y(1 )= -1.2033396041666666d0
  y(N1)= -2.2966603958333334d0
  ! rm_avg
  ! y(1 )=1.5032390970833327d0
  ! y(N1)=-3.3069490279166658d0
  ! y(1)=w(1)*1.d-4/(y(3)-2.d0*y(2))
  ! y(N1)=w(N1)*1.d-4/(y(N1-2)-2.d0*y(N1-1))
endif

end subroutine get_bc
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Set physical boundary conditions for scalar fields p and s
!> in the BiCGSTAB method
! -----------------------------------------------------------------------------
subroutine get_bc_BiCGSTAB_p(N1,NGC,iOmin1,iOmax1,bc_type,y)
integer, intent(in) :: N1, NGC
integer, intent(in) :: iOmin1, iOmax1
character(len=*), intent(in) :: bc_type
double precision, intent(inout) :: y(N1)
integer :: i
double precision :: dpi=dacos(-1.d0)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
y(1 )=y(2)
y(N1)=y(N1-1)
end subroutine get_bc_BiCGSTAB_p
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Set physical boundary conditions for scalar fields p and s
!> in the BiCGSTAB method
! -----------------------------------------------------------------------------
subroutine get_bc_BiCGSTAB_s(N1,NGC,iOmin1,iOmax1,bc_type,y)
integer, intent(in) :: N1, NGC
integer, intent(in) :: iOmin1, iOmax1
character(len=*), intent(in) :: bc_type
double precision, intent(inout) :: y(N1)
integer :: i
double precision :: dpi=dacos(-1.d0)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
y(1 )=y(2)
y(N1)=y(N1-1)
end subroutine get_bc_BiCGSTAB_s
! -----------------------------------------------------------------------------

end module mod_bc
