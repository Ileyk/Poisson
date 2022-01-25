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
!> @warning for non-periodic functions, VITAL to let BCs loose
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

  ! y(1 )= 0.19636568138139013d0 ! 0.058800136d0 ! -1.2033396041666666d0
  ! y(2 )= 0.20518173008201698d0
  ! y(N1-1)= 0.19177611863281882d0
  ! y(N1)= 0.17499802874027051d0 ! -2.2966603958333334d0

  ! y(1)=0.20511991814902028d0
  ! y(N1)=0.20063367205550098d0

  ! y(1 )= 0.08218803305d0
  ! y(N1)= 1.3519784507d0

  y(1 )=2.0034657220268395d0
  y(N1)=1.4982671389865767d0

  ! rm_avg
  ! y(1 )=1.5032390970833327d0
  ! y(N1)=-3.3069490279166658d0
  ! y(1)=w(1)*1.d-4/(y(3)-2.d0*y(2))
  ! y(N1)=w(N1)*1.d-4/(y(N1-2)-2.d0*y(N1-1))
  ! y(1 )=y(2)
  ! y(N1)=y(N1-1)
endif

end subroutine get_bc
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Set physical boundary conditions for scalar fields p
!> in the BiCGSTAB method
!>
!> @warning for non-periodic functions, VITAL to set p BCs to 0
! -----------------------------------------------------------------------------
subroutine get_bc_BiCGSTAB_p(N1,NGC,iOmin1,iOmax1,bc_type,y)
integer, intent(in) :: N1, NGC
integer, intent(in) :: iOmin1, iOmax1
character(len=*), intent(in) :: bc_type
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
  do i=1,NGC
    y(i)=0.d0 ! y(2)
  enddo
  do i=iOmax1+1,N1
    y(i)=0.d0 ! y(N1-1)
  enddo
endif
end subroutine get_bc_BiCGSTAB_p
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Set physical boundary conditions for scalar fields s
!> in the BiCGSTAB method
!>
!> @warning for non-periodic functions, VITAL to set s BCs to 0
! -----------------------------------------------------------------------------
subroutine get_bc_BiCGSTAB_s(N1,NGC,iOmin1,iOmax1,bc_type,y)
integer, intent(in) :: N1, NGC
integer, intent(in) :: iOmin1, iOmax1
character(len=*), intent(in) :: bc_type
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
  do i=1,NGC
    y(i)=0.d0 ! y(2)
  enddo
  do i=iOmax1+1,N1
    y(i)=0.d0 ! y(N1-1)
  enddo
endif
end subroutine get_bc_BiCGSTAB_s
! -----------------------------------------------------------------------------

end module mod_bc
