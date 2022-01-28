module mod_bc

implicit none

public :: get_bc
public :: get_bc_BiCGSTAB_p
public :: get_bc_BiCGSTAB_s
public :: set_GC_to_x0

contains

! -----------------------------------------------------------------------------
!> Set physical boundary conditions for scalar field y (not in corners)
!>
!> @todo add new types of boundary conditions
!> @warning for non-periodic functions, VITAL to let BCs loose
! -----------------------------------------------------------------------------
subroutine get_bc(N,N1,N2,NGC,bc_type,x,y)
use mod_sol
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
abstract interface
  function f_xy (x,y)
     double precision :: f_xy
     double precision, intent (in) :: x, y
  end function f_xy
end interface
procedure (f_xy), pointer :: analytic => null ()
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
integer, intent(in) :: N, N1, N2, NGC
character(len=*), intent(in) :: bc_type
double precision, intent(in) :: x(N1,N2,2)
! double precision, intent(in) :: w(N1,N2)
double precision, intent(inout) :: y(N)
integer :: p, k, i, j
double precision :: dpi=dacos(-1.d0)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if     (bc_type=='periodic') then

  call set_GC_to_periodic(N,N1,N2,NGC,y)

else if (bc_type=='zero') then

  call set_GC_to_x0(N,N1,N2,NGC,0.d0,y)

else if (bc_type=='usr') then

  analytic => analytic_

  ! Bottom
  do p=1,NGC
  do k=(p-1)*N1+NGC+1,p*N1-NGC
    j=int((k-1)/N1)+1
    i=k-(j-1)*N1
    y(k)=analytic(x(i,j,1),x(i,j,2))
  enddo
  enddo
  ! Top
  do p=N2-NGC+1,N2
  do k=(p-1)*N1+NGC+1,p*N1-NGC
    j=int((k-1)/N1)+1
    i=k-(j-1)*N1
    y(k)=analytic(x(i,j,1),x(i,j,2))
  enddo
  enddo
  ! Left
  do p=NGC,N2-NGC-1
  do k=p*N1+1,p*N1+NGC
    j=int((k-1)/N1)+1
    i=k-(j-1)*N1
    y(k)=analytic(x(i,j,1),x(i,j,2))
  enddo
  enddo
  ! Right
  do p=NGC,N2-NGC-1
  do k=(p+1)*N1-NGC+1,(p+1)*N1
    j=int((k-1)/N1)+1
    i=k-(j-1)*N1
    y(k)=analytic(x(i,j,1),x(i,j,2))
  enddo
  enddo

  ! no rm_avg

  ! y(1 )= 0.19636568138139013d0 ! 0.058800136d0 ! -1.2033396041666666d0
  ! y(2 )= 0.20518173008201698d0
  ! y(N1-1)= 0.19177611863281882d0
  ! y(N1)= 0.17499802874027051d0 ! -2.2966603958333334d0

  ! y(1)=0.20511991814902028d0
  ! y(N1)=0.20063367205550098d0

  ! y(1 )= 0.08218803305d0
  ! y(N1)= 1.3519784507d0

  ! y(1 )=2.0034657220268395d0
  ! y(N1)=1.4982671389865767d0

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
subroutine get_bc_BiCGSTAB_p(N,N1,N2,NGC,bc_type,y)
integer, intent(in) :: N, N1, N2, NGC
character(len=*), intent(in) :: bc_type
double precision, intent(inout) :: y(N)
integer :: i
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if     (bc_type=='periodic') then
  call set_GC_to_periodic(N,N1,N2,NGC,y)
else if (bc_type=='zero' .or. bc_type=='usr') then
  call set_GC_to_x0(N,N1,N2,NGC,0.d0,y)
endif
end subroutine get_bc_BiCGSTAB_p
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Set physical boundary conditions for scalar fields s
!> in the BiCGSTAB method
!>
!> @warning for non-periodic functions, VITAL to set s BCs to 0
! -----------------------------------------------------------------------------
subroutine get_bc_BiCGSTAB_s(N,N1,N2,NGC,bc_type,y)
integer, intent(in) :: N, N1, N2, NGC
character(len=*), intent(in) :: bc_type
double precision, intent(inout) :: y(N)
integer :: i
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if     (bc_type=='periodic') then
  call set_GC_to_periodic(N,N1,N2,NGC,y)
else if (bc_type=='zero' .or. bc_type=='usr') then
  call set_GC_to_x0(N,N1,N2,NGC,0.d0,y)
endif
end subroutine get_bc_BiCGSTAB_s
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Set periodic BCs on 1D array y
! -----------------------------------------------------------------------------
subroutine set_GC_to_periodic(N,N1,N2,NGC,y)
integer, intent(in) :: N, N1, N2, NGC
double precision, intent(inout) :: y(N)
integer :: p, k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Bottom
do p=1,NGC
do k=(p-1)*N1+NGC+1,p*N1-NGC
  y(k)=y(k+N1*(N2-2*NGC))
  ! print*, 'bottom', k, k+N1*(N2-2*NGC)
enddo
enddo
! Top
do p=N2-NGC+1,N2
do k=(p-1)*N1+NGC+1,p*N1-NGC
  y(k)=y(k-N1*(N2-2*NGC))
  ! print*, 'top', k, k-N1*(N2-2*NGC)
enddo
enddo
! Left
do p=NGC,N2-NGC-1
do k=p*N1+1,p*N1+NGC
  y(k)=y(k+N1-2*NGC)
  ! print*, 'left', k, k+N1-2*NGC
enddo
enddo
! Right
do p=NGC,N2-NGC-1
do k=(p+1)*N1-NGC+1,(p+1)*N1
  y(k)=y(k-N1+2*NGC)
  ! print*, 'right', k, k-N1+2*NGC
enddo
enddo
end subroutine set_GC_to_periodic
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Set all the ghost cells, including corners, in a 1D array x to value x0
! -----------------------------------------------------------------------------
subroutine set_GC_to_x0(N,N1,N2,NGC,x0,x)
integer, intent(in) :: N, N1, N2, NGC
double precision, intent(in) :: x0
double precision, intent(inout) :: x(N)
integer :: p, k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Bottom (w/ corners)
do p=1,NGC
do k=(p-1)*N1+1,p*N1
  x(k)=x0
enddo
enddo
! Top (w/ corners)
do p=N2-NGC+1,N2
do k=(p-1)*N1+1,p*N1
  x(k)=x0
enddo
enddo
! Left (w/o corners)
do p=NGC,N2-NGC-1
do k=p*N1+1,p*N1+NGC
  x(k)=x0
enddo
enddo
! Right (w/o corners)
do p=NGC,N2-NGC-1
do k=(p+1)*N1-NGC+1,(p+1)*N1
  x(k)=x0
enddo
enddo
end subroutine set_GC_to_x0
! -----------------------------------------------------------------------------

end module mod_bc
