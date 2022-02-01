module mod_bc

implicit none

public :: get_bc
public :: get_bc_BiCGSTAB_p
public :: get_bc_BiCGSTAB_s
public :: set_GC_to_x0

contains

! -----------------------------------------------------------------------------
!> Set physical boundary conditions for scalar field U (not in corners)
!>
!> @todo add new types of boundary conditions
!> @warning for non-periodic functions, VITAL to let BCs loose
! -----------------------------------------------------------------------------
subroutine get_bc(N,N1,N2,NGC,bc_type,x,U)
integer, intent(in) :: N, N1, N2, NGC
character(len=*), dimension(4), intent(in) :: bc_type
double precision, intent(in) :: x(N1,N2,2)
! double precision, intent(in) :: w(N1,N2)
double precision, intent(inout) :: U(N)
double precision :: dpi=dacos(-1.d0)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if        (bc_type(1)=='periodic') then
  call set_GC_to_periodic(N,N1,N2,NGC,'left-right',U)
else
  if      (bc_type(1)=='zero') then
    call set_GC_to_x0(N,N1,N2,NGC,0.d0,'left',U)
  else if (bc_type(1)=='usr') then
    call set_GC_to_profile(N,N1,N2,NGC,'left',x,U)
  else if (bc_type(1)=='cont') then
    call set_GC_to_cont(N,N1,N2,NGC,'left',U)
  endif
  if      (bc_type(2)=='zero') then
    call set_GC_to_x0(N,N1,N2,NGC,0.d0,'right',U)
  else if (bc_type(2)=='usr') then
    call set_GC_to_profile(N,N1,N2,NGC,'right',x,U)
  else if (bc_type(2)=='cont') then
    call set_GC_to_cont(N,N1,N2,NGC,'right',U)
  endif
endif

if        (bc_type(3)=='periodic') then
  call set_GC_to_periodic(N,N1,N2,NGC,'bottom-up',U)
else
  if      (bc_type(3)=='zero') then
    call set_GC_to_x0(N,N1,N2,NGC,0.d0,'bottom',U)
  else if (bc_type(3)=='usr') then
    call set_GC_to_profile(N,N1,N2,NGC,'bottom',x,U)
  else if (bc_type(3)=='cont') then
    call set_GC_to_cont(N,N1,N2,NGC,'bottom',U)
  endif
  if      (bc_type(4)=='zero') then
    call set_GC_to_x0(N,N1,N2,NGC,0.d0,'up',U)
  else if (bc_type(4)=='usr') then
    call set_GC_to_profile(N,N1,N2,NGC,'up',x,U)
  else if (bc_type(4)=='cont') then
    call set_GC_to_cont(N,N1,N2,NGC,'up',U)
  endif
endif

end subroutine get_bc
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Set physical boundary conditions for scalar fields p
!> in the BiCGSTAB method
!>
!> @warning for non-periodic functions, VITAL to set p BCs to 0
! -----------------------------------------------------------------------------
subroutine get_bc_BiCGSTAB_p(N,N1,N2,NGC,bc_type,U)
integer, intent(in) :: N, N1, N2, NGC
character(len=*), dimension(4), intent(in) :: bc_type
double precision, intent(inout) :: U(N)
integer :: i
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if     (bc_type(1)=='periodic') then
  call set_GC_to_periodic(N,N1,N2,NGC,'all',U)
else
  call set_GC_to_x0(N,N1,N2,NGC,0.d0,'all',U)
endif
end subroutine get_bc_BiCGSTAB_p
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Set physical boundary conditions for scalar fields s
!> in the BiCGSTAB method
!>
!> @warning for non-periodic functions, VITAL to set s BCs to 0
! -----------------------------------------------------------------------------
subroutine get_bc_BiCGSTAB_s(N,N1,N2,NGC,bc_type,U)
integer, intent(in) :: N, N1, N2, NGC
character(len=*), dimension(4), intent(in) :: bc_type
double precision, intent(inout) :: U(N)
integer :: i
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if     (bc_type(1)=='periodic') then
  call set_GC_to_periodic(N,N1,N2,NGC,'all',U)
else
  call set_GC_to_x0(N,N1,N2,NGC,0.d0,'all',U)
endif
end subroutine get_bc_BiCGSTAB_s
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Set periodic BCs on 1D array U in direction 1 or 2 (or both for 'all')
!>
!> @note direction is either 'left-right' or 'bottom-up'
! -----------------------------------------------------------------------------
subroutine set_GC_to_periodic(N,N1,N2,NGC,direction,U)
integer, intent(in) :: N, N1, N2, NGC
character(len=*), intent(in) :: direction
double precision, intent(inout) :: U(N)
integer :: p, k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if      (direction=='left-right') then
  ! Left
  do p=NGC,N2-NGC-1
  do k=p*N1+1,p*N1+NGC
    U(k)=U(k+N1-2*NGC)
  enddo
  enddo
  ! Right
  do p=NGC,N2-NGC-1
  do k=(p+1)*N1-NGC+1,(p+1)*N1
    U(k)=U(k-N1+2*NGC)
  enddo
  enddo
else if (direction=='bottom-up') then
  ! Bottom
  do p=1,NGC
  do k=(p-1)*N1+NGC+1,p*N1-NGC
    U(k)=U(k+N1*(N2-2*NGC))
  enddo
  enddo
  ! Up
  do p=N2-NGC+1,N2
  do k=(p-1)*N1+NGC+1,p*N1-NGC
    U(k)=U(k-N1*(N2-2*NGC))
  enddo
  enddo
else if (direction=='all') then
  ! Left
  do p=NGC,N2-NGC-1
  do k=p*N1+1,p*N1+NGC
    U(k)=U(k+N1-2*NGC)
  enddo
  enddo
  ! Right
  do p=NGC,N2-NGC-1
  do k=(p+1)*N1-NGC+1,(p+1)*N1
    U(k)=U(k-N1+2*NGC)
  enddo
  enddo
  ! Bottom
  do p=1,NGC
  do k=(p-1)*N1+NGC+1,p*N1-NGC
    U(k)=U(k+N1*(N2-2*NGC))
  enddo
  enddo
  ! Up
  do p=N2-NGC+1,N2
  do k=(p-1)*N1+NGC+1,p*N1-NGC
    U(k)=U(k-N1*(N2-2*NGC))
  enddo
  enddo
endif

end subroutine set_GC_to_periodic
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Set all the ghost cells on a given side (or all ghost cells including corners)
!> in a 1D array U to value x0
! -----------------------------------------------------------------------------
subroutine set_GC_to_x0(N,N1,N2,NGC,x0,side,U)
integer, intent(in) :: N, N1, N2, NGC
double precision, intent(in) :: x0
character(len=*), intent(in) :: side
double precision, intent(inout) :: U(N)
integer :: p, k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if      (side=='left') then
  do p=NGC,N2-NGC-1
  do k=p*N1+1,p*N1+NGC
    U(k)=x0
  enddo
  enddo
else if (side=='right') then
  do p=NGC,N2-NGC-1
  do k=(p+1)*N1-NGC+1,(p+1)*N1
    U(k)=x0
  enddo
  enddo
else if (side=='bottom') then
  do p=1,NGC
  do k=(p-1)*N1+NGC+1,p*N1-NGC
    U(k)=x0
  enddo
  enddo
else if (side=='up') then
  do p=N2-NGC+1,N2
  do k=(p-1)*N1+NGC+1,p*N1-NGC
    U(k)=x0
  enddo
  enddo
else if (side=='all') then
  ! Left
  do p=NGC,N2-NGC-1
  do k=p*N1+1,p*N1+NGC
    U(k)=x0
  enddo
  enddo
  ! Right
  do p=NGC,N2-NGC-1
  do k=(p+1)*N1-NGC+1,(p+1)*N1
    U(k)=x0
  enddo
  enddo
  ! Bottom
  do p=1,NGC
  do k=(p-1)*N1+1,p*N1
    U(k)=x0
  enddo
  enddo
  ! Up
  do p=N2-NGC+1,N2
  do k=(p-1)*N1+1,p*N1
    U(k)=x0
  enddo
  enddo
endif

end subroutine set_GC_to_x0
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Set all the ghost cells on a given side to a pre-defined profile
!> in a 1D array U
! -----------------------------------------------------------------------------
subroutine set_GC_to_profile(N,N1,N2,NGC,side,x,U)
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
character(len=*), intent(in) :: side
double precision, intent(in) :: x(N1,N2,2)
double precision, intent(inout) :: U(N)
integer :: i, j, p, k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
analytic => analytic_

if      (side=='left') then
  do p=NGC,N2-NGC-1
  do k=p*N1+1,p*N1+NGC
    j=int((k-1)/N1)+1
    i=k-(j-1)*N1
    U(k)=analytic(x(i,j,1),x(i,j,2))
  enddo
  enddo
else if (side=='right') then
  do p=NGC,N2-NGC-1
  do k=(p+1)*N1-NGC+1,(p+1)*N1
    j=int((k-1)/N1)+1
    i=k-(j-1)*N1
    U(k)=analytic(x(i,j,1),x(i,j,2))
  enddo
  enddo
else if (side=='bottom') then
  do p=1,NGC
  do k=(p-1)*N1+NGC+1,p*N1-NGC
    j=int((k-1)/N1)+1
    i=k-(j-1)*N1
    U(k)=analytic(x(i,j,1),x(i,j,2))
  enddo
  enddo
else if (side=='up') then
  do p=N2-NGC+1,N2
  do k=(p-1)*N1+NGC+1,p*N1-NGC
    j=int((k-1)/N1)+1
    i=k-(j-1)*N1
    U(k)=analytic(x(i,j,1),x(i,j,2))
  enddo
  enddo
endif

end subroutine set_GC_to_profile
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Set all the ghost cells on a given side to the value just inside
!> the simulation space
! -----------------------------------------------------------------------------
subroutine set_GC_to_cont(N,N1,N2,NGC,side,U)
integer, intent(in) :: N, N1, N2, NGC
character(len=*), intent(in) :: side
double precision, intent(inout) :: U(N)
integer :: p, k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if      (side=='left') then
  do p=NGC,N2-NGC-1
  do k=p*N1+1,p*N1+NGC
    U(k)=U(p*N1+NGC+1)
  enddo
  enddo
else if (side=='right') then
  do p=NGC,N2-NGC-1
  do k=(p+1)*N1-NGC+1,(p+1)*N1
    U(k)=U(p*N1+N1-NGC)
  enddo
  enddo
else if (side=='bottom') then
  do p=1,NGC
  do k=(p-1)*N1+NGC+1,p*N1-NGC
    U(k)=U(NGC*N1+k-(int((k-1)/N1)+1-1)*N1)
  enddo
  enddo
else if (side=='up') then
  do p=N2-NGC+1,N2
  do k=(p-1)*N1+NGC+1,p*N1-NGC
    U(k)=U(N-(NGC+1)*N1+k-(int((k-1)/N1)+1-1)*N1)
  enddo
  enddo
endif

end subroutine set_GC_to_cont
! -----------------------------------------------------------------------------

end module mod_bc
