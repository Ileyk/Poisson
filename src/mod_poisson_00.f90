module mod_poisson
use mod_io

implicit none

public :: frst_gss
public :: gss_sdl

contains

! -----------------------------------------------------------------------------
!> First guess for potential f
! -----------------------------------------------------------------------------
subroutine frst_gss(N1,min1,max1,x,w,f)
integer, intent(in) :: N1
double precision, intent(in) :: min1, max1
double precision, intent(in) :: x(N1), w(N1)
double precision, allocatable, intent(out) :: f(:)
integer :: i
double precision :: dpi=dacos(-1.d0)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

allocate(f(N1))
do i=1,N1
  f(i)=-(1./40.)*dsin(2.d0*dpi*((x(i)-min1)/(max1-min1)))
enddo

end subroutine frst_gss
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Solver using Gauss-Seidel method
!> Find f given w (e.g. find grav potential given density mass)
! -----------------------------------------------------------------------------
subroutine gss_sdl(N1,min1,max1,x,w,f)
integer, intent(in) :: N1
double precision, intent(in) :: min1, max1
double precision, intent(in) :: x(N1), w(N1)
double precision, intent(inout) :: f(N1)
double precision, parameter :: eps0=1.d-5
integer :: i, it, Nitmax=10000
double precision :: eps, d1, C1, C2, denom
double precision, allocatable :: fth(:), fold(:), c(:), T(:,:)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Uniform grid
d1=x(2)-x(1)

call get_analytic(N1,min1,max1,x,fth)

it=1
eps=1.d99

! Compute vector c based on d1 and rho
! and matrix T which stands for Laplacian operator
allocate(c(N1),T(N1,N1))
T=0.d0
c(1)=-0.5d0*d1**2.d0*w(1)
do i=2,N1-1
  c(i)=-0.5d0*d1**2.d0*w(i)
  T(i,i-1)=0.5d0
  T(i,i+1)=0.5d0
enddo
c(N1)=-0.5d0*d1**2.d0*w(N1)
! These 2 ones correspond to periodic BCs
T(1,N1)=0.5d0
T(N1,1)=0.5d0
! T(1,2)=2./5.
! T(1,3)=-3./5.
! T(2,1)=1./3.
! T(2,3)=-1./9.
! T(3,1)=2./7.
! T(3,2)=-1./7.
! c(1)  =-1./5.
! c(2)  =2./9.
! c(3)  =-3./7.
fold=f

do while (eps>eps0 .and. it<Nitmax)

  f(1)=c(1)+T(1,N1)*fold(N1)       +T(1,2)*fold(2)
  do i=2,N1-1
    f(i)=c(i)+T(i,i-1)*fold(i-1)   +T(i,i+1)*fold(i+1)
  enddo
  f(N1)=c(N1)+T(N1,N1-1)*fold(N1-1)+T(N1,1)*fold(1)

  ! f(1)=c(1)+T(1,2)*fold(2)+T(1,3)*fold(3)
  ! f(2)=c(2)+T(2,1)*fold(1)+T(2,3)*fold(3)
  ! f(3)=c(3)+T(3,1)*fold(1)+T(3,2)*fold(2)

  ! call bcs(N1,f)

  call save_vec(N1,it,'f_',f)

  call get_residual(N1,fold,f,eps)
  print*, it,eps

  ! print*, f(1), f(2), f(3)

  it=it+1
  fold=f

enddo

end subroutine gss_sdl
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Specify analytic solution to be compared to by hand
! -----------------------------------------------------------------------------
subroutine get_analytic(N1,min1,max1,x,fth)
integer, intent(in) :: N1
double precision, intent(in) :: min1, max1
double precision, intent(in) :: x(N1)
double precision, allocatable, intent(out) :: fth(:)
integer :: i
double precision :: dpi=dacos(-1.d0)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

allocate(fth(N1))

do i=1,N1
  fth(i)=-dsin(2.d0*dpi*((x(i)-min1)/(max1-min1)))
enddo

end subroutine get_analytic
! -----------------------------------------------------------------------------

! ! -----------------------------------------------------------------------------
! !> Periodic boundary conditions
! ! -----------------------------------------------------------------------------
! subroutine bcs(N1,f)
! integer, intent(in) :: N1
! double precision, intent(inout) :: f(N1)
! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! f(1) =XXX
! f(N1)=XXX
!
! end subroutine bcs
! ! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Compute residual as relative difference w.r.t. solution @ previous step,
!> based on quadratic norm
! -----------------------------------------------------------------------------
subroutine get_residual(N1,fold,f,eps)
integer, intent(in) :: N1
double precision, intent(in) :: fold(N1), f(N1)
double precision, intent(out) :: eps
integer :: i
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

eps=0.d0
do i=1,N1
  eps=eps+(f(i)-fold(i))**2./max(1.d-10,fold(i)**2.)
enddo
eps=dsqrt(eps)

end subroutine get_residual
! -----------------------------------------------------------------------------

end module mod_poisson
