module mod_laplacian

implicit none

public :: get_laplacian

contains

! -----------------------------------------------------------------------------
!> Define p-th order Laplacian operator w/ pencil made of (p-1)/2
!> neighbors on each side + itself
!>
!> @param pencil is 3 (2 closest neighbours), 5 (4 closest neighbours), etc
!> @see https://en.wikipedia.org/wiki/Finite_difference_coefficient
!> @see http://www.iaeng.org/publication/WCE2017/WCE2017_pp127-129.pdf
!> @warning Laplacian is not defined in ghost cells
!> @warning => 1st and last NGC rows are 0.
!> @warning For now, valid for 2D Cartesian, spherical or any metric
!> @warning WITHOUT non-diagonal (i.e. crossed) terms
!> @todo Account for non-diagonal (i.e. crossed) terms in metric
!> @todo Derive formula for higher order Laplace operator
! -----------------------------------------------------------------------------
subroutine get_laplacian(N,N1,N2,pencil,x,DLU)
use mod_linal
use mod_metric
use mod_csts, only : dpi
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Functions of r and theta
abstract interface
  function f_rt (r,t)
     double precision :: f_rt
     double precision, intent (in) :: r, t
  end function f_rt
end interface
procedure (f_rt), pointer :: det_root => null ()
procedure (f_rt), pointer :: det_root_dr => null ()
procedure (f_rt), pointer :: det_root_dt => null ()
procedure (f_rt), pointer :: rr => null ()
procedure (f_rt), pointer :: rr_dr => null ()
procedure (f_rt), pointer :: tt => null ()
procedure (f_rt), pointer :: tt_dt => null ()
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
integer, intent(in) :: N, N1, N2
integer, intent(in) :: pencil
double precision, intent(in) :: x(N1,N2,2)
double precision, intent(out) :: DLU(N,N)
double precision :: lbd(2*(pencil-1)+1)
double precision :: coeff(2*(pencil-1)+1)
double precision :: r, t
integer :: i, j, p, k, NGC
double precision, parameter :: th0=dpi/2.d0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

det_root => det_root_
det_root_dr => det_root_dr_
det_root_dt => det_root_dt_
rr => rr_
rr_dr => rr_dr_
tt => tt_
tt_dt => tt_dt_

NGC=(pencil-1)/2

DLU=0.d0

do p=NGC,N2-NGC-1
do k=p*N1+NGC+1,p*N1+N1-NGC
  j=int((k-1)/N1)+1
  i=k-(j-1)*N1
  call get_lbd(N1,N2,NGC,pencil,x,i,j,lbd)
  r=x(i,j,1)
  t=x(i,j,2)
  if (pencil==3) then
    ! Neighbours in direction 1
    coeff(1)= ( 2.d0*rr(r,t) - lbd(3) * &
              (det_root_dr(r,t)*rr(r,t)/det_root(r,t)+rr_dr(r,t)) ) / &
              (lbd(1)*(lbd(1)-lbd(3)))
    coeff(3)= ( 2.d0*rr(r,t) - lbd(1) * &
              (det_root_dr(r,t)*rr(r,t)/det_root(r,t)+rr_dr(r,t)) ) / &
              (lbd(3)*(lbd(3)-lbd(1)))
    ! Neighbours in direction 2
    coeff(4)= ( 2.d0*tt(r,t) - lbd(5) * &
              (det_root_dt(r,t)*tt(r,t)/det_root(r,t)+tt_dt(r,t)) ) / &
              (lbd(4)*(lbd(4)-lbd(5)))
    coeff(5)= ( 2.d0*tt(r,t) - lbd(4) * &
              (det_root_dt(r,t)*tt(r,t)/det_root(r,t)+tt_dt(r,t)) ) / &
              (lbd(5)*(lbd(5)-lbd(4)))
    ! Self-coefficient
    coeff(2)=-coeff(1)-coeff(3)-coeff(4)-coeff(5)
  else if (pencil==5) then
    ! a1=-(lbd(1)+lbd(2)+lbd(4)+lbd(5))
    ! a2=lbd(1)*lbd(2)+lbd(1)*lbd(4)+lbd(1)*lbd(5)+lbd(2)*lbd(4)+lbd(2)*lbd(5)+lbd(4)*lbd(5)
    ! coeff(1)=(lbd(1)+a1)/((lbd(1)-lbd(2))*(lbd(1)-lbd(4))*(lbd(1)-lbd(5)))+&
    !          a2/(lbd(1)*(lbd(1)-lbd(2))*(lbd(1)-lbd(4))*(lbd(1)-lbd(5)))
    ! coeff(2)=(lbd(2)+a1)/((lbd(2)-lbd(1))*(lbd(2)-lbd(4))*(lbd(2)-lbd(5)))+&
    !          a2/(lbd(2)*(lbd(2)-lbd(1))*(lbd(2)-lbd(4))*(lbd(2)-lbd(5)))
    ! coeff(3)=a2/(lbd(1)*lbd(2)*lbd(4)*lbd(5))
    ! coeff(4)=(lbd(4)+a1)/((lbd(4)-lbd(1))*(lbd(4)-lbd(2))*(lbd(4)-lbd(5)))+&
    !          a2/(lbd(4)*(lbd(4)-lbd(1))*(lbd(4)-lbd(2))*(lbd(4)-lbd(5)))
    ! coeff(5)=(lbd(5)+a1)/((lbd(5)-lbd(1))*(lbd(5)-lbd(2))*(lbd(5)-lbd(4)))+&
    !          a2/(lbd(5)*(lbd(5)-lbd(1))*(lbd(5)-lbd(2))*(lbd(5)-lbd(4)))
    ! coeff=2.d0*coeff
  endif
  ! 1st direction
  do j=1,pencil
    DLU(k,k+(j-1-NGC))=coeff(j)
  enddo
  ! 2nd direction
  i=pencil+1 ! to shift index
  do j=NGC,1,-1
    DLU(k,k-j*N1)=coeff(i)
    i=i+1
  enddo
  do j=1,NGC
    DLU(k,k+j*N1)=coeff(i)
    i=i+1
  enddo
enddo
enddo



! do i=NGC+1,N1-NGC
!
!   call get_lbd(N1,pencil,x,i,lbd)
!
!   ! To check this formula beforehand, use side-kick code
!   ! laplacian_coeff.f90
!   if (pencil==3) then
!
!     ! 1D Cartesian, spherical or any other metric along direction r
!     coeff(1)= ( rr(x(i)) - 0.5d0*lbd(3) * (det_root_dr(x(i),th0)*rr(x(i))/det_root(x(i),th0)+rr_dr(x(i))) ) / &
!               (lbd(1)*0.5d0*(lbd(1)-lbd(3)))
!     coeff(3)= ( rr(x(i)) - 0.5d0*lbd(1) * (det_root_dr(x(i),th0)*rr(x(i))/det_root(x(i),th0)+rr_dr(x(i))) ) / &
!               (lbd(3)*0.5d0*(lbd(3)-lbd(1)))
!     coeff(2)=-coeff(1)-coeff(3)
!
!   else if (pencil==5) then
!     a1=-(lbd(1)+lbd(2)+lbd(4)+lbd(5))
!     a2=lbd(1)*lbd(2)+lbd(1)*lbd(4)+lbd(1)*lbd(5)+lbd(2)*lbd(4)+lbd(2)*lbd(5)+lbd(4)*lbd(5)
!     coeff(1)=(lbd(1)+a1)/((lbd(1)-lbd(2))*(lbd(1)-lbd(4))*(lbd(1)-lbd(5)))+&
!              a2/(lbd(1)*(lbd(1)-lbd(2))*(lbd(1)-lbd(4))*(lbd(1)-lbd(5)))
!     coeff(2)=(lbd(2)+a1)/((lbd(2)-lbd(1))*(lbd(2)-lbd(4))*(lbd(2)-lbd(5)))+&
!              a2/(lbd(2)*(lbd(2)-lbd(1))*(lbd(2)-lbd(4))*(lbd(2)-lbd(5)))
!     coeff(3)=a2/(lbd(1)*lbd(2)*lbd(4)*lbd(5))
!     coeff(4)=(lbd(4)+a1)/((lbd(4)-lbd(1))*(lbd(4)-lbd(2))*(lbd(4)-lbd(5)))+&
!              a2/(lbd(4)*(lbd(4)-lbd(1))*(lbd(4)-lbd(2))*(lbd(4)-lbd(5)))
!     coeff(5)=(lbd(5)+a1)/((lbd(5)-lbd(1))*(lbd(5)-lbd(2))*(lbd(5)-lbd(4)))+&
!              a2/(lbd(5)*(lbd(5)-lbd(1))*(lbd(5)-lbd(2))*(lbd(5)-lbd(4)))
!     coeff=2.d0*coeff
!   endif
!
!   do j=1,pencil
!     DLU(i,i+(j-1-NGC))=coeff(j)
!   enddo
!
! enddo
!
! do p=NGC,N2-NGC-1
! do k=p*N1+NGC+1,p*N1+N1-NGC
!   ! DLU(k,k-N1)= 4.d0/3.d0
!   ! DLU(k,k+N1)= 4.d0/3.d0
!   ! DLU(k,k-1) = 4.d0/3.d0
!   ! DLU(k,k+1) = 4.d0/3.d0
!   ! DLU(k,k-2*N1)= -1.d0/12.d0
!   ! DLU(k,k+2*N1)= -1.d0/12.d0
!   ! DLU(k,k-2) = -1.d0/12.d0
!   ! DLU(k,k+2) = -1.d0/12.d0
!   ! DLU(k,k)   =-5.d0
!   DLU(k,k-N1)= 1.d0
!   DLU(k,k+N1)= 1.d0
!   DLU(k,k-1) = 1.d0
!   DLU(k,k+1) = 1.d0
!   DLU(k,k)   =-4.d0
! enddo
! enddo
! DLU=DLU/(1.d0/dble(N1-2*NGC))**2.d0

end subroutine get_laplacian
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Compute the algebraic distances x(i)-x(j)
!> for any neighbour j in [|i-(pencil-1)/2;i+(pencil-1)/2|]
!>
!> Valid for any pencil
!> Convention: left --> right --> bottom --> top
!> E.g., for a pencil of 5:
!>          9
!>          |
!>          8
!>          |
!>  1 - 2 - 3 - 4 - 5
!>          |
!>          7
!>          |
!>          6
!> @warning the 1st direction must be treated differently from others
!> @warning because it contains the central value (e.g. 3 above)
! -----------------------------------------------------------------------------
subroutine get_lbd(N1,N2,NGC,pencil,x,i0,j0,lbd)
integer, intent(in) :: N1, N2, NGC
integer, intent(in) :: pencil
double precision, intent(in) :: x(N1,N2,2)
integer, intent(in) :: i0, j0
double precision, intent(out) :: lbd(2*(pencil-1)+1)
integer :: i, j, p
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
lbd=0.d0
p=1
! 1st direction
do i=NGC,1,-1
  lbd(p)=x(i0-i,j0,1)-x(i0,j0,1)
  p=p+1
enddo
p=p+1 ! skip the central one
do i=1,NGC
  lbd(p)=x(i0+i,j0,1)-x(i0,j0,1)
  p=p+1
enddo
! 2nd direction
do j=NGC,1,-1
  lbd(p)=x(i0,j0-j,2)-x(i0,j0,2)
  p=p+1
enddo
do j=1,NGC
  lbd(p)=x(i0,j0+j,2)-x(i0,j0,2)
  p=p+1
enddo

end subroutine get_lbd
! -----------------------------------------------------------------------------

end module mod_laplacian
