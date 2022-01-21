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
! -----------------------------------------------------------------------------
subroutine get_laplacian(N1,pencil,x,DLU)
use mod_linal
integer, intent(in) :: N1
integer, intent(in) :: pencil
double precision, intent(in) :: x(N1)
double precision, intent(out) :: DLU(N1,N1)
double precision :: lbd(pencil)
double precision :: W(pencil,pencil)
double precision :: A(pencil,pencil)
double precision :: V_inv(pencil,pencil)
double precision :: deriv(pencil)
double precision :: coeff(pencil)
double precision :: tmp, a1, a2
integer :: i, j, NGC
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

NGC=(pencil-1)/2

DLU=0.d0

! call deriv_order_n(pencil,2,deriv)

do i=NGC+1,N1-NGC

  call get_lbd(N1,pencil,x,i,lbd)

  ! To check this formula beforehand, use side-kick code
  ! laplacian_coeff.f90
  if (pencil==3) then
      coeff(1)=2.d0/(lbd(1)*(lbd(1)-lbd(3)))
      coeff(2)=2.d0/(lbd(1)*lbd(3))
      coeff(3)=2.d0/(lbd(3)*(lbd(3)-lbd(1)))
  else if (pencil==5) then
    a1=-(lbd(1)+lbd(2)+lbd(4)+lbd(5))
    a2=lbd(1)*lbd(2)+lbd(1)*lbd(4)+lbd(1)*lbd(5)+lbd(2)*lbd(4)+lbd(2)*lbd(5)+lbd(4)*lbd(5)
    coeff(1)=(lbd(1)+a1)/((lbd(1)-lbd(2))*(lbd(1)-lbd(4))*(lbd(1)-lbd(5)))+&
             a2/(lbd(1)*(lbd(1)-lbd(2))*(lbd(1)-lbd(4))*(lbd(1)-lbd(5)))
    coeff(2)=(lbd(2)+a1)/((lbd(2)-lbd(1))*(lbd(2)-lbd(4))*(lbd(2)-lbd(5)))+&
             a2/(lbd(2)*(lbd(2)-lbd(1))*(lbd(2)-lbd(4))*(lbd(2)-lbd(5)))
    coeff(3)=a2/(lbd(1)*lbd(2)*lbd(4)*lbd(5))
    coeff(4)=(lbd(4)+a1)/((lbd(4)-lbd(1))*(lbd(4)-lbd(2))*(lbd(4)-lbd(5)))+&
             a2/(lbd(4)*(lbd(4)-lbd(1))*(lbd(4)-lbd(2))*(lbd(4)-lbd(5)))
    coeff(5)=(lbd(5)+a1)/((lbd(5)-lbd(1))*(lbd(5)-lbd(2))*(lbd(5)-lbd(4)))+&
             a2/(lbd(5)*(lbd(5)-lbd(1))*(lbd(5)-lbd(2))*(lbd(5)-lbd(4)))
    coeff=2.d0*coeff
  endif

  ! call get_A(pencil,lbd,A)
  !
  ! call get_W(pencil,lbd,W)
  !
  ! call mx_x_mx(pencil,W,A,V_inv)
  !
  ! call mx_x_vec(pencil,1,pencil,V_inv,deriv,coeff)

  do j=1,pencil
    DLU(i,i+(j-1-NGC))=coeff(j)
  enddo

enddo

! DLU(1,1)=-2.d0/(x(2)-x(1))**2.d0
! DLU(1,2)= 1.d0/(x(2)-x(1))**2.d0
! DLU(N1,N1-1)=1.d0/(x(2)-x(1))**2.d0
! DLU(N1,N1)=-2.d0/(x(2)-x(1))**2.d0

! DLU(1,2)=1.d0/(x(2)-x(1))**2.d0
! DLU(1,N1-2)= 1.d0/(x(2)-x(1))**2.d0
! DLU(1,N1-1)=-2.d0/(x(2)-x(1))**2.d0
! DLU(N1,3)=1.d0/(x(2)-x(1))**2.d0
! DLU(N1,2)=-2.d0/(x(2)-x(1))**2.d0
! DLU(N1,N1-1)=1.d0/(x(2)-x(1))**2.d0


end subroutine get_laplacian
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Compute the algebraic distances x(i)-x(j)
!> for any neighbour j in [|i-(pencil-1)/2;i+(pencil-1)/2|]
!>
!> Valid for any pencil
! -----------------------------------------------------------------------------
subroutine get_lbd(N1,pencil,x,i0,lbd)
integer, intent(in) :: N1
integer, intent(in) :: pencil
double precision, intent(in) :: x(N1)
integer, intent(in) :: i0
double precision, intent(out) :: lbd(pencil)
integer :: i, j
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

lbd=0.d0
do i=1,pencil
  j=-(pencil-1)/2+(i-1)
  lbd(i)=x(i0+j)-x(i0)
enddo
end subroutine get_lbd

! -----------------------------------------------------------------------------


! ! -----------------------------------------------------------------------------
! !> Get the vector which was multiplied by the inverse of the Vandermonde
! !> matrix will yield the coefficients of the discretized form of the
! !> derivative of order n w/ a given pencil
! !>
! !> @note pencil p has to be such that n < or = 3 + 2*int((n-1)/2)
! !> @param n is order of derivative we want: f' (n=1), f'' (n=2), etc
! ! -----------------------------------------------------------------------------
! subroutine deriv_order_n(pencil,n,deriv)
! integer, intent(in) :: pencil
! integer, intent(in) :: n
! double precision, intent(out) :: deriv(pencil)
! integer :: i
! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! deriv=0.d0
! deriv(n+1)=1.d0
! ! Factorial of n
! do i=2,n
!   deriv(n+1)=deriv(n+1)*i
! enddo
!
! end subroutine deriv_order_n
! ! -----------------------------------------------------------------------------
!
! ! -----------------------------------------------------------------------------
! !> Compute the elements of matrix A given in
! !> http://www.iaeng.org/publication/WCE2017/WCE2017_pp127-129.pdf
! !> such as W*A is the invert of the Vandermonde matrix
! !>
! !> @todo Extend formulas for A(i,1) to any pencil instead of if branching
! ! -----------------------------------------------------------------------------
! subroutine get_A(pencil,lbd,A)
! integer, intent(in) :: pencil
! double precision, intent(in) :: lbd(pencil)
! double precision, intent(out) :: A(pencil,pencil)
! double precision :: tmp
! integer :: i, j, k, l
! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! A=0.d0
!
! if      (pencil==3) then
!
!   A(2,1)=-( lbd(1) + lbd(2) + lbd(3) )
!   A(3,1)= ( lbd(1)*lbd(2) + lbd(1)*lbd(3) + &
!             lbd(2)*lbd(3) )
!
! else if (pencil==5) then
!
!   A(2,1)=-( lbd(1) + lbd(2) + lbd(3) + lbd(4) + lbd(5) )
!   A(3,1)= ( lbd(1)*lbd(2) + lbd(1)*lbd(3) + lbd(1)*lbd(4) + lbd(1)*lbd(5) + &
!             lbd(2)*lbd(3) + lbd(2)*lbd(4) + lbd(2)*lbd(5) + &
!             lbd(3)*lbd(4) + lbd(3)*lbd(5) + &
!             lbd(4)*lbd(5) )
!   A(4,1)=-( lbd(1)*lbd(2)*lbd(3) + lbd(1)*lbd(2)*lbd(4) + lbd(1)*lbd(2)*lbd(5) + &
!             lbd(1)*lbd(3)*lbd(4) + lbd(1)*lbd(3)*lbd(5) + lbd(1)*lbd(4)*lbd(5) + &
!             lbd(2)*lbd(3)*lbd(4) + lbd(2)*lbd(3)*lbd(5) + &
!             lbd(2)*lbd(4)*lbd(5) + &
!             lbd(3)*lbd(4)*lbd(5) )
!   A(5,1)= ( lbd(1)*lbd(2)*lbd(3)*lbd(4) + lbd(1)*lbd(2)*lbd(3)*lbd(5) + lbd(1)*lbd(2)*lbd(4)*lbd(5) + &
!             lbd(2)*lbd(3)*lbd(4)*lbd(5) + lbd(1)*lbd(3)*lbd(4)*lbd(5))
!
! endif
!
! do i=1,pencil
!   do j=2,i-1
!     A(i,j)=A(i-1,j-1)
!   enddo
!   A(i,i)=1.d0
! enddo
!
! end subroutine get_A
! ! -----------------------------------------------------------------------------
!
! ! -----------------------------------------------------------------------------
! !> Compute the elements of matrix W given in
! !> http://www.iaeng.org/publication/WCE2017/WCE2017_pp127-129.pdf
! !> such as W*A is the invert of the Vandermonde matrix
! ! -----------------------------------------------------------------------------
! subroutine get_W(pencil,lbd,W)
! integer, intent(in) :: pencil
! double precision, intent(in) :: lbd(pencil)
! double precision, intent(out) :: W(pencil,pencil)
! double precision :: tmp
! integer :: i, j
! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! do i=1,pencil
!   tmp=1.d0
!   do j=1,pencil
!     if (j/=i) tmp=tmp*(lbd(i)-lbd(j))
!   enddo
!   do j=1,pencil
!     W(i,j)=lbd(i)**(dble(pencil-j))/tmp
!   enddo
! enddo
!
! end subroutine get_W
! ! -----------------------------------------------------------------------------

end module mod_laplacian
