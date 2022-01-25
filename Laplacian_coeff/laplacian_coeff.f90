! -----------------------------------------------------------------------------
!> This code is a side-kick for the main Poisson one, in src/.
!> Its goal is to verify pen-and-paper formulas for derivatives of order n
!> w/ accuracy pencil (i.e. closest neighbors NGC on each side)
!> on a 1D Cartesian but not necessarily uniform grid.
!> The goal is to check the pen-and-paper formula against the much more
!> computationally-demanding method used here and if it matches,
!> implement the pen-and-paper formula in the main Poisson program.
!>
!> Input: pencil, n and the formulas to be checked
!>
!> Coefficients can be found here for a uniform 1D Cartesian grid:
!> @see https://en.wikipedia.org/wiki/Finite_difference_coefficient
!> For a basic test, you can also use this website:
!> @see https://web.media.mit.edu/~crtaylor/calculator.html
!> For video tutorials to adapt to 2D, see
!> @see https://www.youtube.com/watch?v=i0f2DegLBN8&ab_channel=SandipMazumder
!> @see https://www.youtube.com/watch?v=zNBkvYr5CNA&ab_channel=EMPossible
!> @warning I tried to expand it to non-Cartesian by providing a space metric
!> @warning but doing so, we no longer deal w/ a Vandermonde matrix since
!> @warning the combinatory factor has to be put back on the LHS,
!> @warning where it comes from, but only on 1 row
!> @warning Update: no, I don't think it is needed `(see 1.2.1 of my course),
!> @warning we can still work w/ a Vandermonde matrix but on the RHS,
!> @warning coefficients must be multiplied according to the general formula
!> @warning of the Laplacian in curved metric`
! -----------------------------------------------------------------------------
program laplacian_coeff
implicit none
integer, parameter :: pencil=5 ! input
integer, parameter :: n=2 ! input
integer, parameter :: N0=10 ! fiducial, not to be changed
integer :: NGC, N1
double precision, allocatable :: x(:)
double precision, allocatable  :: DLU(:,:)
double precision, allocatable  :: coeff(:,:)
double precision :: lbd(pencil)
double precision :: tmp
double precision :: coeff_i(pencil)
! double precision :: metric_inv, det_root
integer :: i,j,k
! integer, parameter :: r_=1, t_=2, p_=3
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if (pencil/=3 .and. pencil/=5) call mpistop("pencil/=3,5 not implemented yet")

NGC=(pencil-1)/2
N1=N0+NGC*2

allocate(x(N1),DLU(N1,N1),coeff(N1,pencil))
x=0.d0
DLU=0.d0

! Check beforehand that for a uniform grid, we retrieve the coefficients given on
! https://en.wikipedia.org/wiki/Finite_difference_coefficient
call make_grid(N1,NGC,1,N1,1.d0,2.d0,x,'uniform')
call get_operator(N1,pencil,n,x,DLU)
print*, ' '
print*, '---------------'
print*, 'Uniform grid - classic coefficients (see @)'
print*, '---------------'
print*, ' '
i=NGC+1
do j=i-NGC,i+NGC ! Plot only non-zero terms
  print*, DLU(i,j)*(x(2)-x(1))**dble(n)
enddo
print*, ' '

! print*, '---------------'
! print*, ' '

! Then, compute for a fiducial stretched grid
call make_grid(N1,NGC,1,N1,1.d0,2.d0,x,'stretched')
call get_operator(N1,pencil,n,x,DLU)

! do i=NGC+1,N1-NGC
! do j=i-NGC,i+NGC
!   print*, DLU(i,j)
! enddo
! print*, ' '
! enddo

! print*, '---------------'
! print*, ' '

do i=NGC+1,N1-NGC
  call get_lbd(N1,pencil,x,i,lbd)
  call get_coeff(pencil,n,lbd,coeff_i)
  coeff(i,:)=coeff_i
enddo

! Check the match between coeff and the non-zero terms of DLU
tmp=0.d0
do i=NGC+1,N1-NGC
  do k=1,pencil
    j=i-(pencil-1)/2+(k-1)
    tmp=max(tmp,dabs(DLU(i,j)-coeff(i,k)))
    ! print*, i,j,k,DLU(i,j), coeff(i,k), dabs(DLU(i,j)-coeff(i,k))
  enddo
  ! print*, ' '
enddo

if (tmp<1.d-12) then
  print*, 'Congratulations, formula validated'
else
  print*, 'Formula wrong', tmp
endif

deallocate(x,DLU,coeff)

end program laplacian_coeff

! -----------------------------------------------------------------------------
!> Define p-th order operator w/ pencil made of (p-1)/2
!> neighbors on each side + itself
!> for derivative of order n
!>
!> @param pencil is 3 (2 closest neighbours), 5 (4 closest neighbours), etc
!> @see https://en.wikipedia.org/wiki/Finite_difference_coefficient
!> @see http://www.iaeng.org/publication/WCE2017/WCE2017_pp127-129.pdf
!> @warning Laplacian is not defined in ghost cells
!> @warning => 1st and last NGC rows are 0.
! -----------------------------------------------------------------------------
subroutine get_operator(N1,pencil,n,x,DLU)
use mod_linal
integer, intent(in) :: N1
integer, intent(in) :: pencil
integer, intent(in) :: n
double precision, intent(in) :: x(N1)
double precision, intent(out) :: DLU(N1,N1)
double precision :: lbd(pencil)
double precision :: W(pencil,pencil)
double precision :: A(pencil,pencil)
double precision :: V_inv(pencil,pencil)
double precision :: deriv(pencil)
double precision :: coeff(pencil)
double precision :: tmp
integer :: i, j, NGC
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

NGC=(pencil-1)/2

DLU=0.d0

call deriv_order_n(pencil,n,deriv)

do i=NGC+1,N1-NGC

  call get_lbd(N1,pencil,x,i,lbd)

  call get_A(pencil,lbd,A)

  call get_W(pencil,lbd,W)

  call mx_x_mx(pencil,W,A,V_inv)

  call mx_x_vec(pencil,1,pencil,V_inv,deriv,coeff)

  do j=1,pencil
    DLU(i,i+(j-1-NGC))=coeff(j)
  enddo

enddo

end subroutine get_operator
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Get the vector which was multiplied by the inverse of the Vandermonde
!> matrix will yield the coefficients of the discretized form of the
!> derivative of order n w/ a given pencil
!>
!> @note pencil p has to be such that n < or = 3 + 2*int((n-1)/2)
!> @param n is order of derivative we want: f' (n=1), f'' (n=2), etc
! -----------------------------------------------------------------------------
subroutine deriv_order_n(pencil,n,deriv)
integer, intent(in) :: pencil
integer, intent(in) :: n
double precision, intent(out) :: deriv(pencil)
integer :: i
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

deriv=0.d0
deriv(n+1)=1.d0
! Factorial of n
do i=2,n
  deriv(n+1)=deriv(n+1)*i
enddo

end subroutine deriv_order_n
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

! -----------------------------------------------------------------------------
!> Compute the elements of matrix A given in
!> http://www.iaeng.org/publication/WCE2017/WCE2017_pp127-129.pdf
!> such as W*A is the invert of the Vandermonde matrix
!>
!> @todo Extend formulas for A(i,1) to any pencil instead of if branching
! -----------------------------------------------------------------------------
subroutine get_A(pencil,lbd,A)
integer, intent(in) :: pencil
double precision, intent(in) :: lbd(pencil)
double precision, intent(out) :: A(pencil,pencil)
double precision :: tmp
integer :: i, j, k, l
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

A=0.d0

if      (pencil==3) then

  A(2,1)=-( lbd(1) + lbd(2) + lbd(3) )
  A(3,1)= ( lbd(1)*lbd(2) + lbd(1)*lbd(3) + &
            lbd(2)*lbd(3) )

else if (pencil==5) then

  A(2,1)=-( lbd(1) + lbd(2) + lbd(3) + lbd(4) + lbd(5) )
  A(3,1)= ( lbd(1)*lbd(2) + lbd(1)*lbd(3) + lbd(1)*lbd(4) + lbd(1)*lbd(5) + &
            lbd(2)*lbd(3) + lbd(2)*lbd(4) + lbd(2)*lbd(5) + &
            lbd(3)*lbd(4) + lbd(3)*lbd(5) + &
            lbd(4)*lbd(5) )
  A(4,1)=-( lbd(1)*lbd(2)*lbd(3) + lbd(1)*lbd(2)*lbd(4) + lbd(1)*lbd(2)*lbd(5) + &
            lbd(1)*lbd(3)*lbd(4) + lbd(1)*lbd(3)*lbd(5) + lbd(1)*lbd(4)*lbd(5) + &
            lbd(2)*lbd(3)*lbd(4) + lbd(2)*lbd(3)*lbd(5) + &
            lbd(2)*lbd(4)*lbd(5) + &
            lbd(3)*lbd(4)*lbd(5) )
  A(5,1)= ( lbd(1)*lbd(2)*lbd(3)*lbd(4) + lbd(1)*lbd(2)*lbd(3)*lbd(5) + lbd(1)*lbd(2)*lbd(4)*lbd(5) + &
            lbd(2)*lbd(3)*lbd(4)*lbd(5) + lbd(1)*lbd(3)*lbd(4)*lbd(5))

endif

do i=1,pencil
  do j=2,i-1
    A(i,j)=A(i-1,j-1)
  enddo
  A(i,i)=1.d0
enddo

end subroutine get_A
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Compute the elements of matrix W given in
!> http://www.iaeng.org/publication/WCE2017/WCE2017_pp127-129.pdf
!> such as W*A is the invert of the Vandermonde matrix
! -----------------------------------------------------------------------------
subroutine get_W(pencil,lbd,W)
integer, intent(in) :: pencil
double precision, intent(in) :: lbd(pencil)
double precision, intent(out) :: W(pencil,pencil)
double precision :: tmp
integer :: i, j
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

do i=1,pencil
  tmp=1.d0
  do j=1,pencil
    if (j/=i) tmp=tmp*(lbd(i)-lbd(j))
  enddo
  do j=1,pencil
    W(i,j)=lbd(i)**(dble(pencil-j))/tmp
  enddo
enddo

end subroutine get_W
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Make a fiducial 1D Cartesian mesh, just to compute the
!> lambdas and by then, the invert of the Vandermonde matrix
! -----------------------------------------------------------------------------
subroutine make_grid(N1,NGC,iImin1,iImax1,min1,max1,x,type)
integer, intent(in) :: N1
integer, intent(in) :: NGC
integer, intent(in) :: iImin1, iImax1
double precision, intent(in) :: min1, max1
character(len=*), intent(in) :: type
double precision, intent(out) :: x(N1)
integer :: i
double precision :: dx
double precision :: q
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if (type=='uniform') then

  dx=(max1-min1)/dble(N1-2*NGC)
  do i=iImin1,iImax1
    x(i)=min1+dx*(i-0.5-NGC)
  enddo

else if (type=='stretched') then

  q=(max1/min1)**(1.d0/dble(N1-2*NGC))
  ! print*, q
  ! stop
  dx=min1*(q-1.d0)/q
  x(1)=min1-0.5d0*dx
  ! dxold=dx/q
  do i=iImin1+1,iImax1
    x(i)=x(i-1)+dx*0.5d0+dx*q*0.5d0
    dx=dx*q
    ! dxold=dx
  enddo
  ! x(iImax1)=max1+0.5d0*dx*q
  ! stop

endif

end subroutine make_grid
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> HERE IS THE MAIN INPUT, THE FORMULAS DERIVED PEN-AND-PAPER AND SIMPLER
!> THAN THE FULL COMPUTATION
!>
!> @todo Fill more pen-and-paper derived formulas
! -----------------------------------------------------------------------------
subroutine get_coeff(pencil,n,lbd,coeff)
integer, intent(in) :: pencil
integer, intent(in) :: n
double precision, intent(in) :: lbd(pencil)
double precision, intent(out) :: coeff(pencil)
double precision :: a1, a2
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if (pencil==3) then
  if      (n==1) then
    coeff(1)= 1.d0/(lbd(1)-lbd(3))-(lbd(1)+lbd(3))/(lbd(1)*(lbd(1)-lbd(3)))
    coeff(2)=-(lbd(1)+lbd(3))/(lbd(1)*lbd(3))
    coeff(3)= 1.d0/(lbd(3)-lbd(1))-(lbd(1)+lbd(3))/(lbd(3)*(lbd(3)-lbd(1)))
  else if (n==2) then
    coeff(1)=2.d0/(lbd(1)*(lbd(1)-lbd(3)))
    coeff(2)=2.d0/(lbd(1)*lbd(3))
    coeff(3)=2.d0/(lbd(3)*(lbd(3)-lbd(1)))
  endif
else if (pencil==5) then
  if      (n==1) then
    print*, 'TO BE DERIVED'
    stop
  else if (n==2) then
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
endif

end subroutine get_coeff
! -----------------------------------------------------------------------------
