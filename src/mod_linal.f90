!> Linear algebra operations such as matrix products

module mod_linal

implicit none

public :: mx_x_mx
public :: mx_x_vec
public :: dot_pdct

contains

! -----------------------------------------------------------------------------
!> Dot product between 2 vectors of size N1
! -----------------------------------------------------------------------------
subroutine dot_pdct(N,N1,N2,NGC,a,b,ab)
integer, intent(in) :: N, N1, N2, NGC
double precision, intent(in) :: a(N), b(N)
double precision, intent(out) :: ab
integer :: p, k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ab=0.d0
! do p=1,N2-NGC-2
! do k=(p-1+NGC)*N1+NGC+1,(p-1+NGC)*N1+N1-NGC
do p=NGC,N2-NGC-1
do k=p*N1+NGC+1,p*N1+N1-NGC
  ab=ab+a(k)*b(k)
enddo
enddo

end subroutine dot_pdct
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Product between 2 matrix of size N1xN1.
! -----------------------------------------------------------------------------
subroutine mx_x_mx(N1,A,B,AB)
integer, intent(in) :: N1
double precision, intent(in) :: A(N1,N1), B(N1,N1)
double precision, intent(out) :: AB(N1,N1)
integer :: i, j, k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

AB=0.d0
do i=1,N1
  do j=1,N1
    do k=1,N1
      AB(i,j)=AB(i,j)+A(i,k)*B(k,j)
    enddo
  enddo
enddo

end subroutine mx_x_mx
! -----------------------------------------------------------------------------


! -----------------------------------------------------------------------------
!> Product between matrix of size N1xN1 and vector of size N1.
!>
!> @warning Perform the computation only between indices iOmin1 & iOmax1,
!> @warning not in the ghost cells.
! -----------------------------------------------------------------------------
subroutine mx_x_vec(N,N1,N2,NGC,A,x,y)
integer, intent(in) :: N, N1, N2, NGC
double precision, intent(in) :: A(N,N), x(N)
double precision, intent(out) :: y(N)
integer :: i, j, p, k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

y=0.d0
do p=NGC,N2-NGC-1
do k=p*N1+NGC+1,p*N1+N1-NGC
  do j=1,N
    y(k)=y(k)+A(k,j)*x(j)
  enddo
enddo
enddo

end subroutine mx_x_vec
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Product between matrix of size N1xN1 and vector of size N1,
!> w/ a pencil of 3 i.e. a mx which is tridiagonal
! -----------------------------------------------------------------------------
subroutine mx_x_vec_tri_d(N1,A,x,y)
integer, intent(in) :: N1
double precision, intent(in) :: A(N1,N1), x(N1)
double precision, intent(out) :: y(N1)
integer :: i, j
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

y=0.d0
do i=2,N1-1
  do j=i-1,i+1
    y(i)=y(i)+A(i,j)*x(j)
  enddo
enddo

end subroutine mx_x_vec_tri_d
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Product between matrix of size N1xN1 and vector of size N1,
!> w/ a pencil of 3 i.e. a mx which is tridiagonal
!> but has also value in the 2 anti-corners
! -----------------------------------------------------------------------------
subroutine mx_x_vec_tri_d_corners(N1,A,x,y)
integer, intent(in) :: N1
double precision, intent(in) :: A(N1,N1), x(N1)
double precision, intent(out) :: y(N1)
integer :: i, j
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

y=0.d0

y(1)=A(1,N1)*x(N1)+A(1,1)*x(1)+A(1,2)*x(2)
do i=2,N1-1
  do j=i-1,i+1
    y(i)=y(i)+A(i,j)*x(j)
  enddo
enddo
y(N1)=A(N1,N1-1)*x(N1-1)+A(N1,N1)*x(N1)+A(N1,1)*x(1)

end subroutine mx_x_vec_tri_d_corners
! -----------------------------------------------------------------------------


end module mod_linal
