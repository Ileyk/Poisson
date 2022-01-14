module mod_poisson
use mod_io
use mod_params
use mod_csts
use mod_mpi

implicit none

public :: frst_gss
public :: solver
public :: get_operator
public :: get_RHS

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
  f(i)=0.d0 ! -(1.d0/(2.d0*dpi*kwv/(max1-min1))**2.d0)*dsin(2.d0*dpi*kwv*((x(i)-min1)/(max1-min1)))
enddo

end subroutine frst_gss
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Solver using Gauss-Seidel method
!> Find f given w (e.g. find grav potential given density mass)
! -----------------------------------------------------------------------------
subroutine solver(N1,min1,max1,x,w,f)
integer, intent(in) :: N1
double precision, intent(in) :: min1, max1
double precision, intent(in) :: x(N1), w(N1)
double precision, intent(inout) :: f(N1)
double precision, parameter :: eps0=1.d-5
integer :: i, it, Nitmax=5000
double precision :: eps, d1, C1, C2, denom
double precision, allocatable :: fth(:), fold(:), c(:), T(:,:)
double precision :: cfl, dt
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Uniform grid
d1=x(2)-x(1)

call get_analytic(N1,min1,max1,x,fth)
call save_vec(N1,0,'fth_',fth)

it=1
eps=1.d99

cfl=1.d0
dt =cfl*(d1**2.d0/2.d0)

call get_operator(N1,T,'no_time_Jacobi',dt,d1)

call get_RHS(N1,w,c,'no_time_Jacobi',dt,d1)

! allocate(c(N1),T(N1,N1))
! T=0.d0
! c(1)=dt*w(1)
! do i=2,N1-1
!   c(i)=dt*w(i)
!   T(i,i-1)=dt/d1**2.d0
!   T(i,i  )=1.d0-2.d0*dt/d1**2.d0
!   T(i,i+1)=dt/d1**2.d0
! enddo
! c(N1)=dt*w(N1)
! ! These 2 ones correspond to periodic BCs
! T(1,N1)=dt/d1**2.d0
! T(N1,1)=dt/d1**2.d0
! T(1,1)  =1.d0-2.d0*dt/d1**2.d0
! T(N1,N1)=1.d0-2.d0*dt/d1**2.d0

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

  call get_Ax(N1,T,fold,f)
  ! do i=1,N1
  !   print*,f(i),w(i),(f(i)-w(i))/w(i)
  ! enddo
  ! stop
  f=f+c

  ! f(1)=-c(1)+T(1,N1)*fold(N1)+T(1,1)*fold(1)+T(1,2)*fold(2)
  ! do i=2,N1-1
  !   f(i)=-c(i)+T(i,i-1)*fold(i-1)+T(i,i)*fold(i)+T(i,i+1)*fold(i+1)
  ! enddo
  ! f(N1)=-c(N1)+T(N1,N1-1)*fold(N1-1)+T(N1,N1)*fold(N1)+T(N1,1)*fold(1)

  ! f(1)=c(1)+T(1,2)*fold(2)+T(1,3)*fold(3)
  ! f(2)=c(2)+T(2,1)*fold(1)+T(2,3)*fold(3)
  ! f(3)=c(3)+T(3,1)*fold(1)+T(3,2)*fold(2)

  ! call bcs(N1,f)

  call save_vec(N1,it,'f_',f)

  ! To prevent from stopping @ 1st iteration if 0 ICs
  if (it>1) call get_residual(N1,fold,f,eps)
  print*, it,eps, f(1), f(N1)

  ! print*, f(1), f(2), f(3)

  it=it+1
  fold=f

enddo

end subroutine solver
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
  fth(i)=-(1.d0/(2.d0*dpi*kwv/(max1-min1))**2.d0)*dsin(2.d0*dpi*kwv*((x(i)-min1)/(max1-min1)))
  ! fth(i)=-1./x(i)
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
!> based on quadratic norm divided by number of points
! -----------------------------------------------------------------------------
subroutine get_residual(N1,fold,f,eps)
integer, intent(in) :: N1
double precision, intent(in) :: fold(N1), f(N1)
double precision, intent(out) :: eps
integer :: i
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

eps=0.d0
do i=1,N1
  ! print*, i, fold(i), f(i)
  if (fold(i)>1.d-8) eps=eps+(dabs(f(i)-fold(i))/max(1.d-14,fold(i)))**2. ! /max(1.d-14,fold(i)**2.)
enddo
eps=dsqrt(eps)/dble(N1)

end subroutine get_residual
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Compute the operator T to be used in the RHS in the solver,
!> to be multiplied w/ x^n and/or x^{n+1} (for Gauss-Seidel).
!> type can be :
!>  - "no_time_Jacobi" to solve DPhi=rho directly w/ Jacobi
!>    (and then, T stands for Laplacian operator)
!>  - "time_Jacobi" to solve dtPhi=DPhi-rho w/ Jacobi
!>  - "GS" to use Gauss-Seidel
!>
!> We assume periodic boundary conditions
!> N.B.: d1 and dt necessary for "time" methods, not for "no time" ones
! -----------------------------------------------------------------------------
subroutine get_operator(N1,T,type,dt,d1)
integer, intent(in) :: N1
double precision, allocatable, intent(out) :: T(:,:)
character(len=*), intent(in) :: type
double precision, intent(in), optional :: dt, d1
integer :: i, j
! only used for GS => allocatable
double precision, allocatable :: DLU(:,:), D(:,:), L(:,:), U(:,:),  &
                                 DL(:,:), LU(:,:), DL_inv(:,:)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

allocate(T(N1,N1))

T=0.d0

if     (type=='time_Jacobi') then

  T(1,1)  =1.d0-2.d0*dt/d1**2.d0
  T(1,2)  =dt/d1**2.d0
  do i=2,N1-1
    T(i,i-1)=dt/d1**2.d0
    T(i,i  )=1.d0-2.d0*dt/d1**2.d0
    T(i,i+1)=dt/d1**2.d0
  enddo
  T(N1,N1)=1.d0-2.d0*dt/d1**2.d0
  T(N1,N1-1)=dt/d1**2.d0
  ! These 2 ones correspond to periodic BCs
  T(1,N1)=dt/d1**2.d0
  T(N1,1)=dt/d1**2.d0

else if (type=='no_time_Jacobi') then

  allocate(DLU(N1,N1),D(N1,N1),L(N1,N1),U(N1,N1),LU(N1,N1))
  call get_laplacian(N1,DLU)
  call dcmps_laplacian(N1,DLU,D,L,U)
  LU=0.d0
  LU=L+U
  T=(-1./D(1,1))*LU
  deallocate(DLU,D,L,U,LU)

else if (type=='no_time_GS') then

  print*, 'A'
  allocate(DL_inv(N1,N1),DL(N1,N1),DLU(N1,N1),D(N1,N1),L(N1,N1),U(N1,N1))
  call get_laplacian(N1,DLU)
  print*, 'B'
  call dcmps_laplacian(N1,DLU,D,L,U)
  DL=0.
  DL=D+L
  print*, 'C'
  call invert_tgl_inf(N1,DL,DL_inv)
  ! Multiply by the triangular superior part of the Laplacian operator
  ! Reuse DL to save memory
  print*, 'D'
  call get_AB(N1,DL_inv,U,DL)
  print*, 'E'
  T=-DL
  deallocate(DL_inv,DL,DLU,D,L,U)

else if (type=='time_GS') then

  ! XXX

else

  call mpistop('not implemented')

endif

end subroutine get_operator
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Compute right hand side term i.e. the source (e.g. 4 pi rho)
!> N.B.: d1 and dt necessary for "time" methods, not for "no time" ones
! -----------------------------------------------------------------------------
subroutine get_RHS(N1,w,c,type,dt,d1)
integer, intent(in) :: N1
double precision, intent(in) :: w(N1)
double precision, allocatable, intent(out) :: c(:)
character(len=*), intent(in) :: type
double precision, intent(in), optional :: dt, d1
integer :: i
! only used for GS => allocatable
double precision, allocatable :: DL_inv(:,:), D(:,:), L(:,:), U(:,:), DLU(:,:), DL(:,:)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

allocate(c(N1))

c=0.d0

if     (type=='time_Jacobi') then

  do i=1,N1
    c(i)=-dt*w(i)
  enddo

else if (type=='no_time_Jacobi') then

  do i=1,N1
    c(i)=-0.5d0*d1**2.d0*w(i)
  enddo

else if (type=='no_time_GS') then

  allocate(DL_inv(N1,N1),DL(N1,N1),DLU(N1,N1),D(N1,N1),L(N1,N1),U(N1,N1))
  call get_laplacian(N1,DLU)
  call dcmps_laplacian(N1,DLU,D,L,U)
  DL=0.
  DL=D+L
  call invert_tgl_inf(N1,DL,DL_inv)
  call get_Ax(N1,DL_inv,d1**2.d0*w,c)
  deallocate(DL_inv,DL,DLU,D,L,U)

else if (type=='time_GS') then



else

  call mpistop('not implemented')

endif

end subroutine get_RHS
! -----------------------------------------------------------------------------


! -----------------------------------------------------------------------------
!> Compute product of square matrix A by vector X and output Y
! -----------------------------------------------------------------------------
subroutine get_Ax(N1,A,X,Y)
integer, intent(in) :: N1
double precision, intent(in)  :: X(N1)
double precision, intent(in) :: A(N1,N1)
double precision, intent(out) :: Y(N1)
integer :: i, j
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Y=0.d0

do i=1,N1
  do j=1,N1
    Y(i)=Y(i)+A(i,j)*X(j)
  enddo
enddo

end subroutine get_Ax
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Compute product of 2 square matrix A and B
! -----------------------------------------------------------------------------
subroutine get_AB(N1,A,B,AB)
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

end subroutine get_AB
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Compute product of 1 square matrix A w/ only diagonal and
!> sub-diagonal elements + 1 element in lower left corner,
!> w/ 1 square matrix B which is diagonal + triangular inferior
! -----------------------------------------------------------------------------
subroutine get_AB_1(N1,A,B,AB)
integer, intent(in) :: N1
double precision, intent(in) :: A(N1,N1), B(N1,N1)
double precision, intent(out) :: AB(N1,N1)
integer :: i, j, k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

AB=0.d0
do i=1,N1
  do j=1,i
    do k=1,i
      AB(i,j)=AB(i,j)+A(i,k)*B(k,j)
    enddo
  enddo
enddo

end subroutine get_AB_1
! ----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Invert triangular inferior matrix
!> Source: https://math.stackexchange.com/questions/1003801/inverse-of-an-invertible-upper-triangular-matrix-of-order-3
!> Applied to specific case where diagonal elements are all the same,
!> but the triangular inferior is not necessarily a diagonal inferior.
! -----------------------------------------------------------------------------
subroutine invert_tgl_inf(N1,DL,DL_inv)
integer, intent(in) :: N1
double precision, intent(in) :: DL(N1,N1)
double precision, intent(out) :: DL_inv(N1,N1)
double precision :: tmp(N1,N1)
double precision :: pdt(N1,N1) ! successive powers of L
double precision :: L(N1,N1) ! triangular inferior part of DL
double precision :: unity(N1,N1) ! unity matrix
double precision :: diag_inv ! invert of diagonal element in DL
integer :: i, j, k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Extract triangular inferior L from DL
! (i.e. get rid of the diagonal elements).
! Also, construct unity matrix.
L=DL
unity=0.d0
do i=1,N1
  L(i,i)=L(i,i)-DL(i,i)
  unity(i,i)=1.d0
enddo

diag_inv=1./DL(1,1)
DL_inv=unity-diag_inv*L ! iteration 0 and 1
pdt=L
do i=2,N1-1 ! stop @ -2 since L^{n}=0 for n>N1-1
  call get_AB(N1,L,pdt,tmp)
  pdt=tmp
  tmp=(-diag_inv)**dble(i)*pdt
  DL_inv=DL_inv+tmp
enddo
DL_inv=DL_inv*diag_inv

! Sanity check
call get_AB(N1,DL,DL_inv,tmp)
if (any(dabs(tmp-unity)>smalldble)) call mpistop("Mx inversion fucked up")

end subroutine invert_tgl_inf
! -----------------------------------------------------------------------------


! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Define 1st order Laplacian operator
!> To do:
!>  - provide possibility to go higher order
!>  - adapt to ghost cells for non periodic cases
! -----------------------------------------------------------------------------
subroutine get_laplacian(N1,DLU)
integer, intent(in) :: N1
double precision, intent(out) :: DLU(N1,N1)
integer :: i
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

DLU(1 ,N1)= 1.d0
DLU(1 ,1 )=-2.d0
DLU(1 ,2 )= 1.d0
do i=2,N1-1
DLU(i ,i-1)= 1.d0
DLU(i ,i  )=-2.d0
DLU(i ,i+1)= 1.d0
enddo
DLU(N1,1   )= 1.d0
DLU(N1,N1  )=-2.d0
DLU(N1,N1-1)= 1.d0

end subroutine get_laplacian
! -----------------------------------------------------------------------------


! -----------------------------------------------------------------------------
!> Decompose Laplacian operator into diagonal, triangular inferior and
!> triangular superior parts
! -----------------------------------------------------------------------------
subroutine dcmps_laplacian(N1,DLU,D,L,U)
integer, intent(in) :: N1
double precision, intent(in) :: DLU(N1,N1)
double precision, intent(out) :: D(N1,N1), L(N1,N1), U(N1,N1)
integer :: i, j
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

D=0.d0
L=0.d0
U=0.d0

do i=1,N1
  do j=1,i-1
    L(i,j)=DLU(i,j)
  enddo
  D(i,i)=DLU(i,i)
  do j=i+1,N1
    U(i,j)=DLU(i,j)
  enddo
enddo

end subroutine dcmps_laplacian
! -----------------------------------------------------------------------------

end module mod_poisson
