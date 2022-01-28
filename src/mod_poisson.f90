module mod_poisson
use mod_io
use mod_params
use mod_csts
use mod_mpi
use mod_linal

implicit none

public :: solver
! public :: get_operator
! public :: get_RHS

contains

! -----------------------------------------------------------------------------
!>
! -----------------------------------------------------------------------------
subroutine unfold(N,N1,N2,f1,f2)
integer, intent(in) :: N, N1, N2
double precision, intent(in) :: f1(N1,N2)
double precision, intent(out) :: f2(N)
integer :: i, j, k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
do k=1,N
  j=int((k-1)/N1)+1
  i=k-(j-1)*N1
  f2(k)=f1(i,j)
enddo

end subroutine unfold
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!>
! -----------------------------------------------------------------------------
subroutine fold(N,N1,N2,f1,f2)
integer, intent(in) :: N, N1, N2
double precision, intent(in) :: f1(N)
double precision, intent(out) :: f2(N1,N2)
integer :: i, j, k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
do i=1,N1
do j=1,N2
  k=i+(j-1)*N1
  f2(i,j)=f1(k)
enddo
enddo

end subroutine fold
! -----------------------------------------------------------------------------


! -----------------------------------------------------------------------------
!> BiConjugate Gradient STABilized solver
!>
!> @warning It is absolutely necessary to apply the BCs each time
!> @warning the Laplacian is applied (i.e. here immediately after s and p1
!> @warning are computed)
!>
!> @see https://www.codeproject.com/Articles/771337/A-Parallel-Class-Study-Bi-Conjugate-Gradient-Stabi
!> @see https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method#:~:text=In%20numerical%20linear%20algebra%2C%20the,solution%20of%20nonsymmetric%20linear%20systems
!> @see https://utminers.utep.edu/xzeng/2017spring_math5330/MATH_5330_Computational_Methods_of_Linear_Algebra_files/ln07.pdf
!> @see PDF document "NUMERICAL EFFICIENCY OF ITERATIVE SOLVERS FOR
!> @see THE POISSON EQUATION USING COMPUTER CLUSTER" by J. Gościk & J. Gościk
! -----------------------------------------------------------------------------
subroutine solve_BiCGSTAB(N,N1,N2,NGC,iImin1,iImax1,&
  iOmin1,iOmax1,iImin2,iImax2,iOmin2,iOmax2,min1,max1,min2,max2,pencil,eps0,Nitmax,bc_type,grid_type,x,w,f1)
use mod_grid
use mod_bc
use mod_laplacian
integer, intent(in) :: N, N1, N2, NGC
integer, intent(in) :: iImin1, iImax1, iOmin1, iOmax1
integer, intent(in) :: iImin2, iImax2, iOmin2, iOmax2
double precision, intent(in) :: min1, max1
double precision, intent(in) :: min2, max2
integer, intent(in) :: pencil
double precision, intent(in) :: eps0
integer, intent(in) :: Nitmax
character(len=*), intent(in) :: bc_type, grid_type
double precision, intent(in) :: x(N1,N2,2), w(N)
double precision, intent(out) :: f1(N)
double precision :: f1_fold(N1,N2)
double precision :: A(N,N)
double precision :: b(N)
double precision :: f0(N), r0(N), p0(N)
double precision :: Ap(N)
double precision :: r1(N), p1(N)
double precision :: num, denum
double precision :: tmp, eps
double precision :: alpha, omega
double precision :: beta
double precision :: s(N)
double precision :: t(N)
double precision :: r0_star(N)
double precision :: f0_fold(N1,N2)
double precision :: r0_fold(N1,N2)
integer :: i, j, it
double precision :: f_avg
double precision :: dx(N1,N2,2)
double precision :: dpi=dacos(-1.d0)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

call get_laplacian(N,N1,N2,pencil,x,A)

call get_dx(N1,N2,NGC,x,min1,max1,min2,max2,grid_type,dx)

do i=1,N
  b(i)=w(i)
enddo

! 1st guess MUST have the right BCs
! (and 0 average if periodic BCs)
f0=0.d0
! do i=1,N1
! do j=1,N2
!   f0_fold(i,j)=x(i,j,1)**4.d0/12.d0+x(i,j,2)**4.d0/12.d0
! enddo
! enddo
! call rm_avg(N1,N2,NGC,min1,max1,min2,max2,x,dx,f0_fold)
! call unfold(N,N1,N2,f0_fold,f0)

call get_bc(N,N1,N2,NGC,bc_type,x,f0)
call fold(N,N1,N2,f0,f1_fold)
call get_avg(N1,N2,NGC,min1,max1,min2,max2,x,dx,f1_fold,f_avg)
call save_vec(N1,N2,0,'f_',f1_fold-f_avg)

call mx_x_vec(N,N1,N2,NGC,A,f0,Ap)
r0=b-Ap
r0_star=r0
p0  =r0
call get_bc_BiCGSTAB_p(N,N1,N2,NGC,bc_type,p0)
! call get_bc(N,N1,N2,NGC,bc_type,p0)

call fold(N,N1,N2,r0,r0_fold)
call save_vec(N1,N2,0,'r_',r0_fold)

it=1
eps=1.d99
do while (eps>eps0 .and. it<Nitmax)

  call fold(N,N1,N2,r0,r0_fold)
  call save_vec(N1,N2,it,'r_',r0_fold)

  call dot_pdct(N,N1,N2,NGC,r0,r0_star,num)
  call mx_x_vec(N,N1,N2,NGC,A,p0,Ap)
  call dot_pdct(N,N1,N2,NGC,Ap,r0_star,denum)
  alpha=num/denum
  ! print*, alpha
  ! print*, 1.d-2/(4.d0*(dcos(2.d0*dpi*1.d-1)-1.d0))
  ! print*, -1.d0/(2.d0*(2.d0*dpi)**2.d0)
  ! stop

  ! Intermediate check
  ! This step is necessary to prevent having s~0 (which will happen if
  ! it=1, r0=b, p0=b and alpha*Ap=alpha*A*b=A*f) which would lead to
  ! omega~0/0=WTF
  f1=f0+alpha*p0
  call get_residual(N,N1,N2,NGC,f1,b,A,eps)
  if (eps<eps0) then
    call get_bc(N,N1,N2,NGC,bc_type,x,f1)
    call fold(N,N1,N2,f1,f1_fold)
    call rm_avg(N1,N2,NGC,min1,max1,min2,max2,x,dx,f1_fold)
    call save_vec(N1,N2,it,'f_',f1_fold)
    print*, it, eps ! , 'kikou'
    exit
  endif

  s=r0-alpha*Ap
  call get_bc_BiCGSTAB_s(N,N1,N2,NGC,bc_type,s)
  ! call get_bc_BiCGSTAB_s(N1,NGC,iOmin1,iOmax1,bc_type,s)
  ! call get_bc(N,N1,N2,NGC,bc_type,s)

  call mx_x_vec(N,N1,N2,NGC,A,s,t)
  call dot_pdct(N,N1,N2,NGC,t,t,denum)
  call dot_pdct(N,N1,N2,NGC,t,s,num)
  omega=num/denum
  f1=f0+alpha*p0+omega*s
  r1=s-omega*t

  call dot_pdct(N,N1,N2,NGC,r0,r0_star,denum)
  call dot_pdct(N,N1,N2,NGC,r1,r0_star,num)
  beta=(num/denum)*(alpha/omega)
  p1=r1+beta*(p0-omega*Ap)
  call get_bc_BiCGSTAB_p(N,N1,N2,NGC,bc_type,p1)
  ! call get_bc_BiCGSTAB_p(N1,NGC,iOmin1,iOmax1,bc_type,p1)
  ! call get_bc(N,N1,N2,NGC,bc_type,p1)

  ! call rm_avg(N1,NGC,min1,max1,x,dx,f1)
  call get_bc(N,N1,N2,NGC,bc_type,x,f1)

  f0=f1
  r0=r1
  p0=p1

  call fold(N,N1,N2,f1,f1_fold)
  call rm_avg(N1,N2,NGC,min1,max1,min2,max2,x,dx,f1_fold)
  ! print*, 'avg:', f_avg
  call save_vec(N1,N2,it,'f_',f1_fold) ! -sum(f1)/dble(N1))

  call get_residual(N,N1,N2,NGC,f1,b,A,eps)
  print*, it, eps

  it=it+1

enddo

end subroutine solve_BiCGSTAB
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Compute the average of scalar field f,
!> weighted on grid x w/ spacing dx of each cell
!>
!> @warning Works only over the simulation space, not in the ghost cells
! -----------------------------------------------------------------------------
subroutine get_avg(N1,N2,NGC,min1,max1,min2,max2,x,dx,f,f_avg)
integer, intent(in) :: N1
integer, intent(in) :: N2
integer, intent(in) :: NGC
double precision, intent(in) :: min1, max1
double precision, intent(in) :: min2, max2
double precision, intent(in) :: x(N1,N2,2)
double precision, intent(in) :: dx(N1,N2,2)
double precision, intent(in) :: f(N1,N2)
double precision, optional, intent(out) :: f_avg
integer :: i, j
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
f_avg=0.d0
do i=NGC+1,N1-NGC
do j=NGC+1,N2-NGC
  f_avg=f_avg+f(i,j)*dx(i,j,1)*dx(i,j,2)
enddo
enddo
f_avg=f_avg/((max1-min1)*(max2-min2))
! print*, f_avg
end subroutine get_avg
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Removes average value from scalar field f in simulation space,
!> weighted on grid x w/ spacing dx of each cell
!>
!> @note Might be necessary to subtract <w> before (for periodic BCs only)
! -----------------------------------------------------------------------------
subroutine rm_avg(N1,N2,NGC,min1,max1,min2,max2,x,dx,f)
integer, intent(in) :: N1
integer, intent(in) :: N2
integer, intent(in) :: NGC
double precision, intent(in) :: min1, max1
double precision, intent(in) :: min2, max2
double precision, intent(in) :: x(N1,N2,2)
double precision, intent(in) :: dx(N1,N2,2)
double precision, intent(inout) :: f(N1,N2)
double precision :: f_avg
integer :: i
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
call get_avg(N1,N2,NGC,min1,max1,min2,max2,x,dx,f,f_avg)
f=f-f_avg
end subroutine rm_avg
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Solver using Gauss-Seidel method
!> Find f given w (e.g. find grav potential given density mass)
! -----------------------------------------------------------------------------
subroutine solver(N1,N2,NGC,iOmin1,iOmax1,iOmin2,iOmax2,min1,max1,min2,max2,pencil,solver_type,bc_type,grid_type,x,w)
use mod_bc
use mod_grid
use mod_laplacian
integer, intent(in) :: N1, N2, NGC
integer, intent(in) :: iOmin1, iOmax1
integer, intent(in) :: iOmin2, iOmax2
double precision, intent(in) :: min1, max1
double precision, intent(in) :: min2, max2
integer, intent(in) :: pencil
character(len=*), intent(in) :: solver_type, bc_type, grid_type
double precision, intent(in) :: x(N1,N2,2)
double precision, intent(inout) :: w(N1,N2)
double precision :: w_unfold(N1*N2)
double precision :: fth_unfold(N1*N2)
double precision :: f(N1*N2)
double precision :: f_fold(N1,N2)
double precision, parameter :: eps0=1.d-6
integer :: Nitmax=400
double precision, allocatable :: fth(:,:)
double precision :: dx(N1,N2,2)
integer :: N, i
double precision :: num
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

call get_dx(N1,N2,NGC,x,min1,max1,min2,max2,grid_type,dx)

call get_analytic(N1,N2,NGC,min1,max1,min2,max2,x,dx,fth)
call save_vec(N1,N2,0,'fth_',fth)

N=N1*N2

call unfold(N,N1,N2,w,w_unfold)
! open(1,file='test')
! do i=1,N
!   write(1,*) w_unfold(i)
! enddo
! close(1)

call solve_BiCGSTAB(N,N1,N2,NGC,iImin1,iImax1,&
iOmin1,iOmax1,iImin2,iImax2,iOmin2,iOmax2,min1,max1,min2,max2,pencil,eps0,Nitmax,bc_type,grid_type,x,w_unfold,f)

call unfold(N,N1,N2,fth,fth_unfold)
call set_GC_to_x0(N,N1,N2,NGC,0.d0,fth_unfold)

call fold(N,N1,N2,f,f_fold)
call rm_avg(N1,N2,NGC,min1,max1,min2,max2,x,dx,f_fold)
call unfold(N,N1,N2,f_fold,f)
call set_GC_to_x0(N,N1,N2,NGC,0.d0,f)

print*, 'Max gap:', maxval(dabs(fth_unfold-f))

end subroutine solver
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Specify analytic solution to be compared to by hand
! -----------------------------------------------------------------------------
subroutine get_analytic(N1,N2,NGC,min1,max1,min2,max2,x,dx,fth)
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
integer, intent(in) :: N1
integer, intent(in) :: N2
integer, intent(in) :: NGC
double precision, intent(in) :: min1, max1
double precision, intent(in) :: min2, max2
double precision, intent(in) :: x(N1,N2,2)
double precision, intent(in) :: dx(N1,N2,2)
double precision, allocatable, intent(out) :: fth(:,:)
integer :: i, j
double precision :: dpi=dacos(-1.d0)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

analytic => analytic_

allocate(fth(N1,N2))

do i=1,N1
do j=1,N2
  ! fth(i,j)=-(1.d0/(2.d0*(2.d0*dpi*kwv/(max1-min1))**2.d0))*&
  !          dsin(2.d0*dpi*kwv*((x(i,j,1)-min1)/(max1-min1)))*&
  !          dsin(2.d0*dpi*kwv*((x(i,j,2)-min2)/(max2-min2)))
  ! fth(i,j)=-(1.d0/(2.d0*dpi*kwv/(max1-min1))**2.d0)*&
  !          dexp(dsin(2.d0*dpi*kwv*((x(i,j,1)-min1)/(max1-min1))))*&
  !          dexp(dcos(2.d0*dpi*kwv*((x(i,j,2)-min2)/(max2-min2))))
  ! fth(i)=-(1.d0/(2.d0*dpi*kwv/(max1-min1))**2.d0)*dexp(dsin(2.d0*dpi*kwv*((x(i)-min1)/(max1-min1))))
  ! fth(i)=(1.d0/(2.d0*dpi*kwv/(max1-min1))**2.d0)/(2.d0+dsin(2.d0*dpi*kwv*((x(i)-min1)/(max1-min1))))
  fth(i,j)=analytic(x(i,j,1),x(i,j,2)) ! **4.d0/12.d0+x(i,j,2)**4.d0/12.d0
  ! fth(i)=dexp(-((x(i)-1.5d0)/0.4d0)**2.d0)
  ! fth(i)=1.d0/x(i)+1.d0
enddo
enddo
! print*, fth(1,:)
! print*, ' '
! print*, fth(N1,:)
! print*, ' '
! print*, fth(:,1)
! print*, ' '
! print*, fth(:,N2)
! print*, ' '
! stop
call rm_avg(N1,N2,NGC,min1,max1,min2,max2,x,dx,fth)
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
!> Compute residual as norm infinite |·|_{inf}=maxval(|·|)
! -----------------------------------------------------------------------------
subroutine get_residual(N,N1,N2,NGC,f,b,A,eps)
use mod_laplacian
use mod_linal
use mod_bc
integer, intent(in) :: N, N1, N2, NGC
double precision, intent(in) :: f(N)
double precision, intent(in) :: b(N)
double precision, intent(in) :: A(N,N)
double precision, intent(out) :: eps
double precision :: Af(N)
double precision :: res(N), tmp
double precision :: dpi=dacos(-1.d0)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! TO CORRECT XXX

eps=0.d0

call mx_x_vec(N,N1,N2,NGC,A,f,Af)
Af=dabs(Af-b)
! Trick to not pick up ghost cells
call set_GC_to_x0(N,N1,N2,NGC,0.d0,Af)

eps=maxval(Af)
! print*, 'eps', eps

! print*, (f(1)-2.d0*f(2)+f(3))/(1.d0/dble(100))**2.d0
! print*, A(2,1)*f(1)+A(2,2)*f(2)+A(2,3)*f(3)
! print*, eps
! stop

end subroutine get_residual
! -----------------------------------------------------------------------------

! ! -----------------------------------------------------------------------------
! !> Compute residual as relative difference w.r.t. solution @ previous step,
! !> based on quadratic norm divided by number of points
! ! -----------------------------------------------------------------------------
! subroutine get_residual_w_prv(N1,fold,f,eps)
! integer, intent(in) :: N1
! double precision, intent(in) :: fold(N1), f(N1)
! double precision, intent(out) :: eps
! integer :: i
! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! eps=0.d0
! do i=1,N1
!   ! print*, i, fold(i), f(i)
!   if (fold(i)>1.d-8) eps=eps+(dabs(f(i)-fold(i))/max(1.d-14,fold(i)))**2. ! /max(1.d-14,fold(i)**2.)
! enddo
! eps=dsqrt(eps)/dble(N1)
!
! end subroutine get_residual_w_prv
! ! -----------------------------------------------------------------------------
!
! ! -----------------------------------------------------------------------------
! !> Compute the operator T to be used in the RHS in the solver,
! !> to be multiplied w/ x^n and/or x^{n+1} (for Gauss-Seidel).
! !> type can be :
! !>  - "no_time_Jacobi" to solve DPhi=rho directly w/ Jacobi
! !>    (and then, T stands for Laplacian operator)
! !>  - "time_Jacobi" to solve dtPhi=DPhi-rho w/ Jacobi
! !>  - "GS" to use Gauss-Seidel
! !>
! !> We assume periodic boundary conditions
! !> N.B.: d1 and dt necessary for "time" methods, not for "no time" ones
! ! -----------------------------------------------------------------------------
! subroutine get_operator(N1,pencil,T,solver_type,dt,d1,x)
! use mod_laplacian
! integer, intent(in) :: N1
! integer, intent(in) :: pencil
! double precision, allocatable, intent(out) :: T(:,:)
! character(len=*), intent(in) :: solver_type
! double precision, intent(in), optional :: dt, d1
! double precision, intent(in) :: x(N1)
! integer :: i, j
! ! only used for GS => allocatable
! double precision, allocatable :: DLU(:,:), D(:,:), L(:,:), U(:,:),  &
!                                  DL(:,:), LU(:,:), DL_inv(:,:)
! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! allocate(T(N1,N1))
!
! T=0.d0
!
! if     (solver_type=='time_Jacobi') then
!
!   T(1,1)  =1.d0-2.d0*dt/d1**2.d0
!   T(1,2)  =dt/d1**2.d0
!   do i=2,N1-1
!     T(i,i-1)=dt/d1**2.d0
!     T(i,i  )=1.d0-2.d0*dt/d1**2.d0
!     T(i,i+1)=dt/d1**2.d0
!   enddo
!   T(N1,N1)=1.d0-2.d0*dt/d1**2.d0
!   T(N1,N1-1)=dt/d1**2.d0
!   ! These 2 ones correspond to periodic BCs
!   T(1,N1)=dt/d1**2.d0
!   T(N1,1)=dt/d1**2.d0
!
! else if (solver_type=='no_time_Jacobi') then
!
!   allocate(DLU(N1,N1),D(N1,N1),L(N1,N1),U(N1,N1),LU(N1,N1))
!   call get_laplacian(N1,pencil,x,DLU)
!   call dcmps_laplacian(N1,DLU,D,L,U)
!   LU=0.d0
!   LU=L+U
!   T=(-1./D(1,1))*LU
!   deallocate(DLU,D,L,U,LU)
!
! else if (solver_type=='no_time_GS') then
!
!   print*, 'A'
!   allocate(DL_inv(N1,N1),DL(N1,N1),DLU(N1,N1),D(N1,N1),L(N1,N1),U(N1,N1))
!   call get_laplacian(N1,pencil,x,DLU)
!   print*, 'B'
!   call dcmps_laplacian(N1,DLU,D,L,U)
!   DL=0.
!   DL=D+L
!   print*, 'C'
!   call invert_tgl_inf(N1,DL,DL_inv)
!   ! Multiply by the triangular superior part of the Laplacian operator
!   ! Reuse DL to save memory
!   print*, 'D'
!   call get_AB(N1,DL_inv,U,DL)
!   print*, 'E'
!   T=-DL
!   deallocate(DL_inv,DL,DLU,D,L,U)
!
! else if (solver_type=='time_GS') then
!
!   ! XXX
!
! else
!
!   call mpistop('not implemented')
!
! endif
!
! end subroutine get_operator
! ! -----------------------------------------------------------------------------
!
! ! -----------------------------------------------------------------------------
! !> Compute right hand side term i.e. the source (e.g. 4 pi rho)
! !> N.B.: d1 and dt necessary for "time" methods, not for "no time" ones
! ! -----------------------------------------------------------------------------
! subroutine get_RHS(N1,pencil,w,c,solver_type,dt,d1,x)
! use mod_laplacian
! integer, intent(in) :: N1
! integer, intent(in) :: pencil
! double precision, intent(in) :: w(N1)
! double precision, allocatable, intent(out) :: c(:)
! character(len=*), intent(in) :: solver_type
! double precision, intent(in), optional :: dt, d1
! double precision, intent(in) :: x(N1)
! integer :: i
! ! only used for GS => allocatable
! double precision, allocatable :: DL_inv(:,:), D(:,:), L(:,:), U(:,:), DLU(:,:), DL(:,:)
! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! allocate(c(N1))
!
! c=0.d0
!
! if     (solver_type=='time_Jacobi') then
!
!   do i=1,N1
!     c(i)=-dt*w(i)
!   enddo
!
! else if (solver_type=='no_time_Jacobi') then
!
!   do i=1,N1
!     c(i)=-0.5d0*d1**2.d0*w(i)
!   enddo
!
! else if (solver_type=='no_time_GS') then
!
!   allocate(DL_inv(N1,N1),DL(N1,N1),DLU(N1,N1),D(N1,N1),L(N1,N1),U(N1,N1))
!   call get_laplacian(N1,pencil,x,DLU)
!   call dcmps_laplacian(N1,DLU,D,L,U)
!   DL=0.
!   DL=D+L
!   call invert_tgl_inf(N1,DL,DL_inv)
!   call mx_x_vec(N1,iOmin1,iOmax1,DL_inv,d1**2.d0*w,c)
!   deallocate(DL_inv,DL,DLU,D,L,U)
!
! else if (solver_type=='time_GS') then
!
!
!
! else
!
!   call mpistop('not implemented')
!
! endif
!
! end subroutine get_RHS
! ! -----------------------------------------------------------------------------
!
!
! ! ! -----------------------------------------------------------------------------
! ! !> Compute product of square matrix A by vector X and output Y
! ! ! -----------------------------------------------------------------------------
! ! subroutine get_Ax(N1,A,X,Y)
! ! integer, intent(in) :: N1
! ! double precision, intent(in)  :: X(N1)
! ! double precision, intent(in) :: A(N1,N1)
! ! double precision, intent(out) :: Y(N1)
! ! integer :: i, j
! ! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! !
! ! Y=0.d0
! !
! ! do i=1,N1
! !   do j=1,N1
! !     Y(i)=Y(i)+A(i,j)*X(j)
! !   enddo
! ! enddo
! !
! ! end subroutine get_Ax
! ! ! -----------------------------------------------------------------------------
!
! ! -----------------------------------------------------------------------------
! !> Compute product of 2 square matrix A and B
! ! -----------------------------------------------------------------------------
! subroutine get_AB(N1,A,B,AB)
! integer, intent(in) :: N1
! double precision, intent(in) :: A(N1,N1), B(N1,N1)
! double precision, intent(out) :: AB(N1,N1)
! integer :: i, j, k
! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! AB=0.d0
! do i=1,N1
!   do j=1,N1
!     do k=1,N1
!       AB(i,j)=AB(i,j)+A(i,k)*B(k,j)
!     enddo
!   enddo
! enddo
!
! end subroutine get_AB
! ! -----------------------------------------------------------------------------
!
! ! -----------------------------------------------------------------------------
! !> Compute product of 1 square matrix A w/ only diagonal and
! !> sub-diagonal elements + 1 element in lower left corner,
! !> w/ 1 square matrix B which is diagonal + triangular inferior
! ! -----------------------------------------------------------------------------
! subroutine get_AB_1(N1,A,B,AB)
! integer, intent(in) :: N1
! double precision, intent(in) :: A(N1,N1), B(N1,N1)
! double precision, intent(out) :: AB(N1,N1)
! integer :: i, j, k
! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! AB=0.d0
! do i=1,N1
!   do j=1,i
!     do k=1,i
!       AB(i,j)=AB(i,j)+A(i,k)*B(k,j)
!     enddo
!   enddo
! enddo
!
! end subroutine get_AB_1
! ! ----------------------------------------------------------------------------
!
! ! -----------------------------------------------------------------------------
! !> Invert triangular inferior matrix
! !> Source: https://math.stackexchange.com/questions/1003801/inverse-of-an-invertible-upper-triangular-matrix-of-order-3
! !> Applied to specific case where diagonal elements are all the same,
! !> but the triangular inferior is not necessarily a diagonal inferior.
! ! -----------------------------------------------------------------------------
! subroutine invert_tgl_inf(N1,DL,DL_inv)
! integer, intent(in) :: N1
! double precision, intent(in) :: DL(N1,N1)
! double precision, intent(out) :: DL_inv(N1,N1)
! double precision :: tmp(N1,N1)
! double precision :: pdt(N1,N1) ! successive powers of L
! double precision :: L(N1,N1) ! triangular inferior part of DL
! double precision :: unity(N1,N1) ! unity matrix
! double precision :: diag_inv ! invert of diagonal element in DL
! integer :: i, j, k
! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! ! Extract triangular inferior L from DL
! ! (i.e. get rid of the diagonal elements).
! ! Also, construct unity matrix.
! L=DL
! unity=0.d0
! do i=1,N1
!   L(i,i)=L(i,i)-DL(i,i)
!   unity(i,i)=1.d0
! enddo
!
! diag_inv=1./DL(1,1)
! DL_inv=unity-diag_inv*L ! iteration 0 and 1
! pdt=L
! do i=2,N1-1 ! stop @ -2 since L^{n}=0 for n>N1-1
!   call get_AB(N1,L,pdt,tmp)
!   pdt=tmp
!   tmp=(-diag_inv)**dble(i)*pdt
!   DL_inv=DL_inv+tmp
! enddo
! DL_inv=DL_inv*diag_inv
!
! ! Sanity check
! call get_AB(N1,DL,DL_inv,tmp)
! if (any(dabs(tmp-unity)>smalldble)) call mpistop("Mx inversion fucked up")
!
! end subroutine invert_tgl_inf
! ! -----------------------------------------------------------------------------
!
!
! ! -----------------------------------------------------------------------------
!
! ! ! -----------------------------------------------------------------------------
! ! !> Define 1st order Laplacian operator
! ! !>
! ! !> @param pencil is 3 (2 closest neighbours), 5 (4 closest neighbours), etc
! ! !> @see https://en.wikipedia.org/wiki/Finite_difference_coefficient
! ! !> @todo - validate higher order than pencil=3
! ! !> @todo - adapt to ghost cells for non periodic cases
! ! ! -----------------------------------------------------------------------------
! ! subroutine get_laplacian(N1,pencil,x,DLU)
! ! integer, intent(in) :: N1
! ! integer, intent(in) :: pencil
! ! double precision, intent(in) :: x(N1)
! ! double precision, intent(out) :: DLU(N1,N1)
! ! integer :: i, NGC
! ! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! !
! ! NGC=(pencil-1)/2
! !
! ! DLU=0.d0
! !
! ! if      (pencil==3) then
! !   ! DLU(1 ,N1)= 1.d0
! !   DLU(1 ,1 )=0.d0 ! -2.d0
! !   DLU(1 ,2 )=0.d0 ! 1.d0
! !   do i=NGC+1,N1-NGC
! !   DLU(i ,i-1)= (x(i+1)-x(i  ))/((x(i)-x(i-1))*(x(i+1)-x(i))*0.5d0*(x(i+1)-x(i-1)))
! !   DLU(i ,i  )=-(x(i+1)-x(i-1))/((x(i)-x(i-1))*(x(i+1)-x(i))*0.5d0*(x(i+1)-x(i-1)))
! !   DLU(i ,i+1)= (x(i  )-x(i-1))/((x(i)-x(i-1))*(x(i+1)-x(i))*0.5d0*(x(i+1)-x(i-1)))
! !   ! print*, DLU(i ,i-1), DLU(i ,i  ), DLU(i ,i+1)
! !   ! DLU(i ,i-1)= 1.d0
! !   ! DLU(i ,i  )=-2.d0
! !   ! DLU(i ,i+1)= 1.d0
! !   enddo
! !   ! DLU(N1,1   )= 1.d0
! !   DLU(N1,N1-1)=0.d0 ! 1.d0
! !   DLU(N1,N1  )=0.d0 ! -2.d0
! ! else if (pencil==5) then
! !   DLU(1 ,1 )=-5.d0/2.d0
! !   DLU(1 ,2 )= 4.d0/3.d0
! !   DLU(1 ,3 )=-1.d0/12.d0
! !   DLU(2 ,1 )= 4.d0/3.d0
! !   DLU(2 ,2 )=-5.d0/2.d0
! !   DLU(2 ,3 )= 4.d0/3.d0
! !   DLU(2 ,4 )=-1.d0/12.d0
! !   do i=NGC+1,N1-NGC
! !   DLU(i ,i-2)=-1.d0/12.d0
! !   DLU(i ,i-1)= 4.d0/3.d0
! !   DLU(i ,i  )=-5.d0/2.d0
! !   DLU(i ,i+1)= 4.d0/3.d0
! !   DLU(i ,i+2)=-1.d0/12.d0
! !   enddo
! !   DLU(N1-1,N1-3)=-1.d0/12.d0
! !   DLU(N1-1,N1-2)= 4.d0/3.d0
! !   DLU(N1-1,N1-1)=-5.d0/2.d0
! !   DLU(N1-1,N1  )= 4.d0/3.d0
! !   DLU(N1,N1-2)=-1.d0/12.d0
! !   DLU(N1,N1-1)= 4.d0/3.d0
! !   DLU(N1,N1  )=-5.d0/2.d0
! ! endif
! ! ! DLU=DLU/(x(2)-x(1))**2.d0
! !
! ! end subroutine get_laplacian
! ! ! -----------------------------------------------------------------------------
!
! ! ! -----------------------------------------------------------------------------
! ! !> Define p-th order Laplacian operator w/ pencil made of (p-1)/2
! ! !> neighbors on each side + itself
! ! !>
! ! !> @param pencil is 3 (2 closest neighbours), 5 (4 closest neighbours), etc
! ! !> @see https://en.wikipedia.org/wiki/Finite_difference_coefficient
! ! !> @see http://www.iaeng.org/publication/WCE2017/WCE2017_pp127-129.pdf
! ! ! -----------------------------------------------------------------------------
! ! subroutine get_laplacian(N1,pencil,x,DLU)
! ! integer, intent(in) :: N1
! ! integer, intent(in) :: pencil
! ! double precision, intent(in) :: x(N1)
! ! double precision, intent(out) :: DLU(N1,N1)
! ! double precision :: lbd(pencil)
! ! double precision :: W(pencil,pencil)
! ! double precision :: A(pencil,pencil)
! ! double precision :: V_inv(pencil,pencil)
! ! double precision :: out(pencil)
! ! double precision :: coeff(pencil)
! ! double precision :: tmp
! ! integer :: i, j, NGC
! ! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! !
! ! NGC=(pencil-1)/2
! !
! ! DLU=0.d0
! !
! ! if      (pencil==3) then
! !   ! DLU(1 ,1 )=-2.d4
! !   ! DLU(1 ,2 )= 1.d4
! !   do i=NGC+1,N1-NGC
! !
! !     lbd(1)=x(i-1)-x(i)
! !     lbd(2)=0.d0
! !     lbd(3)=x(i+1)-x(i)
! !     A(1,1)=1.d0
! !     A(1,2)=0.d0
! !     A(1,3)=0.d0
! !     A(2,1)=-( lbd(1) + lbd(2) + lbd(3) )
! !     A(2,2)=1.d0
! !     A(2,3)=0.d0
! !     A(3,1)= ( lbd(1)*lbd(2) + lbd(1)*lbd(3) + &
! !               lbd(2)*lbd(3) )
! !     A(3,2)=-( lbd(1) + lbd(2) + lbd(3) )
! !     A(3,3)=1.d0
! !     tmp = (lbd(1) - lbd(2)) * (lbd(1) - lbd(3))
! !     W(1,1)=lbd(1)**2.d0/tmp
! !     W(1,2)=lbd(1)**1.d0/tmp
! !     W(1,3)=lbd(1)**0.d0/tmp
! !     tmp = (lbd(2) - lbd(1)) * (lbd(2) - lbd(3))
! !     W(2,1)=lbd(2)**2.d0/tmp
! !     W(2,2)=lbd(2)**1.d0/tmp
! !     W(2,3)=lbd(2)**0.d0/tmp
! !     tmp = (lbd(3) - lbd(1)) * (lbd(3) - lbd(2))
! !     W(3,1)=lbd(3)**2.d0/tmp
! !     W(3,2)=lbd(3)**1.d0/tmp
! !     W(3,3)=lbd(3)**0.d0/tmp
! !     call mx_x_mx(3,W,A,V_inv)
! !     out(1)=0.d0
! !     out(2)=0.d0
! !     out(3)=2.d0
! !     call mx_x_vec(3,1,3,V_inv,out,coeff)
! !     ! print*, coeff/1.d4
! !     ! stop
! !     ! do i=1,3
! !     !   print*, coeff(i)*9.d0
! !     !   ! do j=1,3
! !     !   !   print*, W(i,j), A(i,j), V_inv(i,j)
! !     !   ! enddo
! !     ! enddo
! !     ! stop
! !     DLU(i ,i-1)= coeff(1) ! (x(i+1)-x(i  ))/((x(i)-x(i-1))*(x(i+1)-x(i))*0.5d0*(x(i+1)-x(i-1)))
! !     DLU(i ,i  )= coeff(2) ! -(x(i+1)-x(i-1))/((x(i)-x(i-1))*(x(i+1)-x(i))*0.5d0*(x(i+1)-x(i-1)))
! !     DLU(i ,i+1)= coeff(3) ! (x(i  )-x(i-1))/((x(i)-x(i-1))*(x(i+1)-x(i))*0.5d0*(x(i+1)-x(i-1)))
! !
! !     ! DLU(i ,i-1)= (x(i+1)-x(i  ))/((x(i)-x(i-1))*(x(i+1)-x(i))*0.5d0*(x(i+1)-x(i-1)))
! !     ! DLU(i ,i  )=-(x(i+1)-x(i-1))/((x(i)-x(i-1))*(x(i+1)-x(i))*0.5d0*(x(i+1)-x(i-1)))
! !     ! DLU(i ,i+1)= (x(i  )-x(i-1))/((x(i)-x(i-1))*(x(i+1)-x(i))*0.5d0*(x(i+1)-x(i-1)))
! !
! !     ! DLU(i ,i-1)= 1.d4
! !     ! DLU(i ,i  )=-2.d4
! !     ! DLU(i ,i+1)= 1.d4
! !
! !     ! print*, DLU(i ,i-1)*1.d-4, DLU(i ,i)*1.d-4, DLU(i ,i+1)*1.d-4, x(i+1)-x(i)
! !
! !   enddo
! !   ! DLU(N1,N1-1)= 1.d-4
! !   ! DLU(N1,N1  )=-2.d-4
! !
! ! else if (pencil==5) then
! !
! !   do i=NGC+1,N1-NGC
! !     lbd(1)=x(i-2)-x(i)
! !     lbd(2)=x(i-1)-x(i)
! !     lbd(3)=0.d0
! !     lbd(4)=x(i+1)-x(i)
! !     lbd(5)=x(i+2)-x(i)
! !     A(1,1)=1.d0
! !     A(1,2)=0.d0
! !     A(1,3)=0.d0
! !     A(1,4)=0.d0
! !     A(1,5)=0.d0
! !     A(2,1)=-( lbd(1) + lbd(2) + lbd(3) + lbd(4) + lbd(5) )
! !     A(2,2)=1.d0
! !     A(2,3)=0.d0
! !     A(2,4)=0.d0
! !     A(2,5)=0.d0
! !     A(3,1)= ( lbd(1)*lbd(2) + lbd(1)*lbd(3) + lbd(1)*lbd(4) + lbd(1)*lbd(5) + &
! !               lbd(2)*lbd(3) + lbd(2)*lbd(4) + lbd(2)*lbd(5) + &
! !               lbd(3)*lbd(4) + lbd(3)*lbd(5) + &
! !               lbd(4)*lbd(5) )
! !     A(3,2)=-( lbd(1) + lbd(2) + lbd(3) + lbd(4) + lbd(5) )
! !     A(3,3)=1.d0
! !     A(3,4)=0.d0
! !     A(3,5)=0.d0
! !     A(4,1)=-( lbd(1)*lbd(2)*lbd(3) + lbd(1)*lbd(2)*lbd(4) + lbd(1)*lbd(2)*lbd(5) + &
! !               lbd(1)*lbd(3)*lbd(4) + lbd(1)*lbd(3)*lbd(5) + &
! !               lbd(2)*lbd(3)*lbd(4) + lbd(2)*lbd(3)*lbd(5) + &
! !               lbd(2)*lbd(4)*lbd(5) + &
! !               lbd(3)*lbd(4)*lbd(5) )
! !     A(4,2)= ( lbd(1)*lbd(2) + lbd(1)*lbd(3) + lbd(1)*lbd(4) + lbd(1)*lbd(5) + &
! !               lbd(2)*lbd(3) + lbd(2)*lbd(4) + lbd(2)*lbd(5) + &
! !               lbd(3)*lbd(4) + lbd(3)*lbd(5) + &
! !               lbd(4)*lbd(5) )
! !     A(4,3)=-( lbd(1) + lbd(2) + lbd(3) + lbd(4) + lbd(5) )
! !     A(4,4)=1.d0
! !     A(4,5)=0.d0
! !     A(5,1)= ( lbd(1)*lbd(2)*lbd(3)*lbd(4) + lbd(1)*lbd(2)*lbd(3)*lbd(5) + &
! !               lbd(2)*lbd(3)*lbd(4)*lbd(5))
! !     A(5,2)=-( lbd(1)*lbd(2)*lbd(3) + lbd(1)*lbd(2)*lbd(4) + lbd(1)*lbd(2)*lbd(5) + &
! !               lbd(1)*lbd(3)*lbd(4) + lbd(1)*lbd(3)*lbd(5) + &
! !               lbd(2)*lbd(3)*lbd(4) + lbd(2)*lbd(3)*lbd(5) + &
! !               lbd(2)*lbd(4)*lbd(5) + &
! !               lbd(3)*lbd(4)*lbd(5) )
! !     A(5,3)= ( lbd(1)*lbd(2) + lbd(1)*lbd(3) + lbd(1)*lbd(4) + lbd(1)*lbd(5) + &
! !               lbd(2)*lbd(3) + lbd(2)*lbd(4) + lbd(2)*lbd(5) + &
! !               lbd(3)*lbd(4) + lbd(3)*lbd(5) + &
! !               lbd(4)*lbd(5) )
! !     A(5,4)=-( lbd(1) + lbd(2) + lbd(3) + lbd(4) + lbd(5) )
! !     A(5,5)=1.d0
! !     tmp = (lbd(1)-lbd(2))*(lbd(1)-lbd(3))*(lbd(1)-lbd(4))*(lbd(1)-lbd(5))
! !     W(1,1)=lbd(1)**4.d0/tmp
! !     W(1,2)=lbd(1)**3.d0/tmp
! !     W(1,3)=lbd(1)**2.d0/tmp
! !     W(1,4)=lbd(1)**1.d0/tmp
! !     W(1,5)=lbd(1)**0.d0/tmp
! !     tmp = (lbd(2)-lbd(1))*(lbd(2)-lbd(3))*(lbd(2)-lbd(4))*(lbd(2)-lbd(5))
! !     W(2,1)=lbd(2)**4.d0/tmp
! !     W(2,2)=lbd(2)**3.d0/tmp
! !     W(2,3)=lbd(2)**2.d0/tmp
! !     W(2,4)=lbd(2)**1.d0/tmp
! !     W(2,5)=lbd(2)**0.d0/tmp
! !     tmp = (lbd(3)-lbd(1))*(lbd(3)-lbd(2))*(lbd(3)-lbd(4))*(lbd(3)-lbd(5))
! !     W(3,1)=lbd(3)**4.d0/tmp
! !     W(3,2)=lbd(3)**3.d0/tmp
! !     W(3,3)=lbd(3)**2.d0/tmp
! !     W(3,4)=lbd(3)**1.d0/tmp
! !     W(3,5)=lbd(3)**0.d0/tmp
! !     tmp = (lbd(4)-lbd(1))*(lbd(4)-lbd(2))*(lbd(4)-lbd(3))*(lbd(4)-lbd(5))
! !     W(4,1)=lbd(4)**4.d0/tmp
! !     W(4,2)=lbd(4)**3.d0/tmp
! !     W(4,3)=lbd(4)**2.d0/tmp
! !     W(4,4)=lbd(4)**1.d0/tmp
! !     W(4,5)=lbd(4)**0.d0/tmp
! !     tmp = (lbd(5)-lbd(1))*(lbd(5)-lbd(2))*(lbd(5)-lbd(3))*(lbd(5)-lbd(4))
! !     W(5,1)=lbd(5)**4.d0/tmp
! !     W(5,2)=lbd(5)**3.d0/tmp
! !     W(5,3)=lbd(5)**2.d0/tmp
! !     W(5,4)=lbd(5)**1.d0/tmp
! !     W(5,5)=lbd(5)**0.d0/tmp
! !     call mx_x_mx(5,W,A,V_inv)
! !     out(1)=0.d0
! !     out(2)=0.d0
! !     out(3)=2.d0
! !     out(4)=0.d0
! !     out(5)=0.d0
! !     call mx_x_vec(5,1,5,V_inv,out,coeff)
! !     DLU(i ,i-2)= coeff(1)
! !     DLU(i ,i-1)= coeff(2)
! !     DLU(i ,i  )= coeff(3)
! !     DLU(i ,i+1)= coeff(4)
! !     DLU(i ,i+2)= coeff(5)
! !   enddo
! !
! ! endif
! !
! ! ! NGC=(pencil-1)/2
! ! !
! ! ! DLU=0.d0
! ! !
! ! ! if      (pencil==3) then
! ! !
! ! !
! ! !
! ! !   ! DLU(1 ,N1)= 1.d0
! ! !   DLU(1 ,1 )=0.d0 ! -2.d0
! ! !   DLU(1 ,2 )=0.d0 ! 1.d0
! ! !   do i=NGC+1,N1-NGC
! ! !   DLU(i ,i-1)= (x(i+1)-x(i  ))/((x(i)-x(i-1))*(x(i+1)-x(i))*0.5d0*(x(i+1)-x(i-1)))
! ! !   DLU(i ,i  )=-(x(i+1)-x(i-1))/((x(i)-x(i-1))*(x(i+1)-x(i))*0.5d0*(x(i+1)-x(i-1)))
! ! !   DLU(i ,i+1)= (x(i  )-x(i-1))/((x(i)-x(i-1))*(x(i+1)-x(i))*0.5d0*(x(i+1)-x(i-1)))
! ! !   ! print*, DLU(i ,i-1), DLU(i ,i  ), DLU(i ,i+1)
! ! !   ! DLU(i ,i-1)= 1.d0
! ! !   ! DLU(i ,i  )=-2.d0
! ! !   ! DLU(i ,i+1)= 1.d0
! ! !   enddo
! ! !   ! DLU(N1,1   )= 1.d0
! ! !   DLU(N1,N1-1)=0.d0 ! 1.d0
! ! !   DLU(N1,N1  )=0.d0 ! -2.d0
! ! ! else if (pencil==5) then
! ! !   DLU(1 ,1 )=-5.d0/2.d0
! ! !   DLU(1 ,2 )= 4.d0/3.d0
! ! !   DLU(1 ,3 )=-1.d0/12.d0
! ! !   DLU(2 ,1 )= 4.d0/3.d0
! ! !   DLU(2 ,2 )=-5.d0/2.d0
! ! !   DLU(2 ,3 )= 4.d0/3.d0
! ! !   DLU(2 ,4 )=-1.d0/12.d0
! ! !   do i=NGC+1,N1-NGC
! ! !   DLU(i ,i-2)=-1.d0/12.d0
! ! !   DLU(i ,i-1)= 4.d0/3.d0
! ! !   DLU(i ,i  )=-5.d0/2.d0
! ! !   DLU(i ,i+1)= 4.d0/3.d0
! ! !   DLU(i ,i+2)=-1.d0/12.d0
! ! !   enddo
! ! !   DLU(N1-1,N1-3)=-1.d0/12.d0
! ! !   DLU(N1-1,N1-2)= 4.d0/3.d0
! ! !   DLU(N1-1,N1-1)=-5.d0/2.d0
! ! !   DLU(N1-1,N1  )= 4.d0/3.d0
! ! !   DLU(N1,N1-2)=-1.d0/12.d0
! ! !   DLU(N1,N1-1)= 4.d0/3.d0
! ! !   DLU(N1,N1  )=-5.d0/2.d0
! ! ! endif
! ! ! ! DLU=DLU/(x(2)-x(1))**2.d0
! !
! ! end subroutine get_laplacian
! ! ! -----------------------------------------------------------------------------
!
! ! -----------------------------------------------------------------------------
! !> Decompose Laplacian operator into diagonal, triangular inferior and
! !> triangular superior parts
! ! -----------------------------------------------------------------------------
! subroutine dcmps_laplacian(N1,DLU,D,L,U)
! integer, intent(in) :: N1
! double precision, intent(in) :: DLU(N1,N1)
! double precision, intent(out) :: D(N1,N1), L(N1,N1), U(N1,N1)
! integer :: i, j
! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! D=0.d0
! L=0.d0
! U=0.d0
!
! do i=1,N1
!   do j=1,i-1
!     L(i,j)=DLU(i,j)
!   enddo
!   D(i,i)=DLU(i,i)
!   do j=i+1,N1
!     U(i,j)=DLU(i,j)
!   enddo
! enddo
!
! end subroutine dcmps_laplacian
! ! -----------------------------------------------------------------------------


end module mod_poisson
