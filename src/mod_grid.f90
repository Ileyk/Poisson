!> Grid (short description)
!>
!> Longer description of the module.
!> This *module* gathers **all subroutines** linked to the grid.
!>
!> ### Example
!> ~~~~~~~~~~~.f90
!> call make_grid(200,-0.5d0,0.5d0,x)
!> ~~~~~~~~~~~
!>
!> @todo
!>  - 2D and 3D
!>  - cylindrical, spherical and more
!>  - adaptive meshes

module mod_grid

implicit none

public :: make_grid
public :: get_dx

contains

! -----------------------------------------------------------------------------
!> Make 1D uniform Cartesian mesh
!>
!> Longer description, maybe even lenghty, which could span several lines
!> and make the reader asleep very quickly. This is similar to
!> mod_poisson::frst_gss().
!> @param N1 number of cells in direction # 1
!> @param NGC number of ghost cells in direction # 1
!> @param iImin1 minimum index including ghost cells
!> @param iImax1 maximum index including ghost cells
!> @param min1 minimum edge of direction # 1
!> @param max1 maximum edge of direction # 1
!> @returns array of positions for cell centers
!> @see mod_init::init_mpi() and mod_poisson::get_analytic()
!> @note PIF
!> @attention PAF
!> @warning POUF
! -----------------------------------------------------------------------------
subroutine make_grid(N1,N2,NGC,iImin1,iImax1,iImin2,iImax2,min1,max1,min2,max2,grid_type,x)
use mod_io
integer, intent(in) :: N1
integer, intent(in) :: N2
integer, intent(in) :: NGC
integer, intent(in) :: iImin1, iImax1
integer, intent(in) :: iImin2, iImax2
double precision, intent(in) :: min1, max1
double precision, intent(in) :: min2, max2
character(len=*), intent(in) :: grid_type
double precision, allocatable, intent(out) :: x(:,:,:)
integer :: i, j
double precision :: dx(2)
double precision :: q ! , dxold
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

allocate(x(N1,N2,2))

if      (grid_type=='uniform') then

  dx(1)=(max1-min1)/dble(N1-2*NGC)
  dx(2)=(max2-min2)/dble(N2-2*NGC)
  do i=iImin1,iImax1
  do j=iImin2,iImax2
    x(i,j,1)=min1+dx(1)*dble(i-0.5d0-NGC)
    x(i,j,2)=min2+dx(2)*dble(j-0.5d0-NGC)
  enddo
  enddo

else if (grid_type=='stretched') then

  q=(max1/min1)**(1.d0/dble(N1-2*NGC))
  ! Lower ghost cells - - -
  dx(1)=min1*(q-1.d0)
  dx(1)=dx(1)/q
  x(NGC,:,1)=min1-0.5d0*dx(1)
  do i=NGC-1,1,-1
    x(i,:,1)=x(i+1,:,1)-0.5d0*dx(1)-0.5d0*dx(1)/q
    dx(1)=dx(1)/q
  enddo
  ! - - -
  dx(1)=min1*(q-1.d0)
  dx(1)=dx(1)/q
  do i=iImin1+1,iImax1
    x(i,:,1)=x(i-1,:,1)+0.5d0*dx(1)+0.5d0*dx(1)*q
    dx(1)=dx(1)*q
  enddo

  q=(max2/min2)**(1.d0/dble(N2-2*NGC))
  ! Lower ghost cells - - -
  dx(2)=min2*(q-1.d0)
  dx(2)=dx(2)/q
  x(:,NGC,2)=min2-0.5d0*dx(2)
  do i=NGC-1,1,-1
    x(:,i,2)=x(:,i+1,2)-0.5d0*dx(2)-0.5d0*dx(2)/q
    dx(2)=dx(2)/q
  enddo
  ! - - -
  dx(2)=min2*(q-1.d0)
  dx(2)=dx(2)/q
  do i=iImin2+1,iImax2
    x(:,i,2)=x(:,i-1,2)+0.5d0*dx(2)+0.5d0*dx(2)*q
    dx(2)=dx(2)*q
  enddo

endif

! print*, x(1), x(2), x(N1-1), x(N1)
! stop

call save_grid(N1,N2,0,'x_',x)

end subroutine make_grid
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Get the spatial step in each cell of whole simulation space
!>
!> @warning Include ghost cells
! -----------------------------------------------------------------------------
subroutine get_dx(N1,N2,NGC,x,min1,max1,min2,max2,grid_type,dx)
integer, intent(in) :: N1
integer, intent(in) :: N2
integer, intent(in) :: NGC
double precision, intent(in) :: x(N1,N2,2)
double precision, intent(in) :: min1, max1
double precision, intent(in) :: min2, max2
character(len=*), intent(in) :: grid_type
double precision, intent(out) :: dx(N1,N2,2)
integer :: i, j
double precision :: q, dx0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if      (grid_type=='uniform') then

  dx(:,:,1)=x(2,1,1)-x(1,1,1)
  dx(:,:,2)=x(1,2,2)-x(1,1,2)

else if (grid_type=='stretched') then

  q=(max1/min1)**(1.d0/dble(N1-2*1))
  ! Lower ghost cells
  dx0=min1*(q-1.d0)
  do i=1,NGC
    dx(i,:,1)=dx0/q**(dble(NGC-i+1))
  enddo
  do i=NGC+1,N1
    dx(i,:,1)=dx(i-1,:,1)*q
  enddo

  q=(max2/min2)**(1.d0/dble(N2-2*1))
  ! Lower ghost cells
  dx0=min2*(q-1.d0)
  do i=1,NGC
    dx(:,i,2)=dx0/q**(dble(NGC-i+1))
  enddo
  do i=NGC+1,N2
    dx(:,i,2)=dx(:,i-1,2)*q
  enddo

endif

end subroutine get_dx
! -----------------------------------------------------------------------------


end module mod_grid
