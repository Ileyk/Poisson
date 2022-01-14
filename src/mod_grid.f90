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
subroutine make_grid(N1,NGC,iImin1,iImax1,min1,max1,x)
integer, intent(in) :: N1
integer, intent(in) :: NGC
integer, intent(in) :: iImin1, iImax1
double precision, intent(in) :: min1, max1
double precision, allocatable, intent(out) :: x(:)
integer :: i
double precision :: dx
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

allocate(x(N1))
dx=(max1-min1)/dble(N1-2*NGC)
do i=iImin1,iImax1
  x(i)=min1+dx*(i-0.5-NGC)
enddo

end subroutine make_grid
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Get the spatial step in each cell of whole simulation space
!>
!> @warning Include ghost cells
! -----------------------------------------------------------------------------
subroutine get_dx(N1,x,dx)
integer, intent(in) :: N1
double precision, intent(in) :: x(N1)
double precision, intent(out) :: dx(N1)
integer :: i
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

dx(1)=x(2)-x(1)
do i=2,N1-1
  dx(i)=x(i+1)-x(i)
enddo
dx(N1)=x(N1)-x(N1-1)
! dx=1.d-2 ! tmp

end subroutine get_dx
! -----------------------------------------------------------------------------


end module mod_grid
