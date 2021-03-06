!> Input parameters
module mod_params
use mod_basic_types

implicit none

public :: init_par

! - - -

! double precision, parameter :: min1=-0.5d0, max1=0.5d0
double precision, parameter :: min1=1.d0, max1=2.d0

integer, parameter :: N1_O=100

character(len=std_len), parameter :: solver_type="BiCGSTAB"

character(len=std_len), parameter :: bc_type="usr"

character(len=std_len), parameter :: grid_type="stretched"

!> Beware, the smaller the pencil, the more stable BiCGSTAB
!> => trick could be to CV 1st w/ pencil=3 and then,
!> keep going w/ pencil=5 (etc if higher order needed)
integer, parameter :: pencil=3

! - - -

double precision :: kwv

integer :: NGC

integer :: N1_I

integer :: iImin1, iOmin1, iOmax1, iImax1

! - - -

contains

! -----------------------------------------------------------------------------
!> Check values of parameter
! -----------------------------------------------------------------------------
subroutine check_par
implicit none

if (pencil/=3 .and. pencil/=5) &
  call mpistop("pencil/=3,5 not implemented yet")

if (grid_type=='stretched' .and. bc_type=='periodic') &
  call mpistop("BCs cannot be periodic if stretched grid")

end subroutine check_par
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Initialize parameters deduced from input ones
! -----------------------------------------------------------------------------
subroutine init_par
implicit none

kwv=max1-min1
NGC=(pencil-1)/2
N1_I=N1_O+2*NGC
iImin1=1
iOmin1=NGC+1
iOmax1=NGC+N1_O
iImax1=N1_I

end subroutine init_par
! -----------------------------------------------------------------------------

end module mod_params
