!> Input parameters
module mod_params
use mod_basic_types
use mod_csts

implicit none

public :: init_par

! - - -

! double precision, parameter :: min1=-0.5d0, max1=0.5d0
double precision, parameter :: min1=1.d0, max1=2.d0
double precision, parameter :: min2=0.d0, max2=dpi

integer, parameter :: N1_O=40
integer, parameter :: N2_O=40

character(len=std_len), parameter :: solver_type="BiCGSTAB"

! WARNING set bc_type in check_par below
character(len=std_len), dimension(4) :: bc_type

character(len=std_len), parameter :: grid_type="uniform"

!> Beware, the smaller the pencil, the more stable BiCGSTAB
!> => trick could be to CV 1st w/ pencil=3 and then,
!> keep going w/ pencil=5 (etc if higher order needed)
integer, parameter :: pencil=3

! - - -

double precision :: kwv

integer :: NGC

integer :: N1_I
integer :: N2_I

integer :: iImin1, iOmin1, iOmax1, iImax1
integer :: iImin2, iOmin2, iOmax2, iImax2

! - - -

contains

! -----------------------------------------------------------------------------
!> Check values of parameter
! -----------------------------------------------------------------------------
subroutine check_par
implicit none

bc_type  = [character(len=std_len) :: "usr","usr","usr","usr"]

if (pencil/=3 .and. pencil/=5) &
  call crash("pencil/=3 not implemented yet")

if (grid_type=='stretched' .and. bc_type(1)=='periodic') &
  call crash("BCs cannot be periodic if stretched grid")

! if (grid_type=='stretched' .and. &
!     ((dabs(min1)<smalldble .or. dabs(max1)<smalldble .or. dabs(min2)<smalldble .or. dabs(max2)<smalldble) .or. &
!      (min1*max1<0.d0 .or. min2*max2<0.d0))) &
!   call crash("Stretched grid incompatible w/ box bounds")

if ((ANY(bc_type(1:2)=='periodic') .and. ANY(bc_type(1:2)/='periodic')).or.&
    (ANY(bc_type(3:4)=='periodic') .and. ANY(bc_type(3:4)/='periodic')))&
  call crash("Beware, both sides should be periodic if periodic BCs are used")

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
N2_I=N2_O+2*NGC
iImin1=1
iOmin1=NGC+1
iOmax1=NGC+N1_O
iImax1=N1_I
iImin2=1
iOmin2=NGC+1
iOmax2=NGC+N2_O
iImax2=N2_I

end subroutine init_par
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Crash after printing an error message
! -----------------------------------------------------------------------------
subroutine crash(message)
character(len=*), intent(in) :: message
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print*, message
stop

end subroutine crash
! -----------------------------------------------------------------------------

end module mod_params
