!> Random number generators
module mod_rdm
implicit none

public :: init_rdm_seed

contains

! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
subroutine init_rdm_seed(id)
!> id of processor
integer, intent(in)  :: id
integer :: i, n, clock
integer, dimension(:), allocatable :: seed
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
call RANDOM_SEED(size = n)
allocate(seed(n))
call SYSTEM_CLOCK(COUNT=clock)
seed = clock + 37 * (/ (i - 1 + id , i = 1, n) /)
call RANDOM_SEED(PUT = seed)
deallocate(seed)
end subroutine init_rdm_seed
! -----------------------------------------------------------------------------

end module mod_rdm
