module mod_share

implicit none

private

public :: share_workload

contains


!=======================================================================
!> Compute the # of particles N_pcl processor id should handled to
!> optimize workload balance given the total # of particles N_pcl_tot in
!> simulation space and the total # of processors NPROC.
!> This subroutine itself is unbalanced since proc #i needs to
!> compute the # of particles handled by the previous processors.
!=======================================================================
subroutine share_workload(id,NPROC,N_pcl_tot,N_pcl)
! ----------------------------------------------------------------------
integer, intent(in)  :: id, NPROC, N_pcl_tot
integer, intent(out) :: N_pcl
integer :: i
! ----------------------------------------------------------------------
! Balanced workload - - - -
N_pcl=N_pcl_tot/NPROC
! Unbalanced workload - - -
! N_pcl=N_pcl_tot-(NPROC-1)
if (id>0) then
  do i=1,id-1
    N_pcl=N_pcl+(N_pcl_tot-N_pcl)/(NPROC-i)
  enddo
  N_pcl=(N_pcl_tot-N_pcl)/(NPROC-id)
endif
end subroutine share_workload
!=======================================================================


end module mod_share
