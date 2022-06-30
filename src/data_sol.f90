module data_sol

  use mainparam,  only  : nvar
  use data_grid,  only  : ncells, nedges, nnodes

  implicit none

  ! !----------------------------------------------------------------------------!
  ! ! data type for solution at cell centers
  ! !----------------------------------------------------------------------------!
  ! type cell_sol
  !   real                          :: r,u,v,p
  ! end type cell_sol
  !
  ! !----------------------------------------------------------------------------!
  ! ! data type for solution at cell vertices
  ! !----------------------------------------------------------------------------!
  ! type node_sol
  !   real                          :: r,u,v,p
  ! end type node_sol

  !----------------------------------------------------------------------------!
  ! allocatable data
  !----------------------------------------------------------------------------!
  ! type(cell_sol),dimension(:),pointer  :: pvar, cvar
  ! type(node_sol),dimension(:),pointer  :: npvar
  real,allocatable    :: pvar(:,:)



contains


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine data_sol_allocate
    implicit none

    ! allocate(pvar(ncells))
    ! allocate(npvar(nnodes))

    allocate(pvar(nvar,ncells))


    return
  end subroutine data_sol_allocate


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine data_sol_close
    implicit none

    deallocate(pvar)
    !deallocate(npvar)

    return
  end subroutine data_sol_close

end module data_sol
