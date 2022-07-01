module data_sol

  use mainparam,  only  : nvar
  use data_grid,  only  : ncells, nedges, nnodes

  implicit none

  ! Primitive variables
  ! ivar=1, density
  ! ivar=2, u-velocity
  ! ivar=3, v-velocity
  ! ivar=4, pressure

  real,allocatable    :: pvar(:,:)



contains


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine data_sol_allocate
    implicit none


    allocate(pvar(nvar,ncells))


    return
  end subroutine data_sol_allocate


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine data_sol_close
    implicit none

    deallocate(pvar)

    return
  end subroutine data_sol_close

end module data_sol
