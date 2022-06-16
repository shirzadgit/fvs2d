module interpolation

  use grid_procs


contains

  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine interpolate_init
    implicit  none



    !--------------------------------------------------------------------------!
    ! setup interpolation opertor for cell centers to cell nodes interpolation
    !--------------------------------------------------------------------------!
    call interpolate_setup_cellcntr2node


    return
  end subroutine interpolate_init


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine interpolate_setup_cellcntr2node
    implicit  none




    return
  end subroutine interpolate_setup_cellcntr2node



  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine interpolate_cellcntr2node
    implicit  none





    return
  end subroutine interpolate_cellcntr2node




end module interpolation
