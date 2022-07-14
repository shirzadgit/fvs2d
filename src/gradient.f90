module gradient

  use data_grid,  only  : ncells
  use input
  use gradient_ggcb
  use gradient_ggnb
  use gradient_lsq

  implicit none

contains

  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine gradient_init
    implicit  none

    !--------------------------------------------------------------------------!
    ! setup gradient opertor at cell centers: Green-Gauss cell-based (ggcb)
    !--------------------------------------------------------------------------!
    if (lgrad_ggcb) then
      !-- setup gradient opertor at cell centers: Green-Gauss cell-based (ggcb)
      call grad_ggcb_init

    elseif (lgrad_ggnb) then
      !-- setup gradient opertor at cell centers: Green-Gauss cell-based (ggcb)
      call grad_ggnb_init

    elseif (lgrad_lsq) then
      !-- setup gradient opertor at cell centers: Green-Gauss cell-based (ggcb)
      call grad_lsq_init
    endif

    return
  end subroutine gradient_init


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine compute_gradient_cellcntr
    implicit  none


    if (lgrad_ggnb) then
      call grad_ggnb
      !call grad_ggnb_exp(fc,dfc)

    elseif (lgrad_ggcb) then
      call grad_ggcb

    elseif (lgrad_lsq) then
      call grad_lsq

    endif

    return
  end subroutine compute_gradient_cellcntr


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine gradient_cellcntr_1var (fc,dfc)
    implicit  none

    real,intent(in)   :: fc(ncells)
    real,intent(out)  :: dfc(ncells,2)

    dfc(:,:)=0.d0

    if (lgrad_ggnb) then
      call grad_ggnb_1var (fc,dfc)
      !call grad_ggnb_exp(fc,dfc)

    elseif (lgrad_ggcb) then
      call grad_ggcb_1var (fc,dfc)

    elseif (lgrad_lsq) then
      call grad_lsq_1var (fc,dfc)

    endif

    return
  end subroutine gradient_cellcntr_1var

end module gradient
