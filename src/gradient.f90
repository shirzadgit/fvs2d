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
  subroutine gradient_cellcntr_test (fc,dfc)
    implicit  none

    real,intent(in)   :: fc(ncells)
    real,intent(out)  :: dfc(ncells,2)

    dfc(:,:)=0.d0

    if (lgrad_ggnb) then
      call grad_ggnb_test (fc,dfc)
      !call grad_ggnb_exp(fc,dfc)

    elseif (lgrad_ggcb) then
      call grad_ggcb_test (fc,dfc)

    elseif (lgrad_lsq) then
      call grad_lsq_test (fc,dfc)
      ! if (lgrad_lsq_fn) then
      !   call grad_lsq_fn(fc,dfc)
      !
      ! elseif (lgrad_lsq_nn) then
      !   call grad_lsq_nn (fc,dfc)
      ! endif
    endif

    return
  end subroutine gradient_cellcntr_test

end module gradient
