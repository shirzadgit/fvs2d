module gradient

  use input
  use grid_procs
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
  subroutine gradient_cellcntr (fc,dfc)
    implicit  none

    real,intent(in)   :: fc(num_cells)
    real,intent(out)  :: dfc(num_cells,2)
    integer           :: i,j, ic,ic1, ie,ie1, in,iv,iv1,iv2, nt
    real              :: af,nxf,nyf
    real,allocatable  :: fv(:)

    dfc(:,:)=0.d0

    if (lgrad_ggnb) then
      call grad_ggnb(fc,dfc)
      call grad_ggnb_exp(fc,dfc)

    elseif (lgrad_ggcb) then
      call grad_ggcb(fc,dfc)

    elseif (lgrad_lsq) then
      call grad_lsq (fc,dfc)

    else
      write(*,*) 'error: Green-Gauss node-based and cell-based are only implemented!'
      stop 'error in gradient_cellcntr'
    endif

    return
  end subroutine gradient_cellcntr

end module gradient
