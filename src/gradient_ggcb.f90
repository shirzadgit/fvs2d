module gradient_ggcb

  use grid_procs
  use interpolation

  implicit none

  private
  integer,allocatable,save  :: grad_ggcb_ptr(:,:)
  real,allocatable,save     :: grad_ggcb_coef0(:,:), grad_ggcb_coefnb(:,:,:)

  private :: setup
  public  :: grad_ggcb_init, grad_ggcb

contains

  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grad_ggcb_init
    implicit  none

    !--------------------------------------------------------------------------!
    ! allocate
    !--------------------------------------------------------------------------!
    allocate( grad_ggcb_ptr(num_cells,num_vert_max), grad_ggcb_coef0(num_cells,2), grad_ggcb_coefnb(num_cells,num_vert_max,2))


    !--------------------------------------------------------------------------!
    ! setup gradient opertor at cell centers: Green-Gauss cell-based (ggcb)
    !--------------------------------------------------------------------------!
    call setup

    return
  end subroutine grad_ggcb_init


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine setup
    implicit  none

    integer     :: i,j, ic,ic1,ic2, ie,ie1,ie2, in,iv,iv1,iv2
    real        :: af,nxf,nyf, xf,yf, xc,yc, xc1,yc1, dx,dy, d0,d1


    grad_ggcb_coef0(:,:)=0.d0
    grad_ggcb_coefnb(:,:,:)=0.d0

    do ic=1,num_cells
      xc=cell_center(ic,1)
      yc=cell_center(ic,2)
      do ie=1,num_vert_cell(ic)
        ie1=cell2edge(ic,ie)

        ic1=edge2cell(ie1,1)
        ic2=edge2cell(ie1,2)

        af =edge_area(ie1)
        nxf=edge_normal(ie1,1) * edge_normal_sign(ic,ie)
        nyf=edge_normal(ie1,2) * edge_normal_sign(ic,ie)

        xf=edge_center(ie1,1)
        yf=edge_center(ie1,2)

        dx=xf-xc
        dy=yf-yc
        d0=dsqrt(dx**2+dy**2)

        xc1=xc
        yc1=yc
        grad_ggcb_ptr(ic,ie)=ic
        if (ic1>0 .and. ic1/=ic) then
          xc1=cell_center(ic1,1)
          yc1=cell_center(ic1,2)
          grad_ggcb_ptr(ic,ie)=ic1
        elseif (ic2>0 .and. ic2/=ic) then
          xc1=cell_center(ic2,1)
          yc1=cell_center(ic2,2)
          grad_ggcb_ptr(ic,ie)=ic2
        endif

        dx=xf-xc1
        dy=yf-yc1
        d1=dsqrt(dx**2+dy**2)

        grad_ggcb_coef0(ic,1)  = grad_ggcb_coef0(ic,1) + d1/(d0+d1)*nxf*af
        grad_ggcb_coef0(ic,2)  = grad_ggcb_coef0(ic,2) + d1/(d0+d1)*nyf*af
        grad_ggcb_coefnb(ic,ie,1) = d0/(d0+d1)*nxf*af
        grad_ggcb_coefnb(ic,ie,2) = d0/(d0+d1)*nyf*af
      enddo
    enddo

    return
  end subroutine setup


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grad_ggcb (fc,dfc)
    implicit  none

    real,intent(in)   :: fc(num_cells)
    real,intent(out)  :: dfc(num_cells,2)
    integer           :: i,j, ic,ic1, ie,ie1, in,iv,iv1,iv2, nt
    real              :: af,nxf,nyf
    real,allocatable  :: fv(:)

    dfc(:,:)=0.d0

    do ic=1,num_cells
      dfc(ic,1)= grad_ggcb_coef0(ic,1) * fc(ic)
      dfc(ic,2)= grad_ggcb_coef0(ic,2) * fc(ic)
      do ie=1,num_vert_cell(ic)
        ic1=grad_ggcb_ptr(ic,ie)
        dfc(ic,1) = dfc(ic,1) + grad_ggcb_coefnb(ic,ie,1)*fc(ic1)
        dfc(ic,2) = dfc(ic,2) + grad_ggcb_coefnb(ic,ie,2)*fc(ic1)
      end do
      dfc(ic,1) = dfc(ic,1)/cell_area(ic)
      dfc(ic,2) = dfc(ic,2)/cell_area(ic)
    end do


    return
  end subroutine grad_ggcb

end module gradient_ggcb
