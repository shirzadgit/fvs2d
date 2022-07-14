module gradient_ggcb

  use mainparam,      only  : nvar
  use data_grid
  use data_solution,  only  : pvar,grad

  implicit none

  private

  type ggcb_type
    integer,dimension(:),pointer  :: ptr
    real,dimension(:),pointer     :: coef0
    real,dimension(:,:),pointer   :: coefnb
  end type ggcb_type

  type(ggcb_type),dimension(:),pointer  :: ggcb

  private :: setup
  public  :: grad_ggcb_init, grad_ggcb, grad_ggcb_1var

contains

  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grad_ggcb_init
    implicit  none

    !--------------------------------------------------------------------------!
    ! allocate
    !--------------------------------------------------------------------------!
    allocate(ggcb(ncells))


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

    integer     :: i,j, ic,ic1,ic2, ie,ie1,ie2, je, in,iv,iv1,iv2
    real        :: af,nxf,nyf, xf,yf, xc,yc, xc1,yc1, dx,dy, d0,d1


    do ic=1,ncells
      xc=cell(ic)%x
      yc=cell(ic)%y

      allocate(ggcb(ic)%ptr(cell(ic)%nvrt))
      allocate(ggcb(ic)%coef0(2))
      allocate(ggcb(ic)%coefnb(cell(ic)%nvrt,2))

      ggcb(ic)%ptr(1:cell(ic)%nvrt)=0
      ggcb(ic)%coef0(1:2)=0.d0
      ggcb(ic)%coefnb(1:cell(ic)%nvrt,1:2)=0.d0

      do ie=1,cell(ic)%nvrt
        je=cell(ic)%edge(ie)

        ic1=edge(je)%c1
        ic2=edge(je)%c2

        af =edge(je)%area
        nxf=edge(je)%nx * cell(ic)%nrmlsign(ie)
        nyf=edge(je)%ny * cell(ic)%nrmlsign(ie)

        xf=edge(je)%x
        yf=edge(je)%y

        dx=xf-xc
        dy=yf-yc
        d0=dsqrt(dx**2+dy**2)

        xc1=xc
        yc1=yc
        ggcb(ic)%ptr(ie)=ic
        if (ic1>0 .and. ic1/=ic) then
          xc1=cell(ic1)%x
          yc1=cell(ic1)%y
          ggcb(ic)%ptr(ie)=ic1
        elseif (ic2>0 .and. ic2/=ic) then
          xc1=cell(ic2)%x
          yc1=cell(ic2)%y
          ggcb(ic)%ptr(ie)=ic2
        endif

        dx=xf-xc1
        dy=yf-yc1
        d1=dsqrt(dx**2+dy**2)

        ggcb(ic)%coef0(1) = ggcb(ic)%coef0(1) + d1/(d0+d1)*nxf*af
        ggcb(ic)%coef0(2) = ggcb(ic)%coef0(2) + d1/(d0+d1)*nyf*af
        ggcb(ic)%coefnb(ie,1) = d0/(d0+d1)*nxf*af
        ggcb(ic)%coefnb(ie,2) = d0/(d0+d1)*nyf*af

      enddo
    enddo

    return
  end subroutine setup


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grad_ggcb

    implicit  none
    integer           :: i,j, ic,ic1, ie,ie1, in,iv,iv1,iv2, ivar

    grad = 0.d0

    do ivar=1,nvar
      do ic=1,ncells
        grad(ivar,ic,1)= ggcb(ic)%coef0(1) * pvar(ivar,ic)
        grad(ivar,ic,2)= ggcb(ic)%coef0(2) * pvar(ivar,ic)
        do ie=1,cell(ic)%nvrt
          ic1=ggcb(ic)%ptr(ie)
          grad(ivar,ic,1) = grad(ivar,ic,1) + ggcb(ic)%coefnb(ie,1)*pvar(ivar,ic1)
          grad(ivar,ic,2) = grad(ivar,ic,2) + ggcb(ic)%coefnb(ie,2)*pvar(ivar,ic1)
        end do
        grad(ivar,ic,1) = grad(ivar,ic,1)/cell(ic)%vol
        grad(ivar,ic,2) = grad(ivar,ic,2)/cell(ic)%vol
      end do
    enddo

    return
  end subroutine grad_ggcb


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grad_ggcb_1var (fc,dfc)
    implicit  none

    real,intent(in)   :: fc(ncells)
    real,intent(out)  :: dfc(ncells,2)
    integer           :: i,j, ic,ic1, ie,ie1, in,iv,iv1,iv2, nt

    dfc(:,:)=0.d0

    do ic=1,ncells
      dfc(ic,1)= ggcb(ic)%coef0(1) * fc(ic)
      dfc(ic,2)= ggcb(ic)%coef0(2) * fc(ic)
      do ie=1,cell(ic)%nvrt
        ic1=ggcb(ic)%ptr(ie)
        dfc(ic,1) = dfc(ic,1) + ggcb(ic)%coefnb(ie,1)*fc(ic1)
        dfc(ic,2) = dfc(ic,2) + ggcb(ic)%coefnb(ie,2)*fc(ic1)
      end do
      dfc(ic,1) = dfc(ic,1)/cell(ic)%vol
      dfc(ic,2) = dfc(ic,2)/cell(ic)%vol
    end do

    return
  end subroutine grad_ggcb_1var

end module gradient_ggcb
