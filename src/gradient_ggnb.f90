module gradient_ggnb

  use mainparam,      only  : nvar
  use data_grid
  use data_solution,  only  : pvar, grad
  use interpolation

  implicit none

  private

  type ggnb_type
    real,dimension(:),pointer     :: coef0, coefnb
    real,dimension(:,:),pointer   :: coefedg
  end type ggnb_type

  type(ggnb_type),dimension(:),pointer  :: ggnb

  private :: setup
  public  :: grad_ggnb_init, grad_ggnb, grad_ggnb_exp, grad_ggnb_1var

contains

  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grad_ggnb_init
    implicit  none

    !--------------------------------------------------------------------------!
    ! allocation
    !--------------------------------------------------------------------------!
    allocate(ggnb(ncells))


    !--------------------------------------------------------------------------!
    ! setup gradient opertor at cell centers: Green-Gauss node-based (ggnb)
    !--------------------------------------------------------------------------!
    call setup


    return
  end subroutine grad_ggnb_init


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine setup
    implicit  none

    integer           :: i,j, ic,jc,ic2, ie,ie1,ie2, je, in,iv,iv1,iv2,jv
    real              :: af,nxf,nyf, xc,yc, xv,yv, dx,dy, w1,w2,wt1,wt2, d,idt
    real,allocatable  :: intrp_idw(:)

    !--------------------------------------------------------------------------!
    ! compute intrep coefficent
    !--------------------------------------------------------------------------!
    allocate(intrp_idw(nnodes))
    do in=1,nnodes
      xv=node(in)%x
      yv=node(in)%y

      idt=0.d0
      do i=1,node(in)%ncells
        ic=node(in)%cell(i)

        xc=cell(ic)%x
        yc=cell(ic)%y

        dx=xc-xv
        dy=yc-yv
        d=dsqrt(dx**2 +dy**2)

        idt=idt + 1.d0/d
      enddo

      intrp_idw(in)=1.d0/idt
    enddo


    !--------------------------------------------------------------------------!
    ! compute gradient coeff for cell
    !--------------------------------------------------------------------------!
    do ic=1,ncells
      allocate(ggnb(ic)%coef0(2))
      ggnb(ic)%coef0(1:2)=0.d0
      do ie=1, cell(ic)%nvrt
        je=cell(ic)%edge(ie)

        af =edge(je)%area
        nxf=edge(je)%nx * cell(ic)%nrmlsign(ie)
        nyf=edge(je)%ny * cell(ic)%nrmlsign(ie)

        iv1=edge(je)%n1
        iv2=edge(je)%n2

        dx = cell(ic)%x - node(iv1)%x
        dy = cell(ic)%y - node(iv1)%y
        w1 = 1.d0/dsqrt(dx**2 + dy**2)

        dx = cell(ic)%x - node(iv2)%x
        dy = cell(ic)%y - node(iv2)%y
        w2 = 1.d0/dsqrt(dx**2 + dy**2)

        wt1 = intrp_idw(iv1)
        wt2 = intrp_idw(iv2)

        ggnb(ic)%coef0(1)= ggnb(ic)%coef0(1) + af*nxf/2.d0 * (w1*wt1 + w2*wt2)
        ggnb(ic)%coef0(2)= ggnb(ic)%coef0(2) + af*nyf/2.d0 * (w1*wt1 + w2*wt2)
      enddo
    enddo


    !--------------------------------------------------------------------------!
    ! compute gradient coeff for neighboring cells
    !--------------------------------------------------------------------------!
    do ic=1,ncells
      i=0
      do in=1,cell(ic)%nvrt
        iv=cell(ic)%node(in)
        i = i + node(iv)%ncells
      enddo
      allocate( ggnb(ic)%coefnb(i) );
    enddo


    do ic=1,ncells
      ggnb(ic)%coefnb(:) =0.d0
      i=0
      do in=1,cell(ic)%nvrt
        iv=cell(ic)%node(in)
        xv=node(iv)%x
        yv=node(iv)%y
        do j=1,node(iv)%ncells
          i=i+1
          jc=node(iv)%cell(j)
          xc=cell(jc)%x
          yc=cell(jc)%y
          dx=xc-xv
          dy=yc-yv
          d =dsqrt(dx**2+dy**2)
          ggnb(ic)%coefnb(i)=intrp_idw(iv)/d
          if (jc==ic) ggnb(ic)%coefnb(i)=0.d0
        enddo
      enddo
    enddo


    !--------------------------------------------------------------------------!
    ! computing edge coeff of each node for a cell
    !--------------------------------------------------------------------------!
    do ic=1,ncells
      allocate( ggnb(ic)%coefedg(cell(ic)%nvrt,2) )
      ggnb(ic)%coefedg(:,:) = 0.d0
      do in=1,cell(ic)%nvrt
        iv=cell(ic)%node(in)
        do ie=1,cell(ic)%nvrt
          je=cell(ic)%edge(ie)
          do j=1,2
            if (j==1) jv=edge(je)%n1
            if (j==2) jv=edge(je)%n2
            if (jv==iv) then
              af =edge(je)%area
              nxf=edge(je)%nx * cell(ic)%nrmlsign(ie)
              nyf=edge(je)%ny * cell(ic)%nrmlsign(ie)

              ggnb(ic)%coefedg(in,1)=ggnb(ic)%coefedg(in,1) + af*nxf/2.d0
              ggnb(ic)%coefedg(in,2)=ggnb(ic)%coefedg(in,2) + af*nyf/2.d0
            endif
          enddo
        enddo
      enddo
    enddo

    return
  end subroutine setup


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grad_ggnb

    implicit  none
    integer           :: i,j, ic,jc, ie,je, in,iv, ivar

    grad = 0.d0

    do ivar=1,nvar
      do ic=1,ncells
        grad(ivar,ic,1)= ggnb(ic)%coef0(1) * pvar(ivar,ic)
        grad(ivar,ic,2)= ggnb(ic)%coef0(2) * pvar(ivar,ic)
        i=0
        do in=1,cell(ic)%nvrt
          iv=cell(ic)%node(in)
          do j=1,node(iv)%ncells
            jc=node(iv)%cell(j)
            i=i+1
            grad(ivar,ic,1) = grad(ivar,ic,1) + ggnb(ic)%coefedg(in,1) * ggnb(ic)%coefnb(i) * pvar(ivar,jc)
            grad(ivar,ic,2) = grad(ivar,ic,2) + ggnb(ic)%coefedg(in,2) * ggnb(ic)%coefnb(i) * pvar(ivar,jc)
          end do
        end do
        grad(ivar,ic,1) = grad(ivar,ic,1)/cell(ic)%vol
        grad(ivar,ic,2) = grad(ivar,ic,2)/cell(ic)%vol
      end do
    enddo

    return
  end subroutine grad_ggnb


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grad_ggnb_1var (fc,dfc)
    implicit  none

    real,intent(in)   :: fc(ncells)
    real,intent(out)  :: dfc(ncells,2)
    integer           :: i,j, ic,jc, ie,je, in,iv

    dfc(:,:)=0.d0

    do ic=1,ncells
      dfc(ic,1)= ggnb(ic)%coef0(1) * fc(ic)
      dfc(ic,2)= ggnb(ic)%coef0(2) * fc(ic)
      i=0
      do in=1,cell(ic)%nvrt
        iv=cell(ic)%node(in)
        do j=1,node(iv)%ncells
          jc=node(iv)%cell(j)
          i=i+1
          dfc(ic,1) = dfc(ic,1) + ggnb(ic)%coefedg(in,1) * ggnb(ic)%coefnb(i) * fc(jc)
          dfc(ic,2) = dfc(ic,2) + ggnb(ic)%coefedg(in,2) * ggnb(ic)%coefnb(i) * fc(jc)
        end do
      end do
      dfc(ic,1) = dfc(ic,1)/cell(ic)%vol
      dfc(ic,2) = dfc(ic,2)/cell(ic)%vol
    end do

    return
  end subroutine grad_ggnb_1var


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grad_ggnb_exp (fc,dfc)
    implicit  none

    real,intent(in)   :: fc(ncells)
    real,intent(out)  :: dfc(ncells,2)
    integer           :: i,j, ic,jc, ie,je, in,iv,iv1,iv2
    real              :: af,nxf,nyf
    real,allocatable  :: fv(:)

    dfc(:,:)=0.d0

    allocate(fv(nnodes))

    call interpolate_cell2node(fc,fv)

    do ic=1,ncells
      do ie=1,cell(ic)%nvrt
        je=cell(ic)%edge(ie)

        af =edge(je)%area
        nxf=edge(je)%nx * cell(ic)%nrmlsign(ie)
        nyf=edge(je)%ny * cell(ic)%nrmlsign(ie)

        iv1=edge(je)%n1
        iv2=edge(je)%n2

        dfc(ic,1) = dfc(ic,1) + nxf*af * (fv(iv1)+fv(iv2))/2.d0
        dfc(ic,2) = dfc(ic,2) + nyf*af * (fv(iv1)+fv(iv2))/2.d0
      enddo
      dfc(ic,1) = dfc(ic,1)/cell(ic)%vol
      dfc(ic,2) = dfc(ic,2)/cell(ic)%vol
    enddo

    deallocate(fv)

    return
  end subroutine grad_ggnb_exp


end module gradient_ggnb
