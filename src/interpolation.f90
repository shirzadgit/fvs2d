module interpolation

  use data_grid
  use grid_procs

  implicit none

  private

  type intrp_type
     real,dimension(:),pointer     :: idw
  end type intrp_type

  type(intrp_type),dimension(:),pointer :: intrp

  private :: cell2node_idw_setup
  public  :: interpolate_init, interpolate_cell2node

contains


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine interpolate_init
    implicit  none


    !--------------------------------------------------------------------------!
    ! allocate
    !--------------------------------------------------------------------------!
    allocate(intrp(nnodes))


    !--------------------------------------------------------------------------!
    ! setup interpolation stencil based on linear inverse distance
    !--------------------------------------------------------------------------!
    call cell2node_idw_setup


    !--------------------------------------------------------------------------!
    ! setup interpolation stencil based on least-squares
    !--------------------------------------------------------------------------!
    !call cell2node_lsq_setup


    return
  end subroutine interpolate_init


  !============================================================================!
  ! Interpolate cell center values (fc) to cell vertices (fv)
  ! scheme: linear inverse distance weighting
  !         _______
  !        \
  !         \
  ! fv(i) =  .      node(i)%intrp_idw(ic)*fc(ic)
  !         /
  !        /_______
  !          ic=nb
  !============================================================================!
  subroutine cell2node_idw_setup
    implicit  none

    integer           :: i,in,ic
    real              :: idt, xc,yc, xv,yv, dx,dy
    real,allocatable  :: d(:)


    node_loop: do in=1,nnodes

      allocate(intrp(in)%idw(node(in)%ncells), d(node(in)%ncells))

      xv=node(in)%x
      yv=node(in)%y

      idt=0.d0
      do i=1,node(in)%ncells
        ic=node(in)%cell(i)

        xc=cell(ic)%x
        yc=cell(ic)%y

        dx=xc-xv
        dy=yc-yv
        d(i)=dsqrt(dx**2 +dy**2)

        idt=idt + 1.d0/d(i)
      enddo

      do i=1,node(in)%ncells
        intrp(in)%idw(i)=1.d0/d(i)/idt
      enddo

      deallocate(d)
    enddo node_loop


    return
  end subroutine cell2node_idw_setup



  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine interpolate_cell2node (fc,fv)
    implicit  none

    real,intent(in)   :: fc(ncells)
    real,intent(out)  :: fv(nnodes)
    integer           :: in,ic,i

    fv(:)=0.d0
    do in=1,nnodes
      do i=1,node(in)%ncells
        ic=node(in)%cell(i)
        fv(in) = fv(in) + intrp(in)%idw(i)*fc(ic)
      enddo
    enddo

    return
  end subroutine interpolate_cell2node

end module interpolation
