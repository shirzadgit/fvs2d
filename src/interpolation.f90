module interpolation

  use grid_procs

  implicit none

  real,allocatable,dimension(:),save  :: interp_nodeweight, interp_cellweight

contains


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine interpolate_init
    implicit  none


    !--------------------------------------------------------------------------!
    ! setup interpolation opertor for cell centers to cell nodes interpolation
    !--------------------------------------------------------------------------!
    call interpolate_cellcntr2node_setup


    !--------------------------------------------------------------------------!
    ! check interpolation
    !--------------------------------------------------------------------------!
    !call interpolate_check


    return
  end subroutine interpolate_init


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine interpolate_cellcntr2node_setup
    implicit  none

    integer   :: i,in,ic,iv,nt
    real      :: d,dt

    nt=size(node2cell,dim=1)
    allocate(interp_nodeweight(num_nodes), interp_cellweight(nt))

    !--------------------------------------------------------------------------!
    ! Interpolate cell center values (fc) to nodes (fv) 
    !                              _______
    !                             \
    !                              \
    ! fv(i) = interp_nodeweight(i)  .       interp_cellweight(ic)*fc(ic)
    !                               /
    !                              /_______
    !                                 ic=nb
    !--------------------------------------------------------------------------!
    do in=1,num_nodes
      dt=0.d0
      do i=node2cell_ptr(in), node2cell_ptr(in) + node2cell_ntot(in) - 1
        ic=node2cell(i)
        do iv=1,num_vert_cell(ic)
          if (in==cell2node(ic,iv)) then
            d=cell_dist2vert(ic,iv)
            dt=dt+1.d0/d
            interp_cellweight(i)=1.0/d
          endif
        enddo
      enddo
      interp_nodeweight(in)=1.d0/dt;
    enddo


    return
  end subroutine interpolate_cellcntr2node_setup



  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine interpolate_cellcntr2node (fc,fv)
    implicit  none

    real,intent(in)   :: fc(num_cells)
    real,intent(out)  :: fv(num_nodes)
    !character(len=*)  :: interp_type
    integer           :: in,ic,i

    fv(:)=0.d0
    do in=1,num_nodes
      do i=node2cell_ptr(in), node2cell_ptr(in) + node2cell_ntot(in) - 1
        ic=node2cell(i)
        fv(in) = fv(in) + interp_cellweight(i)*fc(ic)
      enddo
      fv(in) = interp_nodeweight(in)*fv(in)
    enddo

    return
  end subroutine interpolate_cellcntr2node

end module interpolation
