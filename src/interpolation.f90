module interpolation

  use grid_procs

  real,allocatable,dimension(:),save  :: interp_ggnb_nodeweight, interp_ggnb_cellweight

contains


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine interpolate_init
    implicit  none


    !--------------------------------------------------------------------------!
    ! setup interpolation opertor for cell centers to cell nodes interpolation
    !--------------------------------------------------------------------------!
    call interpolate_cellcntr2node_setup_ggnb


    !--------------------------------------------------------------------------!
    ! check interpolation
    !--------------------------------------------------------------------------!
    call interpolate_check


    return
  end subroutine interpolate_init


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine interpolate_cellcntr2node_setup_ggnb
    implicit  none

    integer   :: i,in,ic,iv,nt
    real      :: d,dt
    
    nt=size(node2cell,dim=1)
    allocate(interp_ggnb_nodeweight(num_nodes), interp_ggnb_cellweight(nt))

    !--------------------------------------------------------------------------!
    ! Green-Gauss Node Based Interpolation
    !                                   _______
    !                                  \
    !                                   \
    ! fv(i) = interp_ggnb_nodeweight(i)  .        interp_ggnb_cellweight(ic)*fc(ic)
    !                                   /
    !                                  /_______
    !                                    ic=nb
    !--------------------------------------------------------------------------!
    do in=1,num_nodes
      dt=0.d0
      do i=node2cell_ptr(in), node2cell_ptr(in) + node2cell_ntot(in) - 1
        ic=node2cell(i)
        do iv=1,num_vert_cell(ic)
          if (in==cell2node(ic,iv)) then
            d=cell_dist2vert(ic,iv)
            dt=dt+1.d0/d
            interp_ggnb_cellweight(i)=1.0/d
          endif
        enddo
      enddo
      interp_ggnb_nodeweight(in)=1.d0/dt;
    enddo


    return
  end subroutine interpolate_cellcntr2node_setup_ggnb



  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine interpolate_cellcntr2node (interp_type,fc,fv)
    implicit  none

    real,intent(in)   :: fc(num_cells)
    real,intent(out)  :: fv(num_nodes)
    character(len=*)  :: interp_type
    integer           :: in,ic,i

    if (interp_type=='ggnb') then
      fv(:)=0.d0
      do in=1,num_nodes
        do i=node2cell_ptr(in), node2cell_ptr(in) + node2cell_ntot(in) - 1
          ic=node2cell(i)
          fv(in) = fv(in) + interp_ggnb_cellweight(i)*fc(ic)
        enddo
        fv(in) = interp_ggnb_nodeweight(in)*fv(in)
      enddo
    else
      write(*,*) 'only green-gauss node-based is implemented'
      stop 'interpolate_cellcntr2node'
    endif

    return
  end subroutine interpolate_cellcntr2node



  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine interpolate_check
    implicit none
    integer             :: ic,ic1,ic2, in,in1,in2, ie,ie1,ie2, i,j,n1,nt,iv,iv1,iv2, num_nearests
    integer,allocatable :: iwrk1d(:),iwrk2d(:,:)
    real,allocatable    :: fc(:),fve(:),fvi(:)
    real                :: vol, r(2), d,dt,xc,yc
    type(kdtree2_result),allocatable	:: fnearest(:)

    allocate(fc(num_cells),fve(num_nodes), fvi(num_nodes))
    do in=1,num_nodes
      xc=cell_vert(in,1)
      yc=cell_vert(in,2)
      fve(in)=sin(1.25d0*xc)*cos(1.25d0*yc)
    enddo

    do ic=1,num_cells
      xc=cell_center(ic,1)
      yc=cell_center(ic,2)
      fc(ic)=sin(1.25d0*xc)*cos(1.25d0*yc)
    enddo

    call interpolate_cellcntr2node('ggnb',fc,fvi)

    open(100,file='interp_new.plt')
    write(100,'(a)') 'TITLE ="grid"'
    write(100,'(a)') 'VARIABLES ="X", "Y", "U_intp", "U"'
    write(100,'(a,i0,a,i0,a)') 'ZONE T="VOL_MIXED",N=', num_nodes, ' E=',num_cells, ' ET=TRIANGLE F=FEPOINT'
    do i=1,num_nodes
      write(100,'(4(e19.8,1x))') cell_vert(i,1),cell_vert(i,2), fvi(i),fve(i)
    enddo
    do i=1,num_cells
        write(100,*) (cell2node(i,j),j=1,3)
    enddo
    close(100)

    return
  end subroutine interpolate_check

end module interpolation
