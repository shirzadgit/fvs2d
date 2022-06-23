module gradient_ggnb

  use grid_procs
  use interpolation

  implicit none

  private
  integer,allocatable,save  :: grad_cell2neighbr_num(:,:), grad_neighbr_ptr(:)
  real,allocatable,save     :: grad_ggnb_coef0(:,:), grad_node_weight(:,:,:), grad_neighbhr_weight(:)

  private :: setup
  public  :: grad_ggnb_init, grad_ggnb, grad_ggnb_exp

contains

  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grad_ggnb_init
    implicit  none

    !--------------------------------------------------------------------------!
    ! allocation
    !--------------------------------------------------------------------------!
    allocate(grad_ggnb_coef0(num_cells,2))
    allocate(grad_cell2neighbr_num(num_cells,num_vert_max))
    allocate(grad_node_weight(num_cells,num_vert_max,2))


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

    integer     :: i,j, ic,ic1,ic2, ie,ie1,ie2, in,iv,iv1,iv2, idim,nt
    real        :: af,nxf,nyf, dx,dy, w1,w2,wt1,wt2



    !--------------------------------------------------------------------------!
    ! compute gradient coeff for cell
    !--------------------------------------------------------------------------!
    grad_ggnb_coef0 (:,:)=0.d0
    do ic=1,num_cells
      do ie=1,num_vert_cell(ic)
        ie1=cell2edge(ic,ie)

        af =edge_area(ie1)
        nxf=edge_normal(ie1,1) * edge_normal_sign(ic,ie)
        nyf=edge_normal(ie1,2) * edge_normal_sign(ic,ie)

        iv1=edge2node(ie1,1)
        iv2=edge2node(ie1,2)

        dx = cell_center(ic,1) - cell_vert(iv1,1)
        dy = cell_center(ic,2) - cell_vert(iv1,2)
        w1 = 1.d0/dsqrt(dx**2 + dy**2)

        dx = cell_center(ic,1) - cell_vert(iv2,1)
        dy = cell_center(ic,2) - cell_vert(iv2,2)
        w2 = 1.d0/dsqrt(dx**2 + dy**2)

        wt1 = interp_nodeweight(iv1)
        wt2 = interp_nodeweight(iv2)

        grad_ggnb_coef0(ic,1)= grad_ggnb_coef0(ic,1) + af*nxf/2.d0 * (w1*wt1 + w2*wt2)
        grad_ggnb_coef0(ic,2)= grad_ggnb_coef0(ic,2) + af*nyf/2.d0 * (w1*wt1 + w2*wt2)
      enddo
    enddo


    !--------------------------------------------------------------------------!
    ! compute gradient coeff for neighboring cells
    !--------------------------------------------------------------------------!
    !--------------------------------------------------------------------------!
    ! step-1: compute array dimension for allocation
    !--------------------------------------------------------------------------!
    nt=0;
    do ic=1,num_cells
      do in=1,num_vert_cell(ic)
        iv=cell2node(ic,in)
        do i=node2cell_ptr(iv), node2cell_ptr(iv) + node2cell_ntot(iv) - 1
          ic1=node2cell(i)
          if (ic/=ic1) nt=nt+1
        end do
      end do
    end do
    idim=nt

    !--------------------------------------------------------------------------!
    ! step 2: determine number of neighbourig cells for each node of a cell, excluding the donor cell
    !--------------------------------------------------------------------------!
    grad_cell2neighbr_num(:,:) = 0
    do ic=1,num_cells
      do in=1,num_vert_cell(ic)
        iv=cell2node(ic,in)

        nt=0
        do i=node2cell_ptr(iv), node2cell_ptr(iv) + node2cell_ntot(iv) - 1
          ic1=node2cell(i)
          if (ic/=ic1) nt=nt+1
        end do
        grad_cell2neighbr_num(ic,in)=nt
      end do
    end do

    do ic=1,num_cells
      do in=1,num_vert_cell(ic)
        iv=cell2node(ic,in)
        if (node2cell_ntot(iv) /=(grad_cell2neighbr_num(ic,in)+1)) then
          write(*,*) 'error in grad_cell2neighbr_num!'
          stop 'gradient_cellcntr_setup_ggnb'
        endif
      end do
    end do

    !--------------------------------------------------------------------------!
    ! step 3: computing weighting coeff for each neighbourig cells for each node of a cell, excluding the donor cell
    ! step 4: assign pointer
    !--------------------------------------------------------------------------!
    allocate(grad_neighbhr_weight(idim), grad_neighbr_ptr(idim))
    grad_neighbhr_weight(:)=0.d0
    grad_neighbr_ptr(:)=0
    nt=0;
    do ic=1,num_cells
      do in=1,num_vert_cell(ic)
        iv=cell2node(ic,in)

        do i=node2cell_ptr(iv), node2cell_ptr(iv) + node2cell_ntot(iv) - 1
          ic1=node2cell(i)
          if (ic/=ic1) then
            nt=nt+1
            grad_neighbhr_weight(nt) = interp_nodeweight(iv) * interp_cellweight(i)
            grad_neighbr_ptr(nt) = ic1
          endif
        end do
      end do
    end do
    if (nt/=idim) then
      write(*,*) 'error in grad_neighbhr_weight!'
      stop 'gradient_cellcntr_setup_ggnb'
    endif

    !--------------------------------------------------------------------------!
    ! step 5: computing weighting coeff of each node for a cell
    !--------------------------------------------------------------------------!
    grad_node_weight(:,:,:) = 0.d0
    do ic=1,num_cells
      do in=1,num_vert_cell(ic)
        iv=cell2node(ic,in)

        do ie=1,num_vert_cell(ic)
          ie1=cell2edge(ic,ie)
          do j=1,2
            iv1=edge2node(ie1,j)
            if (iv1==iv) then
              af =edge_area(ie1)
              nxf=edge_normal(ie1,1) * edge_normal_sign(ic,ie)
              nyf=edge_normal(ie1,2) * edge_normal_sign(ic,ie)

              grad_node_weight(ic,in,1)=grad_node_weight(ic,in,1) + af*nxf/2.d0
              grad_node_weight(ic,in,2)=grad_node_weight(ic,in,2) + af*nyf/2.d0
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
  subroutine grad_ggnb (fc,dfc)
    implicit  none

    real,intent(in)   :: fc(num_cells)
    real,intent(out)  :: dfc(num_cells,2)
    integer           :: i,j, ic,ic1, ie,ie1, in,iv,iv1,iv2, nt
    real              :: af,nxf,nyf
    real,allocatable  :: fv(:)

    dfc(:,:)=0.d0

    nt=0;
    do ic=1,num_cells
      dfc(ic,1)= grad_ggnb_coef0(ic,1) * fc(ic)
      dfc(ic,2)= grad_ggnb_coef0(ic,2) * fc(ic)
      do in=1,num_vert_cell(ic)
        do i=1,grad_cell2neighbr_num(ic,in)
          nt=nt+1
          ic1 = grad_neighbr_ptr(nt)
          dfc(ic,1) = dfc(ic,1) + grad_node_weight(ic,in,1)*grad_neighbhr_weight(nt)*fc(ic1)
          dfc(ic,2) = dfc(ic,2) + grad_node_weight(ic,in,2)*grad_neighbhr_weight(nt)*fc(ic1)
        end do
      end do
      dfc(ic,1) = dfc(ic,1)/cell_area(ic)
      dfc(ic,2) = dfc(ic,2)/cell_area(ic)
    end do

    return
  end subroutine grad_ggnb



  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grad_ggnb_exp (fc,dfc)
    implicit  none

    real,intent(in)   :: fc(num_cells)
    real,intent(out)  :: dfc(num_cells,2)
    integer           :: i,j, ic,ic1, ie,ie1, in,iv,iv1,iv2, nt
    real              :: af,nxf,nyf
    real,allocatable  :: fv(:)

    dfc(:,:)=0.d0

    allocate(fv(num_nodes))

    call interpolate_cellcntr2node(fc,fv)

    do ic=1,num_cells
      do ie=1,num_vert_cell(ic)
        ie1=cell2edge(ic,ie)

        af =edge_area(ie1)
        nxf=edge_normal(ie1,1) * edge_normal_sign(ic,ie)
        nyf=edge_normal(ie1,2) * edge_normal_sign(ic,ie)

        iv1=edge2node(ie1,1)
        iv2=edge2node(ie1,2)

        dfc(ic,1) = dfc(ic,1) + nxf*af * (fv(iv1)+fv(iv2))/2.d0
        dfc(ic,2) = dfc(ic,2) + nyf*af * (fv(iv1)+fv(iv2))/2.d0
      enddo
      dfc(ic,1) = dfc(ic,1)/cell_area(ic)
      dfc(ic,2) = dfc(ic,2)/cell_area(ic)
    enddo

    deallocate(fv)

    return
  end subroutine grad_ggnb_exp


end module gradient_ggnb
