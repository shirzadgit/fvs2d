module gradient

  use grid_procs
  use interpolation

  implicit none

  integer,allocatable,save  :: grad_cell2neighbr_num(:,:), grad_neighbr_ptr(:)
  real,allocatable,save     :: grad_ggnb_coef0(:,:), grad_node_weight(:,:,:), grad_neighbhr_weight(:)

  integer,allocatable,save  :: grad_ggcb_ptr(:,:)
  real,allocatable,save     :: grad_ggcb_coef0(:,:), grad_ggcb_coefnb(:,:,:)

contains

  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine gradient_init
    implicit  none

    !--------------------------------------------------------------------------!
    ! setup gradient opertor at cell centers: Green-Gauss cell-based (ggcb)
    !--------------------------------------------------------------------------!
    call gradient_cellcntr_setup_ggcb


    !--------------------------------------------------------------------------!
    ! setup gradient opertor at cell centers: Green-Gauss node-based (ggnb)
    !--------------------------------------------------------------------------!
    call gradient_cellcntr_setup_ggnb


    return
  end subroutine gradient_init


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine gradient_cellcntr_setup_ggcb
    implicit  none

    integer     :: i,j, ic,ic1,ic2, ie,ie1,ie2, in,iv,iv1,iv2
    real        :: af,nxf,nyf, xf,yf, xc,yc, xc1,yc1, dx,dy, d0,d1

    allocate( grad_ggcb_ptr(num_cells,num_vert_max), grad_ggcb_coef0(num_cells,2), grad_ggcb_coefnb(num_cells,num_vert_max,2))

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
  end subroutine gradient_cellcntr_setup_ggcb

  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine gradient_cellcntr_setup_ggnb
    implicit  none

    integer     :: i,j, ic,ic1,ic2, ie,ie1,ie2, in,iv,iv1,iv2, idim,nt
    real        :: af,nxf,nyf, dx,dy, w1,w2,wt1,wt2

    i=size(node2cell,dim=1)
    allocate(grad_ggnb_coef0(num_cells,2))
    allocate(grad_cell2neighbr_num(num_cells,num_vert_max))
    allocate(grad_node_weight(num_cells,num_vert_max,2))

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
  end subroutine gradient_cellcntr_setup_ggnb



  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine gradient_cellcntr (grad_type,fc,dfc)
    implicit  none

    real,intent(in)   :: fc(num_cells)
    real,intent(out)  :: dfc(num_cells,2)
    character(len=*)  :: grad_type
    integer           :: i,j, ic,ic1, ie,ie1, in,iv,iv1,iv2, nt
    real              :: af,nxf,nyf
    real,allocatable  :: fv(:)

    dfc(:,:)=0.d0
    if (grad_type=='ggnb') then
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

    elseif (grad_type=='ggnb_explicit') then
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

    elseif (grad_type=='ggcb') then
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

    else
      write(*,*) 'error: Green-Gauss node-based and cell-based are only implemented!'
      stop 'error in gradient_cellcntr'
    endif


    return
  end subroutine gradient_cellcntr

end module gradient
