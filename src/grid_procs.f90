module grid_procs

  use mainparam, only : iunit_grid
  use input, only     : file_grid
  use grid_graph
  use kdtree2_module

  implicit none

  integer,parameter                       :: num_vert_max=3
  integer,save                            :: num_cells, num_nodes, num_edges, num_edges_bndry, num_cells_bndry, num_nodes_bndry
  real,allocatable,dimension(:,:),save    :: cell_vert, cell_center, cell_dist2vert, cell_dist2edgecenter
  real,allocatable,dimension(:,:),save    :: edge_center, edge_normal, edge_tangnt, edge_normal_sign
  real,allocatable,dimension(:),save      :: edge_area, cell_area
  integer,allocatable,dimension(:,:),save :: cell2cell, cell2edge, edge2cell, edge2node, cell2node, edgemap
  integer,allocatable,dimension(:),save   :: num_edge_cell, num_vert_edge, num_vert_cell, num_cells_neighbr
  integer,allocatable,dimension(:),save   :: node2nodeptr, node2node, node2edge, node2cell, node2cellptr
  integer,allocatable,dimension(:),save   :: cells_bndry_ptr, edges_bndry_ptr, nodes_bndry_ptr
  type(kdtree2), pointer,save             :: tree_cellcntr, tree_cellvert

contains


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grid_procs_init
    implicit none

    !--------------------------------------------------------------------------!
    ! read grid file
    !--------------------------------------------------------------------------!
    call grid_read


    !--------------------------------------------------------------------------!
    ! prepare grid links
    !--------------------------------------------------------------------------!
    call grid_links


    !--------------------------------------------------------------------------!
    ! compute grid properties (erea, normal, etc.)
    !--------------------------------------------------------------------------!
    call grid_props


    !--------------------------------------------------------------------------!
    ! check grid and computed properties
    !--------------------------------------------------------------------------!
    call grid_check


    return
  end subroutine grid_procs_init


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grid_read
    implicit none

    integer             :: i,istat,nloc1,nloc2,eloc1,eloc2
    real,allocatable    :: rwrk2d(:,:)
    character           :: dchar*400

    !--------------------------------------------------------------------------!
    ! open grid file
    !--------------------------------------------------------------------------!
    open (iunit_grid,file=trim(file_grid),status='unknown',IOSTAT=istat)
    read (iunit_grid,*)
    read (iunit_grid,*)
    read (iunit_grid,'(a)') dchar


    !--------------------------------------------------------------------------!
    ! determine number of nodes and cells
    !--------------------------------------------------------------------------!
    do i=1,len(dchar)-2
      if (dchar(i:i+1)=='N=' .or. dchar(i:i+1)=='n=')		nloc1=i+2
      if (dchar(i:i+1)=='E=' .or. dchar(i:i+1)=='e=')		nloc2=i-1
      if (dchar(i:i+1)=='E=' .or. dchar(i:i+1)=='e=')		eloc1=i+2
      if (dchar(i:i+2)=='ET=' .or. dchar(i:i+2)=='et=')	eloc2=i-1
    enddo
    read (dchar(nloc1:nloc2),*,iostat=istat) num_nodes
    read (dchar(eloc1:eloc2),*,iostat=istat) num_cells


    !--------------------------------------------------------------------------!
    ! allocate and read nodes/vertices coordinates
    !--------------------------------------------------------------------------!
    allocate(cell_vert(num_nodes,2))
    do i=1,num_nodes
      read (iunit_grid,*) cell_vert(i,1),cell_vert(i,2)
    enddo


    !--------------------------------------------------------------------------!
    ! allocate and read cell connectivity
    !--------------------------------------------------------------------------!
    allocate(cell2node(num_cells,num_vert_max))
    do i=1,num_cells
      read (iunit_grid,*) cell2node(i,1),cell2node(i,2),cell2node(i,3)
    enddo


    !--------------------------------------------------------------------------!
    ! allocate and read number of edges/vertices per cell
    !--------------------------------------------------------------------------!
    allocate(num_vert_cell(num_cells))
    num_vert_cell(:)=num_vert_max


    !--------------------------------------------------------------------------!
    ! initialize k-d tree based on cell nodes/vertices
    !--------------------------------------------------------------------------!
    allocate(rwrk2d(2,num_nodes))
    do i=1,num_nodes
      rwrk2d(1:2,i)=cell_vert(i,1:2)
    enddo
    tree_cellvert => kdtree2_create(rwrk2d,sort=.true.,rearrange=.true.)
    deallocate(rwrk2d)


    close(iunit_grid)

    return
  end subroutine grid_read


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grid_links
    implicit none

    integer               :: ndim,i,j
    integer,allocatable   :: iwrk1d(:),iwrk2d(:,:)
    integer,allocatable   :: cell2cell_tmp(:,:), cell2edge_tmp(:,:), cell2node_tmp(:,:)

    !--------------------------------------------------------------------------!
    ! allocation
    !--------------------------------------------------------------------------!
    allocate(cell2cell(num_cells,num_vert_max), cell2edge(num_cells,num_vert_max))

    !--------------------------------------------------------------------------!
    ! find links for edge2node & edge2cell
    !--------------------------------------------------------------------------!
    allocate(cell2node_tmp(num_vert_max,num_cells))
    do i=1,num_cells
      cell2node_tmp(1:num_vert_max,i) = cell2node(i,1:num_vert_max)
    end do

    num_edges=0
    allocate(iwrk2d(6,4*num_cells))
    call findEdges(num_vert_max,num_cells,num_edges,cell2node_tmp,iwrk2d)

    allocate(edgemap(6,num_edges),edge2node(num_edges,2),edge2cell(num_edges,2));
    edgemap(1:6,1:num_edges)=iwrk2d(1:6,1:num_edges)

    do i=1,num_edges
      edge2node(i,1:2)=edgemap(1:2,i)
      edge2cell(i,1:2)=edgemap(3:4,i)
    end do

    !--------------------------------------------------------------------------!
    ! find links for cell2cell and cell2edge
    !--------------------------------------------------------------------------!
    allocate(cell2cell_tmp(num_vert_max,num_cells))
    allocate(cell2edge_tmp(num_vert_max,num_cells))

    call findcell2edgemap(edgemap,cell2cell_tmp,cell2edge_tmp,num_vert_max,num_edges,num_cells)

    do i=1,num_cells
      cell2cell(i,1:num_vert_max)=cell2cell_tmp(1:num_vert_max,i)
      cell2edge(i,1:num_vert_max)=cell2edge_tmp(1:num_vert_max,i)
    end do


    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------
    !allocate(node2nodeptr(num_nodes+1))
    !allocate(node2node(2*num_edges))
    !allocate(node2edge(2*num_edges))
    !call findnodemap(edgemap,lnode2node,lnode2edge,lnode2nodeptr,nnodes,nedges)


    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------
    !allocate(lnode2cellptr(num_nodes+1))
    !call findnode2cellmap(cell2node,iwrk1d,node2cellptr,nvermax,num_nodes,num_cells)
    !ndim=node2cellptr(nnodes+1)-1
    !allocate(node2cell(ndim))
    !node2cell(1:ndim)=iwrk1d(1:ndim)

    !deallocate(iwrk1d, iwrk2d)
    deallocate(cell2cell_tmp, cell2edge_tmp, cell2node_tmp)

    return
  end subroutine grid_links


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grid_props
    implicit none

    integer             :: ic,ic1,ie,ie1,iv,iv1,iv2,in,in1,in2,nt,i,j,dum1,dum2
    integer,allocatable :: iwrk1d(:)
    real                :: tx,ty,xc,yc,dx,dy, vol,xf,af,nx
    real,allocatable    :: rwrk2d(:,:)

    !--------------------------------------------------------------------------!
    ! allocation
    !--------------------------------------------------------------------------!
    allocate(cell_center(num_cells,2), cell_area(num_cells), cell_dist2vert(num_cells,num_vert_max))
    allocate(edge_center(num_edges,2), edge_area(num_edges))
    allocate(edge_normal(num_edges,2), edge_tangnt(num_edges,2), edge_normal_sign(num_cells,num_vert_max) )


    !--------------------------------------------------------------------------!
    ! determine boundary edges
    !--------------------------------------------------------------------------!
    allocate(iwrk1d(num_edges))
    nt=0;
    do ie=1,num_edges
      if (edge2cell(ie,1)*edge2cell(ie,2)==0) then
        nt=nt+1
        if (nt>0) iwrk1d(nt)=ie
      endif
    end do
    num_edges_bndry=nt
    if (nt>0) then
      allocate(edges_bndry_ptr(nt))
      edges_bndry_ptr(1:nt)= iwrk1d(1:nt)
    endif
    deallocate(iwrk1d)


    !--------------------------------------------------------------------------!
    ! determine boundary nodes
    !--------------------------------------------------------------------------!
    allocate(iwrk1d(num_edges*2))
    nt=0;
    do ie=1,num_edges
      if (edge2cell(ie,1)*edge2cell(ie,2)==0) then
        iv1=edge2node(ie,1)
        iv2=edge2node(ie,2)

        if (nt==0) then
          nt=nt+1
          iwrk1d(nt)=iv1; nt=nt+1
          iwrk1d(nt)=iv2
        else
          dum1=1
          dum2=1
          do j=1,nt;
            dum1=min(dum1,abs(iwrk1d(j)-iv1));
            dum2=min(dum2,abs(iwrk1d(j)-iv2));
          enddo
          if (dum1/=0) then
            nt=nt+1
            iwrk1d(nt)=iv1
          endif
          if (dum2/=0) then
            nt=nt+1
            iwrk1d(nt)=iv2
          endif
        endif
      endif
    end do
    num_nodes_bndry=nt
    if (nt>0) then
      allocate(nodes_bndry_ptr(nt))
      nodes_bndry_ptr(1:nt)= iwrk1d(1:nt)
    endif
    deallocate(iwrk1d)


    !--------------------------------------------------------------------------!
    ! determine boundary cells
    !--------------------------------------------------------------------------!
    allocate(iwrk1d(num_cells))
    nt=0;
    do ic=1,num_cells
      if (cell2cell(ic,1)*cell2cell(ic,2)*cell2cell(ic,3)==0) then
        nt=nt+1
        if (nt>0) iwrk1d(nt)=ic;
      endif
    end do
    num_cells_bndry=nt
    if (nt>0) then
      allocate(cells_bndry_ptr(nt))
      cells_bndry_ptr(1:nt)= iwrk1d(1:nt)
    endif
    deallocate(iwrk1d)


    !--------------------------------------------------------------------------!
    ! compute cell center (works for triangle for now)
    !--------------------------------------------------------------------------!
    do ic=1,num_cells
      xc=0.d0;
      yc=0.d0;
      do iv=1,num_vert_cell(ic)
        xc = xc + cell_vert(cell2node(ic,iv),1)
        yc = yc + cell_vert(cell2node(ic,iv),2)
      enddo
      cell_center(ic,1)=xc/dble(num_vert_cell(ic))
      cell_center(ic,2)=yc/dble(num_vert_cell(ic))
    enddo


    !--------------------------------------------------------------------------!
    ! compute edge center and area
    !--------------------------------------------------------------------------!
    do ie=1,num_edges
      iv1=edge2node(ie,1)
      iv2=edge2node(ie,2)
      edge_center(ie,1)=( cell_vert(iv1,1) + cell_vert(iv2,1) ) / 2.d0
      edge_center(ie,2)=( cell_vert(iv1,2) + cell_vert(iv2,2) ) / 2.d0

      dx=cell_vert(iv1,1) - cell_vert(iv2,1)
      dy=cell_vert(iv1,2) - cell_vert(iv2,2)
      edge_area(ie)=dsqrt(dx**2 + dy**2)
    enddo


    !--------------------------------------------------------------------------!
    ! compute edge normal
    !--------------------------------------------------------------------------!
    do ie=1,num_edges
      iv1=edge2node(ie,1)
      iv2=edge2node(ie,2)
      dx=cell_vert(iv2,1) - cell_vert(iv1,1)
      dy=cell_vert(iv2,2) - cell_vert(iv1,2)
      tx=dx/edge_area(ie)
      ty=dy/edge_area(ie)

      edge_normal(ie,1)= ty
      edge_normal(ie,2)=-tx
      edge_tangnt(ie,1)= tx
      edge_tangnt(ie,2)= ty
    end do


    !--------------------------------------------------------------------------!
    ! compute edge normal direction for each cell
    !--------------------------------------------------------------------------!
    do ic=1,num_cells
      do ie=1,num_vert_cell(ic)
        ie1=cell2edge(ic,ie)
        ic1=edge2cell(ie1,1)
        if (ic==ic1) then
          edge_normal_sign(ic,ie)=+1.d0
        else
          edge_normal_sign(ic,ie)=-1.d0
        endif
      enddo
    enddo


    !--------------------------------------------------------------------------!
    ! compute cell volume
    !--------------------------------------------------------------------------!
    do ic=1,num_cells
      vol=0.d0
      do iv=1,3
        ie=cell2edge(ic,iv)
        nx= edge_normal(ie,1) * edge_normal_sign(ic,iv)
        xf= edge_center(ie,1)
        af= edge_area(ie)
        vol = vol +  nx*xf*af
      enddo
        cell_area(ic) = vol
    enddo


    !--------------------------------------------------------------------------!
    ! compute distance between cell center and cell vertices
    !--------------------------------------------------------------------------!
    do ic=1,num_cells
      do iv=1,num_vert_cell(ic)
        in=cell2node(ic,iv)
        dx = cell_vert(in,1) - cell_center(ic,1)
        dy = cell_vert(in,2) - cell_center(ic,2)
        cell_dist2vert(ic,iv)=dsqrt (dx**2 + dy**2)
      enddo
    enddo


    !--------------------------------------------------------------------------!
    ! compute distance between cell center and edge center
    !--------------------------------------------------------------------------!
    do ic=1,num_cells
      do iv=1,num_vert_cell(ic)
        in=cell2node(ic,iv)
        !dx = cell_vert(in,1) - cell_center(ic,1)
        !dy = cell_vert(in,2) - cell_center(ic,2)
        !cell_dist2edgecenter(ic,iv)=dsqrt (dx**2 + dy**2)
      enddo
    enddo


    !--------------------------------------------------------------------------!
    ! determine number of surrounding neighbor cells of each node
    !--------------------------------------------------------------------------!
    allocate(num_cells_neighbr(num_nodes))
    num_cells_neighbr(:) = 0
    do ic=1,num_cells
      do iv=1,num_vert_cell(ic)
        iv1=cell2node(ic,iv)
        num_cells_neighbr(iv1) = num_cells_neighbr(iv1) + 1
      enddo
    enddo


    !--------------------------------------------------------------------------!
    ! initialize k-d tree based on cell centers
    !--------------------------------------------------------------------------!
    allocate(rwrk2d(2,num_cells))
    do i=1,num_cells
      rwrk2d(1:2,i)=cell_center(i,1:2)
    enddo
    tree_cellcntr => kdtree2_create(rwrk2d,sort=.true.,rearrange=.true.)
    deallocate(rwrk2d)



    return
  end subroutine grid_props


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grid_check
    implicit none
    integer             :: ic,in,in1,in2,ie,ie1,ie2,i,j,n1,nt,iv,iv1,iv2, num_nearests
    integer,allocatable :: iwrk1d(:)
    real                :: vol, r(2)
    type(kdtree2_result),allocatable	:: fnearest(:)

    write(*,*)
    write(*,'(a,i0)') 'number of nodes: ',num_nodes
    write(*,'(a,i0)') 'number of cells: ',num_cells
    write(*,'(a,i0)') 'number of edges: ',num_edges
    write(*,'(a,i0)') 'number of bondary nodes: ',num_nodes_bndry
    write(*,'(a,i0)') 'number of bondary edges: ',num_edges_bndry
    write(*,'(a,i0)') 'number of bondary cells: ',num_cells_bndry
    write(*,'(a,en20.11)') 'min cell volume: ',minval(cell_area)
    write(*,'(a,en20.11)') 'max cell volume: ',maxval(cell_area)
    write(*,*)

      !write(*,*)
      !write(*,*)
      !do n=1,5
      !  write(*,*) node2nodeptr(n)
      !  write(*,*) xv(lnode2edge(n)),yv(lnode2edge(n))
      !end

      !n=10
      !write(*,*) lnode2cell(n),lnode2cellptr(n)
      !write(*,*) lcell2node(1,lnode2cellptr(n)),lcell2node(2,lnode2cellptr(n)),lcell2node(3,lnode2cellptr(n))
      !ic=20
      !write(*,*) cell2cell(ic,1), cell2cell(ic,2), cell2cell(ic,3)
      !write(*,*) cell2edge(ic,1), cell2edge(ic,2), cell2edge(ic,3)
      !write(*,*) edge2node(cell2edge(ic,1),1), edge2node(cell2edge(ic,1),2)
      !write(*,*) edge2node(cell2edge(ic,2),1), edge2node(cell2edge(ic,2),2)
      !write(*,*) edge2node(cell2edge(ic,3),1), edge2node(cell2edge(ic,3),2)
      !write(*,*) edge_normal_sign(ic,1), edge_normal_sign(ic,2), edge_normal_sign(ic,3)

      !ie=35;  write(*,*) ie, edge2cell(ie,1), edge2cell(ie,2)
      !ie=4;   write(*,*) ie, edge2cell(ie,1), edge2cell(ie,2)
      !ie=47;  write(*,*) ie, edge2cell(ie,1), edge2cell(ie,2)

    do i=1,num_edges_bndry
        !ie=edges_bndry_ptr(i)
        !write(*,*) cell_vert(edge2node(ie,1),1),cell_vert(edge2node(ie,1),2)
        !write(*,*) cell_vert(edge2node(ie,2),1),cell_vert(edge2node(ie,2),2)
    end do

    do i=1,num_nodes_bndry
      !in=nodes_bndry_ptr(i)
      !write(*,*) cell_vert(in,1),cell_vert(in,2)
    end do

    do i=1,num_edges_bndry
      ie=edges_bndry_ptr(i)
      if (edge2cell(ie,1)>0) then
        ic=edge2cell(ie,1)
      elseif (edge2cell(ie,2)>0) then
        ic=edge2cell(ie,2)
      else
        write(*,*) 'error'
      endif
      !if (ie==cell2edge(ic,1))

      !write(*,*) edge_normal(ie,1),edge_normal(ie,2)
    enddo

    do i=1,1 !num_cells_bndry
      ic=cells_bndry_ptr(i)
      !write(*,*) edge_normal_sign(ic,1),edge_normal_sign(ic,2),edge_normal_sign(ic,3)
      !write(*,*) 'variables="X" "Y"'
      do j=1,3
      !  write(*,*) cell_vert(cell2node(ic,j),1), cell_vert(cell2node(ic,j),2)
      enddo
      j=1
      !write(*,*) cell_vert(cell2node(ic,j),1), cell_vert(cell2node(ic,j),2)
    enddo

    num_nearests=5
    allocate(fnearest(num_nearests))
    do ic=1,1 !num_cells
      do in=1,num_vert_cell(ic)
        iv=cell2node(ic,in)
        write(*,*) cell_vert(iv,1),cell_vert(iv,2)
        r(1:2)=cell_vert(iv,1:2)
        call kdtree2_n_nearest(tp=tree_cellcntr,qv=r,nn=num_nearests,results=fnearest)
        do i=1,num_nearests
          iv1=fnearest(i)%idx
          write(*,*) cell_center(iv1,1),cell_center(iv1,2)
        enddo
      enddo
    enddo


    do i=1,num_nodes_bndry
      in=nodes_bndry_ptr(i)
      !write(*,*) num_cells_neighbr(in)
    enddo



    open(100,file='grid_vol.dat')
    do ic=1,num_cells
      write(100,*) ic,cell_area(ic)
    enddo
    close(100)

    open(100,file='grid_xy.dat')
    do i=1,num_nodes
      write(100,*) cell_vert(i,1),cell_vert(i,2)
    enddo
    close(100)

    open(100,file='grid_con.dat')
    do i=1,num_cells
        write(100,*) (cell2node(i,j),j=1,3)
    enddo
    close(100)



    return
  end subroutine grid_check


end module grid_procs
