module grid_procs

  use mainparam, only : iunit_grid
  use input, only     : file_grid
  use data_type
  use kdtree2_module

  implicit none

  type(kdtree2), pointer,save             :: tree_cellcntr, tree_cellvert

  contains


  !============================================================================!
  !--cell data:
  !  cell(1:ncells)%nvrt    = Number of nodes/vertices of each cell
  !  cell(1:ncells)%node(:) = Pointer to nodes of each cell
  !  cell(1:ncells)%x/y     = cell center x/y-coordinate
  !
  !--node data:
  !  node(1:nnodes)%x     = x-coordinate of the nodes
  !  node(1:nnodes)%y     = y-coordinate of the nodes
  !  node(1:nnodes)%intrp_idw(:)  = intrp coeff based on inverse distance
  !
  !--edge data:
  !  edge(1:nedges)%n1/n2 = edge two end nodes
  !  edge(1:nedges)%c1/c2 = edge left and right cells
  !  edge(1:nedges)%x/y   = edge center x/y-coordinate
  !============================================================================!
  subroutine grid_procs_init
    implicit none

    !--------------------------------------------------------------------------!
    ! read grid file
    !--------------------------------------------------------------------------!
    call grid_read


    !--------------------------------------------------------------------------!
    ! construct grid links and data
    !--------------------------------------------------------------------------!
    call grid_data


    return
  end subroutine grid_procs_init


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grid_read

    implicit none
    integer   :: i,ii,istat
    logical   :: linputfile


    !--------------------------------------------------------------------------!
    ! check grid file
    !--------------------------------------------------------------------------!
    inquire(file=trim(file_grid),exist=linputfile);
    if (.not.linputfile) then
      write(*,*)
      write(*,*) 'linputfile:',linputfile
      write(*,*) 'cannot find "'//trim(file_grid)//'" file!'
      write(*,*) 'error in --> mod:grid_procs, sub:grid_read'
      stop 'program stopped at "grid_read"'
    endif

    !-- open grid file and skip the first line
    open (iunit_grid,file=trim(file_grid),status='unknown',IOSTAT=istat)
    read (iunit_grid,*)
    read (iunit_grid,*) nnodes, ncells_tri, ncells_quad

    !--  allocate node and cell arrays.
    ncells = ncells_tri + ncells_quad
    allocate(node(nnodes))
    allocate(cell(ncells))

    !-- read node/vertices coordinates
    do i=1,nnodes
      read (iunit_grid,*) node(i)%x, node(i)%y
    enddo

    !-- read triangle cells
    do i=1,ncells_tri
      cell(i)%nvrt = 3
      allocate(cell(i)%node(3))
      read (iunit_grid,*) cell(i)%node(1), cell(i)%node(2), cell(i)%node(3)
    enddo

    !-- read quad cells
    do i=1,ncells_quad
      ii = i + ncells_tri
      cell(ii)%nvrt = 4
      allocate(cell(ii)%node(4))
      read (iunit_grid,*) cell(ii)%node(1), cell(ii)%node(2), cell(ii)%node(3), cell(ii)%node(4)
    enddo

    !-- close grid file
    close(iunit_grid)

    return
  end subroutine grid_read


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grid_data

    implicit none
    integer             :: i,j,k,ii,ic,jc,iv,vp
    integer             :: v1,v2,v3,v4, vL,vR, in,im, ie,je
    integer,allocatable :: locedge(:)
    real                :: xc,yc, x1,x2,x3,x4, y1,y2,y3,y4, dx,dy
    real,allocatable    :: tmpr(:,:)
    logical             :: lfound


    !--------------------------------------------------------------------------!
    ! cell centroids
    !--------------------------------------------------------------------------!
    do ic=1,ncells_tri
      xc=0.d0
      yc=0.d0
      do iv=1,3
        vp = cell(ic)%node(iv)
        xc = xc + node(vp)%x
        yc = yc + node(vp)%y
      enddo
      cell(ic)%x = xc/3.d0
      cell(ic)%y = yc/3.d0
    end do

    do i=1,ncells_quad
      ic = ncells_tri + i
      xc=0.d0
      yc=0.d0
      do iv=1,4
        vp = cell(ic)%node(iv)
        xc = xc + node(vp)%x
        yc = yc + node(vp)%y
      enddo
      cell(ic)%x = xc/4.d0
      cell(ic)%y = yc/4.d0
    end do


    !--------------------------------------------------------------------------!
    ! cell volume
    !--------------------------------------------------------------------------!
    do ic=1,ncells_tri
      x1=node(cell(ic)%node(1))%x;  y1=node(cell(ic)%node(1))%y;
      x2=node(cell(ic)%node(2))%x;  y2=node(cell(ic)%node(2))%y;
      x3=node(cell(ic)%node(3))%x;  y3=node(cell(ic)%node(3))%y;
      cell(ic)%vol=tri_area(x1,x2,x3,y1,y2,y3)
    enddo
    do i=1,ncells_quad
      ic=i+ncells_tri
      x1=node(cell(ic)%node(1))%x;  y1=node(cell(ic)%node(1))%y;
      x2=node(cell(ic)%node(2))%x;  y2=node(cell(ic)%node(2))%y;
      x3=node(cell(ic)%node(3))%x;  y3=node(cell(ic)%node(3))%y;
      x4=node(cell(ic)%node(4))%x;  y4=node(cell(ic)%node(4))%y;
      cell(ic)%vol= tri_area(x1,x2,x3,y1,y2,y3) + tri_area(x1,x3,x4,y1,y3,y4)
    enddo


    !--------------------------------------------------------------------------!
    ! build k-d tree based on cell vertices
    !--------------------------------------------------------------------------!
    allocate(tmpr(2,nnodes))
    do i=1,nnodes
      tmpr(1,i)=node(i)%x
      tmpr(2,i)=node(i)%y
    enddo
    tree_cellvert => kdtree2_create(tmpr,sort=.true.,rearrange=.true.)
    deallocate(tmpr)


    !--------------------------------------------------------------------------!
    ! build k-d tree based on cell centers
    !--------------------------------------------------------------------------!
    allocate(tmpr(2,ncells))
    do i=1,ncells
      tmpr(1,i)=cell(i)%x
      tmpr(2,i)=cell(i)%y
    enddo
    tree_cellcntr => kdtree2_create(tmpr,sort=.true.,rearrange=.true.)
    deallocate(tmpr)


    !--------------------------------------------------------------------------!
    ! distribute cell index to nodes
    !--------------------------------------------------------------------------!
    node(1:nnodes)%ncells = 0
    do ic=1,ncells
      do iv=1,cell(ic)%nvrt
        vp = cell(ic)%node(iv)
        node(vp)%ncells = node(vp)%ncells + 1
      enddo
    enddo
    do i=1,nnodes
      allocate(node(i)%cell(node(i)%ncells));
    enddo
    node(1:nnodes)%ncells = 0
    do ic=1,ncells
      do iv=1,cell(ic)%nvrt
        vp = cell(ic)%node(iv)
        node(vp)%ncells = node(vp)%ncells + 1
        node(vp)%cell(node(vp)%ncells) = ic
      enddo
    enddo


    !--------------------------------------------------------------------------!
    ! distribute cell index to nodes
    !--------------------------------------------------------------------------!
    do ic=1,ncells
      cell(ic)%nnghbrs=cell(ic)%nvrt
      allocate(cell(ic)%nghbr(cell(ic)%nvrt));
    enddo

    cells: do ic=1,ncells
      cell_nvrt: do iv=1,cell(ic)%nvrt
        if (iv  < cell(ic)%nvrt) vL = cell(ic)%node(iv+1)
        if (iv == cell(ic)%nvrt) vL = cell(ic)%node(1)
                                 vR = cell(ic)%node(iv)

        !--Loop over the surrounding cells of node vR and find cell neighbor from them
        lfound = .false.
        cells_around_vR : do j = 1, node(vR)%ncells
          jc = node(vR)%cell(j)
          edge_matching : do ii = 1, cell(jc)%nvrt
            if (ii  > 1) v2 = cell(jc)%node(ii-1)
            if (ii == 1) v2 = cell(jc)%node(cell(jc)%nvrt)
                         v1 = cell(jc)%node(ii)

            if (v1==vR .and. v2==vL) then
              lfound = .true.
              im = ii+1
              if (im > cell(jc)%nvrt) im = im - cell(jc)%nvrt
              exit edge_matching
            endif
          end do edge_matching

          if (lfound) exit cells_around_vR
        end do cells_around_vR

        in = iv + 2
        if (in > cell(ic)%nvrt) in = in - cell(ic)%nvrt

        if (lfound) then
          cell(ic)%nghbr(in) = jc
          cell(jc)%nghbr(im) = ic
        else
          cell(ic)%nghbr(in) = 0
        endif
      end do cell_nvrt
    end do cells


    !--------------------------------------------------------------------------!
    ! count number of edges and allocate edge array
    !--------------------------------------------------------------------------!
    nedges=0
    do ic = 1, ncells
      do j=1,cell(ic)%nnghbrs
        if ( cell(ic)%nghbr(j) > ic .or. cell(ic)%nghbr(j)==0 ) nedges = nedges + 1
      enddo
    enddo
    allocate(edge(nedges))


    !--------------------------------------------------------------------------!
    ! assign edge two end nodes (n1/n2), and left and right cells (c1/c2)
    !--------------------------------------------------------------------------!
    nedges=0
    edge(:)%c1 = 0
    edge(:)%c2 = 0
    !--triangle cells
    do ic = 1, ncells_tri
      v1=cell(ic)%node(1)
      v2=cell(ic)%node(2)
      v3=cell(ic)%node(3)

      if ( cell(ic)%nghbr(3) > ic  .or. cell(ic)%nghbr(3)==0 ) then
        nedges = nedges + 1
        edge(nedges)%n1 = v1
        edge(nedges)%n2 = v2
        edge(nedges)%c1 = ic
        edge(nedges)%c2 = cell(ic)%nghbr(3)
      endif

      if ( cell(ic)%nghbr(1) > ic .or. cell(ic)%nghbr(1)==0 ) then
        nedges = nedges + 1
        edge(nedges)%n1 = v2
        edge(nedges)%n2 = v3
        edge(nedges)%c1 = ic
        edge(nedges)%c2 = cell(ic)%nghbr(1)
      endif

      if ( cell(ic)%nghbr(2) > ic .or. cell(ic)%nghbr(2)==0 ) then
        nedges = nedges + 1
        edge(nedges)%n1 = v3
        edge(nedges)%n2 = v1
        edge(nedges)%c1 = ic
        edge(nedges)%c2 = cell(ic)%nghbr(2)
      endif
    enddo

    !--quad cells
    do i = 1, ncells_quad
      ic=ncells_tri + i
      v1=cell(ic)%node(1)
      v2=cell(ic)%node(2)
      v3=cell(ic)%node(3)
      v4=cell(ic)%node(4)

      if ( cell(ic)%nghbr(3) > ic  .or. cell(ic)%nghbr(3)==0 ) then
        nedges = nedges + 1
        edge(nedges)%n1 = v1
        edge(nedges)%n2 = v2
        edge(nedges)%c1 = ic
        edge(nedges)%c2 = cell(ic)%nghbr(3)
      endif

      if ( cell(ic)%nghbr(4) > ic  .or. cell(ic)%nghbr(4)==0 ) then
        nedges = nedges + 1
        edge(nedges)%n1 = v2
        edge(nedges)%n2 = v3
        edge(nedges)%c1 = ic
        edge(nedges)%c2 = cell(ic)%nghbr(4)
      endif

      if ( cell(ic)%nghbr(1) > ic .or. cell(ic)%nghbr(1)==0 ) then
        nedges = nedges + 1
        edge(nedges)%n1 = v3
        edge(nedges)%n2 = v4
        edge(nedges)%c1 = ic
        edge(nedges)%c2 = cell(ic)%nghbr(1)
      endif

      if ( cell(ic)%nghbr(2) > ic .or. cell(ic)%nghbr(2)==0 ) then
        nedges = nedges + 1
        edge(nedges)%n1 = v4
        edge(nedges)%n2 = v1
        edge(nedges)%c1 = ic
        edge(nedges)%c2 = cell(ic)%nghbr(2)
      endif
    enddo


    !--------------------------------------------------------------------------!
    ! assign cell index to edge
    !--------------------------------------------------------------------------!
    do ic=1,ncells
      allocate(cell(ic)%edge(cell(ic)%nvrt))
    enddo

    nedges=0
    !--triangle cells
    allocate(locedge(3));  locedge(1)=2;  locedge(2)=3;  locedge(3)=1;
    do ic = 1, ncells_tri
      if ( cell(ic)%nghbr(3) > ic  .or. cell(ic)%nghbr(3)==0 ) then
        nedges = nedges + 1
        cell(ic)%edge(1)= nedges
      else
        jc=cell(ic)%nghbr(3)
        tlp1: do k=1,3
          if (cell(jc)%nghbr(k)==ic) then
            ie=locedge(k)
            cell(ic)%edge(1)=cell(jc)%edge(ie)
            exit tlp1
          endif
        enddo tlp1
      endif

      if ( cell(ic)%nghbr(1) > ic .or. cell(ic)%nghbr(1)==0 ) then
       nedges = nedges + 1
       cell(ic)%edge(2)= nedges
     else
       jc=cell(ic)%nghbr(1)
       tlp2: do k=1,3
         if (cell(jc)%nghbr(k)==ic) then
           ie=locedge(k)
           cell(ic)%edge(2)=cell(jc)%edge(ie)
           exit tlp2
         endif
       enddo tlp2
      endif

      if ( cell(ic)%nghbr(2) > ic .or. cell(ic)%nghbr(2)==0 ) then
       nedges = nedges + 1
       cell(ic)%edge(3)= nedges
     else
       jc=cell(ic)%nghbr(2)
       tlp3: do k=1,3
         if (cell(jc)%nghbr(k)==ic) then
           ie=locedge(k)
           cell(ic)%edge(3)=cell(jc)%edge(ie)
           exit tlp3
         endif
       enddo tlp3
      endif
    enddo

    !--quad cells
    deallocate(locedge)
    allocate(locedge(4));  locedge(1)=3;  locedge(2)=4;  locedge(3)=1;  locedge(4)=2
    do i = 1, ncells_quad
      ic=ncells_tri + i
      if ( cell(ic)%nghbr(3) > ic  .or. cell(ic)%nghbr(3)==0 ) then
        nedges = nedges + 1
        cell(ic)%edge(1)= nedges
      else
        jc=cell(ic)%nghbr(3)
        qlp1: do k=1,4
          if (cell(jc)%nghbr(k)==ic) then
            ie=locedge(k)
            cell(ic)%edge(1)=cell(jc)%edge(ie)
            exit qlp1
          endif
        enddo qlp1
      endif

      if ( cell(ic)%nghbr(4) > ic  .or. cell(ic)%nghbr(4)==0 ) then
        nedges = nedges + 1
        cell(ic)%edge(2)= nedges
      else
        jc=cell(ic)%nghbr(4)
        qlp2: do k=1,4
          if (cell(jc)%nghbr(k)==ic) then
            ie=locedge(k)
            cell(ic)%edge(2)=cell(jc)%edge(ie)
            exit qlp2
          endif
        enddo qlp2
      endif

      if ( cell(ic)%nghbr(1) > ic .or. cell(ic)%nghbr(1)==0 ) then
        nedges = nedges + 1
        cell(ic)%edge(3)= nedges
      else
        jc=cell(ic)%nghbr(1)
        qlp3: do k=1,4
          if (cell(jc)%nghbr(k)==ic) then
            ie=locedge(k)
            cell(ic)%edge(3)=cell(jc)%edge(ie)
            exit qlp3
          endif
        enddo qlp3
      endif

      if ( cell(ic)%nghbr(2) > ic .or. cell(ic)%nghbr(2)==0 ) then
        nedges = nedges + 1
        cell(ic)%edge(4)= nedges
      else
        jc=cell(ic)%nghbr(2)
        qlp4: do k=1,4
          if (cell(jc)%nghbr(k)==ic) then
            ie=locedge(k)
            cell(ic)%edge(4)=cell(jc)%edge(ie)
            exit qlp4
          endif
        enddo qlp4
      endif
    enddo


    !--------------------------------------------------------------------------!
    ! compute edge center, area, normal, and tangential vectors
    !--------------------------------------------------------------------------!
    do i=1,nedges
      v1=edge(i)%n1
      v2=edge(i)%n2

      dx=node(v2)%x - node(v1)%x
      dy=node(v2)%y - node(v1)%y

      edge(i)%area = dsqrt(dx**2+dy**2)

      edge(i)%x = 0.5d0*(node(v1)%x + node(v2)%x)
      edge(i)%y = 0.5d0*(node(v1)%y + node(v2)%y)

      edge(i)%nx = dy/edge(i)%area
      edge(i)%ny =-dx/edge(i)%area

      edge(i)%tx = dx/edge(i)%area
      edge(i)%ty = dy/edge(i)%area
    enddo


    !--------------------------------------------------------------------------!
    ! compute edge normal direction for each cell
    !--------------------------------------------------------------------------!
    do ic=1,ncells_tri
      allocate(cell(ic)%nrmlsign(3))
      cell(ic)%nrmlsign(1:3)=1.d0
      do ie=1,cell(ic)%nvrt
        je=cell(ic)%edge(ie)
        jc=edge(je)%c1
        if (ic/=jc) cell(ic)%nrmlsign(ie)=-1.d0
      enddo
    enddo

    do i=1,ncells_quad
      ic=i+ncells_tri
      allocate(cell(ic)%nrmlsign(4))
      cell(ic)%nrmlsign(1:4)=1.d0
      do ie=1,cell(ic)%nvrt
        je=cell(ic)%edge(ie)
        jc=edge(je)%c1
        if (ic/=jc) cell(ic)%nrmlsign(ie)=-1.d0
      enddo
    enddo

    return
  end subroutine grid_data


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  function tri_area(x1,x2,x3,y1,y2,y3) result(area)
    implicit none
    real, intent(in) :: x1,x2,x3,y1,y2,y3
    real :: area

    area = 0.5d0*( x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2) )

    return
  end function tri_area


end module grid_procs
