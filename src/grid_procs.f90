module grid_procs

  use mainparam
  use input,      only  : file_grid, file_bc
  use data_grid
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
    ! read grid file
    !--------------------------------------------------------------------------!
    call grid_bc_read


    !--------------------------------------------------------------------------!
    ! construct grid links and data
    !--------------------------------------------------------------------------!
    call grid_data


    !--------------------------------------------------------------------------!
    ! construct grid links and data
    !--------------------------------------------------------------------------!
    call grid_data_verify


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
      write(*,'(a,a,a)') 'cannot find ',trim(file_grid),' file!'
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
  subroutine grid_bc_read

    implicit none
    integer   :: i,ib,istat
    logical   :: linputfile


    !--------------------------------------------------------------------------!
    ! check grid file
    !--------------------------------------------------------------------------!
    inquire(file=trim(file_bc),exist=linputfile);
    if (.not.linputfile) then
      write(*,*)
      write(*,*) 'linputfile:',linputfile
      write(*,'(a,a,a)') 'cannot find ',trim(file_bc),' file!'
      write(*,*) 'error in --> mod:grid_procs, sub:grid_bc_read'
      stop 'program stopped at "grid_bc_read"'
    endif


    !-- open bc file and skip the first line
    open (iunit_bc,file=trim(file_bc),status='unknown',IOSTAT=istat)

    read (iunit_bc,*) nbndries

    allocate(bndry(nbndries))

    do ib=1,nbndries
      read (iunit_bc,*) bndry(ib)%ncells, bndry(ib)%type
    enddo

    do ib=1,nbndries
      allocate( bndry(ib)%cell(bndry(ib)%ncells) )
      do i=1,bndry(ib)%ncells
        read (iunit_bc,*) bndry(ib)%cell(i)
      enddo
    enddo

    close(iunit_bc)

    return
  end subroutine grid_bc_read


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grid_data

    implicit none
    integer             :: i,j,k,ii,ic,jc,iv,vp, ib
    integer             :: v1,v2,v3,v4, vL,vR, in,im, ie,je
    integer             :: idum4(4), tri_v2e(3), qud_v2e(4)
    integer,allocatable :: locedge(:), idum1(:), idum2(:), idum3(:)
    real                :: xc,yc, x1,x2,x3,x4, y1,y2,y3,y4, dx,dy, xf,yf, at
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
    ! compute effective distance, error
    !--------------------------------------------------------------------------!
    at = 0.d0
    do ic=1,ncells
      at = at + cell(ic)%vol
    enddo
    heff1 = sqrt(at/dble(ncells))

    at = 0.d0
    do ic=1,ncells
      at = at + sqrt(cell(ic)%vol)
    enddo
    heff2 = at/dble(ncells)


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
    ! count number pf cells around a node/vertice and allocate node2cell
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


    !--------------------------------------------------------------------------!
    ! build node to cell link
    !--------------------------------------------------------------------------!
    node(1:nnodes)%ncells = 0
    do ic=1,ncells
      do iv=1,cell(ic)%nvrt
        vp = cell(ic)%node(iv)
        node(vp)%ncells = node(vp)%ncells + 1
        node(vp)%cell(node(vp)%ncells) = ic
      enddo
    enddo


    !--------------------------------------------------------------------------!
    ! Find edge-neighbor of each cell
    !            o                               o
    !           . \                            .   \
    !          .   \                          .      \
    !         .  2  \                        .    7   \
    !     v1 .       \ v3                v3 .           \ v2
    !       o---------o                   o--------------o---o
    !     .  \       / \                 / \            /     \
    !    .    \  ic /   \               /   \    ic    /       \
    !   .  5   \   /  8  \             /  9  \        /    12   \
    !  .        \ /       \           /       \      /           \
    ! o----------o---------o         o---------o====o-------------o
    !            v2                           v4    v1
    !..........................................................................
    ! cell(ic)%nghbr = 3              cell(ic)%nghbr = 4
    ! cell(ic)%nghbr(v1) = 8          cell(ic)%nghbr(v1) = 9
    ! cell(ic)%nghbr(v2) = 2          cell(ic)%nghbr(v2) = 0 (boundary)
    ! cell(ic)%nghbr(v3) = 5          cell(ic)%nghbr(v3) = 12
    !                                 cell(ic)%nghbr(v4) = 7
    !..........................................................................
    ! cell(ic)%nghbre(e1) = 5         cell(ic)%nghbre(e1) = 12
    ! cell(ic)%nghbre(e2) = 8         cell(ic)%nghbre(e2) = 7
    ! cell(ic)%nghbre(e3) = 2         cell(ic)%nghbre(e3) = 9
    !                                 cell(ic)%nghbre(e4) = 0 (boundary)
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
    ! assign neighbor cell number of local edge(ie) of cell(ic)
    !--------------------------------------------------------------------------!
    !-- triangle cells
    allocate(locedge(3));  locedge(1)=2;  locedge(2)=3;  locedge(3)=1;
    do ic=1,ncells_tri
      allocate(cell(ic)%nghbre(cell(ic)%nvrt))
      cell(ic)%nghbre(:)=0
      do k=1,cell(ic)%nvrt
        ie=locedge(k)
        cell(ic)%nghbre(ie) = cell(ic)%nghbr(k)
      enddo
    enddo

    !-- quad cells
    deallocate(locedge)
    allocate(locedge(4));  locedge(1)=3;  locedge(2)=4;  locedge(3)=1;  locedge(4)=2
    do i=1,ncells_quad
      ic=i+ncells_tri
      allocate(cell(ic)%nghbre(cell(ic)%nvrt))
      cell(ic)%nghbre(:)=0
      do k=1,cell(ic)%nvrt
        ie=locedge(k)
        cell(ic)%nghbre(ie) = cell(ic)%nghbr(k)
      enddo
    enddo
    deallocate(locedge)


    !--------------------------------------------------------------------------!
    ! count number of edges and allocate edge array
    !--------------------------------------------------------------------------!
    nedges=0
    do ic = 1, ncells
      do j=1,cell(ic)%nnghbrs
        if (cell(ic)%nghbr(j)>ic .or. cell(ic)%nghbr(j)==0) nedges = nedges + 1
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
    tri_v2e(1)=2;  tri_v2e(2)=3;  tri_v2e(3)=1;
    qud_v2e(1)=3;  qud_v2e(2)=4;  qud_v2e(3)=1;  qud_v2e(4)=2

    do ic=1,ncells
      allocate(cell(ic)%edge(cell(ic)%nvrt))
    enddo

    nedges=0
    !--triangle cells
    do ic = 1, ncells_tri
      if ( cell(ic)%nghbr(3) > ic  .or. cell(ic)%nghbr(3)==0 ) then
        nedges = nedges + 1
        cell(ic)%edge(1)= nedges
      else
        jc=cell(ic)%nghbr(3)
        idum4 = 0
        if (cell(jc)%nvrt==3) idum4(1:3) = tri_v2e(1:3)
        if (cell(jc)%nvrt==4) idum4(1:4) = qud_v2e(1:4)
        tlp1: do k=1,cell(jc)%nvrt !3
          if (cell(jc)%nghbr(k)==ic) then
            ie=idum4(k) !locedge(k)
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
       idum4 = 0
       if (cell(jc)%nvrt==3) idum4(1:3) = tri_v2e(1:3)
       if (cell(jc)%nvrt==4) idum4(1:4) = qud_v2e(1:4)
       tlp2: do k=1,cell(jc)%nvrt !3
         if (cell(jc)%nghbr(k)==ic) then
           ie=idum4(k) !locedge(k)
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
       idum4 = 0
       if (cell(jc)%nvrt==3) idum4(1:3) = tri_v2e(1:3)
       if (cell(jc)%nvrt==4) idum4(1:4) = qud_v2e(1:4)
       tlp3: do k=1,cell(jc)%nvrt !3
         if (cell(jc)%nghbr(k)==ic) then
           ie=idum4(k) !locedge(k)
           cell(ic)%edge(3)=cell(jc)%edge(ie)
           exit tlp3
         endif
       enddo tlp3
      endif
    enddo

    !--quad cells
    do i = 1, ncells_quad
      ic=ncells_tri + i
      if ( cell(ic)%nghbr(3) > ic  .or. cell(ic)%nghbr(3)==0 ) then
        nedges = nedges + 1
        cell(ic)%edge(1)= nedges
      else
        jc=cell(ic)%nghbr(3)
        idum4 = 0
        if (cell(jc)%nvrt==3) idum4(1:3) = tri_v2e(1:3)
        if (cell(jc)%nvrt==4) idum4(1:4) = qud_v2e(1:4)
        qlp1: do k=1,cell(jc)%nvrt !4
          if (cell(jc)%nghbr(k)==ic) then
            ie=idum4(k) !locedge(k)
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
        idum4 = 0
        if (cell(jc)%nvrt==3) idum4(1:3) = tri_v2e(1:3)
        if (cell(jc)%nvrt==4) idum4(1:4) = qud_v2e(1:4)
        qlp2: do k=1,cell(jc)%nvrt !4
          if (cell(jc)%nghbr(k)==ic) then
            ie=idum4(k) !locedge(k)
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
        idum4 = 0
        if (cell(jc)%nvrt==3) idum4(1:3) = tri_v2e(1:3)
        if (cell(jc)%nvrt==4) idum4(1:4) = qud_v2e(1:4)
        qlp3: do k=1,cell(jc)%nvrt !4
          if (cell(jc)%nghbr(k)==ic) then
            ie=idum4(k) !locedge(k)
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
        idum4 = 0
        if (cell(jc)%nvrt==3) idum4(1:3) = tri_v2e(1:3)
        if (cell(jc)%nvrt==4) idum4(1:4) = qud_v2e(1:4)
        qlp4: do k=1,cell(jc)%nvrt !4
          if (cell(jc)%nghbr(k)==ic) then
            ie=idum4(k) !locedge(k)
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
    ! build position vector from cell center to cell-edges
    !--------------------------------------------------------------------------!
    do ic=1,ncells
      xc=cell(ic)%x
      yc=cell(ic)%y
      allocate(cell(ic)%pos2edg(cell(ic)%nvrt,2))
      do ie=1,cell(ic)%nvrt
        je=cell(ic)%edge(ie)

        xf=edge(je)%x
        yf=edge(je)%y

        cell(ic)%pos2edg(ie,1) = xf-xc
        cell(ic)%pos2edg(ie,2) = yf-yc
      enddo
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


    !--------------------------------------------------------------------------!
    ! determine interior cells and check #s boundary cells
    !--------------------------------------------------------------------------!
    if (allocated(idum1)) deallocate(idum1)
    if (allocated(idum2)) deallocate(idum2)
    allocate(idum1(ncells), idum2(ncells))
    idum1(:)=0
    idum2(:)=0
    ncells_intr = 0
    ncells_bndr = 0
    do ic=1,ncells
      im=1
      do ie=1,cell(ic)%nvrt
        if (cell(ic)%nghbre(ie)==0) im=0
      enddo
      if (im==1) then
        ncells_intr = ncells_intr + 1
        idum1(ncells_intr) = ic
      else
        ncells_bndr = ncells_bndr + 1
        idum2(ncells_bndr) = ic
      endif
    enddo
    allocate(cell_intr(ncells_intr))
    allocate(cell_bndr(ncells_bndr))
    cell_intr(1:ncells_intr) = idum1(1:ncells_intr)
    cell_bndr(1:ncells_bndr) = idum2(1:ncells_bndr)

    if (ncells_bndr /= sum(bndry(1:nbndries)%ncells)) then
      write(*,*) '#s of boundary cells does not match'
      write(*,*) 'error in mod: grid_procs, sub: grid_data'
      stop
    endif

    !--------------------------------------------------------------------------!
    ! determine interior and boundary edges
    !--------------------------------------------------------------------------!
    if (allocated(idum1)) deallocate(idum1)
    if (allocated(idum2)) deallocate(idum2)
    if (allocated(idum3)) deallocate(idum3)
    allocate(idum1(nedges), idum2(nedges), idum3(nedges))
    idum1(:)=0
    idum2(:)=0
    idum3(:)=0
    nedges_intr = 0
    nedges_bndr = 0
    do ie=1,nedges
      if (edge(ie)%c1 ==0 .and. edge(ie)%c2 >0) then
        nedges_bndr = nedges_bndr + 1
        idum1(nedges_bndr) = ie
        idum3(nedges_bndr) = edge(ie)%c2
      elseif (edge(ie)%c1 >0 .and. edge(ie)%c2 ==0) then
        nedges_bndr = nedges_bndr + 1
        idum1(nedges_bndr) = ie
        idum3(nedges_bndr) = edge(ie)%c1
      elseif  (edge(ie)%c1 >0 .and. edge(ie)%c2 >0) then
        nedges_intr = nedges_intr + 1
        idum2(nedges_intr) = ie
      else
        write(*,*) 'no way !!!'
      endif
    enddo
    allocate(edge_intr(nedges_intr))
    allocate(edge_bndr(nedges_bndr))
    allocate(bedge_cell(nedges_bndr))
    edge_bndr(1:nedges_bndr) = idum1(1:nedges_bndr)
    edge_intr(1:nedges_intr) = idum2(1:nedges_intr)
    bedge_cell(1:nedges_bndr) = idum3(1:nedges_bndr)

    do ib=1,nbndries
      bndry(ib)%nedges = 0
      idum3 = 0
      do i=1,bndry(ib)%ncells
        ic=bndry(ib)%cell(i)
        do ie=1,cell(ic)%nvrt
          je=cell(ic)%edge(ie)
          if (edge(je)%c1 ==ic .and. edge(je)%c2 ==0) then
            bndry(ib)%nedges = bndry(ib)%nedges + 1
            idum3(bndry(ib)%nedges) = je
          elseif (edge(je)%c1 ==0 .and. edge(je)%c2 ==ic) then
            bndry(ib)%nedges = bndry(ib)%nedges + 1
            idum3(bndry(ib)%nedges) = je
          endif
        enddo
      enddo
      allocate(bndry(ib)%edge( bndry(ib)%nedges ))
      bndry(ib)%edge(1:bndry(ib)%nedges) = idum3(1:bndry(ib)%nedges)
    enddo

    if (nedges_bndr /= sum(bndry(1:nbndries)%nedges)) then
      write(*,*) '#s of boundary edges/faces does not match'
      write(*,*) 'error in mod: grid_procs, sub: grid_data'
      stop
    endif

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


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grid_data_verify
    implicit none

    integer   :: i,j, ic,in,iv,ie, jc,je, istat
    real      :: vol,vol1,vol2, xf,nx,af, xc,yc

    !--------------------------------------------------------------------------!
    ! compare computed aread with Green theorem
    !--------------------------------------------------------------------------!
    vol1 = 0.d0
    do ic=1,ncells
      vol1 = vol1 + cell(ic)%vol
    enddo

    vol2= 0.d0
    do ic=1,ncells
      vol=0.d0
      do ie=1,cell(ic)%nvrt
        je=cell(ic)%edge(ie)
        nx= edge(je)%nx * cell(ic)%nrmlsign(ie)
        xf= edge(je)%x
        af= edge(je)%area
        vol = vol +  nx*xf*af
      enddo
      vol2 = vol2 + vol
    enddo


    !--------------------------------------------------------------------------!
    ! write out grid information
    !--------------------------------------------------------------------------!
    if (proc_id==0) then
      open(iunit_log_grid,file=trim(file_log_grid),status='unknown',IOSTAT=istat)

      write(iunit_log_grid,*)
      write(iunit_log_grid,'(a45,i0)') 'number of nodes: ',nnodes
      write(iunit_log_grid,'(a45,i0)') 'number of triangle cells: ',ncells_tri
      write(iunit_log_grid,'(a45,i0)') 'number of quadrilateral cells: ',ncells_quad
      write(iunit_log_grid,'(a45,i0)') 'number of total cells: ',ncells
      write(iunit_log_grid,'(a45,i0)') 'number of total edges: ',nedges
      write(iunit_log_grid,'(a45,i0)') 'number of total bondary edges: ',nedges_bndr
      write(iunit_log_grid,'(a45,i0)') 'number of total bondary cells: ',ncells_bndr

      write(iunit_log_grid,*)
      write(iunit_log_grid,'(a44,en20.11)') " Sum of the cell volumes via numerical cal:", vol1
      write(iunit_log_grid,'(a44,en20.11)') " Sum of the cell volumes via Green theorem:", vol2
      write(iunit_log_grid,'(a45,en20.11)') 'min cell volume: ',minval(cell(1:ncells)%vol)
      write(iunit_log_grid,'(a45,en20.11)') 'max cell volume: ',maxval(cell(1:ncells)%vol)

      write(iunit_log_grid,*)
      write(iunit_log_grid,'(a,en20.11)') " cell effective length, sqrt[sum(vol)/ncells]:", heff1
      write(iunit_log_grid,'(a,en20.11)') " cell effective length, sum[sqrt(vol)]/ncells:", heff2
    endif


    !--------------------------------------------------------------------------!
    ! write out cell normals
    !--------------------------------------------------------------------------!
    if (proc_id==0) then
      open(iunit_log_grid,file='log_grid_cell_normals.plt', status='unknown',IOSTAT=istat)
      write(iunit_log_grid,'(a)') 'VARIABLES ="x", "y", "nx", "ny"'
      do ic=1,ncells
        do ie=1,cell(ic)%nvrt
          je=cell(ic)%edge(ie)
          xc=edge(je)%nx * cell(ic)%nrmlsign(ie)
          yc=edge(je)%ny * cell(ic)%nrmlsign(ie)
          write(iunit_log_grid,*) cell(ic)%x, cell(ic)%y, xc, yc
        enddo
      enddo
      close(iunit_log_grid)
    endif

    return
  end subroutine grid_data_verify

end module grid_procs
