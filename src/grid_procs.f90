module grid_procs

  use mainparam, only : iunit_grid
  use input, only     : file_grid
  use grid_graph

  implicit none

  private
  integer,parameter                       :: nvmax=3
  integer,save                            :: ncells,nnodes,nedges,nbedges
  real,allocatable,dimension(:),save      :: xv,yv, xc,yc, xf,yf, snx,sny, nedgpcell, nverpedge, vol, af
  integer,allocatable,dimension(:,:),save :: lcell2edge, ledge2cell, ledge2node, lcell2node, cellconct, edgemap


  public		:: grid_procs_init

  !================================================================
  ! nvmax:	maximum number of vertices per cell (set to 3, can change later) 
  ! nnodes:	number of nodes/vertices
  ! ncells:	number of cells
  ! nedges:	number of edges
  ! nbedges:	number of boundary edges
  ! nedgpcell:number of edges of each cell
  ! nverpedge:number of nodes of each edge
  ! xv,yv:	x/y coordinates of vertices
  ! xc,yc:	x/y coordinates of cell center
  ! xf,yf:	x/y coordiantes of edge center
  ! snx,sny	edge normal in x/y
  ! af:		edge area
  ! vol:		cell volume
  !================================================================


  contains


  !====================================================================================================================================================
	subroutine grid_procs_init
    implicit none

    !-------------------------------------------------------------------------------
    ! read grid file
    !-------------------------------------------------------------------------------
    call grid_read 
    

    !-------------------------------------------------------------------------------
    ! read grid file
    !-------------------------------------------------------------------------------    
    call grid_link
    
    write(*,*) 'number of nodes: ',nnodes
    write(*,*) 'number of cells: ',ncells
		
    return
  end subroutine grid_procs_init
  
  
  !====================================================================================================================================================
  subroutine grid_read 
    implicit none

    integer					:: i,istat,nloc1,nloc2,eloc1,eloc2
    character				:: dchar*400
	
    !-------------------------------------------------------------------------------
    ! open grid file
    !-------------------------------------------------------------------------------
    open (iunit_grid,file=trim(file_grid),status='unknown',IOSTAT=istat)	
    read (iunit_grid,*)
    read (iunit_grid,*)
    read (iunit_grid,'(a)') dchar

    !-------------------------------------------------------------------------------
    ! determine number of nodes and cells
    !-------------------------------------------------------------------------------	
    do i=1,len(dchar)-2
      if (dchar(i:i+1)=='N=' .or. dchar(i:i+1)=='n=')		nloc1=i+2
      if (dchar(i:i+1)=='E=' .or. dchar(i:i+1)=='e=')		nloc2=i-1
      if (dchar(i:i+1)=='E=' .or. dchar(i:i+1)=='e=')		eloc1=i+2
      if (dchar(i:i+2)=='ET=' .or. dchar(i:i+2)=='et=')	eloc2=i-1	
    enddo			
    read (dchar(nloc1:nloc2),*,iostat=istat) nnodes
    read (dchar(eloc1:eloc2),*,iostat=istat) ncells

    !-------------------------------------------------------------------------------	
    ! allocate and read nodes/vertices coordinates
    !-------------------------------------------------------------------------------	
    allocate(xv(nnodes),yv(nnodes))
    do i=1,nnodes
      read (iunit_grid,*) xv(i),yv(i)
    enddo

    !-------------------------------------------------------------------------------	
    ! allocate and read cell connectivity
    !-------------------------------------------------------------------------------	
    allocate(cellconct(nvmax,ncells))
    do i=1,ncells
      read (iunit_grid,*) cellconct(1,i),cellconct(2,i),cellconct(3,i)
    enddo

    return
  end subroutine grid_read
  
  
  !====================================================================================================================================================
  subroutine grid_link 
    implicit none

    integer,allocatable   :: edgetmp(:,:)

    nedges=0 
    allocate(edgetmp(6,4*ncells))
    call findEdges(nvmax,ncells,nedges,cellconct,edgetmp)
    allocate(edgemap(6,nedges));
    edgemap(1:6,1:nedges)=edgetmp(1:6,1:nedges)
    deallocate(edgetmp)
    !edge(1:4,1:nedges) = [node1,node2,cellLeft,cellRight]

    
    return
  end subroutine grid_link


end module grid_procs
