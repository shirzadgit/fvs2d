module grid_graph
contains
  !
  ! find all edges in a given
  !
  ! edge(1:4,i) = [node1,node2,cellLeft,cellRight]
  ! id = pointer to the next edge that connects to node1
  ! 
  ! iptr(:)  = pointer to the last edge inserted for each node
  !
  ! unstructured grid
  !
  subroutine findEdges(nv,ncells,nedges,cellCon,edgeout)
    !
    implicit none
    !
    ! subroutine arguments
    !
    integer, intent(in) :: nv
    integer,intent(in) :: ncells
    integer, intent(in) :: cellCon(nv,ncells)
    integer, intent(inout) :: nedges
    integer, allocatable, intent(out) :: edgeout(:,:)
    !
    ! local variables
    !
    integer :: i,j,jp1,maxedges,m
    integer, allocatable :: iflag(:),etmp(:,:),iptr(:)
    integer :: eloc(2)
    integer :: nvert,nnodes
    !
    ! begin
    !
    maxedges=ncells*4
    !
    nnodes=0
    do i=1,ncells
       do m=1,nv
          nnodes=max(nnodes,cellCon(m,i))
       enddo
    enddo
    !
    allocate(iptr(nnodes))
    allocate(etmp(7,maxedges))
    allocate(iflag(ncells))
    !
    etmp=0
    !
    iptr=0
    nedges=0
    !
    do i=1,ncells
       nvert=nv
       if (nvert > 3) then
          if (cellCon(3,i)==cellCon(4,i)) nvert=3
       endif
       do j=1,nvert
          jp1=mod(j,nvert)+1
          eloc(1)=cellCon(j,i)
          eloc(2)=cellCon(jp1,i)
          call insert_edge(iptr,etmp,eloc,i,j,nedges,nnodes,maxedges)
       enddo
    enddo
    !
    allocate(edgeout(6,nedges))
    do i=1,nedges
       edgeout(:,i)=etmp(1:6,i)
    enddo
    !
    deallocate(etmp)
    deallocate(iflag)
    deallocate(iptr)
    !
    return
  end subroutine findEdges
  !>
  !> support routine for findEdges
  !>
  subroutine insert_edge(iptr,edge,eloc,cellIndex,edgeIndex,nedges,nnodes,maxedges)
    implicit none
    !
    ! subroutine arguments
    !
    integer, intent(in) :: nnodes,maxedges,eloc(2),cellIndex,edgeIndex
    integer, intent(inout) :: nedges
    integer, intent(inout) :: iptr(nnodes)
    integer, intent(inout) :: edge(7,maxedges)
    !
    ! local variables
    !
    integer :: e1(2),e2(2),ip,te
    !
    ! begin
    !
    e1=eloc
    if (e1(1).gt.e1(2)) then
       te=e1(1)
       e1(1)=e1(2)
       e1(2)=te
    endif
    !
    ip=iptr(e1(1))
    !
    checkloop: do while(ip > 0)
       e2=edge(1:2,ip)
       if (e2(1).gt.e2(2)) then
          te=e2(1)
          e2(1)=e2(2)
          e2(2)=te
       endif
       if (sum(abs(e1-e2))==0) then
          edge(4,ip)=cellIndex
          edge(6,ip)=edgeIndex
          return
       endif
       ip=edge(7,ip)
    enddo checkloop
    !
    nedges=nedges+1
    edge(1:2,nedges)=eloc
    edge(3,nedges)=cellIndex
    edge(5,nedges)=edgeIndex
    edge(7,nedges)=iptr(e1(1))
    iptr(e1(1))=nedges
    !
    return
  end subroutine insert_edge  
  !>
  !> find node2node and node2edge map from edges
  !>
  subroutine findnodemap(edge,node2node,node2edge,node2nodeptr,nnode,nedges)
    implicit none
    integer, intent(in) :: nnode,nedges
    integer, intent(in) :: edge(6,nedges)
    integer, allocatable, intent(out) :: node2node(:)
    integer, allocatable, intent(out) :: node2edge(:)
    integer, allocatable, intent(out) :: node2nodeptr(:)
    !
    integer :: i,nsave,isum
    !
    allocate(node2nodeptr(nnode+1))
    allocate(node2node(2*nedges))
    allocate(node2edge(2*nedges))
    do i=1,nnode+1
      node2nodeptr(i)=0
    enddo
    !
    do i=1,nedges
       node2nodeptr(edge(1,i))=node2nodeptr(edge(1,i))+1
       node2nodeptr(edge(2,i))=node2nodeptr(edge(2,i))+1
    enddo
    !
    isum=1
    do i=1,nnode
       isum=isum+node2nodeptr(i)
       node2nodeptr(i)=isum
    enddo
    node2nodeptr(nnode+1)=isum
    !
    do i=1,nedges
       node2nodeptr(edge(1,i))=node2nodeptr(edge(1,i))-1
       node2node(node2nodeptr(edge(1,i)))=edge(2,i)
       node2edge(node2nodeptr(edge(1,i)))=i
       node2nodeptr(edge(2,i))=node2nodeptr(edge(2,i))-1
       node2node(node2nodeptr(edge(2,i)))=edge(1,i)
       node2edge(node2nodeptr(edge(2,i)))=i
    enddo
    !
  end subroutine findnodemap
  !>
  !> find node2cell map from cell2node
  !>
  subroutine findnode2cellmap(cell2node,node2cell,node2cellptr,nv,nnode,ncells)
    implicit none
    integer, intent(in) :: nv,nnode,ncells
    integer, intent(in) :: cell2node(nv,ncells)
    integer, allocatable, intent(out) :: node2cell(:)
    integer, allocatable, intent(out) :: node2cellptr(:)
    !
    integer :: i,nsave,isum,nvi,j
    integer, dimension(nv) :: iv
    !
    allocate(node2cellptr(nnode+1))
    do i=1,nnode+1
       node2cellptr(i)=0
    enddo
    !
    do i=1,ncells
       !
       ! make sure degenerate elements don't get counted twice
       call unique(cell2node(:,i),iv,nv,nvi)
       !
       do j=1,nvi
          node2cellptr(iv(j))=node2cellptr(iv(j))+1
       enddo
    end do
    !
    isum=1
    do i=1,nnode
       isum=isum+node2cellptr(i)
       node2cellptr(i)=isum
    enddo
    node2cellptr(nnode+1)=isum
    allocate(node2cell(isum-1))
    !
    do i=1,ncells
       call unique(cell2node(:,i),iv,nv,nvi)
       do j=1,nvi
          node2cellptr(iv(j))=node2cellptr(iv(j))-1
          node2cell(node2cellptr(iv(j)))=i
       enddo
    enddo
    !
  end subroutine findnode2cellmap
  !>
  !> find cell2cell and cell2edge map from edges (faces in 3D)
  !>
  subroutine findcell2edgemap(edge,cell2cell,cell2edge,nve,nedges,ncells)
    implicit none
    integer, intent(in) :: nve,nedges,ncells
    integer, intent(in) :: edge(6,nedges)
    integer, allocatable, intent(out) :: cell2cell(:,:)
    integer, allocatable, intent(out) :: cell2edge(:,:)
    !
    integer :: i
    !
    allocate(cell2cell(nve,ncells))
    allocate(cell2edge(nve,ncells))
    cell2cell=0
    cell2edge=0
    !
    do i=1,nedges
       cell2cell(edge(5,i),edge(3,i))=edge(4,i)
       cell2edge(edge(5,i),edge(3,i))=i
       if (edge(4,i) > 0) then
          cell2cell(edge(6,i),edge(4,i))=edge(3,i)
          cell2edge(edge(6,i),edge(4,i))=i
       endif
    enddo
    !
  end subroutine findcell2edgemap    
  !>
  subroutine unique(vec,vec_unique,n,num)
    ! Return only the unique values from vec.
    
    implicit none
    integer, intent(in) :: n
    integer, intent(inout) :: num
    integer,dimension(n),intent(in) :: vec
    integer,dimension(n),intent(out) :: vec_unique
    integer :: i,j
    logical,dimension(n) :: mask
    
    mask = .false.
    
    do i=1,n
       !count the number of occurrences of this element:
       num = count( vec(i)==vec )
        if (num==1) then
          !there is only one, flag it:
          mask(i) = .true.
       else
          !flag this value only if it hasn't already been flagged:
          if (.not. any(vec(i)==vec .and. mask) ) mask(i) = .true.
       end if
    end do
    !return only flagged elements:    
    num=0
    do i=1,n
       if (mask(i)) then
          num=num+1
          vec_unique(num)=vec(i)
       endif
    enddo
  end subroutine unique
  !
  subroutine matchEdges(ngroup,groupcft,edge,nedges)
  !
    implicit none
    integer, intent(inout) :: ngroup
    integer, intent(in)    :: nedges
    integer, allocatable, intent(out) :: groupcft(:)
    integer, intent(inout)    :: edge(6,nedges)
    !
    integer :: i,j,k,nnode,kedge,e1,e2
    integer, allocatable :: etmp(:,:)
    integer, allocatable :: node2edge(:),node2node(:),node2nodeptr(:)
    logical, allocatable :: selected(:)
    !
    if (nedges < 1) then
       ngroup=0
       return
    endif
    !
    nnode=0
    do i=1,nedges
       nnode=max(nnode,edge(1,i))
       nnode=max(nnode,edge(2,i))
    enddo
    !
    call findnodemap(edge,node2node,node2edge,node2nodeptr,nnode,nedges)
    !
    allocate(etmp(6,nedges))
    allocate(selected(nedges))    
    selected=.false.
    !
    kedge=1
    etmp(:,kedge)=edge(:,1)
    e1=etmp(1,kedge)
    selected(1)=.true.
    ngroup=0
    !
    do while(kedge < nedges)
       e2=abs(etmp(2,kedge))
       do j=node2nodeptr(e2),node2nodeptr(e2+1)-1
          k=node2edge(j)
          if (.not.selected(k)) then
             if (edge(1,k)==e2) then
                kedge=kedge+1
                etmp(:,kedge)=edge(:,k)               
             else
                kedge=kedge+1
                etmp(1,kedge)=edge(2,k)
                etmp(2,kedge)=edge(1,k)
                etmp(3,kedge)=edge(4,k)
                etmp(4,kedge)=edge(3,k)
                etmp(5,kedge)=edge(6,k)
                etmp(6,kedge)=edge(5,k)
             endif
             selected(k)=.true.
             exit
          endif
       enddo
       if (etmp(2,kedge)==e1) then
          ngroup=ngroup+1
          etmp(2,kedge)=-e1
          if (kedge < nedges) then
             do j=1,nedges
                if (.not. selected(j)) then
                   kedge=kedge+1
                   etmp(:,kedge)=edge(:,j)
                   e1=etmp(1,kedge)
                   exit
                endif
             enddo
          endif
       endif
    enddo
    allocate(groupcft(ngroup+1))
    groupcft(1)=1
    j=1
    do i=1,nedges
       edge(:,i)=etmp(:,i)
       if (edge(2,i) < 0) then
          edge(2,i)=-edge(2,i)
          j=j+1
          groupcft(j)=i+1
       endif
    enddo
    !
    if (allocated(etmp)) deallocate(etmp)
    if (allocated(node2edge)) deallocate(node2edge)
    if (allocated(node2node)) deallocate(node2node)
    if (allocated(node2nodeptr)) deallocate(node2nodeptr)
    if (allocated(selected)) deallocate(selected)
    !
  end subroutine matchEdges


end module grid_graph
