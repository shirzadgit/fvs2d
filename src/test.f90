module test

  use input
  use data_type
  use grid_procs
  use interpolation
  use gradient

  implicit none

  real,allocatable,save :: fve(:),dfve(:,:), fce(:),dfce(:,:)

contains

  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine test_init
    implicit none

    allocate(fve(nnodes), fce(ncells), dfve(nnodes,2), dfce(ncells,2))

    call test_analytic

    call test_grid

    call test_interpolation

    call test_gradient

    return
  end subroutine test_init



  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine test_grid
    implicit none

    integer             :: ic,ic1,ic2, in,in1,in2, ie,ie1,ie2, je,jc, i,j,n1,nt,iv,iv1,iv2, num_nearests, iptr, vL,vR
    integer,allocatable :: iwrk1d(:),iwrk2d(:,:)
    real,allocatable    :: fc(:),fv(:),fv_interp(:) !!, interp_ggnb_cellweight(:), interp_ggnb_nodeweight(:)
    real                :: vol, r(2), d,dt,xc,yc
    type(kdtree2_result),allocatable	:: fnearest(:)

    write(*,*)
    write(*,'(a,i0)') 'number of nodes: ',nnodes
    write(*,'(a,i0)') 'number of cells-triangle: ',ncells_tri
    write(*,'(a,i0)') 'number of cells-quad: ',ncells_quad
    write(*,'(a,i0)') 'number of cells: ',ncells
    write(*,'(a,i0)') 'number of edges: ',nedges
    ! write(*,'(a,i0)') 'number of bondary nodes: ',num_nodes_bndry
    ! write(*,'(a,i0)') 'number of bondary edges: ',num_edges_bndry
    ! write(*,'(a,i0)') 'number of bondary cells: ',num_cells_bndry
    write(*,'(a,en20.11)') 'min cell volume: ',minval(cell(1:ncells)%vol)
    write(*,'(a,en20.11)') 'max cell volume: ',maxval(cell(1:ncells)%vol)
    write(*,*)

    open(100,file='grid_cell_normals.dat')
    write(100,'(a)') 'VARIABLES ="X", "Y", "N1", "N2"'
    do ic=1,ncells
      do ie=1,cell(ic)%nvrt
        je=cell(ic)%edge(ie)
        xc=edge(je)%nx * cell(ic)%nrmlsign(ie)
        yc=edge(je)%ny * cell(ic)%nrmlsign(ie)
        write(100,*) cell(ic)%x, cell(ic)%y, xc, yc
      enddo
    enddo
    close(100)

    open(100,file='grid_cell_neighbor.plt')
    do ic=1,ncells,4
      write(100,'(a)') 'VARIABLES ="X", "Y"'
      write(100,'(a)') 'zone i=1, j=1, f=point'
      write(100,*) cell(ic)%x,cell(ic)%y
      do i=1,cell(ic)%nnghbrs
        if (i>4) write(*,*) 'wow'
        jc=cell(ic)%nghbr(i)
        !write(*,*) jc
        if (jc>0) then
          write(100,'(a)') ' '
          write(100,'(a)') 'TITLE ="grid"'
          write(100,'(a)') 'VARIABLES ="X", "Y"'
          if (cell(jc)%nvrt==3) write(100,'(a,i0,a)') 'ZONE T="VOL_MIXED",N=', cell(jc)%nvrt, ' E=1, ET=TRIANGLE F=FEPOINT'
          if (cell(jc)%nvrt==4) write(100,'(a)') 'ZONE N=4, ELEMENTS=1, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL'
          do j=1,cell(jc)%nvrt
            write(100,*) node(cell(jc)%node(j))%x, node(cell(jc)%node(j))%y
          enddo
          do j=1,cell(jc)%nvrt
            write(100,*) j
          enddo
        endif
      enddo
    enddo
    close(100)

    ! open(100,file='del.dat')
    ! do ic=1,ncells
    !   do iv=1,cell(ic)%nvrt
    !
    !     if (iv  < cell(ic)%nvrt) vL = cell(ic)%node(iv+1)
    !     if (iv == cell(ic)%nvrt) vL = cell(ic)%node(1)
    !     if (cell(ic)%nghbr(iv)==0) write(100,*) node(vL)%x, node(vL)%y
    !   enddo
    ! enddo
    ! close(100)

    return



    ! open(100,file='neighbor.plt')
    ! nt=0
    ! do ic=1,20
    !   write(100,'(a)') 'VARIABLES ="X", "Y"'
    !   write(100,*) 'zone i=',1, ', j=1, f=point'
    !   write(100,*) cell_center(ic,1),cell_center(ic,2)
    !
    !   write(100,'(a)') 'VARIABLES ="X", "Y"'
    !   write(100,*) 'zone i=',cell_neighbr_nn_num(ic), ', j=1, f=point'
    !   do iptr=nt+1,nt+cell_neighbr_nn_num(ic)
    !     ic1=cell_neighbr_nn_ptr(iptr)
    !     write(100,*) cell_center(ic1,1),cell_center(ic1,2)
    !   enddo
    !   nt=nt+cell_neighbr_nn_num(ic)
    !   write(100,*)
    ! enddo
    ! return
    !
    !
    ! open(100,file='interp.plt')
    ! write(100,'(a)') 'TITLE ="grid"'
    ! write(100,'(a)') 'VARIABLES ="X", "Y", "U_intp", "U"'
    ! write(100,'(a,i0,a,i0,a)') 'ZONE T="VOL_MIXED",N=', num_nodes, ' E=',num_cells, ' ET=TRIANGLE F=FEPOINT'
    ! do i=1,num_nodes
    !   write(100,'(4(e19.8,1x))') cell_vert(i,1),cell_vert(i,2), fv_interp(i),fv(i)
    ! enddo
    ! do i=1,num_cells
    !     write(100,*) (cell2node(i,j),j=1,3)
    ! enddo
    ! close(100)
    !
    !
    ! write(*,*) node2cell_ntot(601)
    ! open(100,file='neighbor.plt')
    ! do in=1,1001,50
    !   write(100,'(a)') 'TITLE ="grid"'
    !   write(100,'(a)') 'VARIABLES ="X", "Y"'
    !   write(100,*) cell_vert(in,1),cell_vert(in,2)
    !   do i=node2cell_ptr(in), node2cell_ptr(in) + node2cell_ntot(in) - 1
    !     ic=node2cell(i)
    !     if (ic>0) then
    !       write(100,'(a)') ' '
    !       write(100,'(a)') 'TITLE ="grid"'
    !       write(100,'(a)') 'VARIABLES ="X", "Y"'
    !       write(100,'(a)') 'ZONE T="VOL_MIXED",N= 3, E=1, ET=TRIANGLE F=FEPOINT'
    !       do j=1,3
    !         write(100,*) cell_vert( cell2node(ic,j) ,1),cell_vert( cell2node(ic,j) ,2)
    !       enddo
    !       write(100,*) 1,2,3 !cell2node(ic,1),cell2node(ic,2),cell2node(ic,3)
    !     endif
    !   enddo
    ! enddo
    ! close(100)
    !
    ! open(100,file='grid_vol.dat')
    ! do ic=1,num_cells
    !   write(100,*) ic,cell_area(ic)
    ! enddo
    ! close(100)
    !
    ! open(100,file='grid_nodes.dat')
    ! write(100,'(a)') 'VARIABLES ="X", "Y"'
    ! do i=1,num_nodes
    !   write(100,*) cell_vert(i,1),cell_vert(i,2)
    ! enddo
    ! close(100)
    !
    ! open(100,file='grid_edge_center.dat')
    ! write(100,'(a)') 'VARIABLES ="X", "Y"'
    ! do i=1,num_edges
    !   write(100,*) edge_center(i,1),edge_center(i,2)
    ! enddo
    ! close(100)
    !
    ! open(100,file='grid_cell_center.dat')
    ! write(100,'(a)') 'VARIABLES ="X", "Y"'
    ! do i=1,num_cells
    !   write(100,*) cell_center(i,1),cell_center(i,2)
    ! enddo
    ! close(100)
    !
    ! open(100,file='grid_cell_normals.dat')
    ! write(100,'(a)') 'VARIABLES ="X", "Y", "N1", "N2"'
    ! do ic=1,num_cells
    !   do ie=1,num_vert_cell(ic)
    !     ie1=cell2edge(ic,ie)
    !     xc=edge_normal(ie1,1) * edge_normal_sign(ic,ie)
    !     yc=edge_normal(ie1,2) * edge_normal_sign(ic,ie)
    !     write(100,*) cell_center(ic,1),cell_center(ic,2),xc,yc
    !   enddo
    ! enddo
    ! close(100)
    !
    ! open(100,file='grid_con.dat')
    ! do i=1,num_cells
    !     write(100,*) (cell2node(i,j),j=1,3)
    ! enddo
    ! close(100)


    return
  end subroutine test_grid


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine test_interpolation
    implicit none

    integer             :: ic,ic1,ic2, in,in1,in2, ie,ie1,ie2, i,j,n1,nt,iv,iv1,iv2, num_nearests
    integer,allocatable :: iwrk1d(:),iwrk2d(:,:)
    real,allocatable    :: fv(:),fw(:,:)
    real                :: vol, r(2), d,dt,xc,yc
    character           :: fout*124

    allocate(fv(nnodes))


    call interpolate_cell2node(fce,fv)

    fout='test_interp.plt'
    allocate(fw(nnodes,1))
    fw(:,1)=fv(:)
    if (ncells_quad>0) then
      call test_tecplot_mixed(fout,1,fw)
    else
      call test_tecplot(fout,1,fw)
    endif
    deallocate(fw)


    return
  end subroutine test_interpolation


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine test_gradient
    implicit none
    integer             :: ic,ic1,ic2, in,in1,in2, ie,ie1,ie2, i,j,n1,nt,iv,iv1,iv2, num_nearests
    integer,allocatable :: iwrk1d(:),iwrk2d(:,:), tmpi(:), tmpi2d(:,:)
    real,allocatable    :: f(:), fxv(:), fyv(:), df(:,:), fw(:,:)
    real                :: vol, r(2), d,dt,xc,yc, nxf,nyf,af, ax,ay
    character           :: fout*124

    allocate(f(ncells),fxv(nnodes), fyv(nnodes), df(ncells,2))


!    call gradient_cellcntr('ggcb',fce,df)
!    call interpolate_cellcntr2node(df(1:num_cells,1),fxv)
!    call interpolate_cellcntr2node(df(1:num_cells,2),fyv)

    call gradient_cellcntr(fce,df)
    call interpolate_cell2node(df(1:ncells,1),fxv)
    call interpolate_cell2node(df(1:ncells,2),fyv)
    fout='test_gradient.plt'
    allocate(fw(nnodes,2))
    fw(:,1)=fxv(:)
    fw(:,2)=fyv(:)
    call test_tecplot(fout,2,fw)


    call interpolate_cell2node(dfce(1:ncells,1),fxv)
    call interpolate_cell2node(dfce(1:ncells,2),fyv)
    fout='test_gradient_exact.plt'
    fw(:,1)=fxv(:)
    fw(:,2)=fyv(:)
    call test_tecplot(fout,2,fw)
    deallocate(fw)


    return
  end subroutine test_gradient


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine test_tecplot (fout,nvar,f)
    implicit none

    integer,intent(in)  :: nvar
    real,intent(in)     :: f(nnodes,nvar)
    character(len=*)    :: fout

    integer             :: i,j,k

    open(100,file=trim(fout))
    write(100,'(a)') 'TITLE ="grid_sol"'
    if (nvar==1) then
      write(100,'(a)') 'VARIABLES ="x", "y", "f1"'
      write(100,'(a,i0,a,i0,a)') 'ZONE T="VOL_MIXED",N=', nnodes, ' E=',ncells, ' ET=TRIANGLE F=FEPOINT'
      do i=1,nnodes
        write(100,'(3(e19.8,1x))') node(i)%x,node(i)%y, f(i,1)
      enddo
    elseif (nvar==2) then
      write(100,'(a)') 'VARIABLES ="x", "y", "f1", "f2"'
      write(100,'(a,i0,a,i0,a)') 'ZONE T="VOL_MIXED",N=', nnodes, ' E=',ncells, ' ET=TRIANGLE F=FEPOINT'
      do i=1,nnodes
        write(100,'(4(e19.8,1x))') node(i)%x,node(i)%y, f(i,1),f(i,2)
      enddo
    else
      write(100,'(a)') 'VARIABLES ="x", "y", "f1", "U"'
      write(100,'(a,i0,a,i0,a)') 'ZONE T="VOL_MIXED",N=', nnodes, ' E=',ncells, ' ET=TRIANGLE F=FEPOINT'
      do i=1,nnodes
        write(100,'(5(e19.8,1x))') node(i)%x,node(i)%y, f(i,1),f(i,2), f(i,3)
      enddo
    endif
    do i=1,ncells_tri
        write(100,*) (cell(i)%node(j),j=1,3)
    enddo
    do i=1,ncells_quad
      k=i+ncells_tri
        write(100,*) (cell(k)%node(j),j=1,4)
    enddo
    close(100)


    return
  end subroutine test_tecplot


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine test_tecplot_mixed (fout,nvar,f)
    implicit none

    integer,intent(in)  :: nvar
    real,intent(in)     :: f(nnodes,nvar)
    character(len=*)    :: fout

    integer             :: i,j,k

    open(100,file=trim(fout))
    write(100,'(a)') 'TITLE ="grid_sol"'
    if (nvar==1) then
      write(100,'(a)') 'VARIABLES ="x", "y", "f1"'
      write(100,'(a,i0,a,i0,a)') 'ZONE NODES=',nnodes, ' ELEMENTS=',ncells, ' DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL'
      do i=1,nnodes
        write(100,'(3(e19.8,1x))') node(i)%x,node(i)%y, f(i,1)
      enddo
    elseif (nvar==2) then
      write(100,'(a)') 'VARIABLES ="x", "y", "f1", "f2"'
      write(100,'(a,i0,a,i0,a)') 'ZONE NODES=',nnodes, ' ELEMENTS=',ncells, ' DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL'
      do i=1,nnodes
        write(100,'(4(e19.8,1x))') node(i)%x,node(i)%y, f(i,1),f(i,2)
      enddo
    else
      write(100,'(a)') 'VARIABLES ="x", "y", "f1", "U"'
      write(100,'(a,i0,a,i0,a)') 'ZONE NODES=',nnodes, ' ELEMENTS=',ncells, ' DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL'
      do i=1,nnodes
        write(100,'(5(e19.8,1x))') node(i)%x,node(i)%y, f(i,1),f(i,2), f(i,3)
      enddo
    endif
    do i=1,ncells_tri
        write(100,*) cell(i)%node(1),cell(i)%node(2),cell(i)%node(3),cell(i)%node(3)
    enddo
    do i=1,ncells_quad
      k=i+ncells_tri
        write(100,*) (cell(k)%node(j),j=1,4)
    enddo
    close(100)


    return
  end subroutine test_tecplot_mixed


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine test_analytic
    implicit none

    integer :: ic,in
    real    :: ax,ay, xc,yc, xmin,xmax, ymin,ymax

    xmin=minval(node(:)%x)
    xmax=maxval(node(:)%x)
    ax=10.d0*acos(-1.d0)/(xmax-xmin)

    ymin=minval(node(:)%y)
    ymax=maxval(node(:)%y)
    ay=10.d0*acos(-1.d0)/(ymax-ymin)

    do ic=1,ncells
      xc=cell(ic)%x
      yc=cell(ic)%y
      fce(ic)=sin(ax*xc)*cos(ay*yc)
      dfce(ic,1)= ax*cos(ax*xc)*cos(ay*yc)
      dfce(ic,2)=-ay*sin(ax*xc)*sin(ay*yc)
    enddo

    do in=1,nnodes
      xc=node(in)%x
      yc=node(in)%y
      fve(in)=sin(ax*xc)*cos(ay*yc)
      dfve(in,1)= ax*cos(ax*xc)*cos(ay*yc)
      dfve(in,2)=-ay*sin(ax*xc)*sin(ay*yc)
    enddo

    return
  end subroutine test_analytic


end module test
