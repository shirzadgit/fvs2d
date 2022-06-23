module test

  use input
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

    allocate(fve(num_nodes), fce(num_cells), dfve(num_nodes,2), dfce(num_cells,2))

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

    integer             :: ic,ic1,ic2, in,in1,in2, ie,ie1,ie2, i,j,n1,nt,iv,iv1,iv2, num_nearests
    integer,allocatable :: iwrk1d(:),iwrk2d(:,:)
    real,allocatable    :: fc(:),fv(:),fv_interp(:) !!, interp_ggnb_cellweight(:), interp_ggnb_nodeweight(:)
    real                :: vol, r(2), d,dt,xc,yc
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
    return

    do i=1,num_edges_bndry
      ie=edges_bndry_ptr(i)
      if (edge2cell(ie,1)>0) then
        ic=edge2cell(ie,1)
      elseif (edge2cell(ie,2)>0) then
        ic=edge2cell(ie,2)
      else
        write(*,*) 'error'
      endif
    enddo


!    do ic=1,num_cells
!      xc=cell_center(ic,1)
!      yc=cell_center(ic,2)
!      fc=sin(1.25d0*xc)*cos(1.25d0*yc)
!
!      xc=cell_vert(cell2node(ic,1),1); yc=cell_vert(cell2node(ic,1),2); f1=sin(1.25d0*xc)*cos(1.25d0*yc)
!      xc=cell_vert(cell2node(ic,2),1); yc=cell_vert(cell2node(ic,2),2); f2=sin(1.25d0*xc)*cos(1.25d0*yc)
!      xc=cell_vert(cell2node(ic,3),1); yc=cell_vert(cell2node(ic,3),2); f3=sin(1.25d0*xc)*cos(1.25d0*yc)
!    enddo


    allocate(fc(num_cells),fv(num_nodes), fv_interp(num_nodes))
    do in=1,num_nodes
      xc=cell_vert(in,1)
      yc=cell_vert(in,2)
      fv(in)=sin(1.25d0*xc)*cos(1.25d0*yc)
    enddo

    do ic=1,num_cells
      xc=cell_center(ic,1)
      yc=cell_center(ic,2)
      fc(ic)=sin(1.25d0*xc)*cos(1.25d0*yc)
    enddo


    !--------------------------------------------------------------------------!
    !
    !--------------------------------------------------------------------------!
    fv_interp(:)=0.d0
    do in=1,num_nodes
      dt=0.d0
      do i=node2cell_ptr(in), node2cell_ptr(in) + node2cell_ntot(in) - 1
        ic=node2cell(i)
        do iv=1,num_vert_cell(ic)
          if (in==cell2node(ic,iv)) then
            d=cell_dist2vert(ic,iv)
            fv_interp(in) = fv_interp(in) + 1.d0/d*fc(ic)
            dt=dt+1.d0/d
          endif
        enddo
      enddo
      fv_interp(in) = fv_interp(in)/dt
    enddo

    !fv_interp(:)=0.d0
    !do in=1,num_nodes
    !  do i=node2cell_ptr(in), node2cell_ptr(in) + node2cell_ntot(in) - 1
    !    ic=node2cell(i)
    !    fv_interp(in) = fv_interp(in) + interp_ggnb_cellweight(i)*fc(ic)
    !  enddo
    !  fv_interp(in) = interp_ggnb_nodeweight(in)*fv_interp(in)
    !enddo


    open(100,file='interp.plt')
    write(100,'(a)') 'TITLE ="grid"'
    write(100,'(a)') 'VARIABLES ="X", "Y", "U_intp", "U"'
    write(100,'(a,i0,a,i0,a)') 'ZONE T="VOL_MIXED",N=', num_nodes, ' E=',num_cells, ' ET=TRIANGLE F=FEPOINT'
    do i=1,num_nodes
      write(100,'(4(e19.8,1x))') cell_vert(i,1),cell_vert(i,2), fv_interp(i),fv(i)
    enddo
    do i=1,num_cells
        write(100,*) (cell2node(i,j),j=1,3)
    enddo
    close(100)


    write(*,*) node2cell_ntot(601)
    open(100,file='neighbor.plt')
    do in=1,1001,50
      write(100,'(a)') 'TITLE ="grid"'
      write(100,'(a)') 'VARIABLES ="X", "Y"'
      write(100,*) cell_vert(in,1),cell_vert(in,2)
      do i=node2cell_ptr(in), node2cell_ptr(in) + node2cell_ntot(in) - 1
        ic=node2cell(i)
        if (ic>0) then
          write(100,'(a)') ' '
          write(100,'(a)') 'TITLE ="grid"'
          write(100,'(a)') 'VARIABLES ="X", "Y"'
          write(100,'(a)') 'ZONE T="VOL_MIXED",N= 3, E=1, ET=TRIANGLE F=FEPOINT'
          do j=1,3
            write(100,*) cell_vert( cell2node(ic,j) ,1),cell_vert( cell2node(ic,j) ,2)
          enddo
          write(100,*) 1,2,3 !cell2node(ic,1),cell2node(ic,2),cell2node(ic,3)
        endif
      enddo
    enddo
    close(100)

    open(100,file='grid_vol.dat')
    do ic=1,num_cells
      write(100,*) ic,cell_area(ic)
    enddo
    close(100)

    open(100,file='grid_nodes.dat')
    write(100,'(a)') 'VARIABLES ="X", "Y"'
    do i=1,num_nodes
      write(100,*) cell_vert(i,1),cell_vert(i,2)
    enddo
    close(100)

    open(100,file='grid_edge_center.dat')
    write(100,'(a)') 'VARIABLES ="X", "Y"'
    do i=1,num_edges
      write(100,*) edge_center(i,1),edge_center(i,2)
    enddo
    close(100)

    open(100,file='grid_cell_center.dat')
    write(100,'(a)') 'VARIABLES ="X", "Y"'
    do i=1,num_cells
      write(100,*) cell_center(i,1),cell_center(i,2)
    enddo
    close(100)

    open(100,file='grid_cell_normals.dat')
    write(100,'(a)') 'VARIABLES ="X", "Y", "N1", "N2"'
    do ic=1,num_cells
      do ie=1,num_vert_cell(ic)
        ie1=cell2edge(ic,ie)
        xc=edge_normal(ie1,1) * edge_normal_sign(ic,ie)
        yc=edge_normal(ie1,2) * edge_normal_sign(ic,ie)
        write(100,*) cell_center(ic,1),cell_center(ic,2),xc,yc
      enddo
    enddo
    close(100)

    open(100,file='grid_con.dat')
    do i=1,num_cells
        write(100,*) (cell2node(i,j),j=1,3)
    enddo
    close(100)


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

    allocate(fv(num_nodes))


    call interpolate_cellcntr2node(fce,fv)

    fout='interp.plt'
    allocate(fw(num_nodes,1))
    fw(:,1)=fv(:)
    call test_tecplot(fout,1,fw)
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

    allocate(f(num_cells),fxv(num_nodes), fyv(num_nodes), df(num_cells,2))


!    call gradient_cellcntr('ggcb',fce,df)
!    call interpolate_cellcntr2node(df(1:num_cells,1),fxv)
!    call interpolate_cellcntr2node(df(1:num_cells,2),fyv)

    call gradient_cellcntr(fce,df)
    call interpolate_cellcntr2node(df(1:num_cells,1),fxv)
    call interpolate_cellcntr2node(df(1:num_cells,2),fyv)
    fout='gradient_lsq.plt'
    allocate(fw(num_nodes,2))
    fw(:,1)=fxv(:)
    fw(:,2)=fyv(:)
    call test_tecplot(fout,2,fw)


    call interpolate_cellcntr2node(dfce(1:num_cells,1),fxv)
    call interpolate_cellcntr2node(dfce(1:num_cells,2),fyv)
    fout='gradient_exact.plt'
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
    real,intent(in)     :: f(num_nodes,nvar)
    character(len=*)    :: fout

    integer             :: i,j

    open(100,file=trim(fout))
    write(100,'(a)') 'TITLE ="grid_sol"'
    if (nvar==1) then
      write(100,'(a)') 'VARIABLES ="x", "y", "f1"'
      write(100,'(a,i0,a,i0,a)') 'ZONE T="VOL_MIXED",N=', num_nodes, ' E=',num_cells, ' ET=TRIANGLE F=FEPOINT'
      do i=1,num_nodes
        write(100,'(3(e19.8,1x))') cell_vert(i,1),cell_vert(i,2), f(i,1)
      enddo
    elseif (nvar==2) then
      write(100,'(a)') 'VARIABLES ="x", "y", "f1", "f2"'
      write(100,'(a,i0,a,i0,a)') 'ZONE T="VOL_MIXED",N=', num_nodes, ' E=',num_cells, ' ET=TRIANGLE F=FEPOINT'
      do i=1,num_nodes
        write(100,'(4(e19.8,1x))') cell_vert(i,1),cell_vert(i,2), f(i,1),f(i,2)
      enddo
    else
      write(100,'(a)') 'VARIABLES ="x", "y", "f1", "U"'
      write(100,'(a,i0,a,i0,a)') 'ZONE T="VOL_MIXED",N=', num_nodes, ' E=',num_cells, ' ET=TRIANGLE F=FEPOINT'
      do i=1,num_nodes
        write(100,'(5(e19.8,1x))') cell_vert(i,1),cell_vert(i,2), f(i,1),f(i,2), f(i,3)
      enddo
    endif
    do i=1,num_cells
        write(100,*) (cell2node(i,j),j=1,3)
    enddo
    close(100)


    return
  end subroutine test_tecplot


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine test_analytic
    implicit none

    integer :: ic,in
    real    :: ax,ay, xc,yc, xmin,xmax, ymin,ymax

    xmin=minval(cell_vert(:,1))
    xmax=maxval(cell_vert(:,1))
    ax=10.d0*acos(-1.d0)/(xmax-xmin)

    ymin=minval(cell_vert(:,2))
    ymax=maxval(cell_vert(:,2))
    ay=10.d0*acos(-1.d0)/(ymax-ymin)

    do ic=1,num_cells
      xc=cell_center(ic,1)
      yc=cell_center(ic,2)
      fce(ic)=sin(ax*xc)*cos(ay*yc)
      dfce(ic,1)= ax*cos(ax*xc)*cos(ay*yc)
      dfce(ic,2)=-ay*sin(ax*xc)*sin(ay*yc)
    enddo

    do in=1,num_nodes
      xc=cell_vert(in,1)
      yc=cell_vert(in,2)
      fve(in)=sin(ax*xc)*cos(ay*yc)
      dfve(in,1)= ax*cos(ax*xc)*cos(ay*yc)
      dfve(in,2)=-ay*sin(ax*xc)*sin(ay*yc)
    enddo

    return
  end subroutine test_analytic


end module test


! do ic=1,num_cells
!   nt=0;
!   do in=1,num_vert_cell(ic)
!     iv=cell2node(ic,in)
!     nt = nt + node2cell_ntot(iv)
!   enddo
!   allocate(tmp(nt))
!
!   tmp(:)=0
!   nt=0
!   do in=1,num_vert_cell(ic)
!     iv=cell2node(ic,in)
!     do i=node2cell_ptr(iv), node2cell_ptr(iv) + node2cell_ntot(iv) - 1
!       ic1=node2cell(i)
!       nt=nt+1
!       tmp(nt)=ic1
!     enddo
!
!     i=0
!     min_val = minval(tmp)-1
!     max_val = maxval(tmp)
!     do while (min_val<max_val)
!         i = i+1
!         min_tmp = minval(tmp, mask=val>min_tmp)
!         unique(i) = min_tmp
!     enddo
!     allocate(final(i), source=unique(1:i))
!     tmp
!
!
!     call unique (tmp)
!
!
!     do iv=1,num_vert_cell(ic)
!       if (in==cell2node(ic,iv)) then
!         d=cell_dist2vert(ic,iv)
!         dt=dt+1.d0/d
!         interp_cellweight(i)=1.0/d
!       endif
!     enddo
!   enddo
!   interp_nodeweight(in)=1.d0/dt;
! enddo
