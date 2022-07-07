module test

  use mainparam,  only  : nvar
  use input
  use data_grid
  use data_solution,  only  : resid
  use grid_procs
  use interpolation
  use gradient
  use residual, only  : compute_residual
  use mms
  use tecplot

  implicit none

  real,allocatable,save :: fve(:),dfve(:,:), fce(:),dfce(:,:)

contains

  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine test_init
    implicit none

    allocate(fve(nnodes), fce(ncells), dfve(nnodes,2), dfce(ncells,2))

    !call verify_vortex
    !return

    call test_analytic
    !
    !call test_tec
    !
    !call test_grid
    !
    !call test_interpolation
    !
    !call test_gradient
    call test_resid

    return
  end subroutine test_init


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine verify_vortex

    use data_solution , only : pvar, cvar
    implicit none

    integer             :: ivar, np, i
    real(kind=4)        :: f(nnodes,nvar), xy(nnodes,2)
    real(kind=8)        :: sol_time, f8(nnodes)
    character(len=*)    :: file_out*100, varinfo*100

    do ivar=1,nvar
      call interpolate_cell2node(cvar(ivar,1:ncells), f8)
      f(1:nnodes,ivar) = f8(1:nnodes)
    enddo

    !--
    do i=1,nnodes
      xy(i,1) = node(i)%x
      xy(i,2) = node(i)%y
    enddo
    file_out = trim('vortex_grid.plt')
    call tecplot_write_grid (trim(file_out), xy)

    !--
    file_out = trim('vortex_solution.plt')
    varinfo  = 'rho u v p'
    np = 4
    sol_time = 1.d0
    call tecplot_write_solution (trim(file_out), trim(varinfo), np, f, sol_time)

    return
  end subroutine verify_vortex


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine test_tec
    implicit none

    integer             :: nedge_max, i,ic,jc, np,j
    logical             :: lGridSolTogether, lSZL
    real(kind=4)        :: xy(nnodes,2), sol(nnodes)
    character           :: file_out*100, varinfo*100
    real                :: sol_time
    character(len=80)    :: new_file1, new_file2

    ! np = 1
    ! varinfo = trim('f')
    ! lGridSolTogether = .false.
    ! lSZL=.false.
    !
    ! !call test_resid
    ! nedge_max=3
    ! if (ncells_quad>0) nedge_max=4
    !
    ! allocate(cell2node_tec(nedge_max,ncells))
    !
    ! do ic=1,ncells_tri
    !   do i=1,cell(ic)%nvrt
    !     cell2node_tec(i,ic) = cell(ic)%node(i)
    !   enddo
    !   cell2node_tec(nedge_max,ic) = cell(ic)%node(3)
    ! enddo
    ! do jc=1,ncells_quad
    !   ic=jc + ncells_tri
    !   do i=1,cell(ic)%nvrt
    !     cell2node_tec(i,ic) = cell(ic)%node(i)
    !   enddo
    ! enddo
    !
    ! call ios2tec_init (nnodes, ncells, nedges, nedge_max, cell2node_tec, np, varinfo, lSZL, lGridSolTogether)

    do i=1,nnodes
      xy(i,1) = node(i)%x
      xy(i,2) = node(i)%y
    enddo

    sol = fve

    file_out = trim('w_grid_solution.plt')

    new_file1 = trim(adjustl(file_out))
    call StripSpaces (new_file1)
    j=len(trim(new_file1))
    new_file2 = new_file1(1:j)

    varinfo  = 'x, y, f'
    np = 1
    sol_time = 1.d0
    call tecplot_write_grid_solution (new_file2, trim(varinfo), np, xy, sol, sol_time)


    file_out = trim('w_grid.plt')
    call tecplot_write_grid (trim(file_out), xy)

    file_out = trim('w_solution.plt')
    varinfo  = 'f'
    np = 1
    sol_time = 1.d0
    call tecplot_write_solution (trim(file_out), trim(varinfo), np, sol, sol_time)






    return
  end subroutine test_tec


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
        !if (cell(ic)%nvrt==4) write(*,*) ic, edge(je)%nx , edge(je)%ny
      enddo
    enddo
    close(100)

    do ic=1,ncells
      do ie=1,cell(ic)%nvrt
        je=cell(ic)%edge(ie)
        if (je<=0) write(*,'(a,4(i0,1x))') 'wow ',ic,cell(ic)%nvrt,ie,je
      enddo
    enddo

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


    open(100,file='grid_interior.plt')
    write(100,*) 'variables = "X" "Y"'
    do i=1,ncells_intr
      ic=cell_intr(i)
      write(100,*) cell(ic)%x,cell(ic)%y
    enddo
    close(100)

    open(100,file='grid_boundarry.plt')
    write(100,*) 'variables = "X" "Y"'
    do i=1,ncells_bndr
      ic=cell_bndr(i)
      write(100,*) cell(ic)%x,cell(ic)%y
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

    fout='log_interp.plt'
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
    real                :: pi,cr0,crs,crx,cry
    real,allocatable    :: f(:),fx(:),fy(:),errx(:),erry(:), df(:,:), fw(:,:)
    real                :: vol, r(2), d,dt,xc,yc, nxf,nyf,af, ax,ay, at, heff
    character           :: fout*124

    allocate(df(ncells,2))


    !--------------------------------------------------------------------------!
    !  manufactured solution and its gradients
    !--------------------------------------------------------------------------!
    allocate( f(ncells), fx(ncells), fy(ncells), errx(ncells), erry(ncells) )
    pi  = acos(-1.d0)
    cr0 =  1.12
    crs =  0.15
    crx =  3.12*pi
    cry =  2.92*pi

    do ic=1,ncells
      xc=cell(ic)%x
      yc=cell(ic)%y
      f (ic) = manufactured_sol(cr0,crs,crx,cry, 0,0, xc,yc)
      fx(ic) = manufactured_sol(cr0,crs,crx,cry, 1,0, xc,yc)
      fy(ic) = manufactured_sol(cr0,crs,crx,cry, 0,1, xc,yc)
    enddo

    !-- compute numerical gradient
    call gradient_cellcntr_test(fce,df)

    !-- compute effective distance, error
    errx = 0.d0
    erry = 0.d0
    do i=1,ncells_intr
      ic = cell_intr(i)
      at = at + cell(ic)%vol
      errx(ic) = abs( dfce(ic,1)-df(ic,1) ) !abs( fx(ic)-df(ic,1) )
      erry(ic) = abs( dfce(ic,2)-df(ic,2) ) !abs( fy(ic)-df(ic,2) )
    enddo
    heff = sqrt(at/dble(ncells_intr))

    write(*,'(a)')  '       heff           err_max           err_l2           err_max           err_l2'
    write(*,'(5(e16.9,1x))') heff,maxval(errx), dsqrt(sum(errx**2))/dble(ncells_intr),  maxval(erry), dsqrt(sum(erry**2))/dble(ncells_intr)
    !write(*,'(3(e16.9,1x))') heff,maxval(erry), dsqrt(sum(erry**2))/dble(ncells_intr)

    !write(*,*) cell(maxloc(errx))%x, cell(maxloc(errx))%y


    allocate(fw(nnodes,2))
    call interpolate_cell2node(df  (1:ncells,1), fw(:,1))
    call interpolate_cell2node(dfce(1:ncells,1), fw(:,2))
    fout='test_gradientx.plt'
    if (ncells_quad>0) then
      call test_tecplot_mixed (fout,2,fw)
    else
      call test_tecplot (fout,2,fw)
    endif
    !
    !
    ! call interpolate_cell2node(dfce(1:ncells,1),fxv)
    ! call interpolate_cell2node(dfce(1:ncells,2),fyv)
    ! fout='test_gradient_exact.plt'
    ! fw(:,1)=fxv(:)
    ! fw(:,2)=fyv(:)
    ! call test_tecplot(fout,2,fw)
    ! deallocate(fw)


    return
  end subroutine test_gradient


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine test_resid

    implicit none

    real,allocatable  :: wrk1(:), wrk2(:,:), wrk3(:)
    real              :: l2(4), linf(4), heff, at
    integer           :: i,ic,ivar
    character         :: fout*124
    logical           :: lexist


    allocate( wrk1(nnodes), wrk3(nnodes), wrk2(nnodes,2))

    call compute_residual


    !-- compute error
    do ivar=1,4
      l2(ivar) = 0.d0
      linf(ivar) = 0.d0
      do i=1,ncells_intr
        ic = cell_intr(i)

        l2(ivar) = l2(ivar) + ( resid(ivar,ic)+mms_source(ivar,ic) )**2
        linf(ivar) = max( linf(ivar), abs( resid(ivar,ic)+mms_source(ivar,ic) ) )
      enddo
      l2(ivar) = dsqrt(l2(ivar))/dble(ncells_intr)
    enddo

    fout='error_resid.plt'
    inquire(file=trim(fout), exist=lexist)
    if (lexist) then
      open(100, file=trim(fout), status="old", position="append", action="write")
    else
      open(100, file=trim(fout), status="new", action="write")
      write(100,'(a)')  'variables = "h<sub>eff" "L<sub>2,rho" "L<sub>2,u"   "L<sub>2,v"  "L<sub>2,e"  "L<sub>inf,rho" "L<sub>inf,u"   "L<sub>inf,v"  "L<sub>inf,e"    '
    end if

    write(100,'(9(e16.9,1x))') heff,(l2(ivar),ivar=1,4),(linf(ivar),ivar=1,4)

    ! write(*,'(a)')  '       heff           err_max           err_l2           err_max           err_l2'
    ! write(*,'(5(e16.9,1x))') heff,maxval(errx), dsqrt(sum(errx**2))/dble(ncells_intr),  maxval(erry), dsqrt(sum(erry**2))/dble(ncells_intr)


    call interpolate_cell2node(resid(1:1,1:ncells),wrk1)
    call interpolate_cell2node(mms_source(1:1,1:ncells),wrk3)
    wrk2(:,1)=abs(wrk1(:)-wrk3(:))
    wrk2(:,2)=wrk3(:)
    fout='test_resid_rho.plt'
    if (ncells_quad==0) call test_tecplot (fout,2,wrk2)
    if (ncells_quad>0)  call test_tecplot_mixed(fout,2,wrk2)


    call interpolate_cell2node(resid(2:2,1:ncells),wrk1)
    call interpolate_cell2node(mms_source(2:2,1:ncells),wrk3)
    wrk2(:,1)=abs(wrk1(:)-wrk3(:))
    wrk2(:,2)=wrk3(:)
    fout='test_resid_u.plt'
    if (ncells_quad==0) call test_tecplot (fout,2,wrk2)
    if (ncells_quad>0)  call test_tecplot_mixed(fout,2,wrk2)

    call interpolate_cell2node(resid(3:3,1:ncells),wrk1)
    call interpolate_cell2node(mms_source(3:3,1:ncells),wrk3)
    wrk2(:,1)=abs(wrk1(:)-wrk3(:))
    wrk2(:,2)=wrk3(:)
    fout='test_resid_v.plt'
    if (ncells_quad==0) call test_tecplot (fout,2,wrk2)
    if (ncells_quad>0)  call test_tecplot_mixed(fout,2,wrk2)

    call interpolate_cell2node(resid(4:4,1:ncells),wrk1)
    call interpolate_cell2node(mms_source(4:4,1:ncells),wrk3)
    wrk2(:,1)=abs(wrk1(:)-wrk3(:))
    wrk2(:,2)=wrk3(:)
    fout='test_resid_e.plt'
    if (ncells_quad==0) call test_tecplot (fout,2,wrk2)
    if (ncells_quad>0)  call test_tecplot_mixed(fout,2,wrk2)

    return
  end subroutine test_resid


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

    ! open(100,file='test_grid.plt')
    ! write(100,'(a)') 'TITLE = "grid"'
    ! write(100,'(a)') 'VARIABLES ="x", "y"'
    ! write(100,'(a)') 'ZONE T="zone 1"'
    ! write(100,'(a)') 'STRANDID=1, SOLUTIONTIME=0'
    ! write(100,'(a,i0,a,i0,a)') 'NODES=',nnodes, ' ELEMENTS=',ncells, ' DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL'
    ! do i=1,nnodes
    !   write(100,'(2(e19.8,1x))') node(i)%x,node(i)%y
    ! enddo
    !
    ! do i=1,ncells_tri
    !   write(100,*) cell(i)%node(1),cell(i)%node(2),cell(i)%node(3),cell(i)%node(3)
    ! enddo
    ! do i=1,ncells_quad
    !   k=i+ncells_tri
    !   write(100,*) (cell(k)%node(j),j=1,4)
    ! enddo
    ! close(100)
    !
    ! open(100,file='test_sol.plt')
    ! write(100,'(a)') 'TITLE = "grid | solution"'
    ! write(100,'(a)') 'VARIABLES ="f"'
    ! write(100,'(a)') 'ZONE T="zone 1"'
    ! write(100,'(a)') 'STRANDID=1, SOLUTIONTIME=1'
    ! write(100,'(a,i0,a,i0,a)') 'NODES=',nnodes, ' ELEMENTS=',ncells, ' DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL'
    ! do i=1,nnodes
    !   write(100,'(1(e19.8,1x))') f(i,1)
    ! enddo
    ! close(100)
    ! return


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
    real    :: ax,ay, xc,yc, xmin,xmax, ymin,ymax, pi

    pi = acos(-1.d0)

    xmin=minval(node(:)%x)
    xmax=maxval(node(:)%x)
    ax=10.d0*pi/(xmax-xmin)

    ymin=minval(node(:)%y)
    ymax=maxval(node(:)%y)
    ay=10.d0*pi/(ymax-ymin)

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


  subroutine StripSpaces(string)
     character(len=*) :: string
     integer :: stringLen
     integer :: last, actual

     stringLen = len (string)
     last = 1
     actual = 1

     do while (actual < stringLen)
         if (string(last:last) == ' ') then
             actual = actual + 1
             string(last:last) = string(actual:actual)
             string(actual:actual) = ' '
         else
             last = last + 1
             if (actual < last) &
                 actual = last
         endif
     end do

     end subroutine

end module test
