module io

  use mainparam
  use input,        only  : dtsave, mach_inf
  use data_grid
  use data_solution
  use tecplot
  use interpolation
  use gradient


  integer,save,private  :: counter=0
  logical,save,private  :: lcp=.false.


contains

  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine io_write_inst (sol_time,ntimes_sub)
    implicit none

    real,intent(in)             :: sol_time
    integer,intent(in)          :: ntimes_sub
    integer                     :: ivar, np, i, ib
    real(kind=4),allocatable    :: f4(:,:),xy4(:,:)
    real,allocatable            :: f8(:)
    character(len=*)            :: file_out*100, varinfo*100, dchar5*5, dchar6*6


    !-- write out grid file only once
    if (counter == 0) then
      allocate(xy4(nnodes,2))
      do i=1,nnodes
        xy4(i,1) = node(i)%x
        xy4(i,2) = node(i)%y
      enddo
      file_out = trim('solution_grid.plt')
      call tecplot_write_grid (trim(file_out), xy4)

      deallocate(xy4)
    endif


    !-- open file for cp
    if (counter==0) then
      do ib=1,nbndries
        if (trim(bndry(ib)%type)=='slip_wall' .or. trim(bndry(ib)%type)=='solid_wall') lcp=.true.
      enddo
      if (lcp) open(iunit_cp,file=trim(file_cp))
    endif


    !-- determine time-step number for output file
    counter = counter + ntimes_sub
    write(dchar5,'(i5.5)') counter
    write(dchar6,'(i6.6)') counter


    !-- convert conservative variables to primitive variables
    call cvar2pvar


    !--
    if (lgrad_limiter) then
      np=nvar + 1
      varinfo  = 'rho u v p phi'
    else
      np=nvar
      varinfo  = 'rho u v p'
    endif


    !-- interpolate cell center solution to to cell nodes
    allocate(f8(nnodes), f4(nnodes,np))
    do ivar=1,nvar
      call interpolate_cell2node(pvar(ivar,1:ncells), f8)
      f4(1:nnodes,ivar) = f8(1:nnodes)
    enddo
    if (lgrad_limiter) then
      call interpolate_cell2node(phi_lim, f8)
      f4(1:nnodes,nvar+1) = f8(1:nnodes)
    endif


    !--write solution
    file_out = trim('solution_it'//dchar6//'.plt')
    !varinfo  = 'rho u v p'
    !np = 4
    call tecplot_write_solution (trim(file_out), trim(varinfo), np, f4, sol_time)


    !-- limiter
    ! allocate(g4(nnodes))
    ! call interpolate_cell2node(phi_lim, f8)
    ! g4(1:nnodes) = f8(1:nnodes)
    ! file_out = trim('limiter_it'//dchar6//'.plt')
    ! varinfo  = 'phi'
    ! np = 1
    ! call tecplot_write_solution (trim(file_out), trim(varinfo), np, g4, sol_time)


    !-- compute and write cp
    if (lcp) call write_cp (sol_time)

    deallocate(f4,f8)

    return
  end subroutine io_write_inst


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine io_write_save (sol_time)
    implicit none

    real,intent(in)             :: sol_time
    integer                     :: ivar, np, i, ib
    real(kind=8),allocatable    :: f(:,:),g(:)
    real(kind=4),allocatable    :: f4(:,:)
    character(len=*)            :: file_out*100, varinfo*100, dchar5*5, dchar6*6


    !-- write out grid file only once
    allocate(f4(nnodes,2))
    do i=1,nnodes
      f4(i,1) = node(i)%x
      f4(i,2) = node(i)%y
    enddo
    file_out = trim('save_grid.plt')
    call tecplot_write_grid (trim(file_out), f4)
    deallocate(f4)


    !-- determine time-step number for output file
    write(dchar6,'(i6.6)') ntimes


    !-- convert conservative variables to primitive variables
    call cvar2pvar


    !-- interpolate cell center solution to to cell nodes
    allocate(g(nnodes), f(nnodes,nvar))
    do ivar=1,nvar
      call interpolate_cell2node(pvar(ivar,1:ncells), g)
      f(1:nnodes,ivar) = g(1:nnodes)
    enddo


    !--write solution
    file_out = trim('save.plt')
    varinfo  = 'rho u v p'
    np = 4
    call tecplot_write_solution8 (trim(file_out), trim(varinfo), np, f, sol_time)


    deallocate(f,g)

    return
  end subroutine io_write_save



  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine write_cp (sol_time)
    implicit none

    real,intent(in) :: sol_time
    integer         :: i,ib,ie,ic
    real            :: xf,yf,af,nx,ny, cp, pw
    real,allocatable:: dp(:,:)



    !-- compute dp/dx and dp/dy
    allocate(dp(ncells,2))
    call gradient_cellcntr_1var(pvar(ip,1:ncells), dp)

    do ib=1,nbndries
      if (trim(bndry(ib)%type)=='slip_wall' .or. trim(bndry(ib)%type)=='slip_wall') then
        write(iunit_cp,'(a)') 'TITLE     = "cp"'
        write(iunit_cp,'(a)') 'VARIABLES = "x" "cp_w" "cp_cell"'
        write(iunit_cp,'(a,i0,a)') 'ZONE I=',bndry(ib)%nedges,' J=1'
        write(iunit_cp,'(a,e16.8)') 'STRANDID=1, SOLUTIONTIME=',sol_time
        do i=1,bndry(ib)%nedges
          ie = bndry(ib)%edge(i)
          ic = bndry(ib)%cell(i)

          xf = edge(ie)%x
          yf = edge(ie)%y
          af = edge(ie)%area
          nx = edge(ie)%nx
          ny = edge(ie)%ny

          pw = pvar(ip,ic) + (edge(ie)%x-cell(ic)%x)*dp(ic,1) + (edge(ie)%y-cell(ic)%y)*dp(ic,2)
          cp = 2.d0/mach_inf**2 * (pw - p_inf)

          write(iunit_cp,'(3(e16.8,1x))') xf,cp, 2.d0/mach_inf**2 * (pvar(ip,ic) - p_inf) !yf,pvar(ip,ic),pvar(iu,ic),pvar(iv,ic),pvar(iu,ic)*nx + pvar(iv,ic)*ny
        enddo
      endif
    enddo

    deallocate(dp)

    return
  end subroutine write_cp

end module io
