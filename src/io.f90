module io

  use mainparam
  use input
  use data_grid
  use data_solution
  !use tecplot
  use interpolation
  use gradient
  use ios_unstrc


  integer,save,private    :: counter=0, icont_cp=0
  logical,save,private    :: lcp=.false.

  logical,save            :: lw_inst(nvar)

  integer,parameter       :: imach_save = 2, imach_cont = 2
  integer,save            :: itape_inst, itape_save, itape_cont
  integer,save            :: iunit_inst, iunit_save, iunit_cont
  character(len=124),save :: fbase_inst, fbase_save, fbase_cont

contains


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine io_setup
    implicit none

    fbase_inst = 'inst'
    fbase_save = 'save'
    fbase_cont = 'cont'

    itape_inst = 11;    iunit_inst = itape_inst + 100;
    itape_save = 12;    iunit_save = itape_save + 100;
    itape_cont = 13;    iunit_cont = itape_cont + 100;

    lw_inst = .false.
    if (lw_rho) lw_inst(ir)=.true.
    if (lw_u  ) lw_inst(iu)=.true.
    if (lw_v  ) lw_inst(iv)=.true.
    if (lw_p  ) lw_inst(ip)=.true.

    return
  end subroutine io_setup


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine io_init
    implicit none

    integer               :: it,mt,mp, minf, ierror, ivar
    integer               :: itimes(maxtime)
    character(len=maxlen) :: inf(maxinf), varinfo(maxinf)
    character             :: dchart*6

    !--------------------------------------------------------------------------!
    ! transient files
    !--------------------------------------------------------------------------!
    itimes(1) = ntstart - 1 + nsubsteps(1)
    do it=2,nsaves
      itimes(it) = nsubsteps(it) + itimes(it-1)
    enddo

    mt = nsaves
    minf = 0

    inf(1) = 'rho'
    inf(2) = 'u'
    inf(3) = 'v'
    inf(4) = 'p'
    it = 0
    do ivar=1,nvar
      if (lw_inst(ivar)) then
        it = it + 1
        varinfo(it) = inf(ivar)
      endif
    enddo
    mp = it

    if (lw_inst_atall) then
      call writecd (itape_inst, iunit_inst, fbase_inst, imach_inst, &
                    mt, nnodes, ncells, mp, itimes, inf, minf, ierror)
    endif



    !--------------------------------------------------------------------------!
    ! final time-step: save file (it will be saved on cell centers)
    !--------------------------------------------------------------------------!
    mt = 1
    mp = 4

    itimes(1) = ntstart - 1 + ntimes
    inf(1) = 'rho'
    inf(2) = 'rhou'
    inf(3) = 'rhov'
    inf(4) = 'rhoE'

    write(dchart,'(i6)') ntimes
    minf = 4
    inf(5) = 'number of time-step computed = '//trim(dchart)
    inf(6) = 'conservatve variables are saved in cell centers'
    inf(7) = '#ncells and #nodes are replaced'
    inf(8) = ' '
    ! call writecd (itape_save, iunit_save, fbase_save, imach_save, &
    !               mt, nnodes, ncells, mp, itimes, inf, minf, ierror)
    call writecd (itape_save, iunit_save, fbase_save, imach_save, &
                  mt, ncells, nnodes, mp, itimes, inf, minf, ierror)

    return
  end subroutine io_init


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine write_inst_ios (sol_time,ntimes_sub)

    implicit none
    real,intent(in)     :: sol_time
    integer,intent(in)  :: ntimes_sub
    integer             :: ivar
    real,allocatable    :: f8(:)

    !-- precheck
    if (.not.lw_inst_atall) return
    if (lw_inst_atall) allocate(f8(nnodes))


    !-- convert conservative variables to primitive variables
    call cvar2pvar

    !-- write out transient ios files
    do ivar=1,nvar
      if (lw_inst(ivar)) then
        call interpolate_cell2node(pvar(ivar,1:ncells), f8)
        call writed(itape_inst, f8)
      endif
    enddo

    !-- deallocate
    deallocate(f8)

    return
  end subroutine write_inst_ios


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine write_save_ios

    implicit none
    integer             :: ivar
    real,allocatable    :: f8(:)

    !-- allocate array for node
    !allocate(f8(nnodes))

    !-- convert conservative variables to primitive variables
    !call cvar2pvar

    !-- write out all conservatve variables in cell centers
    do ivar=1,nvar
      !call interpolate_cell2node(pvar(ivar,1:ncells), f8)
      call writed(itape_save, cvar(ivar,1:ncells))
    enddo

    !-- deallocate
    !deallocate(f8)

    return
  end subroutine write_save_ios


  ! !============================================================================!
  ! !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  ! !============================================================================!
  ! subroutine write_inst_tec (sol_time,ntimes_sub)
  !   implicit none
  !
  !   real,intent(in)             :: sol_time
  !   integer,intent(in)          :: ntimes_sub
  !   integer                     :: ivar, np, i, ib
  !   real(kind=4),allocatable    :: f4(:,:),xy4(:,:)
  !   real,allocatable            :: f8(:)
  !   character(len=*)            :: file_out*100, varinfo*100, dchar5*5, dchar6*6
  !
  !
  !   !-- write out grid file only once
  !   if (counter == 0) then
  !     allocate(xy4(nnodes,2))
  !     do i=1,nnodes
  !       xy4(i,1) = node(i)%x
  !       xy4(i,2) = node(i)%y
  !     enddo
  !     file_out = trim('solution_grid.plt')
  !     call tecplot_write_grid (trim(file_out), xy4)
  !
  !     deallocate(xy4)
  !   endif
  !
  !
  !   !-- open file for cp
  !   if (counter==0) then
  !     do ib=1,nbndries
  !       if (trim(bndry(ib)%type)=='slip_wall' .or. trim(bndry(ib)%type)=='solid_wall') lcp=.true.
  !     enddo
  !     if (lcp) open(iunit_cp,file=trim(file_cp))
  !     if (lcp) then
  !       open(iunit_un,file=trim(file_un))
  !       write(iunit_un,'(a)') 'VARIABLES = "time" "|V<sub>n</sub>|<sub><math>%</math></sub>",   "|V<sub>n</sub>|<sub>2</sub>" , "|V<sub>n</sub>|<sub>1</sub>"'
  !     endif
  !   endif
  !
  !
  !   !-- determine time-step number for output file
  !   counter = counter + ntimes_sub
  !   write(dchar5,'(i5.5)') counter
  !   write(dchar6,'(i6.6)') counter
  !
  !
  !   !-- convert conservative variables to primitive variables
  !   call cvar2pvar
  !
  !
  !   !--
  !   if (lgrad_limiter) then
  !     np=nvar + 1
  !     varinfo  = 'rho u v p phi'
  !   else
  !     np=nvar
  !     varinfo  = 'rho u v p'
  !   endif
  !
  !
  !   !-- interpolate cell center solution to to cell nodes
  !   allocate(f8(nnodes), f4(nnodes,np))
  !   do ivar=1,nvar
  !     call interpolate_cell2node(pvar(ivar,1:ncells), f8)
  !     f4(1:nnodes,ivar) = f8(1:nnodes)
  !   enddo
  !   if (lgrad_limiter) then
  !     call interpolate_cell2node(phi_lim, f8)
  !     f4(1:nnodes,nvar+1) = f8(1:nnodes)
  !   endif
  !
  !
  !   !--write solution
  !   file_out = trim('solution_it'//dchar6//'.plt')
  !   !varinfo  = 'rho u v p'
  !   !np = 4
  !   call tecplot_write_solution (trim(file_out), trim(varinfo), np, f4, sol_time)
  !
  !
  !   !-- limiter
  !   ! allocate(g4(nnodes))
  !   ! call interpolate_cell2node(phi_lim, f8)
  !   ! g4(1:nnodes) = f8(1:nnodes)
  !   ! file_out = trim('limiter_it'//dchar6//'.plt')
  !   ! varinfo  = 'phi'
  !   ! np = 1
  !   ! call tecplot_write_solution (trim(file_out), trim(varinfo), np, g4, sol_time)
  !
  !
  !   !-- compute and write cp
  !   !if (lcp) call write_cp (sol_time)
  !
  !
  !   !--
  !   !call write_un (sol_time)
  !
  !   deallocate(f4,f8)
  !
  !   return
  ! end subroutine write_inst_tec
  !
  !
  ! !============================================================================!
  ! !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  ! !============================================================================!
  ! subroutine write_save_tec (sol_time)
  !   implicit none
  !
  !   real,intent(in)             :: sol_time
  !   integer                     :: ivar, np, i, ib, iprec
  !   real(kind=8),allocatable    :: xy(:,:), sol(:,:), f(:)
  !   character(len=*)            :: file_out*100, varinfo*100, dchar5*5, dchar6*6
  !
  !
  !   !-- 0:single / 1=double precision
  !   iprec = 1
  !
  !   !-- grid data
  !   allocate(xy(nnodes,2))
  !   do i=1,nnodes
  !     xy(i,1) = node(i)%x
  !     xy(i,2) = node(i)%y
  !   enddo
  !
  !
  !   !-- determine time-step number for output file
  !   write(dchar6,'(i6.6)') ntimes
  !
  !
  !   !-- convert conservative variables to primitive variables
  !   call cvar2pvar
  !
  !
  !   !-- interpolate cell center solution to to cell nodes
  !   allocate(f(nnodes), sol(nnodes,nvar))
  !   do ivar=1,nvar
  !     call interpolate_cell2node(pvar(ivar,1:ncells), f)
  !     sol(1:nnodes,ivar) = f(1:nnodes)
  !   enddo
  !
  !
  !   !--write grid solution
  !   file_out = trim('save.plt')
  !   varinfo  = 'x y rho u v p'
  !   np = 4
  !   call tecplot_write_grid_solution (trim(file_out), trim(varinfo), np, iprec, xy, sol, sol_time)
  !
  !
  !   deallocate(xy, sol, f)
  !
  !   return
  ! end subroutine write_save_tec



  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine write_inst_cp_un (sol_time)
    implicit none

    real,intent(in) :: sol_time
    integer         :: i,ib,ie,ic
    real            :: xf,yf,af,nx,ny, cp, cp1, pw, un,uw,vw, un_max, un_l1, un_l2
    real            :: ca,cn, ca1,cn1, cl,cd, cl1,cd1
    real,allocatable:: dp(:,:), du(:,:), dv(:,:)

    !-- open file for cp
    if (icont_cp==0) then
      do ib=1,nbndries
        if (trim(bndry(ib)%type)=='slip_wall' .or. trim(bndry(ib)%type)=='solid_wall') lcp=.true.
      enddo
      if (lcp) open(iunit_cp,file=trim(file_cp))
      if (lcp) then
        open(iunit_un,file=trim(file_un))
        write(iunit_un,'(a)') 'VARIABLES = "time" "|V<sub>n</sub>|<sub><math>%</math></sub>",   "|V<sub>n</sub>|<sub>2</sub>" , "|V<sub>n</sub>|<sub>1</sub>"'

        open(iunit_clcd,file=trim(file_clcd))
        write(iunit_clcd,'(a)') 'VARIABLES = "time" "c<sub>l</sub>",   "c<sub>d</sub>", "c<sub>l1</sub>",   "c<sub>d1</sub>"'
      endif
    endif


    !-- return if no slip_wall
    if (.not.lcp) return


    !-- counter
    icont_cp = icont_cp + 1


    !-- compute dp/dx and dp/dy
    allocate(dp(ncells,2))
    call gradient_cellcntr_1var(pvar(ip,1:ncells), dp)


    !-- compute d(u,v)/dx and d(u,v)/dy
    allocate(du(ncells,2), dv(ncells,2))
    call gradient_cellcntr_1var(pvar(iu,1:ncells), du)
    call gradient_cellcntr_1var(pvar(iv,1:ncells), dv)


    !-- write out cp
    do ib=1,nbndries
      if (trim(bndry(ib)%type)=='slip_wall' .or. trim(bndry(ib)%type)=='solid_wall') then
        write(iunit_cp,'(a)') 'TITLE     = "cp"'
        write(iunit_cp,'(a)') 'VARIABLES = "x" "cp_w" "cp_cell"'
        write(iunit_cp,'(a,i0,a)') 'ZONE I=',bndry(ib)%nedges,' J=1'
        write(iunit_cp,'(a,e16.8)') 'STRANDID=1, SOLUTIONTIME=',sol_time

        un_max = 0.d0
        un_l2  = 0.d0
        un_l1  = 0.d0

        cn = 0.d0
        ca = 0.d0
        cn1= 0.d0
        ca1= 0.d0

        do i=1,bndry(ib)%nedges
          ie = bndry(ib)%edge(i)
          ic = bndry(ib)%cell(i)

          xf = edge(ie)%x
          yf = edge(ie)%y
          af = edge(ie)%area
          nx = edge(ie)%nx
          ny = edge(ie)%ny

          !-- cp
          pw = pvar(ip,ic) + (edge(ie)%x-cell(ic)%x)*dp(ic,1) + (edge(ie)%y-cell(ic)%y)*dp(ic,2)
          cp = 2.d0/mach_inf**2 * (pw - p_inf)
          cp1= 2.d0/mach_inf**2 * (pvar(ip,ic) - p_inf)
          write(iunit_cp,'(3(e16.8,1x))') xf, cp, cp1


          !-- normal and axial force coefficents
          cn = cn + cp *ny*af
          ca = ca + cp *nx*af
          cn1= cn1+ cp1*ny*af
          ca1= ca1+ cp1*nx*af

          !-- Vn
          uw = pvar(iu,ic) + (edge(ie)%x-cell(ic)%x)*du(ic,1) + (edge(ie)%y-cell(ic)%y)*du(ic,2)
          vw = pvar(iv,ic) + (edge(ie)%x-cell(ic)%x)*dv(ic,1) + (edge(ie)%y-cell(ic)%y)*dv(ic,2)
          un = uw*nx + vw*ny
          un_l2 = un_l2 + un*un
          un_l1 = un_l1 + abs(un)
          un_max = max(un_max,abs(un))
        enddo

        un_l2 = sqrt(un_l2/dble(bndry(ib)%nedges))
        un_l1 = un_l1/dble(bndry(ib)%nedges)
        write(iunit_un,'(4(e16.8,1x))')    sol_time, un_max, un_l2, un_l1

        cl = cn *cosd(aoa_inf_deg) - ca *sind(aoa_inf_deg)
        cd = cn *sind(aoa_inf_deg) + ca *cosd(aoa_inf_deg)
        cl1= cn1*cosd(aoa_inf_deg) - ca1*sind(aoa_inf_deg)
        cd1= cn1*sind(aoa_inf_deg) + ca1*cosd(aoa_inf_deg)

        write(iunit_clcd,'(5(e16.8,1x))') sol_time, cl,cd, cl1,cd1
      endif
    enddo

    deallocate(dp, du, dv)

    return
  end subroutine write_inst_cp_un


  ! !============================================================================!
  ! !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  ! !============================================================================!
  ! subroutine write_un (sol_time)
  !   implicit none
  !
  !   real,intent(in) :: sol_time
  !   integer         :: i,ib,ie,ic
  !   real            :: xf,yf,af,nx,ny, un,uw,vw, un_max, un_l1, un_l2
  !   real,allocatable:: du(:,:), dv(:,:)
  !
  !
  !
  !   !-- compute dp/dx and dp/dy
  !   allocate(du(ncells,2), dv(ncells,2))
  !   call gradient_cellcntr_1var(pvar(iu,1:ncells), du)
  !   call gradient_cellcntr_1var(pvar(iv,1:ncells), dv)
  !
  !
  !   do ib=1,nbndries
  !     if (trim(bndry(ib)%type)=='slip_wall' .or. trim(bndry(ib)%type)=='solid_wall') then
  !       ! write(iunit_cp,'(a)') 'TITLE     = "un"'
  !       ! write(iunit_cp,'(a)') 'VARIABLES = "x" "u_n" "|u_n|"'
  !       ! write(iunit_cp,'(a,i0,a)') 'ZONE I=',bndry(ib)%nedges,' J=1'
  !       ! write(iunit_cp,'(a,e16.8)') 'STRANDID=1, SOLUTIONTIME=',sol_time
  !       un_max = 0.d0
  !       un_l2  = 0.d0
  !       un_l1  = 0.d0
  !       do i=1,bndry(ib)%nedges
  !         ie = bndry(ib)%edge(i)
  !         ic = bndry(ib)%cell(i)
  !
  !         xf = edge(ie)%x
  !         yf = edge(ie)%y
  !         af = edge(ie)%area
  !         nx = edge(ie)%nx
  !         ny = edge(ie)%ny
  !
  !         uw = pvar(iu,ic) + (edge(ie)%x-cell(ic)%x)*du(ic,1) + (edge(ie)%y-cell(ic)%y)*du(ic,2)
  !         vw = pvar(iv,ic) + (edge(ie)%x-cell(ic)%x)*dv(ic,1) + (edge(ie)%y-cell(ic)%y)*dv(ic,2)
  !
  !         un = uw*nx + vw*ny
  !
  !         un_l2 = un_l2 + un*un
  !         un_l1 = un_l1 + abs(un)
  !         un_max = max(un_max,abs(un))
  !       enddo
  !       un_l2 = sqrt(un_l2/dble(bndry(ib)%nedges))
  !       un_l1 = un_l1/dble(bndry(ib)%nedges)
  !       write(iunit_un,'(4(e16.8,1x))') sol_time, un_max, un_l2, un_l1
  !     endif
  !   enddo
  !
  !
  !   deallocate(du, dv)
  !
  !   return
  ! end subroutine write_un

end module io
