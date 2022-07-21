module input

  use mainparam
  use mpi

  implicit none

  !-- Grid paramteres
  character,save            :: file_grid*127, file_bc*127

  !-- Flow paramteres
  real,save                 :: rey        !-- Reynolds number
  real,save                 :: mach_inf   !-- Mach number
  real,save                 :: gamma
  real,save                 :: aoa_inf_deg!-- free-stream angle of attack (deg)

  !-- time integration parameters
  real,save                 :: dt                     !-- time-integration time-step
  integer,save              :: ntstart,ntimes,nsaves
  integer,allocatable,save  :: nsubsteps(:)
  real,allocatable,save     :: dtsave(:)              !-- time-step for writing out solution
  logical,save              :: lsteady, lwrite_resid
  real,save                 :: cfl_user

  !-- isentropic vortex
  logical,save              :: lvortex

  !-- Runge-Kutta time integration
  integer,save              :: rk_nstages !-- #s Runge-Kutta stages
  integer,save              :: rk_order   !-- order of accuracy of R-K
  logical,save              :: lSSPRK     !-- Strong stability preserving R-K

  !-- input parameters for gradient
  integer,save              :: grad_cellcntr_imethd
  logical,save              :: lgrad_ggcb, lgrad_ggnb, lgrad_ggnb_exp, lgrad_lsq
  logical,save              :: lgrad_lsq_fn, lgrad_lsq_nn
  real,save                 :: grad_cellcntr_lsq_pow

  !-- gradient limiter scheme
  integer,save              :: grad_limiter_imethd
  logical,save              :: lgrad_limiter
  character,save            :: limiter_type*20

  !-- face re-construction scheme
  integer,save              :: face_reconst_imethd
  logical,save              :: lface_reconst_upwind1st, lface_reconst_upwind2nd, lface_reconst_umuscl
  real,save                 :: umuscl_cst

  !-- inviscid flux discretization scheme
  integer,save              :: flux_inviscd_imethd
  logical,save              :: lflux_inviscd_roe


  !-- instantenous output files
  integer,save              :: imach_inst
  logical,save              :: lw_rho, lw_u, lw_v, lw_p, lw_inst_atall


contains

  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine input_read
    implicit none

    integer         :: istat,ierr,iloc,jloc,i
    logical         :: linputfile, lw_check
    character       :: message*124, dchar_grid*120, grad_cellcntr_lsq_nghbr*2, cmach_inst*2


    !--------------------------------------------------------------------------!
    ! user input file
    !--------------------------------------------------------------------------!
    inquire(file=trim(file_input),exist=linputfile);
    if (.not.linputfile.and.proc_id==0) then
      write(*,*)
      write(*,*) 'linputfile:',linputfile
      write(*,*) 'cannot find "'//trim(file_input)//'" file!'
      write(*,*) 'error in --> mod:input, sub:input_read'
      stop 'program stopped at "input_read"'
    endif

    open(iunit_input,file=trim(file_input),status='unknown',form='formatted')

    read(iunit_input,*)
    read(iunit_input,*)
    read(iunit_input,*)
    read(iunit_input,'(a)') dchar_grid

    read(iunit_input,*) rey
    read(iunit_input,*) mach_inf
    read(iunit_input,*) aoa_inf_deg
    read(iunit_input,*) gamma
    read(iunit_input,*)
    read(iunit_input,*) dt
    read(iunit_input,*) ntimes
    read(iunit_input,*) nsaves
    read(iunit_input,*) ntstart
    read(iunit_input,*) lsteady, cfl_user
    read(iunit_input,*) lvortex

    read(iunit_input,*)
    read(iunit_input,*) lw_rho, lw_u, lw_v, lw_p
    read(iunit_input,*) cmach_inst

    read(iunit_input,*)
    read(iunit_input,*)
    read(iunit_input,*)
    read(iunit_input,*) grad_cellcntr_imethd, grad_cellcntr_lsq_nghbr, grad_cellcntr_lsq_pow;
    read(iunit_input,*) grad_limiter_imethd
    read(iunit_input,*) face_reconst_imethd, umuscl_cst
    read(iunit_input,*) flux_inviscd_imethd

    read(iunit_input,*)
    read(iunit_input,*) rk_nstages
    read(iunit_input,*) rk_order
    read(iunit_input,*) lSSPRK

    close(iunit_input)


    !--
    lwrite_resid = .true.

    !dtsubstep = dt/dble(nsubsteps)
    allocate(dtsave(nsaves), nsubsteps(nsaves))
    if (mod(ntimes,nsaves)==0) then
      nsubsteps(1:nsaves) = ntimes/nsaves
         dtsave(1:nsaves) = dt*dble(nsubsteps(1))
    else
      nsubsteps(1:nsaves-1) = int(ntimes/nsaves) + 1
      nsubsteps(nsaves)     = ntimes - (int(ntimes/nsaves) + 1)*(nsaves-1)
         dtsave(1:nsaves-1) = dt*dble(nsubsteps(1))
         dtsave(nsaves)     = dt*dble(nsubsteps(nsaves))
    endif


    !-- task:
    if (ntstart==0) lvortex=.false.


    !--------------------------------------------------------------------------!
    ! Transient output files
    !--------------------------------------------------------------------------!
    if (lw_rho .or. lw_u .or. lw_v .or. lw_p) lw_inst_atall = .true.

    if (cmach_inst=='s4' .or. cmach_inst=='S4') then
      imach_inst=1
    elseif (cmach_inst=='s8' .or. cmach_inst=='S8') then
      imach_inst=2
    else
      if (proc_id==0) then
        write(*,*)
        write(*,*) 'format for regular output files in physical space:',cmach_inst
        write(*,*) 'format for regular output files must be either s4 or s8'
        write(*,*) 'error in --> mod:input, sr:input_read'
        stop 'error'
      endif
    endif


    !--------------------------------------------------------------------------!
    ! Grid and BC file names
    !--------------------------------------------------------------------------!
    loop1: do i=1,len(dchar_grid)
      if (dchar_grid(i:i)=='-') then
        iloc=i-1
        exit loop1
      endif
    enddo loop1
    jloc=len(trim(dchar_grid(1:iloc)))

    file_grid = dchar_grid(1:jloc)//grid_ext
    file_bc   = dchar_grid(1:jloc)//bc_ext


    !--------------------------------------------------------------------------!
    ! Gradient scheme to compute gradient at cell center
    !--------------------------------------------------------------------------!
    lgrad_ggcb=.false.

    lgrad_ggnb=.false.
    lgrad_ggnb_exp=.false.

    lgrad_lsq=.false.
    lgrad_lsq_fn =.false.
    lgrad_lsq_nn =.false.

    if (grad_cellcntr_imethd==1) then
      lgrad_ggcb=.true.

    elseif (grad_cellcntr_imethd==2) then
      lgrad_ggnb=.true.
      lgrad_ggnb_exp=.true.

    elseif (grad_cellcntr_imethd==3) then
      lgrad_lsq=.true.
      if (trim(grad_cellcntr_lsq_nghbr)=='fn' .or. trim(grad_cellcntr_lsq_nghbr)=='FN') then
        lgrad_lsq_fn=.true.

      elseif (trim(grad_cellcntr_lsq_nghbr)=='nn' .or. trim(grad_cellcntr_lsq_nghbr)=='NN') then
        lgrad_lsq_nn=.true.

      else
        write(*,*) 'check Least-Squares gradient scheme in input file'
        write(*,*) 'error at mod: input,  sub: input_read'
        stop
      endif

    else
      write(*,*) 'check cell-center gradient scheme in input file'
      write(*,*) 'error at mod: input,  sub: input_read'
      stop
    endif


    !--------------------------------------------------------------------------!
    ! Gradient limiter scheme
    !--------------------------------------------------------------------------!
    lgrad_limiter=.false.
    if (grad_limiter_imethd>0) then
      lgrad_limiter=.true.
      if (grad_limiter_imethd==1) then
        limiter_type = 'venk'

      elseif (grad_limiter_imethd==2) then
        limiter_type = 'barth'

      elseif (grad_limiter_imethd==3) then
        limiter_type = 'albada'

      else
        write(*,*) 'check gradient limiter scheme in input file'
        write(*,*) 'error at mod: input,  sub: input_read'
        stop
      endif
    endif


    !--------------------------------------------------------------------------!
    ! Face reconstruction scheme
    !--------------------------------------------------------------------------!
    lface_reconst_upwind1st=.false.
    lface_reconst_upwind2nd=.false.
    lface_reconst_umuscl=.false.

    if (face_reconst_imethd==1) then
      lface_reconst_upwind1st=.true.
      umuscl_cst = 0.d0

    elseif (face_reconst_imethd==2) then
      lface_reconst_upwind2nd=.true.
      umuscl_cst = 0.d0

    elseif (face_reconst_imethd==3) then
      lface_reconst_umuscl=.true.

    else
      write(*,*) 'check face reconstruction scheme in input file'
      write(*,*) 'error at mod: input,  sub: input_read'
      stop
    endif


    !--------------------------------------------------------------------------!
    ! Inviscid flux scheme
    !--------------------------------------------------------------------------!
    lflux_inviscd_roe=.false.
    if (flux_inviscd_imethd==1) then
      lflux_inviscd_roe=.true.

    else
      write(*,*) 'check inviscid flux discretization scheme in input file'
      write(*,*) 'error at mod: input,  sub: input_read'
      stop
    endif


    !--------------------------------------------------------------------------!
    ! Output file
    !--------------------------------------------------------------------------!
    if (proc_id==0) then
      open (iunit_log_input,file=trim(file_log_input),status='unknown',IOSTAT=istat)

      write(iunit_log_input,'(a)') '==========================================================================================================================================='
      write(iunit_log_input,'(a)') '     FVM2D CODE                       '
      write(iunit_log_input,'(a)') '==========================================================================================================================================='
      write(iunit_log_input,'(a35,a)') 'grid file name: ',trim(file_grid)
      write(iunit_log_input,'(a35,a)') 'bc file name: ',trim(file_bc)

      write(iunit_log_input,'(a35,E16.8)')    'Reynolds number: ',rey
      write(iunit_log_input,'(a35,E16.8)')    'Mach number: ',mach_inf
      write(iunit_log_input,'(a35,E16.8)')    'Flow angle (deg): ',aoa_inf_deg

      write(iunit_log_input,'(a)') '-------------------------------------------------------------------------------------------------------------------------------------------'

      write(iunit_log_input,'(a35,E16.8)')    'Time-step: ',dt
      write(iunit_log_input,'(a35,i0)')       '#s of total time-steps: ',ntimes
      write(iunit_log_input,'(a35,i0)')       '#s of output solution files: ',nsaves
      if (mod(ntimes,nsaves)==0) then
        write(iunit_log_input,'(a35,i0)')       'Interval to output solution files: ',int(ntimes/nsaves)
      else
        write(iunit_log_input,'(a35,i0,a3,i0)') 'Interval to output solution files: ',int(ntimes/nsaves)+1,' & ', ntimes-(int(ntimes/nsaves)+1)*(nsaves-1)
      endif
      if (lsteady) then
        write(iunit_log_input,'(a35,a)') 'Steady flow is computed: ','local time-stepping is employed (input dt is ignored)'
        write(iunit_log_input,'(a35,f4.2)') 'Local dt is computed based on CFL=',cfl_user
      else
        write(iunit_log_input,'(a35,i0)')       'Starting time-step: ',ntstart
        write(iunit_log_input,'(a35,e16.8)')    't_inital: ',dble(ntstart-1)*dt
        write(iunit_log_input,'(a35,e16.8)')    't_final : ',dble(ntstart-1)*dt+dble(ntimes)*dt
      endif

      if (ntstart<1) then
        write(iunit_log_input,'(a)') ' primative varialbes are initialized with manufactured solution'
      elseif (ntstart>1)   then
        write(iunit_log_input,'(a)') ' primative varialbes are initialized with continuation files'
        if (lvortex) then
          write(iunit_log_input,'(a)') ' primative varialbes are initialized with isentropic vortex'
          write(iunit_log_input,'(a)') ' Code will read freestream parameters from fvs2d.vortex file'
        endif
      elseif (ntstart==1) then
        if (lvortex) then
          write(iunit_log_input,'(a)') ' primative varialbes are initialized with isentropic vortex'
          write(iunit_log_input,'(a)') ' Code will read freestream parameters from fvs2d.vortex file'
        else
          write(iunit_log_input,'(a)') ' primative varialbes are initialized with freestream values'
        endif
      endif

      write(iunit_log_input,'(a)') '-------------------------------------------------------------------------------------------------------------------------------------------'
      if (lw_inst_atall .and. imach_inst==1) write(iunit_log_input,'(a52)') ' Write out following variables in single precision:'
      if (lw_inst_atall .and. imach_inst==2) write(iunit_log_input,'(a52)') ' Write out following variables in double precision:'
      if (lw_rho) write(iunit_log_input,'(a52)') ' density'
      if (lw_u  ) write(iunit_log_input,'(a52)') ' u-velocity'
      if (lw_v  ) write(iunit_log_input,'(a52)') ' v-velcoity'
      if (lw_p  ) write(iunit_log_input,'(a52)') ' pressure'

      write(iunit_log_input,'(a)') ' '
      write(iunit_log_input,'(a)') '==========================================================================================================================================='
      write(iunit_log_input,'(a)') '     Temporal & Spatial Discretization Schemes      '
      write(iunit_log_input,'(a)') '==========================================================================================================================================='

      !-- Cell-center gradient scheme
      if (lgrad_ggcb) then
        write(iunit_log_input,'(a38,a)') ' Cell-center gradient method:',' Green-Gauss Cell-Base'
      elseif (lgrad_ggnb) then
        write(iunit_log_input,'(a38,a)') ' Cell-center gradient method:',' Green-Gauss Node-Base'
      elseif (lgrad_lsq) then
        if (grad_cellcntr_lsq_pow==0.d0) then
          if (lgrad_lsq_fn) write(iunit_log_input,'(a38,a)') ' Cell-center gradient method:',' Unweigghted Least-Squeres based on face neighbor stencil'
          if (lgrad_lsq_nn) write(iunit_log_input,'(a38,a)') ' Cell-center gradient method:',' Unweigghted Least-Squeres based on node neighbor stencil'
        else
          if (lgrad_lsq_fn) write(iunit_log_input,'(a38,a,f3.1,a)') ' Cell-center gradient method:',' Weigghted (1/d^',grad_cellcntr_lsq_pow,') Least-Squeres based on face neighbor stencil'
          if (lgrad_lsq_nn) write(iunit_log_input,'(a38,a,f3.1,a)') ' Cell-center gradient method:',' Weigghted (1/d^',grad_cellcntr_lsq_pow,') Least-Squeres based on node neighbor stencil'
        endif
      endif

      !-- Gradient limiter
      if (.not.lgrad_limiter) then
        write(iunit_log_input,'(a38,a)') 'Gradient limiter:',' not applied'
      elseif (lgrad_limiter) then
        select case (limiter_type)
        case ('venk')
          write(iunit_log_input,'(a38,a)') 'Gradient limiter:',' Venkatakrishnan'
        case ('barth')
          write(iunit_log_input,'(a38,a)') 'Gradient limiter:',' Barth and Jespersen'
        case ('albada')
          write(iunit_log_input,'(a38,a)') 'Gradient limiter:',' Van Albada'
        case default
          write(iunit_log_input,'(a)')     'Error in limiter type'
          stop 'error!'
        end select
      endif

      !-- Face reconstruction
      if (lface_reconst_upwind1st) then
        write(iunit_log_input,'(a38,a)') ' Face reconstruction method:',' 1st-order upwind'
      elseif (lface_reconst_upwind2nd) then
        write(iunit_log_input,'(a38,a)') ' Face reconstruction method:',' 2nd-order upwind'
      elseif (lface_reconst_umuscl) then
        write(iunit_log_input,'(a38,a,f8.4)') ' Face reconstruction method:',' UMUSCL with Kappa=',umuscl_cst
      endif

      !-- Inviscid flux discretization
      if (lflux_inviscd_roe) then
        write(iunit_log_input,'(a38,a)') ' Inviscid flux discretization scheme:',' Roe'
      endif

      !--Runge-Kutta time-integration
      write(iunit_log_input,'(a)') '-------------------------------------------------------------------------------------------------------------------------------------------'
      if (lSSPRK) then
        write(iunit_log_input,'(a50)') 'Runge-Kutta SSP formualtion is employed'
      else
        write(iunit_log_input,'(a50)') 'Runge-Kutta standard formualtion is employed'
      endif
      write(iunit_log_input,'(a52,i0)') '#s of stages for Runge-Kutta time-integration: ',rk_nstages
      write(iunit_log_input,'(a52,i0)') 'Order of accuracy of Runge-Kutta time-integration: ',rk_order


      write(iunit_log_input,'(a)') '-------------------------------------------------------------------------------------------------------------------------------------------'
      write(iunit_log_input,'(a)') ' Transient output files         '




      close(iunit_log_input)
    endif

    return
  end subroutine input_read

end module input
