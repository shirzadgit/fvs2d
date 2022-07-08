module input

  use mainparam
  use mpi

  implicit none

  !-- Grid paramteres
  character,save    :: file_grid*127, file_bc*127

  !-- Flow paramteres
  real,save         :: rey        !-- Reynolds number
  real,save         :: mach_inf   !-- Mach number
  real,save         :: gamma
  real,save         :: aoa_inf_deg!-- free-stream angle of attack (deg)
  real,save         :: dt         !-- time-step for writing out
  real,save         :: dtsubstep  !-- time-integration time-step
  integer,save      :: ntstart,ntimes,nsubsteps,itsave,ntsave

  !-- isentropic vortex
  logical,save      :: lvortex


  !-- Runge-Kutta time integration
  integer,save      :: rk_nstages !-- #s Runge-Kutta stages
  integer,save      :: rk_order   !-- order of accuracy of R-K
  logical,save      :: lSSPRK     !-- Strong stability preserving R-K

  !-- input parameters for gradient
  integer,save      :: grad_cellcntr_imethd
  logical,save      :: lgrad_ggcb, lgrad_ggnb, lgrad_ggnb_exp, lgrad_lsq
  logical,save      :: lgrad_lsq_fn, lgrad_lsq_nn
  real,save         :: grad_cellcntr_lsq_pow

  !-- gradient limiter scheme
  integer,save      :: grad_limiter_imethd
  logical,save      :: lgrad_limiter

  !-- face re-construction scheme
  integer,save      :: face_reconst_imethd
  logical,save      :: lface_reconst_linear, lface_reconst_umuscl

  !-- inviscid flux discretization scheme
  integer,save      :: flux_inviscd_imethd
  logical,save      :: lflux_inviscd_roe

contains

  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine input_read
    implicit none

    integer         :: istat,ierr,iloc,jloc,i
    logical         :: linputfile
    character       :: message*124, dchar_grid*120, grad_cellcntr_lsq_nghbr*2


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
    read(iunit_input,*) nsubsteps
    read(iunit_input,*) itsave
    read(iunit_input,*) ntstart
    read(iunit_input,*) lvortex

    read(iunit_input,*)
    read(iunit_input,*)
    read(iunit_input,*)
    read(iunit_input,*) grad_cellcntr_imethd, grad_cellcntr_lsq_nghbr, grad_cellcntr_lsq_pow;
    read(iunit_input,*) grad_limiter_imethd
    read(iunit_input,*) face_reconst_imethd
    read(iunit_input,*) flux_inviscd_imethd

    read(iunit_input,*)
    read(iunit_input,*) rk_nstages
    read(iunit_input,*) rk_order
    read(iunit_input,*) lSSPRK

    close(iunit_input)




    dtsubstep = dt/dble(nsubsteps)
    if (ntstart==0) lvortex=.false.

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
      write(*,*) 'scheme not implemented yet!'
      stop 'scheme not implemented yet!'
      ! if (grad_limiter_imethd==1) then
      !
      ! elseif
      !
      !
      ! else
      !   write(*,*) 'check gradient limiter scheme in input file'
      !   write(*,*) 'error at mod: input,  sub: input_read'
      !   stop
      ! endif
    endif


    !--------------------------------------------------------------------------!
    ! Face reconstruction scheme
    !--------------------------------------------------------------------------!
    lface_reconst_linear=.false.
    lface_reconst_umuscl=.false.
    if (face_reconst_imethd==1) then
      lface_reconst_linear=.true.

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
      write(iunit_log_input,'(a35,i0)')       '#s of sub-steps per time-steps: ',nsubsteps
      write(iunit_log_input,'(a35,i0)')       'Interval to output save files: ',itsave
      write(iunit_log_input,'(a35,i0)')       'Starting time-step: ',ntstart
      write(iunit_log_input,'(a35,e16.8)')    't_inital: ',dble(ntstart-1)*dt
      write(iunit_log_input,'(a35,e16.8)')    't_final : ',dble(ntstart-1)*dt+dble(ntimes)*dt

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
          if (lgrad_lsq_fn) write(iunit_log_input,'(a38,a)') ' Cell-center gradient method:',' Unweigghted Least-Squeres based face neighbor stencil'
          if (lgrad_lsq_nn) write(iunit_log_input,'(a38,a)') ' Cell-center gradient method:',' Unweigghted Least-Squeres based node neighbor stencil'
        else
          if (lgrad_lsq_fn) write(iunit_log_input,'(a38,a,f3.1,a)') ' Cell-center gradient method:',' Weigghted (1/d^',grad_cellcntr_lsq_pow,') Least-Squeres based face neighbor stencil'
          if (lgrad_lsq_nn) write(iunit_log_input,'(a38,a,f3.1,a)') ' Cell-center gradient method:',' Weigghted (1/d^',grad_cellcntr_lsq_pow,') Least-Squeres based node neighbor stencil'
        endif
      endif

      !-- Gradient limiter
      if (.not.lgrad_limiter) then
        write(iunit_log_input,'(a38,a)') 'Gradient limiter:',' not applied'
      elseif (lgrad_limiter) then

      endif

      !-- Face reconstruction
      if (lface_reconst_linear) then
        write(iunit_log_input,'(a38,a)') ' Face reconstruction method:',' Linear extrapolation'
      elseif (lface_reconst_umuscl) then
        write(iunit_log_input,'(a38,a)') ' Face reconstruction method:',' UMUSCL'
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


      close(iunit_log_input)
    endif

    return
  end subroutine input_read

end module input
