module input

  use mainparam
  use mpi

  implicit none

  !-- Grid paramteres
  character,save          :: file_grid*127

  !-- input parameters for gradient
  integer,save      :: grad_cellcntr_imethd1, grad_cellcntr_imethd2
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


  !==================================================================================================================================
  subroutine input_read
    implicit none

    integer         :: istat,ierr,iloc,i
    logical         :: linputfile
    character       :: message*124,dchar_grid*200

    !-------------------------------------------------------------------------------
    ! user input file
    !-------------------------------------------------------------------------------
    inquire(file=trim(file_input),exist=linputfile);
    if (.not.linputfile.and.proc_id==0) then
      write(*,*)
      write(*,*) 'linputfile:',linputfile
      write(*,*) 'cannot find "'//trim(file_input)//'" file!'
      write(*,*) 'error in --> mod:input, sub:input_read'
      message=trim('program stopped at "input_read"')
      !call stop(message)
    endif

    open(iunit_input,file=trim(file_input),status='unknown',form='formatted')
    read(iunit_input,*)
    read(iunit_input,*)
    read(iunit_input,*)
    read(iunit_input,'(a)') dchar_grid

    read(iunit_input,*)
    read(iunit_input,*)
    read(iunit_input,*)
    read(iunit_input,*)
    read(iunit_input,*) grad_cellcntr_imethd1, grad_cellcntr_imethd2, grad_cellcntr_lsq_pow
    read(iunit_input,*) grad_limiter_imethd
    read(iunit_input,*) face_reconst_imethd
    read(iunit_input,*) flux_inviscd_imethd
    close(iunit_input)


    !--------------------------------------------------------------------------!
    ! Grid file name
    !--------------------------------------------------------------------------!
    do i=1,len(dchar_grid)
      if (dchar_grid(i:i)=='-') then
        iloc=i-1
        exit
      endif
    enddo
    file_grid(1:iloc)=dchar_grid(1:iloc)


    !--------------------------------------------------------------------------!
    ! Gradient scheme to compute gradient at cell center
    !--------------------------------------------------------------------------!
    lgrad_ggcb=.false.

    lgrad_ggnb=.false.
    lgrad_ggnb_exp=.false.

    lgrad_lsq=.false.
    lgrad_lsq_fn =.false.
    lgrad_lsq_nn =.false.

    if (grad_cellcntr_imethd1==1) then
      lgrad_ggcb=.true.

    elseif (grad_cellcntr_imethd1==2) then
      lgrad_ggnb=.true.
      lgrad_ggnb_exp=.true.

    elseif (grad_cellcntr_imethd1==3) then
      lgrad_lsq=.true.
      if (grad_cellcntr_imethd2==0) then
        lgrad_lsq_fn=.true.

      elseif (grad_cellcntr_imethd2==1) then
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


    !-------------------------------------------------------------------------------
    ! Output file
    !-------------------------------------------------------------------------------
    if (proc_id==0) then
      open(iunit_output,file=trim(file_output),status='unknown',IOSTAT=istat)
      write(iunit_output,'(a)') '==========================================================================================================================================='
      write(iunit_output,'(a)') '     FVM2D CODE                       '
      write(iunit_output,'(a)') '==========================================================================================================================================='
      write(iunit_output,'(a,a)') 'grid file name: ',trim(file_grid)

      write(iunit_output,'(a)') ' '
      write(iunit_output,'(a)') '==========================================================================================================================================='
      write(iunit_output,'(a)') '     Numerics                       '
      write(iunit_output,'(a)') '==========================================================================================================================================='

      !-- Cell-center gradient scheme
      write(iunit_output,'(a,a)') 'Cell-center gradient method: '
      if (lgrad_ggcb) then
        write(iunit_output,'(a)') '   Green-Gauss Cell-Base'
      elseif (lgrad_ggnb) then
        write(iunit_output,'(a)') '   Green-Gauss Node-Base'
      elseif (lgrad_lsq) then
        if (grad_cellcntr_lsq_pow==0.d0) then
          write(iunit_output,'(a)') '   Unweigghted Least-Squeres'
        else
          write(iunit_output,'(a,f3.1)') '   Weigghted Least-Squeres with 1/d^',grad_cellcntr_lsq_pow
        endif
        if (lgrad_lsq_fn) then
          write(iunit_output,'(a)') '   LSQ employs face neighbor stencil'
        elseif (lgrad_lsq_nn) then
          write(iunit_output,'(a)') '   LSQ employs node neighbors stencil'
        endif
      endif

      !-- Gradient limiter
      write(iunit_output,'(a)') ' '
      if (.not.lgrad_limiter) then
        write(iunit_output,'(a)') 'Gradient limiter is not applied'
      elseif (lgrad_limiter) then

      endif

      !-- Face reconstruction
      write(iunit_output,'(a)') ' '
      write(iunit_output,'(a)') 'Face reconstruction method: '
      if (lface_reconst_linear) then
        write(iunit_output,'(a)') '   Linear extrapolation'
      elseif (lface_reconst_umuscl) then
        write(iunit_output,'(a)') '   UMUSCL'
      endif

      !-- Inviscid flux discretization
      write(iunit_output,'(a)') ' '
      write(iunit_output,'(a)') 'Inviscid flux discretization scheme: '
      if (lflux_inviscd_roe) then
        write(iunit_output,'(a)') '   Roe'
      endif

      close(iunit_output)
    endif

    return
  end subroutine input_read

end module input
