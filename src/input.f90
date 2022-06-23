module input

  use mainparam
  use mpi

  implicit none

  !-- Grid paramteres
  character,save          :: file_grid*127

  integer,save            :: grad_cellcntr_imethd1, grad_cellcntr_imethd2
  real,save               :: grad_cellcntr_lsq_pow
  logical,save            :: lgrad_ggcb, lgrad_ggnb, lgrad_ggnb_exp, lgrad_lsq
  logical,save            :: lgrad_lsq_fn1, lgrad_lsq_fn2, lgrad_lsq_nn

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
    lgrad_lsq_fn1=.false.
    lgrad_lsq_fn2=.false.
    lgrad_lsq_nn =.false.

    if (grad_cellcntr_imethd1==1) then
      lgrad_ggcb=.true.

    elseif (grad_cellcntr_imethd1==2) then
      lgrad_ggnb=.true.
      lgrad_ggnb_exp=.true.

    elseif (grad_cellcntr_imethd1==3) then
      lgrad_lsq=.true.
      if (grad_cellcntr_imethd2==0) then
        lgrad_lsq_fn1=.true.

      elseif (grad_cellcntr_imethd2==1) then
        lgrad_lsq_fn2=.true.

      elseif (grad_cellcntr_imethd2==2) then
        lgrad_lsq_nn=.true.
      endif

    else

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
      write(iunit_output,'(a)') 'Gradient method'
      if (lgrad_ggcb) then
        write(iunit_output,'(a,a)') 'Scheme: ','Green-Gauss Cell-Base'
      elseif (lgrad_ggnb) then
        write(iunit_output,'(a,a)') 'gradient scheme: ','Green-Gauss Node-Base'
      elseif (lgrad_lsq) then
        write(iunit_output,'(a,a)') 'gradient scheme: ','Least-Squeres'
        if (grad_cellcntr_lsq_pow==0.d0) then
          write(iunit_output,'(a)') 'Unweigghted LSQ'
        else
          write(iunit_output,'(a,f3.1)') 'Weigghted LSQ with 1/d^',grad_cellcntr_lsq_pow
        endif
        if (lgrad_lsq_fn1) then
          write(iunit_output,'(a)') 'LSQ employs face neighbor stencil'
        elseif (lgrad_lsq_fn2) then
          write(iunit_output,'(a)') 'LSQ employs face neighbors of faceneighbors stencil'
        elseif (lgrad_lsq_nn) then
          write(iunit_output,'(a)') 'LSQ employs node neighbors stencil'
        endif
      endif

      close(iunit_output)
    endif

    return
  end subroutine input_read

end module input
