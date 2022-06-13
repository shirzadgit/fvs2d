module input

  use mainparam
  use mpi

  implicit none

  !-- Grid paramteres
  character,save          :: file_grid*127


  contains


  !==================================================================================================================================
  subroutine input_read
    implicit none

    integer         :: istat,ierr,iloc,i
    logical         :: linputfile
    character       :: message*124,dchar*200

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
    read(iunit_input,'(a)') dchar
    close(iunit_input)

    do i=1,len(dchar)
      if (dchar(i:i)=='-') then
        iloc=i-1
        exit
      endif
    enddo
    file_grid(1:iloc)=dchar(1:iloc)    



    !-------------------------------------------------------------------------------
    ! Output file
    !-------------------------------------------------------------------------------
    if (proc_id==0) then
      open(iunit_output,file=trim(file_output),status='unknown',IOSTAT=istat)
      write(iunit_output,'(a)') '==========================================================================================================================================='
      write(iunit_output,'(a)') '     FVM2D CODE                       '
      write(iunit_output,'(a)') '==========================================================================================================================================='
      write(iunit_output,'(a,a)') 'grid file name: ',trim(file_grid)
      close(iunit_output)
    endif

    return
  end subroutine input_read

end module input
