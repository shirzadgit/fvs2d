module io

  use mainparam
  use data_grid
  use data_solution
  use tecplot
  use interpolation


  integer,save,private  :: counter=0


contains

  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine io_write_inst
    implicit none

    integer                     :: ivar, np, i
    real(kind=4),allocatable    :: f4(:,:),xy4(:,:)
    real,allocatable            :: f8(:)
    real                        :: sol_time
    character(len=*)            :: file_out*100, varinfo*100, dchar*5


    counter = counter + 1
    write(dchar,'(i5.5)') counter

    !--
    if (counter == 1) then
      allocate(xy4(nnodes,2))
      do i=1,nnodes
        xy4(i,1) = node(i)%x
        xy4(i,2) = node(i)%y
      enddo
      file_out = trim('solution_grid.plt')
      call tecplot_write_grid (trim(file_out), xy4)

      deallocate(xy4)
    endif

    allocate(f8(nnodes), f4(nnodes,nvar))
    do ivar=1,nvar
      call interpolate_cell2node(pvar(ivar,1:ncells), f8)
      f4(1:nnodes,ivar) = f8(1:nnodes)
    enddo


    !--
    file_out = trim('solution_it'//dchar(1:5)//'.plt')
    varinfo  = 'rho u v p'
    np = 4
    sol_time = dble(counter)
    call tecplot_write_solution (trim(file_out), trim(varinfo), np, f4, sol_time)

    deallocate(f4,f8)

    return
  end subroutine io_write_inst

end module io
