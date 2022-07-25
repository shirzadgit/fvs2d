module grid_procs

  integer,parameter       :: grid_dim=2
  integer,save            :: nnodes, ncells_tri, ncells_quad, ncells
  real,allocatable        :: node_xy(:,:)
  integer(4),allocatable  :: cell2node(:,:)

  integer,parameter       :: iunit_grid=901

contains


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine read_grid2d (file_grid)

    implicit none

    character,intent(in)    :: file_grid*127
    integer                 :: i,j,istat
    logical                 :: linputfile

    !--------------------------------------------------------------------------!
    ! check grid file
    !--------------------------------------------------------------------------!
    inquire(file=trim(file_grid),exist=linputfile);
    if (.not.linputfile) then
      write(*,*)
      write(*,*) 'linputfile:',linputfile
      write(*,'(a,a,a)') 'cannot find ',trim(file_grid),' file!'
      write(*,*) 'error in --> mod:grid_procs, sub:grid_read'
      stop 'program stopped at "grid_read"'
    endif

    !-- open grid file and skip the first line
    open (iunit_grid,file=trim(file_grid),status='unknown',IOSTAT=istat)
    read (iunit_grid,*)
    read (iunit_grid,*) nnodes, ncells_tri, ncells_quad

    !--  allocate node and cell arrays.
    ncells = ncells_tri + ncells_quad
    allocate(node_xy(nnodes,grid_dim))

    if (ncells_quad>0)  allocate(cell2node(4,ncells))
    if (ncells_quad==0) allocate(cell2node(3,ncells))

    !-- read node/vertices coordinates
    do i=1,nnodes
      read (iunit_grid,*) node_xy(i,1), node_xy(i,2)
    enddo

    !-- read triangle cells
    do i=1,ncells_tri
      read (iunit_grid,*) cell2node(1,i), cell2node(2,i), cell2node(3,i)
    enddo
    if (ncells_quad>0) cell2node(4, 1:ncells_tri) = cell2node(3, 1:ncells_tri)

    !-- read quad cells
    do j=1,ncells_quad
      i=j+ncells_tri
      read (iunit_grid,*) cell2node(1,i), cell2node(2,i), cell2node(3,i), cell2node(4,i)
    enddo

    !-- close grid file
    close(iunit_grid)

    return
  end subroutine read_grid2d

end module grid_procs
