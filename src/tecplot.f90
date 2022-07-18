module tecplot

  use mainparam,  only    : grid_dim
  use data_grid,  only    : ncells, nnodes, ncells_quad, ncells_tri, nedges, cell

  implicit none

  Include "tecio.f90"

  private

  !-- tecplot
  integer,save                      :: debug, IsDouble, IsBlock, ZoneType, StrandId,  unused
  character(len=128),save           :: title, variables, FileName,ZoneTitle
  character,save                    :: ScratchDir*2
  integer,pointer                   :: NullPtr(:)
  real(8),save                      :: SolTime
  !character(len=1),parameter        :: nullchar=char(0)
  integer                           :: ierr
  logical,save                      :: lplt,lszplt, lSZL

  !-- file
  integer,save                      :: nedge_max
  integer(kind=4),allocatable,save  :: cell2node(:,:)

  public  :: tecplot_init
  public  :: tecplot_write_grid_solution
  public  :: tecplot_write_grid
  public  :: tecplot_write_solution, tecplot_write_solution8


contains


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine tecplot_init
    implicit none
    integer     :: i,j,ic,jc

    ! title     - Title of the data set, must be null terminated.
    ! variables - List of variable names (seprated with space or comma), must be null terminated.
    ! FileName  - Name of the file to create, must be null terminated.
    ! ScratchDir- Name of the directory to put the scratch file, must be null terminated.
    ! ZoneTitle - The title of the zone, must be null-terminated.

    debug       = 0 !-- 0=no debugging, 1=debug
    ZoneType    = 0 !-- 0=ORDERED, 1=FELINESEG, 2=FETRIANGLE, 3=FEQUADRILATERAL, 4=FETETRAHEDRON, 5=FEBRICK
    StrandId    = 1 !-- 0=zone is static and not associated with a strand, >0=zone is assigned to a given strand, -1=strand id for this zone is pending
    unused      = 0 !-- ParentZone is no longer used
    IsBlock     = 1 !-- 1= data will be passed in BLOCK (i,j,k,v), 0=POINT format (v,i,j,k)
    IsDouble    = 0 !-- 0= single precision, 1=double precision

    ScratchDir  = '.'//char(0)
    ZoneTitle   = 'zone 1'//char(0)

    !-- FETRIANGLE or FEQUADRILATERAL
    if ( ncells_quad==0 ) then
      nedge_max= 3
      ZoneType = 2
    else
      nedge_max= 4
      ZoneType = 3
    endif


    !-- SZL Technolog
    !-- FileFormat --> 0=plt, 1=szplt (szplt is not complabit with separate grid and solution files)
    lSZL = .false.
    if (lSZL) then
      lszplt=.true.
      !FileFormat = 1
      lplt=.false.
    else
      lplt=.true.
      !FileFormat = 0
      lszplt=.false.
    endif


    !-- cell 2 node link
    allocate(cell2node(nedge_max,ncells))
    do ic=1,ncells_tri
      do i=1,cell(ic)%nvrt
        cell2node(i,ic) = cell(ic)%node(i)
      enddo
      cell2node(nedge_max,ic) = cell(ic)%node(3)
    enddo
    do jc=1,ncells_quad
      ic=jc + ncells_tri
      do i=1,cell(ic)%nvrt
        cell2node(i,ic) = cell(ic)%node(i)
      enddo
    enddo

    return
  end subroutine tecplot_init


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine tecplot_write_grid_solution (file_out, varinfo, np,  xyz, sol, real_time)
    implicit none

    integer,intent(in)          :: np
    real(8),intent(in)			    :: real_time
    real(4),intent(in)          :: sol(nnodes,np), xyz(nnodes,grid_dim)
    character(len=*),intent(in) :: file_out, varinfo
    integer                     :: k,ip,FileFormat, FileType
    character                   :: title1*14

    FileFormat = 0  !-- 0=plt
    FileType = 0    !-- 0=grid & solution, 1=grid, 2=solution

    title1=trim('grid-solution')//char(0)

    variables=trim(varinfo)//char(0)
    !if (lszplt) FileName=trim(trim(file_out))
    !if (lplt)   FileName=trim(trim(file_out)//'.plt')

    !SolTime = SolTime + 1.d0
    SolTime = real_time

    !-- Open the file and write the tecplot datafile header information
    ierr = TECINI142(trim(title1), trim(variables), trim(file_out), trim(ScratchDir), FileFormat, FileType, debug, isDouble)

    !-- Create ordered zone
    ierr = TECZNE142(trim(ZoneTitle),  ZoneType,  nnodes, ncells, 0, 0,0,0,   SolTime, StrandId, unused, IsBlock,    &
                     0,         &   ! NumFaceConnections
                     0,         &   ! FaceNeighborMode
                     0,         &   ! TotalNumFaceNodes
                     0,         &   ! NumConnectedBoundaryFaces
                     0,         &   ! TotalNumBoundaryConnections
                     NullPtr,   &   ! PassiveVarList
                     NullPtr,   &   ! ValueLocation
                     NullPtr,   &   ! ShareVarFromZone
                     0)             ! ShareConnectivityFromZone

    !-- Write Grid Data
    do ip=1,grid_dim
      ierr = TECDAT142(nnodes,xyz(1:nnodes,ip),IsDouble)
    enddo

    !-- Write Solution Data
    do ip=1,np
      ierr = TECDAT142(nnodes,sol(1:nnodes,ip),IsDouble)
    enddo

    !-- Write connectivity
    ierr = TECNOD142(cell2node)

    !-- Close tecplot file
    ierr = TECEND142()

    return
  end subroutine tecplot_write_grid_solution


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine tecplot_write_grid (file_out,xyz)

    implicit none
    real(4),intent(in)          :: xyz(nnodes,grid_dim)
    character(len=*),intent(in) :: file_out
    integer                     :: k,ip,FileFormat, FileType, clen
    real(8)                     :: SolTimeG
    character(len=6)            :: varG

    FileFormat = 0  !-- 0=plt
    FileType = 1    !-- 0=grid & solution, 1=grid, 2=solution

    title=trim('grid')//char(0)
    !if (lszplt) FileName=trim(trim(file_out))
    !if (lplt)   FileName=trim(trim(file_out)//'.plt')

    SolTimeG = 0.d0

    if (grid_dim==1) then
      clen=2
      varG(1:clen)='x'//char(0)
    elseif (grid_dim==2) then
      clen=4
      varG(1:clen)='x y'//char(0)
    elseif (grid_dim==3) then
      clen=6
      varG(1:clen)='x y z'//char(0)
    endif

    !-- Open the file and write the tecplot datafile header information
    ierr = TECINI142('grid'//char(0), trim(varG(1:clen)), trim(file_out)//char(0), '.'//char(0), FileFormat, FileType, debug, isDouble)

    !-- Create ordered zone
    ierr = TECZNE142(trim(ZoneTitle),  ZoneType,  nnodes, ncells, 0, 0,0,0,   SolTimeG, StrandId, unused, IsBlock,    &
                     0,         &   ! NumFaceConnections
                     0,         &   ! FaceNeighborMode
                     0,         &   ! TotalNumFaceNodes
                     0,         &   ! NumConnectedBoundaryFaces
                     0,         &   ! TotalNumBoundaryConnections
                     NullPtr,   &   ! PassiveVarList
                     NullPtr,   &   ! ValueLocation
                     NullPtr,   &   ! ShareVarFromZone
                     0)             ! ShareConnectivityFromZone

    !-- Write Grid Data
    do ip=1,grid_dim
      ierr = TECDAT142(nnodes,xyz(1:nnodes,ip),IsDouble)
    enddo

    !-- Write connectivity
    ierr = TECNOD142(cell2node)

    !-- Close tecplot file
    ierr = TECEND142()

    return
  end subroutine tecplot_write_grid


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine tecplot_write_solution (file_out, varinfo, np, sol, real_time)

    implicit none
    integer,intent(in)          :: np
    real(8),intent(in)			    :: real_time
    real(4),intent(in)          :: sol(nnodes,np)
    character(len=*),intent(in) :: file_out, varinfo
    integer                     :: k,ip,FileFormat, FileType


    FileFormat = 0  !-- 0=plt
    FileType = 2    !-- 0=grid & solution, 1=grid, 2=solution

    title=trim('grid-solution')//char(0)

    variables=trim(varinfo)//char(0)

    SolTime = real_time

    !-- Open the file and write the tecplot datafile header information
    ierr = TECINI142('solution'//char(0), trim(variables), trim(file_out)//char(0), trim(ScratchDir), FileFormat, FileType, debug, IsDouble)

    !-- Create ordered zone
    ierr = TECZNE142(trim(ZoneTitle),  ZoneType,  nnodes, ncells, 0, 0,0,0,   SolTime, StrandId, unused, IsBlock,    &
                     0,         &   ! NumFaceConnections
                     0,         &   ! FaceNeighborMode
                     0,         &   ! TotalNumFaceNodes
                     0,         &   ! NumConnectedBoundaryFaces
                     0,         &   ! TotalNumBoundaryConnections
                     NullPtr,   &   ! PassiveVarList
                     NullPtr,   &   ! ValueLocation
                     NullPtr,   &   ! ShareVarFromZone
                     0)             ! ShareConnectivityFromZone

    !-- Write Solution Data
    do ip=1,np
      ierr = TECDAT142(nnodes,sol(1:nnodes,ip),IsDouble)
    enddo

    !-- Close tecplot file
    ierr = TECEND142()

    return
  end subroutine tecplot_write_solution


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine tecplot_write_solution8 (file_out, varinfo, np, sol, real_time)

    implicit none
    integer,intent(in)          :: np
    real(8),intent(in)			    :: real_time
    real(8),intent(in)          :: sol(nnodes,np)
    character(len=*),intent(in) :: file_out, varinfo
    integer                     :: k,ip,FileFormat, FileType


    FileFormat = 0  !-- 0=plt
    FileType = 2    !-- 0=grid & solution, 1=grid, 2=solution

    title=trim('grid-solution')//char(0)

    variables=trim(varinfo)//char(0)

    SolTime = real_time

    !-- Open the file and write the tecplot datafile header information
    ierr = TECINI142('solution'//char(0), trim(variables), trim(file_out)//char(0), trim(ScratchDir), FileFormat, FileType, debug, 1)

    !-- Create ordered zone
    ierr = TECZNE142(trim(ZoneTitle),  ZoneType,  nnodes, ncells, 0, 0,0,0,   SolTime, StrandId, unused, IsBlock,    &
                     0,         &   ! NumFaceConnections
                     0,         &   ! FaceNeighborMode
                     0,         &   ! TotalNumFaceNodes
                     0,         &   ! NumConnectedBoundaryFaces
                     0,         &   ! TotalNumBoundaryConnections
                     NullPtr,   &   ! PassiveVarList
                     NullPtr,   &   ! ValueLocation
                     NullPtr,   &   ! ShareVarFromZone
                     0)             ! ShareConnectivityFromZone

    !-- Write Solution Data
    do ip=1,np
      ierr = TECDATD142(nnodes,sol(1:nnodes,ip))
    enddo

    !-- Close tecplot file
    ierr = TECEND142()

    return
  end subroutine tecplot_write_solution8

end module tecplot
