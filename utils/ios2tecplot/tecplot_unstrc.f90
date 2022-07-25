module tecplot_unstrc

  implicit none

  Include "tecio.f90"

  private

  !-- tecplot
  integer,save                      :: debug, IsDouble, IsBlock, ZoneType, StrandId,  unused
  character(len=128),save           :: title, variables, FileName, ZoneTitle
  character,save                    :: ScratchDir*2
  integer,pointer                   :: NullPtr(:)
  real(8),save                      :: SolTime
  integer                           :: ierr
  logical,save                      :: lplt,lszplt, lSZL
  integer,save                      :: nvrt_max


  public  :: tecplot_init
  public  :: tecplot_write_grid_solution
  public  :: tecplot_write_grid
  public  :: tecplot_write_solution


contains


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine tecplot_init (ncells_quad)

    implicit none
    integer,intent(in)  :: ncells_quad
    integer             :: i,j,ic,jc

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
    !IsDouble    = 0 !-- 0= single precision, 1=double precision

    ScratchDir  = '.'//char(0)
    ZoneTitle   = 'zone 1'//char(0)

    !-- FETRIANGLE or FEQUADRILATERAL
    if ( ncells_quad==0 ) then
      ZoneType = 2
      nvrt_max = 3
    else
      ZoneType = 3
      nvrt_max = 4
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


    return
  end subroutine tecplot_init


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine tecplot_write_grid_solution (nnodes, ncells, np, grid_dim, imach, file_out, varinfo, cell2node, xy, sol, real_time)
    implicit none

    integer,intent(in)          :: nnodes, ncells, np, grid_dim, imach
    integer(4),intent(in)       :: cell2node(nvrt_max,ncells)
    real(8),intent(in)			    :: real_time
    real(8),intent(in)          :: sol(nnodes,np), xy(nnodes,grid_dim)
    character(len=*),intent(in) :: file_out, varinfo

    real(4),allocatable         :: sol4(:,:), xy4(:,:)
    integer                     :: k,ip,FileFormat, FileType
    character                   :: title1*14

    IsDouble = imach - 1; !1 !-- 0= single precision, 1=double precision

    !if (imach==1) then
    !  IsDouble = 0 !-- 0= single precision, 1=double precision
    !endif

    FileFormat = 0  !-- 0=plt
    FileType = 0    !-- 0=grid & solution, 1=grid, 2=solution

    title1=trim('grid-solution')//char(0)

    variables=trim(varinfo)//char(0)
    !if (lszplt) FileName=trim(trim(file_out))
    !if (lplt)   FileName=trim(trim(file_out)//'.plt')

    !SolTime = SolTime + 1.d0
    SolTime = real_time

    !-- Open the file and write the tecplot datafile header information
    ierr = TECINI142(trim(title1), trim(variables), trim(file_out)//char(0), trim(ScratchDir), FileFormat, FileType, debug, isDouble)

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
    if (IsDouble==1) then
      !-- Write Grid Data
      do ip=1,grid_dim
        ierr = TECDATD142(nnodes,xy(1:nnodes,ip))
      enddo

      !-- Write Solution Data
      do ip=1,np
        ierr = TECDATD142(nnodes,sol(1:nnodes,ip))
      enddo
    else

      allocate( sol4(nnodes,np), xy4(nnodes,grid_dim))
      sol4 = sol
      xy4  = xy

      !-- Write Grid Data
      do ip=1,grid_dim
        ierr = TECDAT142(nnodes,xy4(1:nnodes,ip), IsDouble)
      enddo

      !-- Write Solution Data
      do ip=1,np
        ierr = TECDAT142(nnodes,sol4(1:nnodes,ip), IsDouble)
      enddo

      deallocate(sol4, xy4)

    endif

    !-- Write connectivity
    ierr = TECNOD142(cell2node)

    !-- Close tecplot file
    ierr = TECEND142()

    return
  end subroutine tecplot_write_grid_solution


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine tecplot_write_solution (nnodes, ncells, np, imach, file_out, varinfo, sol, real_time)
    implicit none

    integer,intent(in)          :: nnodes, ncells, np, imach
    real(8),intent(in)			    :: real_time
    real(8),intent(in)          :: sol(nnodes,np)
    character(len=*),intent(in) :: file_out, varinfo

    real(4),allocatable         :: sol4(:,:)
    integer                     :: k,ip,FileFormat, FileType
    character                   :: title1*14

    IsDouble = imach - 1; !1 !-- 0= single precision, 1=double precision

    !if (imach==1) then
    !  IsDouble = 0 !-- 0= single precision, 1=double precision
    !endif

    FileFormat = 0  !-- 0=plt
    FileType = 2    !-- 0=grid & solution, 1=grid, 2=solution

    title1=trim('solution')//char(0)

    variables=trim(varinfo)//char(0)
    !if (lszplt) FileName=trim(trim(file_out))
    !if (lplt)   FileName=trim(trim(file_out)//'.plt')

    !SolTime = SolTime + 1.d0
    SolTime = real_time

    !-- Open the file and write the tecplot datafile header information
    ierr = TECINI142(trim(title1), trim(variables), trim(file_out)//char(0), trim(ScratchDir), FileFormat, FileType, debug, isDouble)

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
    if (IsDouble==1) then
      !-- Write Solution Data
      do ip=1,np
        ierr = TECDATD142(nnodes,sol(1:nnodes,ip))
      enddo
    else

      allocate( sol4(nnodes,np))
      sol4 = sol

      !-- Write Solution Data
      do ip=1,np
        ierr = TECDAT142(nnodes,sol4(1:nnodes,ip), IsDouble)
      enddo

      deallocate(sol4)

    endif

    !-- Close tecplot file
    ierr = TECEND142()

    return
  end subroutine tecplot_write_solution



  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine tecplot_write_grid (nnodes, ncells, grid_dim, imach, file_out, cell2node, xy)
    implicit none

    integer,intent(in)          :: nnodes, ncells, grid_dim, imach
    integer(4),intent(in)       :: cell2node(nvrt_max,ncells)
    real(8),intent(in)          :: xy(nnodes,grid_dim)
    character(len=*),intent(in) :: file_out

    real(4),allocatable         :: xy4(:,:)
    integer                     :: k,ip,FileFormat, FileType
    character                   :: title1*14, variables*10

    IsDouble = imach - 1; !1 !-- 0= single precision, 1=double precision

    !if (imach==1) then
    !  IsDouble = 0 !-- 0= single precision, 1=double precision
    !endif

    FileFormat = 0  !-- 0=plt
    FileType = 1    !-- 0=grid & solution, 1=grid, 2=solution

    title1=trim('grid')//char(0)

    variables='x, y'//char(0)
    !if (lszplt) FileName=trim(trim(file_out))
    !if (lplt)   FileName=trim(trim(file_out)//'.plt')

    !SolTime = SolTime + 1.d0
    SolTime = 0.d0

    !-- Open the file and write the tecplot datafile header information
    ierr = TECINI142(trim(title1), trim(variables), trim(file_out)//char(0), trim(ScratchDir), FileFormat, FileType, debug, isDouble)

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
    if (IsDouble==1) then
      !-- Write Grid Data
      do ip=1,grid_dim
        ierr = TECDATD142(nnodes,xy(1:nnodes,ip))
      enddo
    else

      allocate(  xy4(nnodes,grid_dim))
      xy4  = xy

      !-- Write Grid Data
      do ip=1,grid_dim
        ierr = TECDAT142(nnodes,xy4(1:nnodes,ip), IsDouble)
      enddo

      deallocate(xy4)

    endif

    !-- Write connectivity
    ierr = TECNOD142(cell2node)

    !-- Close tecplot file
    ierr = TECEND142()

    return
  end subroutine tecplot_write_grid

  ! !============================================================================!
  ! !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  ! !============================================================================!
  ! subroutine tecplot_write_grid (file_out,xyz)
  !
  !   implicit none
  !   real(4),intent(in)          :: xyz(nnodes,grid_dim)
  !   character(len=*),intent(in) :: file_out
  !   integer                     :: k,ip,FileFormat, FileType, clen
  !   real(8)                     :: SolTimeG
  !   character(len=6)            :: varG
  !
  !   FileFormat = 0  !-- 0=plt
  !   FileType = 1    !-- 0=grid & solution, 1=grid, 2=solution
  !
  !   title=trim('grid')//char(0)
  !   !if (lszplt) FileName=trim(trim(file_out))
  !   !if (lplt)   FileName=trim(trim(file_out)//'.plt')
  !
  !   SolTimeG = 0.d0
  !
  !   if (grid_dim==1) then
  !     clen=2
  !     varG(1:clen)='x'//char(0)
  !   elseif (grid_dim==2) then
  !     clen=4
  !     varG(1:clen)='x y'//char(0)
  !   elseif (grid_dim==3) then
  !     clen=6
  !     varG(1:clen)='x y z'//char(0)
  !   endif
  !
  !   !-- Open the file and write the tecplot datafile header information
  !   ierr = TECINI142('grid'//char(0), trim(varG(1:clen)), trim(file_out)//char(0), '.'//char(0), FileFormat, FileType, debug, isDouble)
  !
  !   !-- Create ordered zone
  !   ierr = TECZNE142(trim(ZoneTitle),  ZoneType,  nnodes, ncells, 0, 0,0,0,   SolTimeG, StrandId, unused, IsBlock,    &
  !                    0,         &   ! NumFaceConnections
  !                    0,         &   ! FaceNeighborMode
  !                    0,         &   ! TotalNumFaceNodes
  !                    0,         &   ! NumConnectedBoundaryFaces
  !                    0,         &   ! TotalNumBoundaryConnections
  !                    NullPtr,   &   ! PassiveVarList
  !                    NullPtr,   &   ! ValueLocation
  !                    NullPtr,   &   ! ShareVarFromZone
  !                    0)             ! ShareConnectivityFromZone
  !
  !   !-- Write Grid Data
  !   do ip=1,grid_dim
  !     ierr = TECDAT142(nnodes,xyz(1:nnodes,ip),IsDouble)
  !   enddo
  !
  !   !-- Write connectivity
  !   ierr = TECNOD142(cell2node)
  !
  !   !-- Close tecplot file
  !   ierr = TECEND142()
  !
  !   return
  ! end subroutine tecplot_write_grid
  !
  !
  ! !============================================================================!
  ! !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  ! !============================================================================!
  ! subroutine tecplot_write_solution (file_out, varinfo, np, sol, real_time)
  !
  !   implicit none
  !   integer,intent(in)          :: np
  !   real(8),intent(in)			    :: real_time
  !   real(4),intent(in)          :: sol(nnodes,np)
  !   character(len=*),intent(in) :: file_out, varinfo
  !   integer                     :: k,ip,FileFormat, FileType
  !
  !
  !   FileFormat = 0  !-- 0=plt
  !   FileType = 2    !-- 0=grid & solution, 1=grid, 2=solution
  !
  !   title=trim('grid-solution')//char(0)
  !
  !   variables=trim(varinfo)//char(0)
  !
  !   SolTime = real_time
  !
  !   !-- Open the file and write the tecplot datafile header information
  !   ierr = TECINI142('solution'//char(0), trim(variables), trim(file_out)//char(0), trim(ScratchDir), FileFormat, FileType, debug, IsDouble)
  !
  !   !-- Create ordered zone
  !   ierr = TECZNE142(trim(ZoneTitle),  ZoneType,  nnodes, ncells, 0, 0,0,0,   SolTime, StrandId, unused, IsBlock,    &
  !                    0,         &   ! NumFaceConnections
  !                    0,         &   ! FaceNeighborMode
  !                    0,         &   ! TotalNumFaceNodes
  !                    0,         &   ! NumConnectedBoundaryFaces
  !                    0,         &   ! TotalNumBoundaryConnections
  !                    NullPtr,   &   ! PassiveVarList
  !                    NullPtr,   &   ! ValueLocation
  !                    NullPtr,   &   ! ShareVarFromZone
  !                    0)             ! ShareConnectivityFromZone
  !
  !   !-- Write Solution Data
  !   do ip=1,np
  !     ierr = TECDAT142(nnodes,sol(1:nnodes,ip),IsDouble)
  !   enddo
  !
  !   !-- Close tecplot file
  !   ierr = TECEND142()
  !
  !   return
  ! end subroutine tecplot_write_solution
  !
  !
  ! !============================================================================!
  ! !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  ! !============================================================================!
  ! subroutine tecplot_write_solution8 (file_out, varinfo, np, sol, real_time)
  !
  !   implicit none
  !   integer,intent(in)          :: np
  !   real(8),intent(in)			    :: real_time
  !   real(8),intent(in)          :: sol(nnodes,np)
  !   character(len=*),intent(in) :: file_out, varinfo
  !   integer                     :: k,ip,FileFormat, FileType
  !
  !
  !   FileFormat = 0  !-- 0=plt
  !   FileType = 2    !-- 0=grid & solution, 1=grid, 2=solution
  !
  !   title=trim('grid-solution')//char(0)
  !
  !   variables=trim(varinfo)//char(0)
  !
  !   SolTime = real_time
  !
  !   !-- Open the file and write the tecplot datafile header information
  !   ierr = TECINI142('solution'//char(0), trim(variables), trim(file_out)//char(0), trim(ScratchDir), FileFormat, FileType, debug, 1)
  !
  !   !-- Create ordered zone
  !   ierr = TECZNE142(trim(ZoneTitle),  ZoneType,  nnodes, ncells, 0, 0,0,0,   SolTime, StrandId, unused, IsBlock,    &
  !                    0,         &   ! NumFaceConnections
  !                    0,         &   ! FaceNeighborMode
  !                    0,         &   ! TotalNumFaceNodes
  !                    0,         &   ! NumConnectedBoundaryFaces
  !                    0,         &   ! TotalNumBoundaryConnections
  !                    NullPtr,   &   ! PassiveVarList
  !                    NullPtr,   &   ! ValueLocation
  !                    NullPtr,   &   ! ShareVarFromZone
  !                    0)             ! ShareConnectivityFromZone
  !
  !   !-- Write Solution Data
  !   do ip=1,np
  !     ierr = TECDATD142(nnodes,sol(1:nnodes,ip))
  !   enddo
  !
  !   !-- Close tecplot file
  !   ierr = TECEND142()
  !
  !   return
  ! end subroutine tecplot_write_solution8

end module tecplot_unstrc
