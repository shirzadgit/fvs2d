module data_grid

  implicit none

  !----------------------------------------------------------------------------!
  ! grid dimension
  !----------------------------------------------------------------------------!
  integer,save                    :: nnodes
  integer,save                    :: ncells_tri, ncells_quad, ncells
  integer,save                    :: nedges


  !----------------------------------------------------------------------------!
  ! data type for cells
  !----------------------------------------------------------------------------!
  type cell_type
    integer                       :: nvrt, nnghbrs
    integer,dimension(:),pointer  :: node, nghbr, nghbre
    integer,dimension(:),pointer  :: edge
    real                          :: x,y,vol
    real,dimension(:),pointer     :: nrmlsign
    real,dimension(:,:),pointer   :: pos2edg
  end type cell_type


  !----------------------------------------------------------------------------!
  ! data type for edges/faces
  !----------------------------------------------------------------------------!
  type edge_type
    integer                       :: n1,n2, c1,c2
    real                          :: x,y,area,nx,ny,tx,ty
  end type edge_type


  !----------------------------------------------------------------------------!
  ! data type for nodes/vertices
  !----------------------------------------------------------------------------!
  type node_type
    integer                       :: ncells
    integer,dimension(:),pointer  :: cell
    real                          :: x,y
  end type node_type


  !----------------------------------------------------------------------------!
  ! allocatable data
  !----------------------------------------------------------------------------!
  type(node_type),dimension(:),pointer  :: node
  type(edge_type),dimension(:),pointer  :: edge
  type(cell_type),dimension(:),pointer  :: cell


contains


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine data_grid_close
    implicit none

    deallocate(cell)
    deallocate(edge)
    deallocate(node)

    return
  end subroutine data_grid_close

end module data_grid
