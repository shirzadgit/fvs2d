module data_grid

  implicit none

  !----------------------------------------------------------------------------!
  ! grid dimension
  !----------------------------------------------------------------------------!
  integer,save                    :: nnodes
  integer,save                    :: ncells_tri, ncells_quad, ncells
  integer,save                    :: nedges
  integer,save                    :: ncells_intr, ncells_bndr
  integer,allocatable,save        :: cell_intr(:), cell_bndr(:)
  integer,save                    :: nedges_intr, nedges_bndr
  integer,allocatable,save        :: edge_intr(:), edge_bndr(:), bedge_cell(:)

  real,save                       :: heff1, heff2
  integer,save                    :: nbndries


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
  ! data type for boundaries
  !----------------------------------------------------------------------------!
  type bc_type
    integer                       :: ncells, nedges
    integer,dimension(:),pointer  :: cell, edge
    character(len=80)             :: type
  end type bc_type


  !----------------------------------------------------------------------------!
  ! allocatable data
  !----------------------------------------------------------------------------!
  type(node_type),dimension(:),pointer  :: node
  type(edge_type),dimension(:),pointer  :: edge
  type(cell_type),dimension(:),pointer  :: cell
  type(bc_type)  ,dimension(:),pointer  :: bndry


contains


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine data_grid_close
    implicit none

    deallocate(cell)
    deallocate(edge)
    deallocate(node)
    deallocate(bndry)

    return
  end subroutine data_grid_close

end module data_grid
