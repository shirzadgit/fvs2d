module mainparam
  implicit none

  !--number of variables (2d --> nvar=4, rho,u,v,p;   3d --> nvar=4, rho,u,v,w,p)
  integer,parameter       :: nvar = 4, grid_dim = 2

  !-- #'s MPI processors
  integer,save            :: nproc,proc_id

  !-- OpenMP parameters
  integer,save            :: omp_nthreads
  logical,save            :: lOMP

  !-- input file
  character(len=127)      :: file_input='fvs2d.input'
  integer,parameter       :: iunit_input=900

  !-- grid file
  character(len=5)        :: grid_ext='.grid'
  character(len=3)        ::   bc_ext='.bc'
  integer,parameter       :: iunit_grid=901
  integer,parameter       :: iunit_bc=902

  !-- isentropic vortex
  character(len=15)       :: file_vortex='fvs2d.vortex'
  integer,parameter       :: iunit_vortex=903

  !-- cp file
  character(len=15)       :: file_cp='log_cp.plt'
  integer,parameter       :: iunit_cp=904

  !-- Vn file
  character(len=15)       :: file_un='log_un.plt'
  integer,parameter       :: iunit_un=905

  !-- log files
  character(len=127)      :: file_log_input='log.fvs2d'
  integer,parameter       :: iunit_log_input=906

  character(len=127)      :: file_log_grid='log.grid'
  integer,parameter       :: iunit_log_grid=907

end module mainparam
