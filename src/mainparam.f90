module mainparam
  implicit none

  !-- #'s MPI processors  
  integer,save            :: nproc,proc_id  

  !-- OpenMP parameters
  integer,save            :: omp_nthreads
  logical,save            :: lOMP

  !-- input file
  character(len=127)      :: file_input='input.fvm2d'
  integer,parameter       :: iunit_input=900

  !-- grid file
  integer,parameter       :: iunit_grid=901

  !-- log files
  character(len=127)      :: file_output='log.fvm2d'
  integer,parameter       :: iunit_output=902

end module mainparam
