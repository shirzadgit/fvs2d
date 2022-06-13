program fvs2d

  use mainparam
  use input
  use grid_procs
  use mpi
  use omp_lib


  implicit none

  integer           :: ierr,errcode
  integer           :: stat(MPI_STATUS_SIZE)
  
  
  !-------------------------------------------------------------------------------
  ! find process ID and how many processes were started
  !-------------------------------------------------------------------------------
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,proc_id,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  if (nproc/=1) then
    write(*,*)'#s of MPI task:',nproc
    write(*,*)'#s of MPI tasks must be equal to 1!'
    call MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
  endif
  
  
  !--------------------------------------------------------------------------------------------------
  ! number of OMP threads
  !--------------------------------------------------------------------------------------------------
  omp_nthreads=1
  lOMP=.false.
!$omp parallel
  omp_nthreads=omp_get_num_threads()
!$omp end parallel
  if (omp_nthreads>1) lOMP=.true.
  
  
  !-------------------------------------------------------------------------------
  ! read input file
  !-------------------------------------------------------------------------------
  call input_read
  
  
  !-------------------------------------------------------------------------------
  ! static grid pre-processing
  !-------------------------------------------------------------------------------
  call grid_procs_init
  
  
  call MPI_FINALIZE(ierr)
  stop 'o.k.'
  
  
end program fvs2d



