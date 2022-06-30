program fvs2d

  use mainparam
  use input
  use data_grid
  use data_sol
  use grid_procs
  use interpolation
  use gradient
  use mms_euler2d
  use test

  use mpi
  use omp_lib


  implicit none

  integer           :: ierr,errcode
  integer           :: stat(MPI_STATUS_SIZE)
  real              :: t1


  !----------------------------------------------------------------------------!
  ! find process ID and how many processes were started
  !----------------------------------------------------------------------------!
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,proc_id,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  if (nproc/=1) then
    write(*,*)'#s of MPI task:',nproc
    write(*,*)'#s of MPI tasks must be equal to 1!'
    call MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
  endif


  !----------------------------------------------------------------------------!
  ! number of OMP threads
  !----------------------------------------------------------------------------!
  omp_nthreads=1
  lOMP=.false.
!$omp parallel
  omp_nthreads=omp_get_num_threads()
!$omp end parallel
  if (omp_nthreads>1) lOMP=.true.


  !----------------------------------------------------------------------------!
  ! read input file
  !----------------------------------------------------------------------------!
  call input_read


  !----------------------------------------------------------------------------!
  ! static grid pre-processing
  !----------------------------------------------------------------------------!
  call grid_procs_init


  !----------------------------------------------------------------------------!
  ! solution variable allocatation
  !----------------------------------------------------------------------------!
  call data_sol_allocate


  !----------------------------------------------------------------------------!
  ! interpolation setup
  !----------------------------------------------------------------------------!
  call interpolate_init


  !----------------------------------------------------------------------------!
  ! gradient setup
  !----------------------------------------------------------------------------!
  call gradient_init


  !----------------------------------------------------------------------------!
  ! initialize runge-kutta: setup coefficents
  !----------------------------------------------------------------------------!
  !call runge_kutta_init


  !----------------------------------------------------------------------------!
  ! compute source term for euler equation
  !----------------------------------------------------------------------------!
  call mms_euler2d_init


  !----------------------------------------------------------------------------!
  ! initialize flow
  !----------------------------------------------------------------------------!
  !t1 = dble(ntstart-1)*dt
  t1=0.d0
  call initialize(t1)


  !----------------------------------------------------------------------------!
  ! open ios files
  !----------------------------------------------------------------------------!
  !call io_init


  !----------------------------------------------------------------------------!
  ! test operators
  !----------------------------------------------------------------------------!
  call test_init


  call MPI_FINALIZE(ierr)
  stop 'o.k.'


end program fvs2d
