program fvs2d

  use mainparam
  use input
  use data_grid
  use data_solution
  use grid_procs
  use interpolation
  use gradient
  use runge_kutta
  use residual
  use mms
  use initialize
  use tecplot
  use io
  use test

  use mpi
  use omp_lib


  implicit none

  integer           :: ierr,errcode
  integer           :: stat(MPI_STATUS_SIZE)

  integer           :: istep,istep_save,istep_av,it, it_tot
  real              :: t1,t2, t0
  real              :: cput1,cput2,cputot
  real,parameter    :: min=60.d0, hr=3600.d0


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
  call data_solution_init


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
  call runge_kutta_init


  !----------------------------------------------------------------------------!
  ! initialize residual
  !----------------------------------------------------------------------------!
  !call residual_init


  !----------------------------------------------------------------------------!
  ! compute manufactured solutio and source terms
  !----------------------------------------------------------------------------!
  call mms_init


  !----------------------------------------------------------------------------!
  ! initialize flow
  !----------------------------------------------------------------------------!
  call initialize_solution


  !----------------------------------------------------------------------------!
  ! open ios files
  !----------------------------------------------------------------------------!
  !call io_init
  call tecplot_init


  !----------------------------------------------------------------------------!
  ! test operators
  !----------------------------------------------------------------------------!
  !call test_init
  !stop 'ok'

  !----------------------------------------------------------------------------!
  ! start time integration
  !----------------------------------------------------------------------------!
  cput1 = MPI_WTIME()
  istep = 0
  istep_save=0
  istep_av=0
  it_tot=0
  t0 = dble(ntstart-1)*dt
  !do it=ntstart,ntstart+ntimes-1
  do it=1,nsaves
    istep = istep + 1

    !--begin/end times
    t1 = t0   !dble(it-1)*dtsave(it) + t0
    t2 = t1 + dtsave(it)

    !-- integrate from t1 to t2
    call time_integration (t1, nsubsteps(it))

    !-- output the solution at check point for user
    it_tot = it_tot + nsubsteps(it)
    write(*,'(i5,a)') it_tot,' time-steps done'

    !-- write out output files
    call io_write_inst(t2, nsubsteps(it))

    t0 = t2;
  enddo


  !----------------------------------------------------------------------------!
  ! end
  !----------------------------------------------------------------------------!
  call MPI_FINALIZE(ierr)
  stop 'o.k.'


end program fvs2d
