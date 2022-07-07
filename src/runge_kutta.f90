module runge_kutta

  use mainparam,  only  : nvar
  use input,      only  : dtsubstep, nsubsteps, rk_nstages, rk_order, lSSPRK
  use data_grid,  only  : ncells
  use data_solution
  use residual
  use mpi

  implicit none

  private
  real,allocatable,save   :: rk_coef(:),h_rk(:),dts(:),dte(:)
  integer,save            :: istep_rk,icont
  real,save               :: cput_drdt, cput_uvwetc
  real,save               :: h

  public  :: runge_kutta_init,istep_rk,cput_drdt,cput_uvwetc,time_integration

contains

  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine runge_kutta_init

    implicit none

    real              :: ak,bk,ck,dk,ek

    h = dtsubstep
    allocate(h_rk(rk_nstages),rk_coef(rk_nstages),dts(rk_nstages),dte(rk_nstages))

    !--------------------------------------------------------------------------!
    ! coefficent
    !--------------------------------------------------------------------------!
    if (rk_nstages==4) then
      if (rk_order==1) then
        ak=3.60897d0; bk=2.04d0;  ck=0.34206d0;  dk=0.00897d0
      elseif (rk_order==2) then
        ak=0.11d0;    bk=3.92d0;  ck=1.86d0;     dk=0.11d0
      elseif (rk_order==3) then
        ak=0.65;      bk=2.7d0;   ck=2.d0;       dk=0.65d0
      elseif (rk_order==4) then
        ak=1.d0;      bk=2.d0;    ck=2.d0;       dk=1.d0
      endif

      h_rk(1)=h/2.d0;  rk_coef(1)=ak
      h_rk(2)=h/2.d0;  rk_coef(2)=bk
      h_rk(3)=h;       rk_coef(3)=ck
      h_rk(4)=h/6.d0;  rk_coef(4)=dk
      write(*,*) h
    endif

    !--------------------------------------------------------------------------!
    ! modify coefficents according to SSP
    !--------------------------------------------------------------------------!
    if (rk_nstages==4.and.rk_order==2.and.lSSPRK) then
      h_rk(1)=h/3.d0;  rk_coef(1)=1.d0
      h_rk(2)=h/3.d0;  rk_coef(2)=1.d0
      h_rk(3)=h/3.d0;  rk_coef(3)=1.d0
      h_rk(4)=h/4.d0;  rk_coef(4)=1.d0
      dts(1)=0.d0;        dte(1)=h/3.d0
      dts(2)=h/3.d0;      dte(2)=h*2.d0/3.d0
      dts(3)=h*2.d0/3.d0; dte(3)=h
      dts(4)=h;           dte(4)=h
    endif

    icont=0;
    cput_drdt=0.d0
    cput_uvwetc=0.d0

    return
  end subroutine runge_kutta_init


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine time_integration (t1)
    implicit none
    real,intent(in)     :: t1

    if (lSSPRK) then
      call time_integ_SSPRK(t1)
    else
      call time_integ_RK(t1)
    endif

    return
  end subroutine time_integration


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine time_integ_RK (t1)

    implicit none
    real,intent(in)     :: t1
    integer             :: rk,istep
    real                :: t,told,t12,tnew,tstart(rk_nstages),tend(rk_nstages),cput1,cput2
    real,allocatable    :: fcvar(:,:), cvar0(:,:)

    !-- allocation
    allocate(fcvar(nvar,ncells), cvar0(nvar,ncells))


    !-- loop to integrate from t1 to t2
    do istep=1,nsubsteps
      t = t1 + dble(istep-1)*h

      told = t
      t12  = told + 0.5d0*h
      tnew = told + h

      icont=icont+1

      if (rk_nstages==4) then
        tstart(1)=told; tstart(2)=t12;  tstart(3)=t12;  tstart(4)=tnew;
          tend(1)=t12;    tend(2)=t12;    tend(3)=tnew;   tend(4)=tnew;
      endif

      fcvar = 0.d0
      cvar0 = cvar

      !-- loop over R-K stages
      do rk=1,rk_nstages

        cput1=MPI_WTIME()
        call compute_residual;
        cput2=MPI_WTIME(); cput_drdt=cput_drdt + cput2 - cput1; 

         fcvar(1:nvar, 1:ncells) = fcvar(1:nvar, 1:ncells) + rk_coef(rk) * resid(1:nvar, 1:ncells)

        if (rk<rk_nstages) then
          cvar(1:nvar, 1:ncells) = cvar0(1:nvar, 1:ncells) +    h_rk(rk) * resid(1:nvar, 1:ncells)
        elseif (rk==rk_nstages) then
          cvar(1:nvar, 1:ncells) = cvar0(1:nvar, 1:ncells) +    h_rk(rk) * fcvar(1:nvar, 1:ncells)
          !write(*,*) maxval(abs(resid)),minval(abs(resid))
        endif

        !cput1=MPI_WTIME();
        !call uvw_etc (tend(rk))
        !cput2=MPI_WTIME(); cput_uvwetc = cput_uvwetc + cput2 - cput1
      enddo

      !if (lavrun) call avrun
      !call cfl_dfl
    enddo

    !-- deallocate
    deallocate(fcvar, cvar0)

    return
  end subroutine time_integ_RK


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine time_integ_SSPRK(t1)

    implicit none
    real,intent(in)     :: t1
    integer             :: rk,istep
    real                :: t,told,tstart(rk_nstages),tend(rk_nstages),cput1,cput2
    real,allocatable    :: fcvar(:,:), cvar0(:,:)

    !-- allocation
    allocate(fcvar(nvar,ncells), cvar0(nvar,ncells))

    !-- loop to integrate from t1 to t2
    do istep=1,nsubsteps
      t = t1 + dble(istep-1)*h
      told = t
      icont=icont+1

      fcvar = 0.d0
      cvar0 = cvar

      !-- loop over R-K stages
      do rk=1,rk_nstages
        tstart= told + dts(rk)
        tend  = told + dte(rk)

        cput1=MPI_WTIME()
        call compute_residual
        cput2=MPI_WTIME(); cput_drdt=cput_drdt + cput2 - cput1

         cvar(1:nvar, 1:ncells) = cvar0(1:nvar, 1:ncells) + h_rk(rk) * ( rk_coef(rk)*resid(1:nvar, 1:ncells) + fcvar(1:nvar, 1:ncells) )
        fcvar(1:nvar, 1:ncells) = fcvar(1:nvar, 1:ncells) + resid(1:nvar, 1:ncells)

        !cput1=MPI_WTIME();
        !call uvw_etc (tend)
        !cput2=MPI_WTIME(); cput_uvwetc = cput_uvwetc + cput2 - cput1
      enddo

      !if (lavrun) call avrun
      !call cfl_dfl
    enddo

    !-- deallocate
    deallocate(fcvar, cvar0)

    return
  end subroutine time_integ_SSPRK


end module runge_kutta
