module runge_kutta

  use mainparam,      only  : nvar, iunit_res, file_res
  use input,          only  : rk_nstages, rk_order, lSSPRK, lvortex, dt, cfl_user, lsteady, lwrite_resid
  use data_grid,      only  : ncells, cell
  use data_solution,  only  : cvar, resid, ws_nrml, ir,iu,iv,ip
  use residual,       only  : compute_residual
  use mms,            only  : error_isentropic_vortex
  use mpi

  implicit none

  private
  real,allocatable,save   :: rk_coef(:),h_rk(:),dts(:),dte(:), dt_local(:)
  integer,save            :: istep_rk,icont, iter
  real,save               :: cput_drdt, cput_uvwetc

  public  :: runge_kutta_init,istep_rk,cput_drdt,cput_uvwetc,time_integration

contains

  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine runge_kutta_init

    implicit none

    real              :: ak,bk,ck,dk,ek

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

      h_rk(1)=dt/2.d0;  rk_coef(1)=ak
      h_rk(2)=dt/2.d0;  rk_coef(2)=bk
      h_rk(3)=dt;       rk_coef(3)=ck
      h_rk(4)=dt/6.d0;  rk_coef(4)=dk
    endif

    !--------------------------------------------------------------------------!
    ! modify coefficents according to SSP
    !--------------------------------------------------------------------------!
    if (rk_nstages==4.and.rk_order==2.and.lSSPRK) then
      h_rk(1)=dt/3.d0;  rk_coef(1)=1.d0
      h_rk(2)=dt/3.d0;  rk_coef(2)=1.d0
      h_rk(3)=dt/3.d0;  rk_coef(3)=1.d0
      h_rk(4)=dt/4.d0;  rk_coef(4)=1.d0
      dts(1)=0.d0;         dte(1)=dt/3.d0
      dts(2)=dt/3.d0;      dte(2)=dt*2.d0/3.d0
      dts(3)=dt*2.d0/3.d0; dte(3)=dt
      dts(4)=dt;           dte(4)=dt
    endif


    !--------------------------------------------------------------------------!
    ! local time stepping
    !--------------------------------------------------------------------------!
    if (lsteady) then
      allocate(dt_local(ncells))
      dt_local(:) = dt
    endif

    !-- open residual file
    if (lwrite_resid) then
      open(iunit_res, file=trim(file_res))
      write(iunit_res,'(a)') 'variables = "iteration" "|<greek>\r</greek>|<sub>2</sub>",  "|<greek>\r</greek>u|<sub>2</sub>", "|<greek>\r</greek>v|<sub>2</sub>", "|<greek>\r</greek>E|<sub>2</sub>"'
    endif

    iter =0;
    icont=0;
    cput_drdt=0.d0
    cput_uvwetc=0.d0

    return
  end subroutine runge_kutta_init


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine time_integration (t1,ntimes_sub)
    implicit none
    real,intent(in)     :: t1
    integer,intent(in)  :: ntimes_sub

    if (lSSPRK) then
      if (lsteady) then
        call time_integ_SSPRK_steady(t1,ntimes_sub)
      else
        call time_integ_SSPRK(t1,ntimes_sub)
      endif
    else
      if (lsteady) then
        call time_integ_RK_steady(t1,ntimes_sub)
      else
        call time_integ_RK(t1,ntimes_sub)
      endif
    endif

    return
  end subroutine time_integration


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine time_integ_RK (t1,ntimes_sub)

    implicit none
    real,intent(in)     :: t1
    integer,intent(in)  :: ntimes_sub
    integer             :: rk,istep
    real                :: t,told,t12,tnew,tstart(rk_nstages),tend(rk_nstages),cput1,cput2
    real,allocatable    :: fcvar(:,:), cvar0(:,:), diff(:)

    !-- allocation
    allocate(fcvar(nvar,ncells), cvar0(nvar,ncells))


    !-- loop to integrate from t1 to t2
    do istep=1,ntimes_sub
      t = t1 + dble(istep-1)*dt

      told = t
      t12  = told + 0.5d0*dt
      tnew = told + dt

      icont=icont+1

      if (rk_nstages==4) then
        tstart(1)=told; tstart(2)=t12;  tstart(3)=t12;  tstart(4)=tnew;
          tend(1)=t12;    tend(2)=t12;    tend(3)=tnew;   tend(4)=tnew;
      endif

      fcvar = 0.d0
      cvar0 = cvar

      !-- loop over R-K stages
      do rk=1,rk_nstages

        call compute_residual (tstart(rk));

        fcvar(1:nvar, 1:ncells) = fcvar(1:nvar, 1:ncells) + rk_coef(rk) * resid(1:nvar, 1:ncells)

        if (rk<rk_nstages) then
          cvar(1:nvar, 1:ncells) = cvar0(1:nvar, 1:ncells) + h_rk(rk) * resid(1:nvar, 1:ncells)
        elseif (rk==rk_nstages) then
          cvar(1:nvar, 1:ncells) = cvar0(1:nvar, 1:ncells) + h_rk(rk) * fcvar(1:nvar, 1:ncells)
        endif

      enddo

      if (lvortex) call error_isentropic_vortex(tnew)

      !-- write out residual
      if (lwrite_resid) then
        allocate(diff(ncells))
        diff(1:ncells) = abs(cvar(ir,1:ncells) - cvar0(ir,1:ncells));
        write(iunit_res,'(i0,1x,e16.8)',advance='no') icont,sqrt(sum(diff**2)/dble(ncells))

        diff(1:ncells) = abs(cvar(iu,1:ncells) - cvar0(iu,1:ncells));
        write(iunit_res,'(e16.8)',advance='no')             sqrt(sum(diff**2)/dble(ncells))

        diff(1:ncells) = abs(cvar(iv,1:ncells) - cvar0(iv,1:ncells));
        write(iunit_res,'(e16.8)',advance='no')             sqrt(sum(diff**2)/dble(ncells))

        diff(1:ncells) = abs(cvar(4 ,1:ncells) - cvar0(4 ,1:ncells));
        write(iunit_res,'(e16.8)')                          sqrt(sum(diff**2)/dble(ncells))

        deallocate(diff)
      endif
    enddo

    !-- deallocate
    deallocate(fcvar, cvar0)

    return
  end subroutine time_integ_RK


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine time_integ_SSPRK(t1,ntimes_sub)

    implicit none
    real,intent(in)     :: t1
    integer,intent(in)  :: ntimes_sub
    integer             :: rk,istep
    real                :: t,told,tstart,tend,cput1,cput2
    real,allocatable    :: fcvar(:,:), cvar0(:,:), diff(:)

    !-- allocation
    allocate(fcvar(nvar,ncells), cvar0(nvar,ncells))

    !-- loop to integrate from t1 to t2
    do istep=1,ntimes_sub
      t = t1 + dble(istep-1)*dt
      told = t
      icont=icont+1

      fcvar = 0.d0
      cvar0 = cvar

      !-- loop over R-K stages
      do rk=1,rk_nstages
        tstart= told + dts(rk)
        tend  = told + dte(rk)

        call compute_residual(tstart)

         cvar(1:nvar, 1:ncells) = cvar0(1:nvar, 1:ncells) + h_rk(rk) * ( rk_coef(rk)*resid(1:nvar, 1:ncells) + fcvar(1:nvar, 1:ncells) )
        fcvar(1:nvar, 1:ncells) = fcvar(1:nvar, 1:ncells) + resid(1:nvar, 1:ncells)
      enddo

      if (lvortex) call error_isentropic_vortex(tend)

      !-- write out residual
      if (lwrite_resid) then
        allocate(diff(ncells))
        diff(1:ncells) = abs(cvar(ir,1:ncells) - cvar0(ir,1:ncells));
        write(iunit_res,'(i0,1x,e16.8)',advance='no') icont,sqrt(sum(diff**2)/dble(ncells))

        diff(1:ncells) = abs(cvar(iu,1:ncells) - cvar0(iu,1:ncells));
        write(iunit_res,'(e16.8)',advance='no')             sqrt(sum(diff**2)/dble(ncells))

        diff(1:ncells) = abs(cvar(iv,1:ncells) - cvar0(iv,1:ncells));
        write(iunit_res,'(e16.8)',advance='no')             sqrt(sum(diff**2)/dble(ncells))

        diff(1:ncells) = abs(cvar(4 ,1:ncells) - cvar0(4 ,1:ncells));
        write(iunit_res,'(e16.8)')                          sqrt(sum(diff**2)/dble(ncells))

        deallocate(diff)
      endif
    enddo

    !-- deallocate
    deallocate(fcvar, cvar0)

    return
  end subroutine time_integ_SSPRK


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine time_integ_RK_steady (t1,ntimes_sub)

    implicit none
    real,intent(in)     :: t1
    integer,intent(in)  :: ntimes_sub
    integer             :: rk,istep,ic, ivar
    real                :: t,told,t12,tnew,tstart(rk_nstages),tend(rk_nstages),cput1,cput2, cst
    real,allocatable    :: fcvar(:,:), cvar0(:,:), diff(:)

    !-- allocation
    allocate(fcvar(nvar,ncells), cvar0(nvar,ncells), diff(ncells))


    !-- loop to integrate from t1 to t2
    do istep=1,ntimes_sub
      t = t1 + dble(istep-1)*dt

      told = t
      t12  = told + 0.5d0*dt
      tnew = told + dt

      icont=icont+1

      if (rk_nstages==4) then
        tstart(1)=told; tstart(2)=t12;  tstart(3)=t12;  tstart(4)=tnew;
          tend(1)=t12;    tend(2)=t12;    tend(3)=tnew;   tend(4)=tnew;
      endif

      fcvar = 0.d0
      cvar0 = cvar

      iter = iter + 1

      !-- loop over R-K stages
      do rk=1,rk_nstages

        call compute_residual (tstart(rk));
        if (rk==1) call compute_local_time

        fcvar(1:nvar, 1:ncells) = fcvar(1:nvar, 1:ncells) + rk_coef(rk) * resid(1:nvar, 1:ncells)

        cst=1.d0/2.d0
        if (rk==3) cst=1.d0
        if (rk==4) cst=1.d0/6.d0

        if (rk<rk_nstages) then
          do ic=1,ncells
            cvar(1:nvar, ic) = cvar0(1:nvar, ic) + dt_local(ic)*cst * resid(1:nvar, ic)
          enddo
        elseif (rk==rk_nstages) then
          do ic=1,ncells
            cvar(1:nvar, ic) = cvar0(1:nvar, ic) + dt_local(ic)*cst * fcvar(1:nvar, ic)
          enddo
        endif
      enddo

      if (lvortex) call error_isentropic_vortex(tnew)


      !-- write out residual
      if (lwrite_resid) then
        allocate(diff(ncells))
        diff(1:ncells) = abs(cvar(ir,1:ncells) - cvar0(ir,1:ncells));
        write(iunit_res,'(i0,1x,e16.8)',advance='no') iter, sqrt(sum(diff**2)/dble(ncells))

        diff(1:ncells) = abs(cvar(iu,1:ncells) - cvar0(iu,1:ncells));
        write(iunit_res,'(e16.8)',advance='no')             sqrt(sum(diff**2)/dble(ncells))

        diff(1:ncells) = abs(cvar(iv,1:ncells) - cvar0(iv,1:ncells));
        write(iunit_res,'(e16.8)',advance='no')             sqrt(sum(diff**2)/dble(ncells))

        diff(1:ncells) = abs(cvar(4 ,1:ncells) - cvar0(4 ,1:ncells));
        write(iunit_res,'(e16.8)')                          sqrt(sum(diff**2)/dble(ncells))

        deallocate(diff)
      endif
    enddo


    !-- deallocate
    deallocate(fcvar, cvar0)

    return
  end subroutine time_integ_RK_steady


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine time_integ_SSPRK_steady(t1,ntimes_sub)

    implicit none
    real,intent(in)     :: t1
    integer,intent(in)  :: ntimes_sub
    integer             :: rk,istep, ic
    real                :: t,told,tstart,tend,cput1,cput2, cst
    real,allocatable    :: fcvar(:,:), cvar0(:,:), diff(:)

    !-- allocation
    allocate(fcvar(nvar,ncells), cvar0(nvar,ncells))

    !-- loop to integrate from t1 to t2
    do istep=1,ntimes_sub
      t = t1 + dble(istep-1)*dt
      told = t
      icont=icont+1

      fcvar = 0.d0
      cvar0 = cvar

      iter = iter + 1

      !-- loop over R-K stages
      do rk=1,rk_nstages
        tstart= told + dts(rk)
        tend  = told + dte(rk)

        call compute_residual(tstart)
        if (rk==1) call compute_local_time

        cst=1.d0/3.d0
        if (rk==4) cst=1.d0/4.d0
        do ic=1,ncells
          cvar(1:nvar, ic) = cvar0(1:nvar, ic) + dt_local(ic)*cst * ( rk_coef(rk)*resid(1:nvar, ic) + fcvar(1:nvar, ic) )
        enddo

        fcvar(1:nvar, 1:ncells) = fcvar(1:nvar, 1:ncells) + resid(1:nvar, 1:ncells)

      enddo

      if (lvortex) call error_isentropic_vortex(tend)


      !-- write out residual
      if (lwrite_resid) then
        allocate(diff(ncells))
        diff(1:ncells) = abs(cvar(ir,1:ncells) - cvar0(ir,1:ncells));
        write(iunit_res,'(i0,1x,e16.8)',advance='no') iter, sqrt(sum(diff**2)/dble(ncells))

        diff(1:ncells) = abs(cvar(iu,1:ncells) - cvar0(iu,1:ncells));
        write(iunit_res,'(e16.8)',advance='no')             sqrt(sum(diff**2)/dble(ncells))

        diff(1:ncells) = abs(cvar(iv,1:ncells) - cvar0(iv,1:ncells));
        write(iunit_res,'(e16.8)',advance='no')             sqrt(sum(diff**2)/dble(ncells))

        diff(1:ncells) = abs(cvar(4 ,1:ncells) - cvar0(4 ,1:ncells));
        write(iunit_res,'(e16.8)')                          sqrt(sum(diff**2)/dble(ncells))

        deallocate(diff)
      endif
    enddo


    !-- deallocate
    deallocate(fcvar, cvar0)

    return
  end subroutine time_integ_SSPRK_steady


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine compute_local_time

    implicit none
    integer           :: ic

    dt_local(:) = dt

    !-- compute local time step: dt = CFL*volume/sum(0.5*max_wave_speed*face_area).
    cell_loop : do ic = 1, ncells
      dt_local(ic) = cfl_user * cell(ic)%vol/( 0.5d0*ws_nrml(ic) )
    end do cell_loop

    return
  end subroutine compute_local_time

end module runge_kutta
