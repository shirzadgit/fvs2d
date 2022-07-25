module initialize

  use input,        only  : ntstart, dt, lvortex
  use data_grid,    only  : nnodes, ncells, cell, bndry, nbndries
  use data_solution
  use mms,          only  : mms_sol, compute_isentropic_vortex
  use ios_unstrc
  use io,           only  : imach_cont, itape_cont, iunit_cont, fbase_cont

  implicit none


contains


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine initialize_solution

    implicit none
    integer               :: ic, ib, i, ivar
    real                  :: t_initial, pvar_vortex(nvar)
    integer               :: nnodes_r, ncells_r, nt_r, np_r, minf_r, ierr
    integer               :: itimes_r(maxtime)
    character(len=maxlen) :: inf_r(maxinf)


    t_initial = dble(ntstart-1)*dt


    !-- set zero
    pvar = 0.d0
    cvar = 0.d0

    if (ntstart == 1) then
      !-- initialize primative variables by the free stream values
      cell_loop1 : do ic = 1, ncells
        pvar(ir,ic) = rho_inf
        pvar(iu,ic) =   u_inf
        pvar(iv,ic) =   v_inf
        pvar(ip,ic) =   p_inf
      enddo cell_loop1

      !-- initialize with isentropic vortex
      if (lvortex) then
        cell_loop2 : do ic = 1, ncells
          call compute_isentropic_vortex (t_initial, cell(ic)%x, cell(ic)%y, pvar_vortex)
          pvar(ir,ic) = pvar_vortex(ir)
          pvar(iu,ic) = pvar_vortex(iu)
          pvar(iv,ic) = pvar_vortex(iv)
          pvar(ip,ic) = pvar_vortex(ip)
        enddo cell_loop2
      endif

    elseif (ntstart > 1) then
      !-- read continuation files

      !-- read control data file
      call readcd (itape_cont, iunit_cont, fbase_cont, imach_cont, &
                   nt_r, ncells_r, nnodes_r, np_r, itimes_r, inf_r, minf_r, ierr)

      !-- check dimension
      if (ncells_r /= ncells .or. nnodes_r /= nnodes) then
        write(*,*) 'dimension between grid and cont files does not match!'
        write(*,*) 'error in mod:initialize, sub: initialize_solution'
        stop 'error'
      elseif (nt_r>1 .or. np_r/=nvar)  then
        write(*,*) 'nt>1 or #s parameters/=nvar!'
        write(*,*) 'error in mod:initialize, sub: initialize_solution'
        stop 'error'
      endif

      !-- read conservatve variables
      do ivar=1,nvar
        call readd(itape_cont, cvar(ivar,1:ncells), 1, ivar, ierr)
      enddo

    elseif (ntstart < 1) then
      !-- initialize with manufactured solution
      pvar = mms_sol

    endif


    !-- compute conservative variables.
    if (ntstart<=1) call pvar2cvar

    return
  end subroutine initialize_solution


end module initialize
