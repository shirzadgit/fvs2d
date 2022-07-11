module initialize

  use input,        only  : ntstart, dt, lvortex
  use data_grid,    only  : ncells, cell, bndry, nbndries
  use data_solution
  use mms,          only  : mms_sol, compute_isentropic_vortex

  implicit none


contains


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine initialize_solution

    implicit none
    integer     :: ic, ib, i
    real        :: t_initial, pvar_vortex(nvar)


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

      ! do ib=1,nbndries
      !   if (trim(bndry(ib)%type) == 'slip_wall') then
      !     do i=1,bndry(ib)%ncells
      !       ic=bndry(ib)%cell(i)
      !       pvar(iu,ic) =   u_inf*0.7
      !     enddo
      !   endif
      ! enddo

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

    elseif (ntstart < 1) then
      !-- initialize with manufactured solution
      pvar = mms_sol

    endif


    !-- compute conservative variables.
    call pvar2cvar

    return
  end subroutine initialize_solution


end module initialize
