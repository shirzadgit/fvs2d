module data_solution

  use mainparam,  only  : nvar, grid_dim
  use input,      only  : gamma, mach_inf, aoa_inf_deg, rey
  use data_grid,  only  : ncells, nedges, nnodes

  implicit none


  ! ---------------------------------------------------------------------------!
  ! cvar = conservative variables (rho, rho*u, rho*v, rho*E)
  ! pvar =    primitive variables (rho,     u,     v,     p)
  ! ---------------------------------------------------------------------------!

  !-- variables, gradient
  real, allocatable, dimension(:,:)   :: pvar     ! primitive variables at cell centers
  real, allocatable, dimension(:,:)   :: cvar     ! conservatve variables at cell centers
  real, allocatable, dimension(:,:,:) :: grad     ! gradients of pvar at cell centers
  real, allocatable, dimension(:,:)   :: resid    ! residual vector at cell center
  real, allocatable, dimension(:)     :: phi_lim  ! limiter
  real, allocatable, dimension(:)     :: ws_nrml  ! normal wave speed


  !-- variables locations in pvar
  integer,parameter         :: ir = 1             ! pvar(ir) = density
  integer,parameter         :: iu = 2             ! pvar(iu) = u-velocity
  integer,parameter         :: iv = 3             ! pvar(iv) = v-velocity
  integer,parameter         :: ip = grid_dim + 2  ! pvar(ip) = pressure


  !-- free-stream values
  real,save                 :: rho_inf, u_inf, v_inf, p_inf, pvar_inf(nvar)


contains


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine data_solution_init
    implicit none


    !-- allocation
    allocate( pvar(nvar,ncells))
    allocate( cvar(nvar,ncells))
    allocate( grad(nvar,ncells,grid_dim))
    allocate(resid(nvar,ncells))
    allocate(phi_lim(ncells))
    allocate(ws_nrml(ncells))


    !-- set free stream values based on the input Mach number.
    rho_inf = 1.d0
      u_inf = mach_inf*cosd(aoa_inf_deg) !aoa converted from degree to radian
      v_inf = mach_inf*sind(aoa_inf_deg) !aoa converted from degree to radian
      p_inf = 1.d0/gamma

    pvar_inf(ir) = rho_inf
    pvar_inf(iu) = u_inf
    pvar_inf(iv) = v_inf
    pvar_inf(ip) = p_inf

    return
  end subroutine data_solution_init


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine cvar2pvar

    implicit none
    integer   :: ic

    pvar = 0.d0
    do ic=1,ncells
      pvar(ir,ic) = cvar(1,ic)
      pvar(iu,ic) = cvar(2,ic)/cvar(1,ic)
      pvar(iv,ic) = cvar(3,ic)/cvar(1,ic)
      pvar(ip,ic) = (gamma-1.d0)*( cvar(4,ic) - 0.5d0*pvar(1,ic)*(pvar(2,ic)**2 + pvar(3,ic)**2) )
    enddo

    return
  end subroutine cvar2pvar


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine pvar2cvar

    implicit none
    integer   :: ic

    cvar = 0.d0
    do ic=1,ncells
      cvar(1,ic) = pvar(ir,ic)
      cvar(2,ic) = pvar(ir,ic)*pvar(iu,ic)
      cvar(3,ic) = pvar(ir,ic)*pvar(iv,ic)
      cvar(4,ic) = pvar(ip,ic)/(gamma-1.d0) + 0.5d0*pvar(ir,ic)*(pvar(iu,ic)**2+pvar(iv,ic)**2)
    enddo

    return
  end subroutine pvar2cvar


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  function prim2cnsr(w) result(u)

    implicit none
    real, dimension(4), intent(in) :: w ! input
    real, dimension(4)             :: u !output

    u(1) = w(ir)
    u(2) = w(ir)*w(iu)
    u(3) = w(ir)*w(iv)
    u(4) = w(ip)/(gamma-1.d0)+0.5d0*w(ir)*(w(iu)*w(iu)+w(iv)*w(iv))

    return
  end function prim2cnsr


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine data_sol_close
    implicit none

    deallocate(pvar)
    deallocate(cvar)
    deallocate(grad)
    deallocate(resid)
    deallocate(phi_lim)

    return
  end subroutine data_sol_close

end module data_solution
