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


  !-- variables locations in pvar
  integer,parameter         :: ir = 1             ! pvar(ir) = density
  integer,parameter         :: iu = 2             ! pvar(iu) = u-velocity
  integer,parameter         :: iv = 3             ! pvar(iv) = v-velocity
  integer,parameter         :: ip = grid_dim + 2  ! pvar(ip) = pressure


  !-- free-stream values
  real,save                 :: rho_inf, u_inf, v_inf, p_inf


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


    !-- set free stream values based on the input Mach number.
    rho_inf = 1.d0
      u_inf = mach_inf*cosd(aoa_inf_deg) !aoa converted from degree to radian
      v_inf = mach_inf*sind(aoa_inf_deg) !aoa converted from degree to radian
      p_inf = 1.d0/gamma

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
      pvar(ip,ic) = (gamma-1.d0)*( cvar(4,ic) - 0.5d0*cvar(1,ic)*(cvar(2,ic)**2 + cvar(3,ic)**2) )
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