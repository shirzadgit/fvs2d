module residual

  use mainparam,      only  : nvar
  use input
  use data_grid
  use data_solution,  only  : pvar, cvar, resid, grad, cvar2pvar, pvar2cvar, pvar_inf, ir,iu,iv,ip
  use flux_invscid
  use gradient
  use gradient_limiter
  use mms,            only  : mms_compute_euler2d, vortex_inf, compute_isentropic_vortex

  implicit none

  real,save   :: un_max

contains

  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine compute_residual (time)

    implicit none
    real,intent(in)   :: time
    integer           :: i,ic,ie,ivar,je, ieL,ieR, icL,icR, k, idum(4), iloc(1), ib
    real              :: flux(nvar), pfL(nvar),pfR(nvar), ws_max, cfl_max
    real              :: xcL,ycL, xcR,ycR, xc,yc, xf,yf, nxf,nyf, af
    real,allocatable  :: ws_nrml(:), un(:)

    resid(:,:) = 0.d0

    allocate( ws_nrml(ncells) )

    !--------------------------------------------------------------------------!
    ! convert conservatve variables to primative variables
    !--------------------------------------------------------------------------!
    call cvar2pvar


    !--------------------------------------------------------------------------!
    ! compute gradient of primative variables at cell centers
    !--------------------------------------------------------------------------!
    call compute_gradient_cellcntr


    !--------------------------------------------------------------------------!
    ! compute gradient limiter
    !--------------------------------------------------------------------------!
    call compute_gradient_limiter


    !--------------------------------------------------------------------------!
    ! compute left and right inviscid fluxes on interior edges
    !--------------------------------------------------------------------------!
    interior_edges: do i=1,nedges_intr
      ie = edge_intr(i)

      !-- edge properties
       xf = edge(ie)%x
       yf = edge(ie)%y
       af = edge(ie)%area
      nxf = edge(ie)%nx
      nyf = edge(ie)%ny

      !-- left cell index and coordinates
      icL = edge(ie)%c1
      xcL = cell(icL)%x
      ycL = cell(icL)%y

      !-- right cell index and coordinates
      icR = edge(ie)%c2
      xcR = cell(icR)%x
      ycR = cell(icR)%y

      !-- re-construct face-value of primative variables
      pfL(1:nvar) = pvar(1:nvar,icL) + phi_lim(icL) * ((xf-xcL) * grad(1:nvar,icL,1) + (yf-ycL) * grad(1:nvar,icL,2))
      pfR(1:nvar) = pvar(1:nvar,icR) + phi_lim(icR) * ((xf-xcR) * grad(1:nvar,icR,1) + (yf-ycR) * grad(1:nvar,icR,2))

      !-- compute inviscid flux
      call compute_flux_invscid (pfL, pfR, nxf,nyf,  flux, ws_max)

      resid(1:nvar,icL) = resid(1:nvar,icL) + flux(1:nvar)*af
      resid(1:nvar,icR) = resid(1:nvar,icR) - flux(1:nvar)*af

      ws_nrml(icL) = ws_nrml(icL) + ws_max*af
      ws_nrml(icR) = ws_nrml(icR) + ws_max*af
    enddo interior_edges


    !--------------------------------------------------------------------------!
    ! compute left and right fluxes for boundary edges/faces
    !--------------------------------------------------------------------------!
    allocate(un(bndry(2)%nedges)); un=0.d0
    nboundaries: do ib=1,nbndries
      boundary_edges: do i=1,bndry(ib)%nedges

        !-- get global edge number
        ie = bndry(ib)%edge(i)

        !-- edge properties
         xf = edge(ie)%x
         yf = edge(ie)%y
         af = edge(ie)%area
        nxf = edge(ie)%nx
        nyf = edge(ie)%ny

        !-- left cell index and coordinates
        icL= bndry(ib)%cell(i)
        xcL= cell(icL)%x
        ycL= cell(icL)%y

        !-- compute left state
        pfL(1:nvar) = pvar(1:nvar,icL) + phi_lim(icL) * ((xf-xcL) * grad(1:nvar,icL,1) + (yf-ycL) * grad(1:nvar,icL,2))

        !-- compute right state
        call bc_flux (time,xf,yf,nxf,nyf,bndry(ib)%type,pfL, pfR)

        !-- compute un
        if (ib==2) un(i) = pfL(iu)*nxf + pfL(iv)*nyf

        !-- compute inviscid flux
        call compute_flux_invscid (pfL, pfR, nxf,nyf,  flux, ws_max)

        resid(1:nvar,icL) = resid(1:nvar,icL) + flux(1:nvar)*af

        ws_nrml(icL) = ws_nrml(icL) + ws_max*af
      enddo boundary_edges
    enddo nboundaries


    ! dQ/dt = -residual
    do ic=1,ncells
      resid(1:nvar,ic) = -resid(1:nvar,ic)/cell(ic)%vol
    enddo


    cfl_max = 0.d0
    do ic=1,ncells
      cfl_max = max(cfl_max, dt/cell(ic)%vol *0.5d0 * ws_nrml(ic))
    enddo
    !write(*,*) 'cfl_max:',cfl_max
    un_max = maxval(abs(un))
    deallocate(un)

    return
  end subroutine compute_residual


  subroutine bc_flux (time,x,y,nx,ny,bc_type,pfL, pfR)
    implicit none

    real,intent(in)       :: time,x,y,nx,ny,pfL(1:nvar)
    character,intent(in)  :: bc_type*80
    real,intent(out)      :: pfR(nvar)
    real                  :: rhs(nvar), un
    real,parameter        :: slipwall_velcnst = 2.d0

    pfR = 0.d0

    select case(trim(bc_type))
      case('freestream')
        pfR(ir) = pvar_inf(ir)
        pfR(iu) = pvar_inf(iu)
        pfR(iv) = pvar_inf(iv)
        pfR(ip) = pvar_inf(ip)

        !pfR(iu) = u_old + incr
        !u_old = pfR(iu)

      case('slip_wall')
        pfR = pfL

        un = pfL(iu)*nx + pfL(iv)*ny
        pfR(iu) = pfL(iu) - slipwall_velcnst*un*nx
        pfR(iv) = pfL(iv) - slipwall_velcnst*un*ny

      case('solid_wall')


      case('dirichlet')
        if (lvortex) then
          !pfR(1:nvar) = vortex_inf(1:nvar)
          call compute_isentropic_vortex (time,x,y,pfR)
        else
          call mms_compute_euler2d (x,y,pfR,rhs)
        endif

      case default
        write(*,*) "Boundary condition=",trim(bc_type),"  not implemented."
        stop
    end select

    return
  end subroutine bc_flux


end module residual
