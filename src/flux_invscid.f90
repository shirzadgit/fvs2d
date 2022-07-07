module flux_invscid

  use mainparam,  only  : nvar
  use input,      only  : lflux_inviscd_roe, gamma

  implicit none

  private :: flux_invscid_roe


contains

  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine compute_flux_invscid (pvarL,pvarR, nx,ny,  flux, ws_max)

    implicit none
    real,intent(in)   :: pvarL(nvar), pvarR(nvar), nx,ny
    real,intent(out)  :: flux(nvar), ws_max

    if (lflux_inviscd_roe) then
      call flux_invscid_roe (pvarL,pvarR, nx,ny,  flux, ws_max)
    endif

    return
  end subroutine compute_flux_invscid

  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine flux_invscid_roe (pvarL,pvarR, nx,ny,  flux, ws_max)

    implicit none
    real,intent(in)   :: pvarL(nvar), pvarR(nvar), nx,ny
    real,intent(out)  :: flux(nvar), ws_max

    integer           :: i,j
    real              :: tx,ty
    real              :: rhoL,uL,vL,pL, aL,HL, unL,utL, kL
    real              :: rhoR,uR,vR,pR, aR,HR, unR,utR, kR
    real              :: RT,rho,u,v,H,a,un,ut
    real              :: drho,dp,dun,dut,tke
    real              :: Ldu(nvar), ws(nvar), dws(nvar), Rv(nvar,nvar), diss(nvar)
    real              :: fL(nvar),fR(nvar)

    flux(1:nvar) = 0.d0
    ws_max = 0.d0

    !-- tangent vector
    tx = -ny
    ty =  nx

    !--  Left state                   !--  Right state
    rhoL = pvarL(1);                  rhoR = pvarR(1)
      uL = pvarL(2);                    uR = pvarR(2)
      vL = pvarL(3);                    vR = pvarR(3)
      pL = pvarL(4);                    pR = pvarR(4)
     unL = uL*nx+vL*ny;                unR = uR*nx+vR*ny
     utL = uL*tx+vL*ty;                utR = uR*tx+vR*ty
      aL = sqrt(gamma*pL/rhoL);         aR = sqrt(gamma*pR/rhoR)
      kL = 0.5d0*(uL*uL+vL*vL);         kR = 0.5d0*(uR*uR+vR*vR)
      HL = aL*aL/(gamma-1.d0) + kL;     HR = aR*aR/(gamma-1.d0) + kR

    !-- compute Roe averages
      RT = sqrt(rhoR/rhoL)
     rho = RT*rhoL
       u = (uL+RT*uR)/(1.d0+RT)
       v = (vL+RT*vR)/(1.d0+RT)
       H = (HL+RT*HR)/(1.d0+RT)
       a = sqrt( (gamma-1.d0)*(H-0.5d0*(u*u+v*v)) )
      un = u*nx+v*ny
      ut = u*tx+v*ty

    !-- Wave strengths
    drho = rhoR - rhoL
      dp =   pR - pL
     dun =  unR - unL
     dut =  utR - utL

    LdU(1) = (dp - rho*a*dun )/(2.d0*a*a)
    LdU(2) = rho*dut
    LdU(3) = drho - dp/(a*a)
    LdU(4) = (dp + rho*a*dun )/(2.d0*a*a)

    !-- Wave Speed
    ws(1) = abs(un-a)
    ws(2) = abs(un)
    ws(3) = abs(un)
    ws(4) = abs(un+a)

    !-- Harten's entropy fix, ref: J. Comp. Phys., vol 49, pp:357-393, 1983
    dws(1) = 1.d0/5.d0
    if ( ws(1) < dws(1) ) ws(1) = 0.5d0 * ( ws(1)*ws(1)/dws(1)+dws(1) )
    dws(4) = 1.d0/5.d0
    if ( ws(4) < dws(4) ) ws(4) = 0.5d0 * ( ws(4)*ws(4)/dws(4)+dws(4) )

    !-- Right Eigenvectors
    tke=0.5d0*(u*u+v*v)
    Rv(1,1) = 1.d0;       Rv(1,2) = 0.d0;   Rv(1,3) = 1.d0;   Rv(1,4) = 1.d0
    Rv(2,1) = u - a*nx;   Rv(2,2) = tx;     Rv(2,3) = u;      Rv(2,4) = u + a*nx
    Rv(3,1) = v - a*ny;   Rv(3,2) = ty;     Rv(3,3) = v;      Rv(3,4) = v + a*ny
    Rv(4,1) = H - un*a;   Rv(4,2) = ut;     Rv(4,3) = tke;    Rv(4,4) = H + un*a

    !-- Dissipation Term
    diss(1:4) = 0.d0
    do i=1,4
      do j=1,4
        diss(i) = diss(i) + ws(j)*LdU(j)*Rv(i,j)
      end do
    end do

    !-- Compute left & right flux.
    fL(1) = rhoL*unL
    fL(2) = rhoL*unL * uL + pL*nx
    fL(3) = rhoL*unL * vL + pL*ny
    fL(4) = rhoL*unL * HL

    fR(1) = rhoR*unR
    fR(2) = rhoR*unR * uR + pR*nx
    fR(3) = rhoR*unR * vR + pR*ny
    fR(4) = rhoR*unR * HR

    !-- Roe's flux
    flux = 0.5d0 * (fL + fR - 1.d0*diss)

    !-- Normal max wave speed/2
    ws_max = 0.5d0*(abs(un) + a)

    return
  end subroutine flux_invscid_roe


end module flux_invscid
