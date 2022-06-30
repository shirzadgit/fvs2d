module flux_invscid

  ! use mainparam
  ! use data_grid
  ! use data_sol

  implicit none

  real,parameter,private  :: gamma=1.4d0



contains


  ! !============================================================================!
  ! !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  ! !============================================================================!
  ! subroutine face_reconstruction
  !
  !   implicit none
  !   integer           :: je, ieL,ieR, icL,icR, k, idum(4), iloc(1)
  !   real              :: nxf,nyf,af, flux(4), pfL(4),pfR(4)
  !   real,allocatable  :: grad(:,:,:)
  !   real,allocatable,dimension(:,:) :: dr,du,dv,dp, wrk
  !
  !   !--------------------------------------------------------------------------!
  !   ! allocation array for face values
  !   !--------------------------------------------------------------------------!
  !   allocate( fv(ncells) )
  !   do ic=1,ncells
  !     allocate(fv(ic)%r(cell(ic)%nvrt))
  !     allocate(fv(ic)%u(cell(ic)%nvrt))
  !     allocate(fv(ic)%v(cell(ic)%nvrt))
  !     allocate(fv(ic)%p(cell(ic)%nvrt))
  !   enddo
  !
  !
  !   !--------------------------------------------------------------------------!
  !   ! compute gradient of primative variables at cell centers
  !   !--------------------------------------------------------------------------!
  !   allocate(grad(nvar,ncells,2)) !,  du(ncells,2), dv(ncells,2), dp(ncells,2),)
  !   ! ivar=1, density
  !   ! ivar=2, u-velocity
  !   ! ivar=3, v-velocity
  !   ! ivar=4, pressure
  !   do ivar=1,nvar
  !     call gradient_cellcntr(pvar(ivar,1:ncells), dr(ivar,1:ncells,1:2));
  !   enddo
  !   ! call gradient_cellcntr(pvar(1:ncells)%r, dr); !--density gradient
  !   ! call gradient_cellcntr(pvar(1:ncells)%u, du); !--u-velocity gradient
  !   ! call gradient_cellcntr(pvar(1:ncells)%v, dv); !--v-velocity gradient
  !   ! call gradient_cellcntr(pvar(1:ncells)%p, dp); !--pressure gradient
  !
  !
  !   !--------------------------------------------------------------------------!
  !   ! re-construct face-value of primative variables
  !   !--------------------------------------------------------------------------!
  !   if (trim(num_reconstruct_type)=='linear_extrapolation') then
  !
  !     !-- linear extrapolation: u(f)=u(i) + <grad(u), r>
  !     do ic=1,ncells
  !       do ie=1,cell(ic)%nvrt
  !         do ivar=1,nvar
  !           fv(ic,ivar)%r(ie) = pvar(ivar,ic) + cell(ic)%pos2edg(ie,1) * grad(ivar,ic,1) + cell(ic)%pos2edg(ie,2) * grad(ivar,ic,2)
  !         enddo
  !       enddo
  !     enddo
  !     ! do ic=1,ncells
  !     !   do ie=1,cell(ic)%nvrt
  !     !     fv(ic)%r(ie) = pvar(ic)%r + cell(ic)%pos2edg(ie,1) * dr(ic,1) + cell(ic)%pos2edg(ie,2) * dr(ic,2)
  !     !     fv(ic)%u(ie) = pvar(ic)%u + cell(ic)%pos2edg(ie,1) * du(ic,1) + cell(ic)%pos2edg(ie,2) * du(ic,2)
  !     !     fv(ic)%v(ie) = pvar(ic)%v + cell(ic)%pos2edg(ie,1) * dv(ic,1) + cell(ic)%pos2edg(ie,2) * dv(ic,2)
  !     !     fv(ic)%p(ie) = pvar(ic)%p + cell(ic)%pos2edg(ie,1) * dp(ic,1) + cell(ic)%pos2edg(ie,2) * dp(ic,2)
  !     !   enddo
  !     ! enddo
  !
  !     !-- apply limiter
  !
  !
  !   !
  !   else
  !     write(*,*) 'only linear_extrapolation implemented so far'
  !     stop 'error'
  !   endif
  !
  !   return
  ! end subroutine face_reconstruction


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine flux_invscid_roe (pvarL, pvarR, nx,ny,  flux) !, wsn)
    implicit none

    real,intent(in)   :: pvarL(4), pvarR(4), nx,ny
    real,intent(out)  :: flux(4)

    real              :: tx,ty
    real              :: rhoL,uL,vL,pL, aL,HL, unL,utL
    real              :: rhoR,uR,vR,pR, aR,HR, unR,utR
    real              :: fL(4),fR(4)

    flux(4) = 0.d0

    !-- tangent vector
    tx = -ny
    ty =  nx

    !Primitive and other variables.
    !--  Left state
    rhoL = pvarL(1)
      uL = pvarL(2)
      vL = pvarL(3)
      pL = pvarL(4)


    unL = uL*nx+vL*ny
    utL = uL*tx+vL*ty
     aL = sqrt(gamma*pL/rhoL)
     HL = aL*aL/(gamma-1.d0) + 0.5d0*(uL*uL+vL*vL)

    !--  Right state
    rhoR = pvarR(1)
      uR = pvarR(2)
      vR = pvarR(3)
      pR = pvarR(4)

    unR = uR*nx+vR*ny
    utR = uR*tx+vR*ty
     aR = sqrt(gamma*pR/rhoR)
     HR = aR*aR/(gamma-1.d0) + 0.5d0*(uR*uR+vR*vR)


    !-- Compute the flux.
    fL(1) = rhoL*unL
    fL(2) = rhoL*unL * uL + pL*nx
    fL(3) = rhoL*unL * vL + pL*ny
    fL(4) = rhoL*unL * HL

    fR(1) = rhoR*unR
    fR(2) = rhoR*unR * uR + pR*nx
    fR(3) = rhoR*unR * vR + pR*ny
    fR(4) = rhoR*unR * HR

    !--
    flux = 0.5d0 * (fL + fR)

    if (maxval(abs(flux))==0.d0) write(*,*) 'r u kidding me'

    return
  end subroutine flux_invscid_roe


end module flux_invscid
