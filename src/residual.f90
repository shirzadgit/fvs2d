module residual

  use mainparam,  only  : nvar
  use data_grid
  use data_sol,   only  : pvar
  use flux_invscid
  use gradient

  implicit none

  type face_data
    real,dimension(:),pointer :: f
  end type

  type(face_data),dimension(:,:),pointer    :: recnst

  private   :: recnst

contains


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine compute_residual (resid)

    implicit none
    real,intent(out)  :: resid(nvar,ncells)
    integer           :: ic,ie,ivar,je, ieL,ieR, icL,icR, k, idum(4), iloc(1)
    real              :: nxf,nyf,af, flux(nvar), pfL(nvar),pfR(nvar)
    character         :: numeric_frecnst_type*124
    real,allocatable  :: grad(:,:,:)

    ! ivar=1, density
    ! ivar=2, u-velocity
    ! ivar=3, v-velocity
    ! ivar=4, pressure

    numeric_frecnst_type = 'linear_extrapolation'

    resid(:,:) = 0.d0

    !--------------------------------------------------------------------------!
    ! allocation
    !--------------------------------------------------------------------------!
    allocate( recnst(nvar, ncells) )
    do ivar=1,nvar
      do ic=1,ncells
        allocate( recnst(ivar,ic)%f(cell(ic)%nvrt) )
      enddo
    enddo


    !--------------------------------------------------------------------------!
    ! compute gradient of primative variables at cell centers
    !--------------------------------------------------------------------------!
    allocate(grad(nvar,ncells,2))
    do ivar=1,nvar
      call gradient_cellcntr(pvar(ivar,1:ncells), grad(ivar,1:ncells,1:2));
    enddo


    !--------------------------------------------------------------------------!
    ! re-construct face-value of primative variables
    !--------------------------------------------------------------------------!
    if (trim(numeric_frecnst_type)=='linear_extrapolation') then

      !-- linear extrapolation: u(f)=u(i) + <grad(u), r>
      do ivar=1,nvar
        do ic=1,ncells
          do ie=1,cell(ic)%nvrt
            recnst(ivar,ic)%f(ie) = pvar(ivar,ic) + cell(ic)%pos2edg(ie,1) * grad(ivar,ic,1) + cell(ic)%pos2edg(ie,2) * grad(ivar,ic,2)
          enddo
        enddo
      enddo


      !-- apply limiter


    else
      write(*,*) 'only linear_extrapolation implemented so far'
      stop 'error'
    endif


    !--------------------------------------------------------------------------!
    ! compute left and right fluxes based on re-constructed face-value
    !--------------------------------------------------------------------------!
    do icL=1,ncells
      do ieL=1,cell(icL)%nvrt

        je=cell(icL)%edge(ieL)

        af =edge(je)%area
        nxf=edge(je)%nx * cell(icL)%nrmlsign(ieL)
        nyf=edge(je)%ny * cell(icL)%nrmlsign(ieL)

        !-- find cell and edge number of right state
        icR=cell(icL)%nghbre(ieL)
        if (icR>0) then
          idum(:)= 1
          do k=1,cell(icR)%nvrt
            idum(k)=abs( je - cell(icR)%edge(k) )
          enddo
          iloc = minloc(idum)
          ieR=maxval(iloc)
          if (cell(icR)%edge(ieR)/=je) write(*,*) 'no way'
        else
          icR=icL
          ieR=ieL
        endif

        !-- construct left and right state on each edge
        pfL(1:nvar) = recnst(1:nvar,icL)%f(ieL)
        pfR(1:nvar) = recnst(1:nvar,icR)%f(ieR)

        !-- compute flux
        call flux_invscid_roe (pfL, pfR, nxf,nyf,  flux)

        resid(1:nvar,icL) = resid(1:nvar,icL) + flux(1:nvar)*af
      enddo
    enddo

    return
  end subroutine compute_residual


end module residual
