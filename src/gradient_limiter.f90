module gradient_limiter

  use mainparam,      only  : nvar
  use input,          only  : lgrad_limiter, limiter_type, lface_reconst_upwind1st
  use data_solution,  only  : grad, phi_lim, pvar
  use data_grid
  use gradient_lsq

  implicit none


contains


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine compute_gradient_limiter

    implicit none
    integer       :: i,j, ie,je, ic,jc, ivar
    real          :: pmin,pmax, dmax,dmin, diff, phi_edge(4)
    real          :: xc,yc, xf,yf, pf
    real,allocatable  :: phi(:,:)

    !--------------------------------------------------------------------------!
    ! if 1storder upwind --> phi_limiter = 0
    !--------------------------------------------------------------------------!
    if (lface_reconst_upwind1st) then
      phi_lim = 0.d0
      return
    endif


    !-- If no limiter
    if (.not.lgrad_limiter) then
      phi_lim = 1.d0
      return
    endif

    !--------------------------------------------------------------------------!
    ! setup
    !--------------------------------------------------------------------------!
    allocate(phi(nvar,ncells))
    do ivar=1,nvar
      do ic=1,ncells
        xc = cell(ic)%x
        yc = cell(ic)%y

        pmin = pvar(ivar,ic)
        pmax = pvar(ivar,ic)
        do i=1,lsq(ic)%ncells
          jc=lsq(ic)%cell(i)
          pmin = min(pmin, pvar(ivar,jc))
          pmax = max(pmax, pvar(ivar,jc))
        enddo

        phi_edge(:) = 100.d0
        do ie=1,cell(ic)%nvrt
          je=cell(ic)%edge(ie)

          !-- edge properties
          xf = edge(je)%x
          yf = edge(je)%y
          pf = pvar(ivar,ic) + (xf-xc) * grad(ivar,ic,1) + (yf-yc) * grad(ivar,ic,2)

          ! dum1(1:nvar) = pvar(ivar,ic) - pvar(1:nvar,icL)
          ! dumL(1:nvar) = (xf-xc) * grad(ivar,ic,1) + (yf-yc) * grad(ivar,ic,2)
          ! pfL(1:nvar) = pvar(1:nvar,icL) + umuscl_cst/2.d0 * dum1(1:nvar) + (1.d0-umuscl_cst)*dumL(1:nvar)

          dmax = pmax - pvar(ivar,ic)
          dmin = pmin - pvar(ivar,ic)
          diff = pf   - pvar(ivar,ic)
          if (diff>0.d0) then
            phi_edge(ie) = limiter(dmax,diff,cell(ic)%vol)
          elseif (diff<0.d0) then
            phi_edge(ie) = limiter(dmin,diff,cell(ic)%vol)
          else
            phi_edge(ie) = 1.d0
          endif
        enddo

        phi(ivar,ic) = min(1.d0, minval(phi_edge))
      end do
    enddo

    do ic=1,ncells
      phi_lim(ic) = minval(phi(1:nvar,ic))
    enddo

    return
  end subroutine compute_gradient_limiter


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  pure function limiter(a, b, vol)
    implicit none

    real, intent(in)  :: a, b, vol
    real              :: limiter
    real              :: ap, eps2, h, cst
    real,parameter    :: pi=acos(-1.d0)


    select case (limiter_type)
    case ('venk')
      ap = 5.0d0               ! adjustable parameter K
      h = 2.d0*sqrt(vol/pi)    ! h = 2*diamater
      eps2 = (ap*h)**3         ! eps2 = eps^2
      limiter = ( (a**2 + eps2) + 2.d0*b*a )/(a**2 + 2.d0*b**2 + a*b + eps2)

    case ('barth')
      limiter = min(1.d0,a/b)

    case ('albada')
      h = 2.d0*sqrt(vol/pi)
      eps2 = (0.3d0*h)**3
      cst = 0.5d0*( sign(1.0,a*b) + 1.d0 )
      limiter = ( (b**2 + eps2)*a + (a**2 + eps2)*b )/(a**2 + b**2 + 2.d0*eps2)
      limiter = limiter/(b+eps2)

    end select

    !if (limiter<0.d0) write(*,*) 'limiter<0!'

    return
  end function limiter

end module gradient_limiter
