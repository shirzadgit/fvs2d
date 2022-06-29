module mms_euler2d

  use data_type, only   : ncells, cell

  implicit none

  private

  type source_type
    real          :: r,u,v,e
  end type source_type

  type(source_type),dimension(:),pointer  :: source

  real,parameter  :: gamma=1.4d0
  real,save       :: cr0,crs,crx,cry
  real,save       :: cu0,cus,cux,cuy
  real,save       :: cv0,cvs,cvx,cvy
  real,save       :: cp0,cps,cpx,cpy
  real,parameter  :: pi=acos(-1.d0)

   public   :: compute_source_euler, source

contains


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine compute_source_euler

    implicit none
    integer     :: ic
    real        :: xc,yc
    real        :: r,u,v,p, rx,ux,vx,px, ry,uy,vy,py, rH,rHx,rHy

    allocate(source(ncells))

    !--------------------------------------------------------------------------!
    ! define constants
    !--------------------------------------------------------------------------!
    !-- Density    = cr0 + crs*sin(crx*x+cry*y)
    cr0 =  1.12
    crs =  0.15
    crx =  3.12*pi
    cry =  2.92*pi

    !-- u-velocity = cu0 + cus*sin(cux*x+cuy*y)
    cu0 =  1.32
    cus =  0.06
    cux =  2.09*pi
    cuy =  3.12*pi

    !-- v-velocity = cv0 + cvs*sin(cvx*x+cvy*y)
    cv0 =  1.18
    cvs =  0.03
    cvx =  2.15*pi
    cvy =  3.32*pi

    !-- pressure   = cp0 + cps*sin(cpx*x+cpy*y)
    cp0 =  1.62
    cps =  0.31
    cpx =  3.79*pi
    cpy =  2.98*pi

    !--------------------------------------------------------------------------!
    ! compute source terms
    !--------------------------------------------------------------------------!
    do ic=1,ncells
      xc=cell(ic)%x
      yc=cell(ic)%y

      !-- compute primative variables
      r = manufactured_sol(cr0,crs,crx,cry, 0,0, xc,yc)
      u = manufactured_sol(cu0,cus,cux,cuy, 0,0, xc,yc)
      v = manufactured_sol(cv0,cvs,cvx,cvy, 0,0, xc,yc)
      p = manufactured_sol(cp0,cps,cpx,cpy, 0,0, xc,yc)

      !-- compute primative variables derivatives w.r.t x-direction
      rx = manufactured_sol(cr0,crs,crx,cry, 1,0, xc,yc)
      ux = manufactured_sol(cu0,cus,cux,cuy, 1,0, xc,yc)
      vx = manufactured_sol(cv0,cvs,cvx,cvy, 1,0, xc,yc)
      px = manufactured_sol(cp0,cps,cpx,cpy, 1,0, xc,yc)

      !-- compute primative variables derivatives w.r.t y-direction
      ry = manufactured_sol(cr0,crs,crx,cry, 0,1, xc,yc)
      uy = manufactured_sol(cu0,cus,cux,cuy, 0,1, xc,yc)
      vy = manufactured_sol(cv0,cvs,cvx,cvy, 0,1, xc,yc)
      py = manufactured_sol(cp0,cps,cpx,cpy, 0,1, xc,yc)

      !-- compute rho*H=gamma/(gamma-1)*p + rho*(u^2+v^2)/2
      rH = gamma/(gamma-1.d0)*p + r*u*u/2.d0 + r*v*v/2.d0
      rHx= gamma/(gamma-1.d0)*px + rx*(u*u + v*v)/2.d0 + r*(u*ux + v*vx)
      rHy= gamma/(gamma-1.d0)*py + ry*(u*u + v*v)/2.d0 + r*(u*uy + v*vy)

      !-- continuity: source%r = d(ru)/dx + d(rv)/dy
      source(ic)%r =  rx*u + u*rx + ry*v + r*vy

      !-- x-momentum: source%u = d(ruu)/dx + d(ruv)/dy + dp/dx
      source(ic)%u = rx*u*u + 2.d0*r*u*ux + ry*u*v + r*uy*v + r*u*vy + px

      !-- y-momentum: source%v = d(ruv)/dx + d(rvv)/dy + dp/dy
      source(ic)%v = rx*u*v + r*ux*v + r*u*vx + ry*v*v + 2.d0*r*v*vy + py

      !--     energy: source%e =  d(rho*u*H)/dx + (rho*v*H)/dy
      source(ic)%e = u*rHx + ux*rH + v*rHy + vy*rH

    enddo

    return
  end subroutine compute_source_euler



  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  function manufactured_sol(a0,as,ax,ay, nx,ny,x,y) result(fval)

    implicit none
    real, intent(in)    :: a0, as, ax, ay, x, y
    integer, intent(in) :: nx, ny
    real                :: fval

    if (nx < 0 .or. ny < 0) then
      write(*,*) " Invalid input: nx & ny must be greater or equal to zero!"
      stop 'stopped at: mod:mms_euler2d, sub: manufactured_sol'
    endif


    if ( nx+ny == 0 ) then
      fval = a0 + as*sin(ax*x + ay*y)

    elseif ( mod(nx+ny,2) == 0 ) then
      fval = - (ax**nx * ay**ny)*as*sin(ax*x + ay*y)
      if ( mod(nx+ny,4)   == 0 ) fval = -fval

    else
      fval = (ax**nx * ay**ny)*as*cos(ax*x + ay*y)
      if ( mod(nx+ny+1,4) == 0 ) fval = -fval
    endif

    return
  end function manufactured_sol


end module mms_euler2d
