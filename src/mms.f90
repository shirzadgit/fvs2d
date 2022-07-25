module mms

  use mainparam,      only  : nvar, file_vortex, iunit_vortex
  use input,          only  : gamma, lvortex, dt
  use data_grid !,      only  : ncells, cell, node
  use data_solution,  only  : ir,iu,iv,ip

  implicit none

  private

  real,allocatable,dimension(:,:) :: mms_sol, mms_source

  real,save       :: cr0,crs,crx,cry
  real,save       :: cu0,cus,cux,cuy
  real,save       :: cv0,cvs,cvx,cvy
  real,save       :: cp0,cps,cpx,cpy
  real,parameter  :: pi=acos(-1.d0)

  real,save       :: vortex_pos(2), vortex_kappa
  real,save       :: vortex_inf(nvar)

  public  :: mms_init, mms_compute_euler2d, manufactured_sol, mms_source, mms_sol
  public  :: compute_isentropic_vortex, vortex_inf, error_isentropic_vortex

contains


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine mms_init

    implicit none
    integer     :: ic
    real        :: xc,yc
    real        :: r,u,v,p, rx,ux,vx,px, ry,uy,vy,py, rH,rHx,rHy
    real        :: xmax, xmin, ymax, ymin, ax,ay
    logical     :: linputfile


    !--------------------------------------------------------------------------!
    ! isentropic vortex
    !--------------------------------------------------------------------------!
    if (lvortex) then
      inquire(file=trim(file_vortex),exist=linputfile);
      if (.not.linputfile) then
        write(*,*)
        write(*,*) 'linputfile:',linputfile
        write(*,*) 'cannot find "'//trim(file_vortex)//'" file!'
        write(*,*) 'error in --> mod:mms, sub:mms_init'
        stop 'program stopped at "mms_init"'
      endif
      open(iunit_vortex,file=trim(file_vortex),status='unknown',form='formatted')
      read(iunit_vortex,*) vortex_pos(1), vortex_pos(2)
      read(iunit_vortex,*) vortex_kappa
      read(iunit_vortex,*) vortex_inf(ir)
      read(iunit_vortex,*) vortex_inf(iu)
      read(iunit_vortex,*) vortex_inf(iv)
      read(iunit_vortex,*) vortex_inf(ip)
      close(iunit_vortex)
      return
    endif


    allocate(mms_source(nvar,ncells), mms_sol(nvar,ncells))

    !--------------------------------------------------------------------------!
    ! define constants
    !--------------------------------------------------------------------------!
    xmin=minval(node(:)%x)
    xmax=maxval(node(:)%x)
    ax=8.d0*pi/(xmax-xmin)

    ymin=minval(node(:)%y)
    ymax=maxval(node(:)%y)
    ay=7.d0*pi/(ymax-ymin)

    !-- Density    = cr0 + crs*sin(crx*x+cry*y)
    cr0 =  1.12
    crs =  0.15
    crx =  3.12*pi  !8.d0*pi/(xmax-xmin) !
    cry =  2.92*pi  !7.d0*pi/(xmax-xmin) !

    !-- u-velocity = cu0 + cus*sin(cux*x+cuy*y)
    cu0 =  1.32
    cus =  0.06
    cux =  2.09*pi !6.5d0*pi/(xmax-xmin) !
    cuy =  3.12*pi !8.0d0*pi/(xmax-xmin) !

    !-- v-velocity = cv0 + cvs*sin(cvx*x+cvy*y)
    cv0 =  1.18
    cvs =  0.03
    cvx =  2.15*pi !7.2d0*pi/(xmax-xmin) !
    cvy =  3.32*pi !8.2d0*pi/(xmax-xmin) !

    !-- pressure   = cp0 + cps*sin(cpx*x+cpy*y)
    cp0 =  1.62
    cps =  0.31
    cpx =  3.79*pi !7.8d0*pi/(xmax-xmin) !
    cpy =  2.98*pi !6.8d0*pi/(xmax-xmin) !

    !--------------------------------------------------------------------------!
    ! compute manufactured solution and source terms for 2-D Euler
    !--------------------------------------------------------------------------!
    do ic=1,ncells
      xc = cell(ic)%x
      yc = cell(ic)%y

      !-- compute manufactured solution and source terms
      call mms_compute_euler2d (cell(ic)%x, cell(ic)%y, mms_sol(1:nvar,ic), mms_source(1:nvar,ic))

      !-- integrate over volume
      mms_source(1:nvar,ic) = mms_source(1:nvar,ic) !!!* cell(ic)%vol
    enddo

    return
  end subroutine mms_init


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine mms_compute_euler2d (xc,yc,sol,rhs)

    implicit none
    real,intent(in)   :: xc,yc
    real,intent(out)  :: sol(nvar), rhs(nvar)

    integer     :: ic
    real        :: r,u,v,p, rx,ux,vx,px, ry,uy,vy,py, rH,rHx,rHy

    sol = 0.d0
    rhs = 0.d0

    !--------------------------------------------------------------------------!
    ! compute manufactured solution and source terms
    !--------------------------------------------------------------------------!
    !-- compute primative variables
    r = manufactured_sol(cr0,crs,crx,cry, 0,0, xc,yc)
    u = manufactured_sol(cu0,cus,cux,cuy, 0,0, xc,yc)
    v = manufactured_sol(cv0,cvs,cvx,cvy, 0,0, xc,yc)
    p = manufactured_sol(cp0,cps,cpx,cpy, 0,0, xc,yc)

    !-- save manufactured solution
    sol(1) = r
    sol(2) = u
    sol(3) = v
    sol(4) = p

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
    rH = gamma/(gamma-1.d0)*p  + r*u*u/2.d0 + r*v*v/2.d0
    rHx= gamma/(gamma-1.d0)*px + rx*(u*u + v*v)/2.d0 + r*(u*ux + v*vx)
    rHy= gamma/(gamma-1.d0)*py + ry*(u*u + v*v)/2.d0 + r*(u*uy + v*vy)

    !-- continuity: source%r = d(rho*u)/dx + d(rho*v)/dy
    rhs(1) = rx*u + u*rx + ry*v + r*vy

    !-- x-momentum: source%u = d(rho*u*u)/dx + d(rho*u*v)/dy + dp/dx
    rhs(2) = rx*u*u + 2.d0*r*u*ux + ry*u*v + r*uy*v + r*u*vy + px

    !-- y-momentum: source%v = d(rho*u*v)/dx + d(rho*v*v)/dy + dp/dy
    rhs(3) = rx*u*v + r*ux*v + r*u*vx + ry*v*v + 2.d0*r*v*vy + py

    !--     energy: source%e =  d(rho*u*H)/dx + (rho*v*H)/dy
    rhs(4) = u*rHx + ux*rH + v*rHy + vy*rH

    return
  end subroutine mms_compute_euler2d


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


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine compute_isentropic_vortex (t,x,y,pvar_vortex)
    implicit none

    real,intent(in)   :: x,y,t
    real,intent(out)  :: pvar_vortex(nvar)
    real              :: rho_inf, u_inf, v_inf, p_inf, T_inf, K
    real              :: xc,yc, dx,dy,r
    real              :: rho, u, v, p, temp

    !--
    K = vortex_kappa

    !-- free stream values
    rho_inf = vortex_inf(ir)
      u_inf = vortex_inf(iu)
      v_inf = vortex_inf(iv)
      p_inf = vortex_inf(ip)
      T_inf = p_inf/rho_inf

    xc = vortex_pos(1) + u_inf*t
    yc = vortex_pos(2) + v_inf*t

    dx = x - xc
    dy = y - yc
     r = sqrt( dx**2 + dy**2 )

    u = u_inf - K/(2.d0*pi) * dy * exp(0.5d0*(1.d0-r**2))
    v = v_inf + K/(2.d0*pi) * dx * exp(0.5d0*(1.d0-r**2))

    temp =   T_inf - (K/(2.d0*pi))**2 *  (gamma-1.d0)/(2.d0*gamma) * exp(1.d0-r**2)
    !   p = ((rho_inf**gamma)/p_inf * temp**gamma)**(1.d0/(gamma-1.d0))
    ! rho = p/temp

    rho = temp**(1.d0/(gamma-1.d0))
      p = rho**gamma

    ! rho = rho_inf * tem**( 1.d0/(gamma-1.d0)) !Density
    !   p =   p_inf * tem**(gamma/(gamma-1.d0)) !Pressure

    pvar_vortex(ir) = rho
    pvar_vortex(iu) = u
    pvar_vortex(iv) = v
    pvar_vortex(ip) = p


    return
  end subroutine compute_isentropic_vortex


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine error_isentropic_vortex(time)
    use data_solution,  only  : cvar
    implicit none

    real,intent(in) :: time
    integer         :: i,ic
    real            :: dr,dru,drv,dre, pv(nvar), rho,ru,rv,re, err_L2, erho(ncells)
    real            :: err_dr_max , err_dr_L1 , err_dr_L2
    real            :: err_dru_max, err_dru_L1, err_dru_L2
    real            :: err_drv_max, err_drv_L1, err_drv_L2
    real            :: err_dre_max, err_dre_L1, err_dre_L2

    if (time==dt) then
      open(iunit_vortex, file='log_vortex_err.plt')
      write(iunit_vortex,'(a)') 'variables = "t", '
      write(iunit_vortex,'(a)') '"<greek>r</greek><sub>max</sub>",  "<greek>r</greek><sub>L1</sub>", "<greek>r</greek><sub>L2</sub>",'
      write(iunit_vortex,'(a)') '"<greek>r</greek>u<sub>max</sub>",  "<greek>r</greek>u<sub>L1</sub>", "<greek>r</greek>u<sub>L2</sub>",'
      write(iunit_vortex,'(a)') '"<greek>r</greek>v<sub>max</sub>",  "<greek>r</greek>v<sub>L1</sub>", "<greek>r</greek>v<sub>L2</sub>",'
      write(iunit_vortex,'(a)') '"<greek>r</greek>E<sub>max</sub>",  "<greek>r</greek>E<sub>L1</sub>", "<greek>r</greek>E<sub>L2</sub>",'
      write(iunit_vortex,'(a)') '"Q<sub>L2</sub>"'

      open(701, file='log_vortex_err_xy.plt')
      write(701,'(a)') 'variables = "t", "x" "y"'
    endif

    err_dr_max = 0.d0
    err_dru_max= 0.d0
    err_drv_max= 0.d0
    err_dre_max= 0.d0

    err_dr_L1 = 0.d0
    err_dru_L1= 0.d0
    err_drv_L1= 0.d0
    err_dre_L1= 0.d0

    err_dr_L2 = 0.d0
    err_dru_L2= 0.d0
    err_drv_L2= 0.d0
    err_dre_L2= 0.d0

    err_L2 = 0.d0

    erho = 0.d0

    do i=1,ncells_intr
      ic=cell_intr(i)

      call compute_isentropic_vortex (time,cell(ic)%x,cell(ic)%y,pv)

      rho = pv(ir);
       ru = pv(ir)*pv(iu)
       rv = pv(ir)*pv(iv)
       re = pv(ip)/(gamma-1.d0) + 0.5d0*rho*(pv(iu)**2+pv(iv)**2)


      dr  = abs(cvar(ir,ic) - rho)
      dru = abs(cvar(iu,ic) - ru)
      drv = abs(cvar(iv,ic) - rv)
      dre = abs(cvar(4 ,ic) - re)
      erho(ic) = dr;

      err_dr_max = max(err_dr_max, dr)
      err_dr_L1  = err_dr_L1 + dr
      err_dr_L2  = err_dr_L2 + dr*dr

      err_dru_max = max(err_dru_max, dru)
      err_dru_L1  = err_dru_L1 + dru
      err_dru_L2  = err_dru_L2 + dru*dru

      err_drv_max = max(err_drv_max, drv)
      err_drv_L1  = err_drv_L1 + drv
      err_drv_L2  = err_drv_L2 + drv*drv

      err_dre_max = max(err_dre_max, dre)
      err_dre_L1  = err_dre_L1 + dre
      err_dre_L2  = err_dre_L2 + dre*dre

      err_L2  = err_L2 + dr**2 + dru**2 + drv**2 + dre**2
    enddo

    ! write(iunit_vortex, '(14(e16.8,1x))') time, err_dr_max,  err_dr_L1 /dble(ncells_intr), sqrt(err_dr_L2) /dble(ncells_intr), &
    !                                             err_dru_max, err_dru_L1/dble(ncells_intr), sqrt(err_dru_L2)/dble(ncells_intr), &
    !                                             err_drv_max, err_drv_L1/dble(ncells_intr), sqrt(err_drv_L2)/dble(ncells_intr), &
    !                                             err_dre_max, err_dre_L1/dble(ncells_intr), sqrt(err_dre_L2)/dble(ncells_intr), &
    !                                             sqrt(err_L2)/dble(ncells_intr)

    write(iunit_vortex, '(14(e16.8,1x))') time, err_dr_max,  err_dr_L1 /dble(ncells_intr), sqrt(err_dr_L2 /dble(ncells_intr)), &
                                                err_dru_max, err_dru_L1/dble(ncells_intr), sqrt(err_dru_L2/dble(ncells_intr)), &
                                                err_drv_max, err_drv_L1/dble(ncells_intr), sqrt(err_drv_L2/dble(ncells_intr)), &
                                                err_dre_max, err_dre_L1/dble(ncells_intr), sqrt(err_dre_L2/dble(ncells_intr)), &
                                                sqrt(err_L2/dble(ncells_intr))

    write(701,*) time,cell(maxloc(erho))%x, cell(maxloc(erho))%y
    return
  end subroutine error_isentropic_vortex


end module mms
