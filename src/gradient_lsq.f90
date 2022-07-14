module gradient_lsq

  use mainparam,      only  : nvar
  use input,          only  : grad_cellcntr_lsq_pow, lgrad_lsq_fn, lgrad_lsq_nn
  use data_solution,  only  : pvar, grad
  use data_grid
  use grid_procs
  use interpolation

  implicit none

  private
  real,save                 :: lsq_p
  logical,save              :: grad_lsq_lverify

  type lsq_type
    integer                       :: ncells
    integer,dimension(:),pointer  :: cell
    real,dimension(:),pointer     :: w
    real,dimension(:,:),pointer   :: coef
  end type lsq_type

  type(lsq_type),dimension(:),pointer  :: lsq

  private :: setup_fn, setup_nn
  public  :: grad_lsq_init, grad_lsq, grad_lsq_1var !grad_lsq_fn, grad_lsq_nn,
  public  :: lsq

contains

  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grad_lsq_init
    implicit  none

    lsq_p = grad_cellcntr_lsq_pow

    !--------------------------------------------------------------------------!
    ! allocation
    !--------------------------------------------------------------------------!
    allocate(lsq(ncells))


    !--------------------------------------------------------------------------!
    ! setup gradient opertor at cell centers: Green-Gauss node-based (ggnb)
    !--------------------------------------------------------------------------!
    if (lgrad_lsq_fn) then
      !--  setup gradient opertor at cell centers: Green-Gauss node-based (ggnb)
      call setup_fn

    elseif (lgrad_lsq_nn) then
      !-- setup gradient opertor at cell centers: Green-Gauss node-based (ggnb)
      call setup_nn
    endif


    !--------------------------------------------------------------------------!
    ! verify LSQ coefficents via linear function
    !--------------------------------------------------------------------------!
    call grad_lsq_verify

    return
  end subroutine grad_lsq_init


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine setup_fn
    implicit  none

    real                              :: r(2),dum
    type(kdtree2_result),allocatable	:: fnearest(:)
    integer                           :: num_nearests
    integer                           :: i,j,k, ic,jc,in, izb,nt,tmpi(4)
    real                              :: g(2,2),gi(2,2), det,xc,yc,xc1,yc1,w,dis
    real,allocatable                  :: d(:,:),dt(:,:)


    !--------------------------------------------------------------------------!
    ! find additional cells for boundary cells
    ! each cell should have 3 neighboring cells
    !--------------------------------------------------------------------------!
    num_nearests = 8;
    allocate(fnearest(num_nearests))

    do ic=1,ncells
      allocate(lsq(ic)%cell(cell(ic)%nvrt))

      lsq(ic)%ncells=cell(ic)%nvrt

      lsq(ic)%cell(:)=0
      r(1)=cell(ic)%x
      r(2)=cell(ic)%y
      tmpi(:)=0
      izb=0
      do in=1,cell(ic)%nvrt
        jc=cell(ic)%nghbr(in)
        if (jc==0) then
          izb=izb+1
          tmpi(izb)=in
        else
          lsq(ic)%cell(in)=jc
        endif
      enddo
      if (izb>0 .and. izb<3) then
        call kdtree2_n_nearest(tp=tree_cellcntr,qv=r,nn=num_nearests,results=fnearest)
        nt=0
        outer: do i=1,num_nearests
          jc=fnearest(i)%idx
          if (jc/=ic) then
            dum=1.d0
            do in=1,cell(ic)%nvrt
              dum=dum*max(0,abs(cell(ic)%nghbr(in)-jc))
            enddo

            if (dum>0.d0) then
              nt=nt+1
              if (lsq(ic)%cell(tmpi(nt))>0) write(*,*) 'lsqfn: no way'
              lsq(ic)%cell(tmpi(nt))=jc
            endif
            if (nt==izb) exit outer
          endif
        enddo outer
      elseif (izb>2) then
        write(*,*) 'error in gradient_lsq, sub: setup_fn: #s of edges on the boundary>2!'
        stop 'mod: gradient_lsq, sub: setup_fn'
      endif

    enddo


    !--------------------------------------------------------------------------!
    ! compute LSQ coefficents and weights
    !--------------------------------------------------------------------------!
    cell_loop: do ic=1,ncells
      xc=cell(ic)%x
      yc=cell(ic)%y

      !-- allocation
      allocate(lsq(ic)%w(cell(ic)%nvrt))
      allocate(lsq(ic)%coef(2,cell(ic)%nvrt))
      allocate(d(cell(ic)%nvrt,2), dt(2,cell(ic)%nvrt))

      !-- initialize
      lsq(ic)%w(:)=0.d0
      lsq(ic)%coef(:,:)=0.d0
      d (:,:)=0.d0
      dt(:,:)=0.d0

      !-- loop over cell nodes
      do in=1,cell(ic)%nvrt
        jc=lsq(ic)%cell(in)

        xc1=cell(jc)%x
        yc1=cell(jc)%y

        dis=dsqrt( (xc1-xc)**2 + (yc1-yc)**2 );
        w=0.d0
        if (dis>0.d0) then
          w=1.d0/dis**lsq_p
          lsq(ic)%w(in)=w;
        endif

        d(in,1)=w*(xc1-xc);
        d(in,2)=w*(yc1-yc);

        dt(1,in)=d(in,1)
        dt(2,in)=d(in,2)
      enddo

      !-- compute G= d^T * d
      g(:,:)=0.d0
      do i=1,2
        do j=1,2
          do in=1,cell(ic)%nvrt
            g(i,j)= g(i,j) + dt(i,in)*d(in,j)
          enddo
        enddo
      enddo

      !-- compute G^-1
      gi(:,:)=0.d0
      det=g(1,1)*g(2,2)-g(1,2)*g(2,1)
      gi(1,1) = 1.d0/det*g(2,2)
      gi(2,2) = 1.d0/det*g(1,1)
      gi(1,2) =-1.d0/det*g(1,2)
      gi(2,1) =-1.d0/det*g(2,1)

      !-- compute G^-1 * d^T
      do i=1,2
        do in=1,cell(ic)%nvrt
          do k=1,2
            lsq(ic)%coef(i,in)= lsq(ic)%coef(i,in) + gi(i,k)*dt(k,in)
          enddo
        enddo
      enddo

      deallocate(d,dt)
    enddo cell_loop

    do ic=1,ncells
      !if (cell_neighbr_num(ic)==1) grad_lsq_coef(ic,:,:)=0.d0
    enddo

    return
  end subroutine setup_fn


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine setup_nn
    implicit  none

    real,allocatable    :: d(:,:), dt(:,:), wrk(:,:)
    real                :: g(2,2), gi(2,2), det,xc,yc,xc1,yc1,w,dis
    integer             :: i,j,k, ic,jc,in, iv
    integer             :: min_tmp, max_tmp, nt
    integer,allocatable :: tmp(:), tmp_unq(:)

    !--------------------------------------------------------------------------!
    ! find additional cells for boundary cells
    ! each cell should have 3 neighboring cells
    !--------------------------------------------------------------------------!
    do ic=1,ncells

      !-- step 1: get total number of surrounding cells; there will dublicates
      nt=0
      do in=1,cell(ic)%nvrt
        iv=cell(ic)%node(in)
        do i=1,node(iv)%ncells
          jc=node(iv)%cell(i)
          if (ic/=jc) nt=nt+1
        end do
      end do

      !-- step 2: allocate and assign the surrounding cells index
      allocate(tmp(nt), tmp_unq(nt))
      nt=0
      do in=1,cell(ic)%nvrt
        iv=cell(ic)%node(in)
        do i=1,node(iv)%ncells
          jc=node(iv)%cell(i)
          if (ic/=jc) then
            nt=nt+1
            tmp(nt)=jc
          endif
        end do
      end do

      !--- step 3: remove duplicates
      nt=0
      min_tmp = minval(tmp)-1
      max_tmp = maxval(tmp)
      do while (min_tmp<max_tmp)
        nt = nt+1
        min_tmp = minval(tmp, mask=tmp>min_tmp)
        tmp_unq(nt) = min_tmp
      enddo
      tmp(1:nt)=tmp_unq(1:nt)

      lsq(ic)%ncells=nt
      allocate( lsq(ic)%cell(nt) )
      do i=1,nt
        lsq(ic)%cell(i) = tmp(i)
      enddo

      deallocate(tmp,tmp_unq)
    end do


    !--------------------------------------------------------------------------!
    ! Least squares
    !--------------------------------------------------------------------------!
    nt=0
    do ic=1,ncells
      xc=cell(ic)%x
      yc=cell(ic)%y

      !--allocation
      allocate( lsq(ic)%w(lsq(ic)%ncells) )
      allocate( lsq(ic)%coef(2,lsq(ic)%ncells) )
      allocate( d(lsq(ic)%ncells,2), dt(2,lsq(ic)%ncells), wrk(2,lsq(ic)%ncells) )

      !--initialize
      lsq(ic)%w(:)=0.d0
      lsq(ic)%coef(:,:)=0.d0
      d (:,:)=0.d0
      dt(:,:)=0.d0

      do i=1,lsq(ic)%ncells
        jc=lsq(ic)%cell(i)

        xc1=cell(jc)%x
        yc1=cell(jc)%y

        dis=dsqrt( (xc1-xc)**2 + (yc1-yc)**2 );
        w=0.d0
        if (dis>0.d0) then
          w=1.d0/dis**lsq_p
          lsq(ic)%w(i)=w;
        endif

        d(i,1)=w*(xc1-xc);
        d(i,2)=w*(yc1-yc);

        dt(1,i)=d(i,1)
        dt(2,i)=d(i,2)
      enddo

      !-- compute G= d^T * d
      g(:,:)=0.d0
      do i=1,2
        do j=1,2
          do k=1,lsq(ic)%ncells
            g(i,j)= g(i,j) + dt(i,k)*d(k,j)
          enddo
        enddo
      enddo

      !-- compute G^-1
      gi(:,:)=0.d0
      det=g(1,1)*g(2,2)-g(1,2)*g(2,1)
      gi(1,1) = 1.d0/det*g(2,2)
      gi(2,2) = 1.d0/det*g(1,1)
      gi(1,2) =-1.d0/det*g(1,2)
      gi(2,1) =-1.d0/det*g(2,1)

      !-- compute G^-1 * d^T
      wrk(:,:)=0.d0
      do i=1,2
        do j=1,lsq(ic)%ncells
          do k=1,2
            wrk(i,j)= wrk(i,j) + gi(i,k)*dt(k,j)
          enddo
        enddo
      enddo


      do i=1,lsq(ic)%ncells
        lsq(ic)%coef(1:2,i) = wrk(1:2,i)
      enddo

      deallocate(d,dt,wrk)

    enddo

    ! open(100,file='grid_lsqnn.plt')
    ! do ic=1000,1100 !,ncells
    !   write(100,'(a)') 'VARIABLES ="X", "Y"'
    !   write(100,'(a)') 'zone i=1, j=1, f=point'
    !   write(100,*) cell(ic)%x,cell(ic)%y
    !   write(100,'(a)') 'VARIABLES ="X", "Y"'
    !   write(100,'(a,i0,a)') 'zone i=',lsq(ic)%ncells, ', j=1, f=point'
    !   do in=1,lsq(ic)%ncells
    !     jc=lsq(ic)%cell(in)
    !     write(100,*) cell(jc)%x,cell(jc)%y
    !   enddo
    ! end do
    ! close(100)

    return
  end subroutine setup_nn


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grad_lsq

    implicit  none
    integer           :: ic,jc,in,ivar

    grad = 0.d0

    do ivar=1,nvar
      do ic=1,ncells
        do in=1,lsq(ic)%ncells
          jc=lsq(ic)%cell(in)
          grad(ivar,ic,1) = grad(ivar,ic,1) + lsq(ic)%coef(1,in)*(pvar(ivar,jc)-pvar(ivar,ic))*lsq(ic)%w(in)
          grad(ivar,ic,2) = grad(ivar,ic,2) + lsq(ic)%coef(2,in)*(pvar(ivar,jc)-pvar(ivar,ic))*lsq(ic)%w(in)
        enddo
      end do
    enddo

    return
  end subroutine grad_lsq


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grad_lsq_1var (fc,dfc)
    implicit  none

    real,intent(in)   :: fc(ncells)
    real,intent(out)  :: dfc(ncells,2)
    integer           :: ic,jc,in

    do ic=1,ncells
      do in=1,lsq(ic)%ncells !cell(ic)%nvrt
        jc=lsq(ic)%cell(in)
        dfc(ic,1) = dfc(ic,1) + lsq(ic)%coef(1,in)*(fc(jc)-fc(ic))*lsq(ic)%w(in)
        dfc(ic,2) = dfc(ic,2) + lsq(ic)%coef(2,in)*(fc(jc)-fc(ic))*lsq(ic)%w(in)
      enddo
    end do

    return
  end subroutine grad_lsq_1var

  ! !============================================================================!
  ! !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  ! !============================================================================!
  ! subroutine grad_lsq_fn (fc,dfc)
  !   implicit  none
  !
  !   real,intent(in)   :: fc(ncells)
  !   real,intent(out)  :: dfc(ncells,2)
  !   integer           :: ic,jc,in
  !
  !   do ic=1,ncells
  !     do in=1,lsq(ic)%ncells !cell(ic)%nvrt
  !       jc=lsq(ic)%cell(in)
  !       dfc(ic,1) = dfc(ic,1) + lsq(ic)%coef(1,in)*(fc(jc)-fc(ic))*lsq(ic)%w(in)
  !       dfc(ic,2) = dfc(ic,2) + lsq(ic)%coef(2,in)*(fc(jc)-fc(ic))*lsq(ic)%w(in)
  !     enddo
  !   end do
  !
  !   return
  ! end subroutine grad_lsq_fn
  !
  !
  ! !============================================================================!
  ! !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  ! !============================================================================!
  ! subroutine grad_lsq_nn (fc,dfc)
  !   implicit  none
  !
  !   real,intent(in)   :: fc(ncells)
  !   real,intent(out)  :: dfc(ncells,2)
  !   integer           :: ic,jc,in
  !
  !   do ic=1,ncells
  !     do in=1,lsq(ic)%ncells
  !       jc=lsq(ic)%cell(in)
  !       dfc(ic,1) = dfc(ic,1) + lsq(ic)%coef(1,in)*(fc(jc)-fc(ic))*lsq(ic)%w(in)
  !       dfc(ic,2) = dfc(ic,2) + lsq(ic)%coef(2,in)*(fc(jc)-fc(ic))*lsq(ic)%w(in)
  !     enddo
  !   end do
  !
  !   return
  ! end subroutine grad_lsq_nn


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grad_lsq_verify
    implicit  none

    integer           :: ic,jc,in
    real              :: xc,yc, xnb,ynb, diff,dfx,dfy
    real,parameter    :: cstx=2.d0, csty=1.d0

    grad_lsq_lverify = .true.

    do ic=1,ncells
      xc=cell(ic)%x
      yc=cell(ic)%y
      dfx=0.d0
      dfy=0.d0
      do in=1,lsq(ic)%ncells
        jc=lsq(ic)%cell(in)
        xnb=cell(jc)%x
        ynb=cell(jc)%y
        diff = csty*ynb + cstx*xnb - (csty*yc + cstx*xc)
        dfx = dfx + lsq(ic)%coef(1,in)*(diff)*lsq(ic)%w(in)
        dfy = dfy + lsq(ic)%coef(2,in)*(diff)*lsq(ic)%w(in)
      enddo

      if ( abs(dfx-cstx) > 1.0e-10 .or. abs(dfy-csty) > 1.0e-10) then
        write(*,*) " df/dx = ", dfx, " exact df/dx =",cstx
        write(*,*) " df/dy = ", dfy, " exact df/dy =",csty
        grad_lsq_lverify = .false.
      endif
    end do

    if (.not.grad_lsq_lverify) then
      write(*,*) " LSQ coefficients are not correct"
      write(*,*) " Error in mod: gradient_lsq, sub: grad_lsq_verify"
      stop
    !else
    !  write(*,*) " Verified: LSQ coefficients are exact for a linear function."
    endif

    return
  end subroutine grad_lsq_verify

end module gradient_lsq
