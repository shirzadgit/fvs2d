module gradient_lsq

  use input,  only  : grad_cellcntr_lsq_pow, lgrad_lsq_fn1, lgrad_lsq_fn2, lgrad_lsq_nn
  use grid_procs
  use interpolation

  implicit none

  private
  integer,allocatable,save  :: cell_neighbr_lsq_ptr(:,:)
  real,allocatable,save     :: grad_lsq_coef(:,:,:),lsq_w(:,:)
  real,save                 :: lsq_p

  private :: init,setup
  public  :: grad_lsq_init, grad_lsq

contains

  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grad_lsq_init
    implicit  none

    lsq_p = grad_cellcntr_lsq_pow

    !--------------------------------------------------------------------------!
    ! setup gradient opertor at cell centers: Green-Gauss cell-based (ggcb)
    !--------------------------------------------------------------------------!
    call init


    !--------------------------------------------------------------------------!
    ! setup gradient opertor at cell centers: Green-Gauss node-based (ggnb)
    !--------------------------------------------------------------------------!
    call setup


    return
  end subroutine grad_lsq_init


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine init
    implicit  none

    integer               :: ndim,i,j,nt, ie,ie1,ic,ic1,ic2, iv1,iv2
    integer               :: k,jc,jc1,je,je1,je2,jv1,jv2, kv1,kv2
    real                  :: r(2),dum
    integer,allocatable   :: tmpi2d(:,:)
    type(kdtree2_result),allocatable	:: fnearest(:)
    integer                           :: num_nearests

    allocate(cell_neighbr_lsq_ptr(num_cells,num_vert_max))

    allocate(tmpi2d(num_cells,num_vert_max))

    tmpi2d(:,:) = cell_neighbr_ptr(:,:)

    do ic=1,num_cells
      if (cell_neighbr_num(ic)==1) then
        do ie=1,num_vert_cell(ic)
          ic1=cell_neighbr_ptr(ic,ie)
          ie1=cell2edge(ic,ie)
          iv1=edge2node(ie1,1)
          iv2=edge2node(ie1,2)

          if (ic1==ic) then
            do je=1,num_vert_cell(ic)
              je1=cell2edge(ic,je)
              jc1=cell_neighbr_ptr(ic,je)
              if (je1/=ie1 .and. jc1/=ic) then
                !write(*,*) ic,ie,je,jc1
                jv1=edge2node(je1,1)
                jv2=edge2node(je1,2)
                if (jv1==iv1.or.jv1==iv2.or.jv2==iv1.or.jv2==iv2)  then
                  do k=1,num_vert_cell(jc1)
                    kv1=edge2node(cell2edge(jc1,k),1)
                    kv2=edge2node(cell2edge(jc1,k),2)
                    if ( (kv1-iv1)*(kv1-iv2)*(kv2-iv1)*(kv2-iv2)==0 ) then
                      !write(*,*)  ic,jc1,cell_neighbr_ptr(jc1,k) !ic,ie,je,jc1,
                      if (cell_neighbr_ptr(jc1,k)/=jc1.and.cell_neighbr_ptr(jc1,k)/=ic1) tmpi2d(ic,ie)=cell_neighbr_ptr(jc1,k)
                    endif
                  enddo
                endif
              endif
            enddo
          endif
        enddo
      endif
    enddo


    num_nearests=8
    allocate(fnearest(num_nearests))
    do ic=1,num_cells
      r(1)=cell_center(ic,1)
      r(2)=cell_center(ic,2)

      if (cell_neighbr_num(ic)==2) then
        do ie=1,num_vert_cell(ic)
          if (cell_neighbr_ptr(ic,ie)==ic) ie1=ie
        enddo

        call kdtree2_n_nearest(tp=tree_cellcntr,qv=r,nn=num_nearests,results=fnearest)

        outer: do i=1,num_nearests
          ic1=fnearest(i)%idx
          dum=1.d0
          do ie=1,num_vert_cell(ic)
            dum=dum*max(0,abs(cell_neighbr_ptr(ic,ie)-ic1))
          enddo
          if (dum>0.d0) then
            tmpi2d(ic,ie1)=ic1
            !write(*,*) ic,ic1,i,cell_neighbr_ptr(ic,ie)
            exit outer
          endif
        enddo outer
      endif
    enddo



    ! do ic=1,num_cells
    !   if (cell_neighbr_num(ic)==1) then
    !     !write(*,*) ic
    !     do ie=1,num_vert_cell(ic)
    !       !write(*,*) ie,tmpi2d(ic,ie)
    !     enddo
    !   endif
    ! enddo
    ! write(*,*)
    ! do ic=1,num_cells
    !   if (cell_neighbr_num(ic)==2) then
    !     write(*,*) ic
    !     do ie=1,num_vert_cell(ic)
    !       write(*,*) ie,tmpi2d(ic,ie)
    !     enddo
    !   endif
    ! enddo

    cell_neighbr_lsq_ptr(:,:)=tmpi2d(:,:)

    nt=0
    open(100,file='lsq_2neighbor.plt')
    do ic=1,num_cells
      if (cell_neighbr_num(ic)==2 .and. nt<20) then
        nt=nt+1
        !write(100,'(a)') 'TITLE ="grid"'
        write(100,'(a)') 'VARIABLES ="x", "y"'
        write(100,'(a)') 'zone i=4, j=1, f=point'
        write(100,*) cell_center(ic,1),cell_center(ic,2)
        do ie=1,3 !num_vert_cell(ic)
          ic1=cell_neighbr_lsq_ptr(ic,ie)
          write(100,*) cell_center(ic1,1),cell_center(ic1,2)
        enddo
      endif
    enddo
    close(100)

    return
  end subroutine init


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine setup
    implicit  none

    integer     :: i,j,k, ic,ic1,ic2, ie,ie1,ie2, iv,iv1,iv2
    real        :: d(num_vert_max,2), dt(2,num_vert_max), g(2,2), gi(2,2), det,xc,yc,xc1,yc1,w,dis


    allocate(grad_lsq_coef(num_cells,2,num_vert_max), lsq_w(num_cells,num_vert_max))
    grad_lsq_coef(:,:,:) = 0.d0
    lsq_w(:,:)=0.d0

    !--------------------------------------------------------------------------!
    ! Unweighted least squares
    !--------------------------------------------------------------------------!
    grad_lsq_coef(:,:,:)=0.d0
    do ic=1,num_cells
      xc=cell_center(ic,1)
      yc=cell_center(ic,2)
      do ie=1,num_vert_cell(ic)
        ie1=cell2edge(ic,ie)
        ic1=cell_neighbr_lsq_ptr(ic,ie)

        xc1=cell_center(ic1,1)
        yc1=cell_center(ic1,2)

        dis=dsqrt( (xc1-xc)**2 + (yc1-yc)**2 );
        w=0.d0
        if (dis==0.d0) then
          lsq_w(ic,ie)=0.d0
        else
          w=1.d0/dis**lsq_p
          lsq_w(ic,ie)=w;
        endif

        d(ie,1)=w*(xc1-xc);
        d(ie,2)=w*(yc1-yc);

        dt(1,ie)=d(ie,1)
        dt(2,ie)=d(ie,2)
      enddo

      !-- compute G= d^T * d
      g(:,:)=0.d0
      do i=1,2
        do j=1,2
          do ie=1,num_vert_cell(ic)
            g(i,j)= g(i,j) + dt(i,ie)*d(ie,j)
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
        do ie=1,num_vert_cell(ic)
          do k=1,2
            grad_lsq_coef(ic,i,ie)= grad_lsq_coef(ic,i,ie) + gi(i,k)*dt(k,ie)
          enddo
        enddo
      enddo
    enddo

    do ic=1,num_cells
      !if (cell_neighbr_num(ic)==1) grad_lsq_coef(ic,:,:)=0.d0
    enddo


    return
  end subroutine setup


  !============================================================================!
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
  !============================================================================!
  subroutine grad_lsq (fc,dfc)
    implicit  none

    real,intent(in)   :: fc(num_cells)
    real,intent(out)  :: dfc(num_cells,2)
    integer           :: i,j, ic,ic1, ie,ie1, in,iv,iv1,iv2, nt
    real              :: af,nxf,nyf
    real,allocatable  :: fv(:)

    do ic=1,num_cells
      do ie=1,num_vert_cell(ic)
        ic1=cell_neighbr_lsq_ptr(ic,ie)
        dfc(ic,1) = dfc(ic,1) + grad_lsq_coef(ic,1,ie)*(fc(ic1)-fc(ic))*lsq_w(ic,ie)
        dfc(ic,2) = dfc(ic,2) + grad_lsq_coef(ic,2,ie)*(fc(ic1)-fc(ic))*lsq_w(ic,ie)
      enddo
    end do

    return
  end subroutine grad_lsq

end module gradient_lsq
