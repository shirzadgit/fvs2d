program ios2tecplot
  use ios_unstrc
  use grid_procs
  use tecplot_unstrc

  implicit none

  integer               :: it, ip, i, j, mt1, mt2, mt3, clen
  integer               :: nnodes_r, ncells_r, nt, np, ierr, itimes(maxtime), minf, imach
  character(len=maxlen) :: varname(maxinf)
  character             :: fbase*124, fname*127, varinfo*100, fnameplt*127, fnameit5*5, fnameit6*6
  real,allocatable      :: var1(:), var(:,:)
  real                  :: sol_time
  integer,parameter     :: itaped=17, iunitd=117
  logical               :: ltime, lGridSolTogether

  write (*,*)
  write (*,*) '**************************************************************'
  write (*,*) ' read unstructred binary file and write it in Tecplot format  '
  write (*,*) '**************************************************************'

  !--
  write(*,*)
  write(*,'(a)',advance='no') ' enter grid file name (*.grid): '
  read (*,*) fname


  !--
  call read_grid2d(fname)


  !--
  write(*,*)
  write(*,'(a)',advance='no') ' ios data file name: '
  read (*,*) fname
  call mkfname(fname, fbase, imach, ierr)
  call readcd (itaped, iunitd, fbase, imach, nt, nnodes_r, ncells_r, np, itimes, varname, minf, ierr)

  !--
  write(*,*)
  write(*,'(a,i0)') ' #s of time-steps: ',nt
  write(*,'(a,i0)') ' #s of parameters: ',np

  !---
  write(*,*)
  write (*,'(a)',advance='no') ' cut out time interval (T/F): '
  read  (*,*) ltime
  mt1=1
  mt2=nt
  mt3=1
  if (ltime) then
    write(*,*); write(*,'(a)',advance='no') ' mt-min,mt-max,mt-skip: '
    read (*,*)  mt1,mt2,mt3
  endif

  !--
  write (*,*)
  write (*,'(a)',advance='no') ' Tecplot output file without extension (.plt will be added): '
  read  (*,*) fname

  !--
  write (*,*)
  write (*,'(a)',advance='no') ' Grid & Solution in the same file (T/F): '
  read  (*,*) lGridSolTogether

  !--
  if(lGridSolTogether) then
      varinfo(1:4) = 'x,y,'
      clen=5
  else
      clen=1
  endif

  do i=1,np
      varname(i) = trim(adjustl(varname(i)))
      call StripSpaces (varname(i))
      j=len(trim(varname(i)))
      if (i<np) then
         varinfo(clen:clen+j)=trim(varname(i))//','
         clen=clen+j+1
     elseif (i==np) then
         varinfo(clen:clen+j-1)=trim(varname(i))
         clen=clen+j
     endif
  enddo


  !-- allocation
  allocate( var1(nnodes), var(nnodes,np) )


  !-- initialize tecplot
  call tecplot_init (ncells_quad)
  imach = 1


  !--
  if (.not.lGridSolTogether) then
    call tecplot_write_grid (nnodes, ncells, grid_dim, imach, trim(trim(fname)//'_grid.plt'), cell2node, node_xy)
    write(*,*)
    write(*,'(a,a)') ' grid file is written in ', trim(trim(fname)//'_grid.plt')
    write(*,'(a,a)') ' writing solution data in ', trim(trim(fname)//'_itxxxxx.plt')
  else
    write(*,*)
    write(*,'(a,a)') ' writing grid and solution data in ', trim(trim(fname)//'_itxxxxx.plt')
  endif


  !--
  do it=mt1, mt2, mt3
    write(fnameit5, '(i5.5)' ) it
    fnameplt = trim(fname)//'_it'//fnameit5//'.plt'

    sol_time = dble(it)

    do ip=1,np
      call readd(itaped, var1, it, ip, ierr)
      var(1:nnodes,ip) = var1(1:nnodes)
    enddo

    if (lGridSolTogether) then
      call tecplot_write_grid_solution (nnodes, ncells, np, grid_dim, imach, trim(fnameplt), trim(varinfo), cell2node, node_xy, var, sol_time)
    else
      call tecplot_write_solution (nnodes, ncells, np, imach, trim(fnameplt), trim(varinfo), var, sol_time)
    endif
  enddo


end program ios2tecplot


subroutine StripSpaces(string)
   character(len=*) :: string
   integer :: stringLen
   integer :: last, actual

   stringLen = len (string)
   last = 1
   actual = 1

   do while (actual < stringLen)
     if (string(last:last) == ' ') then
       actual = actual + 1
       string(last:last) = string(actual:actual)
       string(actual:actual) = ' '
     else
       last = last + 1
       if (actual < last) actual = last
     endif
   end do

   return
end subroutine
