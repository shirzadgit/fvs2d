! ifort -r8 -O3 -traceback read_gambit.f90 -o read_gambit.x
program gambit

  implicit none

  integer             :: nnodes, ncells_tri, ncells_quad, ncells, i,j, idum1,idum2,idum3,idum33, cell(1000000)
  integer             :: nbcs, ib
  integer,allocatable :: cell2node(:,:), tmp(:)
  real,allocatable    :: x(:),y(:)
  character           :: file_in*124, cdum*11, char39*39

  type bndry_type
    integer                       :: ncells
    integer,dimension(:),pointer  :: cell
    character(len=39)             :: type
  end type bndry_type
  type(bndry_type),dimension(:),pointer :: bndry



  ncells_tri = 0
  ncells_quad= 0

  write(*,*)
  write(*,'(a)',advance='no') 'gambit file name: '
  read (*,*) file_in


  open(100,file=trim(file_in))

  !-- skip 6 lines
  do i=1,6; read(100,*); enddo

  !-- read #s nodes, cells, and boundaries
  read(100,*) nnodes, ncells, idum1, nbcs

  !-- skip 4+nnodes lines
  do i=1,4+nnodes; read(100,*); enddo

  !-- read cell numbers and #triangles and #quads
  do i=1,ncells
    read(100,*) idum1,idum2,idum3
    if (idum3==3) then
      ncells_tri = ncells_tri + 1
      cell(i)=3
    endif
    if (idum3==4) then
      ncells_quad= ncells_quad+ 1
      cell(i)=4
    endif
  enddo
  close (100)

  if (ncells /=(ncells_tri + ncells_quad)) then
    write(*,*) 'error: # total cells /= ncells_tri + ncells_quad !!!'
    stop
  endif


  write(*,*)
  write(*,'(a,i0)') 'nnodes: ',nnodes
  write(*,'(a,i0)') 'ncells triangles: ',ncells_tri
  write(*,'(a,i0)') 'ncells quad: ',ncells_quad
  write(*,'(a,i0)') '# boundries: ',nbcs

  !-- allocate
  allocate(x(nnodes),y(nnodes), cell2node(ncells,4), bndry(nbcs))

  cell2node = 0


  open(100,file=trim(file_in))

  !-- skip 9 lines
  do i=1,9; read(100,*); enddo

  !-- read node coordinates
  do i=1,nnodes
    read(100,*) idum1,x(i),y(i)
  enddo

  !-- skip 2 lines
  do i=1,2; read(100,*); enddo

  !-- read node connectivities
  do i=1,ncells
    if (cell(i)==3) read(100,*) idum1,idum2,idum3,cell2node(i,1),cell2node(i,2),cell2node(i,3)
    if (cell(i)==4) read(100,*) idum1,idum2,idum3,cell2node(i,1),cell2node(i,2),cell2node(i,3),cell2node(i,4)
  enddo

  !-- skip 2 lines
  do i=1,5; read(100,*); enddo

  !-- skip int(ncells/10) lines
  do i=1,int(ncells/10); read(100,*); enddo
  if (ncells>10*int(ncells/10)) read(100,*);

  do ib=1,nbcs
    do i=1,2; read(100,*); enddo

    read(100,'(a39,i1,4x,i4)') char39,idum1,bndry(ib)%ncells
    bndry(ib)%type=trim(adjustl(char39))
    allocate(bndry(ib)%cell( bndry(ib)%ncells ))
    do i=1,bndry(ib)%ncells
      read(100,*) bndry(ib)%cell(i)
    enddo

    write(*,'(a,i0,a,a )') 'boundary(',ib,')%type   = ',(bndry(ib)%type)
    write(*,'(a,i0,a,i0)') 'boundary(',ib,')%ncells = ',bndry(ib)%ncells
  enddo
  close(100)

  !-- write grid file
  do i=1,len(trim(file_in))
    if (file_in(i:i)=='.') j=i-1
  enddo
  file_in = trim(file_in(1:j))

  write(*,*)
  !write(*,'(a)',advance='no') 'output file name (".grid" will be added): '
  write(*,*) 'writing grid file: ',trim(file_in)//".grid"
  !read (*,*) file_in
  open(100,file=trim(file_in)//".grid")
  write(100,*) '#nodes, #cells_tri, #cells_quad'
  write(100,*) nnodes, ncells_tri, ncells_quad
  do i=1,nnodes
    write(100,*) x(i),y(i)
  enddo

  do i=1,ncells
    if (cell(i)==3) write(100,*) cell2node(i,1),cell2node(i,2),cell2node(i,3)
  enddo
  do i=1,ncells
    if (cell(i)==4) write(100,*) cell2node(i,1),cell2node(i,2),cell2node(i,3),cell2node(i,4)
  enddo
  close(100)

  !-- write bc file
  allocate(tmp(ncells))
  j=0
  do i=1,ncells
    if (cell(i)==3) then
      j=j+1;
      tmp(i) = j
    endif
  enddo
  do i=1,ncells
    if (cell(i)==4) then
      j=j+1;
      tmp(i) = j
    endif
  enddo
  if (j/=ncells) write(*,*) 'no way !!!'

  write(*,*) 'writing bc file: ',trim(file_in)//".bc"
  open(100,file=trim(file_in)//".bc")
  write(100,'(i0,18x,a)') nbcs,'!--- number of boundaries'
  do ib=1,nbcs
    write(100,'(i0,3x,a,4x,a,i0,a)') bndry(ib)%ncells,trim(bndry(ib)%type), '!--- boundary ',ib,' : #ncells, type'
  enddo
  do ib=1,nbcs
    do j=1,bndry(ib)%ncells
      i=bndry(ib)%cell(j)
      write(100,'(i0)') tmp(i) !bndry(ib)%cell(i)
    enddo
  enddo
  close(100)


  !-- tecplot file
  write(*,*) 'writing tecplot file: ',trim(file_in)//".plt"
  open(100,file=trim(file_in)//'.plt')
  write(100,*) 'TITLE = "Example: 2D Finite Element Data"'
  write(100,*) 'VARIABLES = "X", "Y"'
  if (ncells_quad>0) then
    write(100,'(a,i0,a,i0,a)') 'ZONE NODES=',nnodes, ' ELEMENTS=',ncells, ' DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL'
  else
    write(100,'(a,i0,a,i0,a)') 'ZONE NODES=',nnodes, ' ELEMENTS=',ncells, ' DATAPACKING=POINT, ZONETYPE=FETRIANGLE'
  endif

  do i=1,nnodes
    write(100,*) x(i),y(i)
  enddo

  if (ncells_quad>0) then
    do i=1,ncells
      if (cell2node(i,4)==0) write(100,*) cell2node(i,1),cell2node(i,2),cell2node(i,3),cell2node(i,3)
    enddo
    do i=1,ncells
      if (cell2node(i,4)/=0) write(100,*) cell2node(i,1),cell2node(i,2),cell2node(i,3),cell2node(i,4)
    enddo
  else
    do i=1,ncells
      write(100,*) cell2node(i,1),cell2node(i,2),cell2node(i,3)
    enddo
  endif
  close(100)

end
