!------------------------------------------------------------------------------!
! module for ios input and output (use instead of IOS library).                !
!
! subroutine: 1) mkfname
!             2) writecd
!             3) writed
!             4) write8
!             5) write4
!             6) readcd
!             7) readd
!             8) read8
!             9) read4
!            10) writecd2d
!            11) writed2d
!            12) write82d
!            13) write42d
!
!------------------------------------------------------------------------------!
module ios_unstrc


implicit none

!--- parameters
   integer                              :: recl = 1 !SGI 4 ,otherwise 1
   integer,parameter                    :: next = 2
   character*3,dimension(2)             :: fext
   data fext/'.s4','.s8'/

!--- variables: integer
   integer,private                      :: i,j
   integer                              :: idot
   integer                              :: mrec

!--- arrays: integer
   integer,dimension(99),private,save     :: iunita
   integer*8,dimension(99),private,save   :: mpa
   integer*8,dimension(99),private,save   :: m1a
   integer*8,dimension(99),private,save   :: m2a
   integer*8,dimension(99),private,save   :: m3a
   integer*8,dimension(99),private,save   :: maa
   integer*8,dimension(99),private,save   :: mba
   integer*8,dimension(99),private,save   :: nta
   integer*8,dimension(99),private,save   :: npa
   integer*8,dimension(99),private,save   :: nza
   integer*8,dimension(99),private,save   :: nya
   integer*8,dimension(99),private,save   :: mreca
   integer*8,dimension(99),private,save   :: minfa
   integer*8,dimension(99),private,save   :: mta
   integer,dimension(99),private,save     :: imacha

   integer,parameter     :: maxtime=4*1024
   integer,parameter     :: maxlen=72
   integer,parameter     :: maxinf=128

contains

!------------------------------------------------------------------------------!
! subroutine: determines file ending.
!
!     purpose
!           determines the filetype corresponding to the extension of
!           filen and returns filename base without extension
!     input variables:
!           fname   :  filename
!           ldef    :  default wordlength in bytes (8 or 4)
!     output variables
!           fbase   :  basic filename (without extension)
!           imach   :  value of imach corresponding to extension
!
!------------------------------------------------------------------------------!
subroutine mkfname(fname, fbase,imach, ldef)

!--- variables
   character(LEN=124),intent(out)        :: fbase
   character(LEN=3)                      :: ex
   character(LEN=127)                    :: fname
   integer                               :: imach,ldef


   do i = 1,124
      fbase(i:i) =' '
   end do
   ex = 'xxx'
   imach = 0

   do i = 1,127
      if (fname(i:i) == '.') then
	  fbase(1:i-1)        = fname(1:i-1)
          do j = i,124
	  fbase(j:j)        = ' '
	  end do
	  ex(1:3)             = fname(i:i+2)
          idot = i
	  exit
      else if (fname(i:i) == ' ') then
          stop 'MKFNAME not an .s4 or .s8 file'
      end if
   end do


   do i = 1,2
      if (ex == fext(i)) then
	 imach = i
      end if
   end do

   if (imach == 0) then
      stop 'MKFNAME extension not allowed'
   end if

end subroutine

!------------------------------------------------------------------------------!
! subroutine: writes .cd file
!
!
!     purpose
!           prepares an output channel. It writes the
!           control data set to the '.cd' file and opens the actual
!           data file according to the user specified format.
!     INPUT parameters:
!     itape    - ordering number to switch between units (1<=itape<=30)
!     iunit    - unit number related to itape (1<=iunit<=999)
!     fbase    - basic filename (without extension for output)
!                 maximum of 11 letters
!     imach    - filetype specification:
!                1 = REAL*4
!                2 = REAL*8
!     mt       - number of timesteps
!     m3,m2,m1 - dimensions of array to be written (specify m3=1 if 2-d)
!                array (m1,m2,m3)
!     mp       - number of parameters to be written
!     itimes   - array (integer(mt)) with timestep numbers (mt<=512)
!     inf      - character*72 array containning information about
!                 parameters in the first mp elements and general
!                 information about the file in the following minf
!                 elements.
!     minf     - number of lines for suplementary information
!------------------------------------------------------------------------------!
subroutine writecd(itape,iunit,fbase,imach,mt,mnodes,mcells,mp,itimes,inf,minf,ierror)

!--- variables
   character(LEN=124),intent(in)         :: fbase
   character(LEN=127)                    :: fname
   integer                               :: itape
   integer                               :: iunit
   integer                               :: imach
   integer                               :: m3,m2,m1, mcells, mnodes
   integer                               :: mt,mp
   integer                               :: minf
   integer                               :: ierror
   integer,dimension(maxtime)            :: itimes
   character(LEN=maxlen),dimension(minf+mp)  :: inf

   m1=mnodes
   m2=1
   m3=1

!--- check input
!--- itape
      if (itape .lt.1) then
        write (*,*) '   ERROR: itape ',itape,' .lt.1 not allowed'
        stop 'WRITECD'
      else if (itape .gt.99) then
        write (*,*) '   ERROR: itape ',itape,' .gt.99 not allowed'
        stop 'WRITECD'
      end if
!---iunit
      if (iunit .lt.1) then
        write (*,*) '   ERROR: iunit ',iunit,' .lt.1 not allowed'
        stop 'WRITECD'
      end if
!---mt
      if (mt .lt.1) then
        write (*,*) '   ERROR: mt ',mt,' .lt.1 not allowed'
        stop 'WRITECD'
      end if
!---m1,m2,m3
      if (m1 .lt.1) then
        write (*,*) '   ERROR: m1 ',m1,' .lt.1 not allowed'
        stop 'WRITECD'
      end if
      if (m2 .lt.1) then
        write (*,*) '   ERROR: m2 ',m2,' .lt.1 not allowed'
        stop 'WRITECD'
      end if
      if (m3 .lt.1) then
        write (*,*) '   ERROR: m3 ',m3,' .lt.1 not allowed'
        stop 'WRITECD'
      end if
!---mp
      if (mp .lt.1) then
        write (*,*) '   ERROR: mp ',mp,' .lt.1 not allowed'
        stop 'WRITECD'
      end if
!---minf
      if (minf .lt.0) then
        write (*,*) '   ERROR: minf ',minf,' .lt.0 not allowed'
        stop 'WRITECD'
      end if

!---store data in common
      iunita(itape)     =     iunit
      m1a   (itape)     =     m1
      m2a   (itape)     =     m2
      m3a   (itape)     =     m3
      maa   (itape)     =     m1 * m2 * m3
      mpa   (itape)     =     mp
      mta   (itape)     =     mt
      minfa(itape)      =     minf
      nta(itape)    =     1
      npa(itape)    =     0
      imacha(itape) = imach

!---rebuild fname
      do i = 1,124
      if (fbase(i:i) == ' ') then
          idot = i
          exit
      end if
      end do
      fname(1:idot+2) = fbase(1:idot-1)//'.cd'
      do i = idot+3,127
         fname(i:i) = ' '
      end do

!--- write .cd file
      open(iunit,file=fname,iostat=ierror)
      rewind(iunit)

!-- write data
!       write(iunit,1010)m3,m2,m1,mp,mt,minf
! 1010  format(5x,'size of array:'/5x,'m3 = ',i5,5x,'m2 = ',i5,5x,  &
!       'm1 = ',i5/5x,'number of parameters = ',i5/5x,		  &
!       'number of timesteps  = ',i5				  &
!       //5x,'Information about file :',3x,'(',i3,2x,'info lines )')

      write(iunit,1010) mnodes, mcells, mp, mt, minf
1010  format(5x,'number of nodes = ',i5/5x, &
                'number of cells = ',i5/5x, &
                'number of parameters = ',i5/5x,	&
                'number of timesteps  = ',i5//5x, &
                'Information about file :',3x,'(',i3,2x,'info lines )')

      write(iunit,1020)(inf(i),i=mp+1,mp+minf)
1020  format(3x,a72)
      write(iunit,*)'     Information about parameters :'
      write(iunit,1030)(inf(i),i=1,mp)
1030  format(3x,a72)
      write(iunit,1040)
1040  format(2x,'Numbers of timesteps :')
      write(iunit,1050)(itimes(i),i=1,min(mt,5000))
1050  format(6(2x,i10))
      close(iunit)

!     calculate recordlenght and open file
!     mreca is recordlenght in numbers of specified type
!     mrec is actual recordlenght in units according to machine

!---rebuild fname
      do i = 1,124
      if (fbase(i:i) == ' ') then
          idot = i
          exit
      end if
      end do
      fname(1:idot+2) = fbase(1:idot-1)//fext(imach)


!--- choose file type
            if (imacha(itape) .eq. 1) then
!             sun: real = 4 bytes, 32 kbytes max. recordlength
!             sgi: real = 1 word, 8 kwords max. recordlength
!              mreca(itape) = min(8000,maa(itape))
              mreca(itape) = maa(itape)
              mrec  = mreca(itape) *recl
            else if (imacha(itape) .eq. 2) then
!             sun: real = 8 bytes, 32 kbytes max. recordlength
!             sgi: real = 2 words, 8 kwords max. recordlength
!              mreca(itape) = min(4000,maa(itape))
              mreca(itape) = maa(itape)
              mrec  = mreca(itape) * 2 *recl
            end if
!--- open .s4 or .s8 file
            open(unit=iunit,recl=mrec,file=fname,    &
     		access='DIRECT',form='unformatted',  &
     		iostat=ierror)

end subroutine


!------------------------------------------------------------------------------!
! subroutine: writes .s4 or .s8 file
!
!     is called when writing actual data.
!     It determines the chosen data format and calls the appropriate
!     writing routine
!------------------------------------------------------------------------------!
      subroutine writed (itape,data)

!--- variables
      integer :: itape
      real*8,dimension(maa(itape))        :: data

!---error check
      if (nta(itape).gt.mta(itape)) then
            write (*,1000) iunita(itape), nta(itape), mta(itape)
            stop 'WRITED'
      end if
 1000 format (3x,'ERROR IN SR WRITED, writing on unit ',i2/	 &
     	    3x,i10,' timesteps are more then allowed'/		 &
     	    3x,i10,' timesteps have been declared when opening')


!---increase counter variables for time and parameter
      npa(itape)  =     npa(itape) + 1
      if (npa(itape).gt.mpa(itape)) then
            npa(itape)  =     1
            nta(itape)  =     nta(itape)+1
      end if

   if (imacha(itape) == 1) then
      call write4(itape,data)
   else if (imacha(itape) ==2) then
      call write8(itape,data)
   else
      stop 'WRITED: imach wrong'
   end if

!      call flush (iunita(itape))

end subroutine

!------------------------------------------------------------------------------!
! subroutine: writes .s8 file
!
!     writes real*8 data to unformatted direct access files.
!
!------------------------------------------------------------------------------!

subroutine write8 (itape,data)

      integer :: itape
      real*8  :: data(maa(itape))
      integer :: adab2,adrec1,ierror
      integer*8 :: adab1


!---calculate adresses
      adab1       =     (nta(itape) - 1) * mpa(itape) * maa(itape)&
                 +     (npa(itape) - 1) * maa(itape) + 1
      adab2 =     adab1 + maa(itape) - 1
      adrec1=     int((adab1-1)/maa(itape)) + 1
!      adpos1=     mod((adab1-1),mreca(itape)) + 1
!      adrec2=     int((adab2-1)/maa(itape)) + 1
!     adpos2=     mod((adab2-1),mreca(itape)) + 1

!---loop over records
    write(iunita(itape),rec=adrec1,err=3000,iostat=ierror) &
       data

      return
!--- error messages
3000  write (*,*) ' SR WRITE8: ERROR WRITING ON UNIT ',iunita(itape)
      write(*,*)'   IOSTAT = ',ierror
      stop 'WRITE8'
end subroutine

!------------------------------------------------------------------------------!
! subroutine: writes .s4 file
!
!     writes real*4 data to unformatted direct access files.
!
!------------------------------------------------------------------------------!

subroutine write4 (itape,data)

      integer :: itape
      real*8 data(maa(itape))
      real*4 data4(maa(itape))
      integer :: adab2,adrec1,ierror
      integer*8 :: adab1


!---calculate adresses
      adab1       =     (nta(itape) - 1) * mpa(itape) * maa(itape)&
                 +     (npa(itape) - 1) * maa(itape) + 1
      adab2 =     adab1 + maa(itape) - 1
      adrec1=     int((adab1-1)/maa(itape)) + 1
!      adpos1=     mod((adab1-1),mreca(itape)) + 1
!      adrec2=     int((adab2-1)/maa(itape)) + 1
!     adpos2=     mod((adab2-1),mreca(itape)) + 1

!---loop over records
    data4=data
    write(iunita(itape),rec=adrec1,err=3000,iostat=ierror) data4

      return
!--- error messages
3000  write (*,*) ' SR WRITE4: ERROR WRITING ON UNIT ',iunita(itape)
      write(*,*)'   IOSTAT = ',ierror
      stop 'WRITE4'
end subroutine



!------------------------------------------------------------------------------!
! subroutine: reads .cd file
!
!     purpose:
!     reads the control data for the specified file
!     and opens the actual data file according to the user specifications
!     INPUT parameters:
!     itape    - ordering number to switch between units (1<=itape<=30)
!     iunit    - unit number related to itape (1<=iunit<=999)
!     fbase    - basic filename (without extension)
!                 maximum of 11 letters
!     imach    - filetype specification:
!                1 = REAL*4
!                2 = REAL*8
!     OUTPUT parameters
!     mt       - number of timesteps
!     m3,m2,m1 - dimensions of array to be read (specify m3=1 if 2-d)
!                array (m1,m2,m3)
!     mp       - number of parameters to be read
!     itimes   - array (integer(mt)) with timestep numbers (mt<=512)
!     inf      - character*72 array containning information about
!                 parameters in the first mp elements and general
!                 information about the file in the following minf
!                 elements.
!     minf     - number of lines for suplementary information
!
!------------------------------------------------------------------------------!


!subroutine readcd(itape,iunit,fbase,imach,mt,m3,m2,m1,mp,itimes,inf,minf,ierror)
subroutine readcd(itape,iunit,fbase,imach,mt,mnodes,mcells,mp,itimes,inf,minf,ierror)

!---variables
   character(LEN=124),intent(in)         :: fbase
   character(LEN=127)                    :: fname
   integer                               :: itape
   integer                               :: iunit
   integer                               :: imach
   integer                               :: m3,m2,m1, mcells, mnodes
   integer                               :: mt,mp
   integer                               :: minf
   integer                               :: ierror
   integer,dimension(maxtime)               :: itimes
   character(LEN=maxlen),dimension(maxinf)      :: inf



!---check input
!---itape
      if (itape .lt.1) then
        write (*,*) '   ERROR: itape ',itape,' .lt.1 not allowed'
        stop 'READCD'
      else if (itape .gt.99) then
        write (*,*) '   ERROR: itape ',itape,' .gt.99 not allowed'
        stop 'READCD'
      end if
!---iunit
      if (iunit .lt.1) then
        write (*,*) '   ERROR: iunit ',iunit,' .lt.1 not allowed'
        stop 'READCD'
      end if

      imacha(itape) = imach

!---rebuild fname
      do i = 1,124
      if (fbase(i:i) == ' ') then
          idot = i
          exit
      end if
      end do
      fname(1:idot+2) = fbase(1:idot-1)//'.cd'
      do i = idot+3,127
         fname(i:i) = ' '
      end do

!---read control data
      !open(iunit,file=trim(fname),err=2000,iostat=ierror,status='old')
      open(iunit,file=trim(trim(fbase)//'.cd'));
      rewind(iunit)

!       read(iunit,1010) m3,m2,m1,mp,mt,minf
! 1010  format(/10x,i5,10x,i5,10x,i5/28x,i5/28x,i5//33x,i3)

      read(iunit,1010) mnodes, mcells, mp, mt, minf
1010  format(23x,i, /23x,i, /28x,i5, /28x,i5, //33x,i3);



      read(iunit,1020)(inf(i),i=mp+1,mp+minf)
1020  format(3x,a72)
      read(iunit,'()')
      read(iunit,1030)(inf(i),i=1,mp)
1030  format(3x,a72)
      read(iunit,'()')
      read(iunit,*)(itimes(i),i=1,min(mt,5000))
      close(iunit)

      m1 = mnodes
      m2 = 1
      m3 = 1

!---store data in common
      iunita(itape)     =     iunit
      maa   (itape)     =     m1 * m2 * m3
      m1a   (itape)     =     m1
      m2a   (itape)     =     m2
      m3a   (itape)     =     m3
      mba   (itape)     =     m1 * m2
      mpa   (itape)     =     mp
      mta   (itape)     =     mt


!     calculate recordlenght and open file
!     mreca is recordlenght in numbers of specified type
!     mrec is actual recordlenght in units according to machine

!---rebuild fname
      do i = 1,124
      if (fbase(i:i) == ' ') then
          idot = i
          exit
      end if
      end do
      fname(1:idot+2) = fbase(1:idot-1)//fext(imach)
      do i = idot+3,127
         fname(i:i) = ' '
      end do


!---unformatted direct access filetype
            if (imacha(itape) .eq. 1) then
!             sun: real = 4 bytes, 32 kbytes max. recordlength
!             sgi: real = 1 word, 8 kwords max. recordlength
              mreca(itape) = maa(itape)
!SGI              mrec  = mreca(itape) * 4
              mrec  = mreca(itape) *recl
            else if (imacha(itape) .eq. 2) then
!             sun: real = 8 bytes, 32 kbytes max. recordlength
!             sgi: real = 2 words, 8 kwords max. recordlength
              mreca(itape) = maa(itape)
!SGI              mrec  = mreca(itape) * 8
              mrec  = mreca(itape) * 2 *recl
            end if

!--- open file
	    open(unit=iunit,err=2000,file=fname,status='old',  &
!                access='DIRECT',form='unformatted/ieee',     &
    	      access='DIRECT',form='unformatted',	       &
    	      iostat=ierror,recl=mrec)
      return

!---error procedure and end
2000  write(*,*)'   ERROR OPENING UNIT ',iunit,', ',fname
      write(*,*)'   IOSTAT = ',ierror
      write(*,*)'   possible cause: file does not exist'
      stop 'READCD'
end subroutine

!------------------------------------------------------------------------------!
! subroutine: reads .s4 and .s8 file
!
!     is called when reading actual data.
!     It determines the chosen data format and calls the appropriate
!     reading routine
!
!------------------------------------------------------------------------------!
      subroutine readd (itape,data,nt,np,ierror)

!---variables
      integer                      :: nt,np,ierror
      real*8,dimension(maa(itape)) :: data
      integer :: itape

!---error check
      if (nt .le. 0) then
        write (*,1000) iunita(itape),nt
        stop 'READD'
      else if (nt .gt. mta(itape) ) then
        write (*,1010) iunita(itape), nt, mta(itape)
        stop 'READD'
      end if
      if (np .le. 0) then
        write (*,1020) iunita(itape),np
        stop 'READD'
      else if (np .gt. mpa(itape) ) then
        write (*,1030) iunita(itape), np, mpa(itape)
        stop 'READD'
      end if
 1000 format (3x,'ERROR: reading from unit ',i2/         &
            3x,'timestep-index',i10,' is less then 1')
 1010 format (3x,'ERROR: reading from unit ',i2/ &
              3x,'timestep-index',i10,' is greater then ',i10)
 1020 format (3x,'ERROR: reading from unit ',i2/ &
            3x,'array number',i3,' is less then 1')
 1030 format (3x,'ERROR: reading from unit ',i2/ &
           3x,'array number',i3,' is greater then ',i3)


   if (imacha(itape) == 1) then
     call read4(itape,data,nt,np,ierror)
   else if (imacha(itape) == 2) then
     call read8(itape,data,nt,np,ierror)
   else
     stop 'READD'
   end if
end subroutine

!------------------------------------------------------------------------------!
! subroutine: reads .s8 file
!
!     reads real*8 data from unformatted direct access files.
!
!------------------------------------------------------------------------------!
subroutine read8 (itape,data,nt,np,ierror)

!---variables
      integer                      :: nt,np,ierror
      real*8,dimension(maa(itape)) :: data
      integer :: adab2,adrec1
      integer*8 :: adab1
      integer :: itape

!---calculate adresses
      adab1       =     (nt - 1) * mpa(itape) * maa(itape)&
                  +     (np - 1) * maa(itape) + 1
      adab2 =     adab1 + maa(itape) - 1
      adrec1=     int((adab1-1)/mreca(itape)) + 1
!      adpos1=     mod((adab1-1),mreca(itape)) + 1
!      adrec2=     int((adab2-1)/mreca(itape)) + 1
!     adpos2=     mod((adab2-1),mreca(itape)) + 1
!      j1    =     1
!      j2    =     mreca(itape) - adpos1 + 1

	    read(iunita(itape),rec=adrec1) data

      return
2000  write (*,*) ' ERROR READING FROM UNIT ',iunita(itape)
      write(*,*)'   IOSTAT = ',ierror
      stop 'READ8'
end subroutine

!------------------------------------------------------------------------------!
! subroutine: reads .s4 file
!
!     reads real*8 data from unformatted direct access files.
!
!------------------------------------------------------------------------------!
subroutine read4 (itape,data,nt,np,ierror)

!---variables
      integer                      :: nt,np,ierror
      real*8,dimension(maa(itape)) :: data
      real*4,dimension(maa(itape)) :: data4
      integer :: adab2,adrec1
      integer*8 :: adab1
      integer :: itape

!---calculate adresses
      adab1       =     (nt - 1) * mpa(itape) * maa(itape)&
                  +     (np - 1) * maa(itape) + 1
      adab2 =     adab1 + maa(itape) - 1
      adrec1=     int((adab1-1)/mreca(itape)) + 1
!      adpos1=     mod((adab1-1),mreca(itape)) + 1
!      adrec2=     int((adab2-1)/mreca(itape)) + 1
!     adpos2=     mod((adab2-1),mreca(itape)) + 1
!      j1    =     1
!      j2    =     mreca(itape) - adpos1 + 1


      read(iunita(itape),rec=adrec1,err=2000,iostat=ierror) data4
      data=data4
      return
2000  write (*,*) ' ERROR READING FROM UNIT ',iunita(itape)
      write(*,*)'   IOSTAT = ',ierror
      stop 'READ8'
end subroutine

!------------------------------------------------------------------------------!
! subroutine: writes .cd file for output in 2D slices
!
!
!     purpose
!           prepares an output channel. It writes the
!           control data set to the '.cd' file and opens the actual
!           data file according to the user specified format.
!     INPUT parameters:
!     itape    - ordering number to switch between units (1<=itape<=30)
!     iunit    - unit number related to itape (1<=iunit<=999)
!     fbase    - basic filename (without extension for output)
!                 maximum of 11 letters
!     imach    - filetype specification:
!                1 = REAL*4
!                2 = REAL*8
!     mt       - number of timesteps
!     m3,m2,m1 - dimensions of array to be written (specify m3=1 if 2-d)
!                array (m1,m2,m3)
!     mp       - number of parameters to be written
!     itimes   - array (integer(mt)) with timestep numbers (mt<=512)
!     inf      - character*72 array containning information about
!                 parameters in the first mp elements and general
!                 information about the file in the following minf
!                 elements.
!     minf     - number of lines for suplementary information
!------------------------------------------------------------------------------!
subroutine writecd2d(itape,iunit,fbase,imach,mt,m3,m2,m1,mp,itimes,inf,minf&
                   ,ierror)

!--- variables
   character(LEN=124),intent(in)         :: fbase
   character(LEN=127)                    :: fname
   integer                               :: itape
   integer                               :: iunit
   integer                               :: imach
   integer                               :: m3,m2,m1
   integer                               :: mt,mp
   integer                               :: minf
   integer                               :: ierror
   integer,dimension(5000)               :: itimes
   character(LEN=72),dimension(minf+mp)  :: inf

!--- check input
!--- itape
      if (itape .lt.1) then
        write (*,*) '   ERROR: itape ',itape,' .lt.1 not allowed'
        stop 'WRITECD'
      else if (itape .gt.99) then
        write (*,*) '   ERROR: itape ',itape,' .gt.99 not allowed'
        stop 'WRITECD'
      end if
!---iunit
      if (iunit .lt.1) then
        write (*,*) '   ERROR: iunit ',iunit,' .lt.1 not allowed'
        stop 'WRITECD'
      end if
!---mt
      if (mt .lt.1) then
        write (*,*) '   ERROR: mt ',mt,' .lt.1 not allowed'
        stop 'WRITECD'
      end if
!---m1,m2,m3
      if (m1 .lt.1) then
        write (*,*) '   ERROR: m1 ',m1,' .lt.1 not allowed'
        stop 'WRITECD'
      end if
      if (m2 .lt.1) then
        write (*,*) '   ERROR: m2 ',m2,' .lt.1 not allowed'
        stop 'WRITECD'
      end if
      if (m3 .lt.1) then
        write (*,*) '   ERROR: m3 ',m3,' .lt.1 not allowed'
        stop 'WRITECD'
      end if
!---mp
      if (mp .lt.1) then
        write (*,*) '   ERROR: mp ',mp,' .lt.1 not allowed'
        stop 'WRITECD'
      end if
!---minf
      if (minf .lt.0) then
        write (*,*) '   ERROR: minf ',minf,' .lt.0 not allowed'
        stop 'WRITECD'
      end if

!---store data in common
      iunita(itape)     =     iunit
      m1a   (itape)     =     m1
      m2a   (itape)     =     m2
      m3a   (itape)     =     m3
      maa   (itape)     =     m1 * m2 * m3
      mpa   (itape)     =     mp
      mta   (itape)     =     mt
      minfa(itape)      =     minf
      nta(itape)    =     1
      npa(itape)    =     1
      nza(itape)    =     0
      imacha(itape) = imach

!---rebuild fname
      do i = 1,124
      if (fbase(i:i) == ' ') then
          idot = i
          exit
      end if
      end do
      fname(1:idot+2) = fbase(1:idot-1)//'.cd'
      do i = idot+3,127
         fname(i:i) = ' '
      end do

!--- write .cd file
      open(iunit,file=fname,iostat=ierror)
      rewind(iunit)

!-- write data
      write(iunit,1010)m3,m2,m1,mp,mt,minf
1010  format(5x,'size of array:'/5x,'m3 = ',i5,5x,'m2 = ',i5,5x,  &
      'm1 = ',i5/5x,'number of parameters = ',i5/5x,		  &
      'number of timesteps  = ',i5				  &
      //5x,'Information about file :',3x,'(',i3,2x,'info lines )')
      write(iunit,1020)(inf(i),i=mp+1,mp+minf)
1020  format(3x,a72)
      write(iunit,*)'     Information about parameters :'
      write(iunit,1030)(inf(i),i=1,mp)
1030  format(3x,a72)
      write(iunit,1040)
1040  format(2x,'Numbers of timesteps :')
      write(iunit,1050)(itimes(i),i=1,min(mt,5000))
1050  format(6(2x,i10))
      close(iunit)

!     calculate recordlenght and open file
!     mreca is recordlenght in numbers of specified type
!     mrec is actual recordlenght in units according to machine

!---rebuild fname
      do i = 1,124
      if (fbase(i:i) == ' ') then
          idot = i
          exit
      end if
      end do
      fname(1:idot+2) = fbase(1:idot-1)//fext(imach)
      do i = idot+3,127
         fname(i:i) = ' '
      end do


!--- choose file type
            if (imacha(itape) .eq. 1) then
!             sun: real = 4 bytes, 32 kbytes max. recordlength
!             sgi: real = 1 word, 8 kwords max. recordlength
!              mreca(itape) = min(8000,maa(itape))
              mreca(itape) = m1a(itape)*m2a(itape)
              mrec  = mreca(itape) *recl
            else if (imacha(itape) .eq. 2) then
!             sun: real = 8 bytes, 32 kbytes max. recordlength
!             sgi: real = 2 words, 8 kwords max. recordlength
!              mreca(itape) = min(4000,maa(itape))
              mreca(itape) = m1a(itape)*m2a(itape)
              mrec  = mreca(itape) * 2 *recl
            end if
!--- open .s4 or .s8 file
            open(unit=iunit,recl=mrec,file=fname,    &
     		access='DIRECT',form='unformatted',  &
     		iostat=ierror)

end subroutine

!------------------------------------------------------------------------------!
! subroutine: writes .s4 or .s8 file (ouptut in 2d slices)
!
!     is called when writing actual data.
!     It determines the chosen data format and calls the appropriate
!     writing routine
!------------------------------------------------------------------------------!
      subroutine writed2d (itape,data)

!--- variables
      integer :: itape
      real*8,dimension(m1a(itape)*m2a(itape))        :: data

!---error check
      if (nta(itape).gt.mta(itape)) then
            write (*,1000) iunita(itape), nta(itape), mta(itape)
            stop 'WRITED'
      end if
 1000 format (3x,'ERROR IN SR WRITED, writing on unit ',i2/	 &
     	    3x,i10,' timesteps are more then allowed'/		 &
     	    3x,i10,' timesteps have been declared when opening')


!---increase counter variables for time and parameter
      nza(itape)   =     nza(itape) + 1
      if (nza(itape).gt.m3a(itape)) then
         nza(itape)     =    1
         npa(itape)  =     npa(itape) + 1
         if (npa(itape).gt.mpa(itape)) then
            npa(itape)  =     1
            nta(itape)  =     nta(itape)+1
         end if
       end if

   if (imacha(itape) == 1) then
      call write42d(itape,data)
   else if (imacha(itape) ==2) then
      call write82d(itape,data)
   else
      stop 'WRITED2D: imach wrong'
   end if

!      call flush (iunita(itape))

end subroutine


!------------------------------------------------------------------------------!
! subroutine: writes .s8 file (output in 2d slices)
!
!     writes real*8 data to unformatted direct access files.
!
!------------------------------------------------------------------------------!

subroutine write82d (itape,data)

      integer :: itape
      real*8  :: data(m1a(itape)*m2a(itape))
      integer :: adab2,adrec1,ierror
      integer*8 :: adab1


!---calculate adresses
      adab1       =     (nta(itape) - 1) * mpa(itape) * maa(itape)&
                  +     (npa(itape) - 1) * maa(itape)             &
		  +     (nza(itape) - 1) * m1a(itape)*m2a(itape) + 1
      adab2 =     adab1 + m1a(itape)*m2a(itape) - 1
      adrec1=     int((adab1-1)/m1a(itape)/m2a(itape)) + 1
!      adpos1=     mod((adab1-1),mreca(itape)) + 1
!      adrec2=     int((adab2-1)/maa(itape)) + 1
!     adpos2=     mod((adab2-1),mreca(itape)) + 1

!---loop over records
    write(iunita(itape),rec=adrec1,err=3000,iostat=ierror) &
       data

      return
!--- error messages
3000  write (*,*) ' SR WRITE8: ERROR WRITING ON UNIT ',iunita(itape)
      write(*,*)'   IOSTAT = ',ierror
      stop 'WRITE82D'
end subroutine

!------------------------------------------------------------------------------!
! subroutine: writes .s4 file (output in 2d slices)
!
!     writes real*4 data to unformatted direct access files.
!
!------------------------------------------------------------------------------!

subroutine write42d (itape,data)

      integer :: itape
      real*8 data(m1a(itape)*m2a(itape))
      real*4 data4(m1a(itape)*m2a(itape))
      integer :: adab2,adrec1,ierror
      integer*8 :: adab1


!---calculate adresses
      adab1       =     (nta(itape) - 1) * mpa(itape) * maa(itape)&
                  +     (npa(itape) - 1) * maa(itape)             &
		  +     (nza(itape) - 1) * m1a(itape)*m2a(itape) + 1
      adab2 =     adab1 + m1a(itape)*m2a(itape) - 1
      adrec1=     int((adab1-1)/(m1a(itape)*m2a(itape))) + 1
!      adpos1=     mod((adab1-1),mreca(itape)) + 1
!      adrec2=     int((adab2-1)/maa(itape)) + 1
!     adpos2=     mod((adab2-1),mreca(itape)) + 1

!---loop over records
    data4=data
    write(iunita(itape),rec=adrec1,err=3000,iostat=ierror) data4

      return
!--- error messages
3000  write (*,*) ' SR WRITE4: ERROR WRITING ON UNIT ',iunita(itape)
      write(*,*)'   IOSTAT = ',ierror
      stop 'WRITE42D'
end subroutine

!------------------------------------------------------------------------------!
! subroutine: reads .cd file (input in 2d slices)
!
!     purpose:
!     reads the control data for the specified file
!     and opens the actual data file according to the user specifications
!     INPUT parameters:
!     itape    - ordering number to switch between units (1<=itape<=30)
!     iunit    - unit number related to itape (1<=iunit<=999)
!     fbase    - basic filename (without extension)
!                 maximum of 11 letters
!     imach    - filetype specification:
!                1 = REAL*4
!                2 = REAL*8
!     OUTPUT parameters
!     mt       - number of timesteps
!     m3,m2,m1 - dimensions of array to be read (specify m3=1 if 2-d)
!                array (m1,m2,m3)
!     mp       - number of parameters to be read
!     itimes   - array (integer(mt)) with timestep numbers (mt<=512)
!     inf      - character*72 array containning information about
!                 parameters in the first mp elements and general
!                 information about the file in the following minf
!                 elements.
!     minf     - number of lines for suplementary information
!
!------------------------------------------------------------------------------!


subroutine readcd2d(itape,iunit,fbase,imach,mt,m3,m2,m1,mp,itimes,inf,minf,&
                    ierror)


!---variables
   character(LEN=124),intent(in)         :: fbase
   character(LEN=127)                    :: fname
   integer                               :: itape
   integer                               :: iunit
   integer                               :: imach
   integer                               :: m3,m2,m1
   integer                               :: mt,mp
   integer                               :: minf
   integer                               :: ierror
   integer,dimension(5000)               :: itimes
   character(LEN=72),dimension(128)      :: inf



!---check input
!---itape
      if (itape .lt.1) then
        write (*,*) '   ERROR: itape ',itape,' .lt.1 not allowed'
        stop 'READCD'
      else if (itape .gt.99) then
        write (*,*) '   ERROR: itape ',itape,' .gt.99 not allowed'
        stop 'READCD'
      end if
!---iunit
      if (iunit .lt.1) then
        write (*,*) '   ERROR: iunit ',iunit,' .lt.1 not allowed'
        stop 'READCD'
      end if

      imacha(itape) = imach

!---rebuild fname
      do i = 1,124
      if (fbase(i:i) == ' ') then
          idot = i
          exit
      end if
      end do
      fname(1:idot+2) = fbase(1:idot-1)//'.cd'
      do i = idot+3,127
         fname(i:i) = ' '
      end do

!---read control data
      open(iunit,file=fname,err=2000,iostat=ierror,status='old')
      rewind(iunit)
      read(iunit,1010)m3,m2,m1,mp,mt,minf
1010  format(/10x,i5,10x,i5,10x,i5/28x,i5/28x,i5//33x,i3)
      read(iunit,1020)(inf(i),i=mp+1,mp+minf)
1020  format(3x,a72)
      read(iunit,'()')
      read(iunit,1030)(inf(i),i=1,mp)
1030  format(3x,a72)
      read(iunit,'()')
      read(iunit,*)(itimes(i),i=1,min(mt,5000))
      close(iunit)


!---store data in common
      iunita(itape)     =     iunit
      maa   (itape)     =     m1 * m2 * m3
      m1a   (itape)     =     m1
      m2a   (itape)     =     m2
      m3a   (itape)     =     m3
      mba   (itape)     =     m1 * m2
      mpa   (itape)     =     mp
      mta   (itape)     =     mt


!     calculate recordlenght and open file
!     mreca is recordlenght in numbers of specified type
!     mrec is actual recordlenght in units according to machine

!---rebuild fname
      do i = 1,124
      if (fbase(i:i) == ' ') then
          idot = i
          exit
      end if
      end do
      fname(1:idot+2) = fbase(1:idot-1)//fext(imach)
      do i = idot+3,127
         fname(i:i) = ' '
      end do


!---unformatted direct access filetype
            if (imacha(itape) .eq. 1) then
!             sun: real = 4 bytes, 32 kbytes max. recordlength
!             sgi: real = 1 word, 8 kwords max. recordlength
              mreca(itape) = m1a(itape)*m2a(itape)
!SGI              mrec  = mreca(itape) * 4
              mrec  = mreca(itape) *recl
            else if (imacha(itape) .eq. 2) then
!             sun: real = 8 bytes, 32 kbytes max. recordlength
!             sgi: real = 2 words, 8 kwords max. recordlength
              mreca(itape) = m1a(itape)*m2a(itape)
!SGI              mrec  = mreca(itape) * 8
              mrec  = mreca(itape) * 2 *recl
            end if

!--- open file
            open(unit=iunit,err=2000,file=fname,status='old',  &
!                access='DIRECT',form='unformatted/ieee',     &
    	      access='DIRECT',form='unformatted',	       &
    	      iostat=ierror,recl=mrec)
      return

!---error procedure and end
2000  write(*,*)'   ERROR OPENING UNIT ',iunit,', ',fname
      write(*,*)'   IOSTAT = ',ierror
      write(*,*)'   possible cause: file does not exist'
      stop 'READCD'
end subroutine

!------------------------------------------------------------------------------!
! subroutine: reads .s4 and .s8 file (input in 2d slices)
!
!     is called when reading actual data.
!     It determines the chosen data format and calls the appropriate
!     reading routine
!
!------------------------------------------------------------------------------!
      subroutine readd2d (itape,data,nt,np,nk,ierror)

!---variables
      integer                      :: nt,np,nk,ierror
      real*8,dimension(m1a(itape)*m2a(itape)) ::data
      integer :: itape

!---error check
      if (nt .le. 0) then
        write (*,1000) iunita(itape),nt
        stop 'READD'
      else if (nt .gt. mta(itape) ) then
        write (*,1010) iunita(itape), nt, mta(itape)
        stop 'READD'
      end if
      if (np .le. 0) then
        write (*,1020) iunita(itape),np
        stop 'READD'
      else if (np .gt. mpa(itape) ) then
        write (*,1030) iunita(itape), np, mpa(itape)
        stop 'READD'
      end if
      if (nk .le. 0) then
        write (*,1020) iunita(itape),nk
        stop 'READD'
      else if (nk .gt. m3a(itape) ) then
        write (*,1030) iunita(itape), nk, m3a(itape)
        stop 'READD'
      end if
 1000 format (3x,'ERROR: reading from unit ',i2/         &
            3x,'timestep-index',i10,' is less then 1')
 1010 format (3x,'ERROR: reading from unit ',i2/ &
              3x,'timestep-index',i10,' is greater then ',i10)
 1020 format (3x,'ERROR: reading from unit ',i2/ &
            3x,'array number',i3,' is less then 1')
 1030 format (3x,'ERROR: reading from unit ',i2/ &
           3x,'array number',i3,' is greater then ',i3)


   if (imacha(itape) == 1) then
     call read42d(itape,data,nt,np,nk,ierror)
   else if (imacha(itape) == 2) then
     call read82d(itape,data,nt,np,nk,ierror)
   else
     stop 'READD'
   end if
end subroutine

!------------------------------------------------------------------------------!
! subroutine: reads .s8 file (input in 2d slices)
!
!     reads real*8 data from unformatted direct access files.
!
!------------------------------------------------------------------------------!
subroutine read82d (itape,data,nt,np,nk,ierror)

!---variables
      integer                      :: nt,np,nk,ierror
      real*8,dimension(m1a(itape)*m2a(itape)) :: data
      integer :: adab2,adrec1
      integer*8 :: adab1
      integer :: itape

!---calculate adresses
      adab1       =     (nt - 1) * mpa(itape) * maa(itape) &
                  +     (np - 1) * maa(itape)              &
		  +     (nk - 1) * m1a(itape)*m2a(itape) + 1
      adab2 =     adab1 + m1a(itape)*m2a(itape) - 1
      adrec1=     int((adab1-1)/mreca(itape)) + 1
!      adpos1=     mod((adab1-1),mreca(itape)) + 1
!      adrec2=     int((adab2-1)/mreca(itape)) + 1
!     adpos2=     mod((adab2-1),mreca(itape)) + 1
!      j1    =     1
!      j2    =     mreca(itape) - adpos1 + 1

            read(iunita(itape),rec=adrec1,err=2000,iostat=ierror) data

      return
2000  write (*,*) ' ERROR READING FROM UNIT ',iunita(itape)
      write(*,*)'   IOSTAT = ',ierror
      stop 'READ8'
end subroutine

!------------------------------------------------------------------------------!
! subroutine: reads .s4 file (input in 2d slices)
!
!     reads real*8 data from unformatted direct access files.
!
!------------------------------------------------------------------------------!
subroutine read42d (itape,data,nt,np,nk,ierror)

!---variables
      integer                      :: nt,np,nk,ierror
      real*8,dimension(m1a(itape)*m2a(itape)) :: data
      real*4,dimension(m1a(itape)*m2a(itape)) :: data4
      integer :: adab2,adrec1
      integer*8 :: adab1
      integer :: itape

!---calculate adresses
      adab1       =     (nt - 1) * mpa(itape) * maa(itape)&
                  +     (np - 1) * maa(itape)               &
		  +     (nk - 1) * m1a(itape)*m2a(itape) + 1
      adab2 =     adab1 + m1a(itape)*m2a(itape) - 1
      adrec1=     int((adab1-1)/mreca(itape)) + 1
!      adpos1=     mod((adab1-1),mreca(itape)) + 1
!      adrec2=     int((adab2-1)/mreca(itape)) + 1
!     adpos2=     mod((adab2-1),mreca(itape)) + 1
!      j1    =     1
!      j2    =     mreca(itape) - adpos1 + 1


      read(iunita(itape),rec=adrec1,err=2000,iostat=ierror) data4
      data=data4
      return
2000  write (*,*) ' ERROR READING FROM UNIT ',iunita(itape)
      write(*,*)'   IOSTAT = ',ierror
      stop 'READ8'
end subroutine

!------------------------------------------------------------------------------!
! subroutine: writes .cd file for output in 2D slices
!
!
!     purpose
!           prepares an output channel. It writes the
!           control data set to the '.cd' file and opens the actual
!           data file according to the user specified format.
!     INPUT parameters:
!     itape    - ordering number to switch between units (1<=itape<=30)
!     iunit    - unit number related to itape (1<=iunit<=999)
!     fbase    - basic filename (without extension for output)
!                 maximum of 11 letters
!     imach    - filetype specification:
!                1 = REAL*4
!                2 = REAL*8
!     mt       - number of timesteps
!     m3,m2,m1 - dimensions of array to be written (specify m3=1 if 2-d)
!                array (m1,m2,m3)
!     mp       - number of parameters to be written
!     itimes   - array (integer(mt)) with timestep numbers (mt<=512)
!     inf      - character*72 array containning information about
!                 parameters in the first mp elements and general
!                 information about the file in the following minf
!                 elements.
!     minf     - number of lines for suplementary information
!------------------------------------------------------------------------------!
subroutine writecd1d(itape,iunit,fbase,imach,mt,m3,m2,m1,mp,itimes,inf,minf&
                   ,ierror)

!--- variables
   character(LEN=124),intent(in)         :: fbase
   character(LEN=127)                    :: fname
   integer                               :: itape
   integer                               :: iunit
   integer                               :: imach
   integer                               :: m3,m2,m1
   integer                               :: mt,mp
   integer                               :: minf
   integer                               :: ierror
   integer,dimension(5000)               :: itimes
   character(LEN=72),dimension(minf+mp)  :: inf

!--- check input
!--- itape
      if (itape .lt.1) then
        write (*,*) '   ERROR: itape ',itape,' .lt.1 not allowed'
        stop 'WRITECD'
      else if (itape .gt.99) then
        write (*,*) '   ERROR: itape ',itape,' .gt.99 not allowed'
        stop 'WRITECD'
      end if
!---iunit
      if (iunit .lt.1) then
        write (*,*) '   ERROR: iunit ',iunit,' .lt.1 not allowed'
        stop 'WRITECD'
      end if
!---mt
      if (mt .lt.1) then
        write (*,*) '   ERROR: mt ',mt,' .lt.1 not allowed'
        stop 'WRITECD'
      end if
!---m1,m2,m3
      if (m1 .lt.1) then
        write (*,*) '   ERROR: m1 ',m1,' .lt.1 not allowed'
        stop 'WRITECD'
      end if
      if (m2 .lt.1) then
        write (*,*) '   ERROR: m2 ',m2,' .lt.1 not allowed'
        stop 'WRITECD'
      end if
      if (m3 .lt.1) then
        write (*,*) '   ERROR: m3 ',m3,' .lt.1 not allowed'
        stop 'WRITECD'
      end if
!---mp
      if (mp .lt.1) then
        write (*,*) '   ERROR: mp ',mp,' .lt.1 not allowed'
        stop 'WRITECD'
      end if
!---minf
      if (minf .lt.0) then
        write (*,*) '   ERROR: minf ',minf,' .lt.0 not allowed'
        stop 'WRITECD'
      end if

!---store data in common
      iunita(itape)     =     iunit
      m1a   (itape)     =     m1
      m2a   (itape)     =     m2
      m3a   (itape)     =     m3
      maa   (itape)     =     m1 * m2 * m3
      mpa   (itape)     =     mp
      mta   (itape)     =     mt
      minfa(itape)      =     minf
      nta(itape)    =     1
      npa(itape)    =     1
      nza(itape)    =     1
      nya(itape)    =     0
      imacha(itape) = imach

!---rebuild fname
      do i = 1,124
      if (fbase(i:i) == ' ') then
          idot = i
          exit
      end if
      end do
      fname(1:idot+2) = fbase(1:idot-1)//'.cd'
      do i = idot+3,127
         fname(i:i) = ' '
      end do

!--- write .cd file
      open(iunit,file=fname,iostat=ierror)
      rewind(iunit)

!-- write data
      write(iunit,1010)m3,m2,m1,mp,mt,minf
1010  format(5x,'size of array:'/5x,'m3 = ',i5,5x,'m2 = ',i5,5x,  &
      'm1 = ',i5/5x,'number of parameters = ',i5/5x,		  &
      'number of timesteps  = ',i5				  &
      //5x,'Information about file :',3x,'(',i3,2x,'info lines )')
      write(iunit,1020)(inf(i),i=mp+1,mp+minf)
1020  format(3x,a72)
      write(iunit,*)'     Information about parameters :'
      write(iunit,1030)(inf(i),i=1,mp)
1030  format(3x,a72)
      write(iunit,1040)
1040  format(2x,'Numbers of timesteps :')
      write(iunit,1050)(itimes(i),i=1,min(mt,5000))
1050  format(6(2x,i10))
      close(iunit)

!     calculate recordlenght and open file
!     mreca is recordlenght in numbers of specified type
!     mrec is actual recordlenght in units according to machine

!---rebuild fname
      do i = 1,124
      if (fbase(i:i) == ' ') then
          idot = i
          exit
      end if
      end do
      fname(1:idot+2) = fbase(1:idot-1)//fext(imach)
      do i = idot+3,127
         fname(i:i) = ' '
      end do


!--- choose file type
            if (imacha(itape) .eq. 1) then
!             sun: real = 4 bytes, 32 kbytes max. recordlength
!             sgi: real = 1 word, 8 kwords max. recordlength
!              mreca(itape) = min(8000,maa(itape))
              mreca(itape) = m1a(itape)
              mrec  = mreca(itape) *recl
            else if (imacha(itape) .eq. 2) then
!             sun: real = 8 bytes, 32 kbytes max. recordlength
!             sgi: real = 2 words, 8 kwords max. recordlength
!              mreca(itape) = min(4000,maa(itape))
              mreca(itape) = m1a(itape)
              mrec  = mreca(itape) * 2 *recl
            end if
!--- open .s4 or .s8 file
            open(unit=iunit,recl=mrec,file=fname,    &
     		access='DIRECT',form='unformatted',  &
     		iostat=ierror)

end subroutine

!------------------------------------------------------------------------------!
! subroutine: writes .s4 or .s8 file (ouptut in 2d slices)
!
!     is called when writing actual data.
!     It determines the chosen data format and calls the appropriate
!     writing routine
!------------------------------------------------------------------------------!
      subroutine writed1d (itape,data)

!--- variables
      integer :: itape
      real*8,dimension(m1a(itape))        :: data

!---error check
      if (nta(itape).gt.mta(itape)) then
            write (*,1000) iunita(itape), nta(itape), mta(itape)
            stop 'WRITED'
      end if
 1000 format (3x,'ERROR IN SR WRITED, writing on unit ',i2/	 &
     	    3x,i10,' timesteps are more then allowed'/		 &
     	    3x,i10,' timesteps have been declared when opening')


!---increase counter variables for time and parameter
      nya(itape)   =     nya(itape) + 1
      if (nya(itape).gt.m2a(itape)) then
         nya(itape)     =    1
         nza(itape)     =    nza(itape) + 1
         if (nza(itape).gt.m3a(itape)) then
            nza(itape)     =    1
            npa(itape)  =     npa(itape) + 1
            if (npa(itape).gt.mpa(itape)) then
               npa(itape)  =     1
               nta(itape)  =     nta(itape)+1
            end if
          end if
       end if

   if (imacha(itape) == 1) then
      call write41d(itape,data)
   else if (imacha(itape) ==2) then
      call write81d(itape,data)
   else
      stop 'WRITED2D: imach wrong'
   end if

!      call flush (iunita(itape))

end subroutine


!------------------------------------------------------------------------------!
! subroutine: writes .s8 file (output in 2d slices)
!
!     writes real*8 data to unformatted direct access files.
!
!------------------------------------------------------------------------------!

subroutine write81d (itape,data)

      integer :: itape
      real*8  :: data(m1a(itape))
      integer :: adab2,adrec1,ierror
      integer*8 :: adab1


!---calculate adresses
      adab1       =     (nta(itape) - 1) * mpa(itape) * maa(itape)&
                  +     (npa(itape) - 1) * maa(itape)             &
		  +     (nza(itape) - 1) * m1a(itape)*m2a(itape)  &
		  +     (nya(itape) - 1) * m1a(itape)  + 1
      adab2 =     adab1 + m1a(itape) - 1
      adrec1=     int((adab1-1)/m1a(itape)) + 1
!      adpos1=     mod((adab1-1),mreca(itape)) + 1
!      adrec2=     int((adab2-1)/maa(itape)) + 1
!     adpos2=     mod((adab2-1),mreca(itape)) + 1

!---loop over records
    write(iunita(itape),rec=adrec1,err=3000,iostat=ierror) &
       data

      return
!--- error messages
3000  write (*,*) ' SR WRITE8: ERROR WRITING ON UNIT ',iunita(itape)
      write(*,*)'   IOSTAT = ',ierror
      stop 'WRITE82D'
end subroutine

!------------------------------------------------------------------------------!
! subroutine: writes .s4 file (output in 2d slices)
!
!     writes real*4 data to unformatted direct access files.
!
!------------------------------------------------------------------------------!

subroutine write41d (itape,data)

      integer :: itape
      real*8 data(m1a(itape))
      real*4 data4(m1a(itape))
      integer :: adab2,adrec1,ierror
      integer*8 :: adab1


!---calculate adresses
      adab1       =     (nta(itape) - 1) * mpa(itape) * maa(itape)&
                  +     (npa(itape) - 1) * maa(itape)             &
		  +     (nza(itape) - 1) * m1a(itape)*m2a(itape)  &
		  +     (nya(itape) - 1) * m1a(itape)  + 1
      adab2 =     adab1 + m1a(itape) - 1
      adrec1=     int((adab1-1)/(m1a(itape))) + 1
!      adpos1=     mod((adab1-1),mreca(itape)) + 1
!      adrec2=     int((adab2-1)/maa(itape)) + 1
!     adpos2=     mod((adab2-1),mreca(itape)) + 1

!---loop over records
    data4=data
    write(iunita(itape),rec=adrec1,err=3000,iostat=ierror) data4

      return
!--- error messages
3000  write (*,*) ' SR WRITE4: ERROR WRITING ON UNIT ',iunita(itape)
      write(*,*)'   IOSTAT = ',ierror
      stop 'WRITE42D'
end subroutine

!------------------------------------------------------------------------------!
! subroutine: writes .s4 or .s8 file (ouptut in 2d slices)
!
!     is called when writing actual data.
!     It determines the chosen data format and calls the appropriate
!     writing routine
!------------------------------------------------------------------------------!
      subroutine writed1dpos (itape,data,nt,np,nk,nj,ierror)

!--- variables
      integer                      :: nt,np,nk,nj,ierror
      integer :: itape
      real*8,dimension(m1a(itape))        :: data

!---error check
      if (nt .le. 0) then
        write (*,1000) iunita(itape),nt
        stop 'READD'
      else if (nt .gt. mta(itape) ) then
        write (*,1010) iunita(itape), nt, mta(itape)
        stop 'READD'
      end if
      if (np .le. 0) then
        write (*,1020) iunita(itape),np
        stop 'READD'
      else if (np .gt. mpa(itape) ) then
        write (*,1030) iunita(itape), np, mpa(itape)
        stop 'READD'
      end if
      if (nk .le. 0) then
        write (*,1020) iunita(itape),nk
        stop 'READD'
      else if (nk .gt. m3a(itape) ) then
        write (*,1030) iunita(itape), nk, m3a(itape)
        stop 'READD'
      end if
      if (nj .le. 0) then
        write (*,1020) iunita(itape),nj
        stop 'READD'
      else if (nj .gt. m2a(itape) ) then
        write (*,1030) iunita(itape), nj, m2a(itape)
        stop 'READD'
      end if
 1000 format (3x,'ERROR: reading from unit ',i2/         &
            3x,'timestep-index',i10,' is less then 1')
 1010 format (3x,'ERROR: reading from unit ',i2/ &
              3x,'timestep-index',i10,' is greater then ',i10)
 1020 format (3x,'ERROR: reading from unit ',i2/ &
            3x,'array number',i3,' is less then 1')
 1030 format (3x,'ERROR: reading from unit ',i2/ &
           3x,'array number',i3,' is greater then ',i3)


   if (imacha(itape) == 1) then
     call write41dpos(itape,data,nt,np,nk,nj,ierror)
   else if (imacha(itape) == 2) then
     call write81dpos(itape,data,nt,np,nk,nj,ierror)
   else
     stop 'WRITED1DPOS'
   end if


end subroutine


!------------------------------------------------------------------------------!
! subroutine: writes .s8 file (output in 2d slices)
!
!     writes real*8 data to unformatted direct access files.
!
!------------------------------------------------------------------------------!

subroutine write81dpos (itape,data,nt,np,nk,nj,ierror)

      integer                      :: nt,np,nk,nj
      integer :: itape
      real*8  :: data(m1a(itape))
      integer :: adab2,adrec1,ierror
      integer*8 :: adab1


!---calculate adresses
!---calculate adresses
      adab1       =     (nt - 1) * mpa(itape) * maa(itape) &
                  +     (np - 1) * maa(itape)              &
		  +     (nk - 1) * m1a(itape)*m2a(itape)   &
		  +     (nj - 1) * m1a(itape) + 1
      adab2 =     adab1 + m1a(itape) - 1
      adrec1=     int((adab1-1)/mreca(itape)) + 1

!---loop over records
    write(iunita(itape),rec=adrec1,err=3000,iostat=ierror) &
       data

      return
!--- error messages
3000  write (*,*) ' SR WRITE8: ERROR WRITING ON UNIT ',iunita(itape)
      write(*,*)'   IOSTAT = ',ierror
      stop 'WRITE81DPOS'
end subroutine

!------------------------------------------------------------------------------!
! subroutine: writes .s4 file (output in 2d slices)
!
!     writes real*4 data to unformatted direct access files.
!
!------------------------------------------------------------------------------!

subroutine write41dpos (itape,data,nt,np,nk,nj,ierror)

      integer                      :: nt,np,nk,nj
      integer :: itape
      real*8 data(m1a(itape))
      real*4 data4(m1a(itape))
      integer :: adab2,adrec1,ierror
      integer*8 :: adab1


!---calculate adresses
      adab1       =     (nt - 1) * mpa(itape) * maa(itape) &
                  +     (np - 1) * maa(itape)              &
		  +     (nk - 1) * m1a(itape)*m2a(itape)   &
		  +     (nj - 1) * m1a(itape) + 1
      adab2 =     adab1 + m1a(itape) - 1
      adrec1=     int((adab1-1)/mreca(itape)) + 1

!---loop over records
    data4=data
    write(iunita(itape),rec=adrec1,err=3000,iostat=ierror) data4

      return
!--- error messages
3000  write (*,*) ' SR WRITE4: ERROR WRITING ON UNIT ',iunita(itape)
      write(*,*)'   IOSTAT = ',ierror
      stop 'WRITE41DPOS'
end subroutine
!------------------------------------------------------------------------------!
! subroutine: reads .cd file (input in 2d slices)
!
!     purpose:
!     reads the control data for the specified file
!     and opens the actual data file according to the user specifications
!     INPUT parameters:
!     itape    - ordering number to switch between units (1<=itape<=30)
!     iunit    - unit number related to itape (1<=iunit<=999)
!     fbase    - basic filename (without extension)
!                 maximum of 11 letters
!     imach    - filetype specification:
!                1 = REAL*4
!                2 = REAL*8
!     OUTPUT parameters
!     mt       - number of timesteps
!     m3,m2,m1 - dimensions of array to be read (specify m3=1 if 2-d)
!                array (m1,m2,m3)
!     mp       - number of parameters to be read
!     itimes   - array (integer(mt)) with timestep numbers (mt<=512)
!     inf      - character*72 array containning information about
!                 parameters in the first mp elements and general
!                 information about the file in the following minf
!                 elements.
!     minf     - number of lines for suplementary information
!
!------------------------------------------------------------------------------!


subroutine readcd1d(itape,iunit,fbase,imach,mt,m3,m2,m1,mp,itimes,inf,minf,&
                    ierror)


!---variables
   character(LEN=124),intent(in)         :: fbase
   character(LEN=127)                    :: fname
   integer                               :: itape
   integer                               :: iunit
   integer                               :: imach
   integer                               :: m3,m2,m1
   integer                               :: mt,mp
   integer                               :: minf
   integer                               :: ierror
   integer,dimension(5000)               :: itimes
   character(LEN=72),dimension(128)      :: inf



!---check input
!---itape
      if (itape .lt.1) then
        write (*,*) '   ERROR: itape ',itape,' .lt.1 not allowed'
        stop 'READCD'
      else if (itape .gt.99) then
        write (*,*) '   ERROR: itape ',itape,' .gt.99 not allowed'
        stop 'READCD'
      end if
!---iunit
      if (iunit .lt.1) then
        write (*,*) '   ERROR: iunit ',iunit,' .lt.1 not allowed'
        stop 'READCD'
      end if

      imacha(itape) = imach

!---rebuild fname
      do i = 1,124
      if (fbase(i:i) == ' ') then
          idot = i
          exit
      end if
      end do
      fname(1:idot+2) = fbase(1:idot-1)//'.cd'
      do i = idot+3,127
         fname(i:i) = ' '
      end do

!---read control data
      open(iunit,file=fname,err=2000,iostat=ierror,status='old')
      rewind(iunit)
      read(iunit,1010)m3,m2,m1,mp,mt,minf
1010  format(/10x,i5,10x,i5,10x,i5/28x,i5/28x,i5//33x,i3)
      read(iunit,1020)(inf(i),i=mp+1,mp+minf)
1020  format(3x,a72)
      read(iunit,'()')
      read(iunit,1030)(inf(i),i=1,mp)
1030  format(3x,a72)
      read(iunit,'()')
      read(iunit,*)(itimes(i),i=1,min(mt,5000))
      close(iunit)


!---store data in common
      iunita(itape)     =     iunit
      maa   (itape)     =     m1 * m2 * m3
      m1a   (itape)     =     m1
      m2a   (itape)     =     m2
      m3a   (itape)     =     m3
      mba   (itape)     =     m1 * m2
      mpa   (itape)     =     mp
      mta   (itape)     =     mt


!     calculate recordlenght and open file
!     mreca is recordlenght in numbers of specified type
!     mrec is actual recordlenght in units according to machine

!---rebuild fname
      do i = 1,124
      if (fbase(i:i) == ' ') then
          idot = i
          exit
      end if
      end do
      fname(1:idot+2) = fbase(1:idot-1)//fext(imach)
      do i = idot+3,127
         fname(i:i) = ' '
      end do


!---unformatted direct access filetype
            if (imacha(itape) .eq. 1) then
!             sun: real = 4 bytes, 32 kbytes max. recordlength
!             sgi: real = 1 word, 8 kwords max. recordlength
              mreca(itape) = m1a(itape)
!SGI              mrec  = mreca(itape) * 4
              mrec  = mreca(itape) *recl
            else if (imacha(itape) .eq. 2) then
!             sun: real = 8 bytes, 32 kbytes max. recordlength
!             sgi: real = 2 words, 8 kwords max. recordlength
              mreca(itape) = m1a(itape)
!SGI              mrec  = mreca(itape) * 8
              mrec  = mreca(itape) * 2 *recl
            end if

!--- open file
            open(unit=iunit,err=2000,file=fname,status='old',  &
!                access='DIRECT',form='unformatted/ieee',     &
    	      access='DIRECT',form='unformatted',	       &
    	      iostat=ierror,recl=mrec)
      return

!---error procedure and end
2000  write(*,*)'   ERROR OPENING UNIT ',iunit,', ',fname
      write(*,*)'   IOSTAT = ',ierror
      write(*,*)'   possible cause: file does not exist'
      stop 'READCD'
end subroutine

!------------------------------------------------------------------------------!
! subroutine: reads .s4 and .s8 file (input in 2d slices)
!
!     is called when reading actual data.
!     It determines the chosen data format and calls the appropriate
!     reading routine
!
!------------------------------------------------------------------------------!
      subroutine readd1d (itape,data,nt,np,nk,nj,ierror)

!---variables
      integer                      :: nt,np,nk,nj,ierror
      real*8,dimension(m1a(itape)) ::data
      integer :: itape

!---error check
      if (nt .le. 0) then
        write (*,1000) iunita(itape),nt
        stop 'READD'
      else if (nt .gt. mta(itape) ) then
        write (*,1010) iunita(itape), nt, mta(itape)
        stop 'READD'
      end if
      if (np .le. 0) then
        write (*,1020) iunita(itape),np
        stop 'READD'
      else if (np .gt. mpa(itape) ) then
        write (*,1030) iunita(itape), np, mpa(itape)
        stop 'READD'
      end if
      if (nk .le. 0) then
        write (*,1020) iunita(itape),nk
        stop 'READD'
      else if (nk .gt. m3a(itape) ) then
        write (*,1030) iunita(itape), nk, m3a(itape)
        stop 'READD'
      end if
      if (nj .le. 0) then
        write (*,1020) iunita(itape),nj
        stop 'READD'
      else if (nj .gt. m2a(itape) ) then
        write (*,1030) iunita(itape), nj, m2a(itape)
        stop 'READD'
      end if
 1000 format (3x,'ERROR: reading from unit ',i2/         &
            3x,'timestep-index',i10,' is less then 1')
 1010 format (3x,'ERROR: reading from unit ',i2/ &
              3x,'timestep-index',i10,' is greater then ',i10)
 1020 format (3x,'ERROR: reading from unit ',i2/ &
            3x,'array number',i3,' is less then 1')
 1030 format (3x,'ERROR: reading from unit ',i2/ &
           3x,'array number',i3,' is greater then ',i3)


   if (imacha(itape) == 1) then
     call read41d(itape,data,nt,np,nk,nj,ierror)
   else if (imacha(itape) == 2) then
     call read81d(itape,data,nt,np,nk,nj,ierror)
   else
     stop 'READD'
   end if
end subroutine

!------------------------------------------------------------------------------!
! subroutine: reads .s8 file (input in 2d slices)
!
!     reads real*8 data from unformatted direct access files.
!
!------------------------------------------------------------------------------!
subroutine read81d (itape,data,nt,np,nk,nj,ierror)

!---variables
      integer                      :: nt,np,nk,nj,ierror
      real*8,dimension(m1a(itape)) :: data
      integer :: adab2,adrec1
      integer*8 :: adab1
      integer :: itape

!---calculate adresses
      adab1       =     (nt - 1) * mpa(itape) * maa(itape) &
                  +     (np - 1) * maa(itape)              &
		  +     (nk - 1) * m1a(itape)*m2a(itape)   &
		  +     (nj - 1) * m1a(itape) + 1
      adab2 =     adab1 + m1a(itape) - 1
      adrec1=     int((adab1-1)/mreca(itape)) + 1
!      adpos1=     mod((adab1-1),mreca(itape)) + 1
!      adrec2=     int((adab2-1)/mreca(itape)) + 1
!     adpos2=     mod((adab2-1),mreca(itape)) + 1
!      j1    =     1
!      j2    =     mreca(itape) - adpos1 + 1

            read(iunita(itape),rec=adrec1,err=2000,iostat=ierror) data

      return
2000  write (*,*) ' ERROR READING FROM UNIT ',iunita(itape)
      write(*,*)'   IOSTAT = ',ierror
      stop 'READ8'
end subroutine

!------------------------------------------------------------------------------!
! subroutine: reads .s4 file (input in 2d slices)
!
!     reads real*8 data from unformatted direct access files.
!
!------------------------------------------------------------------------------!
subroutine read41d (itape,data,nt,np,nk,nj,ierror)

!---variables
      integer                      :: nt,np,nk,nj,ierror
      real*8,dimension(m1a(itape)) :: data
      real*4,dimension(m1a(itape)) :: data4
      integer :: adab2,adrec1
      integer*8 :: adab1
      integer :: itape

!---calculate adresses
      adab1       =     (nt - 1) * mpa(itape) * maa(itape) &
                  +     (np - 1) * maa(itape)              &
		  +     (nk - 1) * m1a(itape)*m2a(itape)   &
		  +     (nj - 1) * m1a(itape) + 1
      adab2 =     adab1 + m1a(itape) - 1
      adrec1=     int((adab1-1)/mreca(itape)) + 1
!      adpos1=     mod((adab1-1),mreca(itape)) + 1
!      adrec2=     int((adab2-1)/mreca(itape)) + 1
!     adpos2=     mod((adab2-1),mreca(itape)) + 1
!      j1    =     1
!      j2    =     mreca(itape) - adpos1 + 1


      read(iunita(itape),rec=adrec1,err=2000,iostat=ierror) data4
      data=data4
      return
2000  write (*,*) ' ERROR READING FROM UNIT ',iunita(itape)
      write(*,*)'   IOSTAT = ',ierror
      stop 'READ8'
end subroutine



end module ios_unstrc
