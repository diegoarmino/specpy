program read_gaus
   integer,parameter          :: nqmatoms=14, nclatoms=0
   integer          :: at_numbers(nqmatoms) ! Atomic numbers of QM atoms.
   double precision :: qmcoords(3,nqmatoms) ! QM atom coordinates
   double precision :: atmass(nqmatoms)
   double precision :: nmodes(3*nqmatoms-6,3*nqmatoms)
   double precision :: dxyzqm(3,nqmatoms)
   double precision :: freq(3*nqmatoms-6)
   call read_gaus_chekpoint(nqmatoms,nclatoms,at_numbers,qmcoords,atmass,nmodes,dxyzqm,freq)

end program read_gaus


   subroutine read_gaus_chekpoint(nqmatoms,nclatoms,at_numbers,qmcoords,atmass,nmodes,dxyzqm,freq)
      implicit none

      integer,          intent(in)  :: nqmatoms, nclatoms
      integer,          intent(in)  :: at_numbers(nqmatoms) ! Atomic numbers of QM atoms.
      double precision, intent(out) :: qmcoords(3,nqmatoms) ! QM atom coordinates
      double precision, intent(out) :: atmass(nqmatoms)
      double precision, intent(out) :: nmodes(3*nqmatoms-6,3*nqmatoms)
      double precision, intent(out) :: dxyzqm(3,nqmatoms)
      double precision, intent(out) :: freq(3*nqmatoms-6)

      double precision              :: excoords(3,nqmatoms) ! QM atom coordinates
      double precision              :: del0(3,nqmatoms) ! QM atom coordinates
      double precision              :: delv(3*nqmatoms) ! QM atom coordinates
      double precision              :: del_nc(3*nqmatoms-6) ! QM atom coordinates
      double precision              :: gnd_coords_v(3*nqmatoms) ! QM atom coordinates
      double precision              :: exc_coords_v(3*nqmatoms) ! QM atom coordinates
      double precision              :: gnd_nc(3*nqmatoms-6)
      double precision              :: exc_nc(3*nqmatoms-6)
      double precision              :: del(3*nqmatoms-6)
      double precision              :: orth(3*nqmatoms-6,3*nqmatoms-6)
      double precision              :: gnc(3*nqmatoms-6)
      double precision              :: grad0(3*nqmatoms-6)
      double precision              :: delta(3*nqmatoms-6)
      double precision              :: nmt(3*nqmatoms)
      double precision              :: gtmp(3*nqmatoms)
      character (len=4)             :: cfx,cfy,cfz
      character (len=17)            :: keycoords
      character (len=19)            :: keymass
      character (len=9)             :: keynmodes
      character (len=18)            :: keygrad
      character (len=6)             :: keyfreq
      integer                       :: estat,ios,exitstat
      integer                       :: i,j,k,r,s,m,nm,at,nat,count
      integer                       :: ndf, nvdf


!     Constants
      double precision, parameter :: hbarjs    =  1.05457187d-34        ! J.s
      double precision, parameter :: ccm       =  2.99792458d10         ! Speed of light cm/s
      double precision, parameter :: A0        =  0.52917710d00         ! Angstrms/bohr
      double precision, parameter :: he2j      =  4.35974394d-18        ! hartree/J
      double precision, parameter :: b2m       =  5.29177249d-11        ! meter/bohr
      double precision, parameter :: kg2da     =  1.66054020d-27        ! Kg/Daltons


!      double precision, parameter :: gfac      = 1.0431997d9            ! Factor convert nondimensionalized Gk
                                                                        ! (units of Hartree) into units of
                                                                        ! Da a0^2 cm-2
      double precision :: gfac
      double precision :: hbarsqrt                           ! h^1/2 bar in Da a0^2 cm-1
      double precision :: b2msq                              ! meter^2 / bohr^2
      double precision :: ccmsq                              ! Speed of light cm^2/s^2
      double precision :: hbar                               ! h bar in Da a0^2 cm-1

!     Allocatable variables
      double precision,ALLOCATABLE  :: nmd(:)
!     External functions
      double precision, external :: dnrm2
      double precision, external :: ddot

      gfac = he2j/(ccm**2 * b2m**2 * kg2da)
      b2msq = b2m**2
      ccmsq = ccm**2
      hbar  = hbarjs/(ccm * b2m**2 * kg2da)
      hbarsqrt = sqrt(hbar)
      write(*,*) 'hbar = ', hbar, 'hbarsqrt = ', hbarsqrt

      ndf=3*nqmatoms
      nvdf=ndf-6

!     Opening checkpoint file.
      open(unit=234, file="freq.fchk", iostat=ios, action="read")
      if ( ios /= 0 ) stop "Error opening fchk file to read."
      open(unit=235, file="td.fchk", iostat=ios, action="read")
      if ( ios /= 0 ) stop "Error opening fchk file to read."
      open(unit=236, file="td-opt.fchk", iostat=ios, action="read")
      if ( ios /= 0 ) stop "Error opening fchk file to read."

!     Reading Cartesian Coordinates.
      read(234,'(A17)') keycoords
      count=1
      do while (trim(keycoords) /= 'Current cartesian')
         read(234,'(A17)') keycoords
         count=count+1
      end do
      write(77,*) 'Number of lines skipped in fchk file is = ', count
      read(234,*) ( ( qmcoords(i,j),i=1,3 ), j=1,nqmatoms )
      qmcoords=qmcoords*a0
      write(77,*) qmcoords
      rewind 234

!     Reading Masses
      read(234,'(A18)') keymass
      count=1
      do while (trim(keymass) /= 'Real atomic weight')
         read(234,'(A18)') keymass
         count=count+1
      end do
      write(77,*) 'Number of lines skipped in fchk file is = ', count
      read(234,*) ( atmass(i),i=1,nqmatoms )
      rewind 234

!     Reading Normal Modes.
      read(234,'(A9)') keynmodes
      count=1
      do while (trim(keynmodes) /= 'Vib-Modes')
         read(234,'(A9)') keynmodes
         count=count+1
      end do

      allocate(nmd(ndf*nvdf))

      write(77,*) 'Number of lines skipped in fchk file is = ', count
      read(234,*) ( nmd(j),j=1,ndf*nvdf )

      k=1
      do m=1,nvdf
         do i=1,ndf
            nmodes(m,i) = nmd(k)
            k=k+1
         end do
      end do

      do nm=1,nvdf
         write(77,*)
         write(77,*) nmodes(nm,:)
      end do

      deallocate(nmd)

 !     UN-MASS-WEIGHT NORMAL MODES
 !     Normal modes come mass weighted from Gaussian.
 !     To be more precise, l = M^-1/2 L
 !     So we have to moltiply them by M^1/2 in order to obtain
 !     pure normal modes.

       do nm=1,nvdf
          do at=1,nqmatoms
             do j=1,3
                k=j+3*(at-1)
                nmodes(nm,k) = nmodes(nm,k)*Sqrt(atmass(at))
             end do
          end do
       end do
       write(77,*) (sqrt(atmass(at)),at=1,nqmatoms)


 !     NORMALIZE MASS-WEIGHTED NORMAL MODES
       do nm=1,nvdf
          nmt=nmodes(nm,:)
!          write(77,'(99D16.8)') nmodes(:,nm)
          write(77,*) 'Not normal'
          write(77,'(99D16.8)') nmt
!          write(77,'(F16.8)') dnrm2(ndf,nmt,1)
          nmt=nmt/dnrm2(ndf,nmt,1)
          write(77,*) 'Yes normal'
          write(77,'(99D16.8)') nmt
!          write(77,'(F16.8)') dnrm2(ndf,nmt,1)
          nmodes(nm,:)=nmt
       end do


       call dgemm('T','N',nvdf,nvdf,ndf,1d0,nmodes,ndf,nmodes,ndf,0d0,orth,nvdf)
       write(77,'(A)') 'DEBUG> Writing normal modes product matrix (should be identity)'
       do i=1,nvdf
          write(77,'(99F15.6)') orth(i,:)
       end do

      rewind 235
      read(235,'(A18)') keygrad
      count=1
      do while (trim(keygrad) /= 'Cartesian Gradient')
         read(235,'(A18)') keygrad
         count=count+1
      end do
      write(77,*) 'Number of lines skipped in fchk file is = ', count
      read(235,*) ( ( dxyzqm(j,i),j=1,3 ), i=1,nqmatoms )

      write(77,*) 'Gradient Vector = \n', dxyzqm

      rewind 234
      read(234,'(A6)') keyfreq
      count=1
      do while (trim(keyfreq) /= 'Vib-E2')
         read(234,'(A6)') keyfreq
         count=count+1
      end do
      write(77,*) 'Number of lines skipped in fchk file is = ', count
      read(234,*) ( freq(j),j=1,nvdf )

      write(77,*) 'Frequencies= ', freq
      write(77,*)
      write(77,*)             '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(77,*)             '           WRITING NORMAL MODES            '
      write(77,*)             '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(77,*)
      do j=7,ndf
         write(77,'(A,i3)')   'MODE =',j-6
         write(77,'(A,F9.2)') 'FREQUENCY (CM-1) =',freq(j)
         write(77,'(A)')      '-------------------------------------------'
         write(77,'(A)')      ' AT       X            Y            Z      '
         write(77,'(A)')      '-------------------------------------------'
         do k=1,nqmatoms
            r=3*(k-1)+1
            s=r+2
            write(77,'(i3,3F13.8)')  k, nmodes(r:s,j)
         end do
         write(77,*)
         write(77,*)
      end do


      !-------------------------------------------------------------------------------
      !     GRADIENT
      !-------------------------------------------------------------------------------
!     Mass weight gradient
      do i=1,nqmatoms
          do j=1,3
              gtmp(3*(i-1)+j)=dxyzqm(j,i)/sqrt(atmass(i))
          end do
      end do

!     Convert Gradient to normal coordinates.
      gtmp=gtmp*1d6
      call dgemv('T',ndf,nvdf,1d0,nmodes,ndf,gtmp,1,0d0,gnc,1)
      grad0=gnc*1d-6   ! Gradient without electric field.

!     Gradient is considered zero if less than 10d-10
      do i=1,nvdf
         if (abs(grad0(i))<1d-10) THEN
            grad0(i) = 0d0
         else
            grad0(i) = hbarsqrt * gfac * grad0(i) / Sqrt(freq(i))
         end if
      end do


      write(77,*) 'Mass-Weighted Gradient in cartesians'
      write(77,*) gtmp
      write(77,*) 'Normal modes L ='
      write(77,*) nmodes
      write(77,*) 'Mass-weighted Gradient in normal coordinates with no electric field.'
      do i=1,nvdf
        write(77,*) grad0(i)
      end do
      !-------------------------------------------------------------------------------
      !     COMPUTING DISPLACEMENTS
      !-------------------------------------------------------------------------------

      delta  = - grad0 / freq
      delta  =   delta / hbar

      write(77,*) 'Displacements in adimensional normal cartesian coordinates'
      do i=1,nvdf
        write(77,'(2f14.6)') delta(i), freq(i)
      end do

!     Reading Cartesian Coordinates of TD excited state optimization.
      read(236,'(A17)') keycoords
      count=1
      do while (trim(keycoords) /= 'Current cartesian')
         read(236,'(A17)') keycoords
         count=count+1
      end do
      write(77,*) 'Number of lines skipped in fchk file is = ', count
      read(236,*) ( ( excoords(i,j),i=1,3 ), j=1,nqmatoms )
      excoords=excoords*a0
      write(77,*) excoords
      rewind 236

      del0 = excoords-qmcoords
      !-------------------------------------------------------------------------------
      !     CONVERT CARTESIAN COORDS INTO ADIMENSIONAL NORMAL COORDS
      !-------------------------------------------------------------------------------
!     Mass weight
      do i=1,nqmatoms
          do j=1,3
              gnd_coords_v(3*(i-1)+j)=qmcoords(j,i)*sqrt(atmass(i))
              exc_coords_v(3*(i-1)+j)=excoords(j,i)*sqrt(atmass(i))
              delv(3*(i-1)+j)=del0(j,i)*sqrt(atmass(i))
          end do
      end do

      write(77,'(A)') 'MASAS ATOMICAS'
      write(77,'(999F8.3)') atmass
      write(77,'(999F8.3)') sqrt(atmass)

!     Convert to normal coordinates.
      call dgemv('T',ndf,nvdf,1d0,nmodes,ndf,gnd_coords_v,1,0d0,gnd_nc,1)
      call dgemv('T',ndf,nvdf,1d0,nmodes,ndf,exc_coords_v,1,0d0,exc_nc,1)
      call dgemv('T',ndf,nvdf,1d0,nmodes,ndf,delv,1,0d0,del_nc,1)

!     Gradient is considered zero if less than 10d-10
      do i=1,nvdf
         gnd_nc(i) = gnd_nc(i) * Sqrt(freq(i))/hbarsqrt
         exc_nc(i) = exc_nc(i) * Sqrt(freq(i))/hbarsqrt
         del_nc(i) = exc_nc(i) * Sqrt(freq(i))/hbarsqrt
      end do


      write(77,*) 'Normal coords of ground state'
      write(77,'(999F8.3)') gnd_nc
      write(77,*) 'Normal coords of excited state'
      write(77,'(999F8.3)') exc_nc

      del=gnd_nc-exc_nc

      write(77,*) 'Displacements in adimensional normal cartesian coordinates'
      write(77,*) 'difference    gradient     frequency (cm-1)'
      do i=1,nvdf
        write(77,'(4f14.6)') del_nc(i), del(i), delta(i),freq(i)
      end do

      close(234)
      close(235)
      close(236)

   end subroutine
