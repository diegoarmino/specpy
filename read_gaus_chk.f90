program read_gaus
   integer,parameter          :: nclatoms=0
   integer          :: nqmatoms,mode
   integer,allocatable           :: at_numbers(:) ! Atomic numbers of QM atoms.
   double precision,allocatable  :: qmcoords(:,:) ! QM atom coordinates
   double precision,allocatable  :: atmass(:)
   double precision,allocatable  :: nmodes(:,:)
   double precision,allocatable  :: dxyzqm(:,:)
   double precision,allocatable  :: freq(:)

   character(100)   :: numchar
   character(100)   :: modechar

   call GET_COMMAND_ARGUMENT(1,numchar)
   read(numchar,*) nqmatoms

   allocate( at_numbers(nqmatoms), qmcoords(3,nqmatoms) , atmass(nqmatoms),  &
             nmodes(3*nqmatoms-6,3*nqmatoms), dxyzqm(3,nqmatoms),            &
             freq(3*nqmatoms-6) )

   call read_gaus_chekpoint(nqmatoms,nclatoms,at_numbers,qmcoords,atmass,nmodes,dxyzqm,freq)

end program read_gaus


   subroutine read_gaus_chekpoint(nqmatoms,nclatoms,at_numbers,qmcoords,atmass,nmodes,dxyzqm,freq)
      implicit none

      integer,          intent(in)  :: nqmatoms, nclatoms
      integer,          intent(in)  :: at_numbers(nqmatoms) ! Atomic numbers of QM atoms.
      double precision, intent(out) :: qmcoords(3,nqmatoms) ! QM atom coordinates
      double precision, intent(out) :: atmass(nqmatoms)
      double precision, intent(out) :: nmodes(3*nqmatoms,3*nqmatoms-6)
      double precision, intent(out) :: dxyzqm(3,nqmatoms)
      double precision, intent(out) :: freq(3*nqmatoms-6)

      double precision              :: gaus_nm(3*nqmatoms,3*nqmatoms-6)
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
      double precision              :: tmp
      character (len=4)             :: cfx,cfy,cfz
      character (len=17)            :: keycoords
      character (len=19)            :: keymass
      character (len=9)             :: keynmodes
      character (len=18)            :: keygrad
      character (len=6)             :: keyfreq
      integer                       :: estat,ios,exitstat
      integer                       :: i,j,k,r,s,m,nm,at,nat,count,xi
      integer                       :: ndf, nvdf


!     Constants
      double precision, parameter :: hbarjs    =  1.05457187d-34        ! J.s
      double precision, parameter :: hjs       =  6.62607015d-34        ! J.s
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

      character(100)                :: modechar
      integer                       :: mode

!     Allocatable variables
      double precision,ALLOCATABLE  :: nmd(:)
!     External functions
      double precision, external :: dnrm2
      double precision, external :: ddot

      call GET_COMMAND_ARGUMENT(2,modechar)
      read(modechar,*) mode

      gfac = he2j/(ccm**2 * b2m**2 * kg2da)
      b2msq = b2m**2
      ccmsq = ccm**2
      hbar  = hbarjs/(ccm * b2m**2 * kg2da)
!      hbar  = hjs/(ccm * b2m**2 * kg2da)
      hbarsqrt = sqrt(hbar)

      ndf=3*nqmatoms
      nvdf=ndf-6

!     Opening checkpoint file.
      open(unit=234, file="freq.fchk", iostat=ios, action="read")
      if ( ios /= 0 ) stop "Error opening freq.fchk file to read."
      open(unit=235, file="td.fchk", iostat=ios, action="read")
      if ( ios /= 0 ) stop "Error opening td.fchk file to read."
      if (mode==1) then
         open(unit=236, file="td-opt.fchk", iostat=ios, action="read")
         if ( ios /= 0 ) stop "Error opening fchk file to read."
      endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     MASSES
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      read(234,'(A18)') keymass
      count=1
      do while (trim(keymass) /= 'Real atomic weight')
         read(234,'(A18)') keymass
         count=count+1
      end do
      write(77,*) 'Number of lines skipped in fchk file is = ', count
      read(234,*) ( atmass(i),i=1,nqmatoms )
      rewind 234

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     FREQUENCIES
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      rewind 234
      read(234,'(A6)') keyfreq
      count=1
      do while (trim(keyfreq) /= 'Vib-E2')
         read(234,'(A6)') keyfreq
         count=count+1
      end do
      write(77,*) 'Number of lines skipped in fchk file is = ', count
      read(234,*) ( freq(j),j=1,nvdf )

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     NORMAL MODES
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      read(234,'(A9)') keynmodes
      count=1
      do while (trim(keynmodes) /= 'Vib-Modes')
         read(234,'(A9)') keynmodes
         count=count+1
      end do

      allocate(nmd(ndf*nvdf))
      write(77,*) 'NMODES Number of lines skipped in fchk file is = ', count
      read(234,*) ( nmd(j),j=1,ndf*nvdf )

      k=1
      do m=1,nvdf
         do i=1,ndf
            nmodes(i,m) = nmd(k)
            k=k+1
         end do
      end do
      deallocate(nmd)

 !     UN-MASS-WEIGHT NORMAL MODES
 !     Normal modes come mass weighted from Gaussian.
 !     To be more precise, l = M^-1/2 L
 !     So we have to moltiply them by M^1/2 in order to obtain
 !     pure normal modes.

       gaus_nm=nmodes
       do nm=1,nvdf
          do at=1,nqmatoms
             do j=1,3
                k=j+3*(at-1)
                nmodes(k,nm) = nmodes(k,nm)*Sqrt(atmass(at))
             end do
          end do
       end do

 !     NORMALIZE NORMAL MODES
       do nm=1,nvdf
          nmt=nmodes(:,nm)
          nmt=nmt/dnrm2(ndf,nmt,1)
          nmodes(:,nm)=nmt
       end do

       call dgemm('T','N',nvdf,nvdf,ndf,1d0,nmodes,ndf,nmodes,ndf,0d0,orth,nvdf)
       write(77,'(A)') 'DEBUG> Writing normal modes product matrix (should be identity)'
       do i=1,nvdf
          write(77,'(99F8.3)') orth(i,:)
       end do

      write(77,*) 'Frequencies= ', freq
      write(77,*)
      write(77,*)             '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(77,*)             '           WRITING NORMAL MODES            '
      write(77,*)             '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(77,*)
      do j=1,nvdf
         write(77,'(A,i3)')   'MODE =',j
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

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     GRADIENTS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!     Read gradient
      rewind 235
      read(235,'(A18)') keygrad
      count=1
      do while (trim(keygrad) /= 'Cartesian Gradient')
         read(235,'(A18)') keygrad
         count=count+1
      end do
      read(235,*) ( ( dxyzqm(j,i),j=1,3 ), i=1,nqmatoms )

      !write(77,*) 'Gradient Vector = \n', dxyzqm

!     Mass weight gradient
      do i=1,nqmatoms
          do j=1,3
              gtmp(3*(i-1)+j)=dxyzqm(j,i)/sqrt(atmass(i))
          end do
      end do

!     Convert gradient to mass-weighted normal coordinates
      grad0=0d0
      do nm = 1,nvdf
         tmp = 0d0
         do xi = 1,ndf
            grad0(nm) = grad0(nm) + nmodes(xi,nm) * gtmp(xi)
         end do
      end do

!     Convert gradient to adimensional coordinates
      !grad0 = hbarsqrt * gfac * grad0 / Sqrt(freq) ! DANGER
      grad0 = gfac * grad0 / Sqrt(freq)

!     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!                   COMPUTING DISPLACEMENTS FROM GRADIENTS
!     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      delta  = - grad0 / (freq * hbar)

!     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!         COMPUTING DISPLACEMENTS FROM EXCITED AND GROUND STATE GEOMETRIES.
!     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      if (mode==1) then
   !     READING GROUND STATE EQUILIBRIUM GEOMETRY
   !     --------------------------------------------------------------------------
         rewind 235
         read(235,'(A17)') keycoords
         count=1
         do while (trim(keycoords) /= 'Current cartesian')
            read(235,'(A17)') keycoords
            count=count+1
         end do
         write(77,*) 'Number of lines skipped in fchk file is = ', count
         read(235,*) ( ( qmcoords(i,j),i=1,3 ), j=1,nqmatoms )
         qmcoords=qmcoords/a0  !DANGER a0 should not be there
         write(77,*) qmcoords
         rewind 235


   !     READING EXCITED STATE EQUILIBRIUM GEOMETRY (TD OPTIMIZATION)
   !     --------------------------------------------------------------------------
         read(236,'(A17)') keycoords
         count=1
         do while (trim(keycoords) /= 'Current cartesian')
            read(236,'(A17)') keycoords
            count=count+1
         end do
         write(77,*) 'Number of lines skipped in fchk file is = ', count
         read(236,*) ( ( excoords(i,j),i=1,3 ), j=1,nqmatoms )
         excoords=excoords/a0  !DANGER a0 should not be there
         write(77,*) excoords
         rewind 236


   !     TRANSLATE TO CENTER OF MASS AND ROTATE TO PRINCIPAL AXES OF INERTIA.
   !     --------------------------------------------------------------------------
         call orient(qmcoords,at_numbers,atmass,ndf,nqmatoms)
         call orient(excoords,at_numbers,atmass,ndf,nqmatoms)


   !     TAKING DIFFERENCE BETWEEN EXCITED AND GROUND STATE GEOMETRIES IN CARTESIAN
   !     --------------------------------------------------------------------------
   !      del0 = excoords-qmcoords
   !      do i=1,nqmatoms
   !          do j=1,3
   !              del0(j,i)         = qmcoords(j,i) - excoords(j,i)
   !          end do
   !      end do

   !     NONDIMENSIONALIZATION OF GEOMETRIES AND DISPLACEMENTS
   !     --------------------------------------------------------------------------

   !     MASS WEIGHT
         do i=1,nqmatoms
             do j=1,3
                 gnd_coords_v(3*(i-1)+j) = qmcoords(j,i) * sqrt(atmass(i))
                 exc_coords_v(3*(i-1)+j) = excoords(j,i) * sqrt(atmass(i))
                 !delv(3*(i-1)+j)         = del0(j,i)     * sqrt(atmass(i))
                 delv(3*(i-1)+j)         = (qmcoords(j,i) - excoords(j,i)) * sqrt(atmass(i))
             end do
         end do

   !     CONVERT TO MASS-WEIGHTED NORMAL COORDINATES
         gnd_nc=0d0
         exc_nc=0d0
         del_nc=0d0
         do nm = 1,nvdf
            do xi = 1,ndf
               gnd_nc(nm) = gnd_nc(nm) + nmodes(xi,nm) * gnd_coords_v(xi)
               exc_nc(nm) = exc_nc(nm) + nmodes(xi,nm) * exc_coords_v(xi)
               del_nc(nm) = del_nc(nm) + nmodes(xi,nm) * delv(xi)
            end do
         end do

   !     CONVERT TO ADIMENSIONAL NORMAL COORDINATES
         gnd_nc = gnd_nc * Sqrt(freq/hbar)
         exc_nc = exc_nc * Sqrt(freq/hbar)
         del_nc = del_nc * Sqrt(freq/hbar)

         write(77,*) 'Coords of ground state'
         write(77,'(999F8.3)') qmcoords
         write(77,*) 'Coords of excited state'
         write(77,'(999F8.3)') excoords
         write(77,*) 'Coords of delta'
         write(77,'(999F8.3)') del0
         write(77,*) 'Normal coords of ground state'
         write(77,'(999F8.3)') gnd_nc
         write(77,*) 'Normal coords of excited state'
         write(77,'(999F8.3)') exc_nc

         del=gnd_nc-exc_nc

         write(77,*) 'Displacements in adimensional normal cartesian coordinates'
         write(77,*) '    LMXg-LMXex     frequency (cm-1)'
         do i=1,nvdf
           !write(77,'(4f14.6)') 2d0*del_nc(i), 2d0*del(i), 2d0*delta(i),freq(i)
           !write(77,'(4f14.6)') del_nc(i), del(i), delta(i)/a0,freq(i) ! DANGER, A0 MAYBE SHOULD NOT BE HERE
           write(77,'(4f14.6)') del_nc(i), freq(i) ! DANGER, A0 MAYBE SHOULD NOT BE HERE
         end do
      end if


        write(77,'(4f14.6)')
        write(77,*) '    gradient     frequency (cm-1)'
      do i=1,nvdf
        write(77,'(4f14.6)') delta(i)/a0, freq(i) ! DANGER, A0 MAYBE SHOULD NOT BE HERE
      end do

      close(234)
      close(235)
      if (mode==1) then
         close(236)
      end if

   end subroutine


  subroutine orient(qmcoords,at_numbers,atmass,ndf,nqmatoms)
        implicit none
   !    ------------------------------------------------------------------
        integer,intent(in)                :: ndf
        integer,intent(in)                :: nqmatoms
        integer,intent(in)                :: at_numbers(nqmatoms) ! Atomic numbers of QM atoms.
        double precision,intent(in)       :: atmass(nqmatoms)
        double precision,intent(inout)    :: qmcoords(3,nqmatoms)
   !    Internal variables
        integer             :: i,j,k,iat
        integer             :: zeros(2),nonz(2)
        integer             :: ipiv(3)
        real*8              :: Mass(ndf)
        real*8              :: X0(3,nqmatoms)
        real*8              :: Xrot(3,nqmatoms)
        real*8              :: Itsr(3,3)
        real*8              :: innp(3,3)
        real*8              :: outp(3,3)
        real*8              :: ai(3)
        real*8              :: atmp
        real*8              :: cutoff
        real*8              :: det
        real*8              :: trp
        real*8              :: totM
        real*8              :: com(3)
        real*8              :: Ieigval(3)

        real*8,external  :: ddot

        real*8,  dimension(:,:), allocatable :: VT,U

        real*8,  dimension(:), allocatable :: WORK,WORK2,IWORK,IWORK2
        integer                            :: LWORK,LIWORK,INFO
!        include "qvmbia_param.f"
   !    ------------------------------------------------------------------

        !unitv=reshape((/ 1d0, 0d0, 0d0, &
        !             &   0d0, 1d0, 0d0, &
        !             &   0d0, 0d0, 1d0 /),(/3,3/))

   !    Converting initial coordinates to AU
        X0=qmcoords!/a0
   !     write(77,'(A)') 'CONVERSION FACTOR a0 -> Angs'
   !     write(77,'(D12.3)') a0

   !    Computing total mass and atomic mass array and center of mass.
        totM=0.0d0
        do i=1,nqmatoms
            do j=1,3
                com(j)=com(j)+atmass(i)*X0(j,i)
            end do
            totM=totM+atmass(i)
        end do
        com = com/totM

   !    Translating to center of mass
        do i=1,nqmatoms
           do j=1,3
              X0(j,i)=(X0(j,i)-com(j))*sqrt(atmass(i))
           end do
        end do


   !    Moment of inertia tensor.
   !    Itsr = Sum_i [ai**T.ai - ai.ai**T]
        Itsr=0.0d0
        do i=1,nqmatoms
           ai=X0(:,i)
           innp=0.0d0
           outp=0.0d0
           innp(1,1) = ddot(3,ai,1,ai,1)
           innp(2,2) = innp(1,1)
           innp(3,3) = innp(1,1)
           call dger(3,3,1d0,ai,1,ai,1,outp,3)
   !        call dsyr('U',3,1d0,ai,1,outp,3)
           Itsr=Itsr+innp-outp
        end do

   !    Symmetrizing inertia tensor.
        do i=1,2
           do j=i+1,3
             Itsr(j,i)=Itsr(i,j)
           end do
        end do

        write(77,'(A)') 'INERTIA TENSOR'
        do i=1,3
           write(77,'(3D11.2)') Itsr(i,:)
        end do

  !     Scaling tensor of inertia to avoid numerical problems.
        Itsr=Itsr*1000d0


  !     DIAGONALIZATION OF THE MOMENT OF INERTIA
        allocate ( work(100), IWORK(100) )
        LWORK=-1
        call dsyevd('V','U',3,Itsr,3,Ieigval,WORK,LWORK,IWORK,LWORK,INFO)
        LWORK=WORK(1)
        LIWORK=IWORK(1)
        if(allocated(WORK2)) deallocate (WORK2,IWORK2)
        allocate (WORK2(LWORK),IWORK2(LIWORK))
        call dsyevd('V','U',3,Itsr,3,Ieigval,WORK2,LWORK,IWORK2,LIWORK,INFO)
        deallocate( WORK, IWORK, WORK2, IWORK2 )


  !     Descaling eigenvalues.
        Ieigval=Ieigval/1000d0

  !     ROTATE TO PRINCIPAL AXES
        Xrot=0d0
        do iat = 1,nqmatoms
          do i = 1,3
            do j = 1,3
              Xrot(i,iat) = Xrot(i,iat) + Itsr(j,i) * X0(j,iat)
            end do
          end do
        end do

        X0 = Xrot

  end subroutine
