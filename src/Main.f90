!
!*********************************************************************************
!** Compute conductivity **
!*********************************************************************************
!
!    
 program Transport
  use PrecisionMod
  use MathConstantsMod
  use ParallelMod
  use LibMod
  use InputMod
  use KpointsMod
  use HamiltonianMod
  use CurrentVertexMod
  use SelfEnMod
  use OpticMod
  implicit none
!
  include 'mpif.h'
!
  character(len=2) :: indx1
  character(len=3) :: indx2
  character(len=4) :: indx4
  character(len=22) :: form
  character(len=100) :: FileName
  integer(i4b) :: nk, kindex, ns, nw, Niwn, Nipn, Nivn, Nv, error, Nsize, n, m
  integer(i4b) :: Norb1, Norb2, Norb3, Norb4, orbindex1, orbindex2
  integer(i4b) :: i, j, ii, jj
  real(dp) :: Omega_n
  complex(dpc), dimension(:,:,:,:), allocatable :: Hamk
  complex(dpc), dimension(:,:,:,:), allocatable :: SelfEn   
  complex(dpc), dimension(:,:,:,:,:), allocatable :: BareVertex
! 
  integer(i4b), dimension(2) :: dir
  real(dp), dimension(:), allocatable :: val, val2
  complex(dpc), dimension(:), allocatable :: vec1, vec2
  complex(dpc), dimension(:,:), allocatable :: mat, identity
  complex(dpc), dimension(:,:,:,:), allocatable :: LocalGFn, LocGFn
  complex(dpc), dimension(:,:), allocatable :: GFn
  real(dp), dimension(:,:,:), allocatable :: Sigma, Sigma1, Sigma2
  complex(dpc), dimension(:,:,:), allocatable :: Pol, Pol1, Pol2  
  complex(dpc), dimension(:,:), allocatable :: Polarization
!
!--variable for parallel run  
  integer(i4b) :: anid, ierr, isender, numsent
  integer(i4b), parameter :: senddatatag = 2013
  integer(i4b), parameter :: returndatatag = 2014
  integer(i4b), dimension(MPI_STATUS_SIZE) :: status
!
  real(dp) :: starttime, endtime
!
  call InitParallel
!  
  call walltime (starttime)
!
  if(MYID.eq.0 .and. NPE.gt.1) then
    write(*,'(/,a,i5,a11)') ' MPI-parallel calculation using ', NPE, ' processors'
  endif
!  
!******************************
!** Read system's parameters ** 
!******************************
!
  write(*,'(/,a15,2x,i4)') 'Reading Inputs:', Myid
  call InputInit()
!
!
!****************************
!** reading the input data **
!****************************
!
  write(*,'(/,a17,2x,i4)') 'Reading K-points:', Myid
  call KpointsInit(casename)
!
  allocate(Hamk( Numorbitals, Numorbitals, NumIrKpoints,NumSpin), stat=error)
  if(error .ne. 0) stop 'Main: error allocating Hamk'
  allocate(BareVertex(3, Numorbitals, Numorbitals, NumIrKpoints, NumSpin), stat=error)
  if(error .ne. 0) stop 'Main: error allocating BareVertex'
!
  allocate(SelfEn(Numorbitals, Numorbitals, 0:NiwnMax, NumSpin), stat=error)
    if(error .ne. 0) stop 'Main: error allocating SelfEn'
!
! read Ham(k)
  write(*,'(a8,2x,i4)') 'Read Hk:', Myid 
  call HkInit(casename, Numspin, Numorbitals, numIrkpoints, Hamk)
!
! Read Bare current vertices
  write(*,'(a8,2x,i4)') 'Read vk:', Myid
  call BareCurrentVertexInit(casename, Numspin, Numorbitals, NumIrkpoints, BareVertex)      
!
! Read Self energy 
  SelfEn = (0.0_dp,0.0_dp)
  if (InteractingSystem) then
    write(*,'(a12,2x,i4)') 'Read SelfEn:', Myid   
    call SelfEnInit2(casename, Numspin, Numorbitals, Niwnmax, InverseTemperature, SelfEn)
  end if 
!
  call barrier()
!
!---------------
!-- Local GFn --
!---------------
!
  write(*,'(a20,2x,i4)') 'calculate Local GFn:', Myid
  allocate(Identity(NumOrbitals,NumOrbitals), stat= error)
  if (error .ne. 0) stop 'Main: unable to allocate mat'
  Identity = (zero,zero)
  do Norb1=1, NumOrbitals
    Identity(Norb1,Norb1) = (one,zero)
  end do
!
    if (Myid .eq. 0) then
      allocate(LocalGFn(NumOrbitals, NumOrbitals, 0:NiwnMax, NumSpin), GFn(NumOrbitals, NumOrbitals), stat= error)
      if (error .ne. 0) stop 'Main: unable to allocate GFn'
      LocalGFn = (zero,zero)
      do ns = 1, NumSpin
        do Niwn = 0, NiwnMax
          Omega_n = (two*Niwn+one)*pi_d/InverseTemperature
          nk = 0
          do kindex = startkp(Myid), endkp(Myid)
            nk = nk+1
            GFn = (xi*Omega_n+ ChemicalPotential)*Identity-Hamk(:,:,kindex,ns)-SelfEn(:,:,niwn,ns)
            call InverseUsingLapack(GFn)
!
            LocalGFn(:,:,niwn,ns) = LocalGFn(:,:,Niwn,ns) + GFn(:,:)/NumIrkpoints
          end do
        end do
      end do
!     recive LocalGFn from computing nodes
      if(NPE .gt. 1) then
        allocate(LocGFn(NumOrbitals, NumOrbitals, 0:NiwnMax, NumSpin), stat= error)
        if (error .ne. 0) stop 'Main: unable to allocate LocGFn'
        Nsize = (NumOrbitals**2)*(NiwnMax+1)*NumSpin
        do Anid =1, NPE-1
          call MPI_RECV(LocGFn, Nsize, MPI_DOUBLE_COMPLEX, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, status, ierr)
          isender = status(MPI_SOURCE)
          LocalGFn = LocalGFn + LocGFn
        end do
	deallocate(LocGFn, stat= error)
        if (error .ne. 0) stop 'Main: unable to allocate LocGFn'
      end if
!     write LocalGFn
      filename = 'GFniwn.dat'
      open(100,file=trim(FileName), status='unknown')
      write(indx2,'(i3.3)') 2*NumSpin*(NumOrbitals**2)
      form = '(f16.10,'//indx2//'f24.16)'
      do Niwn = 0, NiwnMax
        Omega_n = (two*Niwn+one)*pi_d/InverseTemperature
        write(100,form) Omega_n,  (((LocalGFn(Norb1,Norb2,niwn,ns),Norb2=1,NumOrbitals),Norb1=1,NumOrbitals),ns=1,NumSpin)
      end do
      close(100)
      deallocate(LocalGFn, GFn, stat= error)
      if (error .ne. 0) stop 'Main: unable to allocate GFn'
    else  ! computing node
      allocate(LocalGFn(NumOrbitals, NumOrbitals, 0:NiwnMax, NumSpin), GFn(NumOrbitals, NumOrbitals), stat= error)
      if (error .ne. 0) stop 'Main: unable to allocate GFn'
      LocalGFn = (zero,zero)
      do ns = 1, NumSpin
        do Niwn = 0, NiwnMax
          Omega_n = (two*Niwn+one)*pi_d/InverseTemperature
          nk = 0
          do kindex = startkp(Myid), endkp(Myid)
            nk = nk+1
            GFn = (xi*Omega_n+ ChemicalPotential)*Identity-Hamk(:,:,kindex,ns)-SelfEn(:,:,niwn,ns)
            call InverseUsingLapack(GFn)
!
            LocalGFn(:,:,niwn,ns) = LocalGFn(:,:,Niwn,ns) + GFn(:,:)/NumIrkpoints
          end do
        end do
      end do
!     send the LocalGFn to master node
      Nsize = (NumOrbitals**2)*(NiwnMax+1)*NumSpin
      call MPI_SEND(LocalGFn, Nsize, MPI_DOUBLE_COMPLEX, 0, Myid, comm, ierr)
      deallocate(LocalGFn, GFn, stat= error)
      if (error .ne. 0) stop 'Main: unable to allocate GFn'
    end if
!
!
  deallocate(Identity, stat= error)
  if (error .ne. 0) stop 'Main: unable to deallocate Identity'
!
  call barrier()
!   
!
!*********************************************
!** Now actual calculation of the responses **
!*********************************************
!
!
  write(*,'(a19,2x,i4)') 'calculate response:', Myid
!
!-- Calculation of Polarization --
!
    if (MYId .eq. 0) then
!
        allocate( Pol(NivnMin:NivnMax, NumIrKpoints, NumSpin), Pol1(NivnMin:NivnMax, endkp(myid)-startkp(myid)+1,NumSpin), stat=error)
        if(error .ne. 0) stop 'Main: error allocating Pol'
        Pol = (zero,zero)
        Pol1 = (zero,zero)
        call OpticalPolarization (Nchoice, Numspin, NumOrbitals, numIrkpoints, NiwnMax, startkp(myid), endkp(myid), InverseTemperature, &
                                  ChemicalPotential, BareVertex, Hamk, SelfEn, NivnMin, NivnMax, Pol1)
        write(*,'(a8,2x,i4)') 'master: ', myid
        nk = 0
        do kindex = startkp(myid), endkp(myid)
          nk = nk + 1
          do nivn = NivnMin, NivnMax
            Pol(nivn,kindex,:) = Pol1(nivn, nk,:)
          end do
        end do
!
!-- recieving data from compute nodes --
! 
        if(NPE .gt. 1) then
          allocate( Pol2(NivnMin:NivnMax, endkp(NPE-1)-startkp(NPE-1)+1,NumSpin), stat=error)
          if(error .ne. 0) stop 'Main: error allocating Pol2'
          Nsize = (NivnMax-NivnMin+1)*(endkp(NPE-1)-startkp(NPE-1)+1)*NumSpin
          do Anid =1, NPE-1
            Pol2 = (zero,zero)
            call MPI_RECV(Pol2, Nsize, MPI_DOUBLE_COMPLEX, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, status, ierr)
            isender = status(MPI_SOURCE)
            write(*,'(a8,2x,i4)') 'Sender: ', isender
            nk = 0
            do kindex = startkp(isender), endkp(isender)
              nk = nk + 1
              do nivn = NivnMin, NivnMax
                Pol(nivn,kindex,:) = Pol2(nivn, nk,:)
              end do
            end do
          end do
	  deallocate( Pol2, stat=error)
          if(error .ne. 0) stop 'Main: error allocating Pol2'
        end if
        deallocate(Pol1, stat=error)
        if(error .ne. 0) stop 'Main: error deallocating Pol1'
!
    else         ! computing nodes
        allocate( Pol1(NivnMin:NivnMax, endkp(myid)-startkp(myid)+1,NumSpin), stat=error)
        if(error .ne. 0) stop 'Main: error allocating Pol1'
        Pol1 = (zero,zero)
        call OpticalPolarization (Nchoice, Numspin, NumOrbitals, numIrkpoints, NiwnMax, startkp(myid), endkp(myid), InverseTemperature, &
                                  ChemicalPotential, BareVertex, Hamk, SelfEn, NivnMin, NivnMax, Pol1)
!
!     send the Pol1 to master node
        Nsize =  (NivnMax-NivnMin+1)*(endkp(myid)-startkp(myid)+1)*NumSpin
        call MPI_SEND(Pol1, Nsize, MPI_DOUBLE_COMPLEX, 0, myid, comm, ierr)
!
        deallocate(Pol1, stat=error)
        if(error .ne. 0) stop 'Main: error deallocating Pol1'
    end if
!
!
!*****************************
!** Summation over k-points **
!*****************************
!
  if (MYId .eq. 0) then
        allocate( Polarization(NivnMin:NivnMax,NumSpin), stat=error)
        if(error .ne. 0) stop 'Main: error allocating Polarization'
        Polarization = (zero,zero)
!
        do nivn = NivnMin, NivnMax
          do nk = 1, NumIrKpoints
            Polarization(nivn,:) = Polarization(nivn,:) + Pol(nivn,nk,:)
          end do
        end do
        Polarization = Polarization/NumIrKpoints
!        if (NumSpin .eq. 1) Polarization = two*Polarization  ! counting for spin degeneracy
!
        fileName = 'OpticalPolarization.dat'
!
!     write Polarization
      open(100,file=trim(FileName), status='unknown')
        write(indx1,'(i2.2)') 2*NumSpin
        form = '(a32,'//indx1//'f24.16)'
!       Adding the diamagnetic term which is equal to -Paramagnetic term at v_n=0
        if (Nchoice .le. 3) write(100, form) '# diamagnetic contribution is : ', (Polarization(0,ns),ns=1,NumSpin) 
        form = '(f16.10,'//indx1//'f24.16)'
        do nivn = NivnMin, NivnMax
          Omega_n = two*nivn*pi_d/InverseTemperature
          write(100,form) Omega_n, (Polarization(nivn,ns),ns=1,NumSpin)
        end do
!
        deallocate( Polarization, Pol, stat=error)
        if(error .ne. 0) stop 'Main: error deallocating Polarization'
      close(100)
!
  end if
!
!*****************
!** free memory **
!*****************
!
  deallocate(SelfEn, stat=error)
    if(error .ne. 0) stop 'Main: error deallocating SelfEn'
!
  deallocate(BareVertex, stat=error)
  if(error .ne. 0) stop 'Main: error deallocating BareVertex'
!
  deallocate(Hamk, stat=error)
  if(error .ne. 0) stop 'Main: error deallocating Hamk'
!
  call KpointsDestroy()
!
  call InputDestroy()
!
  call barrier()
!
  call walltime (endtime) 
  if (myid .eq. 0) then
    print *, 'time is ', (endtime-starttime)/3600.0, ' hours'
  endif
!
  call CloseParallel
!
 end program Transport

