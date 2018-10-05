!
 module InputMod
  use PrecisionMod
  implicit none
!
  character(len=100) :: casename
  logical(lgt) :: InteractingSystem
  real(dp) :: ChemicalPotential
  integer(i4b) :: NumOrbitals, NumSpin
  real(dp) :: InverseTemperature
  integer(i4b) :: NiwnMax, NivnMin, NivnMax 
!
  integer(i4b) :: NChoice 
!
  integer(i4b), dimension(:,:), allocatable :: PackOrbIndex
  integer(i4b), dimension(:), allocatable :: UnPackOrbIndex1, UnPackOrbIndex2
!
  contains
!------------------------------------------------
!
   subroutine Inputinit()
    implicit none
    integer(i4b) :: n, m, error, Norb1, norb2, i
!
!   read some input parameters
    open(100, file='./Input/Input.dat', status='old')
    read(100,*)  casename
    write(*,'(/,a14,2x)',advance='no') 'Case name is: '
    write(*,*) trim(casename)
    read(100,*) InteractingSystem
    write(*,'(a20, l)') 'Interacting system: ', InteractingSystem
! 
    read(100,*)
    read(100,*) NumOrbitals
    read(100,*) Numspin
    write(*,'(a20,4x,i4)') 'Number of orbitals =', NumOrbitals
    write(*,'(a23,4x,i4)') 'Number of Spin Blocks =', Numspin
    read(100,*)
    read(100,*) chemicalpotential
    write(*,'(a21,4x,f16.8)') 'Chemical potential = ', chemicalpotential
    read(100,*) InverseTemperature
    write(*,'(a21,4x,f16.8)') 'Inverse Temperature =', InverseTemperature
    read(100,*)
    read(100,*) Niwnmax
    write(*,'(a9,4x,i6)') 'NiwnMax =', NiwnMax
    read(100,*) Nivnmax
    if(Nivnmax .gt. Niwnmax/2-1) Nivnmax = Niwnmax/2 - 1 
     NivnMin = 0
     write(*,'(a9,4x,i6)') 'NivnMin =', NivnMin
     write(*,'(a9,4x,i6)') 'NivnMax =', NivnMax
!
    read(100,*)
    read(100,*) Nchoice
    if (NChoice .lt. 1 .and. NChoice .gt. 3) then
      write(*,'(/,4x,a25)') 'The NChoice is undefined!'
      stop
    end if
    write(*,'(a11,4x,i4)') 'NChoice is:', NChoice
!
    close(100)
!
    allocate(PackOrbIndex(NumOrbitals,NumOrbitals), UnPackOrbIndex1(NumOrbitals**2), &
             UnPackOrbIndex2(NumOrbitals**2), stat = error)
    if (error .ne. 0) stop 'InputInit: unable to allocate PackOrbIndex'
    i = 0
    do Norb1 = 1, NumOrbitals
      do Norb2 = 1, NumOrbitals
        i = i + 1
        PackOrbIndex(NOrb1, NOrb2) = i
        UnPackOrbIndex1(i) = NOrb1
        UnPackOrbIndex2(i) = NOrb2
      end do
    end do
!
   end subroutine Inputinit
!------------------------------------------------
!
   subroutine InputDestroy()
    implicit none
    integer(i4b) :: error
!
    deallocate(PackOrbIndex, UnPackOrbIndex1, UnPackOrbIndex2, stat = error)
    if (error .ne. 0) stop 'InputDestroy: unable to deallocate PackOrbIndex'    
!
   end subroutine InputDestroy
!------------------------------------------------
!
  end module InputMod
