!
 module HamiltonianMod
  use PrecisionMod
  implicit none
!
!
  contains
!------------------------------------------------
! this subroutine reads the Hamiltonian from the file
  subroutine HkInit(casename, NumSpin, Numorbitals, numkpoints, Hk)
   use MathConstantsMod,  only : xi
   implicit none
!
   character(len=100) :: casename
   integer(i4b), intent(in) :: NumSpin, Numorbitals, Numkpoints
   complex(dpc), dimension(Numorbitals, Numorbitals, NumKpoints, NumSpin) :: Hk
!
   integer(i4b), parameter :: NUnitIn = 111
   integer(i4b) :: n, m, nk, Nkp, ns
   integer(i4b) :: norb1, norb2, i, error
   real(dp), dimension(:), allocatable  :: val
   character(len=4) :: txt
   character(len=100) :: filename
!
   do ns = 1, NumSpin
     if (NumSpin .eq. 1) then
       filename = './Input/'//trim(casename)//'.Hk'
     else 
       if (ns .eq. 1) then
         filename = './Input/'//trim(casename)//'.Hkup'
       else
         filename = './Input/'//trim(casename)//'.Hkdn'
       end if 
     end if	 	 
!
!     open(NUnitIn, file=trim(filename), form='formatted', status='old', position='REWIND', action='READ')
!     read(NUnitIn,*)
     open(NUnitIn, file=trim(filename), form='unformatted', status='old', position='REWIND', action='READ')
     read(NUnitIn) N
     if(N .ne. Numorbitals) then
       write(*,*), N, NumOrbitals
       stop 'HkInit: inconsistent number of orbital'
     end if
!
     read(NUnitIn) Nk
     if(Nk .ne. Numkpoints) then
       write(*,*) Nk, Numkpoints
       stop 'HkInit:  inconsistent number of k point'
     end if
!
     allocate(val(2*NumOrbitals**2), stat=error)
     if(error .ne. 0) stop 'HkInit: unable to allocate val'
!
     do nk = 1, Numkpoints 
       read(NUnitIn)  val
       i = 1
       do norb1=1,NumOrbitals
         do norb2=1,NumOrbitals
           Hk(Norb1,Norb2, nk, ns) = val(i)+xi*val(i+1)
           i = i+2
         end do
       end do
!!       read(NUnitIn,'(/,2x,a4,i6)') txt, nkp
!       read(NUnitIn) nkp  
!       do n=1,NumOrbitals
!         do m=1,NumOrbitals
!!           read(NUnitIn,'(3x,2i4,2e13.6)') Norb1, Norb2, val
!           read(NUnitIn) Norb1, Norb2, val 
!           Hk(Norb1, Norb2, nkp, ns) = val(1)+xi*val(2) 
!         end do
!       end do
     end do
     deallocate(val, stat=error)
     if(error .ne. 0) stop 'HkReader: unable to deallocate val'
!
     close(NUnitIn)
   end do  
!
  end subroutine HkInit
!
!------------------------------------------------
 end module HamiltonianMod
