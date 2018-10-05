!
 module CurrentVertexMod
  use PrecisionMod
  implicit none
!
!
  contains
!------------------------------------------------
  subroutine BareCurrentVertexInit(casename, Numspin, Numorbitals, numkpoints, BareVertex)
   use MathConstantsMod,  only : xi
   implicit none
!
   character(len=100) :: casename
   integer(i4b), intent(in) :: Numspin, Numorbitals, Numkpoints
   complex(dpc), dimension(3, Numorbitals, NumOrbitals, NumKpoints, NumSpin) :: BareVertex
!
   integer(i4b), parameter :: NUnitIn = 110
   integer(i4b) :: n, m, Nkp, Nk, ns
   integer(i4b) :: norb1, norb2, i, error
   real(dp), dimension(6) :: val
   character(len=4) :: txt
   character(len=100) :: filename
   real(dp) :: factor
!
   do ns = 1, NumSpin
     if (NumSpin .eq. 1) then
       filename = './Input/'//trim(casename)//'.vk'
     else
       if (ns .eq. 1) then
         filename = './Input/'//trim(casename)//'.vk'//'up'
       else
         filename = './Input/'//trim(casename)//'.vk'//'dn'
       end if
     end if
!
!     open(NUnitIn, file=filename, form='formatted', status='old', position='REWIND', action='READ')
!     read(NUnitIn,*)
     open(NUnitIn, file=filename, form='unformatted', status='old', position='REWIND', action='READ')
!
     do nk = 1, Numkpoints 
!       read(NUnitIn,'(/,2x,a4,i6)') txt, nkp
       read(NUnitIn) nkp
       do n=1,NumOrbitals
         do m=1,NumOrbitals
!           read(NUnitIn,'(3x,2i4,6e13.6)') Norb1,Norb2, val
           read(NUnitIn) Norb1,Norb2, val
           BareVertex(1,Norb1, Norb2, nkp,ns) = val(1)+xi*val(2)
           BareVertex(2,Norb1, Norb2, nkp,ns) = val(3)+xi*val(4)
           BareVertex(3,Norb1, Norb2, nkp,ns) = val(5)+xi*val(6) 
!           write(*,'(i8, 2f16.10,2f16.10,2f16.10)') nk, BareVertex(1,Norb1, Norb2, nkp,ns), BareVertex(2,Norb1, Norb2, nkp,ns), BareVertex(3,Norb1, Norb2, nkp,ns)
         end do
       end do
     end do
!
     close(NUnitIn)
   end do
!
  end subroutine BareCurrentVertexInit
!
!------------------------------------------------
 end module CurrentVertexMod
