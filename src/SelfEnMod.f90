!
 module SelfEnMod
  use PrecisionMod
  implicit none
!
!
  contains
!------------------------------------------------
  subroutine SelfEnInit2(casename, Numspin, Numorbitals, Niwnmax, InverseTemperature, SelfEn)
   use MathConstantsMod,  only : xi, pi_d
   implicit none
!
   character(len=100) :: casename
   integer(i4b), intent(in) :: Numspin, Numorbitals, NiwnMax
   real(dp), intent(in) :: InverseTemperature
   complex(dpc), dimension(Numorbitals, NumOrbitals, 0:Niwnmax, NumSpin) :: SelfEn
!
   integer(i4b), parameter :: NUnitIn = 111
   integer(i4b) :: n, m, Ns, Nnd, Nw
   integer(i4b) :: norb1, norb2, i, error
   real(dp) :: beta, omega_n
   real(dp), dimension(2*NumOrbitals**2) :: val
   character(len=100) :: filename
!
   do ns = 1, NumSpin
     SelfEn(:,:,:,Ns) = (0.0_dp,0.0_dp)
!
     if (NumSpin .eq. 1) then
       filename = './Input/'//trim(casename)//'.SelfEnIA'
     else
       if (ns .eq. 1) then
         filename = './Input/'//trim(casename)//'.SelfEnIAup'
       else
         filename = './Input/'//trim(casename)//'.SelfEnIAdn'
       end if
     end if
!
     open(NUnitIn, file=filename, status='old', form='formatted', position='REWIND', action='READ')
!
     do n = 0, Niwnmax 
       read(NUnitIn,*) m, Omega_n, val
       i = 1
       do norb1=1,NumOrbitals
         do norb2=1,NumOrbitals
           SelfEn(Norb1, Norb2,m,Ns) = val(i)+xi*val(i+1) 
	   i = i+2
         end do
       end do
     end do
     close(NUnitIn)
!
     beta = (dble(2*Niwnmax)+1.0_dp)*pi_d/Omega_n
     if(abs(beta - InverseTemperature) .gt. 1.0D-4) stop 'ReadSelfEn:  Inverse temperatures do not match in two files'
!
   end do
!
  end subroutine SelfEnInit2
!
!------------------------------------------------
 end module SelfEnMod
