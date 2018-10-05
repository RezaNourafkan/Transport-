!
 module OpticMod
  use PrecisionMod
  implicit none
!
  contains
!------------------------------------------
!      Adding the contribution from |niwn| > NiwnMax-Nivn. We calculate this part analytically assuming
!      NiwnMax-Nivn is large enough that the Greens function for higher frequency can be described by
!      asymptotic behaviour GFn(iw_m) ~ (iw_m)^{-1}. We ignore the off-diagonal terms.
!      Therefore the tail contrbution is given by
!      -T\sum_{m=-\infty}^{\infty} (iw_m)^{-1} (iw_m+iv_n)^{-1} 
!      + T\sum_{m=-(NiwnMax+1-Nivn)}^{NiwnMax-Nivn} (iw_m)^{-1} (iw_m+iv_n)^{-1}
!      which is equal to
!      (1/4T) \delta_{n,0} + &
!       T\sum_{m=-(NiwnMax+1-Nivn)}^{NiwnMax-Nivn} (iw_m)^{-1} (iw_m+iv_n)^{-1}
!
!      where we used
!      -T\sum_{m=-\infty}^{\infty} (iw_m)^{-1} (iw_m+iv_n)^{-1}  = 0 for v_n \neq 0
!      -T\sum_{m=-\infty}^{\infty} (iw_m)^{-1} (iw_m)^{-1} = -f'(0)
!      -f'(0) = f(0)(1-f(0))/T = 1/4T
!      f(0) = 1/2 
!
  subroutine OpticalPolarization(Nchoice, Numspin, NumOrbitals, NumIrkpoints, NiwnMax, startkp, endkp, InverseTemperature, &
                                 ChemicalPotential, vk, Hamk, SelfEn, NivnMin, NivnMax, pol)
   use MathConstantsMod
   use LibMod
   implicit none
!
   integer(i4b), intent(in) :: Nchoice, Numspin, NumOrbitals, NumIrkpoints, NiwnMax, startkp, endkp, NivnMin, NivnMax
   real(dp), intent(in) :: InverseTemperature, ChemicalPotential
   complex(dpc), dimension(3,NumOrbitals, NumOrbitals, NumIrKpoints, NumSpin), intent(in) :: vk
!
   complex(dpc), dimension(NumOrbitals, NumOrbitals, NumIrkpoints, NumSpin), intent(in) :: Hamk   
   complex(dpc), dimension(NumOrbitals, NumOrbitals, 0:NiwnMax, NumSpin), intent(in) :: SelfEn
   complex(dpc), dimension(NivnMin:NivnMax, endkp-startkp+1,NumSpin), intent(inout) :: pol
!
   integer(i4b) :: Norb1, error, ns, kindex, niwn, Nivn, i, nk
   integer(i4b), dimension(2) :: dir 
!
   real(dp) :: Omega_1, Omega_2, Omega_3, Omega_n
   complex(dpc), dimension(:,:,:,:), allocatable :: GFnvertex
   complex(dpc), dimension(:,:), allocatable :: mat, Identity, GFn, tail, dpol
   complex(dpc) :: para
!
   select case(Nchoice)
    case(1)   ! xx conductivity
      dir(1) = 1
      dir(2) = 1
    case(2)   ! yy conductivity
      dir(1) = 2
      dir(2) = 2
    case(3)   ! zz conductivity
      dir(1) = 3
      dir(2) = 3
    case(4)   ! xy conductivity
      dir(1) = 1
      dir(2) = 2
    case(5)   ! xz conductivity
      dir(1) = 1
      dir(2) = 3
    case(6)   ! yz conductivity
      dir(1) = 2
      dir(2) = 3
    case default
      write(*,*) 'NChoice is not valid.'
      stop
   end select
!
   allocate(Identity(NumOrbitals,NumOrbitals), GFn(NumOrbitals,NumOrbitals), mat(NumOrbitals,NumOrbitals), &
            tail(NumOrbitals,NumOrbitals), dpol(NumOrbitals,NumOrbitals), stat= error)
   if (error .ne. 0) stop 'OpticalPolarization: unable to allocate mat'
   allocate(GFnVertex(2,NumOrbitals,NumOrbitals,-(NiwnMax+1):NiwnMax), stat= error)
   if (error .ne. 0) stop 'OpticalPolarization: unable to allocate GFnVertex'
!
   Identity = (zero,zero)
   do Norb1=1, NumOrbitals
     Identity(Norb1,Norb1) = (one,zero)
   end do
!
   Pol = (zero, zero)
   para = (zero, zero)
!
   do ns = 1, NumSpin
     nk = 0
     do kindex = startkp, endkp
       nk = nk+1
       do Niwn = -(NiwnMax+1), NiwnMax
         if (Niwn .lt. 0) then
           Omega_n = (two*Niwn+one)*pi_d/InverseTemperature
           GFn = (xi*Omega_n+ ChemicalPotential)*Identity-Hamk(:,:,kindex,ns)-transpose(dconjg(SelfEn(:,:,-niwn-1,ns)))
           call InverseUsingLapack(GFn)           
         else
           Omega_n = (two*Niwn+one)*pi_d/InverseTemperature
           GFn = (xi*Omega_n+ ChemicalPotential)*Identity-Hamk(:,:,kindex,ns)-SelfEn(:,:,niwn,ns)
           call InverseUsingLapack(GFn)
         end if
         GFnVertex(1,:,:,Niwn) = matmul( GFn(:,:), vk(dir(1),:,:,kindex,ns) )
         GFnVertex(2,:,:,Niwn) = matmul( GFn(:,:), vk(dir(2),:,:,kindex,ns) )
       end do
!
!      Paramegnetic term       
       mat = matmul( vk(dir(2),:,:,kindex,ns),vk(dir(1),:,:,kindex,ns) )
       do nivn = NivnMin, NivnMax
!
!        Summation over frequencies
         dpol = (zero,zero)
         tail = (zero,zero)
         do Niwn = -(NiwnMax+1-abs(Nivn)), NiwnMax-abs(Nivn)
           Omega_1 = (two*Niwn+one)*pi_D/InverseTemperature
           Omega_2 = (two*(Niwn+Nivn)+one)*pi_D/InverseTemperature
           Omega_3 = (two*(Niwn-Nivn)+one)*pi_D/InverseTemperature         
!           tail = tail - (one/(xi*Omega_1))*(one/(xi*Omega_2))*mat
           tail = tail - ((one/(xi*Omega_1))*(one/(xi*Omega_2))+(one/(xi*Omega_1))*(one/(xi*Omega_3)))*mat/Two
!           dpol = dpol - matmul( GFnVertex(2,:,:,Niwn+Nivn), GFnVertex(1,:,:,Niwn) )
           dpol = dpol - (matmul( GFnVertex(2,:,:,Niwn+Nivn), GFnVertex(1,:,:,Niwn) )+matmul( GFnVertex(2,:,:,Niwn), GFnVertex(1,:,:,Niwn-Nivn) ))/Two 
         end do  ! frequency summation
         if(nivn .eq. 0) then
           tail = tail - (InverseTemperature**2)*mat/4.0_dp
         end if
!        Trace over orbital degree of freedom 
         do i = 1, NumOrbitals
           Pol(nivn,Nk,ns)  =  Pol(nivn,Nk,ns) + (dpol(i,i) - tail(i,i)) / InverseTemperature
         end do
!
!         if (nivn .eq. 1) write(*,'(i4,2x,i6,4x,8f18.12)') nivn,nk,dpol / InverseTemperature
       end do  ! Bosonic Matsubara freq.
       para = para + Pol(0,Nk,ns)
!
     end do    ! Kpoints  
   end do      ! nspin
!
!  Free memory
   deallocate(GFnVertex, Identity, GFn, mat, tail, dpol, stat= error)
   if (error .ne. 0) stop 'OpticalPolarization: unable to deallocate GFnVertex'
!
  end subroutine OpticalPolarization
!---------------------------------------------------
 end module OpticMod
