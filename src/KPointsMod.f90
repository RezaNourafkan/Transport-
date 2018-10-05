!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This is a module which reads the k points from 
! the output of the wien2k package
!
!
! Let us define \mathcal{S}, the set of real 3x3 matrices describing the symmetry operations of the crystal 
! with respect to a real Cartesian basis. Furthernore, {k_i^ir} denotes the irreduciable k-mesh, i.e.
! the reduced mesh when the symmetry operations of \mathcal{S} are exploited, with the corresponding mapping
! m_s(n) = n_s such than n labels a k-point in full mesh and n_s label a k-point in irreduciable k-mesh. 
! If one replaces all vertices m of the full tetrahedral \mathcal{T} by their reduced vertices m_s(m), one formally
! obtains the new tetrahedral mesh \mathcal{T'}_s. Note that \mathcal{T'}_s might include 
! elements that do not correspond to real tetrahedra
! but have e.g. equal nodes when multiple vertices of a element of T \in \mathcal{T} 
! have been mapped onto same k-point in the reduced
! set. After mapping \mathcal{T} --> \mathcal{T}'_s there are in general multiple occurrences of an element. 
! Then the irreducible tetrahedral \mathcal{T}_s
! denotes the reduced symmetrized mesh, i.e. all elements T_s \in \mathcal{T}_s are only considered once 
! and the weigth w_{T_s} accounts for the volume and the multiplicity of T_S.
! 
! For scalar integrand      g = \sum_{T \in \mathcal{T}} g^T = \sum_{T_s \in \mathcal{T}_s} w_{T_s} g^{T_s} 
! for tensorial integrand   g = \sum_{T \in \mathcal{T}} g^T = 
!                               (1/|\mathcal{S}|)\sum_{s \in \mathcal{S}} \sum_{T_s \in \mathcal{T}_s} w_{T_s}  S^+ g^{T_s} S
! Therefore the integrand is only computed for irreducible mesh, but summation is on full k-mesh. Note that in this case
! (1/|\mathcal{S}|)\sum_{s \in \mathcal{S}} w_{T_s} only accounts  for the volume of tetrahedra. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 module KpointsMod
  use PrecisionMod
  implicit none
!
  real(dp)                                    :: UnitCellVolume
  real(dp), dimension(3,3)                    :: ReciprocalVectors
  integer(i4b)                                :: NumIrKpoints
  integer(i4b)                                :: NumIrTetrahedra  
  real(dp)                                    :: TetrahedraVolume               ! All tetrahedra have same volume
  integer(i4b), dimension(:,:), allocatable   :: IrTetrahedra                   ! Different tetrahedra in Irreducible part of BZ
  integer(i4b), dimension(:), allocatable     :: IrTetrahedraMultiplicity
!
  integer(i4b)                            :: NumSymmetryOperations 
  real(dp), dimension(:,:,:), allocatable :: SymmetryOperations
!
  integer(i4b), dimension(:), allocatable :: StartKP, EndKP
!
  contains
!------------------------------------------------
   subroutine KpointsInit(casename)
    use MathConstantsMod
    use ParallelMod, only : NPE
    implicit none
!
    character(len=100) :: casename
!
    character(len=100) :: filename
    integer(i4b), parameter :: NUnitIn = 1111
    logical(lgt) :: i_exist
    integer(i4b) :: n, i, j, error, step
!
!   read-in information about reciprocal basis
    filename = './Input/'//trim(casename)//'.kpoints'
    inquire(file=trim(filename), exist=i_exist)
    if(.not. i_exist) then
      stop 'KpointsReader: file case.kpoints does not exsit'
    end if
    open(NUnitIn, file=trim(filename), position='rewind', action='read')
!
    write(*,'(/,a40)')'Reciprocal unit vectors in 2*pi*bohr^-1:'
    do j = 1, 3
       read(NUnitIn,*) (ReciprocalVectors(j,i), i = 1, 3)
       write(*,'(3f10.6)') ReciprocalVectors(j,1:3)
    enddo
!
!   compute unit cell volume in bohr^3
    UnitCellVolume = ReciprocalVectors(1,1)*(ReciprocalVectors(2,2)*ReciprocalVectors(3,3)-ReciprocalVectors(3,2)*ReciprocalVectors(2,3)) &
                   + ReciprocalVectors(1,2)*(ReciprocalVectors(2,3)*ReciprocalVectors(3,1)-ReciprocalVectors(3,3)*ReciprocalVectors(2,1)) &
                   + ReciprocalVectors(1,3)*(ReciprocalVectors(2,1)*ReciprocalVectors(3,2)-ReciprocalVectors(3,1)*ReciprocalVectors(2,2))
    UnitCellVolume = one/UnitCellVolume*(two*pi_d)**3
    write(*,'(/,a27,3x,f10.6)')'Unit cell volume in bohr^3:', UnitCellVolume
!
   read(NUnitIn,*)  numIrKPoints
   write(*,'(/,a20,4x,i10)')'Number of Irkpoints:', numIrKPoints
!
   read(NUnitIn, *) NumIrtetrahedra
   if (NumIrtetrahedra .ne. 0) then
     read(NUnitIn, *) TetrahedraVolume
     write(*,'(a24,x,i10)')' Number of Irtetrahedra:', NumIrtetrahedra
     write(*,'(a30,x,f10.6)')' Normalized Tetrahedra volume:', TetrahedraVolume   
     read(NUnitIn, *) 
     allocate(IrTetrahedra(4,NumIrtetrahedra), IrTetrahedraMultiplicity(NumIrtetrahedra), stat=error)
     if(error .ne. 0) stop 'KpointsInit: unable to allocate TetrahedraMultiplicity'
     do j = 1, NumIrtetrahedra
       read(NUnitIn, *) IrTetrahedraMultiplicity(j), (IrTetrahedra(i,j), i=1,4)
     end do 
!
     read(NUnitIn, *)
     read(NUnitIn, *) NumSymmetryOperations
     allocate(SymmetryOperations(NumSymmetryOperations,3,3), stat=error)
     if (error .ne. 0) stop 'KpointsInit: unable to allocate SymmetryOperations'
     do i = 1, NumSymmetryOperations
       do j = 1, 3
         read(NUnitIn, *) SymmetryOperations(i,j,1:3)
       end do
       read(NUnitIn, *)
     end do
   end if  
!
   close(NUnitIn)
!
!
!  find startkp and endkp for parallel camputation
   allocate(startKP(0:NPE-1), endKP(0:NPE-1))
   step = numIrKPoints/NPE
   endKP(NPE-1) = numIrKPoints
   startKP(NPE-1) = EndKP(NPE-1) - step + 1
   do n = NPE-2, 0, -1
     endKP(n) = startKP(n+1) - 1
     startKP(n) = EndKP(n) - step + 1
     if(n .eq. 0) StartKP(n) = 1
   end do
!
   write(*,'(/,a17, 2x, i8)') 'Number of kpoints', numIrKPoints
   do n=0, NPE-1
     write(*,'(a29)') 'Process, startkp, endkp, load'
     write(*,'(i3, 2x, 3i8)') n, startKP(n), endKP(n), endKP(n)-startKP(n)+1
   end do
!
  end subroutine KpointsInit
!
!------------------------------------------------ 
   subroutine kPointsDestroy()
    implicit none
! 
    integer(i4b) :: error
    if(allocated(IrTetrahedra)) then
      deallocate(IrTetrahedra, IrTetrahedraMultiplicity, stat=error)
      if(error .ne. 0) stop 'KpointsDestroy: unable to deallocate IrTetrahedraMultiplicity'
      deallocate(SymmetryOperations, stat=error)
      if(error .ne. 0) stop 'KpointsDestroy: unable to deallocate SymmetryOperations'
    end if
!
    deallocate(startKP, endKP, stat=error)
    if(error .ne. 0) stop 'KpointsDestroy: unable to deallocate startKP'
!    
   end subroutine kPointsDestroy   
!------------------------------------------------
 end module KpointsMod    
