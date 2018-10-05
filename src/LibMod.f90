 module LibMod
  use PrecisionMod
  implicit none
!
  contains
!
!----------------------------------------
  subroutine InverseUsingLapack(A)
   use MathConstantsMod, only : zero, one
   implicit none
   complex(dpc), dimension(:,:) :: A
!
   integer(i4b) :: n, lda, lwork, info, i0, j0
   Complex(dpc) :: det,tmp 
   integer(i4b), dimension(:), pointer :: ipiv
   complex(dpc), dimension(:), pointer :: work
   character(len=1) :: uplo = 'U'
!
   n  = size(A, dim=1)
   i0 = lbound(A, dim=1)
   j0 = lbound(A, dim=2)  
   if (n .eq. 1) then
     A(i0,j0)= One/A(i0,j0)
!   else if (n .eq. 2) then 
!     det = A(i0,j0) * A(i0+1, j0+1) - A(i0+1,j0) * A(i0, j0+1) 
!     if (abs(det) < 1.0E-9) then
!       write(*,*) A 
!       stop 'InverseUsingLapack: det(A) is zero'
!     end if 
!     tmp = A(i0,j0)
!     A(i0,j0) = A(i0+1, j0+1)
!     A(i0+1, j0+1) = tmp
!     A(i0+1,j0) = -One*A(i0+1,j0)
!     A(i0, j0+1)= -One*A(i0, j0+1)
!     A = A/det 
   else
     lda = n
     lwork = 2*n
     allocate(ipiv(n), work(lwork))     
     call zgetrf( n, n, a, lda, ipiv, info )
      if (info .ne. 0) then  
       write(*,*) "******* ERROR *******"  
       write(*,*) "from ZGETRF library info = ",info 
       stop  
      end if 

     call zgetri( n, a, lda, ipiv, work, lwork, info )
!!    call zhetrf(uplo, n, a, lda, ipiv, work, lwork, info)
!!     call zhetri(uplo, n, a, lda, ipiv, work, info)
     if (info .ne. 0) then  
       write(*,*) "******* ERROR *******"  
       write(*,*) "from ZGETRI library info = ",info 
       stop  
     end if 
     deallocate(ipiv, work)
   end if
!
   
  end subroutine InverseUsingLapack
!
!----------------------------------------
 end module LibMod
