!
 module ParallelMod
  use PrecisionMod
  implicit none
!
!       NPE       - Number of processing elements
!       MYID      - Number of current processor
  integer(i4b) :: NPE, MYID, comm
!
  contains
   subroutine InitParallel
    implicit none
    include 'mpif.h'
!
    integer(i4b) :: ierr
!
    call MPI_INIT( ierr )
    comm = MPI_COMM_WORLD
    call MPI_COMM_RANK( MPI_COMM_WORLD, MYID, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, NPE, ierr )
    call Barrier 
!
   end subroutine InitParallel
!
   subroutine Barrier
    implicit none
    include 'mpif.h'
    integer(i4b) :: ierr
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)    
   end subroutine barrier
!
   subroutine CloseParallel
    implicit none
    include 'mpif.h'
    integer(i4b) :: ierr
    CALL MPI_FINALIZE(IERR)
   end subroutine CloseParallel
!
   subroutine AbrotParallel
    implicit none
    include 'mpif.h'
    integer(i4b) :: ierr
    CALL MPI_ABORT(MPI_COMM_WORLD, IERR)
   end subroutine AbrotParallel
!
   subroutine Walltime(time)
    implicit none
    include 'mpif.h'
    real(dp) :: time
    time = MPI_WTIME()
   end subroutine walltime
!
 end module ParallelMod
