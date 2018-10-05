!------------------------------------------------
 module PrecisionMod
  implicit none
! Symbolic names for kind types of 4-, 2-, and 1-byte integers:
  integer, parameter :: I4b = selected_int_kind(9)
  integer, parameter :: I2b = selected_int_kind(4)
  integer, parameter :: I1b = selected_int_kind(2)
! Symbolic names for kind types of single- and double-precision reals:
  integer, parameter :: sp = kind(1.0)
  integer, parameter :: dp = kind(1.0D0)
! Symbolic names for kind types of single- and double-precision complex:
  integer, parameter :: spc = kind((1.0,1.0))
  integer, parameter :: dpc = kind((1.0D0,1.0D0))
! Symbolic name for kind type of default logical:
  integer, parameter :: lgt = kind(.true.) 
 end module PrecisionMod 
!------------------------------------------------
