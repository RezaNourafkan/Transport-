!------------------------------------------------
 module MathConstantsMod
  use PrecisionMod
  implicit none
! Frequently used mathematical constants (with precision to spare):
  real(dp), parameter     :: zero  = 0.0_dp
  real(dp), parameter     :: one   = 1.0_dp
  real(dp), parameter     :: two   = 2.0_dp
  complex(dpc), parameter ::  xi   = (0.0_dp,1.0_dp)
  real(sp), parameter     :: sqrt2 = 1.41421356237309504880168872420969807856967_sp

  real(sp), parameter :: pi      = 3.141592653589793238462643383279502884197_sp
  real(sp), parameter :: pio2    = 1.57079632679489661923132169163975144209858_sp
  real(sp), parameter :: twopi   = 6.283185307179586476925286766559005768394_sp   
  real(dp), parameter :: pi_d    = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: pio2_d  = 1.57079632679489661923132169163975144209858_dp
  real(dp), parameter :: twopi_d = 6.283185307179586476925286766559005768394_dp
 
 end module MathConstantsMod 
!------------------------------------------------
!
 module PhysConstantsMod
  use PrecisionMod
  implicit none
  real(dp), parameter :: BohrR = 5.29177249e-11
  real(dp), parameter :: e2Overhbar = 0.000243413479154821    !(1/Ohm)
 end module PhysConstantsMod
