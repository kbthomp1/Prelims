program main

  use kinddefs, only : dp
  use tests,    only : test_s_sum, test_coef_conversion

  implicit none

continue

  call test_s_sum
  !call test_coef_conversion

end program main
