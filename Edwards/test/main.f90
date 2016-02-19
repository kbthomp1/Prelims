program main

  use kinddefs, only : dp
  use tools,    only : get_s, get_s_test, get_kc

  implicit none

  real(dp), dimension(6) :: kf, kc, kb
  real(dp), dimension(7) :: s, c, m_wt

  integer :: i

  real(dp) :: sum_of_s

  real(dp), parameter :: zero   = 0.0_dp
  real(dp), parameter ::  one   = 1.0_dp
  real(dp), parameter ::  two   = 2.0_dp
  real(dp), parameter ::  three = 3.0_dp

  logical :: debug = .false.
continue

! Forward rate constants
  kf = [ 1.9e19_dp, 3.e21_dp, 6.e18_dp, 1.e7_dp, 2.4e24_dp, 1.e5_dp ]

! Concentations  
  do i = 1,7
    c(i) = real(i,dp)
  end do

! Molecular weights
  m_wt(1:7) = [ two, two, three, one, one, two, three ]

! Equilibrium constants
  kc(1:6) = get_kc(c)

! Backward rate constants
  kb(1:6) = kf(1:6)/kc(1:6)

  if (debug) then
    s(1:2) = get_s_test(kf,kb,c)
    s(1) = two*s(1)
    s(2) = s(2)
    sum_of_s = sum(s(1:2))
  else
    s = get_s(kf,kb,c)
    do i = 1, 7
      s = m_wt(i)*s(i)
    end do
    sum_of_s = sum(s(1:7))
  end if

  write(*,*) "CHECK: sum = ",sum_of_s
  if ( sum_of_s /= zero) then
    write(*,*) "FAILED: sum of the source term vector is /= 0"
    stop 1
  else
    write(*,*) "PASSED"
  end if

end program main
