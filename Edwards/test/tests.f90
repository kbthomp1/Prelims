module tests

  use kinddefs, only : dp
  use tools,    only : get_s, get_s_test, get_kc

  implicit none
  private

  public :: test_s_sum
  public :: test_coef_conversion

  real(dp), parameter :: zero   = 0.0_dp
  real(dp), parameter ::  one   = 1.0_dp
  real(dp), parameter ::  two   = 2.0_dp
  real(dp), parameter ::  three = 3.0_dp

  real(dp), dimension(6), parameter :: &
    kf_cgs = [ 1.9e19_dp, 3.e21_dp, 6.e18_dp, 1.e7_dp, 2.4e24_dp, 1.e5_dp ]

  real(dp), dimension(6), parameter :: &
    kf_mks = [ 1.9e17_dp, 3.e19_dp, 6.e17_dp, 1.e6_dp, 2.4e23_dp, 1.e5_dp ]

contains

  subroutine test_coef_conversion

    real(dp), dimension(6) :: kf, kc, kb
    real(dp), dimension(7) :: c, s_cgs, s_mks

    integer :: i

  continue

!   Everything in cm, mol, s
    kf = kf_cgs
    do i = 1,7
      c(i) = real(i,dp) !mol/cm^(2,3)
    end do
    kc(1:6) = get_kc(c)
    kb(1:6) = kf(1:6)/kc(1:6)
    s_cgs = get_s(kf,kb,c) !mol/(cm^2 s)

!   Everything in m, kmol, s
    c(1:3) = 1.e3_dp*c(1:3) !kmol/m^3
    c(4:7) = 10._dp*c(4:7)  !kmol/m^2
    kf = kf_mks
    kc(1:6) = get_kc(c)
    kb(1:6) = kf(1:6)/kc(1:6)
    s_mks = get_s(kf,kb,c) !kmol/(m^2 s)

    write(*,*) "CHECK: diff = ", s_cgs*10._dp - s_mks

  end subroutine test_coef_conversion

  subroutine test_s_sum

    real(dp), dimension(6) :: kf, kc, kb
    real(dp), dimension(7) :: s, c, m_wt
  
    integer :: i
  
    real(dp) :: sum_of_s

    logical :: debug = .false.

  continue

!   Forward rate constants
    kf = kf_cgs

!   Concentations  
    do i = 1,7
      c(i) = real(i,dp)
    end do

!   Molecular weights
    m_wt(1:7) = [ two, two, three, one, one, two, three ]

!   Equilibrium constants
    kc(1:6) = get_kc(c)

!   Backward rate constants
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

  end subroutine test_s_sum

end module tests
