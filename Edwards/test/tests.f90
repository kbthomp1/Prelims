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

  real(dp), parameter ::  tolerance = 1.E-15_dp

  real(dp), dimension(6), parameter :: &
    kf_cgs = [ 1.9e19_dp, 3.e21_dp, 6.e18_dp, 1.e7_dp, 2.4e24_dp, 1.e5_dp ]

  real(dp), dimension(6), parameter :: &
    kf_mks = [ 1.9e17_dp, 3.e19_dp, 6.e17_dp, 1.e6_dp, 2.4e23_dp, 1.e5_dp ]

contains

!======================== TEST_COEF_CONVERSION ==============================80
! Test that kf and kb are correct in both mks and cgs unit basis
!============================================================================80
  subroutine test_coef_conversion

    real(dp), dimension(6) :: kf, kc, kb, convert
    real(dp), dimension(7) :: c, s_cgs, s_mks, diff

    integer :: i

  continue

!   Everything in cm, mol, s
    kf = kf_cgs
    do i = 1,7
      c(i) = real(i,dp) !mol/cm^(2,3)
    end do
    !kc(1:6) = get_kc(c)
    kc = one
    kb(1:6) = kf(1:6)/kc(1:6)
    s_cgs = get_s(kf,kb,c) !mol/(cm^2 s)

!   Everything in m, kmol, s
    c(1:3) = 1.e3_dp*c(1:3) !kmol/m^3
    c(4:7) = 10._dp*c(4:7)  !kmol/m^2
    kf = kf_mks
    !kc(1:6) = get_kc(c)
    kc = one
    kb(1:6) = kf(1:6)/kc(1:6)
    s_mks = get_s(kf,kb,c) !kmol/(m^2 s)

    convert = [ 1.e-2_dp, 1.e-2_dp, 1.e-1_dp, 1.e-1_dp, 1.e-1_dp, one ]

    write(*,*) "CHECK: s_cgs =", s_cgs
    write(*,*) "CHECK: s_mks =", s_mks

    diff = (s_cgs - 10._dp*s_mks)

    write(*,*) "CHECK: diff = ", diff

    do i = 1,6
      if (diff(i) > tolerance) then
        write(*,*) "FAILED: test_coef_conversion"
        write(*,*) "something is wrong in the conversion of cgs to mks."
        stop 1
      end if
    end do

    write(*,*) "PASSED: test_coef_conversion"

  end subroutine test_coef_conversion

!============================ TEST_S_SUM ====================================80
! Test that the sum of the source terms is zero (wdot, actually...)
!============================================================================80
  subroutine test_s_sum

    real(dp), dimension(6) :: kf, kc, kb
    real(dp), dimension(7) :: s, c, m_wt
  
    integer :: i
  
    real(dp) :: sum_of_s

    logical :: debug = .false.

    real(dp), parameter :: m_H = one
    real(dp), parameter :: m_O = 8._dp
    real(dp), parameter :: m_N = 14._dp

  continue

!   Forward rate constants
    kf = kf_cgs

!   Concentations  
    do i = 1,7
      c(i) = real(i,dp)
    end do

!   Molecular weights
    m_wt(1:7) = [ two*m_H, two*m_O, two*m_H + m_O, m_H, m_O, m_O + m_H, &
                  two*m_H + m_O ]

!   Equilibrium constants
    !kc(1:6) = get_kc(c)
    kc(1:6) = one

!   Backward rate constants
    kb(1:6) = kf(1:6)/kc(1:6)

    !debug = .true.
    if (debug) then
      s = zero
      s(1:2) = get_s_test(kf,kb,c)
      s(1) = two*s(1)
      s(2) = s(2)
      sum_of_s = sum(s(1:2))
    else
      s = get_s(kf,kb,c)
      do i = 1, 7
        s(i) = m_wt(i)*s(i)
      end do
      sum_of_s = sum(s(1:7))
    end if

20 format(A,20e13.4)

    write(*,20) "CHECK: M = ",m_wt
    write(*,20) "CHECK: s = ",s

    write(*,20) "CHECK: sum = ",sum_of_s
    if ( sum_of_s /= zero) then
      write(*,*) "FAILED: test_s_sum"
      write(*,*) "sum of the source term vector is /= 0"
      stop 1
    else
      write(*,*) "PASSED: test_s_sum"
    end if

  end subroutine test_s_sum

end module tests
