module tools

  use kinddefs, only : dp

  implicit none
  private

  public :: get_s, get_s_test
  public :: get_kc

  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter ::  one = 1.0_dp
  real(dp), parameter ::  two = 2.0_dp

contains

!=============================== NOMENCLATURE ================================80
! s is the source term vector    -> \dot{s}
! c is the species concentration -> C_i
!
! The array number is formed as
!
! s/c(1) => H_2  (gas)
! s/c(2) => O_2  (gas)
! s/c(3) => H_2O (gas)
! s/c(4) => H    (surface)
! s/c(5) => O    (surface)
! s/c(6) => OH   (surface)
! s/c(7) => H_2O (surface)
!=============================================================================80

  function get_s(kf,kb,c) result (s)

    real(dp), dimension(7) :: s
    real(dp), dimension(:), intent(in) :: kf,kb,c

  continue

    s(1) = (-one)*(kf(1)*c(1)      - kb(1)*c(4)**2)

    s(2) = (-one)*(kf(2)*c(2)      - kb(2)*c(5)**2)

    s(3) = ( one)*(kf(6)*c(7)      - kb(6)*c(3))

    s(4) = ( two)*(kf(1)*c(1)      - kb(1)*c(4)**2)    &
         + (-one)*(kf(3)*c(4)*c(5) - kb(3)*c(6))       &
         + (-one)*(kf(4)*c(4)*c(6) - kb(4)*c(7))

    s(5) = ( two)*(kf(2)*c(2)      - kb(2)*c(5)**2)    &
         + (-one)*(kf(3)*c(4)*c(5) - kb(3)*c(6))       &
         + (-one)*(kf(5)*c(6)**2   - kb(5)*c(7)*c(5))

    s(6) = ( two)*(kf(3)*c(6)      - kb(3)*c(4)*c(5))  &
         + (-one)*(kf(4)*c(4)*c(6) - kb(4)*c(7))       &
         + (-two)*(kf(5)*c(6)**2   - kb(5)*c(7)*c(5))

    s(7) = ( one)*(kf(4)*c(4)*c(6) - kb(4)*c(7))       &
         + ( one)*(kf(5)*c(6)**2   - kb(5)*c(7)*c(5))  &
         + (-one)*(kf(6)*c(7)      - kb(6)*c(3))

  end function get_s

  function get_s_test(kf,kb,c) result (s)

    real(dp), dimension(2) :: s
    real(dp), dimension(:), intent(in) :: kf,kb,c

  continue

    s(1) = (-one)*(kf(1)*c(1) - kb(1)*c(4)**2)
    s(2) = ( two)*(kf(1)*c(1) - kb(1)*c(4)**2)

  end function get_s_test

  function get_kc(c) result(kc)

    real(dp), dimension(6) :: kc
    real(dp), dimension(:), intent(in) :: c

  continue

    kc(1) = c(4)**2/c(1)
    kc(2) = c(5)**2/c(2)
    kc(3) = c(6)/(c(4)*c(5))
    kc(4) = c(7)/(c(4)*c(6)) 
    kc(5) = c(7)*c(5)/c(6)**2
    kc(6) = c(3)/c(7)

  end function get_kc

end module tools
