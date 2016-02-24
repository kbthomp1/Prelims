test_suite solver

  integer, parameter :: dp = selected_real_kind(P=15)

  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter :: one  = 1.0_dp

test check_lhspo

  use test_helper

  real(dp), dimension(:,:),     allocatable :: rhspo
  real(dp), dimension(:,:,:,:), allocatable :: lhspo

  integer :: nvars, neqns, i, j

continue

  nvars = 1
  neqns = 1

  call setup_cube
  allocate(lhspo(nvars,neqns,npoin,npoin))
  allocate(rhspo(nvars,npoin))

  lhspo = get_lhspo(npoin,nelem,nnode,nvars,neqns,inpoel,geoel)

  do i = 1, npoin
    write(*,*) "CHECK: lhspo = ",(lhspo(1,1,i,j), j=1,npoin)
  end do

end test

end test_suite
