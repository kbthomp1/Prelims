test_suite solver

  integer, parameter :: dp = selected_real_kind(P=15)

test check_lhspo

  use test_helper

  !real(dp), dimension(:),   allocatable :: rhspo, phi, exact
  real(dp), dimension(:,:), allocatable :: lhspo

  integer :: ndof, i, j

continue

  call setup_cube
  ndof = 2*grid%npoin
  allocate(lhspo(ndof,ndof))

  lhspo = get_lhspo(grid,ndof)

  do i = 1, ndof
    write(*,20) "CHECK: lhspo = ",(lhspo(i,j), j=1,ndof)
  end do

20 format(A,300g11.4)

end test

end test_suite
