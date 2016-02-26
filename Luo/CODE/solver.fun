test_suite solver

  integer, parameter :: dp = selected_real_kind(P=15)

test check_lhspo

  use test_helper

  real(dp), dimension(:),   allocatable :: rhspo
  real(dp), dimension(:,:), allocatable :: lhspo

  integer :: ndof, i, j

  real(dp) :: uinf, vinf

continue

  call setup_cube
  ndof = grid%npoin + grid%numfac
  uinf = one; vinf = zero
  allocate(lhspo(ndof,ndof))
  allocate(rhspo(ndof))

  lhspo = get_lhspo(grid,ndof)
  rhspo = get_rhspo(grid,ndof,uinf,vinf)
  !lhspo = zero
  !call add_flux_contributions(lhspo,grid,ndof)

  !call add_lift_primal_domain(1,2,5,grid,lhspo)
  !call add_lift_second_domain(1,2,5,grid,lhspo)
  !call add_jump_second_face(1,2,5,grid,lhspo)
  !call add_flux_primal_face(1,2,5,grid,lhspo)

  do i = 1, ndof
    write(*,20) "CHECK: lhspo = ",(lhspo(i,j), j=1,ndof)
  end do
  do i = 1, ndof
    write(*,20) "CHECK: rhspo = ",rhspo(i)
  end do

20 format(A,300g11.4)

end test

end test_suite
