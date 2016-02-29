test_suite flux_functions

  integer, parameter :: dp = selected_real_kind(P=15)

test check_mass_matrix
  use test_helper
  real(dp), dimension(3,3) :: m_inv
  integer :: i,j,ndof
continue
  ndof = 3
  m_inv = get_inv_mass_matrix(two)
  do i = 1, ndof
    write(*,20) "CHECK: M_inv = ",(m_inv(i,j), j=1,ndof)
  end do
20 format(A,300g11.4)
end test

test check_residual

  use test_helper
  use solver, only : add_flux_contributions

  real(dp), dimension(:),   allocatable :: residual, phi

  integer :: ndof, i, j

  real(dp) :: uinf, vinf

continue

  call setup_cube
  ndof = grid%npoin
  uinf = one; vinf = zero
  allocate(residual(ndof))
  allocate(phi(ndof))
  phi = zero
  residual = zero

  !lhspo = get_lhspo(grid,ndof)
  !rhspo = get_rhspo(grid,ndof,uinf,vinf)
  !call add_flux_contributions(lhspo,grid,ndof)

  !call add_lift_primal_domain(1,2,5,grid,lhspo)
  !call add_lift_second_domain(1,2,5,grid,lhspo)
  !call add_jump_second_face(1,2,5,grid,lhspo)
  !call add_flux_primal_face(1,2,5,grid,phi,residual)
  call add_boundary_flux(residual,grid,uinf,vinf)

  do i = 1, ndof
    write(*,20) "CHECK: residual = ",residual(i)
  end do

20 format(A,300g11.4)

end test

end test_suite
