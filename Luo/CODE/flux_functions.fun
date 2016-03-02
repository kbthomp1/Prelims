test_suite flux_functions

  integer, parameter :: dp = selected_real_kind(P=15)

setup
  use test_helper
  logical, save :: init = .true.
  if(init) then
    call setup_cube
    ndof = grid%nnode*grid%nelem
    init = .false.
  end if
end setup

test check_local_lift
  use test_helper
  real(dp), dimension(ndof)        :: phi, residual
  real(dp), dimension(grid%numfac) :: lift
  integer :: i
continue
  !phi = one
  do i = 1, ndof
    phi(i) = real(i,dp)
  end do
  residual = zero
  call compute_local_lift(phi,lift,grid,ndof)
  call add_lift_primal_domain(1,2,5,grid,lift,residual)
  call add_flux_primal_face(1,2,5,grid,phi,lift,residual)
  call wrt(lift,"lift")
  call wrt(residual,"residual")
end test

test check_mass_matrix
  use test_helper
  real(dp), dimension(3,3) :: m_inv
continue
  m_inv = get_mass_matrix_inv(two)
  !call wrt(m_inv,"m_inv")
end test

test check_residual
  use test_helper
  use solver, only : add_flux_contributions
  real(dp), dimension(ndof) :: residual, phi
  real(dp) :: uinf, vinf
continue
  phi = one
  residual = zero
  uinf = one; vinf = zero

  !call add_lift_primal_domain(1,2,5,grid,lhspo)
  !call add_lift_second_domain(1,2,5,grid,lhspo)
  !call add_jump_second_face(1,2,5,grid,lhspo)
  !call add_flux_primal_face(1,2,5,grid,phi,residual)
  !call add_boundary_flux(residual,grid,uinf,vinf)
  !call wrt(residual,"residual")
end test

end test_suite
