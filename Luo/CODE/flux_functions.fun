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
  use solver, only : init_freestream
  real(dp), dimension(ndof)        :: phi, residual
  real(dp), dimension(grid%numfac) :: lift
  real(dp) :: uinf, vinf
continue
  !phi = one
  uinf = one; vinf = zero;
  phi = init_freestream(ndof,grid,uinf,vinf)
  residual = zero
  call compute_local_lift(phi,lift,grid,ndof)
  call add_lift_primal_domain(1,2,9,grid,lift,residual)
  !call add_flux_primal_face(1,2,9,grid,phi,lift,residual)
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

end test_suite
