test_suite solver

  integer, parameter :: dp = selected_real_kind(P=15)

test mass_matrix_inversion

  use test_helper

  real(dp), dimension(:),   allocatable :: residual, phi

  integer :: ndof, i, j

  real(dp) :: uinf, vinf

continue

  call setup_cube
  ndof = grid%npoin
  uinf = one; vinf = zero
  allocate(residual(ndof))
  allocate(phi(ndof))

  call invert_mass_matrix(residual,grid,ndof)

end test

!test check_residual
!
!  use test_helper
!
!  real(dp), dimension(:),   allocatable :: residual, phi
!
!  integer :: ndof, i, j
!
!  real(dp) :: uinf, vinf
!
!continue
!
!  call setup_cube
!  ndof = grid%npoin
!  uinf = one; vinf = zero
!  allocate(residual(ndof))
!  allocate(phi(ndof))
!
!  do i = 1, ndof
!    phi(i) = real(i,dp)
!  end do
!  phi = one
!
!  residual = get_residual(phi,grid,ndof,uinf,vinf)
!
!  !call add_flux_contributions(lhspo,grid,ndof)
!
!  !call add_lift_primal_domain(1,2,5,grid,lhspo)
!  !call add_lift_second_domain(1,2,5,grid,lhspo)
!  !call add_jump_second_face(1,2,5,grid,lhspo)
!  !call add_flux_primal_face(1,2,5,grid,lhspo)
!
!  !do i = 1, ndof
!  !  write(*,20) "CHECK: residual = ",residual(i)
!  !end do
!
!20 format(A,300g11.4)
!
!end test

end test_suite
