test_suite io_helpers

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

test test_restart
  use test_helper
  real(dp), dimension(ndof) :: phi, phi_restart
  integer :: i
continue
  do i = 1,ndof
    phi(i) = real(i,dp)
  end do
  call write_restart(phi,ndof)
  call read_restart(phi_restart,ndof)
  do i = 1,ndof
    assert_equal(phi(i),phi_restart(i))
  end do
end test

end test_suite
