test_suite gridtools

setup
  use test_helper
  logical, save :: init = .true.
  if(init) then
    call setup_cube
    ndof = grid%nnode*grid%nelem
    init = .false.
  end if
end setup

test check_local_global_map
  use test_helper
  integer :: dof
continue
  dof = get_global_dof(3,1,grid)
  assert_equal(dof,3)
end test

end test_suite
