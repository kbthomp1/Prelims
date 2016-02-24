test_suite solver

  integer, parameter :: dp = selected_real_kind(P=15)

  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter :: one  = 1.0_dp

test check_lhspo

  use gridtools, only : basis_function

  integer, parameter :: nnode  = 3
  integer, parameter :: nvars  = 2
  integer, parameter :: neqns  = 2
  integer, parameter :: nelem  = 2
  integer, parameter :: npoin  = 6
  integer, parameter :: ndimn  = 2
  integer, parameter :: nsteps = 3
  
  real(dp), parameter :: tolerance = 1.e-20_dp

  real(dp), dimension(:,:),     allocatable :: coord, geoel, rface
  real(dp), dimension(:,:),     allocatable :: rhspo, phi
  real(dp), dimension(:,:,:,:), allocatable :: lhspo

  integer, dimension(nnode,npoin) :: inpoel

continue

  allocate(inpoel(3,nelem))
  allocate(coord(ndimn,npoin))

  call basis_function(nelem,geoel,inpoel,coord)

  lhspo = get_lhspo(npoin,nelem,nnode,nvars,neqns,inpoel,geoel)

end test

end test_suite
