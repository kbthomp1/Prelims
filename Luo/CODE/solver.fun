test_suite solver

  integer, parameter :: dp = selected_real_kind(P=15)

test check_lhspo

  use test_helper

  real(dp), dimension(:),   allocatable :: rhspo, phi, exact
  real(dp), dimension(:,:), allocatable :: lhspo

  real(dp) :: tolerance = 1.E-5_dp

  integer :: ndof, i, j, nsteps

continue

  call setup_cube
  ndof = 2*npoin
  allocate(lhspo(ndof,ndof))
  allocate(rhspo(ndof))
  allocate(phi(ndof))
  allocate(exact(ndof))

  lhspo = get_lhspo(ndof,nelem,nnode,npoin,inpoel,geoel)

  exact = zero
  exact(1:npoin) = [two,one,one,one]
  !exact(1:npoin) = [rand(1),rand(2),rand(3),rand(4)]
  rhspo = zero

  lhspo = zero
  do i = 1, ndof
    lhspo(i,i) = real(i,dp)
    write(*,20) "CHECK: lhspo = ",(lhspo(i,j), j=1,ndof)
  end do

  do i = 1, ndof
    rhspo(i) = sum(lhspo(i,:)*exact(:))
    write(*,*) "CHECK: rhs = ",rhspo(i)
  end do

  nsteps = 4
  call solve("point jacobi",lhspo,rhspo,phi,ndof,nsteps,tolerance)

  do i = 1, ndof
    assert_equal(phi(i),exact(i))
    write(*,30) "CHECK: phi=exact =>",phi(i),"=",exact(i)
  end do

20 format(A,300g11.4)
30 format(2(a,g11.3))

end test

end test_suite
