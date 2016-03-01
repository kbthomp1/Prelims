test_suite linalg

  integer, parameter :: dp = selected_real_kind(P=15)

test cholesky_solver

  use test_helper

  real(dp), dimension(3)   :: b, x, x_exact
  real(dp), dimension(3,3) :: A
  integer :: i, j

continue

  A = zero
  do i = 1, 3
    A(i,i) = rand(i)
  end do
  x_exact = [rand(1), rand(2), rand(3)]
  !x_exact = [ one, two, 3._dp]
  !A = zero
  !do i = 1,3
  !  A(i,i) = real(i,dp)
  !end do
  
  do i = 1,3
    b(i) = sum(A(i,:)*x_exact)
  end do
  
  call cholesky_decomp(A,3)
  !do i = 1,3
  !  write(*,20) "A = ",(A(i,j),j=1,3)
  !end do
  
  call cholesky_solve(A,b,3)
  x = b

  do i = 1, 3
    assert_equal_within(x(i),x_exact(i),1.E-15_dp)
    write(*,30) "CHECK: phi=exact =>",x(i),"=",x_exact(i)
  end do

20 format(A,300g11.4)
30 format(2(a,g17.7))
end test

test jacobi_solve

  use test_helper

  real(dp), dimension(:),   allocatable :: rhspo, phi, exact
  real(dp), dimension(:,:), allocatable :: lhspo

  real(dp) :: tolerance = 1.E-5_dp

  integer :: ndof, i, j, nsteps

continue

  call setup_cube
  ndof = 2*grid%npoin
  allocate(lhspo(ndof,ndof))
  allocate(rhspo(ndof))
  allocate(phi(ndof))
  allocate(exact(ndof))

  exact = zero
  exact(1:grid%npoin) = [rand(1),rand(2),rand(3),rand(4)]

  lhspo = zero
  do i = 1, ndof
    lhspo(i,i) = real(i,dp)
    !write(*,20) "CHECK: lhspo = ",(lhspo(i,j), j=1,ndof)
  end do

  rhspo = zero
  do i = 1, ndof
    rhspo(i) = sum(lhspo(i,:)*exact(:))
    !write(*,*) "CHECK: rhs = ",rhspo(i)
  end do

  nsteps = 4
  call jacobi(lhspo,rhspo,phi,ndof,nsteps,tolerance)

  write(*,*)
  write(*,*) "CHECK: jacobi test"
  do i = 1, ndof
    assert_equal(phi(i),exact(i))
    write(*,30) "CHECK: phi=exact =>",phi(i),"=",exact(i)
  end do

20 format(A,300g11.4)
30 format(2(a,g11.3))

end test

end test_suite
