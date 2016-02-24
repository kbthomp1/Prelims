test_suite linalg

  integer, parameter :: dp = selected_real_kind(P=15)

  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter :: half = 0.5_dp
  real(dp), parameter :: one  = 1.0_dp

test jacobi_solve

  integer, parameter :: nvars  = 1
  integer, parameter :: neqns  = 1
  integer, parameter :: npoin  = 3
  integer, parameter :: nsteps = 3
  
  real(dp), parameter :: tolerance = 1.e-20_dp
  
  real(dp), dimension(neqns,npoin) :: b
  real(dp), dimension(nvars,npoin) :: x
  
  real(dp), dimension(nvars,neqns,npoin,npoin) :: A

  real(dp), dimension(npoin,npoin) :: identity

  integer :: i, j

continue

  identity = zero
  do i = 1,npoin
    identity(i,i) = one
  end do

  do i = 1,nvars
    do j = 1,neqns
      A(i,j,:,:) = 2.0_dp*identity
    end do
  end do

  do i = 1, neqns
    do j = 1, npoin
      b(i,j) = real(i,dp)*real(j,dp)
      write(*,*) "CHECK: b = ",b(i,j)
    end do
  end do

  call jacobi(A,b,x,npoin,nvars,neqns,nsteps,tolerance)

  do i = 1, nvars
    do j = 1, npoin
      assert_equal(x(i,j),half*b(i,j))
      write(*,"(2(a,g11.3))") "CHECK: x = b => ",x(i,j), "=",b(i,j)
    end do
  end do

end test

end test_suite
