test_suite gridtools

test interior_face

  use test_helper

continue

  call setup_cube

  write(*,*) "CHECK: intface = ",grid%intfac(3:4,5)
  write(*,*) "CHECK: normal  = ",grid%del(:,5)

end test

end test_suite
