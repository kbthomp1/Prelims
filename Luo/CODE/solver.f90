module solver

  use kinddefs,  only : dp
  use gridtools, only : gridtype


  implicit none
  private

  public :: get_lhspo
  public :: get_rhspo
  public :: get_soln
  public :: solve

  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter :: half = 0.5_dp
  real(dp), parameter :: one  = 1.0_dp

contains

!============================= GET_RHSPO ====================================80
! Formulate the RHS load vector with neumann BC's
!============================================================================80
  function get_rhspo(grid,ndof,uinf,vinf) result(rhspo)

    type(gridtype), intent(in) :: grid

    integer, intent(in) :: ndof

    real(dp), intent(in) :: uinf,vinf
  
    real(dp), dimension(ndof) :: rhspo

    real(dp) :: cface

    integer :: ip1, ip2, iface
    
  continue 
  
    rhspo = zero
    
    do iface=1,grid%nface
      if(grid%bcface(3,iface) == 4) then
        
        ip1 = grid%bcface(1,iface)
        ip2 = grid%bcface(2,iface)
        
        cface = 0.5_dp*(uinf*grid%rface(1,iface) + vinf*grid%rface(2,iface))
    
        rhspo(ip1) = rhspo(ip1) + cface
        rhspo(ip2) = rhspo(ip2) + cface
    
      end if
    end do
  
  end function get_rhspo

!================================ GET_LHSPO =================================80
! Form the global stiffness matrix
!============================================================================80

  function get_lhspo(grid,ndof) result (lhspo)

    integer, intent(in) :: ndof

    type(gridtype), intent(in) :: grid

    real(dp), dimension(ndof,ndof) :: lhspo

    real(dp), dimension(3) :: bx, by

    real(dp) :: area, D

    integer :: ielem, ip1, ip2,ip3
    integer :: i, ip, j, jp

    real(dp), parameter :: fact = one/12._dp

  continue

  lhspo(1:ndof,1:ndof) = 0.0

    do ielem=1,grid%nelem

      ip1   = grid%inpoel(1,ielem)
      ip2   = grid%inpoel(2,ielem)
      ip3   = grid%inpoel(3,ielem)

      bx(1) = grid%geoel(1,ielem)
      bx(2) = grid%geoel(2,ielem)
      bx(3) = -(bx(1)+bx(2))
      by(1) = grid%geoel(3,ielem)
      by(2) = grid%geoel(4,ielem)
      by(3) = -(by(1)+by(2))

      D = grid%geoel(5,ielem)
      area = half*D

      do i=1,grid%nnode
        ip = grid%inpoel(i,ielem)
        do j=1,grid%nnode
          jp = grid%inpoel(j,ielem)
          lhspo(ip,jp) = lhspo(ip,jp) + (bx(i)*bx(j) + by(i)*by(j))*area
          lhspo(ip,jp+grid%npoin) = fact*(bx(j) + by(j))*D**2
        end do
        lhspo(ip+grid%npoin,ip+grid%npoin) = one
      end do

    end do

  end function get_lhspo

!============================ GET_SOLN ======================================80
! Calculate velocities Vx, Vy, and Vt by interpolating solution to nodes
!============================================================================80

  subroutine get_soln(Vx,Vy,Vt,phi,grid)

    type(gridtype), intent(in) :: grid

    real(dp), dimension(:),   intent(in)  :: phi
    real(dp), dimension(:),   intent(out) :: Vx, Vy, Vt

    real(dp), dimension(grid%nelem) :: Vx_local, Vy_local
    real(dp), dimension(grid%npoin) :: Vxarea, Vyarea, area

    integer :: ielem, ipoin, ip1, ip2, ip3

    real(dp), parameter :: zero = 0.0_dp

  continue

    Vx_local = zero
    Vy_local = zero
    Vxarea   = zero
    Vyarea   = zero
    area     = zero

    do ielem=1,grid%nelem

      ip1=grid%inpoel(1,ielem)
      ip2=grid%inpoel(2,ielem)
      ip3=grid%inpoel(3,ielem)

      Vx_local(ielem) = (grid%geoel(1,ielem)*(phi(ip1)-phi(ip3))       &
                       + grid%geoel(2,ielem)*(phi(ip2)-phi(ip3)))      &
                       * grid%geoel(5,ielem)

      Vy_local(ielem) = (grid%geoel(3,ielem)*(phi(ip1)-phi(ip3))       &
                       + grid%geoel(4,ielem)*(phi(ip2)-phi(ip3)))      &
                       * grid%geoel(5,ielem)

      Vxarea(ip1) = Vxarea(ip1)+Vx_local(ielem)
      Vxarea(ip2) = Vxarea(ip2)+Vx_local(ielem)
      Vxarea(ip3) = Vxarea(ip3)+Vx_local(ielem)

      Vyarea(ip1) = Vyarea(ip1)+Vy_local(ielem)
      Vyarea(ip2) = Vyarea(ip2)+Vy_local(ielem)
      Vyarea(ip3) = Vyarea(ip3)+Vy_local(ielem)

      area(ip1)  = area(ip1)+grid%geoel(5,ielem)
      area(ip2)  = area(ip2)+grid%geoel(5,ielem)
      area(ip3)  = area(ip3)+grid%geoel(5,ielem)

    end do

    do ipoin=1,grid%npoin

      Vx(ipoin) = Vxarea(ipoin)/area(ipoin)
      Vy(ipoin) = Vyarea(ipoin)/area(ipoin)

      Vt(ipoin) = sqrt(Vx(ipoin)**2 + Vy(ipoin)**2)

    end do
  end subroutine get_soln

!============================== SOLVE =======================================80
! Solve the linear system per user's specified method
!============================================================================80
  subroutine solve(lin_solver,lhspo,rhspo,phi,ndof,nsteps,tolerance)

    use linalg, only : gauss_seidel,jacobi,conjgrad

    character(len=*), intent(in) :: lin_solver

    integer, intent(in) :: nsteps, ndof

    real(dp), dimension(ndof),      intent(in)    :: rhspo
    real(dp), dimension(ndof),      intent(inout) :: phi
    real(dp), dimension(ndof,ndof), intent(in)    :: lhspo

    real(dp), intent(in) :: tolerance

  continue

    select case(lin_solver)
    case("conjugate gradient")
      call conjgrad(lhspo,rhspo,phi,ndof,nsteps)
    case("point jacobi")
      call jacobi(lhspo,rhspo,phi,ndof,nsteps,tolerance)
    case("gauss seidel")
      call gauss_seidel(lhspo,rhspo,phi,ndof,nsteps,tolerance)
    case default
      write(*,*) "Error: linear solver specified not available:",lin_solver
      stop 1
    end select

    !phi(1,:) = x

  end subroutine solve

end module solver
