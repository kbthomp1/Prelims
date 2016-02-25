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

          ! \nabla(B) Basis for lift operator
          lhspo(ip,jp+grid%npoin) = lhspo(ip,jp)                               &
                                  + (bx(i)*bx(j) + by(i)*by(j))*area

          ! P1^2 basis for lift operator
          !lhspo(ip,jp+grid%npoin) = fact*(bx(j) + by(j))*D**2
        end do
        lhspo(ip+grid%npoin,ip+grid%npoin) = one
      end do

    end do

    call add_flux_contributions(lhspo,grid,ndof)

  end function get_lhspo

!===================== ADD_FLUX_CONTRIBUTIONS ===============================80
! Add the flux contributions to the LHS matrix
!============================================================================80
  subroutine add_flux_contributions(lhspo,grid,ndof)

    integer,        intent(in) :: ndof
    type(gridtype), intent(in) :: grid

    real(dp), dimension(ndof,ndof), intent(inout) :: lhspo

    integer :: iface, icell, jcell
    integer :: ip1, ip2, ip3, jp1, jp2, jp3

    real(dp), dimension(2) :: bx_i, by_i, bx_j, by_j

    real(dp), dimension(grid%nnode) :: icell_coefs, jcell_coefs

    real(dp) :: nx, ny, face_area

  continue

    write(*,*) "CHECK: nface, numfac = ",grid%nface,grid%numfac

    interior_faces: do iface = grid%nface+1, grid%nface+grid%numfac

      icell = grid%intfac(1,iface)  ! left cell
      jcell = grid%intfac(2,iface)  ! right cell

      !compute unit normal vector
      nx = grid%del(1,iface)/grid%del(3,iface)
      ny = grid%del(2,iface)/grid%del(3,iface)
      face_area = grid%del(3,iface)

      write(*,*) "CHECK: i/j cell:",icell,jcell
      write(*,*) "CHECK: nx, ny:",nx,ny
      write(*,*) "CHECK: face area: ", face_area

      ip1 = grid%inpoel(1,icell)
      ip2 = grid%inpoel(2,icell)
      ip3 = grid%inpoel(3,icell)

      jp1 = grid%inpoel(1,jcell)
      jp2 = grid%inpoel(2,jcell)
      jp3 = grid%inpoel(3,jcell)

      write(*,*) "CHECK:  left cell pts:",ip1,ip2,ip3
      write(*,*) "CHECK: right cell pts:",jp1,jp2,jp3

      ! Basis functions
      bx_i(1:2) = grid%geoel(1:2,icell)
      by_i(1:2) = grid%geoel(3:4,icell)
      bx_j(1:2) = grid%geoel(1:2,jcell)
      by_j(1:2) = grid%geoel(3:4,jcell)

      write(*,*) "CHECK: bx_i:", bx_i
      write(*,*) "CHECK: bx_j:", bx_j

      ! Left cell contributions to flux
      icell_coefs(1) = bx_i(1)*nx + by_i(1)*ny
      icell_coefs(2) = bx_i(2)*nx + by_i(2)*ny
      icell_coefs(3) = -(bx_i(1) + bx_i(2))*nx   &
                       -(by_i(1) + by_i(2))*ny

      ! Right cell contributions to flux
      jcell_coefs(1) = bx_j(1)*nx + by_j(1)*ny
      jcell_coefs(2) = bx_j(2)*nx + by_j(2)*ny
      jcell_coefs(3) = -(bx_j(1) + bx_j(2))*nx   &
                       -(by_j(1) + by_j(2))*ny

      write(*,"(A,3g12.4)") "CHECK: icell_coefs: ", icell_coefs
      write(*,"(A,3g12.4)") "CHECK: jcell_coefs: ", jcell_coefs

    end do interior_faces

  end subroutine add_flux_contributions

!============================ GET_SOLN ======================================80
! Calculate velocities Vx, Vy, and Vt by interpolating solution to nodes
!
! Note that the gradients of the linear basis functions are constant across
! each element, thus we can use an area weighted average to recover the 
! the gradient at a node.  Since phi = sum(phi_i*B_i)
! 
!   phi_x = sum(phi_i*Bx_i) = phi_1*a_1 + phi_2*a_2 + phi_3*a_3
!   phi_y = sum(phi_i*By_i) = phi_1*b_1 + phi_2*b_2 + phi_3*b_3
! 
! NOTE: since B_1 + B_2 + B_3 = 1, a/b_3 = 1 - a/b_2 - a/b_1 
!============================================================================80

  subroutine get_soln(Vx,Vy,Vt,phi,grid)

    type(gridtype), intent(in) :: grid

    real(dp), dimension(:),   intent(in)  :: phi
    real(dp), dimension(:),   intent(out) :: Vx, Vy, Vt

    real(dp), dimension(grid%npoin) :: Vxarea, Vyarea, area

    integer :: ielem, ipoin, ip1, ip2, ip3

    real(dp) :: local_area, Vx_local, Vy_local

    real(dp), parameter :: zero = 0.0_dp
    real(dp), parameter :: half = 0.5_dp

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

      local_area = half*grid%geoel(5,ielem)

      Vx_local = (grid%geoel(1,ielem)*(phi(ip1)-phi(ip3))       &
                + grid%geoel(2,ielem)*(phi(ip2)-phi(ip3)))      &
                * local_area

      Vy_local = (grid%geoel(3,ielem)*(phi(ip1)-phi(ip3))       &
                + grid%geoel(4,ielem)*(phi(ip2)-phi(ip3)))      &
                * local_area

      Vxarea(ip1) = Vxarea(ip1)+Vx_local
      Vxarea(ip2) = Vxarea(ip2)+Vx_local
      Vxarea(ip3) = Vxarea(ip3)+Vx_local

      Vyarea(ip1) = Vyarea(ip1)+Vy_local
      Vyarea(ip2) = Vyarea(ip2)+Vy_local
      Vyarea(ip3) = Vyarea(ip3)+Vy_local

      area(ip1)  = area(ip1)+local_area
      area(ip2)  = area(ip2)+local_area
      area(ip3)  = area(ip3)+local_area

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
