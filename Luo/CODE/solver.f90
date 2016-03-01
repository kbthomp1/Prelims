module solver

  use kinddefs,  only : dp
  use gridtools, only : gridtype

  implicit none
  private

  public :: get_residual
  public :: get_soln
  public :: solve
  public :: iterate
  public :: init_freestream

  !unit tests
  public :: add_flux_contributions
  public :: invert_mass_matrix

  real(dp), parameter :: my_4th = 0.25_dp

  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter :: half = 0.5_dp
  real(dp), parameter :: one  = 1.0_dp
  real(dp), parameter :: two  = 1.0_dp

contains

!======================= INIT_FREESTREAM ====================================80
! Initialize phi to freestream
!============================================================================80
  function init_freestream(ndof) result(phi)
    integer,                   intent(in)    :: ndof
    real(dp), dimension(ndof) :: phi
    integer :: i
  continue
    do i = 1, ndof
      phi(i) = one
    end do
  end function init_freestream

!============================ ITERATE =======================================80
! Iterate to explicitly advance the solution in pseudo-time
!============================================================================80
  subroutine iterate(phi,grid,ndof)

    use namelist_data, only : uinf, vinf, nsteps, dt

    type(gridtype),            intent(in)    :: grid
    integer,                   intent(in)    :: ndof
    real(dp), dimension(ndof), intent(inout) :: phi
    real(dp), dimension(ndof)  :: residual

    real(dp) :: L2_error

    integer :: timestep, i

  continue

    do timestep = 1, nsteps
      residual = get_residual(phi,grid,ndof,uinf,vinf)
      L2_error = compute_L2(residual,ndof)
      do i = 1, ndof
        phi(i) = phi(i) - dt*residual(i)
      end do

      write(*,*) "CHECK: iter error = ",timestep,L2_error
    end do

  end subroutine

!============================ GET_RESIDUAL ==================================80
! Form the residual
!============================================================================80
  function get_residual(phi,grid,ndof,uinf,vinf) result(residual)
    type(gridtype),            intent(in) :: grid
    integer,                   intent(in) :: ndof
    real(dp),                  intent(in) :: uinf,vinf
    real(dp), dimension(ndof), intent(in) :: phi
    real(dp), dimension(ndof)  :: residual, lift
  continue

    residual = zero    

    call compute_local_lift(phi,lift,grid,ndof)

    call compute_domain_integral(residual,phi,grid,ndof)
    call add_flux_contributions(residual,phi,grid,ndof)
    call add_boundary_flux(residual,grid,uinf,vinf)
    call invert_mass_matrix(residual,grid,ndof)

  end function get_residual

!======================= INVERT_MASS_MATRIX  ================================80
! Form and invert the mass matrix via cholesky solve
!============================================================================80
  subroutine invert_mass_matrix(residual,grid,ndof)

    use linalg,         only : cholesky_decomp, cholesky_solve
    use flux_functions, only : get_mass_matrix

    type(gridtype),         intent(in)    :: grid
    integer,                intent(in)    :: ndof
    real(dp), dimension(:), intent(inout) :: residual

    real(dp), dimension(3,3)       :: m
    real(dp), dimension(ndof,ndof) :: m_global

    real(dp), dimension(:,:), allocatable, save :: m_inv

    real(dp) :: D
    integer :: ip1, ip2, ip3, ip, jp, ielem
    integer :: i, j

  continue

    if (.not. allocated(m_inv)) then
      allocate(m_inv(ndof,ndof))
      m_global = zero

      do ielem = 1, grid%nelem

        ip1   = grid%inpoel(1,ielem)
        ip2   = grid%inpoel(2,ielem)
        ip3   = grid%inpoel(3,ielem)

        D     = grid%geoel(5,ielem)

        m = get_mass_matrix(D)

        do i=1,grid%nnode
          ip = grid%inpoel(i,ielem)
          do j=1,grid%nnode
            jp = grid%inpoel(j,ielem)
            m_global(ip,jp) = m_global(ip,jp) + m(i,j)
          end do
        end do

      end do

      m_inv = m_global
      call cholesky_decomp(m_inv,ndof)
    end if

    call cholesky_solve(m_inv,residual,ndof)

  end subroutine invert_mass_matrix

!======================= ADD_BOUNDARY_FLUX ==================================80
! Add flux from Neumann BC
!============================================================================80
  subroutine add_boundary_flux(residual,grid,uinf,vinf)

    type(gridtype),         intent(in)    :: grid
    real(dp),               intent(in)    :: uinf,vinf
    real(dp), dimension(:), intent(inout) :: residual

    real(dp) :: cface
    integer :: iface, ip1, ip2

  continue

    do iface=1,grid%nface
      if(grid%bcface(3,iface) == 4) then
        
        ip1 = grid%bcface(1,iface)
        ip2 = grid%bcface(2,iface)
        
        cface = 0.5_dp*(uinf*grid%rface(1,iface) + vinf*grid%rface(2,iface))
    
        residual(ip1) = residual(ip1) + cface
        residual(ip2) = residual(ip2) + cface
    
      end if
    end do

  end subroutine add_boundary_flux

!======================== COMPUTE_DOMAIN_INTEGRALS ==========================80
! compute the domain integrals and add them to the residual
!============================================================================80
  subroutine compute_domain_integral(residual,phi,grid,ndof)

    integer,                   intent(in)    :: ndof
    type(gridtype),            intent(in)    :: grid
    real(dp), dimension(ndof), intent(in)    :: phi
    real(dp), dimension(ndof), intent(inout) :: residual

    real(dp), dimension(3) :: bx, by

    real(dp) :: area, D

    integer :: ielem, ip1, ip2,ip3
    integer :: i, ip, j, jp

  continue

    do ielem = 1, grid%nelem

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
          residual(ip) = residual(ip) &
                       + (bx(i)*bx(j)*phi(j) + by(i)*by(j)*phi(j))*area
        end do
      end do

    end do

  end subroutine compute_domain_integral

!===================== COMPUTE_LOCAL_LIFT ===================================80
! Compute the local lifting function at interior faces
!============================================================================80
  subroutine compute_local_lift(phi,lift,grid,ndof)

    integer,        intent(in) :: ndof
    type(gridtype), intent(in) :: grid

    real(dp), dimension(ndof), intent(in)  :: phi
    real(dp), dimension(ndof), intent(out) :: lift

    integer :: iface, icell, jcell
  continue

    interior_faces: do iface = grid%nface+1, grid%nface+grid%numfac

      i = iface - grid%nface
      icell = grid%intfac(1,iface)  ! left cell
      jcell = grid%intfac(2,iface)  ! right cell

      lift(i) = 

    end do interior_faces

  end subroutine compute_local_lift

!===================== ADD_FLUX_CONTRIBUTIONS ===============================80
! Add the flux contributions to the residual 
!============================================================================80
  subroutine add_flux_contributions(residual,phi,grid,ndof)

    use flux_functions, only : add_lift_primal_domain,  &
                               add_lift_second_domain,  & 
                               add_jump_second_face,    &
                               add_flux_primal_face

    integer,        intent(in) :: ndof
    type(gridtype), intent(in) :: grid

    real(dp), dimension(ndof), intent(in)    :: phi
    real(dp), dimension(ndof), intent(inout) :: residual

    integer :: iface, icell, jcell

  continue

    interior_faces: do iface = grid%nface+1, grid%nface+grid%numfac

      icell = grid%intfac(1,iface)  ! left cell
      jcell = grid%intfac(2,iface)  ! right cell

      !call add_lift_primal_domain(icell,jcell,iface,grid,lhspo)
      !call add_lift_second_domain(icell,jcell,iface,grid,lhspo)
      !call add_jump_second_face(icell,jcell,iface,grid,lhspo)
      call add_flux_primal_face(icell,jcell,iface,grid,phi,residual)

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
  subroutine solve(lin_solver,lhs,rhs,phi,ndof,nsteps,tolerance)

    use linalg, only : gauss_seidel,jacobi,conjgrad

    character(len=*), intent(in) :: lin_solver

    integer, intent(in) :: nsteps, ndof

    real(dp), dimension(ndof),      intent(in)    :: rhs
    real(dp), dimension(ndof),      intent(inout) :: phi
    real(dp), dimension(ndof,ndof), intent(in)    :: lhs

    real(dp), intent(in) :: tolerance

  continue

    select case(lin_solver)
    case("conjugate gradient")
      call conjgrad(lhs,rhs,phi,ndof,nsteps)
    case("point jacobi")
      call jacobi(lhs,rhs,phi,ndof,nsteps,tolerance)
    case("gauss seidel")
      call gauss_seidel(lhs,rhs,phi,ndof,nsteps,tolerance)
    case default
      write(*,*) "Error: linear solver specified not available:",lin_solver
      stop 1
    end select

  end subroutine solve

!============================ COMPUTE_L2 ====================================80
! Compute the L2 error norm of the residual
!============================================================================80
  function compute_L2(residual,ndof) result(L2)
    integer, intent(in) :: ndof
    real(dp), dimension(ndof), intent(in) :: residual
    real(dp) :: L2
    integer :: i
  continue
    L2 = zero
    do i = 1, ndof
      L2 = L2 + residual(i)**2
    end do
    L2 = sqrt(L2)
  end function compute_L2

end module solver
