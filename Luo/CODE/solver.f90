module solver

  use kinddefs,  only : dp
  use gridtools, only : gridtype, get_global_dof

  implicit none
  private

  public :: get_residual
  public :: get_soln
  public :: solve
  public :: iterate
  public :: init_freestream

  !unit tests
  public :: compute_domain_integral
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
  function init_freestream(ndof,grid,uinf,vinf) result(phi)
    use flux_functions, only : get_global_dof, get_basis
    type(gridtype),            intent(in)    :: grid
    integer,                   intent(in)    :: ndof
    real(dp),                  intent(in)    :: uinf,vinf
    real(dp), dimension(ndof) :: phi
    real(dp), dimension(3)    :: bx, by, phi_local
    real(dp) :: denom
    integer :: ielem, ip1, ip2, ip3
  continue
    phi = zero
    do ielem= 1, grid%nelem
      call get_basis(bx,by,grid,ielem)
      ip1 = get_global_dof(grid%inpoel(1,ielem),ielem,grid)
      ip2 = get_global_dof(grid%inpoel(2,ielem),ielem,grid)
      ip3 = get_global_dof(grid%inpoel(3,ielem),ielem,grid)

      ! Formulation so that gradient of phi is (uinf,vinf)
      denom = (by(1)*bx(2)-bx(1)*by(2))
      phi(ip1) = -(bx(2)*(by(3)-vinf)+by(2)*uinf-by(2)*bx(3)) / denom
      phi(ip2) = (bx(1)*(by(3)-vinf)+by(1)*uinf-by(1)*bx(3)) / denom
      phi(ip3) = one
      phi_local(1) = phi(ip1)
      phi_local(2) = phi(ip2)
      phi_local(3) = phi(ip3)

    end do
  end function init_freestream

!============================ ITERATE =======================================80
! Iterate to explicitly advance the solution in pseudo-time
!============================================================================80
  subroutine iterate(phi,grid,ndof,tolerance)

    use namelist_data, only : uinf, vinf, nsteps, dt

    type(gridtype),            intent(in)    :: grid
    integer,                   intent(in)    :: ndof
    real(dp),                  intent(in)    :: tolerance
    real(dp), dimension(ndof), intent(inout) :: phi
    real(dp), dimension(ndof)  :: residual

    real(dp) :: L2_error, rel_res, res0

    integer :: timestep

  continue

    write(*,*) "Iteration  L2_error"
    do timestep = 1, nsteps

      ! Integrate explicitly in pseudo time to evolve phi
      call rk_integrate(phi,residual,dt,grid,ndof,uinf,vinf)

      ! Compute the residual L2 error norm for reference
      L2_error = compute_L2(residual,ndof)
      if (timestep == 1) res0 = L2_error + tiny(one)
      rel_res = L2_error/res0

      write(*,11) timestep,rel_res, L2_error

      if (rel_res < tolerance) then
        write(*,10) "Solution converged to tolerance:",tolerance,   &
                    "at iteration:",timestep
        exit
      end if

    end do
10 format(A,x,g10.3,x,A,x,i5)
11 format(i10,2(x,e11.4))
  end subroutine

!============================ RK_INTEGRATE ==================================80
! Explicitly integrate in pseudo time via Runge-Kutta
!============================================================================80
  subroutine rk_integrate(phi,residual,dt,grid,ndof,uinf,vinf)
    type(gridtype),            intent(in)    :: grid
    integer,                   intent(in)    :: ndof
    real(dp),                  intent(in)    :: dt, uinf, vinf
    real(dp), dimension(ndof), intent(out)   :: residual
    real(dp), dimension(ndof), intent(inout) :: phi
    real(dp), dimension(ndof) :: phi0
    real(dp), dimension(5,5)  :: rk_coefs

    real(dp) :: coeff, rk_coef

    integer :: istag, nstag, dof
  continue

    nstag = 5
    ! Runge-Kutta coefficients
    data rk_coefs  /1.0    , 0.0    , 0.0   , 0.0  , 0.0,  &
                    0.5    , 1.0    , 0.0   , 0.0  , 0.0,  &
                    0.1919 , 0.4390 , 1.0   , 0.0  , 0.0,  &
                    0.25   , 0.3333 , 0.5   , 1.0  , 0.0,  &
                    0.25   , 0.1666 , 0.375 , 0.5  , 1.0   /

    ! Save solution
    phi0 = phi
    do istag=1,nstag
      rk_coef = rk_coefs(istag,nstag)
      !phi(1) = one ! Strongly set phi at point for constant
      residual = get_residual(phi,grid,ndof,uinf,vinf)
      do dof=1,ndof
        coeff = rk_coef*dt
        phi(dof) = phi0(dof) - coeff*residual(dof)
      end do
    end do

  end subroutine rk_integrate

!============================ GET_RESIDUAL ==================================80
! Form the residual
!============================================================================80
  function get_residual(phi,grid,ndof,uinf,vinf) result(residual)
    use flux_functions, only : compute_local_lift
    type(gridtype),            intent(in) :: grid
    integer,                   intent(in) :: ndof
    real(dp),                  intent(in) :: uinf,vinf
    real(dp), dimension(ndof), intent(in) :: phi
    real(dp), dimension(ndof)          :: residual
    real(dp), dimension(2,grid%numfac) :: lift
  continue

    residual = zero    

    call compute_local_lift(phi,lift,grid,ndof)
    !lift = zero

    call compute_domain_integral(residual,phi,grid,ndof)
    call add_flux_contributions(residual,phi,lift,grid,ndof,uinf,vinf)
    call invert_mass_matrix(residual,grid,ndof)

  end function get_residual

!======================= INVERT_MASS_MATRIX  ================================80
! Form and invert the mass matrix via cholesky solve
!============================================================================80
  subroutine invert_mass_matrix(residual,grid,ndof)

    use flux_functions, only : get_mass_matrix_inv

    type(gridtype),         intent(in)    :: grid
    integer,                intent(in)    :: ndof
    real(dp), dimension(:), intent(inout) :: residual

    real(dp), dimension(ndof) :: res

    real(dp), dimension(:,:,:), allocatable, save :: m_inv

    real(dp) :: D, row_sum
    integer :: ielem, i, j, k

  continue

    if (.not. allocated(m_inv)) then
      allocate(m_inv(3,3,grid%nelem))
      do ielem = 1, grid%nelem
        D = grid%geoel(5,ielem)
        m_inv(:,:,ielem) = get_mass_matrix_inv(d)
      end do
    end if

    ! matrix vector product of the inverse mass matrix and residual
    do ielem = 1, grid%nelem
      do i = 1, grid%nnode
        row_sum = zero
        do k = 1, grid%nnode
          j = (ielem-1)*grid%nnode + k
          row_sum = row_sum + m_inv(i,k,ielem)*residual(j)
        end do
        j = (ielem-1)*grid%nnode + i
        res(j) = row_sum
      end do
    end do

    residual = res

  end subroutine invert_mass_matrix

!======================== COMPUTE_DOMAIN_INTEGRALS ==========================80
! compute the domain integrals and add them to the residual
!============================================================================80
  subroutine compute_domain_integral(residual,phi,grid,ndof)

    use flux_functions, only : get_basis

    integer,                   intent(in)    :: ndof
    type(gridtype),            intent(in)    :: grid
    real(dp), dimension(ndof), intent(in)    :: phi
    real(dp), dimension(ndof), intent(inout) :: residual

    real(dp), dimension(3) :: bx, by

    real(dp) :: area, D

    integer :: ielem, i, ip, j, jp, global_node

  continue

    do ielem = 1, grid%nelem

      call get_basis(bx,by,grid,ielem)
      D = grid%geoel(5,ielem)
      area = half*D

      do i=1,grid%nnode
        global_node = grid%inpoel(i,ielem)
        ip = get_global_dof(global_node,ielem,grid)
        do j=1,grid%nnode
          global_node = grid%inpoel(j,ielem)
          jp = get_global_dof(global_node,ielem,grid)
          residual(ip) = residual(ip) &
                       + (bx(i)*bx(j)*phi(jp) + by(i)*by(j)*phi(jp))*area
        end do
      end do

    end do

  end subroutine compute_domain_integral

!===================== ADD_FLUX_CONTRIBUTIONS ===============================80
! Add the flux contributions to the residual 
!============================================================================80
  subroutine add_flux_contributions(residual,phi,lift,grid,ndof,uinf,vinf)

    use flux_functions, only : add_lift_primal_domain,  &
                               add_flux_primal_face,    &
                               add_boundary_flux

    integer,        intent(in) :: ndof
    type(gridtype), intent(in) :: grid
    real(dp),       intent(in) :: uinf,vinf

    real(dp), dimension(:),    intent(in)    :: phi
    real(dp), dimension(:,:),  intent(in)    :: lift
    real(dp), dimension(ndof), intent(inout) :: residual

    integer :: iface, icell, jcell

  continue

    interior_faces: do iface = grid%nface+1, grid%nface+grid%numfac

      icell = grid%intfac(1,iface)  ! left cell
      jcell = grid%intfac(2,iface)  ! right cell

      call add_lift_primal_domain(icell,jcell,iface,grid,lift,residual)
      call add_flux_primal_face(icell,jcell,iface,grid,phi,lift,residual)

    end do interior_faces

    call add_boundary_flux(residual,grid,uinf,vinf)

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

  subroutine get_soln(Vx,Vy,Vt,phi,nodal_phi,grid)

    use flux_functions, only : get_basis

    type(gridtype), intent(in) :: grid

    real(dp), dimension(:),   intent(in)  :: phi
    real(dp), dimension(:),   intent(out) :: Vx, Vy, Vt, nodal_phi

    real(dp), dimension(3)          :: phi_local, bx, by
    real(dp), dimension(grid%npoin) :: Vxarea, Vyarea, phiarea,area

    integer :: ielem, ipoin, ip1, ip2, ip3, dof1, dof2, dof3

    real(dp) :: local_area, Vx_local, Vy_local

    real(dp), parameter :: zero = 0.0_dp
    real(dp), parameter :: half = 0.5_dp

  continue

    phiarea  = zero
    Vx_local = zero
    Vy_local = zero
    Vxarea   = zero
    Vyarea   = zero
    area     = zero

    do ielem=1,grid%nelem

      call get_basis(bx,by,grid,ielem)

      ip1=grid%inpoel(1,ielem)
      ip2=grid%inpoel(2,ielem)
      ip3=grid%inpoel(3,ielem)

      dof1=get_global_dof(ip1,ielem,grid)
      dof2=get_global_dof(ip2,ielem,grid)
      dof3=get_global_dof(ip3,ielem,grid)

      local_area = half*grid%geoel(5,ielem)

      phi_local(1) = phi(dof1)*local_area
      phi_local(2) = phi(dof2)*local_area
      phi_local(3) = phi(dof3)*local_area

      Vx_local = sum(phi_local(1:3)*bx(1:3))
      Vy_local = sum(phi_local(1:3)*by(1:3))

      phiarea(ip1) = phiarea(ip1)+phi_local(1)
      phiarea(ip2) = phiarea(ip2)+phi_local(2)
      phiarea(ip3) = phiarea(ip3)+phi_local(3)

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

      nodal_phi(ipoin) = phiarea(ipoin)/area(ipoin)

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
